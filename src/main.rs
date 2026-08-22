// nanofem: a headless 3D finite element electromagnetic field solver in a
// single source file, hard budget 1000 lines of code, std only.
//
// Algorithm: time harmonic (exp(+j w t)) curl-curl equation for the electric
// field, discretized with first order Nedelec (Whitney) edge elements on
// tetrahedra. Boundary conditions: PEC as eliminated unknowns, first order
// absorbing boundary, natural PMC. Lumped rectangular ports are modeled as
// impedance sheets with an impressed surface current (Norton equivalent);
// scattering parameters follow from the port voltages. Materials carry
// relative permittivity, loss tangent and permeability per mesh region.
// Input: a Gmsh .msh version 2.2 ASCII mesh with physical groups plus a deck
// that maps group names to materials, boundaries and ports. Output:
// Touchstone data on stdout. Linear algebra: reverse Cuthill-McKee ordering
// and a complex symmetric skyline LDL^T factorization.

use std::collections::HashMap;

const C0: f64 = 299792458.0;
const ETA0: f64 = 376.73031346177066;

fn die(msg: &str) -> ! {
    eprintln!("nanofem: {}", msg);
    std::process::exit(1);
}

// stdout writer that exits quietly when the pipe is closed (e.g. piped to head)
fn outln(args: std::fmt::Arguments) {
    use std::io::Write;
    if writeln!(std::io::stdout().lock(), "{}", args).is_err() {
        std::process::exit(0);
    }
}

// ---------------------------------------------------------------- complex

#[derive(Clone, Copy)]
struct Cx {
    re: f64,
    im: f64,
}

fn cx(re: f64, im: f64) -> Cx {
    Cx { re, im }
}

impl std::ops::Add for Cx {
    type Output = Cx;
    fn add(self, o: Cx) -> Cx { cx(self.re + o.re, self.im + o.im) }
}

impl std::ops::Sub for Cx {
    type Output = Cx;
    fn sub(self, o: Cx) -> Cx { cx(self.re - o.re, self.im - o.im) }
}

impl std::ops::Mul for Cx {
    type Output = Cx;
    fn mul(self, o: Cx) -> Cx { cx(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
}

impl std::ops::Div for Cx {
    type Output = Cx;
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        cx((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}

impl Cx {
    fn mag(self) -> f64 { self.re.hypot(self.im) }
    fn rs(self, s: f64) -> Cx { cx(self.re * s, self.im * s) }
    // j * s * self
    fn js(self, s: f64) -> Cx { cx(-self.im * s, self.re * s) }
    fn sqrt(self) -> Cx {
        let (r, t) = (self.mag().sqrt(), self.im.atan2(self.re) / 2.0);
        cx(r * t.cos(), r * t.sin())
    }
}

// ---------------------------------------------------------------- vectors

type V3 = [f64; 3];

fn sub3(a: V3, b: V3) -> V3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
fn dot3(a: V3, b: V3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn cross3(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn sc3(a: V3, s: f64) -> V3 { [a[0] * s, a[1] * s, a[2] * s] }

// ---------------------------------------------------------------- mesh

// Gmsh .msh 2.2 ASCII. Elements keep only triangles (type 2) and tetrahedra
// (type 4); the first element tag is the physical group. Group names come
// from $PhysicalNames keyed by (dimension, id), unnamed groups get the id as
// name. Elements without tags land in group "0".
struct Mesh {
    nodes: Vec<V3>,
    tets: Vec<([usize; 4], usize)>,
    tris: Vec<([usize; 3], usize)>,
    names: Vec<String>,
}

fn parse_msh(path: &str) -> Mesh {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("cannot read {}: {}", path, e)));
    let mut m = Mesh { nodes: vec![], tets: vec![], tris: vec![], names: vec![] };
    let mut phys: HashMap<(i64, i64), String> = HashMap::new();
    let mut idmap: HashMap<i64, usize> = HashMap::new();
    let mut groups: HashMap<String, usize> = HashMap::new();
    let ls: Vec<&str> = txt.lines().collect();
    let mut li = 0;
    fn next<'a>(ls: &[&'a str], li: &mut usize) -> &'a str {
        if *li >= ls.len() {
            die("mesh file ends inside a section");
        }
        *li += 1;
        ls[*li - 1]
    }
    let count = |li: &mut usize| -> usize {
        next(&ls, li).trim().parse().unwrap_or_else(|_| die("bad section count"))
    };
    while li < ls.len() {
        let l = next(&ls, &mut li).trim().to_string();
        match l.as_str() {
            "$MeshFormat" => {
                if !next(&ls, &mut li).trim_start().starts_with("2.2") {
                    die("mesh must be Gmsh ASCII v2.2 (export with -format msh22)");
                }
            }
            "$PhysicalNames" => {
                for _ in 0..count(&mut li) {
                    let l = next(&ls, &mut li).to_string();
                    let t: Vec<&str> = l.split_whitespace().collect();
                    let (dim, id): (i64, i64) = (t[0].parse().unwrap_or(-1), t[1].parse().unwrap_or(-1));
                    phys.insert((dim, id), t[2..].join(" ").trim_matches('"').to_string());
                }
            }
            "$Nodes" => {
                for _ in 0..count(&mut li) {
                    let l = next(&ls, &mut li).to_string();
                    let t: Vec<f64> = l.split_whitespace().map(|s| s.parse().unwrap_or(f64::NAN)).collect();
                    idmap.insert(t[0] as i64, m.nodes.len());
                    m.nodes.push([t[1], t[2], t[3]]);
                }
            }
            "$Elements" => {
                for _ in 0..count(&mut li) {
                    let l = next(&ls, &mut li).to_string();
                    let t: Vec<i64> = l.split_whitespace().map(|s| s.parse().unwrap_or(-1)).collect();
                    let (typ, ntag) = (t[1], t[2] as usize);
                    let dim = match typ {
                        2 => 2,
                        4 => 3,
                        _ => continue,
                    };
                    let pid = if ntag >= 1 { t[3] } else { 0 };
                    let name = phys.get(&(dim, pid)).cloned().unwrap_or_else(|| pid.to_string());
                    let gi = *groups.entry(name.clone()).or_insert_with(|| {
                        m.names.push(name);
                        m.names.len() - 1
                    });
                    let nd: Vec<usize> = t[3 + ntag..]
                        .iter()
                        .map(|id| *idmap.get(id).unwrap_or_else(|| die("element references unknown node")))
                        .collect();
                    if dim == 2 {
                        m.tris.push(([nd[0], nd[1], nd[2]], gi));
                    } else {
                        m.tets.push(([nd[0], nd[1], nd[2], nd[3]], gi));
                    }
                }
            }
            _ => {}
        }
    }
    if m.tets.is_empty() {
        die("mesh has no tetrahedra");
    }
    m
}

// ---------------------------------------------------------------- deck

// Deck cards, one per line, * starts a comment:
//   mesh <path>                          (relative to the deck file)
//   mat <group> eps <er> [tand <d>] [mur <mr>]
//   pec <group> [<group> ...]
//   abc <group> [<group> ...]
//   port <n> <group> <jx> <jy> <jz> <z0> (n counting from 1, j = voltage direction)
//   sweep lin <f0> <f1> <npoints>
struct Mat {
    eps: f64,
    tand: f64,
    mur: f64,
}

struct PortDef {
    group: String,
    dir: V3,
    z0: f64,
}

struct Deck {
    mesh: String,
    mats: Vec<(String, Mat)>,
    pec: Vec<String>,
    abc: Vec<String>,
    ports: Vec<PortDef>,
    f0: f64,
    f1: f64,
    nf: usize,
}

fn num(s: &str) -> f64 {
    s.parse().unwrap_or_else(|_| die(&format!("bad number '{}'", s)))
}

fn parse_deck(path: &str) -> Deck {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("cannot read {}: {}", path, e)));
    let mut d = Deck { mesh: String::new(), mats: vec![], pec: vec![], abc: vec![], ports: vec![], f0: 0.0, f1: 0.0, nf: 0 };
    for l in txt.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.is_empty() || t[0].starts_with('*') {
            continue;
        }
        let bad = || -> ! { die(&format!("bad card: {}", l)) };
        match t[0].to_lowercase().as_str() {
            "mesh" if t.len() == 2 => {
                let base = std::path::Path::new(path).parent().unwrap_or(std::path::Path::new(""));
                d.mesh = base.join(t[1]).to_string_lossy().into_owned();
            }
            "mat" if t.len() >= 4 => {
                let mut mt = Mat { eps: 1.0, tand: 0.0, mur: 1.0 };
                for kv in t[2..].chunks(2) {
                    match kv[0].to_lowercase().as_str() {
                        "eps" => mt.eps = num(kv[1]),
                        "tand" => mt.tand = num(kv[1]),
                        "mur" => mt.mur = num(kv[1]),
                        _ => bad(),
                    }
                }
                d.mats.push((t[1].to_string(), mt));
            }
            "pec" if t.len() >= 2 => d.pec.extend(t[1..].iter().map(|s| s.to_string())),
            "abc" if t.len() >= 2 => d.abc.extend(t[1..].iter().map(|s| s.to_string())),
            "port" if t.len() == 7 => {
                if num(t[1]) as usize != d.ports.len() + 1 {
                    die("ports must be numbered 1, 2, ... in order");
                }
                let dir = [num(t[3]), num(t[4]), num(t[5])];
                let n = dot3(dir, dir).sqrt();
                if n == 0.0 {
                    die("port direction must be nonzero");
                }
                d.ports.push(PortDef { group: t[2].to_string(), dir: sc3(dir, 1.0 / n), z0: num(t[6]) });
            }
            "sweep" if t.len() == 5 && t[1] == "lin" => {
                d.f0 = num(t[2]);
                d.f1 = num(t[3]);
                d.nf = num(t[4]) as usize;
            }
            _ => bad(),
        }
    }
    if d.mesh.is_empty() || d.ports.is_empty() || d.nf == 0 {
        die("deck needs at least mesh, one port and a sweep");
    }
    d
}

// ---------------------------------------------------------------- ordering

// Fill reducing ordering by geometric nested dissection: split the dof set
// at the median of the longest bounding box axis, take as separator the
// second half vertices with a neighbor in the first half, order both halves
// recursively, the separator last. Returns perm with perm[old] = new.
fn ndorder(n: usize, keys: &[(u32, u32)], pts: &[V3]) -> Vec<usize> {
    let mut adj: Vec<Vec<u32>> = vec![vec![]; n];
    for &(a, b) in keys {
        if a != b {
            adj[a as usize].push(b);
            adj[b as usize].push(a);
        }
    }
    fn rec(set: &mut Vec<u32>, adj: &[Vec<u32>], pts: &[V3], order: &mut Vec<u32>, mark: &mut [u64], stamp: &mut u64) {
        if set.len() <= 32 {
            order.append(set);
            return;
        }
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for &v in set.iter() {
            for k in 0..3 {
                lo[k] = lo[k].min(pts[v as usize][k]);
                hi[k] = hi[k].max(pts[v as usize][k]);
            }
        }
        let mut ax = 0;
        for k in 1..3 {
            if hi[k] - lo[k] > hi[ax] - lo[ax] {
                ax = k;
            }
        }
        set.sort_unstable_by(|&a, &b| pts[a as usize][ax].total_cmp(&pts[b as usize][ax]));
        let half = set.len() / 2;
        *stamp += 1;
        let sa = *stamp;
        for &v in &set[..half] {
            mark[v as usize] = sa;
        }
        let mut a: Vec<u32> = set[..half].to_vec();
        let (mut b, mut sep) = (vec![], vec![]);
        for &v in &set[half..] {
            if adj[v as usize].iter().any(|&w| mark[w as usize] == sa) {
                sep.push(v);
            } else {
                b.push(v);
            }
        }
        rec(&mut a, adj, pts, order, mark, stamp);
        rec(&mut b, adj, pts, order, mark, stamp);
        order.append(&mut sep);
    }
    let mut order = Vec::with_capacity(n);
    let (mut mark, mut stamp) = (vec![0u64; n], 0u64);
    let mut set: Vec<u32> = (0..n as u32).collect();
    rec(&mut set, &adj, pts, &mut order, &mut mark, &mut stamp);
    let mut perm = vec![0; n];
    for (new, &old) in order.iter().enumerate() {
        perm[old as usize] = new;
    }
    perm
}

// ---------------------------------------------------------------- sparse ldl

// Complex symmetric sparse LDL^T, up looking, no pivoting, after Davis' LDL.
// The permuted upper triangle sits in CSC (ap, ai), unsorted within columns.
// Symbolic pass: elimination tree and exact column counts of L. Numeric
// pass: for column k the reach through the etree gives the pattern of row k
// of L in topological order, a dense accumulator carries the values.
struct Sym {
    parent: Vec<usize>,
    lp: Vec<usize>,
}

fn symbolic(n: usize, ap: &[usize], ai: &[u32]) -> Sym {
    let (mut parent, mut lnz, mut flag) = (vec![usize::MAX; n], vec![0usize; n], vec![usize::MAX; n]);
    for k in 0..n {
        flag[k] = k;
        for p in ap[k]..ap[k + 1] {
            let mut i = ai[p] as usize;
            while i < k && flag[i] != k {
                if parent[i] == usize::MAX {
                    parent[i] = k;
                }
                lnz[i] += 1;
                flag[i] = k;
                i = parent[i];
            }
        }
    }
    let mut lp = vec![0usize; n + 1];
    for k in 0..n {
        lp[k + 1] = lp[k] + lnz[k];
    }
    Sym { parent, lp }
}

struct Work {
    li: Vec<u32>,
    lx: Vec<Cx>,
    d: Vec<Cx>,
    y: Vec<Cx>,
    flag: Vec<usize>,
    pat: Vec<usize>,
    lnz: Vec<usize>,
}

fn work(n: usize, sym: &Sym) -> Work {
    Work {
        li: vec![0; sym.lp[n]],
        lx: vec![cx(0.0, 0.0); sym.lp[n]],
        d: vec![cx(0.0, 0.0); n],
        y: vec![cx(0.0, 0.0); n],
        flag: vec![usize::MAX; n],
        pat: vec![0; n],
        lnz: vec![0; n],
    }
}

fn numeric(n: usize, ap: &[usize], ai: &[u32], ax: &[Cx], sym: &Sym, w: &mut Work) {
    for k in 0..n {
        let mut top = n;
        w.flag[k] = k;
        w.lnz[k] = 0;
        for p in ap[k]..ap[k + 1] {
            let mut i = ai[p] as usize;
            w.y[i] = w.y[i] + ax[p];
            let mut len = 0;
            while i < k && w.flag[i] != k {
                w.pat[len] = i;
                len += 1;
                w.flag[i] = k;
                i = sym.parent[i];
            }
            while len > 0 {
                len -= 1;
                top -= 1;
                w.pat[top] = w.pat[len];
            }
        }
        w.d[k] = w.y[k];
        w.y[k] = cx(0.0, 0.0);
        for t in top..n {
            let i = w.pat[t];
            let yi = w.y[i];
            w.y[i] = cx(0.0, 0.0);
            for p in sym.lp[i]..sym.lp[i] + w.lnz[i] {
                w.y[w.li[p] as usize] = w.y[w.li[p] as usize] - w.lx[p] * yi;
            }
            let lki = yi / w.d[i];
            w.d[k] = w.d[k] - lki * yi;
            let p = sym.lp[i] + w.lnz[i];
            w.li[p] = k as u32;
            w.lx[p] = lki;
            w.lnz[i] += 1;
        }
        if w.d[k].mag() < 1e-300 {
            die("singular system matrix");
        }
    }
}

fn ldsolve(n: usize, sym: &Sym, w: &Work, b: &mut [Cx]) {
    for j in 0..n {
        let xj = b[j];
        for p in sym.lp[j]..sym.lp[j + 1] {
            b[w.li[p] as usize] = b[w.li[p] as usize] - w.lx[p] * xj;
        }
    }
    for j in 0..n {
        b[j] = b[j] / w.d[j];
    }
    for j in (0..n).rev() {
        let mut s = b[j];
        for p in sym.lp[j]..sym.lp[j + 1] {
            s = s - w.lx[p] * b[w.li[p] as usize];
        }
        b[j] = s;
    }
}

// ---------------------------------------------------------------- fem

const TE: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
const SE: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];

struct PortDat {
    z0: f64,
    w: f64,
    h: f64,
    area: f64,
    // dof -> integral of (dir . whitney) over the port face, oriented
    exc: HashMap<usize, f64>,
}

fn simulate(deck: &Deck, mesh: &Mesh) {
    let ng = mesh.names.len();
    let gidx = |n: &str| -> usize {
        mesh.names.iter().position(|x| x == n).unwrap_or_else(|| die(&format!("unknown physical group '{}'", n)))
    };
    // roles: 0 none, 1 pec, 2 abc; ports via pmap
    let mut role = vec![0u8; ng];
    for n in &deck.pec {
        role[gidx(n)] = 1;
    }
    for n in &deck.abc {
        role[gidx(n)] = 2;
    }
    let mut pmap: HashMap<usize, usize> = HashMap::new();
    for (p, pd) in deck.ports.iter().enumerate() {
        pmap.insert(gidx(&pd.group), p);
    }
    // materials per volume group: (complex eps, mur)
    let mut matg = vec![(cx(1.0, 0.0), 1.0); ng];
    for (name, mt) in &deck.mats {
        matg[gidx(name)] = (cx(mt.eps, -mt.eps * mt.tand), mt.mur);
    }
    // global edges, oriented low node -> high node
    let mut emap: HashMap<(u32, u32), u32> = HashMap::new();
    let mut enodes: Vec<(u32, u32)> = vec![];
    for (t, _) in &mesh.tets {
        for (a, b) in TE {
            let k = (t[a].min(t[b]) as u32, t[a].max(t[b]) as u32);
            if !emap.contains_key(&k) {
                emap.insert(k, enodes.len() as u32);
                enodes.push(k);
            }
        }
    }
    let ne = emap.len();
    let edge = |a: usize, b: usize| -> (usize, f64) {
        let k = (a.min(b) as u32, a.max(b) as u32);
        match emap.get(&k) {
            Some(&e) => (e as usize, if a < b { 1.0 } else { -1.0 }),
            None => die("boundary triangle does not match the volume mesh"),
        }
    };
    // PEC: eliminate all edges of pec triangles
    let mut cons = vec![false; ne];
    for (t, g) in &mesh.tris {
        if role[*g] == 1 {
            for (a, b) in SE {
                cons[edge(t[a], t[b]).0] = true;
            }
        }
    }
    let mut dof = vec![usize::MAX; ne];
    let mut ndof = 0;
    for e in 0..ne {
        if !cons[e] {
            dof[e] = ndof;
            ndof += 1;
        }
    }
    // boundary face -> adjacent tet, for the abc material weight
    let mut fmap: HashMap<[u32; 3], usize> = HashMap::new();
    for (ti, (t, _)) in mesh.tets.iter().enumerate() {
        for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut k = [t[f[0]] as u32, t[f[1]] as u32, t[f[2]] as u32];
            k.sort();
            fmap.insert(k, ti);
        }
    }
    // port geometry: total area, extent h along dir, then width w = area / h
    let mut ports: Vec<PortDat> = deck
        .ports
        .iter()
        .map(|pd| PortDat { z0: pd.z0, w: 0.0, h: 0.0, area: 0.0, exc: HashMap::new() })
        .collect();
    let mut ext = vec![(f64::MAX, f64::MIN); ports.len()];
    for (t, g) in &mesh.tris {
        if let Some(&p) = pmap.get(g) {
            let (q0, q1, q2) = (mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]]);
            let n = cross3(sub3(q1, q0), sub3(q2, q0));
            ports[p].area += 0.5 * dot3(n, n).sqrt();
            for n in t {
                let x = dot3(mesh.nodes[*n], deck.ports[p].dir);
                ext[p].0 = ext[p].0.min(x);
                ext[p].1 = ext[p].1.max(x);
            }
        }
    }
    for (p, pt) in ports.iter_mut().enumerate() {
        pt.h = ext[p].1 - ext[p].0;
        if pt.area <= 0.0 || pt.h <= 0.0 {
            die(&format!("port {} has no faces or zero extent along its direction", p + 1));
        }
        pt.w = pt.area / pt.h;
    }
    // assembly: per matrix entry (upper, dof indices) keep the frequency
    // independent parts (curl stiffness s, eps weighted mass m, boundary b)
    // and combine per frequency as s - k0^2 m + j k0 b
    let mut acc: HashMap<(u32, u32), (f64, Cx, Cx)> = HashMap::new();
    let mut add = |di: usize, dj: usize, s: f64, m: Cx, b: Cx| {
        let k = (di.min(dj) as u32, di.max(dj) as u32);
        let e = acc.entry(k).or_insert((0.0, cx(0.0, 0.0), cx(0.0, 0.0)));
        e.0 += s;
        e.1 = e.1 + m;
        e.2 = e.2 + b;
    };
    for (t, g) in &mesh.tets {
        let p0 = mesh.nodes[t[0]];
        let (r1, r2, r3) = (sub3(mesh.nodes[t[1]], p0), sub3(mesh.nodes[t[2]], p0), sub3(mesh.nodes[t[3]], p0));
        let det = dot3(r1, cross3(r2, r3));
        if det == 0.0 {
            die("degenerate tetrahedron");
        }
        let vol = det.abs() / 6.0;
        let g1 = sc3(cross3(r2, r3), 1.0 / det);
        let g2 = sc3(cross3(r3, r1), 1.0 / det);
        let g3 = sc3(cross3(r1, r2), 1.0 / det);
        let gg = [sub3(sub3(sub3([0.0; 3], g1), g2), g3), g1, g2, g3];
        let (epsc, mur) = matg[*g];
        let ii = |p: usize, q: usize| vol * if p == q { 0.1 } else { 0.05 };
        let mut ed = [(0usize, 0.0f64); 6];
        for i in 0..6 {
            ed[i] = edge(t[TE[i].0], t[TE[i].1]);
        }
        for i in 0..6 {
            if cons[ed[i].0] {
                continue;
            }
            for j in i..6 {
                if cons[ed[j].0] {
                    continue;
                }
                let ((a, b), (c, d)) = (TE[i], TE[j]);
                let sij = 4.0 * vol * dot3(cross3(gg[a], gg[b]), cross3(gg[c], gg[d]));
                let mij = dot3(gg[b], gg[d]) * ii(a, c) - dot3(gg[b], gg[c]) * ii(a, d)
                    - dot3(gg[a], gg[d]) * ii(b, c)
                    + dot3(gg[a], gg[c]) * ii(b, d);
                let sg = ed[i].1 * ed[j].1;
                add(dof[ed[i].0], dof[ed[j].0], sij * sg / mur, epsc.rs(mij * sg), cx(0.0, 0.0));
            }
        }
    }
    // boundary faces: abc gets weight sqrt(eps/mur) of the touching tet, a
    // port sheet gets eta0 / Zs with Zs = z0 w / h; both multiply the
    // tangential face mass matrix. Port excitation integrates dir . whitney.
    for (t, g) in &mesh.tris {
        let (isabc, port) = (role[*g] == 2, pmap.get(g));
        if !isabc && port.is_none() {
            continue;
        }
        let (q0, q1, q2) = (mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]]);
        let (u, v) = (sub3(q1, q0), sub3(q2, q0));
        let n = cross3(u, v);
        let n2 = dot3(n, n);
        let at = n2.sqrt() / 2.0;
        let gs1 = sc3(cross3(v, n), 1.0 / n2);
        let gs2 = sc3(cross3(n, u), 1.0 / n2);
        let gs = [sub3(sub3([0.0; 3], gs1), gs2), gs1, gs2];
        let wt = match port {
            Some(&p) => cx(ETA0 * ports[p].h / (ports[p].z0 * ports[p].w), 0.0),
            None => {
                let ti = fmap[&{
                    let mut k = [t[0] as u32, t[1] as u32, t[2] as u32];
                    k.sort();
                    k
                }];
                let (epsc, mur) = matg[mesh.tets[ti].1];
                (epsc.rs(1.0 / mur)).sqrt()
            }
        };
        let i2 = |p: usize, q: usize| at * if p == q { 1.0 / 6.0 } else { 1.0 / 12.0 };
        let mut ed = [(0usize, 0.0f64); 3];
        for i in 0..3 {
            ed[i] = edge(t[SE[i].0], t[SE[i].1]);
        }
        for i in 0..3 {
            if cons[ed[i].0] {
                continue;
            }
            if let Some(&p) = port {
                let (a, b) = SE[i];
                let c = at / 3.0 * dot3(deck.ports[p].dir, sub3(gs[b], gs[a])) * ed[i].1;
                *ports[p].exc.entry(dof[ed[i].0]).or_insert(0.0) += c;
            }
            for j in i..3 {
                if cons[ed[j].0] {
                    continue;
                }
                let ((a, b), (c, d)) = (SE[i], SE[j]);
                let tij = dot3(gs[b], gs[d]) * i2(a, c) - dot3(gs[b], gs[c]) * i2(a, d)
                    - dot3(gs[a], gs[d]) * i2(b, c)
                    + dot3(gs[a], gs[c]) * i2(b, d);
                add(dof[ed[i].0], dof[ed[j].0], 0.0, cx(0.0, 0.0), wt.rs(tij * ed[i].1 * ed[j].1));
            }
        }
    }
    // fill reducing ordering on edge midpoints, then the permuted upper
    // triangle in CSC with a fixed slot per assembled entry
    let keys: Vec<(u32, u32)> = acc.keys().copied().collect();
    let mut pts = vec![[0.0; 3]; ndof];
    for e in 0..ne {
        if dof[e] != usize::MAX {
            let (pa, pb) = (mesh.nodes[enodes[e].0 as usize], mesh.nodes[enodes[e].1 as usize]);
            pts[dof[e]] = sc3([pa[0] + pb[0], pa[1] + pb[1], pa[2] + pb[2]], 0.5);
        }
    }
    let perm = ndorder(ndof, &keys, &pts);
    let mut ap = vec![0usize; ndof + 1];
    for &(a, b) in &keys {
        ap[perm[a as usize].max(perm[b as usize]) + 1] += 1;
    }
    for j in 0..ndof {
        ap[j + 1] += ap[j];
    }
    let mut nxt = ap.clone();
    let mut ai = vec![0u32; keys.len()];
    let ents: Vec<(usize, f64, Cx, Cx)> = acc
        .iter()
        .map(|(&(a, b), &(s, m, bb))| {
            let (i, j) = (perm[a as usize].min(perm[b as usize]), perm[a as usize].max(perm[b as usize]));
            let pos = nxt[j];
            nxt[j] += 1;
            ai[pos] = i as u32;
            (pos, s, m, bb)
        })
        .collect();
    let sym = symbolic(ndof, &ap, &ai);
    eprintln!(
        "nanofem: {} nodes, {} tets, {} edges, {} dofs, {:.2}M nnz A, {:.2}M nnz L",
        mesh.nodes.len(),
        mesh.tets.len(),
        ne,
        ndof,
        ents.len() as f64 / 1e6,
        sym.lp[ndof] as f64 / 1e6
    );
    // frequency sweep: factor once per frequency, one solve per port, port
    // voltage V = -(h/A) integral E . dir, S(q,p) = 2 Vq / sqrt(z0p z0q) - d(q,p)
    let np = ports.len();
    outln(format_args!("# Hz S RI R {}", ports[0].z0));
    let mut ax = vec![cx(0.0, 0.0); ents.len()];
    let mut w = work(ndof, &sym);
    for fi in 0..deck.nf {
        let f = if deck.nf == 1 { deck.f0 } else { deck.f0 + (deck.f1 - deck.f0) * fi as f64 / (deck.nf - 1) as f64 };
        let k0 = 2.0 * std::f64::consts::PI * f / C0;
        for &(pos, s, m, b) in &ents {
            ax[pos] = cx(s, 0.0) + m.rs(-k0 * k0) + b.js(k0);
        }
        numeric(ndof, &ap, &ai, &ax, &sym, &mut w);
        let mut smat = vec![vec![cx(0.0, 0.0); np]; np];
        for p in 0..np {
            let mut rhs = vec![cx(0.0, 0.0); ndof];
            for (&d, &c) in &ports[p].exc {
                rhs[perm[d]] = cx(0.0, -k0 * ETA0 * c / ports[p].w);
            }
            ldsolve(ndof, &sym, &w, &mut rhs);
            for q in 0..np {
                let mut vq = cx(0.0, 0.0);
                for (&d, &c) in &ports[q].exc {
                    vq = vq + rhs[perm[d]].rs(c);
                }
                vq = vq.rs(-ports[q].h / ports[q].area);
                smat[q][p] = vq.rs(2.0 / (ports[p].z0 * ports[q].z0).sqrt());
                if q == p {
                    smat[q][p] = smat[q][p] - cx(1.0, 0.0);
                }
            }
        }
        // touchstone: 2-port data is S11 S21 S12 S22 on one line, larger
        // matrices row by row with the frequency on the first line
        let pair = |s: Cx| format!(" {:.9e} {:.9e}", s.re, s.im);
        if np <= 2 {
            let mut l = format!("{:.9e}", f);
            for p in 0..np {
                for q in 0..np {
                    l += &pair(smat[q][p]);
                }
            }
            outln(format_args!("{}", l));
        } else {
            for q in 0..np {
                let mut l = if q == 0 { format!("{:.9e}", f) } else { String::new() };
                for p in 0..np {
                    l += &pair(smat[q][p]);
                }
                outln(format_args!("{}", l.trim_start()));
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        die("usage: nanofem <deck.nfm>");
    }
    let deck = parse_deck(&args[1]);
    let mesh = parse_msh(&deck.mesh);
    simulate(&deck, &mesh);
}
