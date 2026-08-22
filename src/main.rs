// nanofem: a headless 3D finite element electromagnetic field solver in a
// single source file, hard budget 1000 lines of code, std only.
//
// Algorithm: time harmonic (exp(+j w t)) curl-curl equation for the electric
// field, discretized with first order Nedelec (Whitney) edge elements on
// tetrahedra. Boundary conditions: PEC as eliminated unknowns, first order
// absorbing boundary, PML regions as a complex coordinate stretch, natural
// PMC. Lumped rectangular ports are modeled as impedance sheets with an
// impressed surface current (Norton equivalent); scattering parameters
// follow from the port voltages. Materials carry relative permittivity,
// loss tangent and permeability per mesh region.
// Input: a Gmsh .msh version 2.2 ASCII mesh with physical groups plus a deck
// that maps group names to materials, boundaries and ports. Output:
// Touchstone data on stdout, optionally the E field as legacy VTK. Linear
// algebra: geometric nested dissection ordering and a complex symmetric
// sparse LDL^T factorization, one factorization per frequency in parallel.

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
// a^T diag(w) b for a diagonal complex material tensor
fn cdot3(a: V3, b: V3, w: &[Cx; 3]) -> Cx {
    w[0].rs(a[0] * b[0]) + w[1].rs(a[1] * b[1]) + w[2].rs(a[2] * b[2])
}

// Gradients of the barycentric coordinates of a tetrahedron and its volume.
// The gradients are constant over the element and carry all its geometry.
fn tetgrad(p: [V3; 4]) -> ([V3; 4], f64) {
    let (r1, r2, r3) = (sub3(p[1], p[0]), sub3(p[2], p[0]), sub3(p[3], p[0]));
    let det = dot3(r1, cross3(r2, r3));
    if det == 0.0 {
        die("degenerate tetrahedron");
    }
    let g = [sc3(cross3(r2, r3), 1.0 / det), sc3(cross3(r3, r1), 1.0 / det), sc3(cross3(r1, r2), 1.0 / det)];
    ([sub3(sub3(sub3([0.0; 3], g[0]), g[1]), g[2]), g[0], g[1], g[2]], det.abs() / 6.0)
}

// Mass matrix entry of two Whitney functions W_ab = l_a grad l_b - l_b grad
// l_a, given the gradients and the integrals I(p,q) of l_p l_q over the
// element. The same expression serves the tetrahedron and the face.
fn wmass(g: &[V3], (a, b): (usize, usize), (c, d): (usize, usize), dot: &dyn Fn(V3, V3) -> Cx, i: &dyn Fn(usize, usize) -> f64) -> Cx {
    dot(g[b], g[d]).rs(i(a, c)) - dot(g[b], g[c]).rs(i(a, d)) - dot(g[a], g[d]).rs(i(b, c)) + dot(g[a], g[c]).rs(i(b, d))
}

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
    fn nx<'a>(it: &mut std::str::Lines<'a>) -> &'a str {
        it.next().unwrap_or_else(|| die("mesh file ends inside a section"))
    }
    fn cnt(it: &mut std::str::Lines) -> usize {
        nx(it).trim().parse().unwrap_or_else(|_| die("bad section count"))
    }
    let mut it = txt.lines();
    while let Some(head) = it.next() {
        match head.trim() {
            "$MeshFormat" => {
                if !nx(&mut it).trim_start().starts_with("2.2") {
                    die("mesh must be Gmsh ASCII v2.2 (export with -format msh22)");
                }
            }
            "$PhysicalNames" => {
                for _ in 0..cnt(&mut it) {
                    let t: Vec<&str> = nx(&mut it).split_whitespace().collect();
                    let (dim, id): (i64, i64) = (t[0].parse().unwrap_or(-1), t[1].parse().unwrap_or(-1));
                    phys.insert((dim, id), t[2..].join(" ").trim_matches('"').to_string());
                }
            }
            "$Nodes" => {
                for _ in 0..cnt(&mut it) {
                    let t: Vec<f64> = nx(&mut it).split_whitespace().map(|s| s.parse().unwrap_or(f64::NAN)).collect();
                    idmap.insert(t[0] as i64, m.nodes.len());
                    m.nodes.push([t[1], t[2], t[3]]);
                }
            }
            "$Elements" => {
                for _ in 0..cnt(&mut it) {
                    let t: Vec<i64> = nx(&mut it).split_whitespace().map(|s| s.parse().unwrap_or(-1)).collect();
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
                    let mut nd: Vec<usize> = t[3 + ntag..]
                        .iter()
                        .map(|id| *idmap.get(id).unwrap_or_else(|| die("element references unknown node")))
                        .collect();
                    // vertices ascending: every local edge then runs from
                    // the lower to the higher global node, so the Whitney
                    // functions of neighboring elements agree by
                    // construction and no orientation signs are needed
                    nd.sort();
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
//   pml <group> <ax> <ay> <az>           (imaginary coordinate stretch per axis)
//   port <n> <group> <jx> <jy> <jz> <z0> (n counting from 1, j = voltage direction)
//   sweep lin <f0> <f1> <npoints>
//   field <path.vtk> <f>                 (E field snapshot, port 1 driven)
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
    pmls: Vec<(String, V3)>,
    pec: Vec<String>,
    abc: Vec<String>,
    ports: Vec<PortDef>,
    f0: f64,
    f1: f64,
    nf: usize,
    field: Option<(String, f64)>,
}

fn num(s: &str) -> f64 {
    s.parse().unwrap_or_else(|_| die(&format!("bad number '{}'", s)))
}

fn parse_deck(path: &str) -> Deck {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("cannot read {}: {}", path, e)));
    let mut d = Deck { mesh: String::new(), mats: vec![], pmls: vec![], pec: vec![], abc: vec![], ports: vec![], f0: 0.0, f1: 0.0, nf: 0, field: None };
    let rel = |p: &str| -> String {
        let base = std::path::Path::new(path).parent().unwrap_or(std::path::Path::new(""));
        base.join(p).to_string_lossy().into_owned()
    };
    for l in txt.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.is_empty() || t[0].starts_with('*') {
            continue;
        }
        let bad = || -> ! { die(&format!("bad card: {}", l)) };
        match t[0].to_lowercase().as_str() {
            "mesh" if t.len() == 2 => d.mesh = rel(t[1]),
            "field" if t.len() == 3 => d.field = Some((rel(t[1]), num(t[2]))),
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
            "pml" if t.len() == 5 => d.pmls.push((t[1].to_string(), [num(t[2]), num(t[3]), num(t[4])])),
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

struct Fac {
    li: Vec<u32>,
    lx: Vec<Cx>,
    d: Vec<Cx>,
}

fn fac(n: usize, sym: &Sym) -> Fac {
    Fac { li: vec![0; sym.lp[n]], lx: vec![cx(0.0, 0.0); sym.lp[n]], d: vec![cx(0.0, 0.0); n] }
}

fn numeric(n: usize, ap: &[usize], ai: &[u32], ax: &[Cx], sym: &Sym, w: &mut Fac) {
    let (mut y, mut flag) = (vec![cx(0.0, 0.0); n], vec![usize::MAX; n]);
    let (mut pat, mut lnz) = (vec![0usize; n], vec![0usize; n]);
    for k in 0..n {
        let mut top = n;
        flag[k] = k;
        lnz[k] = 0;
        for p in ap[k]..ap[k + 1] {
            let mut i = ai[p] as usize;
            y[i] = y[i] + ax[p];
            let mut len = 0;
            while i < k && flag[i] != k {
                pat[len] = i;
                len += 1;
                flag[i] = k;
                i = sym.parent[i];
            }
            while len > 0 {
                len -= 1;
                top -= 1;
                pat[top] = pat[len];
            }
        }
        w.d[k] = y[k];
        y[k] = cx(0.0, 0.0);
        for t in top..n {
            let i = pat[t];
            let yi = y[i];
            y[i] = cx(0.0, 0.0);
            for p in sym.lp[i]..sym.lp[i] + lnz[i] {
                y[w.li[p] as usize] = y[w.li[p] as usize] - w.lx[p] * yi;
            }
            let lki = yi / w.d[i];
            w.d[k] = w.d[k] - lki * yi;
            let p = sym.lp[i] + lnz[i];
            w.li[p] = k as u32;
            w.lx[p] = lki;
            lnz[i] += 1;
        }
        if w.d[k].mag() < 1e-300 {
            die("singular system matrix");
        }
    }
}

fn ldsolve(n: usize, sym: &Sym, w: &Fac, b: &mut [Cx]) {
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
    // sheet impedance relative to eta0, the weight of the port face matrix
    zs: f64,
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
    // materials per volume group. A pml region stretches coordinate k by
    // s_k = 1 - j a_k; eps and mu both get the diagonal tensor
    // L_k = s_x s_y s_z / s_k^2, which keeps the wave impedance matched
    // while fields decay along the stretched axes.
    struct Mg {
        wm: [Cx; 3], // mass weight, the eps tensor
        ws: [Cx; 3], // stiffness weight, the 1/mu tensor
        abc: Cx,     // sqrt(eps/mu) of the region, the absorbing weight
    }
    let one = cx(1.0, 0.0);
    let (mut eps, mut mur, mut stretch) = (vec![one; ng], vec![1.0; ng], vec![[0.0; 3]; ng]);
    for (name, mt) in &deck.mats {
        let g = gidx(name);
        eps[g] = cx(mt.eps, -mt.eps * mt.tand);
        mur[g] = mt.mur;
    }
    for (name, a) in &deck.pmls {
        stretch[gidx(name)] = *a;
    }
    let matg: Vec<Mg> = (0..ng)
        .map(|g| {
            let s = [cx(1.0, -stretch[g][0]), cx(1.0, -stretch[g][1]), cx(1.0, -stretch[g][2])];
            let prod = s[0] * s[1] * s[2];
            let lam = |k: usize| prod / (s[k] * s[k]);
            Mg {
                wm: [eps[g] * lam(0), eps[g] * lam(1), eps[g] * lam(2)],
                ws: [(one / lam(0)).rs(1.0 / mur[g]), (one / lam(1)).rs(1.0 / mur[g]), (one / lam(2)).rs(1.0 / mur[g])],
                abc: (eps[g].rs(1.0 / mur[g])).sqrt(),
            }
        })
        .collect();
    // global edges; element vertices are sorted, so a local edge (a, b)
    // always has a < b and needs no orientation sign
    let mut emap: HashMap<(u32, u32), u32> = HashMap::new();
    for (t, _) in &mesh.tets {
        for (a, b) in TE {
            let n = emap.len() as u32;
            emap.entry((t[a] as u32, t[b] as u32)).or_insert(n);
        }
    }
    let ne = emap.len();
    let edge = |a: usize, b: usize| -> usize {
        match emap.get(&(a as u32, b as u32)) {
            Some(&e) => e as usize,
            None => die("boundary triangle does not match the volume mesh"),
        }
    };
    // PEC: eliminate all edges of pec triangles. dof carries the numbering
    // of the remaining unknowns, usize::MAX marks a constrained edge.
    let mut dof = vec![0usize; ne];
    for (t, g) in &mesh.tris {
        if role[*g] == 1 {
            for (a, b) in SE {
                dof[edge(t[a], t[b])] = usize::MAX;
            }
        }
    }
    let mut ndof = 0;
    for e in 0..ne {
        if dof[e] != usize::MAX {
            dof[e] = ndof;
            ndof += 1;
        }
    }
    let free = |e: usize| dof[e] != usize::MAX;
    // boundary face -> adjacent tet, for the abc material weight
    let mut fmap: HashMap<[u32; 3], usize> = HashMap::new();
    for (ti, (t, _)) in mesh.tets.iter().enumerate() {
        for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut k = [t[f[0]] as u32, t[f[1]] as u32, t[f[2]] as u32];
            k.sort();
            fmap.insert(k, ti);
        }
    }
    // Port geometry: the sheet is a rectangle of height h along dir and
    // width w, so its area gives w once h is known. Only w survives: the
    // sheet impedance is z0 w / h and the port voltage is the mean field
    // along dir times h, that is the surface integral divided by w.
    let mut ports: Vec<PortDat> =
        deck.ports.iter().map(|pd| PortDat { z0: pd.z0, w: 0.0, zs: 0.0, exc: HashMap::new() }).collect();
    let mut geo = vec![(0.0, f64::MAX, f64::MIN); ports.len()];
    for (t, g) in &mesh.tris {
        if let Some(&p) = pmap.get(g) {
            let n = cross3(sub3(mesh.nodes[t[1]], mesh.nodes[t[0]]), sub3(mesh.nodes[t[2]], mesh.nodes[t[0]]));
            geo[p].0 += 0.5 * dot3(n, n).sqrt();
            for v in t {
                let x = dot3(mesh.nodes[*v], deck.ports[p].dir);
                geo[p].1 = geo[p].1.min(x);
                geo[p].2 = geo[p].2.max(x);
            }
        }
    }
    for (p, pt) in ports.iter_mut().enumerate() {
        let (area, h) = (geo[p].0, geo[p].2 - geo[p].1);
        if area <= 0.0 || h <= 0.0 {
            die(&format!("port {} has no faces or zero extent along its direction", p + 1));
        }
        pt.w = area / h;
        pt.zs = ETA0 * h / (pt.z0 * pt.w);
    }
    // Assembly: each matrix entry keeps three frequency independent
    // coefficients, one per power of k0. The curl stiffness sits at power
    // 0, the eps weighted mass enters as -k0^2 times its value at power 2,
    // and the boundary term as +j k0 times its value at power 1.
    let mut acc: HashMap<(u32, u32), [Cx; 3]> = HashMap::new();
    let mut add = |di: usize, dj: usize, pow: usize, v: Cx| {
        let e = acc.entry((di.min(dj) as u32, di.max(dj) as u32)).or_insert([cx(0.0, 0.0); 3]);
        e[pow] = e[pow] + v;
    };
    for (t, g) in &mesh.tets {
        let (gg, vol) = tetgrad([mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]], mesh.nodes[t[3]]]);
        let mg = &matg[*g];
        let ii = |p: usize, q: usize| vol * if p == q { 0.1 } else { 0.05 };
        let mut ed = [0usize; 6];
        for i in 0..6 {
            ed[i] = edge(t[TE[i].0], t[TE[i].1]);
        }
        for i in 0..6 {
            if !free(ed[i]) {
                continue;
            }
            for j in i..6 {
                if !free(ed[j]) {
                    continue;
                }
                let ((a, b), (c, d)) = (TE[i], TE[j]);
                let sij = cdot3(cross3(gg[a], gg[b]), cross3(gg[c], gg[d]), &mg.ws).rs(4.0 * vol);
                let mij = wmass(&gg, TE[i], TE[j], &|x, y| cdot3(x, y, &mg.wm), &ii);
                add(dof[ed[i]], dof[ed[j]], 0, sij);
                add(dof[ed[i]], dof[ed[j]], 2, mij.rs(-1.0));
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
            Some(&p) => cx(ports[p].zs, 0.0),
            None => {
                let ti = fmap[&{
                    let mut k = [t[0] as u32, t[1] as u32, t[2] as u32];
                    k.sort();
                    k
                }];
                matg[mesh.tets[ti].1].abc
            }
        };
        let i2 = |p: usize, q: usize| at * if p == q { 1.0 / 6.0 } else { 1.0 / 12.0 };
        let mut ed = [0usize; 3];
        for i in 0..3 {
            ed[i] = edge(t[SE[i].0], t[SE[i].1]);
        }
        for i in 0..3 {
            if !free(ed[i]) {
                continue;
            }
            if let Some(&p) = port {
                let (a, b) = SE[i];
                let c = at / 3.0 * dot3(deck.ports[p].dir, sub3(gs[b], gs[a]));
                *ports[p].exc.entry(dof[ed[i]]).or_insert(0.0) += c;
            }
            for j in i..3 {
                if !free(ed[j]) {
                    continue;
                }
                let tij = wmass(&gs, SE[i], SE[j], &|x, y| cx(dot3(x, y), 0.0), &i2);
                add(dof[ed[i]], dof[ed[j]], 1, (wt * tij).js(1.0));
            }
        }
    }
    // fill reducing ordering on edge midpoints, then the permuted upper
    // triangle in CSC with a fixed slot per assembled entry
    let keys: Vec<(u32, u32)> = acc.keys().copied().collect();
    let mut pts = vec![[0.0; 3]; ndof];
    for (&(a, b), &e) in &emap {
        if dof[e as usize] != usize::MAX {
            let (pa, pb) = (mesh.nodes[a as usize], mesh.nodes[b as usize]);
            pts[dof[e as usize]] = sc3([pa[0] + pb[0], pa[1] + pb[1], pa[2] + pb[2]], 0.5);
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
    let ents: Vec<(usize, [Cx; 3])> = acc
        .iter()
        .map(|(&(a, b), &c)| {
            let (i, j) = (perm[a as usize].min(perm[b as usize]), perm[a as usize].max(perm[b as usize]));
            let pos = nxt[j];
            nxt[j] += 1;
            ai[pos] = i as u32;
            (pos, c)
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
    // voltage V = -(h/A) integral E . dir, S(q,p) = 2 Vq / sqrt(z0p z0q) - d(q,p).
    // Frequencies run in parallel, thread count capped by memory (4 GB for
    // the factorizations), lines are printed in order afterwards.
    let np = ports.len();
    outln(format_args!("# Hz S RI R {}", ports[0].z0));
    let freqs: Vec<f64> = (0..deck.nf)
        .map(|fi| if deck.nf == 1 { deck.f0 } else { deck.f0 + (deck.f1 - deck.f0) * fi as f64 / (deck.nf - 1) as f64 })
        .collect();
    let bytes = 24 * sym.lp[ndof] + 16 * ents.len() + 64 * ndof;
    let nt = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1)
        .min(deck.nf)
        .min((4_000_000_000 / bytes).max(1));
    let lines = std::sync::Mutex::new(vec![String::new(); deck.nf]);
    std::thread::scope(|sc| {
        for t in 0..nt {
            let (ports, ents, ap, ai, sym, perm, freqs, lines) = (&ports, &ents, &ap, &ai, &sym, &perm, &freqs, &lines);
            sc.spawn(move || {
                let mut ax = vec![cx(0.0, 0.0); ents.len()];
                let mut w = fac(ndof, sym);
                for fi in (t..freqs.len()).step_by(nt) {
                    let k0 = 2.0 * std::f64::consts::PI * freqs[fi] / C0;
                    for &(pos, c) in ents.iter() {
                        ax[pos] = c[0] + (c[1] + c[2].rs(k0)).rs(k0);
                    }
                    numeric(ndof, ap, ai, &ax, sym, &mut w);
                    let mut smat = vec![vec![cx(0.0, 0.0); np]; np];
                    for p in 0..np {
                        let mut rhs = vec![cx(0.0, 0.0); ndof];
                        for (&d, &c) in &ports[p].exc {
                            rhs[perm[d]] = cx(0.0, -k0 * ETA0 * c / ports[p].w);
                        }
                        ldsolve(ndof, sym, &w, &mut rhs);
                        for q in 0..np {
                            let mut vq = cx(0.0, 0.0);
                            for (&d, &c) in &ports[q].exc {
                                vq = vq + rhs[perm[d]].rs(c);
                            }
                            vq = vq.rs(-1.0 / ports[q].w);
                            smat[q][p] = vq.rs(2.0 / (ports[p].z0 * ports[q].z0).sqrt());
                            if q == p {
                                smat[q][p] = smat[q][p] - cx(1.0, 0.0);
                            }
                        }
                    }
                    // touchstone: 2-port data is S11 S21 S12 S22 on one
                    // line, larger matrices row by row
                    let pair = |s: Cx| format!(" {:.9e} {:.9e}", s.re, s.im);
                    let mut out = format!("{:.9e}", freqs[fi]);
                    if np <= 2 {
                        for p in 0..np {
                            for q in 0..np {
                                out += &pair(smat[q][p]);
                            }
                        }
                    } else {
                        for q in 0..np {
                            for p in 0..np {
                                out += &pair(smat[q][p]);
                            }
                            if q + 1 < np {
                                out += "\n";
                            }
                        }
                    }
                    lines.lock().unwrap()[fi] = out;
                }
            });
        }
    });
    for l in lines.into_inner().unwrap() {
        outln(format_args!("{}", l));
    }
    // field snapshot: solve once more at the requested frequency with port 1
    // driven and write the E field per tet (centroid value, where the
    // whitney functions reduce to (g_b - g_a) / 4) as legacy VTK cell data
    if let Some((path, ff)) = &deck.field {
        let k0 = 2.0 * std::f64::consts::PI * ff / C0;
        let mut ax = vec![cx(0.0, 0.0); ents.len()];
        for &(pos, c) in &ents {
            ax[pos] = c[0] + (c[1] + c[2].rs(k0)).rs(k0);
        }
        let mut w = fac(ndof, &sym);
        numeric(ndof, &ap, &ai, &ax, &sym, &mut w);
        let mut rhs = vec![cx(0.0, 0.0); ndof];
        for (&d, &c) in &ports[0].exc {
            rhs[perm[d]] = cx(0.0, -k0 * ETA0 * c / ports[0].w);
        }
        ldsolve(ndof, &sym, &w, &mut rhs);
        let mut ef = Vec::with_capacity(mesh.tets.len());
        for (t, _) in &mesh.tets {
            let (gg, _) = tetgrad([mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]], mesh.nodes[t[3]]]);
            let mut e = [cx(0.0, 0.0); 3];
            for (a, b) in TE {
                let ei = edge(t[a], t[b]);
                if free(ei) {
                    let x = rhs[perm[dof[ei]]].rs(0.25);
                    let d = sub3(gg[b], gg[a]);
                    for k in 0..3 {
                        e[k] = e[k] + x.rs(d[k]);
                    }
                }
            }
            ef.push(e);
        }
        let mut o = format!(
            "# vtk DataFile Version 3.0\nnanofem E field at {} Hz, port 1 driven\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS {} double\n",
            ff,
            mesh.nodes.len()
        );
        for p in &mesh.nodes {
            o += &format!("{} {} {}\n", p[0], p[1], p[2]);
        }
        o += &format!("CELLS {} {}\n", mesh.tets.len(), 5 * mesh.tets.len());
        for (t, _) in &mesh.tets {
            o += &format!("4 {} {} {} {}\n", t[0], t[1], t[2], t[3]);
        }
        o += &format!("CELL_TYPES {}\n", mesh.tets.len());
        o += &"10\n".repeat(mesh.tets.len());
        o += &format!("CELL_DATA {}\n", mesh.tets.len());
        for (head, comp) in [
            ("VECTORS Ere double", 0),
            ("VECTORS Eim double", 1),
            ("SCALARS Emag double\nLOOKUP_TABLE default", 2),
        ] {
            o += &format!("{}\n", head);
            for e in &ef {
                match comp {
                    0 => o += &format!("{:.6e} {:.6e} {:.6e}\n", e[0].re, e[1].re, e[2].re),
                    1 => o += &format!("{:.6e} {:.6e} {:.6e}\n", e[0].im, e[1].im, e[2].im),
                    _ => o += &format!("{:.6e}\n", (e[0].mag().powi(2) + e[1].mag().powi(2) + e[2].mag().powi(2)).sqrt()),
                }
            }
        }
        std::fs::write(path, o).unwrap_or_else(|e| die(&format!("cannot write {}: {}", path, e)));
        eprintln!("nanofem: field written to {}", path);
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
