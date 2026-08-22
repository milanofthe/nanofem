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

fn add3(a: V3, b: V3) -> V3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
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

// -------------------------------------------------------------- quadrature

// Fully symmetric rules on the reference simplex, in barycentric
// coordinates, with weights summing to one. The tetrahedron rule is exact
// to degree 5 and the triangle rule to degree 4, which covers the mass
// matrices of both element orders.
fn qtet() -> Vec<(f64, [f64; 4])> {
    let mut q = vec![];
    for (w, a) in [(0.0734930431163619, 0.0927352503108912), (0.1126879257180162, 0.3108859192633005)] {
        for k in 0..4 {
            let mut p = [a; 4];
            p[k] = 1.0 - 3.0 * a;
            q.push((w, p));
        }
    }
    let (w, a) = (0.0425460207770812, 0.0455037041256497);
    for i in 0..3 {
        for j in i + 1..4 {
            let mut p = [0.5 - a; 4];
            p[i] = a;
            p[j] = a;
            q.push((w, p));
        }
    }
    q
}

fn qtri() -> Vec<(f64, [f64; 3])> {
    let mut q = vec![];
    for (w, a) in [(0.223381589678011, 0.445948490915965), (0.109951743655322, 0.091576213509771)] {
        for k in 0..3 {
            let mut p = [a; 3];
            p[k] = 1.0 - 2.0 * a;
            q.push((w, p));
        }
    }
    q
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
                    // vertices ascending: local edges and faces then follow
                    // the global numbering, which fixes the orientation of
                    // every basis function without a sign convention
                    if dim == 2 {
                        let mut k = [nd[0], nd[1], nd[2]];
                        k.sort();
                        m.tris.push((k, gi));
                    } else {
                        let mut k = [nd[0], nd[1], nd[2], nd[3]];
                        k.sort();
                        m.tets.push((k, gi));
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
//   order <1|2>                          (element order, default 1)
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
    order: usize,
}

fn num(s: &str) -> f64 {
    s.parse().unwrap_or_else(|_| die(&format!("bad number '{}'", s)))
}

fn parse_deck(path: &str) -> Deck {
    let txt = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("cannot read {}: {}", path, e)));
    let mut d = Deck { mesh: String::new(), mats: vec![], pmls: vec![], pec: vec![], abc: vec![], ports: vec![], f0: 0.0, f1: 0.0, nf: 0, field: None, order: 1 };
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
            "order" if t.len() == 2 => {
                d.order = num(t[1]) as usize;
                if d.order < 1 || d.order > 2 {
                    die("order must be 1 or 2");
                }
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

// ------------------------------------------------------------------ basis

const TE: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
const SE: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];
const TF: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
// largest local basis size, reached by the second order tetrahedron
const NB: usize = 20;

// Value and curl of the local H(curl) basis of a tetrahedron at the point
// with barycentric coordinates l, given the gradients g of those
// coordinates. Order 1 is the six Whitney edge functions
// W_ab = l_a grad l_b - l_b grad l_a, whose curl 2 grad l_a x grad l_b is
// constant over the element. Order 2 appends the six curl free edge
// gradients grad(l_a l_b) and then, per face abc, the pair l_a W_bc and
// l_b W_ca; the third such product is minus their sum, so two are
// independent. Element vertices arrive in ascending global order, so every
// function is pinned by the global numbering and no signs are needed.
fn tbasis(order: usize, g: &[V3; 4], l: [f64; 4], out: &mut Vec<(V3, V3)>) {
    let w = |a: usize, b: usize| sub3(sc3(g[b], l[a]), sc3(g[a], l[b]));
    let cw = |a: usize, b: usize| sc3(cross3(g[a], g[b]), 2.0);
    out.clear();
    for (a, b) in TE {
        out.push((w(a, b), cw(a, b)));
    }
    if order < 2 {
        return;
    }
    for (a, b) in TE {
        out.push((add3(sc3(g[b], l[a]), sc3(g[a], l[b])), [0.0; 3]));
    }
    for f in TF {
        for (x, y, z) in [(f[0], f[1], f[2]), (f[1], f[2], f[0])] {
            // curl(l_x W_yz) = grad l_x x W_yz + l_x curl W_yz
            out.push((sc3(w(y, z), l[x]), add3(cross3(g[x], w(y, z)), sc3(cw(y, z), l[x]))));
        }
    }
}

// The same basis restricted to a boundary triangle, where only the value
// matters. g holds the gradients of the surface barycentric coordinates,
// whose tangential parts equal those of the volume gradients.
fn sbasis(order: usize, g: &[V3; 3], l: [f64; 3], out: &mut Vec<V3>) {
    let w = |a: usize, b: usize| sub3(sc3(g[b], l[a]), sc3(g[a], l[b]));
    out.clear();
    for (a, b) in SE {
        out.push(w(a, b));
    }
    if order < 2 {
        return;
    }
    for (a, b) in SE {
        out.push(add3(sc3(g[b], l[a]), sc3(g[a], l[b])));
    }
    for (x, y, z) in [(0, 1, 2), (1, 2, 0)] {
        out.push(sc3(w(y, z), l[x]));
    }
}

// ---------------------------------------------------------------- fem

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
    // materials per volume group. A pml region stretches coordinate k by
    // s_k = 1 - j a_k; eps and mu both get the diagonal tensor
    // L_k = s_x s_y s_z / s_k^2, which keeps the wave impedance matched
    // while fields decay along the stretched axes.
    struct Mg {
        epsc: Cx,
        mur: f64,
        wm: [Cx; 3], // mass weight, the eps tensor
        ws: [Cx; 3], // stiffness weight, the 1/mu tensor
    }
    let one = cx(1.0, 0.0);
    let mut matg: Vec<Mg> = (0..ng).map(|_| Mg { epsc: one, mur: 1.0, wm: [one; 3], ws: [one; 3] }).collect();
    for (name, mt) in &deck.mats {
        let g = &mut matg[gidx(name)];
        g.epsc = cx(mt.eps, -mt.eps * mt.tand);
        g.mur = mt.mur;
    }
    let mut stretch = vec![[0.0; 3]; ng];
    for (name, a) in &deck.pmls {
        stretch[gidx(name)] = *a;
    }
    for g in 0..ng {
        let s = [cx(1.0, -stretch[g][0]), cx(1.0, -stretch[g][1]), cx(1.0, -stretch[g][2])];
        let prod = s[0] * s[1] * s[2];
        for k in 0..3 {
            let lam = prod / (s[k] * s[k]);
            matg[g].wm[k] = matg[g].epsc * lam;
            matg[g].ws[k] = (one / lam).rs(1.0 / matg[g].mur);
        }
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
    let edge = |a: usize, b: usize| -> usize {
        match emap.get(&(a.min(b) as u32, a.max(b) as u32)) {
            Some(&e) => e as usize,
            None => die("boundary triangle does not match the volume mesh"),
        }
    };
    // global faces, numbered for the order 2 face unknowns and carrying one
    // adjacent tetrahedron for the abc material weight
    let mut fmap: HashMap<[u32; 3], (usize, usize)> = HashMap::new();
    let mut fnodes: Vec<[u32; 3]> = vec![];
    for (ti, (t, _)) in mesh.tets.iter().enumerate() {
        for f in TF {
            let k = [t[f[0]] as u32, t[f[1]] as u32, t[f[2]] as u32];
            let n = fnodes.len();
            let e = fmap.entry(k).or_insert_with(|| {
                fnodes.push(k);
                (n, ti)
            });
            e.1 = ti;
        }
    }
    let face = |t: &[usize]| -> usize {
        match fmap.get(&[t[0] as u32, t[1] as u32, t[2] as u32]) {
            Some(&(f, _)) => f,
            None => die("boundary triangle does not match the volume mesh"),
        }
    };
    // Unknowns sit in slots: p per edge, then nfd per face. usize::MAX
    // marks a slot eliminated by a pec boundary, the rest are numbered.
    let (p, nf) = (deck.order, fnodes.len());
    let nfd = 2 * (p - 1);
    let eslot = |e: usize, k: usize| e * p + k;
    let fslot = |f: usize, k: usize| ne * p + f * nfd + k;
    let mut dof = vec![0usize; ne * p + nf * nfd];
    for (t, g) in &mesh.tris {
        if role[*g] == 1 {
            for (a, b) in SE {
                for k in 0..p {
                    dof[eslot(edge(t[a], t[b]), k)] = usize::MAX;
                }
            }
            for k in 0..nfd {
                dof[fslot(face(t), k)] = usize::MAX;
            }
        }
    }
    let mut ndof = 0;
    for s in dof.iter_mut() {
        if *s != usize::MAX {
            *s = ndof;
            ndof += 1;
        }
    }
    // Local unknowns of an element, in the order tbasis and sbasis produce
    // their functions: edge whitney, edge gradients, then the face pairs.
    let tdofs = |t: &[usize; 4]| -> [usize; NB] {
        let mut ld = [usize::MAX; NB];
        for i in 0..6 {
            let e = edge(t[TE[i].0], t[TE[i].1]);
            for k in 0..p {
                ld[i + 6 * k] = dof[eslot(e, k)];
            }
        }
        for (j, f) in TF.iter().enumerate() {
            for k in 0..nfd {
                ld[12 + 2 * j + k] = dof[fslot(face(&[t[f[0]], t[f[1]], t[f[2]]]), k)];
            }
        }
        ld
    };
    let sdofs = |t: &[usize; 3]| -> [usize; NB] {
        let mut ld = [usize::MAX; NB];
        for i in 0..3 {
            let e = edge(t[SE[i].0], t[SE[i].1]);
            for k in 0..p {
                ld[i + 3 * k] = dof[eslot(e, k)];
            }
        }
        for k in 0..nfd {
            ld[6 + k] = dof[fslot(face(t), k)];
        }
        ld
    };
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
    let mut acc: HashMap<(u32, u32), (Cx, Cx, Cx)> = HashMap::new();
    let mut add = |di: usize, dj: usize, s: Cx, m: Cx, b: Cx| {
        let k = (di.min(dj) as u32, di.max(dj) as u32);
        let e = acc.entry(k).or_insert((cx(0.0, 0.0), cx(0.0, 0.0), cx(0.0, 0.0)));
        e.0 = e.0 + s;
        e.1 = e.1 + m;
        e.2 = e.2 + b;
    };
    let (qt, qs) = (qtet(), qtri());
    let zero = cx(0.0, 0.0);
    let (mut bf, mut sf) = (vec![], vec![]);
    let (mut se, mut me) = ([[zero; NB]; NB], [[zero; NB]; NB]);
    for (t, g) in &mesh.tets {
        let (gg, vol) = tetgrad([mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]], mesh.nodes[t[3]]]);
        let mg = &matg[*g];
        let ld = tdofs(t);
        let nb = 6 * p + 4 * nfd;
        for r in se.iter_mut().chain(me.iter_mut()) {
            r.fill(zero);
        }
        for &(w, l) in &qt {
            tbasis(p, &gg, l, &mut bf);
            for i in 0..nb {
                for j in i..nb {
                    se[i][j] = se[i][j] + cdot3(bf[i].1, bf[j].1, &mg.ws).rs(w * vol);
                    me[i][j] = me[i][j] + cdot3(bf[i].0, bf[j].0, &mg.wm).rs(w * vol);
                }
            }
        }
        for i in 0..nb {
            for j in i..nb {
                if ld[i] != usize::MAX && ld[j] != usize::MAX {
                    add(ld[i], ld[j], se[i][j], me[i][j], zero);
                }
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
                let ti = fmap[&[t[0] as u32, t[1] as u32, t[2] as u32]].1;
                let mgf = &matg[mesh.tets[ti].1];
                (mgf.epsc.rs(1.0 / mgf.mur)).sqrt()
            }
        };
        let ld = sdofs(t);
        let nb = 3 * p + nfd;
        for r in me.iter_mut() {
            r.fill(zero);
        }
        let mut ex = [0.0f64; NB];
        for &(w, l) in &qs {
            sbasis(p, &gs, l, &mut sf);
            for i in 0..nb {
                if let Some(&pt) = port {
                    ex[i] += w * at * dot3(deck.ports[pt].dir, sf[i]);
                }
                for j in i..nb {
                    me[i][j] = me[i][j] + cx(dot3(sf[i], sf[j]) * w * at, 0.0);
                }
            }
        }
        for i in 0..nb {
            if ld[i] == usize::MAX {
                continue;
            }
            if let Some(&pt) = port {
                *ports[pt].exc.entry(ld[i]).or_insert(0.0) += ex[i];
            }
            for j in i..nb {
                if ld[j] != usize::MAX {
                    add(ld[i], ld[j], zero, zero, wt * me[i][j]);
                }
            }
        }
    }
    // fill reducing ordering on edge midpoints, then the permuted upper
    // triangle in CSC with a fixed slot per assembled entry
    let keys: Vec<(u32, u32)> = acc.keys().copied().collect();
    let mut pts = vec![[0.0; 3]; ndof];
    for (s, &d) in dof.iter().enumerate() {
        if d == usize::MAX {
            continue;
        }
        // slot to its geometric location: edge midpoint or face centroid
        pts[d] = if s < ne * p {
            let e = enodes[s / p];
            sc3(add3(mesh.nodes[e.0 as usize], mesh.nodes[e.1 as usize]), 0.5)
        } else {
            let f = fnodes[(s - ne * p) / nfd];
            let mut c = [0.0; 3];
            for v in f {
                c = add3(c, mesh.nodes[v as usize]);
            }
            sc3(c, 1.0 / 3.0)
        };
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
    let ents: Vec<(usize, Cx, Cx, Cx)> = acc
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
                let mut w = work(ndof, sym);
                for fi in (t..freqs.len()).step_by(nt) {
                    let k0 = 2.0 * std::f64::consts::PI * freqs[fi] / C0;
                    for &(pos, s, m, b) in ents.iter() {
                        ax[pos] = s + m.rs(-k0 * k0) + b.js(k0);
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
                            vq = vq.rs(-ports[q].h / ports[q].area);
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
        for &(pos, s, m, b) in &ents {
            ax[pos] = s + m.rs(-k0 * k0) + b.js(k0);
        }
        let mut w = work(ndof, &sym);
        numeric(ndof, &ap, &ai, &ax, &sym, &mut w);
        let mut rhs = vec![cx(0.0, 0.0); ndof];
        for (&d, &c) in &ports[0].exc {
            rhs[perm[d]] = cx(0.0, -k0 * ETA0 * c / ports[0].w);
        }
        ldsolve(ndof, &sym, &w, &mut rhs);
        let mut ef = Vec::with_capacity(mesh.tets.len());
        let mut bf = vec![];
        for (t, _) in &mesh.tets {
            let (gg, _) = tetgrad([mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]], mesh.nodes[t[3]]]);
            tbasis(p, &gg, [0.25; 4], &mut bf);
            let ld = tdofs(t);
            let mut e = [cx(0.0, 0.0); 3];
            for i in 0..bf.len() {
                if ld[i] != usize::MAX {
                    let x = rhs[perm[ld[i]]];
                    for k in 0..3 {
                        e[k] = e[k] + x.rs(bf[i].0[k]);
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
        o += &format!("CELL_DATA {}\nVECTORS Ere double\n", mesh.tets.len());
        for e in &ef {
            o += &format!("{:.6e} {:.6e} {:.6e}\n", e[0].re, e[1].re, e[2].re);
        }
        o += "VECTORS Eim double\n";
        for e in &ef {
            o += &format!("{:.6e} {:.6e} {:.6e}\n", e[0].im, e[1].im, e[2].im);
        }
        o += "SCALARS Emag double\nLOOKUP_TABLE default\n";
        for e in &ef {
            let m = (e[0].mag().powi(2) + e[1].mag().powi(2) + e[2].mag().powi(2)).sqrt();
            o += &format!("{:.6e}\n", m);
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
