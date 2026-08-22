use std::collections::HashMap;
use std::process::Command;

const ETA0: f64 = 376.73031346177066;
const C0: f64 = 299792458.0;

fn tmp(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

// Structured box mesh in Gmsh 2.2 ASCII, each grid cell split into 6
// tetrahedra (Kuhn split, consistent across cells). Boundary triangles are
// classified by their centroid: sname returns Some(physical name) to emit
// the triangle, None leaves the face natural. vname names the region of a
// tet by its cell center.
fn box_msh(
    name: &str,
    n: [usize; 3],
    l: [f64; 3],
    sname: &dyn Fn(f64, f64, f64) -> Option<&'static str>,
    vname: &dyn Fn(f64, f64, f64) -> &'static str,
) -> std::path::PathBuf {
    let (nx, ny, nz) = (n[0], n[1], n[2]);
    let d = [l[0] / nx as f64, l[1] / ny as f64, l[2] / nz as f64];
    let nid = |v: [usize; 3]| 1 + v[0] + v[1] * (nx + 1) + v[2] * (nx + 1) * (ny + 1);
    let pos = |id: usize| -> [f64; 3] {
        let q = id - 1;
        [(q % (nx + 1)) as f64 * d[0], ((q / (nx + 1)) % (ny + 1)) as f64 * d[1], (q / ((nx + 1) * (ny + 1))) as f64 * d[2]]
    };
    const PERMS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut tets: Vec<([usize; 4], &str)> = vec![];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let nm = vname((i as f64 + 0.5) * d[0], (j as f64 + 0.5) * d[1], (k as f64 + 0.5) * d[2]);
                for p in PERMS {
                    let mut v = [[i, j, k]; 4];
                    for s in 0..3 {
                        v[s + 1] = v[s];
                        v[s + 1][p[s]] += 1;
                    }
                    tets.push(([nid(v[0]), nid(v[1]), nid(v[2]), nid(v[3])], nm));
                }
            }
        }
    }
    let mut count: HashMap<[usize; 3], u32> = HashMap::new();
    for (t, _) in &tets {
        for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut key = [t[f[0]], t[f[1]], t[f[2]]];
            key.sort();
            *count.entry(key).or_insert(0) += 1;
        }
    }
    let mut tris: Vec<([usize; 3], &str)> = vec![];
    for (key, c) in &count {
        if *c != 1 {
            continue;
        }
        let (p0, p1, p2) = (pos(key[0]), pos(key[1]), pos(key[2]));
        let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
        if let Some(nm) = sname(cen[0], cen[1], cen[2]) {
            tris.push((*key, nm));
        }
    }
    let mut ids: HashMap<(i32, &str), usize> = HashMap::new();
    for (_, nm) in &tris {
        let n = ids.len() + 1;
        ids.entry((2, nm)).or_insert(n);
    }
    for (_, nm) in &tets {
        let n = ids.len() + 1;
        ids.entry((3, nm)).or_insert(n);
    }
    let mut s = String::from("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n");
    s += &format!("{}\n", ids.len());
    let mut sorted: Vec<(&(i32, &str), &usize)> = ids.iter().collect();
    sorted.sort_by_key(|(_, id)| **id);
    for ((dim, nm), id) in sorted {
        s += &format!("{} {} \"{}\"\n", dim, id, nm);
    }
    s += "$EndPhysicalNames\n$Nodes\n";
    let nn = (nx + 1) * (ny + 1) * (nz + 1);
    s += &format!("{}\n", nn);
    for id in 1..=nn {
        let p = pos(id);
        s += &format!("{} {} {} {}\n", id, p[0], p[1], p[2]);
    }
    s += "$EndNodes\n$Elements\n";
    s += &format!("{}\n", tris.len() + tets.len());
    let mut eid = 0;
    for (t, nm) in &tris {
        eid += 1;
        s += &format!("{} 2 2 {} 1 {} {} {}\n", eid, ids[&(2, *nm)], t[0], t[1], t[2]);
    }
    for (t, nm) in &tets {
        eid += 1;
        s += &format!("{} 4 2 {} 1 {} {} {} {}\n", eid, ids[&(3, *nm)], t[0], t[1], t[2], t[3]);
    }
    s += "$EndElements\n";
    let path = tmp(name);
    std::fs::write(&path, s).unwrap();
    path
}

// Runs nanofem on a deck, returns (freq, values) rows of the touchstone data
fn run(deck_name: &str, deck: &str) -> Vec<(f64, Vec<f64>)> {
    let path = tmp(deck_name);
    std::fs::write(&path, deck).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_nanofem")).arg(&path).output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#') && !l.starts_with('!'))
        .map(|l| {
            let v: Vec<f64> = l.split_whitespace().map(|s| s.parse().unwrap()).collect();
            (v[0], v[1..].to_vec())
        })
        .collect()
}

fn mag(re: f64, im: f64) -> f64 {
    re.hypot(im)
}

// wrapped phase difference a - b in (-pi, pi]
fn ang_diff(a: f64, b: f64) -> f64 {
    let mut d = a - b;
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d <= -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

// Parallel plate TEM line along x: PEC top and bottom, natural (PMC) side
// walls, one port at each end. Line impedance eta h / w = eta for a square
// cross section.
fn tem_mesh(name: &str, far: &'static str) -> std::path::PathBuf {
    let (lx, lyz) = (0.12, 0.01);
    box_msh(
        name,
        [12, 1, 1],
        [lx, lyz, lyz],
        &move |x, _y, z| {
            if z < 1e-12 || z > lyz - 1e-12 {
                Some("pec")
            } else if x < 1e-12 {
                Some("p1")
            } else if x > lx - 1e-12 {
                Some(far)
            } else {
                None
            }
        },
        &|_, _, _| "air",
    )
}

#[test]
fn tem_matched() {
    let mesh = tem_mesh("tem.msh", "p2");
    let deck = format!(
        "* tem line\nmesh {}\npec pec\nport 1 p1 0 0 1 {}\nport 2 p2 0 0 1 {}\nsweep lin 1e9 1e9 1\n",
        mesh.display(),
        ETA0,
        ETA0
    );
    let rows = run("tem.nfm", &deck);
    let v = &rows[0].1;
    let (s11, s21, s12) = (mag(v[0], v[1]), mag(v[2], v[3]), mag(v[4], v[5]));
    assert!(s11 < 0.05, "S11 = {}", s11);
    assert!((s21 - 1.0).abs() < 0.02, "S21 = {}", s21);
    assert!((s12 - s21).abs() < 1e-9 && (v[3] - v[5]).abs() < 1e-9, "not reciprocal");
    let k0l = 2.0 * std::f64::consts::PI * 1e9 / C0 * 0.12;
    let d = ang_diff(v[3].atan2(v[2]), -k0l);
    assert!(d.abs() < 0.05, "S21 phase error {} rad", d);
}

#[test]
fn tem_dielectric_loss() {
    let (lx, lyz) = (0.12, 0.01);
    let mesh = box_msh(
        "diel.msh",
        [12, 1, 1],
        [lx, lyz, lyz],
        &move |x, _y, z| {
            if z < 1e-12 || z > lyz - 1e-12 {
                Some("pec")
            } else if x < 1e-12 {
                Some("p1")
            } else if x > lx - 1e-12 {
                Some("p2")
            } else {
                None
            }
        },
        &|_, _, _| "diel",
    );
    let z0 = ETA0 / 1.5;
    let deck = format!(
        "mesh {}\nmat diel eps 2.25 tand 0.02\npec pec\nport 1 p1 0 0 1 {}\nport 2 p2 0 0 1 {}\nsweep lin 1e9 1e9 1\n",
        mesh.display(),
        z0,
        z0
    );
    let rows = run("diel.nfm", &deck);
    let v = &rows[0].1;
    let k0 = 2.0 * std::f64::consts::PI * 1e9 / C0;
    assert!(mag(v[0], v[1]) < 0.05, "S11 = {}", mag(v[0], v[1]));
    let expect = (-1.5 * k0 * 0.01 * 0.12).exp();
    assert!((mag(v[2], v[3]) - expect).abs() < 0.01, "S21 = {} expected {}", mag(v[2], v[3]), expect);
    let d = ang_diff(v[3].atan2(v[2]), -1.5 * k0 * 0.12);
    assert!(d.abs() < 0.08, "S21 phase error {} rad", d);
}

#[test]
fn shorted_line() {
    let mesh = tem_mesh("short.msh", "pec");
    let deck = format!("mesh {}\npec pec\nport 1 p1 0 0 1 {}\nsweep lin 1e9 1e9 1\n", mesh.display(), ETA0);
    let rows = run("short.nfm", &deck);
    let v = &rows[0].1;
    let s11 = mag(v[0], v[1]);
    assert!((s11 - 1.0).abs() < 0.005, "S11 = {}, lossless line should reflect all", s11);
    let k0l = 2.0 * std::f64::consts::PI * 1e9 / C0 * 0.12;
    let d = ang_diff(v[1].atan2(v[0]), std::f64::consts::PI - 2.0 * k0l);
    assert!(d.abs() < 0.1, "S11 phase error {} rad", d);
}

#[test]
fn mismatched_port() {
    let mesh = tem_mesh("mis.msh", "p2");
    let deck = format!(
        "mesh {}\npec pec\nport 1 p1 0 0 1 {}\nport 2 p2 0 0 1 {}\nsweep lin 1e9 1e9 1\n",
        mesh.display(),
        2.0 * ETA0,
        ETA0
    );
    let rows = run("mis.nfm", &deck);
    let v = &rows[0].1;
    // line of impedance eta0 matched at the far end, seen through a 2 eta0
    // reference: S11 = (1 - 2) / (1 + 2) = -1/3, no phase
    assert!((v[0] + 1.0 / 3.0).abs() < 0.03 && v[1].abs() < 0.03, "S11 = {} + j {}", v[0], v[1]);
    assert!((v[2] - v[4]).abs() < 1e-9 && (v[3] - v[5]).abs() < 1e-9, "not reciprocal");
}

#[test]
fn abc_absorbs() {
    let mesh = tem_mesh("abc.msh", "abc");
    let deck = format!("mesh {}\npec pec\nabc abc\nport 1 p1 0 0 1 {}\nsweep lin 1e9 1e9 1\n", mesh.display(), ETA0);
    let rows = run("abc.nfm", &deck);
    let s11 = mag(rows[0].1[0], rows[0].1[1]);
    assert!(s11 < 0.05, "S11 = {}, abc should absorb the TEM wave", s11);
}

#[test]
fn pml_absorbs() {
    // TEM line running into a 3 cell PML slab stretched along x, natural
    // (PMC) wall behind it: the wave should not come back
    let (lx, lyz) = (0.15, 0.01);
    let mesh = box_msh(
        "pml.msh",
        [15, 1, 1],
        [lx, lyz, lyz],
        &move |x, _y, z| {
            if z < 1e-12 || z > lyz - 1e-12 {
                Some("pec")
            } else if x < 1e-12 {
                Some("p1")
            } else {
                None
            }
        },
        &|x, _, _| if x > 0.12 { "pml" } else { "air" },
    );
    let deck = format!("mesh {}\npec pec\npml pml 3 0 0\nport 1 p1 0 0 1 {}\nsweep lin 1e9 1e9 1\n", mesh.display(), ETA0);
    let rows = run("pml.nfm", &deck);
    let s11 = mag(rows[0].1[0], rows[0].1[1]);
    assert!(s11 < 0.08, "S11 = {}, pml should absorb the TEM wave", s11);
}

#[test]
fn cavity_resonance() {
    // PEC box 0.1 x 0.05 x 0.1, aperture port in the x = 0 wall. The 101
    // mode resonates at c0 / 2 * sqrt(2) / 0.1 = 2.12 GHz and shows up as a
    // peak of the port input impedance.
    let (a, b, c) = (0.1, 0.05, 0.1);
    let mesh = box_msh(
        "cav.msh",
        [8, 3, 8],
        [a, b, c],
        &move |x, y, z| {
            if x < 1e-12 && (y - b / 2.0).abs() < b / 6.0 && (z - c / 2.0).abs() < c / 8.0 + 1e-12 {
                Some("feed")
            } else {
                Some("pec")
            }
        },
        &|_, _, _| "air",
    );
    let deck = format!(
        "mesh {}\npec pec\nport 1 feed 0 1 0 50\nsweep lin 1.95e9 2.35e9 29\n",
        mesh.display()
    );
    let rows = run("cav.nfm", &deck);
    let mut best = (0.0, 0.0);
    for (f, v) in &rows {
        // Z = z0 (1 + S) / (1 - S)
        let (nr, ni, dr, di) = (1.0 + v[0], v[1], 1.0 - v[0], -v[1]);
        let z = 50.0 * mag(nr, ni) / mag(dr, di);
        if z > best.1 {
            best = (*f, z);
        }
    }
    let f101 = C0 / 2.0 * (2.0f64).sqrt() / 0.1;
    assert!((best.0 - f101).abs() < 0.08e9, "resonance at {} GHz, expected {} GHz", best.0 / 1e9, f101 / 1e9);
}

// A conductive filling makes the line lossy through eps - j sigma / w eps0,
// so the propagation constant is k0 sqrt(eps_c) and S21 follows in closed
// form. Checks that conduction current lands in the right frequency slot.
#[test]
fn substrate_conductivity() {
    let (lx, lyz) = (0.12, 0.01);
    let mesh = box_msh(
        "sig.msh",
        [12, 1, 1],
        [lx, lyz, lyz],
        &move |x, _y, z| {
            if z < 1e-12 || z > lyz - 1e-12 {
                Some("pec")
            } else if x < 1e-12 {
                Some("p1")
            } else if x > lx - 1e-12 {
                Some("p2")
            } else {
                None
            }
        },
        &|_, _, _| "lossy",
    );
    let (f, sigma, er) = (1e9, 0.01, 1.0);
    let (w, eps0) = (2.0 * std::f64::consts::PI * f, 8.8541878128e-12);
    // eps_c = er - j sigma / (w eps0), gamma = j k0 sqrt(eps_c)
    let (re, im) = (er, -sigma / (w * eps0));
    let r = (re * re + im * im).sqrt().sqrt();
    let th = im.atan2(re) / 2.0;
    let (nr, ni) = (r * th.cos(), r * th.sin());
    let k0 = w / C0;
    let expect = (k0 * ni * lx).exp();
    let z0 = ETA0 / nr;
    let deck = format!(
        "mesh {}\nmat lossy eps {} sigma {}\npec pec\nport 1 p1 0 0 1 {}\nport 2 p2 0 0 1 {}\nsweep lin {} {} 1\n",
        mesh.display(),
        er,
        sigma,
        z0,
        z0,
        f,
        f
    );
    let v = run("sig.nfm", &deck)[0].1.clone();
    let got = mag(v[2], v[3]);
    assert!((got - expect).abs() < 0.01, "S21 = {}, analytic {}", got, expect);
    assert!(got < 0.95, "conductivity produced no measurable loss");
}

// Lossy plates instead of PEC: the parallel plate line then attenuates by
// Rs / (eta h) per unit length with Rs the surface resistance. Checks the
// sqrt(k0) slot, which is where the skin effect lives.
#[test]
fn conductor_loss() {
    let (lx, lyz) = (0.12, 0.01);
    let mesh = box_msh(
        "met.msh",
        [12, 1, 1],
        [lx, lyz, lyz],
        &move |x, _y, z| {
            if z < 1e-12 || z > lyz - 1e-12 {
                Some("cu")
            } else if x < 1e-12 {
                Some("p1")
            } else if x > lx - 1e-12 {
                Some("p2")
            } else {
                None
            }
        },
        &|_, _, _| "air",
    );
    let (f, sigma) = (1e9, 1e4);
    let (w, mu0) = (2.0 * std::f64::consts::PI * f, 4e-7 * std::f64::consts::PI);
    let rs = (w * mu0 / (2.0 * sigma)).sqrt();
    let expect = (-rs / (ETA0 * lyz) * lx).exp();
    let deck = format!(
        "mesh {}\nmetal cu {}\nport 1 p1 0 0 1 {}\nport 2 p2 0 0 1 {}\nsweep lin {} {} 1\n",
        mesh.display(),
        sigma,
        ETA0,
        ETA0,
        f,
        f
    );
    let v = run("met.nfm", &deck)[0].1.clone();
    let got = mag(v[2], v[3]);
    assert!((got - expect).abs() < 0.01, "S21 = {}, analytic {}", got, expect);
    assert!(got < 0.99, "metal loss produced no measurable attenuation");
}

// Runs nanofem and returns stderr, expecting it to refuse the input. A
// panic message counts as a failure: bad input must produce a diagnostic,
// never a backtrace.
fn refuses(name: &str, deck: &str) -> String {
    let path = tmp(name);
    std::fs::write(&path, deck).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_nanofem")).arg(&path).output().unwrap();
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(!err.contains("panicked"), "{} panicked: {}", name, err);
    assert!(!out.status.success(), "{} was accepted, stdout: {}", name, String::from_utf8_lossy(&out.stdout));
    err
}

// Values that are meaningless out of range, and roles put on a group of the
// wrong dimension. The latter used to run happily and return a plausible
// answer for a completely different model.
#[test]
fn refuses_bad_decks() {
    let mesh = tem_mesh("valid.msh", "p2");
    let m = format!("mesh {}\n", mesh.display());
    let cases: [(&str, &str, &str); 8] = [
        ("z0neg", "pec pec\nport 1 p1 0 0 1 -50\nsweep lin 1e9 1e9 1\n", "impedance must be positive"),
        ("z0zero", "pec pec\nport 1 p1 0 0 1 0\nsweep lin 1e9 1e9 1\n", "impedance must be positive"),
        ("mur0", "mat air eps 1 mur 0\npec pec\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", "mur must be positive"),
        ("epsneg", "mat air eps -2\npec pec\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", "eps must be positive"),
        ("tandneg", "mat air eps 1 tand -0.1\npec pec\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", "tand must not be negative"),
        ("freq0", "pec pec\nport 1 p1 0 0 1 50\nsweep lin 0 0 1\n", "frequency must be positive"),
        ("pecvol", "pec air\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", "no boundary triangles"),
        ("matsurf", "mat pec eps 2\npec pec\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", "no tetrahedra"),
    ];
    for (name, body, want) in cases {
        let err = refuses(&format!("{}.nfm", name), &format!("{}{}", m, body));
        assert!(err.contains(want), "{}: expected '{}', got '{}'", name, want, err.trim());
    }
}

// Truncated mesh lines must be reported, not indexed past the end.
#[test]
fn refuses_bad_meshes() {
    let head = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    let cases: [(&str, &str, &str); 4] = [
        ("shortnode", "$Nodes\n1\n1 0 0\n$EndNodes\n", "bad node line"),
        ("shortelem", "$Nodes\n1\n1 0 0 0\n$EndNodes\n$Elements\n1\n1 4\n$EndElements\n", "bad element line"),
        ("fewnodes", "$Nodes\n1\n1 0 0 0\n$EndNodes\n$Elements\n1\n1 4 2 1 1 1 2\n$EndElements\n", "fewer nodes"),
        ("shortphys", "$PhysicalNames\n1\n3\n$EndPhysicalNames\n", "$PhysicalNames entry"),
    ];
    for (name, body, want) in cases {
        let mp = tmp(&format!("{}.msh", name));
        std::fs::write(&mp, format!("{}{}", head, body)).unwrap();
        let deck = format!("mesh {}\npec pec\nport 1 p1 0 0 1 50\nsweep lin 1e9 1e9 1\n", mp.display());
        let err = refuses(&format!("{}.nfm", name), &deck);
        assert!(err.contains(want), "{}: expected '{}', got '{}'", name, want, err.trim());
    }
}

// The budget covers the solver in src/main.rs alone. This test file and the
// models in models/ are outside it.
#[test]
fn loc_budget() {
    let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/main.rs")).unwrap();
    let n = src
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with("//")
        })
        .count();
    assert!(n <= 1000, "src/main.rs has {} LOC, budget is 1000", n);
}
