//! Test-only loader for `fixtures/<name>/arrays.npz` (golden data from the
//! Python /flow pipeline).  Not part of the library build: included by
//! integration tests via `#[path = "../src/fixtures.rs"] mod fixtures;`
//! (npyz + zip[deflate] are dev-dependencies; npyz's own `npz` feature
//! pulls bzip2-sys/getrandom, which need a C toolchain we do not have).
#![allow(dead_code)]
use std::path::PathBuf;

pub struct Arr<T> {
    pub shape: Vec<usize>,
    pub data: Vec<T>,
}

pub struct Fixture {
    pub name: String,
    pub dir: PathBuf,
}

pub fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("fixtures")
}

impl Fixture {
    pub fn open(name: &str) -> Fixture {
        Fixture { name: name.to_string(), dir: fixtures_dir().join(name) }
    }
    fn read<T: npyz::Deserialize>(&self, key: &str) -> Arr<T> {
        let file = std::fs::File::open(self.dir.join("arrays.npz")).expect("open npz");
        let mut ar = zip::ZipArchive::new(file).expect("zip");
        let mut ent = ar.by_name(&format!("{key}.npy")).unwrap_or_else(|_| panic!("missing {key}"));
        let mut bytes = Vec::new();
        std::io::Read::read_to_end(&mut ent, &mut bytes).expect("inflate");
        let f = npyz::NpyFile::new(&bytes[..]).expect("npy header");
        let shape = f.shape().iter().map(|&s| s as usize).collect();
        Arr { shape, data: f.into_vec::<T>().expect("decode") }
    }
    pub fn f64(&self, key: &str) -> Arr<f64> {
        self.read(key)
    }
    pub fn i64(&self, key: &str) -> Arr<i64> {
        self.read(key)
    }
    pub fn json(&self) -> String {
        std::fs::read_to_string(self.dir.join("fixture.json")).unwrap()
    }
}
