//! Ricci flow core (Rust -> WASM). Stage 1 modules land here: mesh, flow (BDF),
//! lm (Levenberg-Marquardt), alexandrov (convex realization).
use wasm_bindgen::prelude::*;

pub mod flow;
pub mod mesh;
pub mod sparse;

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

pub mod hull;
pub mod lm;
pub mod quality;
pub mod realize;
pub mod seedmesh;
