// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

#![allow(dead_code)]
//! ========================================================================================
//! ITME v23.0 (`v2.3.0`) Native Rust Simulation — Standalone Single-Entrypoint Wrapper
//! Deduplicated via `#[path = "src/..."]` against the canonical `src/` modules so the
//! zero-dependency `rustc` build stays 100% synchronized with the Cargo library & CLI.
//!
//! Full Cargo workspace (with native `eframe` GUI + integration tests):
//!     cargo run --release
//!     cargo run --release -- --headless --export-newick net.nhx --export-phylo-svg net.svg
//!
//! Zero-dependency standalone `rustc` compilation (headless CLI mode):
//!     rustc --edition 2021 -O itme_sim_standalone.rs -o itme_sim && ./itme_sim --headless
//! ========================================================================================

#[path = "src/model.rs"]
pub mod model;
#[path = "src/spatial.rs"]
pub mod spatial;
#[path = "src/phylogeny.rs"]
pub mod phylogeny;
#[path = "src/regimes.rs"]
pub mod regimes;
#[path = "src/organisms.rs"]
pub mod organisms;
#[path = "src/geometry.rs"]
pub mod geometry;
#[path = "src/engine.rs"]
pub mod engine;
#[path = "src/console.rs"]
pub mod console;
#[cfg(feature = "gui")]
#[path = "src/gui.rs"]
pub mod gui;

use crate::console::{run_headless, ConsoleOptions};
use crate::organisms::OrganismKind;
use crate::regimes::RegimeId;
use std::env;
use std::process;

fn print_help() {
    println!(
        "\
ITME v23.0 (v2.3.0) Continuous 4.25 Ga Hadean-to-Gaia Multi-Epoch Grand Trajectory & Full ETI Engine
Standalone zero-dependency CLI wrapper sharing canonical `src/` simulation modules.

USAGE:
    itme_sim --headless [OPTIONS]
    itme_sim --grand-trajectory [--steps <N>] [--export-csv <FILE>]

OPTIONS:
    --headless, --console         Run in terminal console mode (no display server required)
    --grand-trajectory            Run the 8-epoch 4.25 Ga Hadean-to-Gaia Grand Trajectory
    --gui                         Launch native GUI (when compiled with `--cfg feature=\"gui\"`)
    --organism <SLUG>             Select `rna`, `hominid`, `ai`, or `all` [default: all]
    --regime <CODE>               Select `all` (13 regimes) or a specific code [default: all]
    --steps <N>                   Simulation ticks per regime in headless mode [default: 24]
    --seed <U64>                  Deterministic Xoshiro256** PRNG seed [default: 20250813]
    --no-ascii-3d                 Omit inline 3D ASCII & spatial map projections in console output
    --export-csv <FILE>           Export 67-column step-by-step telemetry (or grand trajectory) CSV
    --export-svg <FILE>           Export 3D/2D projected structural stage atlas to SVG
    --export-newick <FILE>        Export Extended Newick (eNewick #H + MRO/APICO/KLEPTO NHX) networks
    --export-phylo-svg <FILE>     Export anastomosing SVG phylogram
    -h, --help                    Print this help message"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut headless = false;
    let mut force_gui = false;
    let mut opts = ConsoleOptions::default();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print_help();
                return;
            }
            "--headless" | "--console" => {
                headless = true;
            }
            "--grand-trajectory" => {
                headless = true;
                opts.grand_trajectory = true;
            }
            "--gui" => {
                force_gui = true;
            }
            "--no-ascii-3d" => {
                opts.show_ascii_3d = false;
            }
            "--organism" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --organism requires a value (`rna`, `hominid`, `ai`, or `all`).");
                    process::exit(2);
                };
                if val.eq_ignore_ascii_case("all") {
                    opts.organisms = OrganismKind::ALL.to_vec();
                } else if let Some(org) = OrganismKind::from_slug(val) {
                    opts.organisms = vec![org];
                } else {
                    eprintln!("Error: unknown organism '{val}'. Expected `rna`, `hominid`, `ai`, or `all`.");
                    process::exit(2);
                }
            }
            "--regime" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --regime requires a regime code or `all`.");
                    process::exit(2);
                };
                if val.eq_ignore_ascii_case("all") {
                    opts.regimes = RegimeId::ALL.to_vec();
                } else if let Some(reg) = RegimeId::from_slug(val) {
                    opts.regimes = vec![reg];
                } else {
                    eprintln!("Error: unknown regime '{val}'. Run `--help` for valid codes.");
                    process::exit(2);
                }
            }
            "--steps" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --steps requires an integer.");
                    process::exit(2);
                };
                opts.steps_per_regime = val.parse::<usize>().unwrap_or_else(|_| {
                    eprintln!("Error: invalid --steps integer '{val}'.");
                    process::exit(2);
                });
            }
            "--seed" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --seed requires an unsigned integer.");
                    process::exit(2);
                };
                opts.seed = val.parse::<u64>().unwrap_or_else(|_| {
                    eprintln!("Error: invalid --seed integer '{val}'.");
                    process::exit(2);
                });
            }
            "--export-csv" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --export-csv requires a file path.");
                    process::exit(2);
                };
                opts.csv_path = Some(val.clone());
            }
            "--export-svg" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --export-svg requires a file path.");
                    process::exit(2);
                };
                opts.svg_path = Some(val.clone());
            }
            "--export-newick" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --export-newick requires a file path.");
                    process::exit(2);
                };
                opts.newick_path = Some(val.clone());
            }
            "--export-phylo-svg" => {
                i += 1;
                let Some(val) = args.get(i) else {
                    eprintln!("Error: --export-phylo-svg requires a file path.");
                    process::exit(2);
                };
                opts.phylo_svg_path = Some(val.clone());
            }
            other => {
                eprintln!("Unknown flag '{other}'. Run with `--help` for usage.");
                process::exit(2);
            }
        }
        i += 1;
    }

    if headless && !force_gui {
        if let Err(err) = run_headless(&opts) {
            eprintln!("Headless simulation error: {err}");
            process::exit(1);
        }
        return;
    }

    #[cfg(feature = "gui")]
    {
        if let Err(err) = crate::gui::run_gui(opts.seed) {
            eprintln!(
                "Notice: Native GUI could not initialize ({err}). Falling back to headless console mode...\n"
            );
            if let Err(e2) = run_headless(&opts) {
                eprintln!("Headless fallback failed: {e2}");
                process::exit(1);
            }
        }
    }

    #[cfg(not(feature = "gui"))]
    {
        let _ = force_gui;
        println!("Standalone zero-dependency binary built without `gui` feature; executing headless console mode.");
        if let Err(err) = run_headless(&opts) {
            eprintln!("Headless simulation error: {err}");
            process::exit(1);
        }
    }
}
