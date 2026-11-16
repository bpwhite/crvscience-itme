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

use std::fs;
use std::path::PathBuf;

use itme_evolution_sim::console::{run_grand_trajectory_headless, run_headless, ConsoleOptions};
use itme_evolution_sim::engine::{SimulationMode, SimulationWorld};
use itme_evolution_sim::geometry::{
    build_structure_scene, render_ascii_3d, render_stage_atlas_svg, Camera3D,
};
use itme_evolution_sim::model::TrophicGuild;
use itme_evolution_sim::organisms::OrganismKind;
use itme_evolution_sim::phylogeny::PhylogenyTracker;
use itme_evolution_sim::regimes::RegimeId;
use itme_evolution_sim::spatial::BiomeZone;

struct TempDirGuard(PathBuf);

impl TempDirGuard {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "crvscience_itme_{}_{}_{:?}",
            tag,
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("create temp test dir");
        Self(dir)
    }
}

impl Drop for TempDirGuard {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn headless_and_grand_trajectory_export_valid_telemetry() {
    let tmp = TempDirGuard::new("headless");
    let csv_path = tmp.0.join("telemetry.csv");
    let atlas_svg_path = tmp.0.join("atlas.svg");
    let newick_path = tmp.0.join("network.enewick");
    let phylo_svg_path = tmp.0.join("phylo.svg");

    let opts = ConsoleOptions {
        organisms: vec![OrganismKind::RnaDnaProtein],
        regimes: vec![RegimeId::R00Static, RegimeId::R30FullyDarwinian],
        steps_per_regime: 6,
        seed: 2026_1113,
        show_ascii_3d: false,
        csv_path: Some(csv_path.to_string_lossy().into_owned()),
        svg_path: Some(atlas_svg_path.to_string_lossy().into_owned()),
        newick_path: Some(newick_path.to_string_lossy().into_owned()),
        phylo_svg_path: Some(phylo_svg_path.to_string_lossy().into_owned()),
        grand_trajectory: false,
    };
    run_headless(&opts).expect("run_headless should succeed");

    let csv_raw = fs::read_to_string(&csv_path).expect("read csv output");
    assert!(csv_raw.starts_with("organism,regime_code,regime_name,step,"));
    assert!(csv_raw.lines().count() >= 1 + 2 * 6);

    let atlas_raw = fs::read_to_string(&atlas_svg_path).expect("read atlas svg");
    assert!(atlas_raw.starts_with("<svg ") && atlas_raw.contains("</svg>"));

    let newick_raw = fs::read_to_string(&newick_path).expect("read newick");
    assert!(newick_raw.contains(';'));

    let phylo_svg_raw = fs::read_to_string(&phylo_svg_path).expect("read phylo svg");
    assert!(phylo_svg_raw.contains("<svg") && phylo_svg_raw.contains("</svg>"));

    let grand_csv = tmp.0.join("grand.csv");
    let grand_opts = ConsoleOptions {
        steps_per_regime: 12,
        seed: 99,
        csv_path: Some(grand_csv.to_string_lossy().into_owned()),
        grand_trajectory: true,
        ..ConsoleOptions::default()
    };
    run_grand_trajectory_headless(&grand_opts).expect("grand trajectory should succeed");
    let grand_raw = fs::read_to_string(&grand_csv).expect("read grand csv");
    assert!(grand_raw.starts_with("epoch,age_ga,tick_end,"));
    assert!(grand_raw.lines().count() >= 8);
}

#[test]
fn geometry_scene_ascii_3d_and_svg_atlas_render_for_all_organisms() {
    let cam = Camera3D::default();
    let mut panels = Vec::new();

    for &org in &[
        OrganismKind::RnaDnaProtein,
        OrganismKind::HominidPopulation,
        OrganismKind::ArtificialIntelligence,
    ] {
        let mut world = SimulationWorld::new(
            org,
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            8,
            42,
        );
        for _ in 0..4 {
            world.step();
        }
        let unit = world.best_unit().expect("population should have best unit");
        let scene = build_structure_scene(unit, RegimeId::R35Symbiogenetic, world.population.len(), 0.5);
        assert!(!scene.nodes_3d.is_empty(), "{org:?} should produce 3D nodes");
        assert!(!scene.nodes_2d.is_empty(), "{org:?} should produce 2D nodes");

        let ascii = render_ascii_3d(&scene, &cam, 48, 16);
        assert!(ascii.contains('+') && ascii.lines().count() >= 18);

        panels.push((
            org.short_title().to_string(),
            format!("nodes={}", scene.nodes_3d.len()),
            scene,
        ));
    }

    let svg = render_stage_atlas_svg("ITME Structural Geometry Verification & <Test>", &panels);
    assert!(svg.starts_with("<svg "));
    assert!(svg.contains("&amp; &lt;Test&gt;"));
    assert!(svg.contains("</svg>"));
}

#[test]
fn phylogeny_tracker_records_horizontal_transfers_and_reticulations() {
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (1, 1, TrophicGuild::Producer, BiomeZone::Alpha, 1.0),
            (2, 2, TrophicGuild::Syntroph, BiomeZone::Beta, 1.25),
        ],
    );

    tracker.record_birth(
        3,
        1,
        Some(2),
        4,
        1,
        TrophicGuild::Syntroph,
        BiomeZone::Alpha,
        1,
        TrophicGuild::Producer,
        0.35,
        0.25,
        1.40,
    );

    let summary = tracker.compute_summary(4);
    assert_eq!(summary.horizontal_transfers, 1);
    assert_eq!(summary.active_reticulations, 1);
    assert!(summary.total_speciation_splits >= 1);

    let enewick = tracker.to_newick(4, true);
    assert!(enewick.contains("#LGT3") && enewick.contains(';'));

    let svg = tracker.to_svg_phylogram("HGT Reticulation Check", 4);
    assert!(svg.contains("<svg") && svg.contains("</svg>"));
}
