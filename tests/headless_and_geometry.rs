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

    let grand_csv = tmp.0.join("nested_csv/sub/grand.csv");
    let grand_svg = tmp.0.join("nested_svg/sub/grand_atlas.svg");
    let grand_newick = tmp.0.join("nested_nhx/sub/grand.enewick");
    let grand_phylo_svg = tmp.0.join("nested_phylo/sub/grand_phylo.svg");
    let grand_opts = ConsoleOptions {
        steps_per_regime: 12,
        seed: 99,
        csv_path: Some(grand_csv.to_string_lossy().into_owned()),
        svg_path: Some(grand_svg.to_string_lossy().into_owned()),
        newick_path: Some(grand_newick.to_string_lossy().into_owned()),
        phylo_svg_path: Some(grand_phylo_svg.to_string_lossy().into_owned()),
        grand_trajectory: true,
        ..ConsoleOptions::default()
    };
    run_grand_trajectory_headless(&grand_opts).expect("grand trajectory should succeed");
    let grand_raw = fs::read_to_string(&grand_csv).expect("read grand csv");
    assert!(grand_raw.starts_with("epoch,age_ga,tick_end,"));
    assert!(grand_raw.lines().count() >= 8);
    let grand_svg_raw = fs::read_to_string(&grand_svg).expect("read grand atlas svg");
    assert!(grand_svg_raw.starts_with("<svg ") && grand_svg_raw.contains("</svg>"));
    let grand_nhx_raw = fs::read_to_string(&grand_newick).expect("read grand newick");
    assert!(grand_nhx_raw.contains(';'));
    let grand_phylo_raw = fs::read_to_string(&grand_phylo_svg).expect("read grand phylo svg");
    assert!(grand_phylo_raw.contains("<svg") && grand_phylo_raw.contains("</svg>"));
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

#[test]
fn biosphere_coupling_vs_domain_sandbox_isolation() {
    for &org in &OrganismKind::ALL {
        let isolated = RegimeId::R35Symbiogenetic.domain_preset_controls(org);
        assert_eq!(
            isolated.biogeochemical_feedback_gain, 0.0,
            "{org:?} isolated preset must zero biogeochemical_feedback_gain"
        );
        assert_eq!(
            isolated.planetary_stewardship_coupling, 0.0,
            "{org:?} isolated preset must zero planetary_stewardship_coupling"
        );
        assert_eq!(
            isolated.extractive_burn_intensity, 0.0,
            "{org:?} isolated preset must zero extractive_burn_intensity"
        );

        let mut world = SimulationWorld::new(
            org,
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            12,
            2026_0505,
        );
        assert!(!world.coupled_biosphere_mode);
        world.controls = isolated;
        let o2_before = world.biosphere.atmospheric_o2;
        let co2_before = world.biosphere.atmospheric_co2_ppm;
        let gibbs_before = world.biosphere.chemical_disequilibrium_bits;
        for _ in 0..8 {
            world.step();
        }
        assert!(
            (world.biosphere.atmospheric_o2 - o2_before).abs() < 1e-12
                && (world.biosphere.atmospheric_co2_ppm - co2_before).abs() < 1e-12
                && (world.biosphere.chemical_disequilibrium_bits - gibbs_before).abs() < 1e-12,
            "{org:?} isolated sandbox must not mutate planetary biosphere reservoirs"
        );
    }

    let mut coupled = SimulationWorld::new(
        OrganismKind::ArtificialIntelligence,
        RegimeId::R35Symbiogenetic,
        SimulationMode::OrganicSandbox,
        12,
        2026_0505,
    );
    coupled.coupled_biosphere_mode = true;
    coupled.controls = RegimeId::R35Symbiogenetic.preset_controls();
    coupled.controls.biogeochemical_feedback_gain = 0.85;
    coupled.controls.planetary_stewardship_coupling = 0.88;
    let o2_start = coupled.biosphere.atmospheric_o2;
    let co2_start = coupled.biosphere.atmospheric_co2_ppm;
    for _ in 0..10 {
        coupled.step();
    }
    assert!(
        (coupled.biosphere.atmospheric_o2 - o2_start).abs() > 1e-4
            || (coupled.biosphere.atmospheric_co2_ppm - co2_start).abs() > 1e-2,
        "Coupled PlanetaryBiosphere mode must actively step O2 and CO2 reservoirs"
    );
    let m = coupled.latest_metrics();
    assert!((m.planetary_o2_level - coupled.biosphere.atmospheric_o2).abs() < 1e-12);
    assert!((m.planetary_co2_ppm - coupled.biosphere.atmospheric_co2_ppm).abs() < 1e-12);
}

#[test]
fn grand_trajectory_uses_unclamped_raw_mi_and_li() {
    let mut world = SimulationWorld::new(
        OrganismKind::ArtificialIntelligence,
        RegimeId::R35Symbiogenetic,
        SimulationMode::OrganicSandbox,
        12,
        20280101,
    );
    let epochs = world.run_continuous_grand_trajectory(10);
    assert_eq!(epochs.len(), 8);

    // Verify Epoch 8 reflects the exact final recorded metrics without any synthetic floor clamp:
    let final_m = world.latest_metrics();
    let final_ep = epochs.last().unwrap();
    assert!(
        (final_ep.heritable_mutual_info_bits - final_m.heritable_mutual_info_bits).abs() < 1e-12,
        "Grand trajectory MI ({}) must match genuine raw latest_metrics ({})",
        final_ep.heritable_mutual_info_bits,
        final_m.heritable_mutual_info_bits
    );
    assert!(
        (final_ep.informational_leverage_li - final_m.leverage_li).abs() < 1e-12,
        "Grand trajectory L_I ({}) must match genuine raw latest_metrics ({})",
        final_ep.informational_leverage_li,
        final_m.leverage_li
    );
}

#[test]
fn shared_csv_and_stage_atlas_svg_exports_are_well_formed() {
    let mut world = SimulationWorld::new(
        OrganismKind::RnaDnaProtein,
        RegimeId::R30FullyDarwinian,
        SimulationMode::OrganicSandbox,
        12,
        777,
    );
    for _ in 0..5 {
        world.step();
    }
    let hist_csv = world.export_history_csv();
    let mut lines = hist_csv.lines();
    let hdr = lines.next().expect("CSV header present");
    assert!(hdr.starts_with("organism,regime_code,regime_name,step,"));
    assert_eq!(hdr.split(',').count(), 67, "Telemetry CSV must have 67 columns");
    assert_eq!(
        lines.count(),
        world.history.len(),
        "Telemetry CSV must emit one row per history record"
    );

    let epochs = world.run_continuous_grand_trajectory(6);
    let grand_csv = SimulationWorld::export_grand_trajectory_csv(&epochs);
    let mut glines = grand_csv.lines();
    let ghdr = glines.next().expect("Grand CSV header present");
    assert!(ghdr.starts_with("epoch,age_ga,tick_end,"));
    assert_eq!(glines.count(), 8, "Grand trajectory CSV must contain 8 epoch rows");

    let atlas_svg = world.build_stage_atlas_svg();
    assert!(atlas_svg.starts_with("<svg ") && atlas_svg.contains("</svg>"));
    assert!(
        atlas_svg.contains("R(0.0)") && atlas_svg.contains("R(3.5-SYM)"),
        "Stage Atlas SVG must include all 13 regime panels"
    );
}
