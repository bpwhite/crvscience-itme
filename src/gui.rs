//! Native desktop GUI mode built with `eframe` / `egui` (`0.29`) — ITME v12.0 (`v1.2.0`).
//!
//! Features live real-time **Anaerobic Mitochondrial Reduction (`Mito->Hydro->MitoS->MRO->0`) + Apicoplasts + Kleptoplasty**:
//! - `MRO Reduce!`, `Apico Reduce!`, `Klepto Steal!`, `Symbio 1°!`, `Serial 2°/3°!`, `EGT Step!`, and `NM Reduce!` controls,
//! - `O2_level (Hypoxia)`, `mu_mro (Anaerobic MRO Reduction)`, and `suf_rescue (Bacterial Cytosolic SUF)` sliders,
//! - Teal-gold `Hydrogenosome`s, gold-rimmed `ISC` `Mitosome`s, and cyan `SUF` `Amiteukaryote`s on the 2D grid & phylogram!

use crate::engine::{SimulationMode, SimulationWorld};
use crate::geometry::{Camera3D, StructureScene};
use crate::model::{AbiogenesisStage, AgricultureState, ChiralityCodeStage, CulturalState, EusocialState, ExosomaticState, GenomeSubstrateStage, MachineSynthesisState, MroClass, MulticellularState, OrderPhase, PlanetaryState, TrophicGuild};
use crate::organisms::OrganismKind;
use crate::regimes::RegimeId;
use eframe::{egui, App};
use std::fs;
use std::time::{Duration, Instant};

const SPECIES_PALETTE: [[u8; 3]; 8] = [
    [90, 235, 255],
    [255, 168, 76],
    [130, 245, 150],
    [255, 108, 178],
    [195, 145, 255],
    [255, 232, 98],
    [95, 165, 255],
    [255, 118, 108],
];

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum ViewportLayout {
    SpatialAndPhylogeny,
    SpatialGridAnd3D,
    PhylogenyFull,
    Perspective3D,
    Blueprint2D,
    TriOrganism3D,
}

pub struct ItmeGuiApp {
    worlds: [SimulationWorld; 3],
    active_org_idx: usize,
    viewport_layout: ViewportLayout,
    camera: Camera3D,
    auto_rotate: bool,
    show_labels: bool,
    paused: bool,
    steps_per_second: f32,
    last_step_instant: Instant,
    export_status: Option<String>,
}

impl ItmeGuiApp {
    pub fn new(cc: &eframe::CreationContext<'_>, seed: u64) -> Self {
        let mut visuals = egui::Visuals::dark();
        visuals.panel_fill = egui::Color32::from_rgb(13, 18, 32);
        visuals.window_fill = egui::Color32::from_rgb(16, 22, 38);
        cc.egui_ctx.set_visuals(visuals);

        let worlds = [
            SimulationWorld::new(
                OrganismKind::RnaDnaProtein,
                RegimeId::R00Static,
                SimulationMode::GuidedCurriculum,
                32,
                seed,
            ),
            SimulationWorld::new(
                OrganismKind::HominidPopulation,
                RegimeId::R00Static,
                SimulationMode::GuidedCurriculum,
                32,
                seed + 101,
            ),
            SimulationWorld::new(
                OrganismKind::ArtificialIntelligence,
                RegimeId::R00Static,
                SimulationMode::GuidedCurriculum,
                32,
                seed + 202,
            ),
        ];

        Self {
            worlds,
            active_org_idx: 0,
            viewport_layout: ViewportLayout::SpatialAndPhylogeny,
            camera: Camera3D::default(),
            auto_rotate: true,
            show_labels: true,
            paused: false,
            steps_per_second: 8.0,
            last_step_instant: Instant::now(),
            export_status: None,
        }
    }

    fn active_world(&self) -> &SimulationWorld {
        &self.worlds[self.active_org_idx]
    }

    fn active_world_mut(&mut self) -> &mut SimulationWorld {
        &mut self.worlds[self.active_org_idx]
    }

    fn step_all_worlds(&mut self) {
        for w in &mut self.worlds {
            w.step();
        }
    }

    fn jump_all_or_active_to_regime(&mut self, regime: RegimeId) {
        if self.viewport_layout == ViewportLayout::TriOrganism3D {
            for w in &mut self.worlds {
                let preserve = regime.index() >= 2;
                w.apply_regime_preset(regime, preserve);
                w.step();
            }
        } else {
            let preserve = regime.index() >= 2;
            let w = self.active_world_mut();
            w.apply_regime_preset(regime, preserve);
            w.step();
        }
    }

    fn export_snapshot_files(&mut self) {
        let w = self.active_world();
        let slug = w.organism.slug();
        let step = w.step_count;
        let nhx_file = format!("itme_{slug}_step{step}.nhx");
        let svg_file = format!("itme_{slug}_step{step}_phylogram.svg");
        let nhx_res = fs::write(&nhx_file, w.export_newick(true));
        let svg_res = fs::write(&svg_file, w.export_phylogram_svg());
        self.export_status = Some(match (nhx_res, svg_res) {
            (Ok(()), Ok(())) => format!("Saved {nhx_file} & {svg_file}"),
            (Err(e), _) | (_, Err(e)) => format!("Export error: {e}"),
        });
    }
}

impl App for ItmeGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if !self.paused {
            let interval = Duration::from_secs_f32((1.0 / self.steps_per_second).max(0.015));
            if self.last_step_instant.elapsed() >= interval {
                if self.viewport_layout == ViewportLayout::TriOrganism3D {
                    self.step_all_worlds();
                } else {
                    self.active_world_mut().step();
                }
                self.last_step_instant = Instant::now();
            }
            if self.auto_rotate {
                self.camera.yaw += 0.008;
            }
            ctx.request_repaint_after(Duration::from_millis(16));
        } else if self.auto_rotate {
            self.camera.yaw += 0.006;
            ctx.request_repaint_after(Duration::from_millis(16));
        }

        self.draw_top_toolbar(ctx);
        self.draw_left_controls(ctx);
        self.draw_right_inspector(ctx);
        self.draw_bottom_telemetry_plots(ctx);
        self.draw_central_viewport(ctx);
    }
}

impl ItmeGuiApp {
    fn draw_top_toolbar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("itme_top_bar").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.heading(
                    egui::RichText::new(
                        "ITME v23.0 — Continuous 4.25 Ga Grand Trajectory (Homochirality -> SGC -> T-DNA -> LUCA -> LECA -> Archives -> Gaia-2.0)",
                    )
                    .color(egui::Color32::from_rgb(130, 245, 205))
                    .strong(),
                );
                ui.separator();

                let play_label = if self.paused { "Play" } else { "Pause" };
                if ui.button(play_label).clicked() {
                    self.paused = !self.paused;
                    self.last_step_instant = Instant::now();
                }
                if ui.button("Step").clicked() {
                    if self.viewport_layout == ViewportLayout::TriOrganism3D {
                        self.step_all_worlds();
                    } else {
                        self.active_world_mut().step();
                    }
                }
                if ui.button("4.25 Ga -> Gaia-2.0 Grand Run!").clicked() {
                    let epochs = self.active_world_mut().run_continuous_grand_trajectory(10);
                    if let Some(last_ep) = epochs.last() {
                        self.export_status = Some(format!(
                            "Completed 8-Epoch 4.25 Ga Grand Trajectory -> {}: L_I = {:.2}, MI = {:.1}b, O2 = {:.2}, Exergy = {:.2}!",
                            last_ep.epoch.code(),
                            last_ep.informational_leverage_li,
                            last_ep.heritable_mutual_info_bits,
                            last_ep.planetary_o2,
                            last_ep.planetary_exergy_budget
                        ));
                    }
                }
                if ui.button("RNR RNA->T-DNA + Virosphere!").clicked() {
                    let (rn, td, kb) =
                        self.active_world_mut().trigger_rna_to_dna_virosphere_pulse();
                    self.export_status = Some(format!(
                        "RNR & Virosphere: +{rn} U-DNA Transitions, +{td} Megabase T-DNA/UDG Genomes (Cap = {kb:.0} kb)!"
                    ));
                }
                if ui.button("Homochiral D-RNA/L-AA + SGC!").clicked() {
                    let (hm, sgc, mp) =
                        self.active_world_mut().trigger_homochiral_sgc_pulse();
                    self.export_status = Some(format!(
                        "Homochirality & SGC: +{hm} D-RNA Symmetry Breaks, +{sgc} Universal 20-AA SGC Lineages (MP = {:.1}%)!",
                        mp * 100.0
                    ));
                }
                if ui.button("Hadean Abiogenesis -> LUCA!").clicked() {
                    let (rf, lc, eig) =
                        self.active_world_mut().trigger_hadean_abiogenesis_pulse();
                    self.export_status = Some(format!(
                        "Hadean Abiogenesis: +{rf} RAF Protocells, +{lc} LUCA Lineages (Eigen L_max = {eig} nt)!"
                    ));
                }
                if ui.button("Gaia-2.0 Planetary Awakening!").clicked() {
                    let (fm, aw, ex) =
                        self.active_world_mut().trigger_planetary_awakening_pulse();
                    self.export_status = Some(format!(
                        "Gaia-2.0 Homeostasis: +{fm} Formal Symbionts, +{aw} Planetary Stewards, Planetary Exergy Budget = {ex:.2}!"
                    ));
                }
                if ui.button("Formal AI Symbiont + Proof Check!").clicked() {
                    let (ms, fm, prf) =
                        self.active_world_mut().trigger_machine_synthesis_pulse();
                    self.export_status = Some(if ms > 0 || fm > 0 || prf > 0 {
                        format!(
                            "Machine-Verified Synthesis: +{ms} Computational Symbionts, +{fm} Formal-Evaluator Symbioses, {prf} proofs checked!"
                        )
                    } else {
                        "All lineages already at MachineVerifiedSymbiosis (Tails >= 96%)!".to_string()
                    });
                }
                if ui.button("Written Archive + Peer Verify!").clicked() {
                    let (ar, inst, aud) =
                        self.active_world_mut().trigger_exosomatic_archive_pulse();
                    self.export_status = Some(if ar > 0 || inst > 0 || aud > 0 {
                        format!(
                            "Exosomatic Epistemology: +{ar} Written Archives, +{inst} Peer-Verified Institutions, {aud} replication audits!"
                        )
                    } else {
                        "All lineages already at PeerVerifiedInstitution (40% cognitive offloading)!".to_string()
                    });
                }
                if ui.button("Symbolic Culture + Deme Recomb!").clicked() {
                    let (rt, sym, rec) =
                        self.active_world_mut().trigger_cultural_ratchet_pulse();
                    self.export_status = Some(if rt > 0 || sym > 0 || rec > 0 {
                        format!(
                            "Cumulative Culture: +{rt} Ratchet, +{sym} Symbolic Teaching lineages, +{rec} Cross-Deme Recombinations!"
                        )
                    } else {
                        "Set Teaching Fidelity >= 0.88 and Partial Connectivity (0.15..0.72) for Cross-Deme Symbolic Culture!".to_string()
                    });
                }
                if ui.button("Agro-Crop + Actino Shield!").clicked() {
                    let (ag, ob, sh) =
                        self.active_world_mut().trigger_agriculture_domestication_pulse();
                    self.export_status = Some(if ag > 0 || ob > 0 {
                        format!(
                            "Agricultural Ectosymbiosis: +{ag} farming colonies, +{ob} clonal Gongylidia/LeafCutter monocultures, {sh} Pseudonocardia shields!"
                        )
                    } else {
                        "All superorganisms already at 8n LeafCutterAgroEcosystem!".to_string()
                    });
                }
                if ui.button("Eusocial Superorg!").clicked() {
                    let (eu, sup, pol) =
                        self.active_world_mut().trigger_superorganism_transition_pulse();
                    self.export_status = Some(if sup > 0 || eu > 0 {
                        format!(
                            "Superorganism ETI: +{eu} Eusocial colonies, +{sup} Obligate Superorganisms, {pol} rebel eggs policed!"
                        )
                    } else {
                        "Raise Queen Monogamy >= 0.75 (Monogamy Window) on Germ-Soma organisms to cross Point of No Return!".to_string()
                    });
                }
                if ui.button("Multicell + Germ-Soma!").clicked() {
                    let (cl, gs, ap) =
                        self.active_world_mut().trigger_multicellular_transition_pulse();
                    self.export_status = Some(if gs > 0 || cl > 0 {
                        format!(
                            "Multicellular ETI: +{cl} Clonal clusters, +{gs} Germ-Soma organisms, {ap} apoptotic culls!"
                        )
                    } else {
                        "Need bioenergetic Holobionts to nucleate Clonal Multicellularity & Germ-Soma!".to_string()
                    });
                }
                if ui.button("MRO Reduce!").clicked() {
                    let (hy, mi, am) = self.active_world_mut().trigger_mro_reduction_pulse();
                    self.export_status = Some(if am > 0 {
                        format!(
                            "Cytosolic SUF active! {am} mitosome(s) shed completely -> Monocercomonoides MRO->0!"
                        )
                    } else if hy > 0 || mi > 0 {
                        format!(
                            "Anaerobic MRO step: +{hy} Hydrogenosomes (mtDNA=0, H2+ISC), +{mi} Mitosomes (ISC=1)!"
                        )
                    } else {
                        "Need MRO-bearing Holobionts, or raise Cytosolic SUF >= 0.65 to pass Mitosome ISC floor.".to_string()
                    });
                }
                if ui.button("Apico Reduce!").clicked() {
                    let (conv, shed, lost) =
                        self.active_world_mut().trigger_apicoplast_reduction_pulse();
                    self.export_status = Some(if lost > 0 {
                        format!(
                            "IPP rescue active! {lost} parasite(s) shed apicoplast completely (Cryptosporidium Apico->0)!"
                        )
                    } else if conv > 0 || shed > 0 {
                        format!(
                            "Holo->Apico: {conv} converted, {shed} shed Heme/FASII (retaining MEP/IPP + SUF)!"
                        )
                    } else {
                        "Need extant Holobionts, or raise Host IPP Rescue >= 0.65 to pass MEP/SUF floor.".to_string()
                    });
                }
                if ui.button("Klepto Steal!").clicked() {
                    let (st, hg, fx) = self.active_world_mut().trigger_klepto_steal_pulse();
                    self.export_status = Some(if fx > 0 {
                        format!(
                            "Tipping point crossed! {fx} kleptoplastic predator(s) fixed permanent plastids (+{hg} HGT)!"
                        )
                    } else if st > 0 || hg > 0 {
                        format!("Stole {st} kleptoplasts & recruited +{hg} shopping-bag HGT loci!")
                    } else {
                        "Need coexisting non-apicoplast PredatorParasite & Producer/Holobiont prey.".to_string()
                    });
                }
                if ui.button("Symbio 1°!").clicked() {
                    let fused = self.active_world_mut().trigger_endosymbiotic_pulse();
                    self.export_status = Some(if fused {
                        "Fused Producer + Syntroph -> 1° Holobiont (#H, Aerobic Mito)!".to_string()
                    } else {
                        "Need coexisting Producer & Syntroph for 1° fusion.".to_string()
                    });
                }
                if ui.button("Serial 2°/3°!").clicked() {
                    let tier_opt = self.active_world_mut().trigger_serial_endosymbiosis_pulse();
                    self.export_status = Some(match tier_opt {
                        Some(t) => format!(
                            "Engulfed eukaryotic Holobiont -> {t}° Complex Chimera (4 membranes + Nucleomorph)!"
                        ),
                        None => "Need an extant 1° or 2° Holobiont + host to trigger serial engulfment.".to_string(),
                    });
                }
                if ui.button("Export NHX+SVG").clicked() {
                    self.export_snapshot_files();
                }
                ui.add(
                    egui::Slider::new(&mut self.steps_per_second, 1.0..=30.0)
                        .text("ticks/s")
                        .logarithmic(false),
                );

                ui.separator();
                ui.label("Organism:");
                for (idx, org) in OrganismKind::ALL.iter().enumerate() {
                    if ui
                        .selectable_label(self.active_org_idx == idx, org.short_title())
                        .clicked()
                    {
                        self.active_org_idx = idx;
                    }
                }

                ui.separator();
                ui.label("View:");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialAndPhylogeny,
                    "Spatial + Phylogeny",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialGridAnd3D,
                    "Spatial + 3D",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::PhylogenyFull,
                    "Full Phylogram",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::Perspective3D,
                    "3D Structure",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::Blueprint2D,
                    "2D Blueprint",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::TriOrganism3D,
                    "All 3 Organisms",
                );

                ui.separator();
                ui.checkbox(&mut self.auto_rotate, "Auto-Orbit");
                ui.checkbox(&mut self.show_labels, "Labels");
            });
            ui.add_space(4.0);
        });
    }

    fn draw_left_controls(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("itme_left_controls")
            .default_width(345.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Simulation Mode");
                    let curr_mode = self.active_world().mode;
                    let mut selected_mode = curr_mode;
                    ui.horizontal(|ui| {
                        ui.selectable_value(
                            &mut selected_mode,
                            SimulationMode::GuidedCurriculum,
                            "13-Regime Tour",
                        );
                        ui.selectable_value(
                            &mut selected_mode,
                            SimulationMode::OrganicSandbox,
                            "Autonomous Sandbox",
                        );
                    });
                    if selected_mode != curr_mode {
                        for w in &mut self.worlds {
                            w.mode = selected_mode;
                        }
                    }

                    if selected_mode == SimulationMode::GuidedCurriculum {
                        let mut sps = self.active_world().steps_per_stage;
                        if ui
                            .add(egui::Slider::new(&mut sps, 12..=120).text("Ticks / Regime"))
                            .changed()
                        {
                            for w in &mut self.worlds {
                                w.steps_per_stage = sps;
                            }
                        }
                    } else if ui.button("Stochastic Phase Jump").clicked() {
                        self.active_world_mut().trigger_network_transition();
                    }

                    ui.separator();
                    ui.heading("13 Emergent Regimes");
                    let active_classified = self.active_world().classified_regime;
                    let mut clicked_regime = None;
                    for regime in RegimeId::ALL {
                        let is_curr = regime == active_classified;
                        let txt = format!("{}  {}", regime.code(), regime.name());
                        let rich = if is_curr {
                            egui::RichText::new(txt)
                                .strong()
                                .color(egui::Color32::from_rgb(105, 240, 185))
                        } else {
                            egui::RichText::new(txt)
                        };
                        if ui.selectable_label(is_curr, rich).clicked() {
                            clicked_regime = Some(regime);
                        }
                    }
                    if let Some(reg) = clicked_regime {
                        self.jump_all_or_active_to_regime(reg);
                    }

                    ui.separator();
                    ui.heading("Anaerobic MRO, Apicoplast & Symbiogenesis");
                    let theta_f = self.active_world().latest_metrics().critical_fidelity_theta;
                    {
                        let w = self.active_world_mut();
                        let mut changed = false;
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.rnr_thiyl_radical_activity,
                                    0.0..=1.00,
                                )
                                .text("k_rnr (Class III RNR 2'-OH->2'-H Radical)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.thymidylate_udg_repair_efficiency,
                                    0.0..=1.00,
                                )
                                .text("e_udg (ThyX U->T + UDG Deamination Repair)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.virosphere_arms_race_coupling,
                                    0.0..=1.00,
                                )
                                .text("v_sym (Virosphere Pol Transfer & Immunity)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.magnetite_ciss_spin_bias,
                                    0.0..=1.00,
                                )
                                .text("b_ciss (Magnetite Fe3O4 CISS Chiral Bias)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.minihelix_chiral_selectivity,
                                    0.0..=1.00,
                                )
                                .text("s_ac (Minihelix 3'-NCCA L-AA Steric Filter)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.communal_code_hgt_rate,
                                    0.0..=1.00,
                                )
                                .text("h_sgc (Communal Progenote SGC HGT Rate)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.wet_dry_cycle_amplitude,
                                    0.0..=1.00,
                                )
                                .text("A_wd (Hot-Spring Wet-Dry Cycles)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.hydrothermal_delta_ph,
                                    0.0..=4.50,
                                )
                                .text("Delta_pH (Alkaline Vent FeNiS Gradient)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.citrate_mg_chelation,
                                    0.0..=1.00,
                                )
                                .text("c_Mg (Citrate-Mg2+ Vesicle Protection)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.ribozyme_polymerase_fidelity,
                                    0.0..=1.00,
                                )
                                .text("f_ribo (Ribozyme Fidelity vs Eigen Cliff)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.biogeochemical_feedback_gain,
                                    0.0..=1.00,
                                )
                                .text("k_gaia (Planetary Biogeochem Coupling)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.planetary_stewardship_coupling,
                                    0.0..=1.00,
                                )
                                .text("s_gaia2 (Homeostatic Stewardship >=0.65)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.extractive_burn_intensity,
                                    0.0..=1.00,
                                )
                                .text("b_burn (Extractive Carbon Burn)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.synthetic_generation_rate,
                                    0.0..=1.00,
                                )
                                .text("g_synth (Synthetic Sampler Throughput)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.formal_verifier_strictness,
                                    0.0..=1.00,
                                )
                                .text("v_formal (Proof-Evaluator Strictness >=0.72)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.empirical_grounding_floor,
                                    0.0..=0.60,
                                )
                                .text("gamma_emp (Empirical Grounding Anchor >=0.12)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.archival_externalization_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_arch (Exosomatic Archive Drive)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.peer_verification_rigor,
                                    0.0..=1.00,
                                )
                                .text("r_peer (Peer-Verification Rigor >=0.70)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.unverified_publish_pressure,
                                    0.0..=1.00,
                                )
                                .text("p_fdr (Unverified Publish Pressure)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.symbolic_teaching_fidelity,
                                    0.0..=1.00,
                                )
                                .text("f_teach (Symbolic Teaching Fidelity >=0.88)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.inter_deme_connectivity,
                                    0.0..=1.00,
                                )
                                .text("c_deme (Deme Connectivity — 0.35 Opt.)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.cultural_innovation_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_innov (Cultural Innovation Drive)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.crop_domestication_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_crop (Gongylidia Domest. Drive)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.agro_pathogen_pressure,
                                    0.0..=1.00,
                                )
                                .text("p_escov (Escovopsis Blight Press.)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.defensive_symbiont_potency,
                                    0.0..=1.00,
                                )
                                .text("d_pseudo (Pseudonocardia Shield)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.queen_monogamy_index,
                                    0.0..=1.00,
                                )
                                .text("m_queen (Monogamy Window >=0.75)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.eusocial_caste_canalization_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_caste (JH Sterile Caste Canal.)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.worker_policing_rate,
                                    0.0..=1.00,
                                )
                                .text("p_police (Worker Kin-Policing)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.multicellular_adhesion_rate,
                                    0.0..=0.90,
                                )
                                .text("alpha_adh (Clonal Staying-Together)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.germ_soma_differentiation_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_gs (Germ-Soma regA Sequest.)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.apoptosis_policing_rate,
                                    0.0..=1.00,
                                )
                                .text("p_apop (Mito Apoptosis Policing)"),
                            )
                            .changed();
                        changed |= ui
                            .checkbox(
                                &mut w.controls.single_cell_bottleneck,
                                "Zygotic N=1 Bottleneck (r -> 1.0)",
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.ambient_oxygen_level, 0.0..=1.00)
                                    .text("O2_level (Hypoxia < 0.45)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.anaerobic_mro_reduction_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_mro (Mito->Hydro->MitoS)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.cytosolic_suf_hgt_rescue,
                                    0.0..=1.00,
                                )
                                .text("suf_rescue (Cytosolic SUF >=0.65)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.apicoplast_reduction_rate,
                                    0.0..=0.90,
                                )
                                .text("mu_apico (Holo->Apico Reduction)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut w.controls.host_metabolite_rescue,
                                    0.0..=1.00,
                                )
                                .text("ipp_rescue (Host IPP Scavenge)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.kleptoplasty_rate, 0.0..=0.90)
                                    .text("gamma_kp (Kleptoplast Theft)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.klepto_hgt_rate, 0.0..=0.90)
                                    .text("mu_hgt (Shopping-Bag HGT)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.endosymbiosis_rate, 0.0..=0.80)
                                    .text("gamma_1 (1° Endosymbiosis)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.cross_feeding_rate, 0.0..=0.75)
                                    .text("beta_cf (H2 Cross-Feed)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.bath_exergy_flux, -0.6..=1.5)
                                    .text("J_in (Exergy Flux)"),
                            )
                            .changed();

                        ui.separator();
                        ui.heading("Thermodynamic & Evolutionary Controls");
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.landauer_bit_cost, 0.0..=0.008)
                                    .text("k_B T ln2 (Cost)"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut w.controls.replication_fidelity, 0.40..=1.0)
                                    .text("F_0 (Base Fid.)"),
                            )
                            .changed();

                        let f_eff = w.controls.effective_fidelity();
                        ui.label(
                            egui::RichText::new(format!(
                                "F_eff = {:.3} vs Eigen theta_F = {:.3}",
                                f_eff, theta_f
                            ))
                            .small()
                            .color(if f_eff >= theta_f {
                                egui::Color32::from_rgb(100, 235, 165)
                            } else {
                                egui::Color32::from_rgb(255, 145, 105)
                            }),
                        );
                        changed |= ui
                            .checkbox(
                                &mut w.controls.artificial_selection,
                                "Meta-Engineered Objective (4.0)",
                            )
                            .changed();

                        if changed {
                            w.mode = SimulationMode::OrganicSandbox;
                        }
                    }
                });
            });
    }

    fn draw_right_inspector(&mut self, ctx: &egui::Context) {
        egui::SidePanel::right("itme_right_inspector")
            .default_width(370.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    let w = self.active_world();
                    let m = w.latest_metrics();
                    let reg = w.classified_regime;

                    ui.heading("Emergent Classified Regime");
                    ui.label(
                        egui::RichText::new(reg.full_title())
                            .size(16.0)
                            .strong()
                            .color(egui::Color32::from_rgb(255, 214, 102)),
                    );
                    ui.label(egui::RichText::new(reg.core_concept()).italics());
                    ui.monospace(format!("Emergent {}", m.diagnostic));

                    let phase_col = match m.phase {
                        OrderPhase::Adaptive => egui::Color32::from_rgb(95, 235, 165),
                        OrderPhase::Stasis => egui::Color32::from_rgb(135, 195, 255),
                        OrderPhase::Regressive => egui::Color32::from_rgb(255, 110, 110),
                    };
                    ui.colored_label(phase_col, format!("Phase: {}", m.phase.label()));

                    if let Some(ref status) = self.export_status {
                        ui.colored_label(egui::Color32::from_rgb(130, 245, 205), status);
                    }

                    ui.separator();
                    ui.heading("Anaerobic MRO, Apicoplast & Symbio Telemetry");
                    egui::Grid::new("itme_v12_grid")
                        .striped(true)
                        .num_columns(2)
                        .show(ui, |ui| {
                            ui.label("MROs (Aero / Hydro / MitoS / 0)");
                            ui.monospace(format!(
                                "{} / {} / {} / {} (Paths:{:.1})",
                                m.aerobic_mito_count,
                                m.hydrogenosome_count,
                                m.mitosome_count,
                                m.amiteukaryote_count,
                                m.mean_mro_pathways
                            ));
                            ui.end_row();

                            ui.label("MRO Transitions");
                            ui.monospace(format!(
                                "Mito->Hy:{} | Hy->Mi:{} | MRO->0:{}",
                                m.hydrogenosome_transition_events,
                                m.mitosome_transition_events,
                                m.mro_complete_loss_events
                            ));
                            ui.end_row();

                            ui.label("Apicoplast Parasites / Paths");
                            ui.monospace(format!(
                                "{} / {:.1} paths (Holo->Ap:{}, Ap->0:{})",
                                m.apicoplast_parasite_count,
                                m.mean_apicoplast_pathways,
                                m.apicoplast_transition_events,
                                m.apicoplast_loss_events
                            ));
                            ui.end_row();

                            ui.label("Guilds (P / Syn / Hol / X)");
                            ui.monospace(format!(
                                "{} / {} / {} / {}",
                                m.producer_count,
                                m.syntroph_count,
                                m.holobiont_count,
                                m.predator_count
                            ));
                            ui.end_row();

                            ui.label("Order O / O_pop");
                            ui.monospace(format!(
                                "{:.2} / {:.2}",
                                m.order_o, m.ecosystem_order_pop
                            ));
                            ui.end_row();
                        });

                    ui.separator();
                    ui.heading("Emergent Capabilities");
                    ui.horizontal_wrapped(|ui| {
                        flag_chip(ui, "Driven", m.flags.exergy_driven);
                        flag_chip(ui, "Replicate", m.flags.replication);
                        flag_chip(ui, "Select", m.flags.selection);
                        flag_chip(ui, "1°Holo", m.flags.endosymbiotic_fusion);
                        flag_chip(ui, "2°/3°Serial", m.flags.serial_endosymbiosis);
                        flag_chip(ui, "Hydrogenosome", m.flags.hydrogenosome_fermentation);
                        flag_chip(ui, "RNR_dNTP", m.flags.rnr_deoxyribose_synthesis);
                        flag_chip(ui, "ThyX_UDG_BER", m.flags.thymine_udg_deamination_repair);
                        flag_chip(ui, "ViroCoevolved", m.flags.virosphere_dna_coevolution);
                        flag_chip(ui, "Homochiral_DL", m.flags.homochiral_d_rna_l_aa);
                        flag_chip(ui, "MinihelixAcyl", m.flags.stereoselective_aminoacylation);
                        flag_chip(ui, "UniversalSGC", m.flags.error_minimized_universal_sgc);
                        flag_chip(ui, "RAF_Closure", m.flags.prebiotic_raf_closure);
                        flag_chip(ui, "RiboHypercycle", m.flags.ribozyme_hypercycle_active);
                        flag_chip(ui, "DarwinianLUCA", m.flags.darwinian_threshold_luca);
                        flag_chip(ui, "GOE_O2", m.flags.oxygenic_goe_regulation);
                        flag_chip(ui, "NutrientLoop", m.flags.syntrophic_nutrient_closure);
                        flag_chip(ui, "Gaia2Awakened", m.flags.planetary_homeostatic_awakening);
                        flag_chip(ui, "CompSymbiont", m.flags.computational_symbiont);
                        flag_chip(ui, "FormalProof", m.flags.formal_proof_verification);
                        flag_chip(ui, "CollapseGuard", m.flags.model_collapse_safeguard);
                        flag_chip(ui, "ExoOffload", m.flags.exosomatic_memory_offloading);
                        flag_chip(ui, "WrittenArch", m.flags.written_archival_storage);
                        flag_chip(ui, "PeerVerified", m.flags.institutional_peer_verification);
                        flag_chip(ui, "DualInherit", m.flags.dual_gene_culture_inheritance);
                        flag_chip(ui, "CultRatchet", m.flags.cumulative_cultural_ratchet);
                        flag_chip(ui, "SymbolicTeach", m.flags.symbolic_pedagogical_teaching);
                        flag_chip(ui, "AgriDomest", m.flags.agricultural_domestication);
                        flag_chip(ui, "ClonalCrop", m.flags.clonal_crop_monoculture);
                        flag_chip(ui, "ActinoShield", m.flags.defensive_microbiome_symbiosis);
                        flag_chip(ui, "Eusocial", m.flags.eusocial_colony);
                        flag_chip(ui, "SuperorgCastes", m.flags.sterile_worker_castes);
                        flag_chip(ui, "KinPolicing", m.flags.worker_kin_policing);
                        flag_chip(ui, "ClonalMC", m.flags.clonal_multicellularity);
                        flag_chip(ui, "GermSoma", m.flags.germ_soma_division);
                        flag_chip(ui, "Apoptosis", m.flags.apoptotic_cheater_policing);
                        flag_chip(ui, "Mitosome(ISC)", m.flags.mitosome_isc_retention);
                        flag_chip(ui, "MRO->0(SUF)", m.flags.amiteukaryote_mro_loss);
                        flag_chip(ui, "Apicoplast", m.flags.apicoplast_reduction);
                        flag_chip(ui, "Apico->0", m.flags.apicoplast_loss);
                    });

                    ui.separator();
                    ui.heading("Live Extended Newick (eNewick #H + MRO NHX)");
                    let nhx_str = w.export_newick(true);
                    let snippet = if nhx_str.chars().count() > 240 {
                        let head: String = nhx_str.chars().take(240).collect();
                        format!("{head}...")
                    } else {
                        nhx_str
                    };
                    ui.code(snippet);
                });
            });
    }

    fn draw_bottom_telemetry_plots(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("itme_bottom_plots")
            .min_height(150.0)
            .default_height(160.0)
            .resizable(true)
            .show(ctx, |ui| {
                let w = self.active_world();
                ui.horizontal(|ui| {
                    ui.strong(format!("ITME v12.0 Trajectories — {}", w.organism.title()));
                    ui.separator();
                    ui.colored_label(egui::Color32::from_rgb(80, 220, 255), "Order (O)");
                    ui.colored_label(egui::Color32::from_rgb(6, 214, 160), "Hydrogenosomes");
                    ui.colored_label(egui::Color32::from_rgb(255, 209, 102), "Mitosomes (ISC=1)");
                });

                let avail = ui.available_size();
                let (_resp, painter) = ui.allocate_painter(avail, egui::Sense::hover());
                let rect = painter.clip_rect().shrink(4.0);
                paint_framed_box(
                    &painter,
                    rect,
                    6.0,
                    egui::Color32::from_rgb(10, 14, 26),
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(38, 52, 84)),
                );

                if w.history.len() < 2 {
                    return;
                }

                let left_plot = egui::Rect::from_min_max(
                    rect.min + egui::vec2(8.0, 8.0),
                    egui::pos2(rect.center().x - 6.0, rect.max.y - 8.0),
                );
                let right_plot = egui::Rect::from_min_max(
                    egui::pos2(rect.center().x + 6.0, rect.min.y + 8.0),
                    rect.max - egui::vec2(8.0, 8.0),
                );

                let max_val = w
                    .history
                    .iter()
                    .map(|m| m.order_o.max(1.4))
                    .fold(1.4_f64, f64::max) as f32;
                draw_series(
                    &painter,
                    left_plot,
                    &w.history,
                    |m| m.order_o as f32,
                    0.0,
                    max_val,
                    egui::Color32::from_rgb(80, 220, 255),
                );
                painter.text(
                    left_plot.left_top() + egui::vec2(6.0, 4.0),
                    egui::Align2::LEFT_TOP,
                    format!("Order O [0 .. {:.2}]", max_val),
                    egui::FontId::monospace(10.5),
                    egui::Color32::from_rgb(175, 190, 220),
                );

                let max_right = 10.0_f32;
                draw_series(
                    &painter,
                    right_plot,
                    &w.history,
                    |m| m.hydrogenosome_count as f32,
                    0.0,
                    max_right,
                    egui::Color32::from_rgb(6, 214, 160),
                );
                draw_series(
                    &painter,
                    right_plot,
                    &w.history,
                    |m| m.mitosome_count as f32,
                    0.0,
                    max_right,
                    egui::Color32::from_rgb(255, 209, 102),
                );
                painter.text(
                    right_plot.left_top() + egui::vec2(6.0, 4.0),
                    egui::Align2::LEFT_TOP,
                    "Hydrogenosome (H2+ISC) & Mitosome (ISC-Only) Lineage Counts",
                    egui::FontId::monospace(10.5),
                    egui::Color32::from_rgb(175, 190, 220),
                );
            });
    }

    fn draw_central_viewport(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let avail = ui.available_size();
            let (response, painter) = ui.allocate_painter(avail, egui::Sense::click_and_drag());

            if response.dragged() {
                let d = response.drag_delta();
                self.camera.yaw += (d.x as f64) * 0.011;
                self.camera.pitch = (self.camera.pitch + (d.y as f64) * 0.011).clamp(-1.35, 1.35);
            }
            if response.hovered() {
                let scroll = ui.input(|i| i.smooth_scroll_delta.y);
                if scroll.abs() > 0.1 {
                    self.camera.zoom =
                        (self.camera.zoom * (1.0 + (scroll as f64) * 0.0014)).clamp(0.45, 2.8);
                }
            }

            let full_rect = response.rect;
            paint_framed_box(
                &painter,
                full_rect,
                8.0,
                egui::Color32::from_rgb(9, 13, 24),
                egui::Stroke::new(1.2, egui::Color32::from_rgb(36, 52, 88)),
            );

            match self.viewport_layout {
                ViewportLayout::SpatialAndPhylogeny => {
                    let w = self.active_world();
                    let mid_x = full_rect.left() + full_rect.width() * 0.46;
                    let left_r = egui::Rect::from_min_max(
                        full_rect.min + egui::vec2(6.0, 6.0),
                        egui::pos2(mid_x - 4.0, full_rect.max.y - 6.0),
                    );
                    let right_r = egui::Rect::from_min_max(
                        egui::pos2(mid_x + 4.0, full_rect.min.y + 6.0),
                        full_rect.max - egui::vec2(6.0, 6.0),
                    );
                    paint_framed_box(
                        &painter,
                        left_r,
                        6.0,
                        egui::Color32::from_rgb(11, 17, 30),
                        egui::Stroke::new(1.0, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_framed_box(
                        &painter,
                        right_r,
                        6.0,
                        egui::Color32::from_rgb(11, 17, 30),
                        egui::Stroke::new(1.0, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_spatial_grid(&painter, left_r.shrink(4.0), w, self.show_labels);
                    paint_phylogram_tree(&painter, right_r.shrink(4.0), w, self.show_labels);
                }
                ViewportLayout::SpatialGridAnd3D => {
                    let w = self.active_world();
                    let scene = w.current_scene();
                    let mid_x = full_rect.left() + full_rect.width() * 0.48;
                    let left_r = egui::Rect::from_min_max(
                        full_rect.min + egui::vec2(6.0, 6.0),
                        egui::pos2(mid_x - 4.0, full_rect.max.y - 6.0),
                    );
                    let right_r = egui::Rect::from_min_max(
                        egui::pos2(mid_x + 4.0, full_rect.min.y + 6.0),
                        full_rect.max - egui::vec2(6.0, 6.0),
                    );
                    paint_framed_box(
                        &painter,
                        left_r,
                        6.0,
                        egui::Color32::from_rgb(11, 17, 30),
                        egui::Stroke::new(1.0, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_spatial_grid(&painter, left_r.shrink(4.0), w, self.show_labels);
                    paint_scene_3d(&painter, right_r, &scene, &self.camera, self.show_labels);
                }
                ViewportLayout::PhylogenyFull => {
                    let w = self.active_world();
                    paint_phylogram_tree(&painter, full_rect.shrink(8.0), w, self.show_labels);
                }
                ViewportLayout::Perspective3D => {
                    let scene = self.active_world().current_scene();
                    paint_scene_3d(
                        &painter,
                        full_rect.shrink(6.0),
                        &scene,
                        &self.camera,
                        self.show_labels,
                    );
                }
                ViewportLayout::Blueprint2D => {
                    let scene = self.active_world().current_scene();
                    paint_scene_2d(&painter, full_rect.shrink(6.0), &scene, self.show_labels);
                }
                ViewportLayout::TriOrganism3D => {
                    let third = full_rect.width() / 3.0;
                    for i in 0..3 {
                        let sub = egui::Rect::from_min_max(
                            egui::pos2(
                                full_rect.left() + third * (i as f32) + 4.0,
                                full_rect.top() + 4.0,
                            ),
                            egui::pos2(
                                full_rect.left() + third * ((i + 1) as f32) - 4.0,
                                full_rect.bottom() - 4.0,
                            ),
                        );
                        paint_framed_box(
                            &painter,
                            sub,
                            6.0,
                            egui::Color32::from_rgb(11, 16, 30),
                            egui::Stroke::new(1.0, egui::Color32::from_rgb(38, 54, 90)),
                        );
                        let scene = self.worlds[i].current_scene();
                        paint_scene_3d(
                            &painter,
                            sub.shrink(4.0),
                            &scene,
                            &self.camera,
                            self.show_labels,
                        );
                    }
                }
            }
        });
    }
}

fn guild_color(guild: TrophicGuild) -> egui::Color32 {
    match guild {
        TrophicGuild::Producer => egui::Color32::from_rgb(76, 201, 240),
        TrophicGuild::Syntroph => egui::Color32::from_rgb(80, 227, 137),
        TrophicGuild::PredatorParasite => egui::Color32::from_rgb(255, 92, 138),
        TrophicGuild::Holobiont => egui::Color32::from_rgb(199, 125, 255),
    }
}

fn paint_phylogram_tree(
    painter: &egui::Painter,
    rect: egui::Rect,
    world: &SimulationWorld,
    show_labels: bool,
) {
    let layout = world.current_phylogram_layout();
    let p_sum = world.phylo_summary();

    painter.text(
        rect.left_top() + egui::vec2(8.0, 4.0),
        egui::Align2::LEFT_TOP,
        format!(
            "Network — MRCA N{} | Mito->Hy:{} | Hy->MitoS:{} | MRO->0:{} | #H:{}",
            p_sum.mrca_id,
            p_sum.hydrogenosome_transition_events,
            p_sum.mitosome_transition_events,
            p_sum.mro_complete_loss_events,
            p_sum.endosymbiotic_fusions
        ),
        egui::FontId::proportional(13.0),
        egui::Color32::from_rgb(230, 240, 255),
    );
    painter.text(
        rect.left_top() + egui::vec2(8.0, 22.0),
        egui::Align2::LEFT_TOP,
        "Teal=Hydrogenosome(H2+ISC) | Gold=Mitosome(ISC=1) | Cyan=MRO->0(SUF) | Violet=AeroMito",
        egui::FontId::monospace(9.5),
        egui::Color32::from_rgb(142, 168, 212),
    );

    let plot = egui::Rect::from_min_max(
        rect.left_top() + egui::vec2(18.0, 44.0),
        rect.right_bottom() - egui::vec2(if show_labels { 118.0 } else { 18.0 }, 24.0),
    );
    if plot.width() < 20.0 || plot.height() < 20.0 {
        return;
    }

    let to_pt = |xn: f64, yn: f64| -> egui::Pos2 {
        egui::pos2(
            plot.left() + (xn as f32) * plot.width(),
            plot.top() + (yn as f32) * plot.height(),
        )
    };

    for k in 0..=4 {
        let f = (k as f32) / 4.0;
        let x = plot.left() + f * plot.width();
        painter.line_segment(
            [egui::pos2(x, plot.top()), egui::pos2(x, plot.bottom())],
            egui::Stroke::new(0.8, egui::Color32::from_rgba_unmultiplied(55, 76, 120, 80)),
        );
        let t_val =
            layout.min_step + (((layout.max_step - layout.min_step) as f32) * f).round() as usize;
        painter.text(
            egui::pos2(x, plot.bottom() + 4.0),
            egui::Align2::CENTER_TOP,
            format!("t={t_val}"),
            egui::FontId::monospace(9.5),
            egui::Color32::from_rgb(135, 155, 195),
        );
    }

    for br in &layout.branches {
        let p_parent = to_pt(br.parent_x, br.parent_y);
        let p_corner = to_pt(br.parent_x, br.child_y);
        let p_child = to_pt(br.child_x, br.child_y);
        let col = match br.child_mro_class {
            MroClass::Hydrogenosome => egui::Color32::from_rgb(6, 214, 160),
            MroClass::Mitosome => egui::Color32::from_rgb(255, 209, 102),
            MroClass::Amiteukaryote => egui::Color32::from_rgb(56, 189, 248),
            _ if br.child_guild == TrophicGuild::Holobiont => match br.child_endo_tier {
                3..=u8::MAX => egui::Color32::from_rgb(247, 37, 133),
                2 => egui::Color32::from_rgb(255, 159, 28),
                _ => guild_color(br.child_guild),
            },
            _ if br.child_has_apicoplast => egui::Color32::from_rgb(224, 102, 255),
            _ => guild_color(br.child_guild),
        };
        let sw = if br.child_mro_class.is_reduced_anaerobic_mro()
            || br.child_endo_tier >= 2
            || br.child_has_apicoplast
        {
            2.8
        } else if br.child_guild == TrophicGuild::Holobiont {
            2.4
        } else {
            1.9
        };

        painter.line_segment([p_parent, p_corner], egui::Stroke::new(sw, col));
        painter.line_segment([p_corner, p_child], egui::Stroke::new(sw, col));

        let mid = egui::pos2((p_corner.x + p_child.x) * 0.5, p_child.y);
        if br.is_hydrogenosome_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(6, 214, 160));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8, egui::Color32::from_rgb(255, 209, 102)),
            );
        } else if br.is_mitosome_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(255, 209, 102));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8, egui::Color32::from_rgb(255, 92, 138)),
            );
        } else if br.is_mro_complete_loss {
            painter.circle_filled(mid, 5.6, egui::Color32::from_rgb(56, 189, 248));
            painter.circle_stroke(
                mid,
                5.6,
                egui::Stroke::new(1.7, egui::Color32::from_rgb(255, 209, 102)),
            );
        } else if br.is_apicoplast_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(255, 92, 138));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8, egui::Color32::from_rgb(199, 125, 255)),
            );
        }
    }

    for ret in &layout.reticulations {
        let p1 = to_pt(ret.from_x, ret.from_y);
        let p2 = to_pt(ret.to_x, ret.to_y);
        if ret.is_endosymbiosis {
            let arc_col = if ret.is_klepto_tipping {
                egui::Color32::from_rgb(163, 255, 87)
            } else {
                match ret.endo_tier {
                    3..=u8::MAX => egui::Color32::from_rgb(247, 37, 133),
                    2 => egui::Color32::from_rgb(255, 159, 28),
                    _ => egui::Color32::from_rgb(80, 227, 137),
                }
            };
            let segs = 12usize;
            let span = (p2.x - p1.x).abs().max(14.0);
            let c1 = egui::pos2(p1.x + span * 0.55, p1.y);
            let c2 = egui::pos2(p2.x - span * 0.35, p2.y);
            let mut prev = p1;
            for s in 1..=segs {
                let t = (s as f32) / (segs as f32);
                let u = 1.0 - t;
                let bx = u * u * u * p1.x
                    + 3.0 * u * u * t * c1.x
                    + 3.0 * u * t * t * c2.x
                    + t * t * t * p2.x;
                let by = u * u * u * p1.y
                    + 3.0 * u * u * t * c1.y
                    + 3.0 * u * t * t * c2.y
                    + t * t * t * p2.y;
                let next_p = egui::pos2(bx, by);
                if s % 2 == 1 {
                    painter.line_segment([prev, next_p], egui::Stroke::new(2.2, arc_col));
                }
                prev = next_p;
            }
            painter.circle_filled(p2, 5.8, arc_col);
        }
    }

    let extant_count = layout.nodes.iter().filter(|n| n.extant).count();
    let label_stride = if extant_count > 20 { 2 } else { 1 };
    let mut ext_seen = 0usize;
    for nd in &layout.nodes {
        let pt = to_pt(nd.x_norm, nd.y_norm);
        if nd.extant {
            let col = match nd.mro_class {
                MroClass::Hydrogenosome => egui::Color32::from_rgb(6, 214, 160),
                MroClass::Mitosome => egui::Color32::from_rgb(255, 209, 102),
                MroClass::Amiteukaryote => egui::Color32::from_rgb(56, 189, 248),
                _ if nd.has_apicoplast => egui::Color32::from_rgb(224, 102, 255),
                _ => guild_color(nd.guild),
            };
            painter.circle_filled(pt, 4.2, col);
            if show_labels && ext_seen % label_stride == 0 {
                painter.text(
                    pt + egui::vec2(6.0, 0.0),
                    egui::Align2::LEFT_CENTER,
                    &nd.label,
                    egui::FontId::monospace(9.5),
                    col,
                );
            }
            ext_seen += 1;
        } else {
            painter.circle_filled(pt, 2.2, egui::Color32::from_rgb(160, 182, 225));
        }
    }
}

fn flag_chip(ui: &mut egui::Ui, text: &str, active: bool) {
    let col = if active {
        egui::Color32::from_rgb(80, 225, 160)
    } else {
        egui::Color32::from_rgb(100, 112, 138)
    };
    let prefix = if active { "[on] " } else { "[--] " };
    ui.colored_label(col, format!("{}{}", prefix, text));
}

fn paint_framed_box(
    painter: &egui::Painter,
    rect: egui::Rect,
    rounding: f32,
    fill: egui::Color32,
    stroke: egui::Stroke,
) {
    painter.rect_filled(rect, rounding, stroke.color);
    painter.rect_filled(rect.shrink(stroke.width), rounding, fill);
}

fn paint_spatial_grid(
    painter: &egui::Painter,
    rect: egui::Rect,
    world: &SimulationWorld,
    show_labels: bool,
) {
    let m = world.latest_metrics();
    painter.text(
        rect.left_top() + egui::vec2(8.0, 4.0),
        egui::Align2::LEFT_TOP,
        format!(
            "Spatial Grid — AeroMito:{} Hydro:{} MitoS:{} MRO0:{}",
            m.aerobic_mito_count, m.hydrogenosome_count, m.mitosome_count, m.amiteukaryote_count
        ),
        egui::FontId::proportional(13.0),
        egui::Color32::from_rgb(230, 240, 255),
    );
    painter.text(
        rect.left_top() + egui::vec2(8.0, 22.0),
        egui::Align2::LEFT_TOP,
        "TealRing=Hydrogenosome | GoldCore=Mitosome(ISC) | CyanRing=Amiteukaryote(SUF)",
        egui::FontId::monospace(9.5),
        egui::Color32::from_rgb(142, 168, 212),
    );

    let grid_area = egui::Rect::from_min_max(
        rect.left_top() + egui::vec2(10.0, 42.0),
        rect.right_bottom() - egui::vec2(10.0, 14.0),
    );
    let side = grid_area.width().min(grid_area.height());
    let origin = egui::pos2(
        grid_area.center().x - side * 0.5,
        grid_area.center().y - side * 0.5,
    );
    let cw = side / (world.grid.width as f32);
    let ch = side / (world.grid.height as f32);

    for gy in 0..world.grid.height {
        for gx in 0..world.grid.width {
            let e_norm = (world.grid.local_exergy(gx, gy) / 1.65).clamp(0.0, 1.0) as f32;
            let b_norm = (world.grid.local_byproduct(gx, gy) / 1.10).clamp(0.0, 1.0) as f32;
            let n_norm = (world.grid.local_niche(gx, gy) / 1.40).clamp(0.0, 1.0) as f32;

            let r =
                (14.0 + 34.0 * e_norm + 150.0 * n_norm + 78.0 * b_norm).clamp(0.0, 255.0) as u8;
            let g =
                (20.0 + 128.0 * e_norm + 105.0 * n_norm + 45.0 * b_norm).clamp(0.0, 255.0) as u8;
            let b =
                (36.0 + 115.0 * e_norm + 20.0 * n_norm + 135.0 * b_norm).clamp(0.0, 255.0) as u8;

            let c_rect = egui::Rect::from_min_size(
                egui::pos2(origin.x + (gx as f32) * cw, origin.y + (gy as f32) * ch),
                egui::vec2((cw - 1.2).max(1.0), (ch - 1.2).max(1.0)),
            );
            painter.rect_filled(c_rect, 2.0, egui::Color32::from_rgb(r, g, b));
        }
    }

    let best_id = world.best_unit().map(|u| u.id);
    for (idx, unit) in world.population.iter().enumerate() {
        let jitter_x = (((idx % 3) as f32) - 1.0) * cw * 0.14;
        let jitter_y = ((((idx / 3) % 3) as f32) - 1.0) * ch * 0.14;
        let center = egui::pos2(
            origin.x + ((unit.grid_x as f32) + 0.5) * cw + jitter_x,
            origin.y + ((unit.grid_y as f32) + 0.5) * ch + jitter_y,
        );
        let [sr, sg, sb] = SPECIES_PALETTE[unit.species_id % SPECIES_PALETTE.len()];
        let is_champ = best_id == Some(unit.id);
        let rad = if is_champ {
            (cw * 0.35).clamp(4.2, 9.8)
        } else {
            (cw * 0.26).clamp(3.2, 7.4)
        };
        let fill = egui::Color32::from_rgb(sr, sg, sb);

        match unit.guild {
            TrophicGuild::Producer => {
                painter.circle_filled(center, rad, fill);
            }
            TrophicGuild::Syntroph => {
                let pts = vec![
                    egui::pos2(center.x, center.y - rad * 1.15),
                    egui::pos2(center.x + rad * 1.15, center.y),
                    egui::pos2(center.x, center.y + rad * 1.15),
                    egui::pos2(center.x - rad * 1.15, center.y),
                ];
                painter.add(egui::Shape::convex_polygon(
                    pts,
                    fill,
                    egui::Stroke::new(1.3, egui::Color32::from_rgb(80, 250, 170)),
                ));
            }
            TrophicGuild::Holobiont => {
                let host_col = match unit.mro_class {
                    MroClass::Hydrogenosome => egui::Color32::from_rgb(6, 214, 160),
                    MroClass::Mitosome => egui::Color32::from_rgb(255, 209, 102),
                    MroClass::Amiteukaryote => egui::Color32::from_rgb(56, 189, 248),
                    _ => match unit.endo_tier {
                        3..=u8::MAX => egui::Color32::from_rgb(247, 37, 133),
                        2 => egui::Color32::from_rgb(255, 159, 28),
                        _ => egui::Color32::from_rgb(199, 125, 255),
                    },
                };
                let mc_scale = match unit.multicell_state {
                    MulticellularState::GermSomaMulticellular => 1.55,
                    MulticellularState::UndifferentiatedClonal => 1.38,
                    _ => 1.25,
                };
                painter.circle_filled(center, rad * mc_scale, host_col);
                if unit.genome_substrate_stage != GenomeSubstrateStage::CoevolvedCellularVirosphere {
                    painter.circle_stroke(
                        center,
                        rad * 2.92,
                        egui::Stroke::new(2.2, egui::Color32::from_rgb(34, 211, 238)),
                    );
                } else if unit.chirality_code_stage != ChiralityCodeStage::UniversalErrorMinimizedSgc {
                    painter.circle_stroke(
                        center,
                        rad * 2.75,
                        egui::Stroke::new(2.2, egui::Color32::from_rgb(192, 132, 252)),
                    );
                } else if unit.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage {
                    painter.circle_stroke(
                        center,
                        rad * 2.45,
                        egui::Stroke::new(2.1, egui::Color32::from_rgb(251, 191, 36)),
                    );
                } else if unit.planetary_epoch == PlanetaryState::HomeostaticAwakening {
                    painter.circle_stroke(
                        center,
                        rad * 3.04,
                        egui::Stroke::new(2.6, egui::Color32::from_rgb(52, 211, 153)),
                    );
                } else if unit.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis {
                    painter.circle_stroke(
                        center,
                        rad * 2.82,
                        egui::Stroke::new(2.6, egui::Color32::from_rgb(139, 92, 246)),
                    );
                } else if unit.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                    painter.circle_stroke(
                        center,
                        rad * 2.58,
                        egui::Stroke::new(2.5, egui::Color32::from_rgb(20, 184, 166)),
                    );
                } else if unit.cultural_state == CulturalState::SymbolicCumulativeCulture {
                    painter.circle_stroke(
                        center,
                        rad * 2.34,
                        egui::Stroke::new(2.4, egui::Color32::from_rgb(250, 204, 21)),
                    );
                } else if unit.agriculture_state == AgricultureState::LeafCutterAgroEcosystem {
                    painter.circle_stroke(
                        center,
                        rad * 2.12,
                        egui::Stroke::new(2.4, egui::Color32::from_rgb(163, 230, 53)),
                    );
                } else if unit.eusocial_state == EusocialState::ObligateSuperorganism {
                    painter.circle_stroke(
                        center,
                        rad * 1.88,
                        egui::Stroke::new(2.2, egui::Color32::from_rgb(245, 158, 11)),
                    );
                } else if unit.multicell_state == MulticellularState::GermSomaMulticellular {
                    painter.circle_stroke(
                        center,
                        rad * 1.62,
                        egui::Stroke::new(1.8, egui::Color32::from_rgb(56, 189, 248)),
                    );
                }
                if unit.mro_class != MroClass::Amiteukaryote {
                    let core_col = if unit.mro_class == MroClass::Mitosome {
                        egui::Color32::from_rgb(255, 92, 138)
                    } else {
                        egui::Color32::from_rgb(255, 209, 102)
                    };
                    painter.circle_filled(center, rad * 0.44, core_col);
                }
            }
            TrophicGuild::PredatorParasite => {
                let pts = vec![
                    egui::pos2(center.x, center.y - rad * 1.25),
                    egui::pos2(center.x + rad * 1.15, center.y + rad * 0.95),
                    egui::pos2(center.x - rad * 1.15, center.y + rad * 0.95),
                ];
                let rim = if unit.has_apicoplast() {
                    egui::Color32::from_rgb(199, 125, 255)
                } else if unit.has_active_kleptoplast() {
                    egui::Color32::from_rgb(163, 255, 87)
                } else {
                    egui::Color32::from_rgb(255, 215, 120)
                };
                painter.add(egui::Shape::convex_polygon(
                    pts,
                    egui::Color32::from_rgb(255, 86, 122),
                    egui::Stroke::new(1.5, rim),
                ));
                if unit.has_apicoplast() {
                    painter.circle_filled(
                        center + egui::vec2(0.0, rad * 0.20),
                        rad * 0.44,
                        egui::Color32::from_rgb(255, 209, 102),
                    );
                }
            }
        }

        if show_labels && is_champ {
            painter.text(
                center + egui::vec2(0.0, -rad - 3.5),
                egui::Align2::CENTER_BOTTOM,
                format!("{} {}", unit.guild.short_code(), unit.mro_class.short_code()),
                egui::FontId::monospace(9.5),
                egui::Color32::WHITE,
            );
        }
    }
}

fn draw_series<F>(
    painter: &egui::Painter,
    rect: egui::Rect,
    data: &[crate::model::ItmeMetrics],
    extract: F,
    y_min: f32,
    y_max: f32,
    color: egui::Color32,
) where
    F: Fn(&crate::model::ItmeMetrics) -> f32,
{
    if data.len() < 2 {
        return;
    }
    let span = (y_max - y_min).max(1e-5);
    let n = (data.len() - 1) as f32;
    let mut prev = None;
    for (i, item) in data.iter().enumerate() {
        let xf = rect.left() + rect.width() * ((i as f32) / n);
        let norm_y = ((extract(item) - y_min) / span).clamp(0.0, 1.0);
        let yf = rect.bottom() - rect.height() * norm_y;
        let pt = egui::pos2(xf, yf);
        if let Some(p0) = prev {
            painter.line_segment([p0, pt], egui::Stroke::new(1.6, color));
        }
        prev = Some(pt);
    }
}

fn paint_scene_3d(
    painter: &egui::Painter,
    rect: egui::Rect,
    scene: &StructureScene,
    camera: &Camera3D,
    show_labels: bool,
) {
    let w = rect.width() as f64;
    let h = rect.height() as f64;

    let proj: Vec<(egui::Pos2, f64, f64)> = scene
        .nodes_3d
        .iter()
        .map(|n| {
            let (sx, sy, depth, persp) = camera.project(n.pos, w, h);
            (
                egui::pos2(rect.left() + sx as f32, rect.top() + sy as f32),
                depth,
                persp,
            )
        })
        .collect();

    for e in &scene.edges_3d {
        if e.a >= proj.len() || e.b >= proj.len() {
            continue;
        }
        let (p1, _, _) = proj[e.a];
        let (p2, _, _) = proj[e.b];
        if !rect.contains(p1) && !rect.contains(p2) {
            continue;
        }
        let [r, g, b] = e.color;
        let col = egui::Color32::from_rgba_unmultiplied(r, g, b, if e.dashed { 145 } else { 215 });
        painter.line_segment([p1, p2], egui::Stroke::new(e.width as f32, col));
    }

    let mut order: Vec<usize> = (0..scene.nodes_3d.len()).collect();
    order.sort_by(|&a, &b| proj[b].1.total_cmp(&proj[a].1));

    for idx in order {
        let node = &scene.nodes_3d[idx];
        let (center, _, persp) = proj[idx];
        if !rect.contains(center) {
            continue;
        }
        let rad = ((node.radius * persp) as f32).clamp(2.5, 22.0);
        let [r, g, b] = node.color;
        painter.circle_filled(center, rad, egui::Color32::from_rgb(r, g, b));
        if show_labels && !node.label.is_empty() && rad >= 4.5 {
            painter.text(
                center + egui::vec2(rad + 3.0, -2.0),
                egui::Align2::LEFT_CENTER,
                &node.label,
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(222, 232, 252),
            );
        }
    }
}

fn paint_scene_2d(
    painter: &egui::Painter,
    rect: egui::Rect,
    scene: &StructureScene,
    show_labels: bool,
) {
    let plot_r = rect.shrink2(egui::vec2(20.0, 34.0));
    let scale = plot_r.width().min(plot_r.height()) * 0.46;
    let center = plot_r.center();

    let to_screen = |nx: f64, ny: f64| -> egui::Pos2 {
        egui::pos2(center.x + (nx as f32) * scale, center.y - (ny as f32) * scale)
    };

    for e in &scene.edges_2d {
        if e.a < scene.nodes_2d.len() && e.b < scene.nodes_2d.len() {
            let p1 = to_screen(scene.nodes_2d[e.a].x, scene.nodes_2d[e.a].y);
            let p2 = to_screen(scene.nodes_2d[e.b].x, scene.nodes_2d[e.b].y);
            let [r, g, b] = e.color;
            painter.line_segment(
                [p1, p2],
                egui::Stroke::new(
                    e.width as f32,
                    egui::Color32::from_rgba_unmultiplied(r, g, b, 195),
                ),
            );
        }
    }

    for n in &scene.nodes_2d {
        let p = to_screen(n.x, n.y);
        let rad = ((n.radius as f32) * scale).clamp(4.0, 12.0);
        let [r, g, b] = n.color;
        painter.circle_filled(p, rad, egui::Color32::from_rgb(r, g, b));
        if show_labels && !n.label.is_empty() {
            painter.text(
                p + egui::vec2(0.0, -rad - 6.0),
                egui::Align2::CENTER_BOTTOM,
                &n.label,
                egui::FontId::monospace(9.5),
                egui::Color32::from_rgb(210, 224, 248),
            );
        }
    }
}

pub fn run_gui(seed: u64) -> Result<(), String> {
    let viewport = egui::ViewportBuilder::default()
        .with_title(
            "ITME v12.0 (v1.2.0) — Hydrogenosomes, Mitosomes (ISC), Apicoplasts & Symbiogenesis",
        )
        .with_inner_size(egui::vec2(1600.0, 940.0))
        .with_min_inner_size(egui::vec2(980.0, 640.0))
        .with_resizable(true);

    let options = eframe::NativeOptions {
        viewport,
        ..Default::default()
    };

    eframe::run_native(
        "ITME v12.0 Native Rust Simulation",
        options,
        Box::new(move |cc| Ok(Box::new(ItmeGuiApp::new(cc, seed)))),
    )
    .map_err(|e| format!("Failed to launch native GUI: {e}"))
}
