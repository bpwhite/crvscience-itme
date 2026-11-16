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

//! Native desktop GUI mode built with `eframe` / `egui` (`0.29`) — ITME v23.0.
//!
//! Starts at an isolated 3-Domain Main Menu launcher (`RNA/DNA/Protein`, `Hominid Populations`,
//! and `Artificial Intelligences`), where launching any substrate opens a single dedicated
//! `SimulationWorld` with domain-specific regime titles, grouped sliders, and telemetry.

use crate::engine::{SimulationMode, SimulationWorld};
use crate::geometry::{Camera3D, StructureScene};
use crate::model::{
    AbiogenesisStage, AgricultureState, ChiralityCodeStage, CulturalState, EusocialState,
    ExosomaticState, GenomeSubstrateStage, GrandTrajectoryEpochSummary, MachineSynthesisState,
    MroClass, MulticellularState, OrderPhase, PlanetaryState, TrophicGuild,
};
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
    GrandTrajectory,
}

struct MainMenuState {
    mode: SimulationMode,
    steps_per_stage: usize,
    initial_regime: RegimeId,
    seed_text: String,
    seed: u64,
}

impl MainMenuState {
    fn new(seed: u64) -> Self {
        Self {
            mode: SimulationMode::GuidedCurriculum,
            steps_per_stage: 32,
            initial_regime: RegimeId::R00Static,
            seed_text: seed.to_string(),
            seed,
        }
    }

    fn parsed_seed(&self) -> u64 {
        self.seed_text
            .trim()
            .parse::<u64>()
            .unwrap_or(self.seed)
    }
}

struct DomainSession {
    organism: OrganismKind,
    world: SimulationWorld,
    viewport_layout: ViewportLayout,
    camera: Camera3D,
    auto_rotate: bool,
    show_labels: bool,
    paused: bool,
    steps_per_second: f32,
    last_step_instant: Instant,
    export_status: Option<String>,
    grand_epochs: Vec<GrandTrajectoryEpochSummary>,
    wants_main_menu: bool,
}

enum AppScreen {
    MainMenu(MainMenuState),
    Simulation(Box<DomainSession>),
}

pub struct ItmeGuiApp {
    screen: AppScreen,
    default_seed: u64,
}

impl ItmeGuiApp {
    pub fn new(cc: &eframe::CreationContext<'_>, seed: u64) -> Self {
        Self::new_with_ctx(&cc.egui_ctx, seed)
    }

    fn new_with_ctx(ctx: &egui::Context, seed: u64) -> Self {
        let mut visuals = egui::Visuals::dark();
        visuals.panel_fill = egui::Color32::from_rgb(13, 18, 32);
        visuals.window_fill = egui::Color32::from_rgb(16, 22, 38);
        visuals.extreme_bg_color = egui::Color32::from_rgb(9, 13, 24);
        ctx.set_visuals(visuals);

        Self {
            screen: AppScreen::MainMenu(MainMenuState::new(seed)),
            default_seed: seed,
        }
    }

    fn draw_panels(&mut self, ctx: &egui::Context) {
        let mut next_screen: Option<AppScreen> = None;
        match &mut self.screen {
            AppScreen::MainMenu(menu) => {
                if let Some(session) = Self::draw_main_menu(ctx, menu) {
                    next_screen = Some(AppScreen::Simulation(Box::new(session)));
                }
            }
            AppScreen::Simulation(session) => {
                session.draw_panels(ctx);
                if session.wants_main_menu {
                    let seed = self.default_seed;
                    next_screen = Some(AppScreen::MainMenu(MainMenuState::new(seed)));
                }
            }
        }
        if let Some(scr) = next_screen {
            self.screen = scr;
        }
    }

    fn domain_card_meta(org: OrganismKind) -> (&'static str, &'static str, egui::Color32, &'static [&'static str], &'static str) {
        match org {
            OrganismKind::RnaDnaProtein => (
                "Domain 1 • Molecular Evolution",
                "Abiogenesis, Coding, & Organelle Symbiogenesis",
                egui::Color32::from_rgb(56, 189, 248),
                &[
                    "Hadean FeNiS delta-pH & wet-dry protocell closure",
                    "Fe3O4 CISS homochiral D-RNA / L-AA & SGC emergence",
                    "Class III RNR RNA -> T-DNA + UDG repair virosphere",
                    "1°/2°/3° endosymbiosis, Hydrogenosomes, Mitosomes (ISC), MRO->0",
                    "Clonal multicellularity, germ-soma sequestration & apoptosis",
                ],
                "Simulates prebiotic chemistry, Eigen's fidelity cliff, endosymbiont gene transfer, and reductive organelle trajectories.",
            ),
            OrganismKind::HominidPopulation => (
                "Domain 2 • Hominid Populations",
                "Biocultural Niche Construction, Pedagogy, & Written Archives",
                egui::Color32::from_rgb(251, 191, 36),
                &[
                    "Base-camp hearth provisioning & foraging band fissioning",
                    "Inter-band seasonal watershed & lithic guild competition",
                    "Cumulative symbolic teaching & cross-deme recombination",
                    "Agro-ectosymbiosis, cultivar domestication & blight shields",
                    "External written archives & peer-verified replication audits",
                ],
                "Simulates kin-structured bands, behavioral plasticity, cumulative cultural ratchets, subsistence mutualisms, and archival institutions.",
            ),
            OrganismKind::ArtificialIntelligence => (
                "Domain 3 • Artificial Intelligences",
                "Compute Allocation, Recursive Distillation, & Proof Verifiers",
                egui::Color32::from_rgb(167, 139, 250),
                &[
                    "Cluster load-balancing & stateless container checkpoint forking",
                    "Shared-KV / compute-quota contention & online test-time RL",
                    "Adversarial population training (PBT) & teacher distillation",
                    "Composable tool enclaves & deterministic verifier loops",
                    "Formal proof-verifier strictness & anti-collapse empiric guard",
                ],
                "Simulates accelerator compute budgets, synthetic tail-diversity decay (Shumailov cliff), teacher distillation, and formal verifiers.",
            ),
        }
    }

    fn draw_main_menu(ctx: &egui::Context, menu: &mut MainMenuState) -> Option<DomainSession> {
        let mut launched: Option<DomainSession> = None;

        egui::TopBottomPanel::top("itme_menu_top").show(ctx, |ui| {
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.heading(
                    egui::RichText::new("ITME v23.0 — Information-Theoretic & Thermodynamic Evolutionary Simulator")
                        .size(20.0)
                        .strong()
                        .color(egui::Color32::from_rgb(130, 245, 205)),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        egui::RichText::new("CRV Science Network • GPL-3.0")
                            .small()
                            .color(egui::Color32::from_rgb(148, 163, 184)),
                    );
                });
            });
            ui.label(
                egui::RichText::new(
                    "Select an isolated evolutionary substrate below. Each organism launches a separate single-domain SimulationWorld with its own regime sequence, physical controls, and telemetry.",
                )
                .color(egui::Color32::from_rgb(186, 198, 218)),
            );
            ui.add_space(6.0);
        });

        egui::TopBottomPanel::bottom("itme_menu_footer").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                ui.strong("Session Launch Options:");
                ui.separator();
                ui.label("Mode:");
                ui.selectable_value(&mut menu.mode, SimulationMode::GuidedCurriculum, "13-Regime Guided Tour");
                ui.selectable_value(&mut menu.mode, SimulationMode::OrganicSandbox, "Autonomous Sandbox");
                ui.separator();
                ui.label("Initial Regime:");
                ui.selectable_value(&mut menu.initial_regime, RegimeId::R00Static, "R(0.0) Basal");
                ui.selectable_value(&mut menu.initial_regime, RegimeId::R30FullyDarwinian, "R(3.0) Darwinian");
                ui.selectable_value(&mut menu.initial_regime, RegimeId::R30SoftInheritance, "R(3.0) Soft/Distill");
                ui.selectable_value(&mut menu.initial_regime, RegimeId::R35Symbiogenetic, "R(3.5) Integrated");
                ui.separator();
                ui.label("Seed:");
                ui.add(egui::TextEdit::singleline(&mut menu.seed_text).desired_width(90.0));
                if menu.mode == SimulationMode::GuidedCurriculum {
                    ui.separator();
                    ui.add(egui::Slider::new(&mut menu.steps_per_stage, 12..=80).text("Ticks/Regime"));
                }
            });
            ui.add_space(4.0);
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.add_space(10.0);
            ui.columns(3, |cols| {
                for (idx, &org) in OrganismKind::ALL.iter().enumerate() {
                    let (badge, subtitle, accent, highlights, dynamics_note) = Self::domain_card_meta(org);
                    let col_ui = &mut cols[idx];
                    let inner_w = (col_ui.available_width() - 36.0).max(220.0);
                    egui::Frame::none()
                        .fill(egui::Color32::from_rgb(19, 27, 46))
                        .stroke(egui::Stroke::new(1.5_f32, accent))
                        .rounding(10.0)
                        .inner_margin(egui::Margin::same(16.0))
                        .show(col_ui, |card| {
                            card.set_min_width(inner_w);
                            card.set_max_width(inner_w);
                            card.set_min_height(460.0);

                            card.label(
                                egui::RichText::new(badge)
                                    .small()
                                    .strong()
                                    .color(accent),
                            );
                            card.add_space(2.0);
                            card.heading(
                                egui::RichText::new(org.title())
                                    .size(18.0)
                                    .strong()
                                    .color(egui::Color32::from_rgb(241, 245, 249)),
                            );
                            card.add(
                                egui::Label::new(
                                    egui::RichText::new(subtitle)
                                        .small()
                                        .italics()
                                        .color(egui::Color32::from_rgb(148, 163, 184)),
                                )
                                .wrap(),
                            );
                            card.add_space(8.0);
                            card.separator();
                            card.add_space(6.0);

                            card.strong(egui::RichText::new("Symbol Alphabet").color(accent));
                            card.add(
                                egui::Label::new(
                                    egui::RichText::new(org.alphabet_description())
                                        .small()
                                        .monospace()
                                        .color(egui::Color32::from_rgb(203, 213, 225)),
                                )
                                .wrap(),
                            );
                            card.add_space(8.0);

                            card.strong(egui::RichText::new("Domain Regimes & Phenomena").color(accent));
                            for h in highlights {
                                card.horizontal(|row| {
                                    row.label(egui::RichText::new("•").strong().color(accent));
                                    row.add(
                                        egui::Label::new(
                                            egui::RichText::new(*h)
                                                .small()
                                                .color(egui::Color32::from_rgb(215, 224, 238)),
                                        )
                                        .wrap(),
                                    );
                                });
                            }

                            card.add_space(8.0);
                            card.strong(egui::RichText::new("Independent Dynamics").color(accent));
                            card.add(
                                egui::Label::new(
                                    egui::RichText::new(dynamics_note)
                                        .small()
                                        .color(egui::Color32::from_rgb(162, 178, 204)),
                                )
                                .wrap(),
                            );

                            card.add_space(14.0);
                            card.strong(
                                egui::RichText::new(format!(
                                    "Preview Regime {}: {}",
                                    menu.initial_regime.code(),
                                    menu.initial_regime.domain_name(org)
                                ))
                                .small()
                                .color(egui::Color32::from_rgb(226, 232, 240)),
                            );
                            card.add_space(8.0);

                            let launch_btn = egui::Button::new(
                                egui::RichText::new(format!("Launch {} Simulation ->", org.short_title()))
                                    .size(14.0)
                                    .strong()
                                    .color(egui::Color32::from_rgb(10, 15, 28)),
                            )
                            .fill(accent);
                            if card.add_sized([inner_w, 34.0], launch_btn).clicked() {
                                let seed = menu.parsed_seed();
                                launched = Some(DomainSession::new(
                                    org,
                                    menu.initial_regime,
                                    menu.mode,
                                    menu.steps_per_stage,
                                    seed,
                                ));
                            }
                        });
                }
            });
        });

        launched
    }
}

impl App for ItmeGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if let AppScreen::Simulation(session) = &mut self.screen {
            if !session.paused {
                let interval = Duration::from_secs_f32((1.0 / session.steps_per_second).max(0.015));
                if session.last_step_instant.elapsed() >= interval {
                    session.step_once();
                    session.last_step_instant = Instant::now();
                }
                if session.auto_rotate {
                    session.camera.yaw += 0.008;
                }
                ctx.request_repaint_after(Duration::from_millis(16));
            } else if session.auto_rotate {
                session.camera.yaw += 0.006;
                ctx.request_repaint_after(Duration::from_millis(16));
            }
        }

        self.draw_panels(ctx);
    }
}

impl DomainSession {
    fn new(
        organism: OrganismKind,
        initial_regime: RegimeId,
        mode: SimulationMode,
        steps_per_stage: usize,
        seed: u64,
    ) -> Self {
        let mut world = SimulationWorld::new(
            organism,
            initial_regime,
            mode,
            steps_per_stage,
            seed,
        );
        world.controls = initial_regime.domain_preset_controls(organism);
        world.step();

        Self {
            organism,
            world,
            viewport_layout: ViewportLayout::SpatialAndPhylogeny,
            camera: Camera3D::default(),
            auto_rotate: true,
            show_labels: true,
            paused: false,
            steps_per_second: 8.0,
            last_step_instant: Instant::now(),
            export_status: None,
            grand_epochs: Vec::new(),
            wants_main_menu: false,
        }
    }

    fn step_once(&mut self) {
        let prev_reg = self.world.target_regime;
        self.world.step();
        if self.world.mode == SimulationMode::GuidedCurriculum && self.world.target_regime != prev_reg {
            self.world.controls = self.world.target_regime.domain_preset_controls(self.organism);
        }
    }

    fn jump_to_regime(&mut self, regime: RegimeId) {
        let preserve = regime.index() >= 2;
        self.world.apply_regime_preset(regime, preserve);
        self.world.controls = regime.domain_preset_controls(self.organism);
        self.world.step();
    }

    fn export_snapshot_files(&mut self) {
        let slug = self.organism.slug();
        let step = self.world.step_count;
        let nhx_file = format!("itme_{slug}_step{step}.nhx");
        let svg_file = format!("itme_{slug}_step{step}_phylogram.svg");
        let nhx_res = fs::write(&nhx_file, self.world.export_newick(true));
        let svg_res = fs::write(&svg_file, self.world.export_phylogram_svg());
        self.export_status = Some(match (nhx_res, svg_res) {
            (Ok(()), Ok(())) => format!("Saved {nhx_file} & {svg_file}"),
            (Err(e), _) | (_, Err(e)) => format!("Export error: {e}"),
        });
    }

    fn draw_panels(&mut self, ctx: &egui::Context) {
        self.draw_top_toolbar(ctx);
        self.draw_left_controls(ctx);
        self.draw_right_inspector(ctx);
        self.draw_bottom_telemetry_plots(ctx);
        self.draw_central_viewport(ctx);
    }

    fn draw_top_toolbar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("itme_top_bar").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                if ui
                    .button(egui::RichText::new("<- Main Menu").strong())
                    .clicked()
                {
                    self.wants_main_menu = true;
                }
                ui.separator();

                let (_, _, accent, _, _) = ItmeGuiApp::domain_card_meta(self.organism);
                ui.label(
                    egui::RichText::new(format!("ITME v23.0 • {}", self.organism.title()))
                        .color(accent)
                        .strong(),
                );
                ui.separator();

                let play_label = if self.paused { "Play" } else { "Pause" };
                if ui.button(play_label).clicked() {
                    self.paused = !self.paused;
                    self.last_step_instant = Instant::now();
                }
                if ui.button("Step").clicked() {
                    self.step_once();
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
                ui.label("View:");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialAndPhylogeny,
                    "Spatial + Network",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialGridAnd3D,
                    "Spatial + 3D",
                );
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::PhylogenyFull,
                    "Full Lineage",
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
                if self.organism == OrganismKind::RnaDnaProtein {
                    ui.selectable_value(
                        &mut self.viewport_layout,
                        ViewportLayout::GrandTrajectory,
                        "4.25 Ga Timeline",
                    );
                }

                ui.separator();
                ui.checkbox(&mut self.auto_rotate, "Auto-Orbit");
                ui.checkbox(&mut self.show_labels, "Labels");
            });

            ui.add_space(2.0);
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    egui::RichText::new("Domain Milestones:")
                        .small()
                        .strong()
                        .color(egui::Color32::from_rgb(148, 163, 184)),
                );
                match self.organism {
                    OrganismKind::RnaDnaProtein => {
                        if ui.button("4.25 Ga Grand Trajectory").clicked() {
                            let epochs = self.world.run_continuous_grand_trajectory(10);
                            self.grand_epochs = epochs.clone();
                            self.viewport_layout = ViewportLayout::GrandTrajectory;
                            if let Some(last_ep) = epochs.last() {
                                self.export_status = Some(format!(
                                    "8-Epoch Trajectory -> {}: L_I={:.2}, MI={:.1}b, O2={:.2}",
                                    last_ep.epoch.code(),
                                    last_ep.informational_leverage_li,
                                    last_ep.heritable_mutual_info_bits,
                                    last_ep.planetary_o2
                                ));
                            }
                        }
                        if ui.button("Hadean Abiogenesis -> LUCA").clicked() {
                            let (rf, lc, eig) = self.world.trigger_hadean_abiogenesis_pulse();
                            self.export_status = Some(format!(
                                "Hadean Abiogenesis: +{rf} RAF Protocells, +{lc} LUCA (Eigen L_max={eig} nt)"
                            ));
                        }
                        if ui.button("Homochiral D/L + SGC").clicked() {
                            let (hm, sgc, mp) = self.world.trigger_homochiral_sgc_pulse();
                            self.export_status = Some(format!(
                                "Homochirality & SGC: +{hm} D-RNA breaks, +{sgc} 20-AA SGC (MP={:.1}%)",
                                mp * 100.0
                            ));
                        }
                        if ui.button("RNR RNA->T-DNA").clicked() {
                            let (rn, td, kb) = self.world.trigger_rna_to_dna_virosphere_pulse();
                            self.export_status = Some(format!(
                                "RNR & UDG Repair: +{rn} U-DNA, +{td} T-DNA Genomes (Cap={kb:.0} kb)"
                            ));
                        }
                        if ui.button("1°/2° Symbio + MRO").clicked() {
                            let fused = self.world.trigger_endosymbiotic_pulse();
                            let ser = self.world.trigger_serial_endosymbiosis_pulse();
                            let (hy, mi, am) = self.world.trigger_mro_reduction_pulse();
                            self.export_status = Some(format!(
                                "Symbiosis & MRO: 1°={fused}, Serial={ser:?}, +Hy:{hy} +MitoS:{mi} +SUF(MRO->0):{am}"
                            ));
                        }
                    }
                    OrganismKind::HominidPopulation => {
                        if ui.button("Kin Division & Lithic Guilds").clicked() {
                            self.jump_to_regime(RegimeId::R30FullyDarwinian);
                            self.export_status = Some(
                                "Promoted foraging bands to kin-structured lithic guilds R(3.0-DAR)."
                                    .to_string(),
                            );
                        }
                        if ui.button("Symbolic Teaching + Ratchet").clicked() {
                            self.world.controls.symbolic_teaching_fidelity =
                                self.world.controls.symbolic_teaching_fidelity.max(0.90);
                            self.world.controls.inter_deme_connectivity =
                                self.world.controls.inter_deme_connectivity.clamp(0.25, 0.55);
                            let (rt, sym, rec) = self.world.trigger_cultural_ratchet_pulse();
                            self.export_status = Some(format!(
                                "Symbolic Culture: +{rt} ratchets, +{sym} pedagogical lineages, +{rec} deme recombinations"
                            ));
                        }
                        if ui.button("Agro-Ectosymbiosis & Crops").clicked() {
                            let (ag, ob, sh) =
                                self.world.trigger_agriculture_domestication_pulse();
                            self.export_status = Some(format!(
                                "Cultivar Care: +{ag} farming demes, +{ob} clonal crops, {sh} defensive shields"
                            ));
                        }
                        if ui.button("Written Archive + Peer Audit").clicked() {
                            let (ar, inst, aud) = self.world.trigger_exosomatic_archive_pulse();
                            self.export_status = Some(format!(
                                "Written Archives: +{ar} archives, +{inst} peer institutions, {aud} audits"
                            ));
                        }
                    }
                    OrganismKind::ArtificialIntelligence => {
                        if ui.button("Container Checkpoint Scale").clicked() {
                            self.jump_to_regime(RegimeId::R15Replicative);
                            self.export_status = Some(
                                "Scaled stateless container weights across cluster R(1.5-RP)."
                                    .to_string(),
                            );
                        }
                        if ui.button("Adversarial Arena (PBT)").clicked() {
                            self.jump_to_regime(RegimeId::R30FullyDarwinian);
                            self.export_status = Some(
                                "Activated population-based training & adversarial arena R(3.0-DAR)."
                                    .to_string(),
                            );
                        }
                        if ui.button("Formal Proof Verifier").clicked() {
                            let (ms, fm, prf) = self.world.trigger_machine_synthesis_pulse();
                            self.export_status = Some(format!(
                                "Formal Verifiers: +{ms} compute symbionts, +{fm} verified loops, {prf} proofs checked"
                            ));
                        }
                        if ui.button("Anti-Collapse Safeguard").clicked() {
                            self.world.controls.empirical_grounding_floor =
                                self.world.controls.empirical_grounding_floor.max(0.36);
                            self.world.controls.formal_verifier_strictness =
                                self.world.controls.formal_verifier_strictness.max(0.92);
                            self.world.controls.unverified_publish_pressure =
                                self.world.controls.unverified_publish_pressure.min(0.18);
                            for _ in 0..4 {
                                self.step_once();
                            }
                            let m = self.world.latest_metrics();
                            self.export_status = Some(format!(
                                "Anti-Collapse Guard active: empiric floor >= 0.36 | TailDiv = {:.2}, Theorems = {:.1}b",
                                m.mean_tail_diversity, m.mean_verified_theorem_bits
                            ));
                        }
                    }
                }
            });
            ui.add_space(4.0);
        });
    }

    fn draw_left_controls(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("itme_left_controls")
            .default_width(340.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Simulation Mode");
                    let mut selected_mode = self.world.mode;
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
                    self.world.mode = selected_mode;

                    if selected_mode == SimulationMode::GuidedCurriculum {
                        ui.add(
                            egui::Slider::new(&mut self.world.steps_per_stage, 12..=120)
                                .text("Ticks / Regime"),
                        );
                    } else if ui.button("Stochastic Phase Jump").clicked() {
                        self.world.trigger_network_transition();
                    }

                    ui.separator();
                    ui.heading(format!("{} Regimes (13)", self.organism.short_title()));
                    let active_classified = self.world.classified_regime;
                    let mut clicked_regime = None;
                    for regime in RegimeId::ALL {
                        let is_curr = regime == active_classified;
                        let txt = regime.domain_title(self.organism);
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
                        self.jump_to_regime(reg);
                    }

                    ui.separator();
                    let theta_f = self.world.latest_metrics().critical_fidelity_theta;
                    let w = &mut self.world;
                    let mut changed = false;

                    ui.heading("Energetics & Transmission");
                    changed |= ui
                        .add(
                            egui::Slider::new(&mut w.controls.bath_exergy_flux, -0.6..=1.5)
                                .text("J_in (Supply / Exergy Flux)"),
                        )
                        .changed();
                    changed |= ui
                        .add(
                            egui::Slider::new(&mut w.controls.replication_fidelity, 0.40..=1.0)
                                .text("F_0 (Copy / Update Fidelity)"),
                        )
                        .changed();
                    changed |= ui
                        .add(
                            egui::Slider::new(&mut w.controls.landauer_bit_cost, 0.0..=0.008)
                                .text("k_B T ln2 (Maintenance Cost)"),
                        )
                        .changed();

                    ui.add_space(4.0);
                    match self.organism {
                        OrganismKind::RnaDnaProtein => {
                            ui.heading("Prebiotic & Genetic Code");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.wet_dry_cycle_amplitude,
                                        0.0..=1.0,
                                    )
                                    .text("A_wd (Wet-Dry Condensation)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.hydrothermal_delta_ph,
                                        0.0..=4.5,
                                    )
                                    .text("Delta_pH (Vent FeNiS Gradient)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.magnetite_ciss_spin_bias,
                                        0.0..=1.0,
                                    )
                                    .text("b_ciss (Fe3O4 Chiral Spin Bias)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.rnr_thiyl_radical_activity,
                                        0.0..=1.0,
                                    )
                                    .text("k_rnr (Class III RNR 2'-OH->2'-H)"),
                                )
                                .changed();

                            ui.add_space(4.0);
                            ui.heading("Symbiogenesis & MRO Continuum");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.endosymbiosis_rate,
                                        0.0..=0.80,
                                    )
                                    .text("gamma_1 (1° Endosymbiosis)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.ambient_oxygen_level,
                                        0.0..=1.0,
                                    )
                                    .text("O2_level (Hypoxia < 0.45)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.cytosolic_suf_hgt_rescue,
                                        0.0..=1.0,
                                    )
                                    .text("suf_rescue (Cytosolic SUF >=0.65)"),
                                )
                                .changed();
                        }
                        OrganismKind::HominidPopulation => {
                            ui.heading("Cultural Transmission & Demes");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.symbolic_teaching_fidelity,
                                        0.0..=1.0,
                                    )
                                    .text("f_teach (Symbolic Pedagogy)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.inter_deme_connectivity,
                                        0.0..=1.0,
                                    )
                                    .text("c_deme (Inter-Band Contact)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.cultural_innovation_rate,
                                        0.0..=0.90,
                                    )
                                    .text("mu_innov (Innovation Drive)"),
                                )
                                .changed();

                            ui.add_space(4.0);
                            ui.heading("Subsistence & Written Archives");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.crop_domestication_rate,
                                        0.0..=0.90,
                                    )
                                    .text("mu_crop (Cultivar Domestication)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.archival_externalization_rate,
                                        0.0..=0.90,
                                    )
                                    .text("mu_arch (Written Archival Drive)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.peer_verification_rigor,
                                        0.0..=1.0,
                                    )
                                    .text("r_peer (Peer Replication Rigor)"),
                                )
                                .changed();
                        }
                        OrganismKind::ArtificialIntelligence => {
                            ui.heading("Training, Plasticity & Distillation");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.synthetic_generation_rate,
                                        0.0..=1.0,
                                    )
                                    .text("g_synth (Synthetic Rollout Rate)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.memory_plasticity,
                                        0.0..=1.0,
                                    )
                                    .text("M_plast (In-Context / LoRA)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.soft_inheritance_coupling,
                                        0.0..=1.0,
                                    )
                                    .text("alpha_distill (Teacher Distillation)"),
                                )
                                .changed();

                            ui.add_space(4.0);
                            ui.heading("Formal Verification & Anti-Collapse");
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.formal_verifier_strictness,
                                        0.0..=1.0,
                                    )
                                    .text("v_formal (Proof Verifier >=0.72)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.empirical_grounding_floor,
                                        0.0..=0.60,
                                    )
                                    .text("gamma_emp (Empirical Anchor >=0.12)"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut w.controls.unverified_publish_pressure,
                                        0.0..=1.0,
                                    )
                                    .text("p_noise (Unverified Publish Rate)"),
                                )
                                .changed();
                        }
                    }

                    let f_eff = w.controls.effective_fidelity();
                    ui.add_space(4.0);
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
                    if changed {
                        w.mode = SimulationMode::OrganicSandbox;
                    }
                });
            });
    }

    fn draw_right_inspector(&mut self, ctx: &egui::Context) {
        egui::SidePanel::right("itme_right_inspector")
            .default_width(365.0)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    let w = &self.world;
                    let m = w.latest_metrics();
                    let reg = w.classified_regime;

                    ui.heading("Emergent Domain Regime");
                    ui.label(
                        egui::RichText::new(reg.domain_title(self.organism))
                            .size(15.5)
                            .strong()
                            .color(egui::Color32::from_rgb(255, 214, 102)),
                    );
                    ui.add(
                        egui::Label::new(
                            egui::RichText::new(reg.domain_concept(self.organism)).italics(),
                        )
                        .wrap(),
                    );
                    ui.monospace(format!("Diagnostic {}", m.diagnostic));

                    let phase_col = match m.phase {
                        OrderPhase::Adaptive => egui::Color32::from_rgb(95, 235, 165),
                        OrderPhase::Stasis => egui::Color32::from_rgb(135, 195, 255),
                        OrderPhase::Regressive => egui::Color32::from_rgb(255, 110, 110),
                    };
                    ui.colored_label(phase_col, format!("Phase: {}", m.phase.label()));

                    if let Some(ref status) = self.export_status {
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(status)
                                    .small()
                                    .color(egui::Color32::from_rgb(130, 245, 205)),
                            )
                            .wrap(),
                        );
                    }

                    ui.separator();
                    ui.heading(format!("{} Telemetry", self.organism.short_title()));
                    egui::Grid::new("itme_domain_metrics_grid")
                        .striped(true)
                        .num_columns(2)
                        .show(ui, |ui| {
                            ui.label("Order O / O_pop");
                            ui.monospace(format!("{:.2} / {:.2}", m.order_o, m.ecosystem_order_pop));
                            ui.end_row();

                            ui.label("MI (bits) / Adami rho");
                            ui.monospace(format!(
                                "{:.2}b / {:.3}",
                                m.heritable_mutual_info_bits, m.adami_conservation_rho
                            ));
                            ui.end_row();

                            match self.organism {
                                OrganismKind::RnaDnaProtein => {
                                    ui.label("Oligomer / Eigen L_max");
                                    ui.monospace(format!(
                                        "{:.1} nt / {} nt",
                                        m.mean_rna_oligomer_nt, m.eigen_max_length_nt
                                    ));
                                    ui.end_row();

                                    ui.label("Genome Cap / D-ee");
                                    ui.monospace(format!(
                                        "{:.1} kb / {:.0}%",
                                        m.genome_capacity_kb,
                                        m.ribose_ee_d * 100.0
                                    ));
                                    ui.end_row();

                                    ui.label("MROs (Aero/Hy/Mi/0)");
                                    ui.monospace(format!(
                                        "{}/{}/{}/{} (P:{:.1})",
                                        m.aerobic_mito_count,
                                        m.hydrogenosome_count,
                                        m.mitosome_count,
                                        m.amiteukaryote_count,
                                        m.mean_mro_pathways
                                    ));
                                    ui.end_row();

                                    ui.label("Apicoplast Lineages");
                                    ui.monospace(format!(
                                        "{} ({:.1} paths)",
                                        m.apicoplast_parasite_count, m.mean_apicoplast_pathways
                                    ));
                                    ui.end_row();
                                }
                                OrganismKind::HominidPopulation => {
                                    ui.label("Population / Species");
                                    ui.monospace(format!("{} / {} clades", m.population_size, m.species_count));
                                    ui.end_row();

                                    ui.label("Deme Divergence F_ST");
                                    ui.monospace(format!(
                                        "{:.3} (Rec: {})",
                                        m.fst_fixation_index, m.cross_deme_recombinations
                                    ));
                                    ui.end_row();

                                    ui.label("Cultural Rep / Teach");
                                    ui.monospace(format!(
                                        "{:.1}b / {:.2}",
                                        m.mean_cultural_repertoire_bits,
                                        m.mean_neural_teaching_capacity
                                    ));
                                    ui.end_row();

                                    ui.label("Written Archive / FDR");
                                    ui.monospace(format!(
                                        "{:.1}b / {:.1}%",
                                        m.mean_exosomatic_archive_bits,
                                        m.mean_archive_fdr_load * 100.0
                                    ));
                                    ui.end_row();
                                }
                                OrganismKind::ArtificialIntelligence => {
                                    ui.label("Active Workers / Symbionts");
                                    ui.monospace(format!(
                                        "{} / {} ({} verified)",
                                        m.population_size,
                                        m.computational_symbiont_count,
                                        m.formal_symbio_count
                                    ));
                                    ui.end_row();

                                    ui.label("Verified Theorem Bits");
                                    ui.monospace(format!("{:.2} bits", m.mean_verified_theorem_bits));
                                    ui.end_row();

                                    ui.label("Tail Diversity (Anti-Collapse)");
                                    ui.monospace(format!("{:.3}", m.mean_tail_diversity));
                                    ui.end_row();

                                    ui.label("Proofs Checked / Collapses");
                                    ui.monospace(format!(
                                        "{} / {}",
                                        m.formal_proofs_verified, m.recursive_model_collapses
                                    ));
                                    ui.end_row();
                                }
                            }
                        });

                    ui.separator();
                    ui.heading("Domain Capabilities");
                    ui.horizontal_wrapped(|ui| {
                        flag_chip(ui, "Driven", m.flags.exergy_driven);
                        flag_chip(ui, "Replicate", m.flags.replication);
                        flag_chip(ui, "Select", m.flags.selection);
                        match self.organism {
                            OrganismKind::RnaDnaProtein => {
                                flag_chip(ui, "RAF_Closure", m.flags.prebiotic_raf_closure);
                                flag_chip(ui, "RiboHypercycle", m.flags.ribozyme_hypercycle_active);
                                flag_chip(ui, "Homochiral_DL", m.flags.homochiral_d_rna_l_aa);
                                flag_chip(ui, "UniversalSGC", m.flags.error_minimized_universal_sgc);
                                flag_chip(ui, "RNR_dNTP", m.flags.rnr_deoxyribose_synthesis);
                                flag_chip(ui, "ThyX_UDG", m.flags.thymine_udg_deamination_repair);
                                flag_chip(ui, "DarwinianLUCA", m.flags.darwinian_threshold_luca);
                                flag_chip(ui, "1°Holo", m.flags.endosymbiotic_fusion);
                                flag_chip(ui, "2°/3°Serial", m.flags.serial_endosymbiosis);
                                flag_chip(ui, "Hydrogenosome", m.flags.hydrogenosome_fermentation);
                                flag_chip(ui, "Mitosome(ISC)", m.flags.mitosome_isc_retention);
                                flag_chip(ui, "MRO->0(SUF)", m.flags.amiteukaryote_mro_loss);
                                flag_chip(ui, "ClonalMC", m.flags.clonal_multicellularity);
                                flag_chip(ui, "GermSoma", m.flags.germ_soma_division);
                            }
                            OrganismKind::HominidPopulation => {
                                flag_chip(ui, "Plasticity", m.flags.memory);
                                flag_chip(ui, "NicheBuild", m.flags.niche_construction);
                                flag_chip(ui, "DualInherit", m.flags.dual_gene_culture_inheritance);
                                flag_chip(ui, "CultRatchet", m.flags.cumulative_cultural_ratchet);
                                flag_chip(ui, "SymbolicTeach", m.flags.symbolic_pedagogical_teaching);
                                flag_chip(ui, "AgriDomest", m.flags.agricultural_domestication);
                                flag_chip(ui, "ClonalCrop", m.flags.clonal_crop_monoculture);
                                flag_chip(ui, "ExoOffload", m.flags.exosomatic_memory_offloading);
                                flag_chip(ui, "WrittenArch", m.flags.written_archival_storage);
                                flag_chip(ui, "PeerVerified", m.flags.institutional_peer_verification);
                            }
                            OrganismKind::ArtificialIntelligence => {
                                flag_chip(ui, "WeightPlastic", m.flags.memory);
                                flag_chip(ui, "Distillation", m.flags.soft_inheritance);
                                flag_chip(ui, "VectorArchive", m.flags.written_archival_storage);
                                flag_chip(ui, "PeerAudit", m.flags.institutional_peer_verification);
                                flag_chip(ui, "CompSymbiont", m.flags.computational_symbiont);
                                flag_chip(ui, "FormalProof", m.flags.formal_proof_verification);
                                flag_chip(ui, "CollapseGuard", m.flags.model_collapse_safeguard);
                            }
                        }
                    });

                    ui.separator();
                    ui.heading("Live Extended Newick (NHX)");
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
            .min_height(145.0)
            .default_height(155.0)
            .resizable(true)
            .show(ctx, |ui| {
                let w = &self.world;
                ui.horizontal(|ui| {
                    ui.strong(format!("ITME v23.0 Trajectories — {}", self.organism.title()));
                    ui.separator();
                    ui.colored_label(egui::Color32::from_rgb(80, 220, 255), "Order (O)");
                    match self.organism {
                        OrganismKind::RnaDnaProtein => {
                            ui.colored_label(egui::Color32::from_rgb(6, 214, 160), "Hydrogenosomes");
                            ui.colored_label(egui::Color32::from_rgb(255, 209, 102), "Mitosomes (ISC=1)");
                        }
                        OrganismKind::HominidPopulation => {
                            ui.colored_label(egui::Color32::from_rgb(6, 214, 160), "Cultural Repertoire (bits)");
                            ui.colored_label(egui::Color32::from_rgb(255, 209, 102), "Written Archive (bits)");
                        }
                        OrganismKind::ArtificialIntelligence => {
                            ui.colored_label(egui::Color32::from_rgb(6, 214, 160), "Verified Theorem Bits");
                            ui.colored_label(egui::Color32::from_rgb(199, 125, 255), "Tail Diversity (x10)");
                        }
                    }
                });

                let avail = ui.available_size();
                let (_resp, painter) = ui.allocate_painter(avail, egui::Sense::hover());
                let rect = painter.clip_rect().shrink(4.0);
                paint_framed_box(
                    &painter,
                    rect,
                    6.0,
                    egui::Color32::from_rgb(10, 14, 26),
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(38, 52, 84)),
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

                match self.organism {
                    OrganismKind::RnaDnaProtein => {
                        let max_right = w
                            .history
                            .iter()
                            .map(|m| m.hydrogenosome_count.max(m.mitosome_count) as f32)
                            .fold(10.0_f32, f32::max)
                            * 1.15;
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
                            format!("Hydrogenosomes & Mitosomes [0 .. {:.0}]", max_right),
                            egui::FontId::monospace(10.5),
                            egui::Color32::from_rgb(175, 190, 220),
                        );
                    }
                    OrganismKind::HominidPopulation => {
                        let max_right = w
                            .history
                            .iter()
                            .map(|m| m.mean_cultural_repertoire_bits.max(m.mean_exosomatic_archive_bits) as f32)
                            .fold(8.0_f32, f32::max)
                            * 1.15;
                        draw_series(
                            &painter,
                            right_plot,
                            &w.history,
                            |m| m.mean_cultural_repertoire_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(6, 214, 160),
                        );
                        draw_series(
                            &painter,
                            right_plot,
                            &w.history,
                            |m| m.mean_exosomatic_archive_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(255, 209, 102),
                        );
                        painter.text(
                            right_plot.left_top() + egui::vec2(6.0, 4.0),
                            egui::Align2::LEFT_TOP,
                            format!("Cultural & Written Archive Bits [0 .. {:.1}]", max_right),
                            egui::FontId::monospace(10.5),
                            egui::Color32::from_rgb(175, 190, 220),
                        );
                    }
                    OrganismKind::ArtificialIntelligence => {
                        let max_right = w
                            .history
                            .iter()
                            .map(|m| (m.mean_verified_theorem_bits as f32).max(m.mean_tail_diversity as f32 * 10.0))
                            .fold(10.0_f32, f32::max)
                            * 1.15;
                        draw_series(
                            &painter,
                            right_plot,
                            &w.history,
                            |m| m.mean_verified_theorem_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(6, 214, 160),
                        );
                        draw_series(
                            &painter,
                            right_plot,
                            &w.history,
                            |m| (m.mean_tail_diversity as f32) * 10.0,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(199, 125, 255),
                        );
                        painter.text(
                            right_plot.left_top() + egui::vec2(6.0, 4.0),
                            egui::Align2::LEFT_TOP,
                            format!("Verified Theorem Bits & TailDiv [0 .. {:.1}]", max_right),
                            egui::FontId::monospace(10.5),
                            egui::Color32::from_rgb(175, 190, 220),
                        );
                    }
                }
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
                egui::Stroke::new(1.2_f32, egui::Color32::from_rgb(36, 52, 88)),
            );

            match self.viewport_layout {
                ViewportLayout::SpatialAndPhylogeny => {
                    let w = &self.world;
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
                        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_framed_box(
                        &painter,
                        right_r,
                        6.0,
                        egui::Color32::from_rgb(11, 17, 30),
                        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_spatial_grid(&painter, left_r.shrink(4.0), w, self.show_labels);
                    paint_phylogram_tree(&painter, right_r.shrink(4.0), w, self.show_labels);
                }
                ViewportLayout::SpatialGridAnd3D => {
                    let w = &self.world;
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
                        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(42, 60, 98)),
                    );
                    paint_spatial_grid(&painter, left_r.shrink(4.0), w, self.show_labels);
                    paint_scene_3d(&painter, right_r, &scene, &self.camera, self.show_labels);
                }
                ViewportLayout::PhylogenyFull => {
                    paint_phylogram_tree(&painter, full_rect.shrink(8.0), &self.world, self.show_labels);
                }
                ViewportLayout::Perspective3D => {
                    let scene = self.world.current_scene();
                    paint_scene_3d(
                        &painter,
                        full_rect.shrink(6.0),
                        &scene,
                        &self.camera,
                        self.show_labels,
                    );
                }
                ViewportLayout::Blueprint2D => {
                    let scene = self.world.current_scene();
                    paint_scene_2d(&painter, full_rect.shrink(6.0), &scene, self.show_labels);
                }
                ViewportLayout::GrandTrajectory => {
                    self.paint_grand_trajectory_viewport(&painter, full_rect);
                }
            }
        });
    }

    fn paint_grand_trajectory_viewport(&self, painter: &egui::Painter, full_rect: egui::Rect) {
        let panel_rect = full_rect.shrink(4.0);
        paint_framed_box(
            painter,
            panel_rect,
            6.0,
            egui::Color32::from_rgb(11, 16, 30),
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(38, 54, 90)),
        );
        let painter = &painter.with_clip_rect(panel_rect.shrink(2.0));

        if self.grand_epochs.is_empty() {
            painter.text(
                panel_rect.center(),
                egui::Align2::CENTER_CENTER,
                "Click '4.25 Ga -> Gaia-2.0 Grand Run!' in the top bar to execute the 8-Epoch Evolutionary Trajectory.",
                egui::FontId::proportional(14.0),
                egui::Color32::from_rgb(148, 163, 184),
            );
            return;
        }

        let last = self.grand_epochs.last().unwrap();
        painter.text(
            panel_rect.left_top() + egui::vec2(12.0, 9.0),
            egui::Align2::LEFT_TOP,
            "4.25 Ga -> Gaia-2.0 Continuous Trajectory",
            egui::FontId::proportional(13.0),
            egui::Color32::from_rgb(230, 240, 255),
        );
        painter.text(
            panel_rect.right_top() + egui::vec2(-12.0, 10.0),
            egui::Align2::RIGHT_TOP,
            format!(
                "Final: Ch:{}/5 | MI:{:.1}b | L_I:{:.2} | O2:{:.3} | CO2:{:.0} ppm",
                last.heritable_channels_active,
                last.heritable_mutual_info_bits,
                last.informational_leverage_li,
                last.planetary_o2,
                last.planetary_co2_ppm
            ),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(125, 211, 252),
        );

        let content = egui::Rect::from_min_max(
            egui::pos2(panel_rect.left() + 8.0, panel_rect.top() + 34.0),
            egui::pos2(panel_rect.right() - 8.0, panel_rect.bottom() - 8.0),
        );
        let charts_h = (content.height() * 0.50).max(200.0);
        let charts_rect = egui::Rect::from_min_max(
            content.left_top(),
            egui::pos2(content.right(), content.top() + charts_h),
        );
        let table_rect = egui::Rect::from_min_max(
            egui::pos2(content.left(), charts_rect.bottom() + 8.0),
            content.right_bottom(),
        );

        let mid_x = charts_rect.center().x;
        let left_chart = egui::Rect::from_min_max(
            charts_rect.left_top(),
            egui::pos2(mid_x - 5.0, charts_rect.bottom()),
        );
        let right_chart = egui::Rect::from_min_max(
            egui::pos2(mid_x + 5.0, charts_rect.top()),
            charts_rect.right_bottom(),
        );

        let n = self.grand_epochs.len();
        let draw_epoch_axis = |p: &egui::Painter, plot: egui::Rect| {
            for (idx, ep) in self.grand_epochs.iter().enumerate() {
                let frac = if n <= 1 {
                    0.5
                } else {
                    idx as f32 / (n - 1) as f32
                };
                let x = plot.left() + frac * plot.width();
                p.line_segment(
                    [egui::pos2(x, plot.top()), egui::pos2(x, plot.bottom())],
                    egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(30, 41, 68)),
                );
                let age_lbl = if ep.geologic_age_ga >= 0.5 {
                    format!("E{}\n{:.1}G", idx + 1, ep.geologic_age_ga)
                } else if ep.geologic_age_ga > 0.0 {
                    format!("E{}\n{:.0}M", idx + 1, ep.geologic_age_ga * 1000.0)
                } else {
                    format!("E{}\n0G", idx + 1)
                };
                p.text(
                    egui::pos2(x, plot.bottom() + 4.0),
                    egui::Align2::CENTER_TOP,
                    age_lbl,
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(148, 163, 184),
                );
            }
        };

        // Panel A: Atmospheric O2 & log10(CO2 ppm)
        paint_framed_box(
            painter,
            left_chart,
            5.0,
            egui::Color32::from_rgb(15, 22, 39),
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(45, 65, 105)),
        );
        {
            let p = painter.with_clip_rect(left_chart);
            p.text(
                left_chart.left_top() + egui::vec2(10.0, 6.0),
                egui::Align2::LEFT_TOP,
                "Panel A | Atmosphere: Free O2 (Cyan) & log10(CO2) (Amber)",
                egui::FontId::proportional(10.5),
                egui::Color32::from_rgb(226, 232, 240),
            );
            let plot = egui::Rect::from_min_max(
                egui::pos2(left_chart.left() + 46.0, left_chart.top() + 28.0),
                egui::pos2(left_chart.right() - 44.0, left_chart.bottom() - 32.0),
            );
            let max_o2 = self
                .grand_epochs
                .iter()
                .map(|e| e.planetary_o2)
                .fold(0.25_f64, f64::max)
                * 1.18;
            let log_co2s: Vec<f64> = self
                .grand_epochs
                .iter()
                .map(|e| e.planetary_co2_ppm.max(1.0).log10())
                .collect();
            let min_lc = log_co2s.iter().copied().fold(f64::INFINITY, f64::min).min(2.0);
            let max_lc = log_co2s
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                .max(min_lc + 1.0);
            for k in 0..=4 {
                let frac = k as f32 / 4.0;
                let y = plot.bottom() - frac * plot.height();
                p.line_segment(
                    [egui::pos2(plot.left(), y), egui::pos2(plot.right(), y)],
                    egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(30, 41, 68)),
                );
                p.text(
                    egui::pos2(plot.left() - 6.0, y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{:.2}", max_o2 * frac as f64),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(76, 201, 240),
                );
                p.text(
                    egui::pos2(plot.right() + 6.0, y),
                    egui::Align2::LEFT_CENTER,
                    format!("{:.1}", min_lc + (max_lc - min_lc) * frac as f64),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(244, 162, 97),
                );
            }
            draw_epoch_axis(&p, plot);

            let o2_pts: Vec<egui::Pos2> = self
                .grand_epochs
                .iter()
                .enumerate()
                .map(|(i, ep)| {
                    let fx = if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 };
                    let fy = (ep.planetary_o2 / max_o2).clamp(0.0, 1.0) as f32;
                    egui::pos2(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                })
                .collect();
            let co2_pts: Vec<egui::Pos2> = log_co2s
                .iter()
                .enumerate()
                .map(|(i, &lc)| {
                    let fx = if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 };
                    let fy = ((lc - min_lc) / (max_lc - min_lc)).clamp(0.0, 1.0) as f32;
                    egui::pos2(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                })
                .collect();
            p.add(egui::Shape::line(
                co2_pts.clone(),
                egui::Stroke::new(2.2_f32, egui::Color32::from_rgb(244, 162, 97)),
            ));
            p.add(egui::Shape::line(
                o2_pts.clone(),
                egui::Stroke::new(2.4_f32, egui::Color32::from_rgb(76, 201, 240)),
            ));
            for (i, pt) in o2_pts.iter().enumerate() {
                p.circle_filled(*pt, 3.8, egui::Color32::from_rgb(76, 201, 240));
                if i == 3 {
                    p.text(
                        *pt + egui::vec2(0.0, -9.0),
                        egui::Align2::CENTER_BOTTOM,
                        format!("O2:{:.2}", self.grand_epochs[i].planetary_o2),
                        egui::FontId::monospace(8.5),
                        egui::Color32::from_rgb(186, 230, 253),
                    );
                } else if i == n - 1 {
                    p.text(
                        *pt + egui::vec2(-6.0, -7.0),
                        egui::Align2::RIGHT_BOTTOM,
                        format!("O2:{:.2}", self.grand_epochs[i].planetary_o2),
                        egui::FontId::monospace(8.5),
                        egui::Color32::from_rgb(186, 230, 253),
                    );
                }
            }
            for (i, pt) in co2_pts.iter().enumerate() {
                p.circle_filled(*pt, 3.5, egui::Color32::from_rgb(244, 162, 97));
                if i == 0 || i == n - 1 {
                    p.text(
                        *pt + egui::vec2(if i == 0 { 16.0 } else { -12.0 }, 8.0),
                        egui::Align2::CENTER_TOP,
                        format!("{:.0}ppm", self.grand_epochs[i].planetary_co2_ppm),
                        egui::FontId::monospace(8.5),
                        egui::Color32::from_rgb(253, 186, 116),
                    );
                }
            }
        }

        // Panel B: 5-Channel Heritable Mutual Information, Leverage L_I & Active Channels
        paint_framed_box(
            painter,
            right_chart,
            5.0,
            egui::Color32::from_rgb(15, 22, 39),
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(45, 65, 105)),
        );
        {
            let p = painter.with_clip_rect(right_chart);
            p.text(
                right_chart.left_top() + egui::vec2(10.0, 6.0),
                egui::Align2::LEFT_TOP,
                "Panel B | MI (Green), L_I (Violet) & Active Ch (Gold)",
                egui::FontId::proportional(10.5),
                egui::Color32::from_rgb(226, 232, 240),
            );
            let plot = egui::Rect::from_min_max(
                egui::pos2(right_chart.left() + 46.0, right_chart.top() + 28.0),
                egui::pos2(right_chart.right() - 44.0, right_chart.bottom() - 32.0),
            );
            let max_mi = self
                .grand_epochs
                .iter()
                .map(|e| e.heritable_mutual_info_bits)
                .fold(10.0_f64, f64::max)
                * 1.24;
            let max_li = self
                .grand_epochs
                .iter()
                .map(|e| e.informational_leverage_li)
                .fold(1.0_f64, f64::max)
                * 1.24;
            for k in 0..=4 {
                let frac = k as f32 / 4.0;
                let y = plot.bottom() - frac * plot.height();
                p.line_segment(
                    [egui::pos2(plot.left(), y), egui::pos2(plot.right(), y)],
                    egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(30, 41, 68)),
                );
                p.text(
                    egui::pos2(plot.left() - 6.0, y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{:.0}b", max_mi * frac as f64),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(80, 227, 137),
                );
                p.text(
                    egui::pos2(plot.right() + 6.0, y),
                    egui::Align2::LEFT_CENTER,
                    format!("{:.1}", max_li * frac as f64),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(199, 125, 255),
                );
            }
            draw_epoch_axis(&p, plot);

            let bar_w = (plot.width() / n.max(1) as f32 * 0.42).clamp(16.0, 34.0);
            for (i, ep) in self.grand_epochs.iter().enumerate() {
                let fx = if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 };
                let x = plot.left() + fx * plot.width();
                let ch_frac = (ep.heritable_channels_active as f32 / 5.8).clamp(0.05, 0.92);
                let bar_top = plot.bottom() - ch_frac * plot.height();
                let bar = egui::Rect::from_min_max(
                    egui::pos2(x - bar_w * 0.5, bar_top),
                    egui::pos2(x + bar_w * 0.5, plot.bottom()),
                );
                p.rect_filled(bar, 2.0, egui::Color32::from_rgba_premultiplied(250, 204, 21, 30));
                p.text(
                    egui::pos2(x, plot.bottom() - 10.0),
                    egui::Align2::CENTER_CENTER,
                    format!("{}ch", ep.heritable_channels_active),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(253, 224, 71),
                );
            }

            let mi_pts: Vec<egui::Pos2> = self
                .grand_epochs
                .iter()
                .enumerate()
                .map(|(i, ep)| {
                    let fx = if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 };
                    let fy = (ep.heritable_mutual_info_bits / max_mi).clamp(0.0, 1.0) as f32;
                    egui::pos2(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                })
                .collect();
            let li_pts: Vec<egui::Pos2> = self
                .grand_epochs
                .iter()
                .enumerate()
                .map(|(i, ep)| {
                    let fx = if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 };
                    let fy = (ep.informational_leverage_li / max_li).clamp(0.0, 1.0) as f32;
                    egui::pos2(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                })
                .collect();
            p.add(egui::Shape::line(
                li_pts.clone(),
                egui::Stroke::new(2.2_f32, egui::Color32::from_rgb(199, 125, 255)),
            ));
            p.add(egui::Shape::line(
                mi_pts.clone(),
                egui::Stroke::new(2.5_f32, egui::Color32::from_rgb(80, 227, 137)),
            ));
            for (i, pt) in mi_pts.iter().enumerate() {
                p.circle_filled(*pt, 3.8, egui::Color32::from_rgb(80, 227, 137));
                if i == 2 {
                    p.text(
                        *pt + egui::vec2(-8.0, -24.0),
                        egui::Align2::CENTER_BOTTOM,
                        format!("{:.1}b", self.grand_epochs[i].heritable_mutual_info_bits),
                        egui::FontId::monospace(8.5),
                        egui::Color32::from_rgb(167, 243, 208),
                    );
                } else if i == n - 1 {
                    p.text(
                        *pt + egui::vec2(-8.0, -8.0),
                        egui::Align2::RIGHT_BOTTOM,
                        format!("{:.1}b", self.grand_epochs[i].heritable_mutual_info_bits),
                        egui::FontId::monospace(8.5),
                        egui::Color32::from_rgb(167, 243, 208),
                    );
                }
            }
            for pt in &li_pts {
                p.circle_filled(*pt, 3.4, egui::Color32::from_rgb(199, 125, 255));
            }
        }

        // Bottom Table: 8-Epoch Major Transition Telemetry Ledger
        paint_framed_box(
            painter,
            table_rect,
            5.0,
            egui::Color32::from_rgb(13, 19, 36),
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(45, 65, 105)),
        );
        let tp = painter.with_clip_rect(table_rect);
        let w = table_rect.width();
        let x_cols = [
            table_rect.left() + 10.0,
            table_rect.left() + w * 0.155,
            table_rect.left() + w * 0.235,
            table_rect.left() + w * 0.280,
            table_rect.left() + w * 0.325,
            table_rect.left() + w * 0.390,
            table_rect.left() + w * 0.445,
            table_rect.left() + w * 0.515,
            table_rect.left() + w * 0.575,
            table_rect.left() + w * 0.635,
            table_rect.left() + w * 0.700,
        ];
        let headers = [
            "Epoch",
            "Age",
            "Ch",
            "Pop",
            "MI (b)",
            "L_I",
            "v_O",
            "O2",
            "CO2",
            "Exergy",
            "Landmark Transition",
        ];
        let hdr_y = table_rect.top() + 9.0;
        for (x, h) in x_cols.iter().zip(headers.iter()) {
            tp.text(
                egui::pos2(*x, hdr_y),
                egui::Align2::LEFT_TOP,
                *h,
                egui::FontId::monospace(9.2),
                egui::Color32::from_rgb(125, 211, 252),
            );
        }
        tp.line_segment(
            [
                egui::pos2(table_rect.left() + 8.0, hdr_y + 15.0),
                egui::pos2(table_rect.right() - 8.0, hdr_y + 15.0),
            ],
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(45, 65, 105)),
        );

        let body_top = hdr_y + 19.0;
        let row_h = ((table_rect.bottom() - body_top - 6.0) / n.max(1) as f32).clamp(18.0, 28.0);
        let short_epoch_labels = [
            "E1 Homochiral",
            "E2 SGC + RAF",
            "E3 T-DNA / LUCA",
            "E4 GOE / LECA",
            "E5 Germ-Soma",
            "E6 Agri Superorg",
            "E7 Archives + AI",
            "E8 Gaia-2.0",
        ];
        for (i, ep) in self.grand_epochs.iter().enumerate() {
            let y0 = body_top + i as f32 * row_h;
            if i % 2 == 1 {
                tp.rect_filled(
                    egui::Rect::from_min_max(
                        egui::pos2(table_rect.left() + 6.0, y0),
                        egui::pos2(table_rect.right() - 6.0, y0 + row_h - 1.0),
                    ),
                    2.0,
                    egui::Color32::from_rgba_premultiplied(30, 41, 70, 80),
                );
            }
            let yc = y0 + row_h * 0.5;
            let ep_lbl = short_epoch_labels
                .get(i)
                .copied()
                .unwrap_or_else(|| ep.epoch.code());
            let mut landmark: String = ep.landmark_summary.chars().take(34).collect();
            if ep.landmark_summary.chars().count() > 34 {
                landmark.push_str("...");
            }
            let cells = [
                ep_lbl.to_string(),
                format!("{:.3}", ep.geologic_age_ga),
                format!("{}/5", ep.heritable_channels_active),
                format!("{}", ep.population_size),
                format!("{:.1}", ep.heritable_mutual_info_bits),
                format!("{:.2}", ep.informational_leverage_li),
                format!("{:+.3}", ep.order_velocity_vo),
                format!("{:.3}", ep.planetary_o2),
                format!("{:.0}", ep.planetary_co2_ppm),
                format!("{:.2}", ep.planetary_exergy_budget),
                landmark,
            ];
            for (c_idx, (x, cell)) in x_cols.iter().zip(cells.iter()).enumerate() {
                let col = match c_idx {
                    0 => egui::Color32::from_rgb(226, 232, 240),
                    2 => egui::Color32::from_rgb(253, 224, 71),
                    4 => egui::Color32::from_rgb(130, 245, 150),
                    7 => egui::Color32::from_rgb(76, 201, 240),
                    10 => egui::Color32::from_rgb(203, 213, 225),
                    _ => egui::Color32::from_rgb(186, 198, 218),
                };
                tp.text(
                    egui::pos2(*x, yc),
                    egui::Align2::LEFT_CENTER,
                    cell,
                    egui::FontId::monospace(8.8),
                    col,
                );
            }
        }
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
    let painter = &painter.with_clip_rect(rect);
    let layout = world.current_phylogram_layout();
    let p_sum = world.phylo_summary();

    let (net_title, net_legend) = match world.organism {
        OrganismKind::RnaDnaProtein => (
            format!(
                "Network — N{} | Hy:{} Mit:{} MRO0:{} | #H:{}",
                p_sum.mrca_id,
                p_sum.hydrogenosome_transition_events,
                p_sum.mitosome_transition_events,
                p_sum.mro_complete_loss_events,
                p_sum.endosymbiotic_fusions
            ),
            "Teal=Hydro(H2) | Gold=MitoS | Cyan=MRO0 | Violet=Aero",
        ),
        OrganismKind::HominidPopulation => (
            format!(
                "Cultural & Genealogical Lineage — MRCA:N{} | Deme Recomb:{} | Taxa:{}",
                p_sum.mrca_id,
                world.total_cross_deme_recombinations,
                world.last_species_count
            ),
            "Cyan=Forager | Green=Guild | Gold=Archive | Violet=Synergy",
        ),
        OrganismKind::ArtificialIntelligence => (
            format!(
                "Model & Distillation Graph — Root:N{} | Proofs:{} | Lineages:{}",
                p_sum.mrca_id,
                world.total_formal_proofs_verified,
                world.population.len()
            ),
            "Cyan=BaseWorker | Green=Verified | Gold=Distilled | Violet=Enclave",
        ),
    };
    painter.text(
        rect.left_top() + egui::vec2(8.0, 4.0),
        egui::Align2::LEFT_TOP,
        net_title,
        egui::FontId::proportional(12.5),
        egui::Color32::from_rgb(230, 240, 255),
    );
    painter.text(
        rect.left_top() + egui::vec2(8.0, 21.0),
        egui::Align2::LEFT_TOP,
        net_legend,
        egui::FontId::monospace(9.0),
        egui::Color32::from_rgb(142, 168, 212),
    );

    let plot = egui::Rect::from_min_max(
        rect.left_top() + egui::vec2(18.0, 44.0),
        rect.right_bottom() - egui::vec2(if show_labels { 148.0 } else { 18.0 }, 24.0),
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
            egui::Stroke::new(0.8_f32, egui::Color32::from_rgba_unmultiplied(55, 76, 120, 80)),
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
            1.9_f32
        };

        painter.line_segment([p_parent, p_corner], egui::Stroke::new(sw, col));
        painter.line_segment([p_corner, p_child], egui::Stroke::new(sw, col));

        let mid = egui::pos2((p_corner.x + p_child.x) * 0.5, p_child.y);
        if br.is_hydrogenosome_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(6, 214, 160));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(255, 209, 102)),
            );
        } else if br.is_mitosome_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(255, 209, 102));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(255, 92, 138)),
            );
        } else if br.is_mro_complete_loss {
            painter.circle_filled(mid, 5.6, egui::Color32::from_rgb(56, 189, 248));
            painter.circle_stroke(
                mid,
                5.6,
                egui::Stroke::new(1.7_f32, egui::Color32::from_rgb(255, 209, 102)),
            );
        } else if br.is_apicoplast_transition {
            painter.circle_filled(mid, 5.8, egui::Color32::from_rgb(255, 92, 138));
            painter.circle_stroke(
                mid,
                5.8,
                egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(199, 125, 255)),
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
                    painter.line_segment([prev, next_p], egui::Stroke::new(2.2_f32, arc_col));
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
                    egui::FontId::monospace(8.5),
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
    let painter = &painter.with_clip_rect(rect);
    let m = world.latest_metrics();
    let (sp_title, sp_legend) = match world.organism {
        OrganismKind::RnaDnaProtein => (
            format!(
                "Spatial — Aero:{} Hy:{} Mit:{} MRO0:{}",
                m.aerobic_mito_count, m.hydrogenosome_count, m.mitosome_count, m.amiteukaryote_count
            ),
            "Teal=Hydro | Gold=MitoS(ISC) | Cyan=MRO0(SUF)",
        ),
        OrganismKind::HominidPopulation => (
            format!(
                "Territorial Catchments — Pop:{} | Cult:{:.1}b | F_ST:{:.2}",
                m.population_size, m.mean_cultural_repertoire_bits, m.fst_fixation_index
            ),
            "Teal=Hearth | Gold=Archive | Green=Cultivar | Blue=Band",
        ),
        OrganismKind::ArtificialIntelligence => (
            format!(
                "Compute Cluster Grid — Agents:{} | Proofs:{} | TailDiv:{:.2}",
                m.population_size, m.formal_proofs_verified, m.mean_tail_diversity
            ),
            "Cyan=Worker | Green=Verifier | Violet=Enclave | Gold=Cached",
        ),
    };
    painter.text(
        rect.left_top() + egui::vec2(8.0, 4.0),
        egui::Align2::LEFT_TOP,
        sp_title,
        egui::FontId::proportional(12.5),
        egui::Color32::from_rgb(230, 240, 255),
    );
    painter.text(
        rect.left_top() + egui::vec2(8.0, 21.0),
        egui::Align2::LEFT_TOP,
        sp_legend,
        egui::FontId::monospace(9.0),
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
                    egui::Stroke::new(1.3_f32, egui::Color32::from_rgb(80, 250, 170)),
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
                        egui::Stroke::new(2.2_f32, egui::Color32::from_rgb(34, 211, 238)),
                    );
                } else if unit.chirality_code_stage != ChiralityCodeStage::UniversalErrorMinimizedSgc {
                    painter.circle_stroke(
                        center,
                        rad * 2.75,
                        egui::Stroke::new(2.2_f32, egui::Color32::from_rgb(192, 132, 252)),
                    );
                } else if unit.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage {
                    painter.circle_stroke(
                        center,
                        rad * 2.45,
                        egui::Stroke::new(2.1_f32, egui::Color32::from_rgb(251, 191, 36)),
                    );
                } else if unit.planetary_epoch == PlanetaryState::HomeostaticAwakening {
                    painter.circle_stroke(
                        center,
                        rad * 3.04,
                        egui::Stroke::new(2.6_f32, egui::Color32::from_rgb(52, 211, 153)),
                    );
                } else if unit.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis {
                    painter.circle_stroke(
                        center,
                        rad * 2.82,
                        egui::Stroke::new(2.6_f32, egui::Color32::from_rgb(139, 92, 246)),
                    );
                } else if unit.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                    painter.circle_stroke(
                        center,
                        rad * 2.58,
                        egui::Stroke::new(2.5_f32, egui::Color32::from_rgb(20, 184, 166)),
                    );
                } else if unit.cultural_state == CulturalState::SymbolicCumulativeCulture {
                    painter.circle_stroke(
                        center,
                        rad * 2.34,
                        egui::Stroke::new(2.4_f32, egui::Color32::from_rgb(250, 204, 21)),
                    );
                } else if unit.agriculture_state == AgricultureState::LeafCutterAgroEcosystem {
                    painter.circle_stroke(
                        center,
                        rad * 2.12,
                        egui::Stroke::new(2.4_f32, egui::Color32::from_rgb(163, 230, 53)),
                    );
                } else if unit.eusocial_state == EusocialState::ObligateSuperorganism {
                    painter.circle_stroke(
                        center,
                        rad * 1.88,
                        egui::Stroke::new(2.2_f32, egui::Color32::from_rgb(245, 158, 11)),
                    );
                } else if unit.multicell_state == MulticellularState::GermSomaMulticellular {
                    painter.circle_stroke(
                        center,
                        rad * 1.62,
                        egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(56, 189, 248)),
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
                    egui::Stroke::new(1.5_f32, rim),
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
            let champ_lbl = match world.organism {
                OrganismKind::RnaDnaProtein => {
                    format!("{} {}", unit.guild.short_code(), unit.mro_class.short_code())
                }
                OrganismKind::HominidPopulation => format!("Band {}", unit.species_id + 1),
                OrganismKind::ArtificialIntelligence => format!("Agent {}", unit.species_id + 1),
            };
            painter.text(
                center + egui::vec2(0.0, -rad - 3.5),
                egui::Align2::CENTER_BOTTOM,
                champ_lbl,
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
            painter.line_segment([p0, pt], egui::Stroke::new(1.6_f32, color));
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
            "ITME v23.0 — Domain-Isolated Biosphere, Hominid & AI Evolutionary Simulator",
        )
        .with_inner_size(egui::vec2(1600.0, 940.0))
        .with_min_inner_size(egui::vec2(980.0, 640.0))
        .with_resizable(true);

    let options = eframe::NativeOptions {
        viewport,
        ..Default::default()
    };

    eframe::run_native(
        "ITME v23.0 Native Rust Simulation",
        options,
        Box::new(move |cc| Ok(Box::new(ItmeGuiApp::new(cc, seed)))),
    )
    .map_err(|e| format!("Failed to launch native GUI: {e}"))
}

pub fn run_screenshot(out_path: &str, seed: u64, view_idx: u8) -> Result<(), String> {
    let ctx = egui::Context::default();
    let mut app = ItmeGuiApp::new_with_ctx(&ctx, seed);

    match view_idx {
        1 => {
            app.screen = AppScreen::MainMenu(MainMenuState {
                mode: SimulationMode::GuidedCurriculum,
                steps_per_stage: 32,
                initial_regime: RegimeId::R35Symbiogenetic,
                seed_text: seed.to_string(),
                seed,
            });
        }
        2 => {
            let mut session = DomainSession::new(
                OrganismKind::RnaDnaProtein,
                RegimeId::R20AutoadaptiveReplicative,
                SimulationMode::OrganicSandbox,
                32,
                seed,
            );
            for _ in 0..10 {
                session.step_once();
            }
            session.jump_to_regime(RegimeId::R30FullyDarwinian);
            for _ in 0..14 {
                session.step_once();
            }
            session.jump_to_regime(RegimeId::R35Symbiogenetic);
            session.world.controls.cytosolic_suf_hgt_rescue = 0.72;
            for _ in 0..2 {
                let _ = session.world.trigger_endosymbiotic_pulse();
            }
            for _ in 0..3 {
                let _ = session.world.trigger_mro_reduction_pulse();
            }
            let _ = session.world.trigger_serial_endosymbiosis_pulse();
            for _ in 0..20 {
                session.step_once();
            }
            session.viewport_layout = ViewportLayout::SpatialAndPhylogeny;
            session.paused = true;
            app.screen = AppScreen::Simulation(Box::new(session));
        }
        3 => {
            let mut session = DomainSession::new(
                OrganismKind::ArtificialIntelligence,
                RegimeId::R30SoftInheritance,
                SimulationMode::OrganicSandbox,
                32,
                seed + 202,
            );
            for _ in 0..12 {
                session.step_once();
            }
            session.jump_to_regime(RegimeId::R30PseudoSoftInheritance);
            let _ = session.world.trigger_machine_synthesis_pulse();
            for _ in 0..22 {
                session.step_once();
            }
            session.viewport_layout = ViewportLayout::SpatialGridAnd3D;
            session.camera.yaw = 0.54;
            session.camera.pitch = 0.30;
            session.camera.zoom = 1.25;
            session.paused = true;
            app.screen = AppScreen::Simulation(Box::new(session));
        }
        other => {
            return Err(format!(
                "Unsupported --screenshot-view {other} (expected 1 = Main Menu, 2 = Molecular Domain, 3 = AI Domain)"
            ));
        }
    }
    render_offscreen_png(&ctx, out_path, |c| app.draw_panels(c))
}

fn block_on_wgpu<F: std::future::Future>(fut: F) -> F::Output {
    struct ThreadWaker(std::thread::Thread);
    impl std::task::Wake for ThreadWaker {
        fn wake(self: std::sync::Arc<Self>) {
            self.0.unpark();
        }
        fn wake_by_ref(self: &std::sync::Arc<Self>) {
            self.0.unpark();
        }
    }
    let waker = std::task::Waker::from(std::sync::Arc::new(ThreadWaker(std::thread::current())));
    let mut fut = std::pin::pin!(fut);
    let mut cx = std::task::Context::from_waker(&waker);
    loop {
        match fut.as_mut().poll(&mut cx) {
            std::task::Poll::Ready(out) => return out,
            std::task::Poll::Pending => std::thread::park(),
        }
    }
}

fn render_offscreen_png(
    ctx: &egui::Context,
    out_path: &str,
    mut draw_ui: impl FnMut(&egui::Context),
) -> Result<(), String> {
    use eframe::{egui_wgpu, wgpu};
    const WIDTH: u32 = 1600;
    const HEIGHT: u32 = 1000;
    const BYTES_PER_ROW: u32 = WIDTH * 4;

    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::PRIMARY | wgpu::Backends::GL,
        ..Default::default()
    });
    let adapter = block_on_wgpu(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .or_else(|| {
        block_on_wgpu(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::None,
            compatible_surface: None,
            force_fallback_adapter: true,
        }))
    })
    .ok_or_else(|| "No wgpu adapter available for offscreen screenshot".to_string())?;

    let mut limits = wgpu::Limits::default();
    limits.max_texture_dimension_2d = adapter.limits().max_texture_dimension_2d.min(8192);
    let (device, queue) = block_on_wgpu(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("itme_offscreen_device"),
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: Default::default(),
        },
        None,
    ))
    .map_err(|e| format!("wgpu request_device failed: {e}"))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("itme_offscreen_target"),
        size: wgpu::Extent3d {
            width: WIDTH,
            height: HEIGHT,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("itme_offscreen_staging"),
        size: u64::from(BYTES_PER_ROW) * u64::from(HEIGHT),
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let mut renderer = egui_wgpu::Renderer::new(&device, format, None, 1, false);
    let screen_desc = egui_wgpu::ScreenDescriptor {
        size_in_pixels: [WIDTH, HEIGHT],
        pixels_per_point: 1.0,
    };

    for frame_idx in 0..4usize {
        let mut raw_input = egui::RawInput::default();
        raw_input.screen_rect = Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(WIDTH as f32, HEIGHT as f32),
        ));
        raw_input
            .viewports
            .entry(egui::ViewportId::ROOT)
            .and_modify(|v| v.native_pixels_per_point = Some(1.0));
        raw_input.time = Some(frame_idx as f64 / 60.0);
        ctx.set_pixels_per_point(1.0);

        let full_output = ctx.run(raw_input, &mut draw_ui);
        let paint_jobs = ctx.tessellate(full_output.shapes, full_output.pixels_per_point);

        for (id, delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, delta);
        }

        let mut prep_encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("itme_offscreen_prep"),
        });
        let user_cmds =
            renderer.update_buffers(&device, &queue, &mut prep_encoder, &paint_jobs, &screen_desc);
        queue.submit(user_cmds.into_iter().chain(std::iter::once(prep_encoder.finish())));

        let mut pass_encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("itme_offscreen_pass_encoder"),
        });
        {
            let mut rpass = pass_encoder
                .begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("itme_offscreen_pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color {
                                r: 0.05,
                                g: 0.07,
                                b: 0.12,
                                a: 1.0,
                            }),
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    timestamp_writes: None,
                    occlusion_query_set: None,
                })
                .forget_lifetime();
            renderer.render(&mut rpass, &paint_jobs, &screen_desc);
        }

        if frame_idx == 3 {
            pass_encoder.copy_texture_to_buffer(
                wgpu::ImageCopyTexture {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::ImageCopyBuffer {
                    buffer: &staging,
                    layout: wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(BYTES_PER_ROW),
                        rows_per_image: Some(HEIGHT),
                    },
                },
                wgpu::Extent3d {
                    width: WIDTH,
                    height: HEIGHT,
                    depth_or_array_layers: 1,
                },
            );
        }

        queue.submit(std::iter::once(pass_encoder.finish()));
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
        device.poll(wgpu::Maintain::Wait);
    }

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::sync_channel(1);
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .map_err(|e| format!("GPU readback channel failed: {e}"))?
        .map_err(|e| format!("GPU buffer map failed: {e}"))?;

    let rgba = slice.get_mapped_range().to_vec();
    staging.unmap();
    write_rgba_png(std::path::Path::new(out_path), WIDTH, HEIGHT, &rgba)
}

fn write_rgba_png(path: &std::path::Path, width: u32, height: u32, rgba: &[u8]) -> Result<(), String> {
    if rgba.len() != (width as usize) * (height as usize) * 4 {
        return Err("RGBA buffer length mismatch".to_string());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|e| format!("mkdir failed: {e}"))?;
        }
    }
    fn crc32(bytes: &[u8]) -> u32 {
        let mut c = !0u32;
        for &b in bytes {
            c ^= u32::from(b);
            for _ in 0..8 {
                c = (c >> 1) ^ (0xEDB8_8320 & (0u32.wrapping_sub(c & 1)));
            }
        }
        !c
    }
    fn push_chunk(out: &mut Vec<u8>, tag: &[u8; 4], data: &[u8]) {
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        let start = out.len();
        out.extend_from_slice(tag);
        out.extend_from_slice(data);
        out.extend_from_slice(&crc32(&out[start..]).to_be_bytes());
    }
    let row_bytes = (width as usize) * 4;
    let mut raw = Vec::with_capacity((row_bytes + 1) * (height as usize));
    for row in rgba.chunks_exact(row_bytes) {
        raw.push(0);
        raw.extend_from_slice(row);
    }
    let (mut s1, mut s2) = (1u32, 0u32);
    for &b in &raw {
        s1 = (s1 + u32::from(b)) % 65521;
        s2 = (s2 + s1) % 65521;
    }
    let mut zlib = Vec::with_capacity(raw.len() + 64);
    zlib.extend_from_slice(&[0x78, 0x01]);
    let mut off = 0usize;
    while off < raw.len() {
        let end = (off + u16::MAX as usize).min(raw.len());
        let block = &raw[off..end];
        let len = block.len() as u16;
        zlib.push(u8::from(end == raw.len()));
        zlib.extend_from_slice(&len.to_le_bytes());
        zlib.extend_from_slice(&(!len).to_le_bytes());
        zlib.extend_from_slice(block);
        off = end;
    }
    zlib.extend_from_slice(&((s2 << 16) | s1).to_be_bytes());

    let mut png = Vec::with_capacity(zlib.len() + 64);
    png.extend_from_slice(&[137, 80, 78, 71, 13, 10, 26, 10]);
    let mut ihdr = [0u8; 13];
    ihdr[0..4].copy_from_slice(&width.to_be_bytes());
    ihdr[4..8].copy_from_slice(&height.to_be_bytes());
    ihdr[8] = 8;
    ihdr[9] = 6;
    push_chunk(&mut png, b"IHDR", &ihdr);
    push_chunk(&mut png, b"IDAT", &zlib);
    push_chunk(&mut png, b"IEND", &[]);
    fs::write(path, &png).map_err(|e| format!("Failed writing PNG '{}': {e}", path.display()))
}
