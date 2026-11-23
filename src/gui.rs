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
//! Starts at a 4-Card Main Menu launcher (3 isolated domain sandboxes: `RNA/DNA/Protein`,
//! `Hominid Populations`, and `Artificial Intelligences`, plus the coupled `4.25 Ga PlanetaryBiosphere`
//! engine), where launching any card opens a dedicated `SimulationWorld` with tailored controls and telemetry.

use crate::engine::{SimulationMode, SimulationWorld};
use crate::geometry::{Camera3D, StructureScene};
use crate::model::{
    AbiogenesisStage, AgricultureState, ChiralityCodeStage, CulturalState, EusocialState,
    ExosomaticState, GenomeSubstrateStage, GrandTrajectoryEpochSummary, MachineSynthesisState,
    MroClass, MulticellularState, OrderPhase, PlanetaryBiosphere, PlanetaryState, TrophicGuild,
};
use crate::organisms::OrganismKind;
use crate::regimes::RegimeId;
use eframe::{egui, App};
use std::fs;
use std::io::BufWriter;
use std::sync::Arc;
use std::time::{Duration, Instant};

fn add_param_slider(
    ui: &mut egui::Ui,
    value: &mut f64,
    range: std::ops::RangeInclusive<f64>,
    label: &str,
    tooltip: &str,
) -> bool {
    ui.add(egui::Slider::new(value, range).text(label))
        .on_hover_text(tooltip)
        .changed()
}

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

    fn try_parse_seed(&self) -> Result<u64, &'static str> {
        let trimmed = self.seed_text.trim();
        if trimmed.is_empty() {
            return Err("Seed cannot be empty (enter decimal or 0x-prefixed hex u64).");
        }
        if let Some(hex) = trimmed.strip_prefix("0x").or_else(|| trimmed.strip_prefix("0X")) {
            if hex.is_empty() {
                return Err("Hex seed after 0x cannot be empty.");
            }
            u64::from_str_radix(hex, 16)
                .map_err(|_| "Invalid hexadecimal seed (expected 0x followed by hex digits).")
        } else {
            trimmed
                .parse::<u64>()
                .map_err(|_| "Invalid decimal seed (expected unsigned integer or 0x hex).")
        }
    }

    fn parsed_seed(&self) -> u64 {
        self.try_parse_seed().unwrap_or(self.seed)
    }
}

struct DomainSession {
    organism: OrganismKind,
    coupled_biosphere: bool,
    world: SimulationWorld,
    viewport_layout: ViewportLayout,
    camera: Camera3D,
    auto_rotate: bool,
    show_labels: bool,
    show_milestones: bool,
    paused: bool,
    steps_per_second: f32,
    last_step_instant: Instant,
    export_status: Option<(String, bool)>,
    grand_epochs: Vec<GrandTrajectoryEpochSummary>,
    wants_main_menu: bool,
    wants_about: bool,
}

enum AppScreen {
    MainMenu(MainMenuState),
    Simulation(Box<DomainSession>),
}

pub struct ItmeGuiApp {
    screen: AppScreen,
    default_seed: u64,
    show_about: bool,
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
            show_about: false,
        }
    }

    fn draw_panels(&mut self, ctx: &egui::Context) {
        let esc_pressed = ctx.input(|i| i.key_pressed(egui::Key::Escape));
        if self.show_about && esc_pressed {
            self.show_about = false;
        } else if !ctx.wants_keyboard_input() {
            if ctx.input(|i| i.key_pressed(egui::Key::F1)) {
                self.show_about = !self.show_about;
            } else if !self.show_about {
                if let AppScreen::Simulation(session) = &mut self.screen {
                    session.handle_global_hotkeys(ctx, esc_pressed);
                }
            }
        }

        let mut next_screen: Option<AppScreen> = None;
        match &mut self.screen {
            AppScreen::MainMenu(menu) => {
                if let Some(session) = Self::draw_main_menu(ctx, menu, &mut self.show_about) {
                    next_screen = Some(AppScreen::Simulation(Box::new(session)));
                }
            }
            AppScreen::Simulation(session) => {
                session.draw_panels(ctx);
                if session.wants_about {
                    session.wants_about = false;
                    self.show_about = !self.show_about;
                }
                if session.wants_main_menu {
                    let seed = self.default_seed;
                    next_screen = Some(AppScreen::MainMenu(MainMenuState::new(seed)));
                }
            }
        }
        if let Some(scr) = next_screen {
            self.screen = scr;
        }
        self.draw_about_window(ctx);
    }

    fn draw_about_window(&mut self, ctx: &egui::Context) {
        if !self.show_about {
            return;
        }
        egui::Window::new("About ITME")
            .open(&mut self.show_about)
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                ui.set_min_width(380.0);
                ui.heading(
                    egui::RichText::new("ITME — Information-Theoretic Model of Evolution")
                        .color(egui::Color32::from_rgb(130, 245, 205))
                        .strong(),
                );
                ui.label(
                    egui::RichText::new(concat!("Version ", env!("CARGO_PKG_VERSION"), " (ITME v23.0)"))
                        .strong(),
                );
                ui.add_space(6.0);
                ui.label("Copyright CRV Science Network");
                ui.label("License: GNU General Public License v3.0 or later (GPL-3.0-or-later)");
                ui.add_space(6.0);
                ui.add(
                    egui::Label::new(
                        egui::RichText::new(
                            "Continuous 4.25 Ga multi-epoch evolutionary, symbiogenetic, biocultural, formal AI, and planetary biosphere simulator.",
                        )
                        .small()
                        .color(egui::Color32::from_rgb(186, 198, 218)),
                    )
                    .wrap(),
                );
                ui.add_space(8.0);
                ui.separator();
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.hyperlink_to("CRV Science Network", "https://crvscience.com");
                    ui.separator();
                    ui.hyperlink_to("CRV Science Apps Hub", "https://crvscience.com/apps");
                });
            });
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

    fn biosphere_card_meta() -> (&'static str, &'static str, &'static str, egui::Color32, &'static [&'static str], &'static str) {
        (
            "Domain 4 • Coupled Planetary Biosphere",
            "4.25 Ga PlanetaryBiosphere",
            "Hadean-to-Gaia-2.0 Biogeochemical & 5-Channel Coupling",
            egui::Color32::from_rgb(52, 211, 153),
            &[
                "Hadean anoxic FeNiS cradle -> GOE oxygenation (O2 >= 0.25)",
                "Live O2, CO2 ppm, N/P pool & Gibbs disequilibrium feedbacks",
                "Continuous 8-Epoch 4.25 Ga -> Gaia-2.0 Grand Trajectory",
                "Extractive CO2 overshoot vs verified planetary stewardship",
                "Coupled biological, cultural, archival & AI channels (5/5)",
            ],
            "Couples autotrophic, syntrophic, superorganismal, and machine-verified guilds to live atmospheric O2, CO2, N/P, and Gibbs exergy reservoirs.",
        )
    }

    fn draw_main_menu(
        ctx: &egui::Context,
        menu: &mut MainMenuState,
        show_about: &mut bool,
    ) -> Option<DomainSession> {
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
                    if ui
                        .button("About")
                        .on_hover_text("About ITME & CRV Science Network (F1)")
                        .clicked()
                    {
                        *show_about = !*show_about;
                    }
                    ui.separator();
                    ui.label(
                        egui::RichText::new("Copyright CRV Science Network")
                            .small()
                            .color(egui::Color32::from_rgb(168, 184, 208)),
                    );
                });
            });
            ui.label(
                egui::RichText::new(
                    "Select an isolated single-domain sandbox (Domains 1-3) or launch the coupled 4.25 Ga PlanetaryBiosphere engine (Domain 4) with live biogeochemical feedbacks and the 8-Epoch Grand Trajectory.",
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
                let seed_res = menu.try_parse_seed();
                let seed_resp = ui.add(
                    egui::TextEdit::singleline(&mut menu.seed_text).desired_width(90.0),
                );
                if let Err(err_hint) = seed_res {
                    seed_resp.on_hover_text(format!("{err_hint} Fallback seed: {}.", menu.seed));
                    ui.label(
                        egui::RichText::new(format!("Invalid seed (using {})", menu.seed))
                            .small()
                            .color(egui::Color32::from_rgb(248, 113, 113)),
                    )
                    .on_hover_text(err_hint);
                } else {
                    seed_resp.on_hover_text(
                        "Deterministic Xoshiro256** PRNG seed (decimal or 0x-prefixed hex).",
                    );
                }
                if menu.mode == SimulationMode::GuidedCurriculum {
                    ui.separator();
                    ui.add(egui::Slider::new(&mut menu.steps_per_stage, 12..=80).text("Ticks/Regime"));
                }
            });
            ui.add_space(4.0);
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.add_space(8.0);
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.columns(4, |cols| {
                    for (idx, &org) in OrganismKind::ALL.iter().enumerate() {
                        let (badge, subtitle, accent, highlights, dynamics_note) = Self::domain_card_meta(org);
                        let col_ui = &mut cols[idx];
                        let inner_w = (col_ui.available_width() - 26.0).max(170.0);
                        egui::Frame::none()
                            .fill(egui::Color32::from_rgb(19, 27, 46))
                            .stroke(egui::Stroke::new(1.5_f32, accent))
                            .rounding(10.0)
                            .inner_margin(egui::Margin::same(12.0))
                            .show(col_ui, |card| {
                                card.set_min_width(inner_w);
                                card.set_max_width(inner_w);
                                card.set_min_height(430.0);

                                card.label(
                                    egui::RichText::new(badge)
                                        .small()
                                        .strong()
                                        .color(accent),
                                );
                                card.add_space(2.0);
                                card.heading(
                                    egui::RichText::new(org.title())
                                        .size(16.5)
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
                                card.add_space(6.0);
                                card.separator();
                                card.add_space(4.0);

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
                                card.add_space(6.0);

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

                                card.add_space(6.0);
                                card.strong(egui::RichText::new("Isolated Sandbox").color(accent));
                                card.add(
                                    egui::Label::new(
                                        egui::RichText::new(dynamics_note)
                                            .small()
                                            .color(egui::Color32::from_rgb(162, 178, 204)),
                                    )
                                    .wrap(),
                                );

                                card.add_space(10.0);
                                card.strong(
                                    egui::RichText::new(format!(
                                        "Preview {}: {}",
                                        menu.initial_regime.code(),
                                        menu.initial_regime.domain_name(org)
                                    ))
                                    .small()
                                    .color(egui::Color32::from_rgb(226, 232, 240)),
                                );
                                card.add_space(6.0);

                                let launch_btn = egui::Button::new(
                                    egui::RichText::new(format!("Launch {} ->", org.short_title()))
                                        .size(13.5)
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

                    // Card 4: Coupled 4.25 Ga PlanetaryBiosphere
                    let (badge, title, subtitle, accent, highlights, dynamics_note) =
                        Self::biosphere_card_meta();
                    let col_ui = &mut cols[3];
                    let inner_w = (col_ui.available_width() - 26.0).max(170.0);
                    egui::Frame::none()
                        .fill(egui::Color32::from_rgb(16, 29, 44))
                        .stroke(egui::Stroke::new(1.5_f32, accent))
                        .rounding(10.0)
                        .inner_margin(egui::Margin::same(12.0))
                        .show(col_ui, |card| {
                            card.set_min_width(inner_w);
                            card.set_max_width(inner_w);
                            card.set_min_height(430.0);

                            card.label(
                                egui::RichText::new(badge)
                                    .small()
                                    .strong()
                                    .color(accent),
                            );
                            card.add_space(2.0);
                            card.heading(
                                egui::RichText::new(title)
                                    .size(16.5)
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
                            card.add_space(6.0);
                            card.separator();
                            card.add_space(4.0);

                            card.strong(egui::RichText::new("Planetary Reservoirs").color(accent));
                            card.add(
                                egui::Label::new(
                                    egui::RichText::new("O2 | CO2 ppm | N/P | Gibbs Bits | Exergy")
                                        .small()
                                        .monospace()
                                        .color(egui::Color32::from_rgb(203, 213, 225)),
                                )
                                .wrap(),
                            );
                            card.add_space(6.0);

                            card.strong(egui::RichText::new("Coupled Trajectory & Gaia-2.0").color(accent));
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

                            card.add_space(6.0);
                            card.strong(egui::RichText::new("Coupled Biogeochemistry").color(accent));
                            card.add(
                                egui::Label::new(
                                    egui::RichText::new(dynamics_note)
                                        .small()
                                        .color(egui::Color32::from_rgb(162, 178, 204)),
                                )
                                .wrap(),
                            );

                            card.add_space(10.0);
                            card.strong(
                                egui::RichText::new(format!(
                                    "Preview {}: {}",
                                    menu.initial_regime.code(),
                                    menu.initial_regime.name()
                                ))
                                .small()
                                .color(egui::Color32::from_rgb(226, 232, 240)),
                            );
                            card.add_space(6.0);

                            let launch_btn = egui::Button::new(
                                egui::RichText::new("Launch 4.25 Ga Biosphere ->")
                                    .size(13.5)
                                    .strong()
                                    .color(egui::Color32::from_rgb(10, 15, 28)),
                            )
                            .fill(accent);
                            if card.add_sized([inner_w, 34.0], launch_btn).clicked() {
                                let seed = menu.parsed_seed();
                                launched = Some(DomainSession::new_coupled_biosphere(
                                    menu.initial_regime,
                                    menu.mode,
                                    menu.steps_per_stage,
                                    seed,
                                ));
                            }
                        });
                });
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
        world.coupled_biosphere_mode = false;
        world.controls = initial_regime.domain_preset_controls(organism);
        world.step();

        Self {
            organism,
            coupled_biosphere: false,
            world,
            viewport_layout: ViewportLayout::SpatialAndPhylogeny,
            camera: Camera3D::default(),
            auto_rotate: true,
            show_labels: true,
            show_milestones: true,
            paused: false,
            steps_per_second: 8.0,
            last_step_instant: Instant::now(),
            export_status: None,
            grand_epochs: Vec::new(),
            wants_main_menu: false,
            wants_about: false,
        }
    }

    fn new_coupled_biosphere(
        initial_regime: RegimeId,
        mode: SimulationMode,
        steps_per_stage: usize,
        seed: u64,
    ) -> Self {
        let mut world = SimulationWorld::new(
            OrganismKind::ArtificialIntelligence,
            initial_regime,
            mode,
            steps_per_stage,
            seed,
        );
        world.coupled_biosphere_mode = true;
        world.controls = initial_regime.preset_controls();
        if initial_regime >= RegimeId::R30FullyDarwinian {
            world.controls.biogeochemical_feedback_gain =
                world.controls.biogeochemical_feedback_gain.max(0.65);
            world.controls.planetary_stewardship_coupling =
                world.controls.planetary_stewardship_coupling.max(0.75);
        }
        world.biosphere = PlanetaryBiosphere::new_for_controls(&world.controls);
        world.step();

        Self {
            organism: OrganismKind::ArtificialIntelligence,
            coupled_biosphere: true,
            world,
            viewport_layout: ViewportLayout::SpatialAndPhylogeny,
            camera: Camera3D::default(),
            auto_rotate: true,
            show_labels: true,
            show_milestones: true,
            paused: false,
            steps_per_second: 8.0,
            last_step_instant: Instant::now(),
            export_status: None,
            grand_epochs: Vec::new(),
            wants_main_menu: false,
            wants_about: false,
        }
    }

    fn step_once(&mut self) {
        let prev_reg = self.world.target_regime;
        self.world.step();
        if self.world.mode == SimulationMode::GuidedCurriculum && self.world.target_regime != prev_reg {
            if self.coupled_biosphere {
                self.world.controls = self.world.target_regime.preset_controls();
                if self.world.target_regime >= RegimeId::R30FullyDarwinian {
                    self.world.controls.biogeochemical_feedback_gain =
                        self.world.controls.biogeochemical_feedback_gain.max(0.65);
                }
            } else {
                self.world.controls = self.world.target_regime.domain_preset_controls(self.organism);
            }
        }
    }

    fn jump_to_regime(&mut self, regime: RegimeId) {
        let preserve = regime.index() >= 2;
        self.world.apply_regime_preset(regime, preserve);
        if self.coupled_biosphere {
            self.world.controls = regime.preset_controls();
            if regime >= RegimeId::R30FullyDarwinian {
                self.world.controls.biogeochemical_feedback_gain =
                    self.world.controls.biogeochemical_feedback_gain.max(0.65);
                self.world.controls.planetary_stewardship_coupling =
                    self.world.controls.planetary_stewardship_coupling.max(0.75);
            }
        } else {
            self.world.controls = regime.domain_preset_controls(self.organism);
        }
        self.world.step();
    }

    fn save_text_with_dialog(
        &mut self,
        default_name: &str,
        filter_name: &str,
        extension: &str,
        contents: &str,
    ) {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(default_name)
            .add_filter(filter_name, &[extension])
            .save_file()
        {
            match fs::write(&path, contents) {
                Ok(()) => {
                    self.export_status = Some((format!("Saved {}", path.display()), false));
                }
                Err(e) => {
                    self.export_status =
                        Some((format!("Export error ({}): {e}", path.display()), true));
                }
            }
        } else {
            self.export_status = Some(("Export cancelled.".to_string(), false));
        }
    }

    fn handle_global_hotkeys(&mut self, ctx: &egui::Context, esc_pressed: bool) {
        if esc_pressed {
            self.wants_main_menu = true;
            return;
        }
        ctx.input(|i| {
            if i.key_pressed(egui::Key::Space) {
                self.paused = !self.paused;
                self.last_step_instant = Instant::now();
            }
            if i.key_pressed(egui::Key::Period) {
                self.step_once();
            }
            if i.key_pressed(egui::Key::R) {
                self.camera = Camera3D::default();
            }
            if i.key_pressed(egui::Key::Num1) {
                self.viewport_layout = ViewportLayout::SpatialAndPhylogeny;
            } else if i.key_pressed(egui::Key::Num2) {
                self.viewport_layout = ViewportLayout::SpatialGridAnd3D;
            } else if i.key_pressed(egui::Key::Num3) {
                self.viewport_layout = ViewportLayout::PhylogenyFull;
            } else if i.key_pressed(egui::Key::Num4) {
                self.viewport_layout = ViewportLayout::Perspective3D;
            } else if i.key_pressed(egui::Key::Num5) {
                self.viewport_layout = ViewportLayout::Blueprint2D;
            } else if i.key_pressed(egui::Key::Num6) && self.coupled_biosphere {
                self.viewport_layout = ViewportLayout::GrandTrajectory;
            }
        });
    }

    fn step_to_ga(&self, step: usize) -> f64 {
        if !self.grand_epochs.is_empty() {
            let first = &self.grand_epochs[0];
            if step <= first.tick_end {
                return first.geologic_age_ga;
            }
            for win in self.grand_epochs.windows(2) {
                let a = &win[0];
                let b = &win[1];
                if step <= b.tick_end {
                    let span = (b.tick_end.saturating_sub(a.tick_end)).max(1) as f64;
                    let frac = ((step.saturating_sub(a.tick_end)) as f64 / span).clamp(0.0, 1.0);
                    return a.geologic_age_ga + (b.geologic_age_ga - a.geologic_age_ga) * frac;
                }
            }
            return self.grand_epochs.last().map_or(0.0, |e| e.geologic_age_ga);
        }
        let total_ticks = if self.world.mode == SimulationMode::GuidedCurriculum {
            (self.world.steps_per_stage.max(12) * 13) as f64
        } else {
            240.0
        };
        (4.25 * (1.0 - (step as f64) / total_ticks)).clamp(0.0, 4.25)
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
                    .on_hover_text("Return to the Domain Launcher Main Menu (Esc)")
                    .clicked()
                {
                    self.wants_main_menu = true;
                }
                ui.separator();

                let (title_text, accent) = if self.coupled_biosphere {
                    (
                        "ITME v23.0 • 4.25 Ga PlanetaryBiosphere".to_string(),
                        egui::Color32::from_rgb(52, 211, 153),
                    )
                } else {
                    let (_, _, accent, _, _) = ItmeGuiApp::domain_card_meta(self.organism);
                    (format!("ITME v23.0 • {}", self.organism.title()), accent)
                };
                ui.label(egui::RichText::new(title_text).color(accent).strong());
                ui.separator();

                let play_label = if self.paused { "Play" } else { "Pause" };
                if ui
                    .button(play_label)
                    .on_hover_text("Toggle simulation playback (Space)")
                    .clicked()
                {
                    self.paused = !self.paused;
                    self.last_step_instant = Instant::now();
                }
                if ui
                    .button("Step")
                    .on_hover_text("Advance simulation by one discrete tick (.)")
                    .clicked()
                {
                    self.step_once();
                }
                let slug = if self.coupled_biosphere {
                    "planetary_biosphere"
                } else {
                    self.organism.slug()
                };
                let step = self.world.step_count;
                let export_resp = ui.menu_button("Export...", |ui| {
                    if ui
                        .button("Telemetry CSV (.csv)")
                        .on_hover_text("Export step-by-step telemetry (or 8-Epoch Grand Trajectory) to CSV.")
                        .clicked()
                    {
                        let csv_data = if self.coupled_biosphere
                            && !self.grand_epochs.is_empty()
                            && self.viewport_layout == ViewportLayout::GrandTrajectory
                        {
                            SimulationWorld::export_grand_trajectory_csv(&self.grand_epochs)
                        } else {
                            self.world.export_history_csv()
                        };
                        let fname = format!("itme_{slug}_step{step}.csv");
                        self.save_text_with_dialog(&fname, "CSV Telemetry", "csv", &csv_data);
                        ui.close_menu();
                    }
                    if ui
                        .button("Stage Atlas SVG (.svg)")
                        .on_hover_text("Export 13-regime 3D/2D structural stage atlas as vector SVG.")
                        .clicked()
                    {
                        let atlas_svg = self.world.build_stage_atlas_svg();
                        let fname = format!("itme_{slug}_stage_atlas.svg");
                        self.save_text_with_dialog(&fname, "SVG Vector Graphics", "svg", &atlas_svg);
                        ui.close_menu();
                    }
                    if ui
                        .button("Extended Newick (.nhx)")
                        .on_hover_text("Export Extended Newick (eNewick #H + MRO/APICO/KLEPTO NHX) lineage network.")
                        .clicked()
                    {
                        let nhx = self.world.export_newick(true);
                        let fname = format!("itme_{slug}_step{step}.nhx");
                        self.save_text_with_dialog(&fname, "Extended Newick", "nhx", &nhx);
                        ui.close_menu();
                    }
                    if ui
                        .button("Phylogram SVG (.svg)")
                        .on_hover_text("Export reticulate anastomosing evolutionary phylogram as vector SVG.")
                        .clicked()
                    {
                        let phylo_svg = self.world.export_phylogram_svg();
                        let fname = format!("itme_{slug}_step{step}_phylogram.svg");
                        self.save_text_with_dialog(&fname, "SVG Phylogram", "svg", &phylo_svg);
                        ui.close_menu();
                    }
                });
                export_resp
                    .response
                    .on_hover_text("Export telemetry CSV, 13-regime Stage Atlas SVG, Extended Newick (NHX), or Phylogram SVG.");
                ui.add(
                    egui::Slider::new(&mut self.steps_per_second, 1.0..=30.0)
                        .text("ticks/s")
                        .logarithmic(false),
                )
                .on_hover_text("Simulation playback speed in discrete ticks per second.");

                ui.separator();
                ui.label("View:");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialAndPhylogeny,
                    "Spatial + Network",
                )
                .on_hover_text("Split view: 2D spatial reaction-diffusion grid + reticulate phylogram (1)");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::SpatialGridAnd3D,
                    "Spatial + 3D",
                )
                .on_hover_text("Split view: 2D spatial grid + interactive 3D structure (2)");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::PhylogenyFull,
                    "Full Lineage",
                )
                .on_hover_text("Full-width anastomosing evolutionary lineage tree (3)");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::Perspective3D,
                    "3D Structure",
                )
                .on_hover_text("Interactive 3D structure — drag to orbit, scroll to zoom, R to reset camera (4)");
                ui.selectable_value(
                    &mut self.viewport_layout,
                    ViewportLayout::Blueprint2D,
                    "2D Blueprint",
                )
                .on_hover_text("2D schematic structural blueprint with directed edges (5)");
                if self.coupled_biosphere {
                    ui.selectable_value(
                        &mut self.viewport_layout,
                        ViewportLayout::GrandTrajectory,
                        "4.25 Ga Timeline",
                    )
                    .on_hover_text("Continuous 8-Epoch 4.25 Ga -> Gaia-2.0 Grand Trajectory plots & ledger (6)");
                }

                ui.separator();
                ui.checkbox(&mut self.auto_rotate, "Auto-Orbit")
                    .on_hover_text("Automatically orbit the 3D structural camera.");
                ui.checkbox(&mut self.show_labels, "Labels")
                    .on_hover_text("Show/hide pill labels across spatial, lineage, 3D, and 2D viewports.");
                let ms_lbl = if self.show_milestones {
                    "Milestones [-]"
                } else {
                    "Milestones [+]"
                };
                ui.toggle_value(&mut self.show_milestones, ms_lbl)
                    .on_hover_text("Show or collapse quick milestone trigger buttons to maximize viewport canvas height.");

                ui.separator();
                ui.label(
                    egui::RichText::new("Copyright CRV Science Network")
                        .small()
                        .color(egui::Color32::from_rgb(168, 184, 208)),
                );
                if ui
                    .button("About")
                    .on_hover_text("About ITME & CRV Science Network (F1)")
                    .clicked()
                {
                    self.wants_about = true;
                }
            });

            if self.show_milestones {
                ui.add_space(2.0);
                ui.horizontal_wrapped(|ui| {
                    ui.label(
                        egui::RichText::new(if self.coupled_biosphere {
                            "Biosphere Milestones:"
                        } else {
                            "Domain Milestones:"
                        })
                        .small()
                        .strong()
                        .color(egui::Color32::from_rgb(148, 163, 184)),
                    );
                    if self.coupled_biosphere {
                        if ui
                            .button("4.25 Ga Grand Trajectory")
                            .on_hover_text("Run unbroken 8-Epoch continuous simulation from 4.25 Ga to Gaia-2.0.")
                            .clicked()
                        {
                            let epochs = self.world.run_continuous_grand_trajectory(10);
                            self.grand_epochs = epochs.clone();
                            self.viewport_layout = ViewportLayout::GrandTrajectory;
                            if let Some(last_ep) = epochs.last() {
                                self.export_status = Some((
                                    format!(
                                        "8-Epoch Trajectory -> {}: L_I={:.2}, MI={:.1}b, O2 (PAL)={:.2}",
                                        last_ep.epoch.code(),
                                        last_ep.informational_leverage_li,
                                        last_ep.heritable_mutual_info_bits,
                                        last_ep.planetary_o2
                                    ),
                                    false,
                                ));
                            }
                        }
                        if ui
                            .button("Gaia-2.0 Awakening")
                            .on_hover_text("Activate verified planetary stewardship co-regulating O2 (PAL Index), CO2, and exergy.")
                            .clicked()
                        {
                            self.world.controls.biogeochemical_feedback_gain =
                                self.world.controls.biogeochemical_feedback_gain.max(0.85);
                            self.world.controls.planetary_stewardship_coupling =
                                self.world.controls.planetary_stewardship_coupling.max(0.92);
                            self.world.controls.extractive_burn_intensity =
                                self.world.controls.extractive_burn_intensity.min(0.18);
                            let (_, aw_n, ex_b) = self.world.trigger_planetary_awakening_pulse();
                            self.export_status = Some((
                                format!(
                                    "Gaia-2.0 Awakening: {aw_n} co-regulated units | O2 (PAL)={:.2}, CO2={:.0}ppm, Exergy={ex_b:.2}",
                                    self.world.biosphere.atmospheric_o2,
                                    self.world.biosphere.atmospheric_co2_ppm
                                ),
                                false,
                            ));
                        }
                        if ui
                            .button("Extractive CO2 Overshoot")
                            .on_hover_text("Simulate unguided high-burn extractive stress pushing CO2 and boundary transgression load.")
                            .clicked()
                        {
                            self.world.controls.biogeochemical_feedback_gain = 0.85;
                            self.world.controls.planetary_stewardship_coupling = 0.08;
                            self.world.controls.extractive_burn_intensity = 0.92;
                            for _ in 0..8 {
                                self.step_once();
                            }
                            self.export_status = Some((
                                format!(
                                    "Boundary Stress: State={}, CO2={:.0}ppm, Load={:.2}",
                                    self.world.biosphere.state.short_code(),
                                    self.world.biosphere.atmospheric_co2_ppm,
                                    self.world.biosphere.boundary_transgression_load
                                ),
                                false,
                            ));
                        }
                    } else {
                        match self.organism {
                            OrganismKind::RnaDnaProtein => {
                                if ui
                                    .button("Hadean Abiogenesis -> LUCA")
                                    .on_hover_text("Seed citrate-chelated RAF protocells and step toward Darwinian LUCA.")
                                    .clicked()
                                {
                                    let (rf, lc, eig) = self.world.trigger_hadean_abiogenesis_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Hadean Abiogenesis: +{rf} RAF Protocells, +{lc} LUCA (Eigen L_max={eig} nt)"
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Homochiral D/L + SGC")
                                    .on_hover_text("Trigger Fe3O4 CISS D-RNA symmetry breaking, L-aminoacylation, and 20-AA SGC crystallization.")
                                    .clicked()
                                {
                                    let (hm, sgc, mp) = self.world.trigger_homochiral_sgc_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Homochirality & SGC: +{hm} D-RNA breaks, +{sgc} 20-AA SGC (MP={:.1}%)",
                                            mp * 100.0
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("RNR RNA->T-DNA")
                                    .on_hover_text("Activate Class III RNR deoxyribose reduction and ThyX/UDG cytosine-deamination repair.")
                                    .clicked()
                                {
                                    let (rn, td, kb) = self.world.trigger_rna_to_dna_virosphere_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "RNR & UDG Repair: +{rn} U-DNA, +{td} T-DNA Genomes (Cap={kb:.0} kb)"
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("1°/2° Symbio + MRO")
                                    .on_hover_text("Trigger 1°/2° endosymbiotic fusion and anaerobic MRO streamlining (Mito->Hydro->MitoS->MRO->0).")
                                    .clicked()
                                {
                                    let fused = self.world.trigger_endosymbiotic_pulse();
                                    let ser = self.world.trigger_serial_endosymbiosis_pulse();
                                    let (hy, mi, am) = self.world.trigger_mro_reduction_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Symbiosis & MRO: 1°={fused}, Serial={ser:?}, +Hy:{hy} +MitoS:{mi} +SUF(MRO->0):{am}"
                                        ),
                                        false,
                                    ));
                                }
                            }
                            OrganismKind::HominidPopulation => {
                                if ui
                                    .button("Kin Division & Lithic Guilds")
                                    .on_hover_text("Promote foraging bands into kin-structured lithic guilds in R(3.0-DAR).")
                                    .clicked()
                                {
                                    self.jump_to_regime(RegimeId::R30FullyDarwinian);
                                    self.export_status = Some((
                                        "Promoted foraging bands to kin-structured lithic guilds R(3.0-DAR)."
                                            .to_string(),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Symbolic Teaching + Ratchet")
                                    .on_hover_text("Elevate symbolic pedagogical fidelity and partial deme connectivity to drive cultural ratchets.")
                                    .clicked()
                                {
                                    self.world.controls.symbolic_teaching_fidelity =
                                        self.world.controls.symbolic_teaching_fidelity.max(0.90);
                                    self.world.controls.inter_deme_connectivity =
                                        self.world.controls.inter_deme_connectivity.clamp(0.25, 0.55);
                                    let (rt, sym, rec) = self.world.trigger_cultural_ratchet_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Symbolic Culture: +{rt} ratchets, +{sym} pedagogical lineages, +{rec} deme recombinations"
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Agro-Ectosymbiosis & Crops")
                                    .on_hover_text("Trigger clonal cultivar domestication and defensive antimicrobial shields.")
                                    .clicked()
                                {
                                    let (ag, ob, sh) =
                                        self.world.trigger_agriculture_domestication_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Cultivar Care: +{ag} farming demes, +{ob} clonal crops, {sh} defensive shields"
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Written Archive + Peer Audit")
                                    .on_hover_text("Externalize cultural knowledge into written archives with institutional replication audits.")
                                    .clicked()
                                {
                                    let (ar, inst, aud) = self.world.trigger_exosomatic_archive_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Written Archives: +{ar} archives, +{inst} peer institutions, {aud} audits"
                                        ),
                                        false,
                                    ));
                                }
                            }
                            OrganismKind::ArtificialIntelligence => {
                                if ui
                                    .button("Container Checkpoint Scale")
                                    .on_hover_text("Jump to R(1.5-RP) stateless container and checkpoint cloning.")
                                    .clicked()
                                {
                                    self.jump_to_regime(RegimeId::R15Replicative);
                                    self.export_status = Some((
                                        "Scaled stateless container weights across cluster R(1.5-RP)."
                                            .to_string(),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Adversarial Arena (PBT)")
                                    .on_hover_text("Activate population-based training and adversarial policy arenas in R(3.0-DAR).")
                                    .clicked()
                                {
                                    self.jump_to_regime(RegimeId::R30FullyDarwinian);
                                    self.export_status = Some((
                                        "Activated population-based training & adversarial arena R(3.0-DAR)."
                                            .to_string(),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Formal Proof Verifier")
                                    .on_hover_text("Couple generative rollouts to deterministic formal proof-checker kernels.")
                                    .clicked()
                                {
                                    let (ms, fm, prf) = self.world.trigger_machine_synthesis_pulse();
                                    self.export_status = Some((
                                        format!(
                                            "Formal Verifiers: +{ms} compute symbionts, +{fm} verified loops, {prf} proofs checked"
                                        ),
                                        false,
                                    ));
                                }
                                if ui
                                    .button("Anti-Collapse Safeguard")
                                    .on_hover_text("Enforce empirical grounding floor (>= 0.36) and strict formal verification (>= 0.92) to halt model collapse.")
                                    .clicked()
                                {
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
                                    self.export_status = Some((
                                        format!(
                                            "Anti-Collapse Guard active: empiric floor >= 0.36 | TailDiv = {:.2}, Theorems = {:.1}b",
                                            m.mean_tail_diversity, m.mean_verified_theorem_bits
                                        ),
                                        false,
                                    ));
                                }
                            }
                        }
                    }
                });
            }
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
                        )
                        .on_hover_text("Automatically cycle sequentially across all 13 emergent evolutionary regimes.");
                        ui.selectable_value(
                            &mut selected_mode,
                            SimulationMode::OrganicSandbox,
                            "Autonomous Sandbox",
                        )
                        .on_hover_text("Free-running ecological and evolutionary sandbox driven by live parameters.");
                    });
                    self.world.mode = selected_mode;

                    if selected_mode == SimulationMode::GuidedCurriculum {
                        ui.add(
                            egui::Slider::new(&mut self.world.steps_per_stage, 12..=120)
                                .text("Ticks / Regime"),
                        )
                        .on_hover_text("Number of ticks spent in each regime before transitioning during Guided Tour mode.");
                    } else if ui
                        .button("Stochastic Phase Jump")
                        .on_hover_text("Trigger a stochastic jump along outgoing evolutionary regime transition pathways.")
                        .clicked()
                    {
                        self.world.trigger_network_transition();
                    }

                    ui.separator();
                    ui.heading(if self.coupled_biosphere {
                        "Coupled Biosphere Regimes (13)".to_string()
                    } else {
                        format!("{} Regimes (13)", self.organism.short_title())
                    });
                    let active_classified = self.world.classified_regime;
                    let mut clicked_regime = None;
                    for regime in RegimeId::ALL {
                        let is_curr = regime == active_classified;
                        let txt = if self.coupled_biosphere {
                            regime.full_title()
                        } else {
                            regime.domain_title(self.organism)
                        };
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
                    changed |= add_param_slider(
                        ui,
                        &mut w.controls.bath_exergy_flux,
                        -0.6..=1.5,
                        "J_in (Supply / Exergy Flux)",
                        "Thermodynamic free-energy / exergy inflow driving dissipative self-organization.",
                    );
                    changed |= add_param_slider(
                        ui,
                        &mut w.controls.replication_fidelity,
                        0.40..=1.0,
                        "F_0 (Copy / Update Fidelity)",
                        "Baseline per-symbol replication or state-update fidelity before maintenance decay.",
                    );
                    changed |= add_param_slider(
                        ui,
                        &mut w.controls.landauer_bit_cost,
                        0.0..=0.008,
                        "k_B T ln2 (Maintenance Cost)",
                        "Landauer thermodynamic maintenance cost per heritable bit.",
                    );

                    ui.add_space(4.0);
                    match self.organism {
                        OrganismKind::RnaDnaProtein => {
                            ui.heading("Prebiotic & Genetic Code");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.wet_dry_cycle_amplitude,
                                0.0..=1.0,
                                "A_wd (Wet-Dry Condensation)",
                                "Tidal/geothermal wet-dry cycling amplitude driving RNA oligomer condensation.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.hydrothermal_delta_ph,
                                0.0..=4.5,
                                "Delta_pH (Vent FeNiS Gradient)",
                                "Alkaline hydrothermal FeNiS proton-motive pH gradient powering Hadean protocells.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.magnetite_ciss_spin_bias,
                                0.0..=1.0,
                                "b_ciss (Fe3O4 Chiral Spin Bias)",
                                "Magnetite chiral-induced spin selectivity (CISS) bias breaking D-RNA/L-AA symmetry.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.rnr_thiyl_radical_activity,
                                0.0..=1.0,
                                "k_rnr (Class III RNR 2'-OH->2'-H)",
                                "Anaerobic Class III ribonucleotide reductase activity converting RNA to T-DNA.",
                            );

                            ui.add_space(4.0);
                            ui.heading("Symbiogenesis & MRO Continuum");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.endosymbiosis_rate,
                                0.0..=0.80,
                                "gamma_1 (1° Endosymbiosis)",
                                "Primary endosymbiotic fusion rate integrating bacterial symbionts.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.ambient_oxygen_level,
                                0.0..=1.0,
                                "O2_level (Hypoxia < 0.45)",
                                "Ambient dissolved O2; hypoxia (< 0.45) favors hydrogenosome/mitosome reduction.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.cytosolic_suf_hgt_rescue,
                                0.0..=1.0,
                                "suf_rescue (Cytosolic SUF >=0.65)",
                                "Horizontal transfer of bacterial cytosolic SUF Fe-S assembly enabling MRO->0.",
                            );
                        }
                        OrganismKind::HominidPopulation => {
                            ui.heading("Cultural Transmission & Demes");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.symbolic_teaching_fidelity,
                                0.0..=1.0,
                                "f_teach (Symbolic Pedagogy)",
                                "Fidelity of active symbolic teaching and cross-generational skill transmission.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.inter_deme_connectivity,
                                0.0..=1.0,
                                "c_deme (Inter-Band Contact)",
                                "Migration and seasonal exchange connectivity across foraging bands and demes.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.cultural_innovation_rate,
                                0.0..=0.90,
                                "mu_innov (Innovation Drive)",
                                "Rate of novel lithic, ecological, or symbolic cultural innovation.",
                            );

                            ui.add_space(4.0);
                            ui.heading("Subsistence & Written Archives");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.crop_domestication_rate,
                                0.0..=0.90,
                                "mu_crop (Cultivar Domestication)",
                                "Co-evolutionary agro-ectosymbiosis and cultivar selection intensity.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.archival_externalization_rate,
                                0.0..=0.90,
                                "mu_arch (Written Archival Drive)",
                                "Externalization of symbolic memory into durable written records.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.peer_verification_rigor,
                                0.0..=1.0,
                                "r_peer (Peer Replication Rigor)",
                                "Independent verification rigor filtering archival drift and false claims.",
                            );
                        }
                        OrganismKind::ArtificialIntelligence => {
                            ui.heading("Training, Plasticity & Distillation");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.synthetic_generation_rate,
                                0.0..=1.0,
                                "g_synth (Synthetic Rollout Rate)",
                                "Fraction of training distribution drawn from self-generated synthetic rollouts.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.memory_plasticity,
                                0.0..=1.0,
                                "M_plast (In-Context / LoRA)",
                                "Test-time and in-context plasticity during deployment.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.soft_inheritance_coupling,
                                0.0..=1.0,
                                "alpha_distill (Teacher Distillation)",
                                "Strength of teacher-to-student behavioral and representation distillation.",
                            );

                            ui.add_space(4.0);
                            ui.heading("Formal Verification & Anti-Collapse");
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.formal_verifier_strictness,
                                0.0..=1.0,
                                "v_formal (Proof Verifier >=0.72)",
                                "Strictness of formal theorem verifier gating self-improvement.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.empirical_grounding_floor,
                                0.0..=0.60,
                                "gamma_emp (Empirical Anchor >=0.12)",
                                "Minimum real-world empirical grounding floor preventing tail collapse.",
                            );
                            changed |= add_param_slider(
                                ui,
                                &mut w.controls.unverified_publish_pressure,
                                0.0..=1.0,
                                "p_noise (Unverified Publish Rate)",
                                "Rate at which unverified synthetic outputs contaminate the corpus.",
                            );
                        }
                    }

                    let f_eff = w.controls.effective_fidelity();
                    if self.coupled_biosphere {
                        ui.add_space(6.0);
                        ui.separator();
                        ui.heading("Planetary Biogeochemistry & Gaia-2.0");
                        changed |= add_param_slider(
                            ui,
                            &mut w.controls.biogeochemical_feedback_gain,
                            0.0..=1.0,
                            "g_bio (Biosphere Coupling)",
                            "Strength of two-way feedback between organismal metabolism and planetary reservoirs.",
                        );
                        changed |= add_param_slider(
                            ui,
                            &mut w.controls.planetary_stewardship_coupling,
                            0.0..=1.0,
                            "s_gaia (Gaia-2.0 Stewardship)",
                            "Intentional planetary homeostatic regulation stabilizing O2, CO2, and exergy.",
                        );
                        changed |= add_param_slider(
                            ui,
                            &mut w.controls.extractive_burn_intensity,
                            0.0..=1.0,
                            "b_ext (Extractive Burn)",
                            "Unregulated extractive carbon oxidation intensity driving CO2 overshoot.",
                        );
                        if add_param_slider(
                            ui,
                            &mut w.controls.ambient_oxygen_level,
                            0.02..=0.94,
                            "O2 (PAL Index Target)",
                            "Target atmospheric oxygen PAL (Present Atmospheric Level) index [0.02, 0.94] for coupled planetary biosphere runs.",
                        ) {
                            w.biosphere.atmospheric_o2 = w.controls.ambient_oxygen_level;
                            changed = true;
                        }
                    }

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

                    ui.heading(if self.coupled_biosphere {
                        "Emergent Biosphere Regime"
                    } else {
                        "Emergent Domain Regime"
                    });
                    ui.label(
                        egui::RichText::new(if self.coupled_biosphere {
                            reg.full_title()
                        } else {
                            reg.domain_title(self.organism)
                        })
                        .size(15.5)
                        .strong()
                        .color(egui::Color32::from_rgb(255, 214, 102)),
                    );
                    ui.add(
                        egui::Label::new(
                            egui::RichText::new(if self.coupled_biosphere {
                                reg.core_concept()
                            } else {
                                reg.domain_concept(self.organism)
                            })
                            .italics(),
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

                    if let Some((ref status, is_err)) = self.export_status {
                        let status_col = if is_err {
                            egui::Color32::from_rgb(255, 110, 110)
                        } else {
                            egui::Color32::from_rgb(130, 245, 205)
                        };
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(status)
                                    .small()
                                    .color(status_col),
                            )
                            .wrap(),
                        );
                    }

                    if self.coupled_biosphere {
                        ui.separator();
                        ui.heading("Planetary Biosphere Gauges");
                        let b = &w.biosphere;
                        let (st_col, st_bg) = match b.state {
                            PlanetaryState::AnoxicArchean => (
                                egui::Color32::from_rgb(251, 191, 36),
                                egui::Color32::from_rgb(45, 34, 16),
                            ),
                            PlanetaryState::OxygenatedPhanerozoic => (
                                egui::Color32::from_rgb(76, 201, 240),
                                egui::Color32::from_rgb(16, 40, 56),
                            ),
                            PlanetaryState::AnthropoceneDisequilibrium => (
                                egui::Color32::from_rgb(255, 110, 110),
                                egui::Color32::from_rgb(56, 20, 26),
                            ),
                            PlanetaryState::HomeostaticAwakening => (
                                egui::Color32::from_rgb(52, 211, 153),
                                egui::Color32::from_rgb(14, 46, 36),
                            ),
                        };
                        egui::Frame::none()
                            .fill(st_bg)
                            .stroke(egui::Stroke::new(1.0_f32, st_col))
                            .rounding(5.0)
                            .inner_margin(egui::Margin::symmetric(8.0, 4.0))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.label(
                                        egui::RichText::new("PlanetaryState:")
                                            .small()
                                            .color(egui::Color32::from_rgb(203, 213, 225)),
                                    );
                                    ui.label(
                                        egui::RichText::new(b.state.short_code())
                                            .small()
                                            .strong()
                                            .color(st_col),
                                    );
                                });
                            });
                        ui.add_space(4.0);
                        draw_biosphere_gauge(
                            ui,
                            "O2 (PAL Index)",
                            &format!("{:.3} (GOE>=0.25)", b.atmospheric_o2),
                            b.atmospheric_o2 as f32,
                            0.0,
                            1.0,
                            Some(0.25),
                            egui::Color32::from_rgb(76, 201, 240),
                        );
                        draw_biosphere_gauge(
                            ui,
                            "CO2 ppm",
                            &format!("{:.0} ppm (Safe<=420)", b.atmospheric_co2_ppm),
                            b.atmospheric_co2_ppm as f32,
                            220.0,
                            1800.0,
                            Some(420.0),
                            if b.atmospheric_co2_ppm > 420.0 {
                                egui::Color32::from_rgb(255, 125, 105)
                            } else {
                                egui::Color32::from_rgb(244, 162, 97)
                            },
                        );
                        draw_biosphere_gauge(
                            ui,
                            "N/P Pool",
                            &format!("{:.3}", b.nutrient_pool_np),
                            b.nutrient_pool_np as f32,
                            0.0,
                            1.0,
                            Some(0.35),
                            egui::Color32::from_rgb(80, 227, 137),
                        );
                        draw_biosphere_gauge(
                            ui,
                            "Gibbs Diseq",
                            &format!("{:.2} bits", b.chemical_disequilibrium_bits),
                            b.chemical_disequilibrium_bits as f32,
                            0.0,
                            26.0,
                            None,
                            egui::Color32::from_rgb(167, 139, 250),
                        );
                        ui.monospace(format!(
                            "Exergy: {:.2} | Load: {:.2}",
                            b.planetary_exergy_budget, b.boundary_transgression_load
                        ));
                    }

                    ui.separator();
                    ui.heading(if self.coupled_biosphere {
                        "Biosphere Telemetry".to_string()
                    } else {
                        format!("{} Telemetry", self.organism.short_title())
                    });
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
            .min_height(155.0)
            .default_height(168.0)
            .resizable(true)
            .show(ctx, |ui| {
                let w = &self.world;
                ui.horizontal(|ui| {
                    let plot_title = if self.coupled_biosphere {
                        "ITME v23.0 Trajectories — 4.25 Ga PlanetaryBiosphere".to_string()
                    } else {
                        format!("ITME v23.0 Trajectories — {}", self.organism.title())
                    };
                    ui.strong(plot_title);
                    ui.separator();
                    ui.colored_label(egui::Color32::from_rgb(80, 220, 255), "Order (O)");
                    if self.coupled_biosphere {
                        ui.colored_label(egui::Color32::from_rgb(76, 201, 240), "O2 (PAL Index)");
                        ui.colored_label(egui::Color32::from_rgb(80, 227, 137), "N/P Pool");
                        ui.colored_label(egui::Color32::from_rgb(167, 139, 250), "Gibbs/26");
                    } else {
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
                    }
                });

                let avail = ui.available_size();
                let (resp, painter) = ui.allocate_painter(avail, egui::Sense::hover());
                let hover_pos = resp.hover_pos();
                let rect = painter.clip_rect().shrink(4.0);
                paint_framed_box(
                    &painter,
                    rect,
                    6.0,
                    egui::Color32::from_rgb(10, 14, 26),
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(38, 52, 84)),
                );

                if w.history.len() < 2 {
                    painter.text(
                        rect.center() - egui::vec2(0.0, 9.0),
                        egui::Align2::CENTER_CENTER,
                        "Telemetry plots activate after 2+ simulation ticks.",
                        egui::FontId::proportional(13.0),
                        egui::Color32::from_rgb(186, 206, 236),
                    );
                    painter.text(
                        rect.center() + egui::vec2(0.0, 11.0),
                        egui::Align2::CENTER_CENTER,
                        "Press Play (Space) to stream live trajectories or Step (.) to advance one tick.",
                        egui::FontId::monospace(10.5),
                        egui::Color32::from_rgb(125, 211, 252),
                    );
                    return;
                }

                let left_plot = egui::Rect::from_min_max(
                    rect.min + egui::vec2(8.0, 6.0),
                    egui::pos2(rect.center().x - 6.0, rect.max.y - 6.0),
                );
                let right_plot = egui::Rect::from_min_max(
                    egui::pos2(rect.center().x + 6.0, rect.min.y + 6.0),
                    rect.max - egui::vec2(8.0, 6.0),
                );

                let make_inner = |outer: egui::Rect| -> egui::Rect {
                    egui::Rect::from_min_max(
                        egui::pos2(outer.left() + 42.0, outer.top() + 18.0),
                        egui::pos2(outer.right() - 10.0, outer.bottom() - 16.0),
                    )
                };
                let inner_left = make_inner(left_plot);
                let inner_right = make_inner(right_plot);

                let paint_grid_axes = |inner: egui::Rect, y_max: f32, y_decimals: usize| {
                    if inner.width() < 24.0 || inner.height() < 16.0 {
                        return;
                    }
                    let grid_col = egui::Color32::from_rgba_unmultiplied(44, 62, 102, 130);
                    let lbl_col = egui::Color32::from_rgb(138, 158, 196);
                    let n_hist = w.history.len();
                    for k in 0..=4 {
                        let frac = (k as f32) / 4.0;
                        let y = inner.bottom() - frac * inner.height();
                        painter.line_segment(
                            [egui::pos2(inner.left(), y), egui::pos2(inner.right(), y)],
                            egui::Stroke::new(0.8_f32, grid_col),
                        );
                        let y_val = y_max * frac;
                        let y_lbl = match y_decimals {
                            0 => format!("{y_val:.0}"),
                            1 => format!("{y_val:.1}"),
                            _ => format!("{y_val:.2}"),
                        };
                        painter.text(
                            egui::pos2(inner.left() - 4.0, y),
                            egui::Align2::RIGHT_CENTER,
                            y_lbl,
                            egui::FontId::monospace(8.5),
                            lbl_col,
                        );

                        let x = inner.left() + frac * inner.width();
                        painter.line_segment(
                            [egui::pos2(x, inner.top()), egui::pos2(x, inner.bottom())],
                            egui::Stroke::new(0.8_f32, grid_col),
                        );
                        let h_idx = ((frac * ((n_hist - 1) as f32)).round() as usize).min(n_hist - 1);
                        let ga = self.step_to_ga(w.history[h_idx].step);
                        let align = if k == 0 {
                            egui::Align2::LEFT_TOP
                        } else if k == 4 {
                            egui::Align2::RIGHT_TOP
                        } else {
                            egui::Align2::CENTER_TOP
                        };
                        painter.text(
                            egui::pos2(x, inner.bottom() + 2.0),
                            align,
                            format!("{ga:.2}Ga"),
                            egui::FontId::monospace(8.5),
                            lbl_col,
                        );
                    }
                    painter.rect_stroke(
                        inner,
                        0.0,
                        egui::Stroke::new(0.9_f32, egui::Color32::from_rgb(44, 62, 100)),
                    );
                };

                let max_val = w
                    .history
                    .iter()
                    .map(|m| m.order_o.max(1.4))
                    .fold(1.4_f64, f64::max) as f32;
                paint_grid_axes(inner_left, max_val, 2);
                draw_series(
                    &painter,
                    inner_left,
                    &w.history,
                    |m| m.order_o as f32,
                    0.0,
                    max_val,
                    egui::Color32::from_rgb(80, 220, 255),
                );
                painter.text(
                    egui::pos2(inner_left.left(), left_plot.top() + 2.0),
                    egui::Align2::LEFT_TOP,
                    format!("Order O [0 .. {:.2}]", max_val),
                    egui::FontId::monospace(10.0),
                    egui::Color32::from_rgb(175, 190, 220),
                );

                let right_y_max: f32;
                if self.coupled_biosphere {
                    right_y_max = 1.0;
                    paint_grid_axes(inner_right, right_y_max, 2);
                    draw_series(
                        &painter,
                        inner_right,
                        &w.history,
                        |m| m.planetary_o2_level as f32,
                        0.0,
                        1.0,
                        egui::Color32::from_rgb(76, 201, 240),
                    );
                    draw_series(
                        &painter,
                        inner_right,
                        &w.history,
                        |m| m.planetary_nutrient_pool as f32,
                        0.0,
                        1.0,
                        egui::Color32::from_rgb(80, 227, 137),
                    );
                    draw_series(
                        &painter,
                        inner_right,
                        &w.history,
                        |m| (m.chemical_disequilibrium_bits as f32) / 26.0,
                        0.0,
                        1.0,
                        egui::Color32::from_rgb(167, 139, 250),
                    );
                    painter.text(
                        egui::pos2(inner_right.left(), right_plot.top() + 2.0),
                        egui::Align2::LEFT_TOP,
                        format!(
                            "O2 (PAL)={:.2} | N/P={:.2} | Gibbs={:.1}b | CO2={:.0}ppm",
                            w.biosphere.atmospheric_o2,
                            w.biosphere.nutrient_pool_np,
                            w.biosphere.chemical_disequilibrium_bits,
                            w.biosphere.atmospheric_co2_ppm
                        ),
                        egui::FontId::monospace(10.0),
                        egui::Color32::from_rgb(175, 190, 220),
                    );
                } else {
                    match self.organism {
                    OrganismKind::RnaDnaProtein => {
                        let max_right = w
                            .history
                            .iter()
                            .map(|m| m.hydrogenosome_count.max(m.mitosome_count) as f32)
                            .fold(10.0_f32, f32::max)
                            * 1.15;
                        right_y_max = max_right;
                        paint_grid_axes(inner_right, max_right, 0);
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| m.hydrogenosome_count as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(6, 214, 160),
                        );
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| m.mitosome_count as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(255, 209, 102),
                        );
                        painter.text(
                            egui::pos2(inner_right.left(), right_plot.top() + 2.0),
                            egui::Align2::LEFT_TOP,
                            format!("Hydrogenosomes & Mitosomes [0 .. {:.0}]", max_right),
                            egui::FontId::monospace(10.0),
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
                        right_y_max = max_right;
                        paint_grid_axes(inner_right, max_right, 1);
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| m.mean_cultural_repertoire_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(6, 214, 160),
                        );
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| m.mean_exosomatic_archive_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(255, 209, 102),
                        );
                        painter.text(
                            egui::pos2(inner_right.left(), right_plot.top() + 2.0),
                            egui::Align2::LEFT_TOP,
                            format!("Cultural & Written Archive Bits [0 .. {:.1}]", max_right),
                            egui::FontId::monospace(10.0),
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
                        right_y_max = max_right;
                        paint_grid_axes(inner_right, max_right, 1);
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| m.mean_verified_theorem_bits as f32,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(6, 214, 160),
                        );
                        draw_series(
                            &painter,
                            inner_right,
                            &w.history,
                            |m| (m.mean_tail_diversity as f32) * 10.0,
                            0.0,
                            max_right,
                            egui::Color32::from_rgb(199, 125, 255),
                        );
                        painter.text(
                            egui::pos2(inner_right.left(), right_plot.top() + 2.0),
                            egui::Align2::LEFT_TOP,
                            format!("Verified Theorem Bits & TailDiv [0 .. {:.1}]", max_right),
                            egui::FontId::monospace(10.0),
                            egui::Color32::from_rgb(175, 190, 220),
                        );
                    }
                    }
                }

                // Linked hover crosshair + Ga/series value readout across both telemetry plots
                let hovered_plot = hover_pos.and_then(|hp| {
                    if inner_left.width() > 4.0 && inner_left.contains(hp) {
                        Some((true, inner_left, hp))
                    } else if inner_right.width() > 4.0 && inner_right.contains(hp) {
                        Some((false, inner_right, hp))
                    } else {
                        None
                    }
                });
                if let Some((is_left, active_inner, hp)) = hovered_plot {
                    let n = w.history.len();
                    let t_x = ((hp.x - active_inner.left()) / active_inner.width()).clamp(0.0, 1.0);
                    let idx = ((t_x * ((n - 1) as f32)).round() as usize).min(n - 1);
                    let x_frac = (idx as f32) / ((n - 1) as f32);
                    let m_at = &w.history[idx];
                    let ga_at = self.step_to_ga(m_at.step);

                    let cross_col = egui::Color32::from_rgba_unmultiplied(226, 232, 240, 145);
                    for plot_box in [inner_left, inner_right] {
                        if plot_box.width() > 4.0 && plot_box.height() > 4.0 {
                            let sx = plot_box.left() + x_frac * plot_box.width();
                            painter.line_segment(
                                [egui::pos2(sx, plot_box.top()), egui::pos2(sx, plot_box.bottom())],
                                egui::Stroke::new(1.0_f32, cross_col),
                            );
                        }
                    }
                    painter.line_segment(
                        [
                            egui::pos2(active_inner.left(), hp.y),
                            egui::pos2(active_inner.right(), hp.y),
                        ],
                        egui::Stroke::new(
                            0.9_f32,
                            egui::Color32::from_rgba_unmultiplied(148, 163, 184, 110),
                        ),
                    );

                    let plot_pt = |box_r: egui::Rect, val: f32, ymax: f32| -> egui::Pos2 {
                        let sx = box_r.left() + x_frac * box_r.width();
                        let ny = (val / ymax.max(1e-5)).clamp(0.0, 1.0);
                        egui::pos2(sx, box_r.bottom() - ny * box_r.height())
                    };
                    let o_pt = plot_pt(inner_left, m_at.order_o as f32, max_val);
                    painter.circle_filled(o_pt, 3.3, egui::Color32::from_rgb(80, 220, 255));
                    painter.circle_stroke(o_pt, 3.8, egui::Stroke::new(1.0_f32, egui::Color32::WHITE));

                    let y_frac = ((active_inner.bottom() - hp.y) / active_inner.height()).clamp(0.0, 1.0);
                    let cursor_y = y_frac * if is_left { max_val } else { right_y_max };

                    let l1 = format!("Tick {} | {:.3} Ga", m_at.step, ga_at);
                    let (l2, l3, accent) = if is_left {
                        (
                            format!(
                                "Order O: {:.3} (L_I={:.2}, S_T-={:.2})",
                                m_at.order_o, m_at.leverage_li, m_at.stability_st
                            ),
                            format!("O_pop: {:.2} | Cursor Y: {:.3}", m_at.ecosystem_order_pop, cursor_y),
                            egui::Color32::from_rgb(80, 220, 255),
                        )
                    } else if self.coupled_biosphere {
                        let p_o2 = plot_pt(inner_right, m_at.planetary_o2_level as f32, 1.0);
                        let p_np = plot_pt(inner_right, m_at.planetary_nutrient_pool as f32, 1.0);
                        let p_gb = plot_pt(inner_right, (m_at.chemical_disequilibrium_bits as f32) / 26.0, 1.0);
                        painter.circle_filled(p_o2, 3.1, egui::Color32::from_rgb(76, 201, 240));
                        painter.circle_filled(p_np, 3.1, egui::Color32::from_rgb(80, 227, 137));
                        painter.circle_filled(p_gb, 3.1, egui::Color32::from_rgb(167, 139, 250));
                        (
                            format!(
                                "O2 (PAL): {:.3} | N/P: {:.3}",
                                m_at.planetary_o2_level, m_at.planetary_nutrient_pool
                            ),
                            format!(
                                "Gibbs: {:.2}b | CO2: {:.0}ppm | Y: {:.3}",
                                m_at.chemical_disequilibrium_bits, m_at.planetary_co2_ppm, cursor_y
                            ),
                            egui::Color32::from_rgb(76, 201, 240),
                        )
                    } else {
                        match self.organism {
                            OrganismKind::RnaDnaProtein => {
                                let p1 = plot_pt(inner_right, m_at.hydrogenosome_count as f32, right_y_max);
                                let p2 = plot_pt(inner_right, m_at.mitosome_count as f32, right_y_max);
                                painter.circle_filled(p1, 3.1, egui::Color32::from_rgb(6, 214, 160));
                                painter.circle_filled(p2, 3.1, egui::Color32::from_rgb(255, 209, 102));
                                (
                                    format!(
                                        "Hydro: {} | MitoS: {} | MRO->0: {}",
                                        m_at.hydrogenosome_count, m_at.mitosome_count, m_at.amiteukaryote_count
                                    ),
                                    format!("Order O: {:.3} | Cursor Y: {:.1}", m_at.order_o, cursor_y),
                                    egui::Color32::from_rgb(6, 214, 160),
                                )
                            }
                            OrganismKind::HominidPopulation => {
                                let p1 = plot_pt(inner_right, m_at.mean_cultural_repertoire_bits as f32, right_y_max);
                                let p2 = plot_pt(inner_right, m_at.mean_exosomatic_archive_bits as f32, right_y_max);
                                painter.circle_filled(p1, 3.1, egui::Color32::from_rgb(6, 214, 160));
                                painter.circle_filled(p2, 3.1, egui::Color32::from_rgb(255, 209, 102));
                                (
                                    format!(
                                        "Culture: {:.2}b | Archive: {:.2}b",
                                        m_at.mean_cultural_repertoire_bits, m_at.mean_exosomatic_archive_bits
                                    ),
                                    format!("Order O: {:.3} | Cursor Y: {:.2}b", m_at.order_o, cursor_y),
                                    egui::Color32::from_rgb(6, 214, 160),
                                )
                            }
                            OrganismKind::ArtificialIntelligence => {
                                let p1 = plot_pt(inner_right, m_at.mean_verified_theorem_bits as f32, right_y_max);
                                let p2 = plot_pt(inner_right, (m_at.mean_tail_diversity as f32) * 10.0, right_y_max);
                                painter.circle_filled(p1, 3.1, egui::Color32::from_rgb(6, 214, 160));
                                painter.circle_filled(p2, 3.1, egui::Color32::from_rgb(199, 125, 255));
                                (
                                    format!(
                                        "Theorems: {:.2}b | TailDiv: {:.3}",
                                        m_at.mean_verified_theorem_bits, m_at.mean_tail_diversity
                                    ),
                                    format!("Order O: {:.3} | Cursor Y: {:.2}", m_at.order_o, cursor_y),
                                    egui::Color32::from_rgb(199, 125, 255),
                                )
                            }
                        }
                    };
                    paint_hover_tooltip(&painter, rect, hp, &[&l1, &l2, &l3], accent);
                }
            });
    }

    fn draw_central_viewport(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let avail = ui.available_size();
            let (response, painter) = ui.allocate_painter(avail, egui::Sense::click_and_drag());
            let hover_pos = response.hover_pos();

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
                    let left_sub = left_r.shrink(4.0);
                    paint_spatial_grid(
                        &painter,
                        left_sub,
                        w,
                        self.show_labels,
                        hover_pos.filter(|p| left_sub.contains(*p)),
                    );
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
                    let left_sub = left_r.shrink(4.0);
                    paint_spatial_grid(
                        &painter,
                        left_sub,
                        w,
                        self.show_labels,
                        hover_pos.filter(|p| left_sub.contains(*p)),
                    );
                    paint_scene_3d(
                        &painter,
                        right_r,
                        &scene,
                        &self.camera,
                        self.show_labels,
                        hover_pos.filter(|p| right_r.contains(*p)),
                    );
                }
                ViewportLayout::PhylogenyFull => {
                    paint_phylogram_tree(&painter, full_rect.shrink(8.0), &self.world, self.show_labels);
                }
                ViewportLayout::Perspective3D => {
                    let scene = self.world.current_scene();
                    let sub_r = full_rect.shrink(6.0);
                    paint_scene_3d(
                        &painter,
                        sub_r,
                        &scene,
                        &self.camera,
                        self.show_labels,
                        hover_pos.filter(|p| sub_r.contains(*p)),
                    );
                }
                ViewportLayout::Blueprint2D => {
                    let scene = self.world.current_scene();
                    let sub_r = full_rect.shrink(6.0);
                    paint_scene_2d(
                        &painter,
                        sub_r,
                        &scene,
                        self.show_labels,
                        hover_pos.filter(|p| sub_r.contains(*p)),
                    );
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
                "Click '4.25 Ga Grand Trajectory' in the top bar to execute the 8-Epoch Evolutionary Trajectory.",
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
                "Final: Ch:{}/5 | MI:{:.1}b | L_I:{:.2} | O2 (PAL):{:.3} | CO2:{:.0} ppm",
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
                "Panel A | Atmosphere: O2 PAL Index (Cyan) & log10(CO2) (Amber)",
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
            "O2(PAL)",
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

fn draw_biosphere_gauge(
    ui: &mut egui::Ui,
    label: &str,
    value_text: &str,
    val: f32,
    min_val: f32,
    max_val: f32,
    threshold: Option<f32>,
    bar_color: egui::Color32,
) {
    let span = (max_val - min_val).max(1e-5);
    let frac = ((val - min_val) / span).clamp(0.0, 1.0);
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(label)
                .small()
                .strong()
                .color(egui::Color32::from_rgb(203, 213, 225)),
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(value_text)
                    .small()
                    .monospace()
                    .color(bar_color),
            );
        });
    });
    let bar_w = ui.available_width().max(60.0);
    let (bar_rect, _) = ui.allocate_exact_size(egui::vec2(bar_w, 8.0), egui::Sense::hover());
    let p = ui.painter();
    p.rect_filled(bar_rect, 3.0, egui::Color32::from_rgb(16, 24, 44));
    p.rect_stroke(
        bar_rect,
        3.0,
        egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(44, 62, 100)),
    );
    let fill_w = (bar_rect.width() * frac).max(if frac > 0.0 { 2.0 } else { 0.0 });
    let fill_rect =
        egui::Rect::from_min_size(bar_rect.min, egui::vec2(fill_w, bar_rect.height()));
    p.rect_filled(fill_rect, 3.0, bar_color);
    if let Some(th) = threshold {
        let th_frac = ((th - min_val) / span).clamp(0.0, 1.0);
        let tx = bar_rect.left() + bar_rect.width() * th_frac;
        p.line_segment(
            [
                egui::pos2(tx, bar_rect.top() - 1.5),
                egui::pos2(tx, bar_rect.bottom() + 1.5),
            ],
            egui::Stroke::new(1.4_f32, egui::Color32::from_rgb(248, 250, 252)),
        );
    }
    ui.add_space(3.0);
}

pub fn rect_overlaps_any(candidate: egui::Rect, occupied: &[egui::Rect]) -> bool {
    occupied.iter().any(|r| r.intersects(candidate))
}

fn clamp_pill_rect(rect: egui::Rect, bounds: egui::Rect) -> egui::Rect {
    let min_x = bounds.left();
    let max_x = (bounds.right() - rect.width()).max(min_x);
    let min_y = bounds.top();
    let max_y = (bounds.bottom() - rect.height()).max(min_y);
    let dx = rect.left().clamp(min_x, max_x) - rect.left();
    let dy = rect.top().clamp(min_y, max_y) - rect.top();
    rect.translate(egui::vec2(dx, dy))
}

#[allow(clippy::too_many_arguments)]
fn paint_pill_label(
    painter: &egui::Painter,
    bounds: egui::Rect,
    anchor: egui::Pos2,
    align: egui::Align2,
    text: &str,
    font: egui::FontId,
    text_col: egui::Color32,
    border_col: egui::Color32,
    occupied: &mut Vec<egui::Rect>,
    force: bool,
) -> bool {
    if text.is_empty() || bounds.width() < 24.0 || bounds.height() < 16.0 {
        return false;
    }
    let galley = painter.layout_no_wrap(text.to_string(), font, text_col);
    let pad = egui::vec2(5.0, 2.0);
    let pill_size = galley.size() + pad * 2.0;
    let safe_bounds = bounds.shrink(3.0);

    let make_rect = |pos: egui::Pos2| -> egui::Rect {
        let min_x = match align.x() {
            egui::Align::Min => pos.x,
            egui::Align::Center => pos.x - pill_size.x * 0.5,
            egui::Align::Max => pos.x - pill_size.x,
        };
        let min_y = match align.y() {
            egui::Align::Min => pos.y,
            egui::Align::Center => pos.y - pill_size.y * 0.5,
            egui::Align::Max => pos.y - pill_size.y,
        };
        clamp_pill_rect(
            egui::Rect::from_min_size(egui::pos2(min_x, min_y), pill_size),
            safe_bounds,
        )
    };

    let base_rect = make_rect(anchor);
    let chosen_rect = if force || !rect_overlaps_any(base_rect.expand(1.5), occupied) {
        Some(base_rect)
    } else {
        let step = pill_size.y + 4.0;
        let up_rect = make_rect(anchor + egui::vec2(0.0, -step));
        if !rect_overlaps_any(up_rect.expand(1.5), occupied) {
            Some(up_rect)
        } else {
            let down_rect = make_rect(anchor + egui::vec2(0.0, step));
            if !rect_overlaps_any(down_rect.expand(1.5), occupied) {
                Some(down_rect)
            } else {
                None
            }
        }
    };

    let Some(pill_rect) = chosen_rect else {
        return false;
    };

    occupied.push(pill_rect);
    painter.rect_filled(
        pill_rect,
        4.0,
        egui::Color32::from_rgba_unmultiplied(10, 15, 28, 224),
    );
    painter.rect_stroke(pill_rect, 4.0, egui::Stroke::new(0.9_f32, border_col));
    painter.galley(pill_rect.min + pad, galley, text_col);
    true
}

fn paint_hover_tooltip(
    painter: &egui::Painter,
    bounds: egui::Rect,
    pos: egui::Pos2,
    lines: &[&str],
    accent: egui::Color32,
) {
    if lines.is_empty() {
        return;
    }
    let font = egui::FontId::monospace(10.0);
    let mut galleys = Vec::with_capacity(lines.len());
    let mut max_w = 0.0_f32;
    let mut total_h = 0.0_f32;
    for (idx, line) in lines.iter().enumerate() {
        let col = if idx == 0 {
            accent
        } else {
            egui::Color32::from_rgb(226, 232, 240)
        };
        let g = painter.layout_no_wrap((*line).to_string(), font.clone(), col);
        max_w = max_w.max(g.size().x);
        total_h += g.size().y + 2.0;
        galleys.push((g, col));
    }
    let pad = egui::vec2(8.0, 6.0);
    let box_size = egui::vec2(max_w + pad.x * 2.0, total_h + pad.y * 2.0 - 2.0);
    let raw_rect = egui::Rect::from_min_size(pos + egui::vec2(14.0, 14.0), box_size);
    let tip_rect = clamp_pill_rect(raw_rect, bounds.shrink(4.0));
    painter.rect_filled(
        tip_rect,
        6.0,
        egui::Color32::from_rgba_unmultiplied(8, 13, 24, 242),
    );
    painter.rect_stroke(tip_rect, 6.0, egui::Stroke::new(1.2_f32, accent));
    let mut y = tip_rect.top() + pad.y;
    for (g, col) in galleys {
        let h = g.size().y;
        painter.galley(egui::pos2(tip_rect.left() + pad.x, y), g, col);
        y += h + 2.0;
    }
}

fn paint_spatial_grid(
    painter: &egui::Painter,
    rect: egui::Rect,
    world: &SimulationWorld,
    show_labels: bool,
    hover_pos: Option<egui::Pos2>,
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
    struct SpatialLabelCandidate {
        anchor: egui::Pos2,
        text: String,
        border_col: egui::Color32,
        is_champ: bool,
    }
    let mut label_candidates: Vec<SpatialLabelCandidate> = Vec::new();
    let mut hovered_unit: Option<(usize, egui::Pos2, f32, f32)> = None;
    let mut seen_label_keys: Vec<(usize, u8)> = Vec::new();

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

        if let Some(hp) = hover_pos {
            let d = center.distance(hp);
            let hit_r = (rad + 7.0).max(11.0);
            if d <= hit_r && hovered_unit.map_or(true, |(_, _, _, best_d)| d < best_d) {
                hovered_unit = Some((idx, center, rad, d));
            }
        }

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

        if show_labels {
            let key = (unit.species_id, unit.guild as u8);
            if is_champ || (!seen_label_keys.contains(&key) && seen_label_keys.len() < 6) {
                if !is_champ {
                    seen_label_keys.push(key);
                }
                let lbl = match world.organism {
                    OrganismKind::RnaDnaProtein => {
                        let base = format!("{} {}", unit.guild.short_code(), unit.mro_class.short_code());
                        if is_champ {
                            format!("* {base}")
                        } else {
                            base
                        }
                    }
                    OrganismKind::HominidPopulation => {
                        if is_champ {
                            format!("* Band {}", unit.species_id + 1)
                        } else {
                            format!("Band {}", unit.species_id + 1)
                        }
                    }
                    OrganismKind::ArtificialIntelligence => {
                        if is_champ {
                            format!("* Agent {}", unit.species_id + 1)
                        } else {
                            format!("Agent {}", unit.species_id + 1)
                        }
                    }
                };
                label_candidates.push(SpatialLabelCandidate {
                    anchor: center + egui::vec2(0.0, -rad - 4.0),
                    text: lbl,
                    border_col: if is_champ {
                        egui::Color32::from_rgb(251, 191, 36)
                    } else {
                        fill
                    },
                    is_champ,
                });
            }
        }
    }

    // Post-glyph label pass with dark pill backdrops and AABB overlap culling
    if show_labels {
        label_candidates.sort_by_key(|c| !c.is_champ);
        let mut occupied: Vec<egui::Rect> = Vec::new();
        for cand in label_candidates {
            paint_pill_label(
                painter,
                rect,
                cand.anchor,
                egui::Align2::CENTER_BOTTOM,
                &cand.text,
                egui::FontId::monospace(9.5),
                if cand.is_champ {
                    egui::Color32::from_rgb(250, 252, 255)
                } else {
                    egui::Color32::from_rgb(215, 226, 245)
                },
                cand.border_col,
                &mut occupied,
                cand.is_champ,
            );
        }
    }

    // Hover highlight and rich unit / cell tooltip
    if let Some((u_idx, center, rad, _)) = hovered_unit {
        let u = &world.population[u_idx];
        painter.circle_stroke(
            center,
            rad + 4.0,
            egui::Stroke::new(2.0_f32, egui::Color32::from_rgb(251, 191, 36)),
        );
        let gx = u.grid_x;
        let gy = u.grid_y;
        let lev = u.complexity_c / u.structure_s.max(2.0);
        let l1 = format!(
            "Unit #{} | Sp {} | {} | MRO:{}",
            u.id,
            u.species_id + 1,
            u.guild.short_code(),
            u.mro_class.short_code()
        );
        let l2 = format!(
            "S_eff:{:.1}b  C:{:.1}  L_I:{:.2}  Fit:{:.2}",
            u.structure_s, u.complexity_c, lev, u.fitness
        );
        let l3 = format!(
            "Cell ({},{}) E:{:.2} B:{:.2} N:{:.2} | Endo:{}°",
            gx,
            gy,
            world.grid.local_exergy(gx, gy),
            world.grid.local_byproduct(gx, gy),
            world.grid.local_niche(gx, gy),
            u.endo_tier
        );
        if let Some(hp) = hover_pos {
            paint_hover_tooltip(
                painter,
                rect,
                hp,
                &[&l1, &l2, &l3],
                egui::Color32::from_rgb(251, 191, 36),
            );
        }
    } else if let Some(hp) = hover_pos {
        let grid_box = egui::Rect::from_min_size(origin, egui::vec2(side, side));
        if grid_box.contains(hp) && cw > 0.0 && ch > 0.0 {
            let gx = (((hp.x - origin.x) / cw).floor() as usize).min(world.grid.width - 1);
            let gy = (((hp.y - origin.y) / ch).floor() as usize).min(world.grid.height - 1);
            let biome = world.grid.biome_at(gx, gy);
            let l1 = format!("Grid Cell ({gx}, {gy}) [{biome:?}]");
            let l2 = format!(
                "Exergy:{:.2}  Byproduct:{:.2}  Niche:{:.2}",
                world.grid.local_exergy(gx, gy),
                world.grid.local_byproduct(gx, gy),
                world.grid.local_niche(gx, gy)
            );
            paint_hover_tooltip(
                painter,
                rect,
                hp,
                &[&l1, &l2],
                egui::Color32::from_rgb(76, 201, 240),
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
    hover_pos: Option<egui::Pos2>,
) {
    let painter = &painter.with_clip_rect(rect);
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

    let fog_for_depth = |d: f64| -> f32 {
        ((d - (camera.distance - 2.1)) / 4.2).clamp(0.0, 1.0) as f32
    };
    let apply_fog = |c: [u8; 3], fog_t: f32, alpha: u8| -> egui::Color32 {
        let mix = (fog_t * 0.45).clamp(0.0, 0.65);
        let r = (f32::from(c[0]) * (1.0 - mix) + 9.0 * mix).clamp(0.0, 255.0) as u8;
        let g = (f32::from(c[1]) * (1.0 - mix) + 13.0 * mix).clamp(0.0, 255.0) as u8;
        let b = (f32::from(c[2]) * (1.0 - mix) + 24.0 * mix).clamp(0.0, 255.0) as u8;
        egui::Color32::from_rgba_unmultiplied(r, g, b, alpha)
    };
    let shade_rgb = |c: [u8; 3], fog_t: f32, light_scale: f32, spec_mix: f32| -> egui::Color32 {
        let mix = (fog_t * 0.45).clamp(0.0, 0.65);
        let br = f32::from(c[0]) * (1.0 - mix) + 9.0 * mix;
        let bg = f32::from(c[1]) * (1.0 - mix) + 13.0 * mix;
        let bb = f32::from(c[2]) * (1.0 - mix) + 24.0 * mix;
        let r = (br * light_scale * (1.0 - spec_mix) + 248.0 * spec_mix).clamp(0.0, 255.0) as u8;
        let g = (bg * light_scale * (1.0 - spec_mix) + 252.0 * spec_mix).clamp(0.0, 255.0) as u8;
        let b = (bb * light_scale * (1.0 - spec_mix) + 255.0 * spec_mix).clamp(0.0, 255.0) as u8;
        egui::Color32::from_rgb(r, g, b)
    };

    enum DrawItem3D {
        Segment {
            p1: egui::Pos2,
            p2: egui::Pos2,
            depth: f64,
            width: f32,
            color: egui::Color32,
        },
        Node {
            idx: usize,
            depth: f64,
        },
    }
    impl DrawItem3D {
        fn depth(&self) -> f64 {
            match self {
                Self::Segment { depth, .. } => *depth,
                Self::Node { depth, .. } => *depth,
            }
        }
        fn kind_rank(&self) -> u8 {
            match self {
                Self::Segment { .. } => 0,
                Self::Node { .. } => 1,
            }
        }
    }

    let mut draw_items: Vec<DrawItem3D> =
        Vec::with_capacity(scene.nodes_3d.len() + scene.edges_3d.len() * 4);

    let view_expand = rect.expand(48.0);
    for e in &scene.edges_3d {
        if e.a >= proj.len() || e.b >= proj.len() {
            continue;
        }
        let (p1, d1, _) = proj[e.a];
        let (p2, d2, _) = proj[e.b];
        if !view_expand.intersects(egui::Rect::from_two_pos(p1, p2)) {
            continue;
        }
        let len_px = p1.distance(p2);
        if len_px < 0.5 {
            continue;
        }
        let alpha = if e.dashed { 185 } else { 225 };
        if e.dashed {
            let dash_px = 6.0_f32;
            let gap_px = 4.5_f32;
            let period = dash_px + gap_px;
            let mut s = 0.0_f32;
            while s < len_px {
                let t0 = (s / len_px).clamp(0.0, 1.0);
                let t1 = ((s + dash_px) / len_px).clamp(0.0, 1.0);
                if t1 > t0 + 1e-4 {
                    let q1 = p1.lerp(p2, t0);
                    let q2 = p1.lerp(p2, t1);
                    let d_mid = d1 + (d2 - d1) * (0.5 * f64::from(t0 + t1));
                    let fog_t = fog_for_depth(d_mid);
                    draw_items.push(DrawItem3D::Segment {
                        p1: q1,
                        p2: q2,
                        depth: d_mid,
                        width: e.width as f32,
                        color: apply_fog(e.color, fog_t, alpha),
                    });
                }
                s += period;
            }
        } else {
            let segs = if (d1 - d2).abs() > 0.35 || len_px > 60.0 {
                4
            } else {
                1
            };
            for k in 0..segs {
                let t0 = (k as f32) / (segs as f32);
                let t1 = ((k + 1) as f32) / (segs as f32);
                let q1 = p1.lerp(p2, t0);
                let q2 = p1.lerp(p2, t1);
                let d_mid = d1 + (d2 - d1) * (0.5 * f64::from(t0 + t1));
                let fog_t = fog_for_depth(d_mid);
                draw_items.push(DrawItem3D::Segment {
                    p1: q1,
                    p2: q2,
                    depth: d_mid,
                    width: e.width as f32,
                    color: apply_fog(e.color, fog_t, alpha),
                });
            }
        }
    }

    for (idx, &(_, depth, _)) in proj.iter().enumerate() {
        draw_items.push(DrawItem3D::Node { idx, depth });
    }

    // Back-to-front depth-interleaved painter's pass across edges and nodes
    draw_items.sort_by(|a, b| {
        b.depth()
            .total_cmp(&a.depth())
            .then_with(|| a.kind_rank().cmp(&b.kind_rank()))
    });

    struct Label3D {
        center: egui::Pos2,
        rad: f32,
        depth: f64,
        text: String,
        col: egui::Color32,
    }
    let mut labels_3d: Vec<Label3D> = Vec::new();
    let mut hovered_node: Option<(usize, egui::Pos2, f32, f32, f64)> = None;

    for item in draw_items {
        match item {
            DrawItem3D::Segment {
                p1,
                p2,
                depth: _,
                width,
                color,
            } => {
                painter.line_segment([p1, p2], egui::Stroke::new(width, color));
            }
            DrawItem3D::Node { idx, depth } => {
                let node = &scene.nodes_3d[idx];
                let (center, _, persp) = proj[idx];
                if !rect.expand(24.0).contains(center) {
                    continue;
                }
                let rad = ((node.radius * persp) as f32).clamp(2.5, 22.0);
                let fog_t = fog_for_depth(depth);
                let glow_f = (node.glow as f32).clamp(0.0, 1.0);

                if glow_f > 0.05 {
                    let outer_a =
                        ((20.0 + 36.0 * glow_f) * (1.0 - 0.35 * fog_t)).clamp(4.0, 72.0) as u8;
                    let inner_a =
                        ((42.0 + 52.0 * glow_f) * (1.0 - 0.30 * fog_t)).clamp(8.0, 110.0) as u8;
                    painter.circle_filled(
                        center,
                        rad * (1.42 + 0.58 * glow_f),
                        apply_fog(node.color, fog_t, outer_a),
                    );
                    painter.circle_filled(
                        center,
                        rad * (1.16 + 0.24 * glow_f),
                        apply_fog(node.color, fog_t, inner_a),
                    );
                }

                // 3D sphere shading: shadow rim + lit hemisphere + top-left specular highlight
                let rim_col = shade_rgb(node.color, fog_t, 0.58, 0.0);
                let body_col = shade_rgb(node.color, fog_t, 1.0, 0.0);
                let spec_mix = 0.62 * (1.0 - 0.55 * fog_t);
                let spec_col = shade_rgb(node.color, fog_t, 1.08, spec_mix);

                painter.circle_filled(center, rad, rim_col);
                painter.circle_filled(
                    center + egui::vec2(-rad * 0.13, -rad * 0.13),
                    rad * 0.85,
                    body_col,
                );
                if rad >= 3.2 {
                    painter.circle_filled(
                        center + egui::vec2(-rad * 0.30, -rad * 0.30),
                        (rad * 0.32).max(1.0),
                        spec_col,
                    );
                }
                painter.circle_stroke(
                    center,
                    rad,
                    egui::Stroke::new(0.7_f32, egui::Color32::from_rgb(8, 12, 22)),
                );

                if let Some(hp) = hover_pos {
                    let d = center.distance(hp);
                    let hit_r = (rad + 7.0).max(11.0);
                    if d <= hit_r
                        && hovered_node.map_or(true, |(_, _, _, best_d, best_z)| {
                            depth < best_z - 0.15 || ((depth - best_z).abs() <= 0.15 && d < best_d)
                        })
                    {
                        hovered_node = Some((idx, center, rad, d, depth));
                    }
                }

                if show_labels && !node.label.is_empty() && rad >= 4.2 {
                    labels_3d.push(Label3D {
                        center,
                        rad,
                        depth,
                        text: node.label.clone(),
                        col: body_col,
                    });
                }
            }
        }
    }

    // Post-glyph label pass: closest / largest foreground nodes first, with pill backdrops & AABB overlap culling
    if show_labels {
        labels_3d.sort_by(|a, b| a.depth.total_cmp(&b.depth).then_with(|| b.rad.total_cmp(&a.rad)));
        let mut occupied: Vec<egui::Rect> = Vec::new();
        for lbl in labels_3d {
            paint_pill_label(
                painter,
                rect,
                lbl.center + egui::vec2(lbl.rad + 4.0, 0.0),
                egui::Align2::LEFT_CENTER,
                &lbl.text,
                egui::FontId::monospace(9.5),
                egui::Color32::from_rgb(232, 240, 255),
                lbl.col,
                &mut occupied,
                false,
            );
        }
    }

    if let Some((n_idx, center, rad, _, _)) = hovered_node {
        let n = &scene.nodes_3d[n_idx];
        painter.circle_stroke(
            center,
            rad + 4.0,
            egui::Stroke::new(2.0_f32, egui::Color32::from_rgb(251, 191, 36)),
        );
        let l1 = if n.label.is_empty() {
            format!("3D Node #{n_idx}")
        } else {
            format!("{} (Node #{n_idx})", n.label)
        };
        let l2 = format!(
            "pos=({:.2}, {:.2}, {:.2})  r={:.2}  glow={:.2}",
            n.pos.x, n.pos.y, n.pos.z, n.radius, n.glow
        );
        if let Some(hp) = hover_pos {
            paint_hover_tooltip(
                painter,
                rect,
                hp,
                &[&l1, &l2],
                egui::Color32::from_rgb(251, 191, 36),
            );
        }
    }
}

fn paint_scene_2d(
    painter: &egui::Painter,
    rect: egui::Rect,
    scene: &StructureScene,
    show_labels: bool,
    hover_pos: Option<egui::Pos2>,
) {
    let painter = &painter.with_clip_rect(rect);
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
            let col = egui::Color32::from_rgba_unmultiplied(r, g, b, 195);
            let stroke = egui::Stroke::new(e.width as f32, col);
            let delta = p2 - p1;
            let len = delta.length();
            if e.directed && len > 6.0 {
                let u = delta / len;
                let perp = egui::vec2(-u.y, u.x);
                let rad_b = ((scene.nodes_2d[e.b].radius as f32) * scale).clamp(4.0, 12.0);
                let inset = (rad_b + 1.5).min(len * 0.45);
                let tip = p2 - u * inset;
                let arrow_len = (7.0 + (e.width as f32) * 1.1).min((len - inset).max(4.0) * 0.65);
                let half_w = arrow_len * 0.44;
                let shaft_end = tip - u * (arrow_len * 0.72);
                painter.line_segment([p1, shaft_end], stroke);
                let w1 = tip - u * arrow_len + perp * half_w;
                let w2 = tip - u * arrow_len - perp * half_w;
                painter.add(egui::Shape::convex_polygon(
                    vec![tip, w1, w2],
                    col,
                    egui::Stroke::NONE,
                ));
            } else {
                painter.line_segment([p1, p2], stroke);
            }
        }
    }

    struct Label2D {
        pos: egui::Pos2,
        rad: f32,
        text: String,
        col: egui::Color32,
    }
    let mut labels_2d: Vec<Label2D> = Vec::new();
    let mut hovered_node: Option<(usize, egui::Pos2, f32, f32)> = None;

    for (idx, n) in scene.nodes_2d.iter().enumerate() {
        let p = to_screen(n.x, n.y);
        let rad = ((n.radius as f32) * scale).clamp(4.0, 12.0);
        let [r, g, b] = n.color;
        let col = egui::Color32::from_rgb(r, g, b);
        painter.circle_filled(p, rad, col);
        painter.circle_stroke(
            p,
            rad,
            egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(10, 14, 26)),
        );

        if let Some(hp) = hover_pos {
            let d = p.distance(hp);
            let hit_r = (rad + 7.0).max(11.0);
            if d <= hit_r && hovered_node.map_or(true, |(_, _, _, best_d)| d < best_d) {
                hovered_node = Some((idx, p, rad, d));
            }
        }

        if show_labels && !n.label.is_empty() {
            labels_2d.push(Label2D {
                pos: p,
                rad,
                text: n.label.clone(),
                col,
            });
        }
    }

    if show_labels {
        let mut occupied: Vec<egui::Rect> = Vec::new();
        for lbl in labels_2d {
            paint_pill_label(
                painter,
                rect,
                lbl.pos + egui::vec2(0.0, -lbl.rad - 5.0),
                egui::Align2::CENTER_BOTTOM,
                &lbl.text,
                egui::FontId::monospace(9.5),
                egui::Color32::from_rgb(228, 238, 255),
                lbl.col,
                &mut occupied,
                false,
            );
        }
    }

    if let Some((n_idx, p, rad, _)) = hovered_node {
        let n = &scene.nodes_2d[n_idx];
        painter.circle_stroke(
            p,
            rad + 4.0,
            egui::Stroke::new(2.0_f32, egui::Color32::from_rgb(251, 191, 36)),
        );
        let l1 = if n.label.is_empty() {
            format!("2D Node #{n_idx}")
        } else {
            format!("{} (2D Node #{n_idx})", n.label)
        };
        let l2 = format!("xy=({:.2}, {:.2})  r={:.3}", n.x, n.y, n.radius);
        if let Some(hp) = hover_pos {
            paint_hover_tooltip(
                painter,
                rect,
                hp,
                &[&l1, &l2],
                egui::Color32::from_rgb(251, 191, 36),
            );
        }
    }
}

pub fn itme_icon_data() -> egui::IconData {
    const SIZE: usize = 64;
    let mut rgba = vec![0u8; SIZE * SIZE * 4];
    let cx = (SIZE - 1) as f32 * 0.5;
    let cy = (SIZE - 1) as f32 * 0.5;
    let r_max = SIZE as f32 * 0.5;
    let corner = SIZE as f32 * 0.22;
    let half = SIZE as f32 * 0.48;

    let blend = |rgba: &mut [u8], x: usize, y: usize, r: f32, g: f32, b: f32, a: f32| {
        if a <= 0.0 {
            return;
        }
        let idx = (y * SIZE + x) * 4;
        let sa = (a / 255.0).clamp(0.0, 1.0);
        let da = f32::from(rgba[idx + 3]) / 255.0;
        let out_a = sa + da * (1.0 - sa);
        if out_a <= 1e-6 {
            return;
        }
        rgba[idx] = ((r * sa + f32::from(rgba[idx]) * da * (1.0 - sa)) / out_a)
            .clamp(0.0, 255.0) as u8;
        rgba[idx + 1] = ((g * sa + f32::from(rgba[idx + 1]) * da * (1.0 - sa)) / out_a)
            .clamp(0.0, 255.0) as u8;
        rgba[idx + 2] = ((b * sa + f32::from(rgba[idx + 2]) * da * (1.0 - sa)) / out_a)
            .clamp(0.0, 255.0) as u8;
        rgba[idx + 3] = (out_a * 255.0).clamp(0.0, 255.0) as u8;
    };
    let dist_seg = |px: f32, py: f32, x1: f32, y1: f32, x2: f32, y2: f32| -> f32 {
        let (vx, vy) = (x2 - x1, y2 - y1);
        let (wx, wy) = (px - x1, py - y1);
        let l2 = vx * vx + vy * vy;
        let t = if l2 < 1e-9 {
            0.0
        } else {
            ((wx * vx + wy * vy) / l2).clamp(0.0, 1.0)
        };
        (px - (x1 + t * vx)).hypot(py - (y1 + t * vy))
    };

    let root = (0.0_f32, 0.56_f32);
    let mid = (0.0_f32, 0.10_f32);
    let left_b = (-0.36_f32, -0.22_f32);
    let right_b = (0.36_f32, -0.22_f32);
    let top_c = (0.0_f32, -0.46_f32);
    let edges = [
        (root, mid, (56.0_f32, 189.0_f32, 248.0_f32)),
        (mid, left_b, (52.0, 211.0, 153.0)),
        (mid, right_b, (167.0, 139.0, 250.0)),
        (left_b, top_c, (52.0, 211.0, 153.0)),
        (right_b, top_c, (251.0, 191.0, 36.0)),
        (left_b, right_b, (244.0, 114.0, 182.0)),
    ];
    let nodes = [
        (root, 0.13_f32, (56.0_f32, 189.0_f32, 248.0_f32)),
        (mid, 0.14, (251.0, 191.0, 36.0)),
        (left_b, 0.13, (6.0, 214.0, 160.0)),
        (right_b, 0.13, (167.0, 139.0, 250.0)),
        (top_c, 0.15, (130.0, 245.0, 205.0)),
    ];

    for y in 0..SIZE {
        for x in 0..SIZE {
            let dx = (x as f32 - cx).abs() - (half - corner);
            let dy = (y as f32 - cy).abs() - (half - corner);
            let sdf_box = dx.max(0.0).hypot(dy.max(0.0)) + dx.max(dy).min(0.0) - corner;
            if sdf_box > 1.2 {
                continue;
            }
            let bg_alpha = (1.0 - sdf_box).clamp(0.0, 1.0) * 255.0;
            let nx = (x as f32 - cx) / r_max;
            let ny = (y as f32 - cy) / r_max;
            let rad = nx.hypot(ny);
            let glow = (1.0 - rad * 1.1).max(0.0);
            let idx = (y * SIZE + x) * 4;
            rgba[idx] = (11.0 + 14.0 * glow).clamp(0.0, 255.0) as u8;
            rgba[idx + 1] = (18.0 + 28.0 * glow).clamp(0.0, 255.0) as u8;
            rgba[idx + 2] = (34.0 + 48.0 * glow).clamp(0.0, 255.0) as u8;
            rgba[idx + 3] = bg_alpha.clamp(0.0, 255.0) as u8;

            if sdf_box > -(SIZE as f32) * 0.04 && sdf_box <= 0.5 {
                blend(&mut rgba, x, y, 56.0, 189.0, 248.0, 130.0);
            }
            let ring_d = (rad - 0.72).abs();
            let ring_w = 0.055_f32.max(1.4 / r_max);
            if ring_d < ring_w {
                let t = (nx + 1.0) * 0.5;
                let a = (1.0 - ring_d / ring_w).max(0.0) * 235.0;
                blend(
                    &mut rgba,
                    x,
                    y,
                    56.0 * (1.0 - t) + 52.0 * t,
                    189.0 * (1.0 - t) + 211.0 * t,
                    248.0 * (1.0 - t) + 153.0 * t,
                    a,
                );
            }
            let ew = 0.055_f32.max(1.5 / r_max);
            for &(p1, p2, col) in &edges {
                let d = dist_seg(nx, ny, p1.0, p1.1, p2.0, p2.1);
                if d < ew {
                    let a = (1.0 - d / ew).max(0.0) * 230.0;
                    blend(&mut rgba, x, y, col.0, col.1, col.2, a);
                }
            }
            let aa = 0.04_f32.max(1.2 / r_max);
            for &(np, nr, col) in &nodes {
                let d = (nx - np.0).hypot(ny - np.1);
                if d < nr + aa {
                    let a = ((nr + aa - d) / aa).clamp(0.0, 1.0) * 255.0;
                    blend(&mut rgba, x, y, col.0, col.1, col.2, a);
                    if d < nr * 0.45 {
                        blend(&mut rgba, x, y, 245.0, 250.0, 255.0, 180.0);
                    }
                }
            }
        }
    }

    egui::IconData {
        rgba,
        width: SIZE as u32,
        height: SIZE as u32,
    }
}

pub fn run_gui(seed: u64) -> Result<(), String> {
    let viewport = egui::ViewportBuilder::default()
        .with_title(
            "ITME v23.0 — Domain-Isolated Biosphere, Hominid & AI Evolutionary Simulator",
        )
        .with_icon(Arc::new(itme_icon_data()))
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
        4 => {
            let mut session = DomainSession::new_coupled_biosphere(
                RegimeId::R35Symbiogenetic,
                SimulationMode::OrganicSandbox,
                32,
                seed + 404,
            );
            let epochs = session.world.run_continuous_grand_trajectory(10);
            session.grand_epochs = epochs;
            for _ in 0..8 {
                session.step_once();
            }
            session.viewport_layout = ViewportLayout::GrandTrajectory;
            session.paused = true;
            app.screen = AppScreen::Simulation(Box::new(session));
        }
        other => {
            return Err(format!(
                "Unsupported --screenshot-view {other} (expected 1 = Main Menu, 2 = Molecular Domain, 3 = AI Domain, 4 = Coupled PlanetaryBiosphere)"
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
    let file = fs::File::create(path)
        .map_err(|e| format!("Failed creating PNG '{}': {e}", path.display()))?;
    let mut encoder = png::Encoder::new(BufWriter::new(file), width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_compression(png::Compression::High);
    let mut writer = encoder
        .write_header()
        .map_err(|e| format!("Failed writing PNG header '{}': {e}", path.display()))?;
    writer
        .write_image_data(rgba)
        .map_err(|e| format!("Failed writing PNG data '{}': {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pill_label_overlap_culling_and_clamping_prevent_collisions() {
        let bounds = egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(320.0, 200.0));
        let out_of_bounds = egui::Rect::from_min_size(egui::pos2(305.0, -12.0), egui::vec2(60.0, 18.0));
        let clamped = clamp_pill_rect(out_of_bounds, bounds);
        assert!(clamped.left() >= bounds.left());
        assert!(clamped.right() <= bounds.right() + 1e-4);
        assert!(clamped.top() >= bounds.top());
        assert!(clamped.bottom() <= bounds.bottom() + 1e-4);

        let ctx = egui::Context::default();
        let mut raw = egui::RawInput::default();
        raw.screen_rect = Some(bounds);
        let mut placed_count = 0usize;
        let _ = ctx.run(raw, |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let painter = ui.painter();
                let mut occupied = Vec::new();
                // Place multiple identical-anchor pills: first succeeds, up/down nudges succeed, 4th+ is culled.
                for idx in 0..5 {
                    let ok = paint_pill_label(
                        painter,
                        bounds,
                        egui::pos2(160.0, 100.0),
                        egui::Align2::CENTER_CENTER,
                        &format!("NodeLabel_{idx}"),
                        egui::FontId::monospace(10.0),
                        egui::Color32::WHITE,
                        egui::Color32::LIGHT_BLUE,
                        &mut occupied,
                        false,
                    );
                    if ok {
                        placed_count += 1;
                    }
                }
                assert_eq!(occupied.len(), placed_count);
                assert!(placed_count <= 3, "only base + up + down slots fit without AABB overlap");
                for i in 0..occupied.len() {
                    assert!(bounds.contains_rect(occupied[i]));
                    for j in (i + 1)..occupied.len() {
                        assert!(
                            !occupied[i].intersects(occupied[j]),
                            "placed pills {i} and {j} must not overlap"
                        );
                    }
                }
            });
        });
    }

    #[test]
    fn main_menu_seed_parsing_handles_decimal_hex_and_invalid() {
        let mut menu = MainMenuState::new(20250813);
        assert_eq!(menu.try_parse_seed(), Ok(20250813));
        assert_eq!(menu.parsed_seed(), 20250813);

        menu.seed_text = "  424242  ".to_string();
        assert_eq!(menu.try_parse_seed(), Ok(424242));
        assert_eq!(menu.parsed_seed(), 424242);

        menu.seed_text = "0x1A2B_deadbeef".replace('_', "");
        assert_eq!(menu.try_parse_seed(), Ok(0x1A2B_deadbeef));
        assert_eq!(menu.parsed_seed(), 0x1A2B_deadbeef);

        menu.seed_text = "0XFF".to_string();
        assert_eq!(menu.try_parse_seed(), Ok(255));

        menu.seed_text = "0x".to_string();
        assert!(menu.try_parse_seed().is_err());
        assert_eq!(menu.parsed_seed(), 20250813);

        menu.seed_text = "not_a_number".to_string();
        assert!(menu.try_parse_seed().is_err());
        assert_eq!(menu.parsed_seed(), 20250813);

        menu.seed_text = "   ".to_string();
        assert!(menu.try_parse_seed().is_err());
        assert_eq!(menu.parsed_seed(), 20250813);
    }

    #[test]
    fn gui_f1_and_escape_toggle_about_modal_headlessly() {
        let ctx = egui::Context::default();
        let mut app = ItmeGuiApp::new_with_ctx(&ctx, 20250813);
        assert!(!app.show_about);

        let run_key = |app: &mut ItmeGuiApp, key: egui::Key| {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1440.0, 900.0),
            ));
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            });
            let _ = ctx.run(raw, |c| app.draw_panels(c));
        };

        run_key(&mut app, egui::Key::F1);
        assert!(app.show_about, "F1 should open About ITME modal");

        run_key(&mut app, egui::Key::Escape);
        assert!(!app.show_about, "Escape should close About ITME modal");

        run_key(&mut app, egui::Key::F1);
        assert!(app.show_about);

        run_key(&mut app, egui::Key::F1);
        assert!(!app.show_about);
    }

    #[test]
    fn gui_global_hotkeys_control_playback_camera_layouts_and_escape_with_input_guard() {
        let ctx = egui::Context::default();
        let mut app = ItmeGuiApp::new_with_ctx(&ctx, 20250813);

        let session = DomainSession::new_coupled_biosphere(
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            32,
            20250813,
        );
        app.screen = AppScreen::Simulation(Box::new(session));

        // Verify hotkeys are suppressed while a TextEdit has keyboard focus (`ctx.wants_keyboard_input()`)
        let mut focus_buf = String::from("seed_edit");
        let step_initial = if let AppScreen::Simulation(ref s) = app.screen {
            s.world.step_count
        } else {
            0
        };
        for key in [egui::Key::Space, egui::Key::Period, egui::Key::F1, egui::Key::Escape] {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1440.0, 900.0),
            ));
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            });
            let _ = ctx.run(raw, |c| {
                egui::Area::new(egui::Id::new("test_focus_area")).show(c, |ui| {
                    let resp = ui.text_edit_singleline(&mut focus_buf);
                    resp.request_focus();
                });
                assert!(c.wants_keyboard_input());
                app.draw_panels(c);
            });
        }
        assert!(!app.show_about, "F1 must be ignored while text input has focus");
        if let AppScreen::Simulation(ref s) = app.screen {
            assert!(!s.paused, "Space must be ignored while text input has focus");
            assert_eq!(s.world.step_count, step_initial, "'.' must be ignored while text input has focus");
        } else {
            panic!("Esc must be ignored while text input has focus");
        }

        // Flush one frame without the temporary TextEdit so egui clears text-input focus
        let mut clear_raw = egui::RawInput::default();
        clear_raw.screen_rect = Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(1440.0, 900.0),
        ));
        let _ = ctx.run(clear_raw, |c| app.draw_panels(c));

        let send_key = |app: &mut ItmeGuiApp, key: egui::Key| {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1440.0, 900.0),
            ));
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            });
            let _ = ctx.run(raw, |c| app.draw_panels(c));
        };

        // Space toggles paused
        send_key(&mut app, egui::Key::Space);
        if let AppScreen::Simulation(ref s) = app.screen {
            assert!(s.paused, "Space should pause running simulation");
        } else {
            panic!("Expected Simulation screen");
        }

        // '.' advances 1 step
        let step_before = if let AppScreen::Simulation(ref s) = app.screen {
            s.world.step_count
        } else {
            0
        };
        send_key(&mut app, egui::Key::Period);
        if let AppScreen::Simulation(ref s) = app.screen {
            assert_eq!(s.world.step_count, step_before + 1, "'.' should advance one step");
        }

        // 'R' resets modified camera to Camera3D::default()
        if let AppScreen::Simulation(ref mut s) = app.screen {
            s.camera.yaw = 2.4;
            s.camera.pitch = -0.7;
            s.camera.zoom = 2.1;
        }
        send_key(&mut app, egui::Key::R);
        let def_cam = Camera3D::default();
        if let AppScreen::Simulation(ref s) = app.screen {
            assert!((s.camera.yaw - def_cam.yaw).abs() < 1e-12);
            assert!((s.camera.pitch - def_cam.pitch).abs() < 1e-12);
            assert!((s.camera.zoom - def_cam.zoom).abs() < 1e-12);
        }

        // '1'..'6' switches viewport layouts
        let expected = [
            (egui::Key::Num2, ViewportLayout::SpatialGridAnd3D),
            (egui::Key::Num3, ViewportLayout::PhylogenyFull),
            (egui::Key::Num4, ViewportLayout::Perspective3D),
            (egui::Key::Num5, ViewportLayout::Blueprint2D),
            (egui::Key::Num6, ViewportLayout::GrandTrajectory),
            (egui::Key::Num1, ViewportLayout::SpatialAndPhylogeny),
        ];
        for (key, layout) in expected {
            send_key(&mut app, key);
            if let AppScreen::Simulation(ref s) = app.screen {
                assert_eq!(s.viewport_layout, layout);
            }
        }

        // Esc closes About modal first without leaving simulation, then Esc returns to MainMenu
        send_key(&mut app, egui::Key::F1);
        assert!(app.show_about);
        send_key(&mut app, egui::Key::Escape);
        assert!(!app.show_about);
        assert!(matches!(app.screen, AppScreen::Simulation(_)));

        send_key(&mut app, egui::Key::Escape);
        assert!(matches!(app.screen, AppScreen::MainMenu(_)));
    }

    #[test]
    fn gui_telemetry_onboarding_crosshair_and_3d_2d_render_headlessly() {
        let ctx = egui::Context::default();
        let mut session = DomainSession::new_coupled_biosphere(
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            32,
            202607,
        );
        // Simulate fresh world before 2nd history sample (onboarding hint path)
        session.world.history.truncate(1);
        session.paused = true;
        session.show_milestones = false;
        let mut app = ItmeGuiApp::new_with_ctx(&ctx, 202607);
        app.screen = AppScreen::Simulation(Box::new(session));

        let run_frame = |app: &mut ItmeGuiApp, pointer: Option<egui::Pos2>| {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1440.0, 900.0),
            ));
            if let Some(pos) = pointer {
                raw.events.push(egui::Event::PointerMoved(pos));
            }
            let _ = ctx.run(raw, |c| app.draw_panels(c));
        };

        run_frame(&mut app, None);

        // Advance steps and test hover crosshair on bottom telemetry plots + 3D/2D renderers
        if let AppScreen::Simulation(ref mut s) = app.screen {
            for _ in 0..8 {
                s.step_once();
            }
            s.viewport_layout = ViewportLayout::Perspective3D;
        }
        run_frame(&mut app, Some(egui::pos2(500.0, 820.0)));

        if let AppScreen::Simulation(ref mut s) = app.screen {
            s.viewport_layout = ViewportLayout::Blueprint2D;
        }
        run_frame(&mut app, Some(egui::pos2(1020.0, 820.0)));
    }
}
