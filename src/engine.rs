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

//! ITME v12.0 (`v1.2.0`) Anaerobic Mitochondrial Reduction (`Mito->Hydro->MitoS->MRO->0`),
//! Reductive Evolution of Photosynthetic Holobionts into Apicoplasts (`Holo->Apico` & `Apico->0`),
//! Mixotrophic Kleptoplasty, Shopping-Bag `HGT` Tipping Points, Serial `1°/2°/3°` Endosymbiosis,
//! Periplastidial Nucleomorph Reduction & Reticulate Phylogeny Engine (`SimulationWorld`).

use crate::geometry::{build_structure_scene, render_stage_atlas_svg, StructureScene};
use crate::model::{
    compute_critical_fidelity, compute_food_web_entropy, compute_heritable_mutual_info,
    compute_thermodynamic_hazard, AbiogenesisStage, AgricultureState, ChiralCodeReactor,
    ChiralityCodeStage, CulturalState, EnvironmentControls, EusocialState, ExosomaticState,
    GenomeSubstrateStage, GrandEpochStage, GrandTrajectoryEpochSummary, HadeanGeochemicalBath,
    ItmeMetrics, MachineSynthesisState, MroClass, MulticellularState, OrderPhase,
    PlanetaryBiosphere, PlanetaryState, PriceDecomposition, RnrVirosphereReactor, Rng,
    TelemetryWindow, TrophicGuild,
};
use crate::organisms::{ModeledUnit, OrganismKind};
use crate::phylogeny::{ExtantPhyloState, PhyloSummary, PhylogenyTracker, PhylogramLayout};
use crate::regimes::{classify_emergent_state, RegimeId};
use crate::spatial::{
    analyze_spatial_speciation, normalized_sequence_divergence, BiomeZone, SpatialGrid,
    DEFAULT_GRID_DIM,
};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SimulationMode {
    GuidedCurriculum,
    OrganicSandbox,
}

#[derive(Clone, Debug)]
pub struct SimulationWorld {
    pub organism: OrganismKind,
    pub mode: SimulationMode,
    pub target_regime: RegimeId,
    pub classified_regime: RegimeId,
    pub controls: EnvironmentControls,
    pub grid: SpatialGrid,
    pub phylogeny: PhylogenyTracker,
    pub population: Vec<ModeledUnit>,
    pub step_count: usize,
    pub stage_step: usize,
    pub steps_per_stage: usize,
    pub next_unit_id: u64,
    pub rng: Rng,
    pub telemetry_window: TelemetryWindow,
    pub env_phase: f64,
    pub smoothed_birth_rate: f64,
    pub smoothed_plastic_signal: f64,
    pub smoothed_assimilation: f64,
    pub smoothed_recombination: f64,
    pub smoothed_mut_diversity: f64,
    pub last_observed_fidelity: Option<f64>,
    pub last_price: PriceDecomposition,
    pub last_adami_rho: f64,
    pub last_species_count: usize,
    pub last_fst: f64,
    pub last_trophic_flux: f64,
    pub last_cross_feed_flux: f64,
    pub total_apoptotic_culls: usize,
    pub total_neoplastic_collapses: usize,
    pub total_policed_eggs: usize,
    pub total_colony_anarchy_collapses: usize,
    pub total_antimicrobial_shields: usize,
    pub total_garden_blight_collapses: usize,
    pub total_cross_deme_recombinations: usize,
    pub total_cultural_drift_collapses: usize,
    pub total_peer_retraction_audits: usize,
    pub total_archive_corruption_collapses: usize,
    pub total_formal_proofs_verified: usize,
    pub total_recursive_model_collapses: usize,
    pub biosphere: PlanetaryBiosphere,
    pub coupled_biosphere_mode: bool,
    pub total_planetary_burnout_collapses: usize,
    pub hadean_bath: HadeanGeochemicalBath,
    pub total_eigen_error_catastrophes: usize,
    pub total_mg_vesicle_collapses: usize,
    pub chiral_reactor: ChiralCodeReactor,
    pub total_racemic_cross_inhibitions: usize,
    pub total_mistranslation_babel_collapses: usize,
    pub rnr_virosphere: RnrVirosphereReactor,
    pub total_cytosine_deamination_meltdowns: usize,
    pub total_lytic_viral_collapses: usize,
    pub history: Vec<ItmeMetrics>,
    pub regime_history: Vec<RegimeId>,
    pub phase_time: f64,
}

impl SimulationWorld {
    pub fn new(
        organism: OrganismKind,
        initial_regime: RegimeId,
        mode: SimulationMode,
        steps_per_stage: usize,
        seed: u64,
    ) -> Self {
        let mut rng = Rng::from_seed(seed ^ ((organism as u64 + 1) * 0x517C_C1B7_2722_0A95));
        let controls = initial_regime.preset_controls();
        let grid = SpatialGrid::new(
            DEFAULT_GRID_DIM,
            DEFAULT_GRID_DIM,
            controls.bath_exergy_flux,
        );
        let mut world = Self {
            organism,
            mode,
            target_regime: initial_regime,
            classified_regime: initial_regime,
            controls: controls.clone(),
            grid,
            phylogeny: PhylogenyTracker::new(),
            population: Vec::new(),
            step_count: 0,
            stage_step: 0,
            steps_per_stage: steps_per_stage.max(8),
            next_unit_id: 1,
            rng: rng.clone(),
            telemetry_window: TelemetryWindow::new(12),
            env_phase: 0.4,
            smoothed_birth_rate: 0.0,
            smoothed_plastic_signal: 0.0,
            smoothed_assimilation: 0.0,
            smoothed_recombination: 0.0,
            smoothed_mut_diversity: 0.0,
            last_observed_fidelity: None,
            last_price: PriceDecomposition::default(),
            last_adami_rho: 1.0,
            last_species_count: 1,
            last_fst: 0.0,
            last_trophic_flux: 0.0,
            last_cross_feed_flux: 0.0,
            total_apoptotic_culls: 0,
            total_neoplastic_collapses: 0,
            total_policed_eggs: 0,
            total_colony_anarchy_collapses: 0,
            total_antimicrobial_shields: 0,
            total_garden_blight_collapses: 0,
            total_cross_deme_recombinations: 0,
            total_cultural_drift_collapses: 0,
            total_peer_retraction_audits: 0,
            total_archive_corruption_collapses: 0,
            total_formal_proofs_verified: 0,
            total_recursive_model_collapses: 0,
            biosphere: PlanetaryBiosphere::new_for_controls(&controls),
            coupled_biosphere_mode: false,
            total_planetary_burnout_collapses: 0,
            hadean_bath: HadeanGeochemicalBath::new_for_controls(&controls),
            total_eigen_error_catastrophes: 0,
            total_mg_vesicle_collapses: 0,
            chiral_reactor: ChiralCodeReactor::new_for_controls(&controls),
            total_racemic_cross_inhibitions: 0,
            total_mistranslation_babel_collapses: 0,
            rnr_virosphere: RnrVirosphereReactor::new_for_controls(&controls),
            total_cytosine_deamination_meltdowns: 0,
            total_lytic_viral_collapses: 0,
            history: Vec::with_capacity(780),
            regime_history: Vec::with_capacity(780),
            phase_time: 0.0,
        };
        world.populate_for_regime(initial_regime, &mut rng);
        world.rng = rng;
        world.seed_observable_estimators();
        world.reset_phylogeny_from_population();
        world.record_metrics();
        world
    }

    fn alloc_id(&mut self) -> u64 {
        let id = self.next_unit_id;
        self.next_unit_id += 1;
        id
    }

    fn reset_phylogeny_from_population(&mut self) {
        let founders: Vec<(u64, usize, TrophicGuild, BiomeZone, f64)> = self
            .population
            .iter()
            .map(|u| {
                let b = self.grid.biome_at(u.grid_x, u.grid_y);
                let lev = u.complexity_c / u.structure_s.max(2.0);
                (u.id, u.species_id, u.guild, b, lev)
            })
            .collect();
        self.phylogeny.reset_with_founders(self.step_count, &founders);
    }

    fn seed_observable_estimators(&mut self) {
        let ex_on = self.controls.bath_exergy_flux > 0.03;
        self.smoothed_birth_rate = if ex_on {
            self.controls.replication_rate * 0.75
        } else {
            0.0
        };
        self.smoothed_plastic_signal = if ex_on {
            self.controls.memory_plasticity * 0.45
        } else {
            0.0
        };
        self.smoothed_assimilation = if ex_on && self.controls.replication_rate >= 0.08 {
            self.controls.soft_inheritance_coupling * 0.45
        } else {
            0.0
        };
        self.smoothed_recombination = if ex_on && self.controls.replication_rate >= 0.08 {
            (self.controls.recombination_rate
                + 0.45 * self.controls.endosymbiosis_rate
                + 0.30 * self.controls.serial_endosymbiosis_rate
                + 0.30 * self.controls.klepto_hgt_rate)
                * 0.45
        } else {
            0.0
        };
        self.smoothed_mut_diversity = if ex_on && self.controls.replication_rate >= 0.08 {
            self.controls.heritable_variation
        } else {
            0.0
        };
        self.last_price = PriceDecomposition::default();
        self.last_trophic_flux = 0.0;
        self.last_cross_feed_flux = 0.0;
    }

    pub fn trigger_environmental_shock(&mut self) {
        self.env_phase += 1.55;
        for e in &mut self.grid.exergy {
            *e *= 0.45;
        }
    }

    pub fn trigger_apex_predator_cull(&mut self) {
        self.population
            .retain(|u| u.guild != TrophicGuild::PredatorParasite);
        if self.population.is_empty() {
            let mut rng = self.rng.clone();
            self.populate_for_regime(self.target_regime, &mut rng);
            self.rng = rng;
            self.reset_phylogeny_from_population();
        }
    }

    /// Manually triggers one step of **Anaerobic Mitochondrial Reduction (`Mito -> Hydro -> MitoS -> MRO->0`)**
    /// across extant MRO-bearing lineages, recording phylogenetic transition badges (`Mito→Hydro`,
    /// `Hydro→MitoS`, `MRO→0` when `cytosolic_suf_hgt_rescue >= 0.65`) in `PhylogenyTracker`!
    /// Returns `(hydrogenosomes_formed, mitosomes_formed, amiteukaryote_mro_losses)`.
    /// Manually triggers one step of **Clonal Multicellularity & Soma-Germline Differentiation**
    /// (`Unicellular Holobiont -> UndifferentiatedClonal -> GermSomaMulticellular`) plus a pulse
    /// of **mitochondrial apoptotic cheater policing**, recording `Holo→Clonal` and `Clonal→GermSoma`
    /// transitions in `PhylogenyTracker`!
    /// Manually advances `GermSomaMulticellular` lineages (promoting unicellular/clonal holobionts first
    /// if needed) along the **Eusocial Superorganism ETI** (`Solitary -> Subsocial -> FacultativeEusocial ->
    /// ObligateSuperorganism` when `queen_monogamy_index >= 0.75`) and executes a pulse of **worker kin-policing**!
    /// Manually advances eusocial superorganisms (promoting holobionts to `GermSoma -> ObligateSuperorganism`
    /// first if needed when `queen_monogamy_index >= 0.75`) along the **Agricultural Ectosymbiosis Continuum**
    /// (`HunterForager -> LowerFacultativeFarming -> HigherDomesticatedCrop -> LeafCutterAgroEcosystem`) and
    /// runs cuticular *Pseudonocardia* antimicrobial suppression against *Escovopsis*!
    /// Manually advances symbolic cultural populations along **Externalized Exosomatic Memory &
    /// Institutional Peer-Verification** (`Endosomatic -> Notation -> WrittenArch -> VerifiedInst`)
    /// and executes institutional replication audits purging false discoveries!
    /// Manually advances peer-verified archival institutions along **Autonomous Machine-Verified
    /// Exosomatic Synthesis** (`ManualArchival -> AlgorithmicMediation -> MachineVerifiedSymbiosis`)
    /// and runs formal proof verification + empirical-anchor safeguards!
    /// Manually co-evolves autotrophic, syntrophic, superorganismal, and machine-verified guilds
    /// into **Gaia-2.0 Planetary Homeostatic Awakening**, closing `O2 / CO2 / N-P` cycles and
    /// recording `Archean -> GOE` and `Techno -> Gaia2Homeostasis` landmarks!
    /// Triggers **De-Novo Hadean Abiogenesis (`~4.2 Ga`)**, seeding prebiotic monomer pools in
    /// coupled hot-spring wet-dry + submarine alkaline `Fe(Ni)S` cradles and stepping across
    /// `RAF Protocell -> Ribozyme Hypercycle -> Communal Progenote -> LUCA`!
    /// Triggers **Prebiotic Homochiral Symmetry Breaking (`Fe3O4` CISS + Frank), Stereoselective
    /// `3'-NCCA` Minihelix `L-Aminoacylation`, and Collective `HGT` Standard Genetic Code (`SGC`) Crystallization**!
    /// Triggers **Anaerobic Radical-SAM `RNR` `RNA -> U-DNA -> T-DNA` Transition, `ThyX + UDG`
    /// Cytosine Deamination Repair, and Viral Polymerase Domestication** into Megabase Genomes!
    /// Executes the **Unbroken 8-Epoch 4.25 Ga Hadean-to-Gaia Continuous Grand Trajectory**
    /// on a single persistent `SimulationWorld` and `PhylogenyTracker`!
    pub fn run_continuous_grand_trajectory(
        &mut self,
        ticks_per_epoch: usize,
    ) -> Vec<GrandTrajectoryEpochSummary> {
        let t_ep = ticks_per_epoch.max(4);
        let mut summaries = Vec::with_capacity(8);

        // Configure anoxic Hadean starting atmosphere & prebiotic reactors (4.25 Ga):
        self.coupled_biosphere_mode = true;
        self.controls = RegimeId::R35Symbiogenetic.preset_controls();
        self.controls.ambient_oxygen_level = 0.06;
        self.controls.biogeochemical_feedback_gain = 0.0;
        self.controls.chiral_symmetry_breaking_mode = true;
        self.controls.hadean_abiogenesis_mode = true;
        self.controls.rna_to_dna_virosphere_mode = true;
        self.controls.magnetite_ciss_spin_bias = 0.94;
        self.controls.minihelix_chiral_selectivity = 0.90;
        self.controls.communal_code_hgt_rate = 0.30;
        self.controls.wet_dry_cycle_amplitude = 0.62;
        self.controls.hydrothermal_delta_ph = 2.80;
        self.controls.citrate_mg_chelation = 0.88;
        self.controls.ribozyme_polymerase_fidelity = 0.62;
        self.controls.rnr_thiyl_radical_activity = 0.22;
        self.controls.thymidylate_udg_repair_efficiency = 0.50;
        self.controls.virosphere_arms_race_coupling = 0.72;
        self.controls.lytic_viral_burst_pressure = 0.18;

        self.biosphere = PlanetaryBiosphere::new_for_controls(&self.controls);
        self.chiral_reactor = ChiralCodeReactor::new_for_controls(&self.controls);
        self.hadean_bath = HadeanGeochemicalBath::new_for_controls(&self.controls);
        self.rnr_virosphere = RnrVirosphereReactor::new_for_controls(&self.controls);

        for u in &mut self.population {
            u.seed_hadean_prebiotic_pool(&self.hadean_bath);
            u.seed_racemic_prebiotic_pool(&self.chiral_reactor);
            u.seed_labile_rna_cell(&self.rnr_virosphere);
        }

        for &ep in &GrandEpochStage::ALL {
            let landmark = match ep {
                GrandEpochStage::Epoch1HadeanHomochirality => {
                    // Break racemic symmetry on magnetite Fe3O4 CISS + 3'-NCCA L-aminoacylation:
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!(
                        "Magnetite CISS ee_D={:.0}% + Minihelix ee_L={:.0}%",
                        self.chiral_reactor.ribose_ee_d * 100.0,
                        self.chiral_reactor.amino_ee_l * 100.0
                    )
                }
                GrandEpochStage::Epoch2SgcAndRafProtocells => {
                    self.controls.communal_code_hgt_rate = 0.92;
                    self.controls.wet_dry_cycle_amplitude = 0.88;
                    self.controls.hydrothermal_delta_ph = 3.45;
                    self.controls.ribozyme_polymerase_fidelity = 0.84;
                    let (_, sgc_n, mp) = self.trigger_homochiral_sgc_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Universal 20-AA SGC (n={sgc_n}, MP={:.1}%) + RAF Closure", mp * 100.0)
                }
                GrandEpochStage::Epoch3RnaToDnaAndLuca => {
                    self.controls.ribozyme_polymerase_fidelity = 0.96;
                    self.controls.rnr_thiyl_radical_activity = 0.94;
                    self.controls.thymidylate_udg_repair_efficiency = 0.94;
                    self.controls.virosphere_arms_race_coupling = 0.90;
                    let (_, luca_n, _) = self.trigger_hadean_abiogenesis_pulse();
                    let (_, tdna_n, kb) = self.trigger_rna_to_dna_virosphere_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Darwinian LUCA (n={luca_n}) + RNR/UDG T-DNA ({kb:.0} kb, n={tdna_n})")
                }
                GrandEpochStage::Epoch4GoeAndSymbiogenesis => {
                    self.controls.chiral_symmetry_breaking_mode = false;
                    self.controls.hadean_abiogenesis_mode = false;
                    self.controls.rna_to_dna_virosphere_mode = false;
                    self.controls.biogeochemical_feedback_gain = 0.72;
                    for _ in 0..8 {
                        let _ = self.trigger_endosymbiotic_pulse();
                    }
                    let _ = self.trigger_serial_endosymbiosis_pulse();
                    for _ in 0..4 {
                        let _ = self.trigger_egt_transfer_pulse();
                    }
                    for _ in 0..8 {
                        self.biosphere
                            .step_biogeochemical_feedbacks(1.25, 1.15, 0.45, 0.0, 0.0, &self.controls);
                    }
                    self.controls.ambient_oxygen_level = self.biosphere.atmospheric_o2;
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!(
                        "GOE Oxygenation (O2={:.2}) + Asgard-Mito/Plastid LECA Holobionts",
                        self.biosphere.atmospheric_o2
                    )
                }
                GrandEpochStage::Epoch5ClonalMulticellularity => {
                    self.controls.single_cell_bottleneck = true;
                    self.controls.apoptosis_policing_rate = 0.88;
                    let _ = self.trigger_multicellular_transition_pulse();
                    let (_, gs_n, _) = self.trigger_multicellular_transition_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Clonal r=1.0 Zygote + Germ-Soma Multicellularity (n={gs_n}) + Apoptosis")
                }
                GrandEpochStage::Epoch6SuperorganismAgriculture => {
                    self.controls.queen_monogamy_index = 0.92;
                    self.controls.worker_policing_rate = 0.90;
                    self.controls.defensive_symbiont_potency = 0.90;
                    let _ = self.trigger_superorganism_transition_pulse();
                    let (_, sup_n, _) = self.trigger_superorganism_transition_pulse();
                    let _ = self.trigger_agriculture_domestication_pulse();
                    let _ = self.trigger_agriculture_domestication_pulse();
                    let (_, ag_n, _) = self.trigger_agriculture_domestication_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Kin-Policed Superorganisms (n={sup_n}) + 8n Clonal Crop Symbiosis (n={ag_n})")
                }
                GrandEpochStage::Epoch7SymbolicExosomaticArchives => {
                    self.controls.symbolic_teaching_fidelity = 0.94;
                    self.controls.inter_deme_connectivity = 0.36;
                    self.controls.peer_verification_rigor = 0.90;
                    self.controls.formal_verifier_strictness = 0.92;
                    self.controls.empirical_grounding_floor = 0.35;
                    self.controls.synthetic_generation_rate = 0.72;
                    for _ in 0..3 {
                        let _ = self.trigger_cultural_ratchet_pulse();
                    }
                    for _ in 0..3 {
                        let _ = self.trigger_exosomatic_archive_pulse();
                    }
                    let _ = self.trigger_machine_synthesis_pulse();
                    let (_, fm_n, _) = self.trigger_machine_synthesis_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Symbolic Culture + Verified Written Archives + Formal AI Symbiosis (n={fm_n})")
                }
                GrandEpochStage::Epoch8Gaia2PlanetaryHomeostasis => {
                    self.controls.biogeochemical_feedback_gain = 0.85;
                    self.controls.planetary_stewardship_coupling = 0.92;
                    self.controls.extractive_burn_intensity = 0.16;
                    let (_, aw_n, ex_b) = self.trigger_planetary_awakening_pulse();
                    for _ in 0..t_ep {
                        self.step();
                    }
                    format!("Gaia-2.0 Homeostatic Awakening (n={aw_n}, Exergy={ex_b:.2})")
                }
            };

            let m = self.latest_metrics().clone();

            summaries.push(GrandTrajectoryEpochSummary {
                epoch: ep,
                tick_end: self.step_count,
                geologic_age_ga: ep.geologic_age_ga(),
                population_size: m.population_size,
                heritable_channels_active: ep.active_heritable_channels(),
                heritable_mutual_info_bits: m.heritable_mutual_info_bits,
                informational_leverage_li: m.leverage_li,
                order_velocity_vo: m.order_velocity_vo,
                planetary_o2: self.biosphere.atmospheric_o2,
                planetary_co2_ppm: self.biosphere.atmospheric_co2_ppm,
                planetary_exergy_budget: self.biosphere.planetary_exergy_budget,
                landmark_summary: landmark,
            });
        }

        summaries
    }

    pub fn trigger_rna_to_dna_virosphere_pulse(&mut self) -> (usize, usize, f64) {
        self.controls.rna_to_dna_virosphere_mode = true;
        self.controls.rnr_thiyl_radical_activity =
            self.controls.rnr_thiyl_radical_activity.max(0.92);
        self.controls.thymidylate_udg_repair_efficiency =
            self.controls.thymidylate_udg_repair_efficiency.max(0.94);
        self.controls.virosphere_arms_race_coupling =
            self.controls.virosphere_arms_race_coupling.max(0.88);
        self.controls.lytic_viral_burst_pressure =
            self.controls.lytic_viral_burst_pressure.min(0.28);
        self.rnr_virosphere = RnrVirosphereReactor::new_for_controls(&self.controls);

        let mut rnr_cnt = 0usize;
        let mut tdna_cnt = 0usize;
        let mut transitions = Vec::new();
        for u in &mut self.population {
            let old_id = u.id;
            u.seed_labile_rna_cell(&self.rnr_virosphere);
            let mut crossed_rnr = false;
            for _ in 0..8 {
                self.rnr_virosphere.step_rnr_and_virosphere_cycle(&self.controls);
                let (st, _, _) = u
                    .step_rnr_dna_and_virosphere_coevolution(&self.rnr_virosphere, &mut self.rng);
                if st.has_deoxyribose_rnr() && !crossed_rnr {
                    crossed_rnr = true;
                    rnr_cnt += 1;
                    let mid_id = self.next_unit_id;
                    self.next_unit_id += 1;
                    let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                    let lev = (u.complexity_c / u.structure_s.max(2.0)).max(1.8);
                    transitions.push((
                        mid_id,
                        old_id,
                        u.species_id,
                        u.guild,
                        biome,
                        GenomeSubstrateStage::ViralModifiedUDna,
                        u.chromosome_capacity_kb,
                        lev,
                    ));
                    u.id = mid_id;
                }
            }
            if u.genome_substrate_stage.has_thymine_udg_repair() {
                tdna_cnt += 1;
                let parent_mid = u.id;
                let final_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = final_id;
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                let lev = (u.complexity_c / u.structure_s.max(2.0)).max(3.4);
                transitions.push((
                    final_id,
                    parent_mid,
                    u.species_id,
                    u.guild,
                    biome,
                    u.genome_substrate_stage,
                    u.chromosome_capacity_kb,
                    lev,
                ));
            }
        }
        for (nid, pid, sp, g, bm, stg, kb, lv) in transitions {
            self.phylogeny.record_rna_to_dna_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                stg,
                kb,
                lv,
            );
        }
        self.record_metrics();
        (rnr_cnt, tdna_cnt, self.rnr_virosphere.stable_genome_capacity_kb)
    }

    pub fn trigger_homochiral_sgc_pulse(&mut self) -> (usize, usize, f64) {
        self.controls.chiral_symmetry_breaking_mode = true;
        self.controls.magnetite_ciss_spin_bias = self.controls.magnetite_ciss_spin_bias.max(0.92);
        self.controls.minihelix_chiral_selectivity =
            self.controls.minihelix_chiral_selectivity.max(0.90);
        self.controls.communal_code_hgt_rate = self.controls.communal_code_hgt_rate.max(0.88);
        self.chiral_reactor = ChiralCodeReactor::new_for_controls(&self.controls);

        let mut homo_cnt = 0usize;
        let mut sgc_cnt = 0usize;
        let mut transitions = Vec::new();
        for u in &mut self.population {
            let old_id = u.id;
            u.seed_racemic_prebiotic_pool(&self.chiral_reactor);
            let mut crossed_homo = false;
            for _ in 0..8 {
                self.chiral_reactor.step_chiral_and_code_dynamics(&self.controls);
                let (st, _, _) =
                    u.step_homochiral_code_crystallization(&self.chiral_reactor, &mut self.rng);
                if st.is_homochiral_rna() && !crossed_homo {
                    crossed_homo = true;
                    homo_cnt += 1;
                    let mid_id = self.next_unit_id;
                    self.next_unit_id += 1;
                    let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                    let lev = (u.complexity_c / u.structure_s.max(2.0)).max(1.4);
                    transitions.push((
                        mid_id,
                        old_id,
                        u.species_id,
                        u.guild,
                        biome,
                        ChiralityCodeStage::HomochiralRnaSymmetryBroken,
                        u.genetic_code_mp,
                        lev,
                    ));
                    u.id = mid_id;
                }
            }
            if u.chirality_code_stage == ChiralityCodeStage::UniversalErrorMinimizedSgc {
                sgc_cnt += 1;
                let parent_mid = u.id;
                let final_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = final_id;
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                let lev = (u.complexity_c / u.structure_s.max(2.0)).max(2.8);
                transitions.push((
                    final_id,
                    parent_mid,
                    u.species_id,
                    u.guild,
                    biome,
                    ChiralityCodeStage::UniversalErrorMinimizedSgc,
                    u.genetic_code_mp,
                    lev,
                ));
            }
        }
        for (nid, pid, sp, g, bm, stg, mp, lv) in transitions {
            self.phylogeny.record_chirality_code_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                stg,
                mp,
                lv,
            );
        }
        self.record_metrics();
        (homo_cnt, sgc_cnt, self.chiral_reactor.sgc_error_minimization_mp)
    }

    pub fn trigger_hadean_abiogenesis_pulse(&mut self) -> (usize, usize, usize) {
        self.controls.hadean_abiogenesis_mode = true;
        self.controls.wet_dry_cycle_amplitude = self.controls.wet_dry_cycle_amplitude.max(0.84);
        self.controls.hydrothermal_delta_ph = self.controls.hydrothermal_delta_ph.max(3.35);
        self.controls.citrate_mg_chelation = self.controls.citrate_mg_chelation.max(0.86);
        self.controls.ribozyme_polymerase_fidelity =
            self.controls.ribozyme_polymerase_fidelity.max(0.92);
        self.hadean_bath = HadeanGeochemicalBath::new_for_controls(&self.controls);

        let mut raf_cnt = 0usize;
        let mut luca_cnt = 0usize;
        let mut transitions = Vec::new();
        for u in &mut self.population {
            let old_id = u.id;
            u.seed_hadean_prebiotic_pool(&self.hadean_bath);
            let mut crossed_raf = false;
            for _ in 0..10 {
                let ribo_f = if u.abiogenesis_stage.has_ribozyme_heredity() {
                    1.0
                } else {
                    0.35
                };
                self.hadean_bath.step_geochemical_cycle(
                    &self.controls,
                    u.raf_closure_score,
                    ribo_f,
                );
                let (st, _, _) = u.step_hadean_abiogenesis(&self.hadean_bath, &mut self.rng);
                if st.has_raf_closure() && !crossed_raf {
                    crossed_raf = true;
                    raf_cnt += 1;
                    let mid_id = self.next_unit_id;
                    self.next_unit_id += 1;
                    let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                    let lev = (u.complexity_c / u.structure_s.max(2.0)).max(1.2);
                    transitions.push((
                        mid_id,
                        old_id,
                        u.species_id,
                        u.guild,
                        biome,
                        AbiogenesisStage::AutocatalyticRafProtocell,
                        u.raf_closure_score,
                        lev,
                    ));
                    u.id = mid_id;
                }
            }
            if u.abiogenesis_stage == AbiogenesisStage::LucaCellularLineage {
                luca_cnt += 1;
                let parent_mid = u.id;
                let final_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = final_id;
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                let lev = (u.complexity_c / u.structure_s.max(2.0)).max(2.4);
                transitions.push((
                    final_id,
                    parent_mid,
                    u.species_id,
                    u.guild,
                    biome,
                    AbiogenesisStage::LucaCellularLineage,
                    u.raf_closure_score,
                    lev,
                ));
            }
        }
        for (nid, pid, sp, g, bm, stg, raf_s, lv) in transitions {
            self.phylogeny.record_abiogenesis_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                stg,
                raf_s,
                lv,
            );
        }
        self.hadean_bath.darwinian_threshold_crossed = luca_cnt > 0;
        self.record_metrics();
        (raf_cnt, luca_cnt, self.hadean_bath.eigen_max_length_nt)
    }

    pub fn trigger_planetary_awakening_pulse(&mut self) -> (usize, usize, f64) {
        self.controls.biogeochemical_feedback_gain =
            self.controls.biogeochemical_feedback_gain.max(0.75);
        self.controls.planetary_stewardship_coupling =
            self.controls.planetary_stewardship_coupling.max(0.88);
        let (_, fm, _) = self.trigger_machine_synthesis_pulse();

        let pop_f = (self.population.len().max(1)) as f64;
        if self.controls.rna_to_dna_virosphere_mode && !self.population.is_empty() {
            self.rnr_virosphere
                .step_rnr_and_virosphere_cycle(&self.controls);
            for u in &mut self.population {
                let (_, deam_bad, lytic_bad) = u
                    .step_rnr_dna_and_virosphere_coevolution(&self.rnr_virosphere, &mut self.rng);
                if deam_bad {
                    self.total_cytosine_deamination_meltdowns += 1;
                }
                if lytic_bad {
                    self.total_lytic_viral_collapses += 1;
                }
            }
        }
        if self.controls.chiral_symmetry_breaking_mode && !self.population.is_empty() {
            self.chiral_reactor
                .step_chiral_and_code_dynamics(&self.controls);
            for u in &mut self.population {
                let (_, cross_inh, babel_bad) =
                    u.step_homochiral_code_crystallization(&self.chiral_reactor, &mut self.rng);
                if cross_inh {
                    self.total_racemic_cross_inhibitions += 1;
                }
                if babel_bad {
                    self.total_mistranslation_babel_collapses += 1;
                }
            }
        }
        if self.controls.hadean_abiogenesis_mode && !self.population.is_empty() {
            let mean_raf = self
                .population
                .iter()
                .map(|u| u.raf_closure_score)
                .sum::<f64>()
                / pop_f;
            let ribo_frac = self
                .population
                .iter()
                .filter(|u| u.abiogenesis_stage.has_ribozyme_heredity())
                .count() as f64
                / pop_f;
            self.hadean_bath
                .step_geochemical_cycle(&self.controls, mean_raf, ribo_frac);
            for u in &mut self.population {
                let (_, eigen_bad, mg_bad) =
                    u.step_hadean_abiogenesis(&self.hadean_bath, &mut self.rng);
                if eigen_bad {
                    self.total_eigen_error_catastrophes += 1;
                }
                if mg_bad {
                    self.total_mg_vesicle_collapses += 1;
                }
            }
            self.hadean_bath.darwinian_threshold_crossed = self
                .population
                .iter()
                .any(|u| u.abiogenesis_stage.crossed_darwinian_threshold());
        }
        if self.controls.biogeochemical_feedback_gain > 0.001 && !self.population.is_empty() {
            let (mut s_o2, mut s_rec, mut s_weath, mut s_tech, mut s_stew) =
                (0.0, 0.0, 0.0, 0.0, 0.0);
            for u in &mut self.population {
                let (ao, sr, sw, tb, vs) =
                    u.compute_planetary_guild_fluxes(&self.controls, &self.biosphere);
                s_o2 += ao;
                s_rec += sr;
                s_weath += sw;
                s_tech += tb;
                s_stew += vs;
            }
            let (_, burnout) = self.biosphere.step_biogeochemical_feedbacks(
                s_o2 / pop_f,
                s_rec / pop_f,
                s_weath / pop_f,
                s_tech / pop_f,
                s_stew / pop_f,
                &self.controls,
            );
            self.controls.ambient_oxygen_level = self.biosphere.atmospheric_o2;
            if burnout {
                for u in &mut self.population {
                    if u.apply_planetary_boundary_stress(&self.biosphere) {
                        self.total_planetary_burnout_collapses += 1;
                    }
                }
            }
        } else {
            self.biosphere.atmospheric_o2 = self.controls.ambient_oxygen_level;
            self.biosphere.planetary_exergy_budget = self.controls.bath_exergy_flux.max(0.25);
        }
        let (mut sum_o2, mut sum_rec, mut sum_weath, mut sum_tech, mut sum_stew) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for u in &mut self.population {
            let (ao, sr, sw, tb, vs) =
                u.compute_planetary_guild_fluxes(&self.controls, &self.biosphere);
            sum_o2 += ao;
            sum_rec += sr;
            sum_weath += sw;
            sum_tech += tb;
            sum_stew += vs;
        }
        for _ in 0..8 {
            self.biosphere.step_biogeochemical_feedbacks(
                (sum_o2 / pop_f).max(0.85),
                (sum_rec / pop_f).max(0.80),
                (sum_weath / pop_f).max(0.55),
                sum_tech / pop_f,
                (sum_stew / pop_f).max(0.95),
                &self.controls,
            );
        }
        self.controls.ambient_oxygen_level = self.biosphere.atmospheric_o2;

        let mut awakened = 0usize;
        let mut transitions = Vec::new();
        let pop_n = self.population.len().max(1);
        for u in &mut self.population {
            if u.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis
                && self.biosphere.state == PlanetaryState::HomeostaticAwakening
            {
                let old_id = u.id;
                u.planetary_epoch = PlanetaryState::HomeostaticAwakening;
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le * self.biosphere.planetary_exergy_budget,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    u.eusocial_state,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    u.agriculture_state,
                    u.cultivar_polyploidy,
                    u.gongylidia_yield,
                    u.cultural_state,
                    u.cultural_ratchet_tier,
                    u.cultural_repertoire_bits(),
                    u.exosomatic_state,
                    u.exosomatic_archive_bits(),
                    u.machine_synthesis_state,
                    u.verified_theorem_bits(),
                    u.distribution_tail_diversity,
                    self.biosphere.state,
                    self.biosphere.chemical_disequilibrium_bits,
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                awakened += 1;
            }
        }
        for (
            nid,
            pid,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            agstate,
            ploidy,
            gongyl,
            cult_st,
            rtier,
            cbits,
            exo_st,
            abits,
            ms_st,
            tbits,
            tails,
            pl_st,
            diseq,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_planetary_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                agstate,
                ploidy,
                gongyl,
                cult_st,
                rtier,
                cbits,
                exo_st,
                abits,
                ms_st,
                tbits,
                tails,
                pl_st,
                diseq,
                egt_l,
                org_rem,
                lv,
            );
        }
        self.record_metrics();
        (fm, awakened, self.biosphere.planetary_exergy_budget)
    }

    pub fn trigger_machine_synthesis_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let fid = self.controls.symbolic_teaching_fidelity.max(0.92);
        let conn = self.controls.inter_deme_connectivity.max(0.25);
        let rigor = self.controls.peer_verification_rigor.max(0.85);
        let ver_s = self.controls.formal_verifier_strictness.max(0.88);
        let emp_g = self.controls.empirical_grounding_floor.max(0.25);
        let gen_r = self.controls.synthetic_generation_rate.max(0.65);

        let mut symbio_formed = 0usize;
        let mut formal_formed = 0usize;
        let mut proofs_checked = 0usize;

        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            EusocialState,
            u16,
            f64,
            AgricultureState,
            u8,
            f64,
            CulturalState,
            u8,
            f64,
            ExosomaticState,
            f64,
            MachineSynthesisState,
            f64,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild == TrophicGuild::PredatorParasite && !u.has_apicoplast_flag {
                continue;
            }
            while u.cultural_state != CulturalState::SymbolicCumulativeCulture {
                let (ch, _) = u.step_cumulative_culture_transition(fid, conn, &mut rng);
                if !ch {
                    break;
                }
            }
            while u.exosomatic_state != ExosomaticState::PeerVerifiedInstitution {
                let (ch, _) = u.step_exosomatic_archive_transition(rigor, &mut rng);
                if !ch {
                    break;
                }
            }
            let old_id = u.id;
            let (changed, new_ms) = u.step_machine_synthesis_transition(ver_s, emp_g, &mut rng);
            let (prf, collapsed) =
                u.step_recursive_synthesis_and_collapse_safeguards(gen_r, ver_s, emp_g);
            proofs_checked += prf as usize;
            if collapsed {
                self.total_recursive_model_collapses += 1;
            }
            if changed && new_ms.has_computational_symbiont() {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    u.eusocial_state,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    u.agriculture_state,
                    u.cultivar_polyploidy,
                    u.gongylidia_yield,
                    u.cultural_state,
                    u.cultural_ratchet_tier,
                    u.cultural_repertoire_bits(),
                    u.exosomatic_state,
                    u.exosomatic_archive_bits(),
                    new_ms,
                    u.verified_theorem_bits(),
                    u.distribution_tail_diversity,
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                symbio_formed += 1;
                if new_ms.is_formally_verified() {
                    formal_formed += 1;
                }
            }
        }

        self.total_formal_proofs_verified += proofs_checked;
        self.rng = rng;
        for (
            nid,
            pid,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            agstate,
            ploidy,
            gongyl,
            cult_st,
            rtier,
            cbits,
            exo_st,
            abits,
            ms_st,
            tbits,
            tails,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_machine_synthesis_transition(
                nid,
                pid,
                None,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                agstate,
                ploidy,
                gongyl,
                cult_st,
                rtier,
                cbits,
                exo_st,
                abits,
                ms_st,
                tbits,
                tails,
                egt_l,
                org_rem,
                lv,
            );
        }
        if symbio_formed > 0 || formal_formed > 0 || proofs_checked > 0 {
            self.record_metrics();
        }
        (symbio_formed, formal_formed, proofs_checked)
    }

    pub fn trigger_exosomatic_archive_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let fid = self.controls.symbolic_teaching_fidelity.max(0.92);
        let conn = self.controls.inter_deme_connectivity.max(0.25);
        let rigor = self.controls.peer_verification_rigor.max(0.82);
        let pub_p = self.controls.unverified_publish_pressure;
        let mut archives_formed = 0usize;
        let mut institutions_formed = 0usize;
        let mut audits_executed = 0usize;

        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            EusocialState,
            u16,
            f64,
            AgricultureState,
            u8,
            f64,
            CulturalState,
            u8,
            f64,
            ExosomaticState,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild == TrophicGuild::PredatorParasite && !u.has_apicoplast_flag {
                continue;
            }
            while u.cultural_state != CulturalState::SymbolicCumulativeCulture {
                let (ch, _) = u.step_cumulative_culture_transition(fid, conn, &mut rng);
                if !ch {
                    break;
                }
            }
            let old_id = u.id;
            let (changed, new_exo) = u.step_exosomatic_archive_transition(rigor, &mut rng);
            let (aud, corrupt) = u.step_archive_verification_and_fdr_dynamics(pub_p, rigor);
            audits_executed += aud as usize;
            if corrupt {
                self.total_archive_corruption_collapses += 1;
            }
            if changed && new_exo.is_archival() {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    u.eusocial_state,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    u.agriculture_state,
                    u.cultivar_polyploidy,
                    u.gongylidia_yield,
                    u.cultural_state,
                    u.cultural_ratchet_tier,
                    u.cultural_repertoire_bits(),
                    new_exo,
                    u.exosomatic_archive_bits(),
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                if new_exo == ExosomaticState::WrittenArchive {
                    archives_formed += 1;
                } else if new_exo == ExosomaticState::PeerVerifiedInstitution {
                    institutions_formed += 1;
                }
            }
        }

        self.total_peer_retraction_audits += audits_executed;
        self.rng = rng;
        for (
            nid,
            pid,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            agstate,
            ploidy,
            gongyl,
            cult_st,
            rtier,
            cbits,
            exo_st,
            abits,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_exosomatic_transition(
                nid,
                pid,
                None,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                agstate,
                ploidy,
                gongyl,
                cult_st,
                rtier,
                cbits,
                exo_st,
                abits,
                egt_l,
                org_rem,
                lv,
            );
        }
        if archives_formed > 0 || institutions_formed > 0 || audits_executed > 0 {
            self.record_metrics();
        }
        (archives_formed, institutions_formed, audits_executed)
    }

    pub fn trigger_cultural_ratchet_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let fid = self.controls.symbolic_teaching_fidelity.max(0.90);
        let conn = self.controls.inter_deme_connectivity;
        let mut ratchet_formed = 0usize;
        let mut symbolic_formed = 0usize;
        let mut cross_deme_recombs = 0usize;

        let deme_teachers: Vec<ModeledUnit> = self
            .population
            .iter()
            .filter(|u| u.cultural_state.is_cultural())
            .cloned()
            .collect();

        let mut transitions: Vec<(
            u64,
            u64,
            Option<u64>,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            EusocialState,
            u16,
            f64,
            AgricultureState,
            u8,
            f64,
            CulturalState,
            u8,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild == TrophicGuild::PredatorParasite && !u.has_apicoplast_flag {
                continue;
            }
            let old_id = u.id;
            let (changed, new_cult) = u.step_cumulative_culture_transition(fid, conn.max(0.15), &mut rng);
            let mut teacher_donor_id: Option<u64> = None;

            if u.cultural_state.is_cumulative_ratchet() {
                if let Some(teacher) = deme_teachers.iter().find(|t| t.deme_id != u.deme_id) {
                    let added = u.fuse_cross_deme_cultural_recombination(teacher, conn, fid);
                    if added > 0 {
                        cross_deme_recombs += 1;
                        teacher_donor_id = Some(teacher.id);
                    }
                }
            }

            if changed || teacher_donor_id.is_some() {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    teacher_donor_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    u.eusocial_state,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    u.agriculture_state,
                    u.cultivar_polyploidy,
                    u.gongylidia_yield,
                    u.cultural_state,
                    u.cultural_ratchet_tier,
                    u.cultural_repertoire_bits(),
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                if changed {
                    if new_cult == CulturalState::ImitativeRatchet {
                        ratchet_formed += 1;
                    } else if new_cult == CulturalState::SymbolicCumulativeCulture {
                        symbolic_formed += 1;
                    }
                }
            }
        }

        self.total_cross_deme_recombinations += cross_deme_recombs;
        self.rng = rng;
        for (
            nid,
            pid,
            tdonor,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            agstate,
            ploidy,
            gongyl,
            cult_st,
            rtier,
            cbits,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_cultural_transition(
                nid,
                pid,
                tdonor,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                agstate,
                ploidy,
                gongyl,
                cult_st,
                rtier,
                cbits,
                egt_l,
                org_rem,
                lv,
            );
        }
        if ratchet_formed > 0 || symbolic_formed > 0 || cross_deme_recombs > 0 {
            self.record_metrics();
        }
        (ratchet_formed, symbolic_formed, cross_deme_recombs)
    }

    pub fn trigger_agriculture_domestication_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let monogamy = self.controls.queen_monogamy_index.max(0.85);
        let bottleneck = self.controls.single_cell_bottleneck;
        let actino = self.controls.defensive_symbiont_potency.max(0.75);
        let escov_p = self.controls.agro_pathogen_pressure;
        let mut ag_formed = 0usize;
        let mut obligate_crops = 0usize;
        let mut shields_fired = 0usize;
        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            EusocialState,
            u16,
            f64,
            AgricultureState,
            u8,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild != TrophicGuild::Holobiont {
                continue;
            }
            while u.multicell_state != MulticellularState::GermSomaMulticellular {
                let (ch, _) = u.step_clonal_multicellularity(bottleneck, &mut rng);
                if !ch {
                    break;
                }
            }
            while u.eusocial_state != EusocialState::ObligateSuperorganism {
                let (ch, _) = u.step_eusocial_superorganism(monogamy, &mut rng);
                if !ch {
                    break;
                }
            }
            let old_id = u.id;
            let (changed, new_ag) = u.step_agricultural_domestication(actino, &mut rng);
            let (sh, blight) = u.step_agro_pathogen_and_defensive_microbiome(escov_p, actino);
            shields_fired += sh as usize;
            if blight {
                self.total_garden_blight_collapses += 1;
            }
            if changed && new_ag.is_farming() {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    u.eusocial_state,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    new_ag,
                    u.cultivar_polyploidy,
                    u.gongylidia_yield,
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                ag_formed += 1;
                if new_ag.is_obligate_monoculture() {
                    obligate_crops += 1;
                }
            }
        }

        self.total_antimicrobial_shields += shields_fired;
        self.rng = rng;
        for (
            nid,
            pid,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            agstate,
            ploidy,
            gongyl,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_agriculture_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                agstate,
                ploidy,
                gongyl,
                egt_l,
                org_rem,
                lv,
            );
        }
        if ag_formed > 0 || obligate_crops > 0 || shields_fired > 0 {
            self.record_metrics();
        }
        (ag_formed, obligate_crops, shields_fired)
    }

    pub fn trigger_superorganism_transition_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let monogamy = self.controls.queen_monogamy_index;
        let police_rate = self.controls.worker_policing_rate.max(0.65);
        let bottleneck = self.controls.single_cell_bottleneck;
        let mut eusoc_formed = 0usize;
        let mut superorg_formed = 0usize;
        let mut policed_eggs = 0usize;
        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            EusocialState,
            u16,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild != TrophicGuild::Holobiont {
                continue;
            }
            if u.multicell_state != MulticellularState::GermSomaMulticellular {
                u.step_clonal_multicellularity(bottleneck, &mut rng);
                if u.multicell_state != MulticellularState::GermSomaMulticellular {
                    u.step_clonal_multicellularity(bottleneck, &mut rng);
                }
            }
            let old_id = u.id;
            let (changed, new_eu) = u.step_eusocial_superorganism(monogamy, &mut rng);
            let (pol, anarchy) = u.step_colony_kin_policing(monogamy, police_rate);
            policed_eggs += pol as usize;
            if anarchy {
                self.total_colony_anarchy_collapses += 1;
            }
            if changed && new_eu.is_eusocial() {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    u.multicell_state,
                    u.cell_count,
                    u.soma_fraction,
                    new_eu,
                    u.colony_size,
                    u.total_sterile_caste_fraction(),
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                match new_eu {
                    EusocialState::FacultativeEusocial => eusoc_formed += 1,
                    EusocialState::ObligateSuperorganism => superorg_formed += 1,
                    _ => {}
                }
            }
        }

        self.total_policed_eggs += policed_eggs;
        self.rng = rng;
        for (
            nid,
            pid,
            sp,
            g,
            bm,
            tier,
            mc,
            mp,
            mtdna,
            mstate,
            cells,
            soma,
            eustate,
            col_sz,
            caste_f,
            egt_l,
            org_rem,
            lv,
        ) in transitions
        {
            self.phylogeny.record_superorganism_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                eustate,
                col_sz,
                caste_f,
                egt_l,
                org_rem,
                lv,
            );
        }
        if eusoc_formed > 0 || superorg_formed > 0 || policed_eggs > 0 {
            self.record_metrics();
        }
        (eusoc_formed, superorg_formed, policed_eggs)
    }

    pub fn trigger_multicellular_transition_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let bottleneck = self.controls.single_cell_bottleneck;
        let apop_rate = self.controls.apoptosis_policing_rate.max(0.65);
        let mut clonal_formed = 0usize;
        let mut gs_formed = 0usize;
        let mut apop_culls = 0usize;
        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            MulticellularState,
            u16,
            f64,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild != TrophicGuild::Holobiont {
                continue;
            }
            let old_id = u.id;
            let (changed, new_mc) = u.step_clonal_multicellularity(bottleneck, &mut rng);
            let (culls, collapsed) = u.step_multicellular_homeostasis(apop_rate, bottleneck);
            apop_culls += culls as usize;
            if collapsed {
                self.total_neoplastic_collapses += 1;
            }
            if changed {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let new_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_id;
                let lev = u.complexity_c / u.structure_s.max(2.0);
                transitions.push((
                    new_id,
                    old_id,
                    u.species_id,
                    u.guild,
                    biome,
                    u.endo_tier,
                    u.mro_class,
                    u.mro_pathway_count(),
                    u.mro_mtdna_loci,
                    new_mc,
                    u.cell_count,
                    u.soma_fraction,
                    u.egt_transferred_loci,
                    u.endosymbiont_sequence.len(),
                    lev,
                ));
                match new_mc {
                    MulticellularState::UndifferentiatedClonal
                    | MulticellularState::AggregativeChimerism => clonal_formed += 1,
                    MulticellularState::GermSomaMulticellular => gs_formed += 1,
                    MulticellularState::Unicellular => {}
                }
            }
        }

        self.total_apoptotic_culls += apop_culls;
        self.rng = rng;
        for (nid, pid, sp, g, bm, tier, mc, mp, mtdna, mstate, cells, soma, egt_l, org_rem, lv) in
            transitions
        {
            self.phylogeny.record_multicellular_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                mstate,
                cells,
                soma,
                egt_l,
                org_rem,
                lv,
            );
        }
        if clonal_formed > 0 || gs_formed > 0 || apop_culls > 0 {
            self.record_metrics();
        }
        (clonal_formed, gs_formed, apop_culls)
    }

    pub fn trigger_mro_reduction_pulse(&mut self) -> (usize, usize, usize) {
        let pop_n = self.population.len().max(1);
        let suf_rescue = self.controls.cytosolic_suf_hgt_rescue;
        let mut hydro_formed = 0usize;
        let mut mitos_formed = 0usize;
        let mut mro_lost = 0usize;
        let mut transitions: Vec<(
            u64,
            u64,
            usize,
            TrophicGuild,
            BiomeZone,
            u8,
            MroClass,
            u8,
            usize,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.mro_class == MroClass::None || u.mro_class == MroClass::Amiteukaryote {
                continue;
            }
            let old_id = u.id;
            let (changed, new_cls, _) = u.step_anaerobic_mro_reduction(suf_rescue, false);
            if changed {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                if matches!(
                    new_cls,
                    MroClass::Hydrogenosome | MroClass::Mitosome | MroClass::Amiteukaryote
                ) {
                    let new_id = self.next_unit_id;
                    self.next_unit_id += 1;
                    u.id = new_id;
                    let lev = u.complexity_c / u.structure_s.max(2.0);
                    transitions.push((
                        new_id,
                        old_id,
                        u.species_id,
                        u.guild,
                        biome,
                        u.endo_tier,
                        new_cls,
                        u.mro_pathway_count(),
                        u.mro_mtdna_loci,
                        u.egt_transferred_loci,
                        u.endosymbiont_sequence.len(),
                        lev,
                    ));
                }
                match new_cls {
                    MroClass::Hydrogenosome => hydro_formed += 1,
                    MroClass::Mitosome => mitos_formed += 1,
                    MroClass::Amiteukaryote => mro_lost += 1,
                    _ => {}
                }
            }
        }

        for (nid, pid, sp, g, bm, tier, mc, mp, mtdna, egt_l, org_rem, lv) in transitions {
            self.phylogeny.record_mro_transition(
                nid,
                pid,
                self.step_count,
                sp,
                g,
                bm,
                tier,
                mc,
                mp,
                mtdna,
                egt_l,
                org_rem,
                lv,
            );
        }
        if hydro_formed > 0 || mitos_formed > 0 || mro_lost > 0 {
            self.record_metrics();
        }
        (hydro_formed, mitos_formed, mro_lost)
    }

    pub fn trigger_apicoplast_reduction_pulse(&mut self) -> (usize, usize, usize) {
        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let rescue = self.controls.host_metabolite_rescue;
        let mut pathways_shed = 0usize;
        let mut crypto_losses = 0usize;

        for u in &mut self.population {
            if u.has_apicoplast() {
                let (shed, lost_all, _) = u.streamline_apicoplast_pathway(rescue, false);
                if shed {
                    pathways_shed += 1;
                }
                if lost_all {
                    crypto_losses += 1;
                }
                if shed || lost_all {
                    let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                    let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                    let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                    let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                    u.recompute_phenome(
                        &self.controls,
                        le,
                        lb,
                        ln,
                        0.28,
                        biome,
                        self.last_adami_rho,
                        pop_n,
                        self.env_phase,
                    );
                }
            }
        }

        let mut converted = 0usize;
        if let Some(holo_idx) = self
            .population
            .iter()
            .enumerate()
            .filter(|(_, u)| u.guild == TrophicGuild::Holobiont)
            .max_by_key(|(_, u)| (u.endo_tier, (u.fitness * 1000.0) as i64))
            .map(|(i, _)| i)
        {
            let old_holo_id = self.population[holo_idx].id;
            if self.population[holo_idx].convert_holobiont_to_apicoplast_parasite(&mut rng) {
                let new_apico_id = self.next_unit_id;
                self.next_unit_id += 1;
                let gx = self.population[holo_idx].grid_x;
                let gy = self.population[holo_idx].grid_y;
                let biome = self.grid.biome_at(gx, gy);
                let le = self.grid.local_exergy(gx, gy);
                let lb = self.grid.local_byproduct(gx, gy);
                let ln = self.grid.local_niche(gx, gy);
                self.population[holo_idx].id = new_apico_id;
                self.population[holo_idx].recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.28,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                let u_ref = &self.population[holo_idx];
                let lev = u_ref.complexity_c / u_ref.structure_s.max(2.0);
                self.phylogeny.record_apicoplast_transition(
                    new_apico_id,
                    old_holo_id,
                    self.step_count,
                    u_ref.species_id,
                    biome,
                    u_ref.endo_tier,
                    u_ref.egt_transferred_loci,
                    u_ref.endosymbiont_sequence.len(),
                    u_ref.apicoplast_pathway_count(),
                    lev,
                );
                converted = 1;
            }
        }

        self.rng = rng;
        if converted > 0 || pathways_shed > 0 || crypto_losses > 0 {
            self.record_metrics();
        }
        (converted, pathways_shed, crypto_losses)
    }

    pub fn trigger_klepto_steal_pulse(&mut self) -> (usize, usize, usize) {
        let prey_candidates: Vec<ModeledUnit> = self
            .population
            .iter()
            .filter(|u| u.guild != TrophicGuild::PredatorParasite)
            .cloned()
            .collect();
        if prey_candidates.is_empty() {
            return (0, 0, 0);
        }

        let mut rng = self.rng.clone();
        let pop_n = self.population.len().max(1);
        let mut stolen = 0usize;
        let mut hgt_gained = 0usize;
        let mut fixed = 0usize;
        let mut tipping_transitions: Vec<(
            u64,
            u64,
            Option<u64>,
            usize,
            BiomeZone,
            u8,
            f64,
            usize,
            usize,
            usize,
            f64,
        )> = Vec::new();

        for u in &mut self.population {
            if u.guild != TrophicGuild::PredatorParasite
                || u.has_apicoplast()
                || u.apicoplast_eliminated
            {
                continue;
            }
            let prey = &prey_candidates[rng.gen_index(prey_candidates.len())];
            let old_pred_id = u.id;
            let (did_steal, did_hgt, crossed) =
                u.steal_kleptoplast_from_prey(prey, true, true, true, &mut rng);
            if did_steal {
                stolen += 1;
            }
            if did_hgt {
                hgt_gained += 1;
            }
            let le = self.grid.local_exergy(u.grid_x, u.grid_y);
            let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
            let ln = self.grid.local_niche(u.grid_x, u.grid_y);
            let biome = self.grid.biome_at(u.grid_x, u.grid_y);
            u.recompute_phenome(
                &self.controls,
                le,
                lb,
                ln,
                0.30,
                biome,
                self.last_adami_rho,
                pop_n,
                self.env_phase,
            );
            if crossed {
                let new_holo_id = self.next_unit_id;
                self.next_unit_id += 1;
                u.id = new_holo_id;
                let gamma = (u.endosymbiont_sequence.len() as f64)
                    / ((u.sequence.len() + u.endosymbiont_sequence.len()).max(1) as f64);
                let lev = u.complexity_c / u.structure_s.max(2.0);
                tipping_transitions.push((
                    new_holo_id,
                    old_pred_id,
                    u.endosymbiont_origin_id,
                    u.species_id,
                    biome,
                    u.endo_tier,
                    gamma,
                    u.klepto_hgt_loci,
                    u.endosymbiont_sequence.len(),
                    u.nucleomorph_sequence.len(),
                    lev,
                ));
                fixed += 1;
            }
        }
        self.rng = rng;
        for (hid, pid, did, sp, bm, tier, gm, hgt_l, org_rem, nm_l, lv) in tipping_transitions {
            self.phylogeny.record_klepto_tipping_fixation(
                hid,
                pid,
                did,
                self.step_count,
                sp,
                bm,
                tier,
                gm,
                hgt_l,
                org_rem,
                nm_l,
                lv,
            );
        }
        if stolen > 0 || hgt_gained > 0 || fixed > 0 {
            self.record_metrics();
        }
        (stolen, hgt_gained, fixed)
    }

    pub fn trigger_endosymbiotic_pulse(&mut self) -> bool {
        let prod_idx = self
            .population
            .iter()
            .enumerate()
            .filter(|(_, u)| u.guild == TrophicGuild::Producer)
            .max_by(|(_, a), (_, b)| a.fitness.total_cmp(&b.fitness))
            .map(|(i, _)| i);
        let syn_idx = self
            .population
            .iter()
            .enumerate()
            .filter(|(_, u)| u.guild == TrophicGuild::Syntroph)
            .max_by(|(_, a), (_, b)| a.fitness.total_cmp(&b.fitness))
            .map(|(i, _)| i);

        let (Some(pi), Some(si)) = (prod_idx, syn_idx) else {
            return false;
        };
        if pi == si {
            return false;
        }
        let host = self.population[pi].clone();
        let sym = self.population[si].clone();
        let hid = self.alloc_id();
        let biome = self.grid.biome_at(host.grid_x, host.grid_y);
        let le = self.grid.local_exergy(host.grid_x, host.grid_y);
        let lb = self.grid.local_byproduct(host.grid_x, host.grid_y);
        let ln = self.grid.local_niche(host.grid_x, host.grid_y);
        let pop_n = self.population.len();
        let (holo, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
            &host,
            &sym,
            hid,
            biome,
            le,
            lb,
            ln,
            &self.controls,
            self.last_adami_rho,
            pop_n,
            self.env_phase,
        );
        let lev = holo.complexity_c / holo.structure_s.max(2.0);
        self.phylogeny.record_endosymbiotic_fusion(
            holo.id,
            host.id,
            sym.id,
            self.step_count,
            holo.species_id,
            biome,
            gamma,
            holo.endosymbiont_sequence.len(),
            lev,
        );
        let hi = pi.max(si);
        let lo = pi.min(si);
        self.population.remove(hi);
        self.population.remove(lo);
        self.population.push(holo);
        self.record_metrics();
        true
    }

    pub fn trigger_serial_endosymbiosis_pulse(&mut self) -> Option<u8> {
        let engulfed_idx = self
            .population
            .iter()
            .enumerate()
            .filter(|(_, u)| u.guild == TrophicGuild::Holobiont && u.endo_tier <= 2)
            .max_by_key(|(_, u)| (u.endo_tier, (u.fitness * 1000.0) as i64))
            .map(|(i, _)| i)?;

        let host_idx = self
            .population
            .iter()
            .enumerate()
            .filter(|(i, u)| {
                *i != engulfed_idx
                    && matches!(
                        u.guild,
                        TrophicGuild::Holobiont
                            | TrophicGuild::PredatorParasite
                            | TrophicGuild::Producer
                    )
            })
            .max_by_key(|(_, u)| {
                let holo_pref = u8::from(u.guild == TrophicGuild::Holobiont);
                (holo_pref, (u.fitness * 1000.0) as i64)
            })
            .map(|(i, _)| i)?;

        let host = self.population[host_idx].clone();
        let engulfed = self.population[engulfed_idx].clone();
        let cid = self.alloc_id();
        let biome = self.grid.biome_at(host.grid_x, host.grid_y);
        let le = self.grid.local_exergy(host.grid_x, host.grid_y);
        let lb = self.grid.local_byproduct(host.grid_x, host.grid_y);
        let ln = self.grid.local_niche(host.grid_x, host.grid_y);
        let pop_n = self.population.len();

        let (chimera, gamma) = ModeledUnit::fuse_serial_holobiont(
            &host,
            &engulfed,
            cid,
            biome,
            le,
            lb,
            ln,
            &self.controls,
            self.last_adami_rho,
            pop_n,
            self.env_phase,
        );
        let tier = chimera.endo_tier;
        let lev = chimera.complexity_c / chimera.structure_s.max(2.0);
        self.phylogeny.record_serial_endosymbiotic_fusion(
            chimera.id,
            host.id,
            engulfed.id,
            self.step_count,
            chimera.species_id,
            biome,
            tier,
            gamma,
            chimera.endosymbiont_sequence.len(),
            chimera.nucleomorph_sequence.len(),
            chimera.egt_transferred_loci,
            lev,
        );

        let hi = host_idx.max(engulfed_idx);
        let lo = host_idx.min(engulfed_idx);
        self.population.remove(hi);
        self.population.remove(lo);
        self.population.push(chimera);
        self.record_metrics();
        Some(tier)
    }

    pub fn trigger_egt_transfer_pulse(&mut self) -> usize {
        let mut rng = self.rng.clone();
        let mut moved = 0usize;
        let pop_n = self.population.len().max(1);
        for u in &mut self.population {
            if u.guild == TrophicGuild::Holobiont
                && u.migrate_endosymbiont_locus_to_host(false, &mut rng)
            {
                let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                u.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.30,
                    biome,
                    self.last_adami_rho,
                    pop_n,
                    self.env_phase,
                );
                moved += 1;
            }
        }
        self.rng = rng;
        if moved > 0 {
            self.record_metrics();
        }
        moved
    }

    pub fn trigger_nucleomorph_reduction_pulse(&mut self, allow_final_sweep: bool) -> (usize, usize) {
        let mut rng = self.rng.clone();
        let mut loci_moved = 0usize;
        let mut nm_swept = 0usize;
        let pop_n = self.population.len().max(1);
        for u in &mut self.population {
            if u.endo_tier >= 2 && !u.nucleomorph_sequence.is_empty() {
                let (m, swept) = u.migrate_nucleomorph_locus_to_host(allow_final_sweep, &mut rng);
                if m > 0 {
                    loci_moved += m;
                    if swept {
                        nm_swept += 1;
                    }
                    let le = self.grid.local_exergy(u.grid_x, u.grid_y);
                    let lb = self.grid.local_byproduct(u.grid_x, u.grid_y);
                    let ln = self.grid.local_niche(u.grid_x, u.grid_y);
                    let biome = self.grid.biome_at(u.grid_x, u.grid_y);
                    u.recompute_phenome(
                        &self.controls,
                        le,
                        lb,
                        ln,
                        0.32,
                        biome,
                        self.last_adami_rho,
                        pop_n,
                        self.env_phase,
                    );
                }
            }
        }
        self.rng = rng;
        if loci_moved > 0 {
            self.record_metrics();
        }
        (loci_moved, nm_swept)
    }

    pub fn apply_regime_preset(&mut self, regime: RegimeId, preserve_champion: bool) {
        let coming_from_meltdown = matches!(
            self.target_regime,
            RegimeId::R00Static | RegimeId::R05Decaying | RegimeId::R20ProtoDarwinian
        );
        self.target_regime = regime;
        self.controls = regime.preset_controls();
        if self.coupled_biosphere_mode && !preserve_champion {
            self.biosphere = PlanetaryBiosphere::new_for_controls(&self.controls);
        }
        self.stage_step = 0;
        self.telemetry_window.reset();
        self.last_observed_fidelity = None;
        self.seed_observable_estimators();

        let saved_niche = if preserve_champion && regime >= RegimeId::R20AutoadaptiveReplicative {
            Some((self.grid.niche.clone(), self.grid.byproduct.clone()))
        } else {
            None
        };
        self.grid = SpatialGrid::new(
            DEFAULT_GRID_DIM,
            DEFAULT_GRID_DIM,
            self.controls.bath_exergy_flux,
        );
        if let Some((niche_vec, byp_vec)) = saved_niche {
            if niche_vec.len() == self.grid.niche.len() {
                self.grid.niche = niche_vec;
                self.grid.byproduct = byp_vec;
            }
        }

        let champion = if preserve_champion && !coming_from_meltdown {
            self.population
                .iter()
                .filter(|u| u.guild == TrophicGuild::Producer)
                .max_by(|a, b| a.fitness.total_cmp(&b.fitness))
                .or_else(|| self.best_unit())
                .cloned()
        } else {
            None
        };

        let mut rng = self.rng.clone();
        if let Some(mut champ) = champion {
            champ.integrity = 1.0;
            champ.grid_x = 4;
            champ.grid_y = 4;
            if regime < RegimeId::R35Symbiogenetic && self.controls.endosymbiosis_rate < 0.04 {
                champ.endosymbiont_sequence.clear();
                champ.nucleomorph_sequence.clear();
                champ.endosymbiont_origin_id = None;
                champ.egt_transferred_loci = 0;
                champ.nucleomorph_transferred_loci = 0;
                champ.nucleomorph_eliminated = false;
                champ.endo_tier = 0;
                champ.organelle_membranes = 0;
                champ.kleptoplast_sequence.clear();
                champ.kleptoplast_integrity = 0.0;
                champ.kleptoplast_donor_id = None;
                champ.kleptoplast_donor_tier = 0;
                champ.klepto_hgt_loci = 0;
                champ.klepto_donor_mask = 0;
                champ.has_apicoplast_flag = false;
                champ.apicoplast_pathways_mask = 0;
                champ.apicoplast_eliminated = false;
                champ.mro_class = MroClass::None;
                champ.mro_pathways_mask = 0;
                champ.mro_mtdna_loci = 0;
            }
            if regime == RegimeId::R00Static {
                champ = ModeledUnit::new_founder(
                    self.alloc_id(),
                    self.organism,
                    regime,
                    4,
                    4,
                    BiomeZone::Alpha,
                    self.env_phase,
                    &mut rng,
                );
            } else {
                champ.wiring = ModeledUnit::build_wiring(self.organism, &champ.sequence, regime);
            }

            let mem_retain = if self.controls.soft_inheritance_coupling >= 0.10 {
                (self.controls.soft_inheritance_coupling * 0.55).clamp(0.0, 0.55)
            } else {
                0.0
            };
            for m in &mut champ.lifetime_memory {
                *m *= mem_retain;
            }
            champ.age = 1;

            let initial_pop = match regime {
                RegimeId::R00Static | RegimeId::R15Autoadaptive => 1,
                RegimeId::R05Decaying | RegimeId::R10AdaptiveDissipation => 3,
                RegimeId::R15Replicative => 4,
                _ => 16,
            };

            self.population.clear();
            let loc_e = self.grid.local_exergy(4, 4);
            let loc_b = self.grid.local_byproduct(4, 4);
            let loc_n = self.grid.local_niche(4, 4);
            champ.recompute_phenome(
                &self.controls,
                loc_e,
                loc_b,
                loc_n,
                0.0,
                BiomeZone::Alpha,
                1.0,
                initial_pop,
                self.env_phase,
            );
            let champ_id = champ.id;
            let champ_sp = champ.species_id;
            let champ_guild = champ.guild;
            self.population.push(champ.clone());

            let multi_trophic_stage = regime >= RegimeId::R30FullyDarwinian
                && (self.controls.trophic_coupling >= 0.05
                    || self.controls.cross_feeding_rate >= 0.05);

            for k in 1..initial_pop {
                let cid = self.alloc_id();
                let (gx, gy) = if regime >= RegimeId::R30FullyDarwinian && k % 2 == 1 {
                    (11 + (k % 3), 11 + ((k / 3) % 3))
                } else {
                    (3 + (k % 3), 3 + ((k / 3) % 3))
                };
                let biome = self.grid.biome_at(gx, gy);
                let le = self.grid.local_exergy(gx, gy);
                let lb = self.grid.local_byproduct(gx, gy);
                let ln = self.grid.local_niche(gx, gy);

                let unit_added = if multi_trophic_stage && (k == 4 || k == 9 || k == 13) {
                    let mut syn = ModeledUnit::new_founder_with_role(
                        cid,
                        self.organism,
                        regime,
                        TrophicGuild::Syntroph,
                        gx,
                        gy,
                        biome,
                        self.env_phase,
                        &mut rng,
                    );
                    syn.age = 1;
                    syn
                } else if multi_trophic_stage && (k == 7 || k == 14) {
                    let mut pred = ModeledUnit::new_founder_with_role(
                        cid,
                        self.organism,
                        regime,
                        TrophicGuild::PredatorParasite,
                        gx,
                        gy,
                        biome,
                        self.env_phase,
                        &mut rng,
                    );
                    pred.age = 1;
                    pred
                } else if regime >= RegimeId::R20ProtoDarwinian {
                    let (mut child, _, _, _) = champ.replicate_with_partner(
                        None,
                        cid,
                        gx,
                        gy,
                        biome,
                        le,
                        lb,
                        ln,
                        &self.controls,
                        0.92,
                        initial_pop,
                        self.env_phase,
                        &mut rng,
                    );
                    child.age = 1;
                    child
                } else {
                    let mut clone_u = champ.clone();
                    clone_u.id = cid;
                    clone_u.grid_x = gx;
                    clone_u.grid_y = gy;
                    clone_u.age = 1;
                    clone_u.recompute_phenome(
                        &self.controls,
                        le,
                        lb,
                        ln,
                        0.0,
                        biome,
                        1.0,
                        initial_pop,
                        self.env_phase,
                    );
                    clone_u
                };

                let div = normalized_sequence_divergence(&champ.sequence, &unit_added.sequence);
                let lev = unit_added.complexity_c / unit_added.structure_s.max(2.0);
                self.phylogeny.record_birth_full(
                    unit_added.id,
                    champ_id,
                    None,
                    self.step_count,
                    unit_added.species_id,
                    unit_added.guild,
                    biome,
                    champ_sp,
                    champ_guild,
                    div,
                    self.controls.speciation_threshold,
                    lev,
                    unit_added.egt_transferred_loci,
                    unit_added.endosymbiont_sequence.len(),
                    unit_added.endo_tier,
                    unit_added.nucleomorph_sequence.len(),
                    unit_added.nucleomorph_eliminated,
                );
                self.population.push(unit_added);
            }
        } else {
            self.populate_for_regime(regime, &mut rng);
            self.reset_phylogeny_from_population();
        }
        self.rng = rng;
    }

    fn populate_for_regime(&mut self, regime: RegimeId, rng: &mut Rng) {
        self.population.clear();
        let count = match regime {
            RegimeId::R00Static | RegimeId::R15Autoadaptive => 1,
            RegimeId::R05Decaying | RegimeId::R10AdaptiveDissipation => 3,
            RegimeId::R15Replicative => 4,
            _ => 16,
        };

        let founder_id = self.alloc_id();
        let mut founder = ModeledUnit::new_founder(
            founder_id,
            self.organism,
            regime,
            4,
            4,
            BiomeZone::Alpha,
            self.env_phase,
            rng,
        );
        let le0 = self.grid.local_exergy(4, 4);
        let lb0 = self.grid.local_byproduct(4, 4);
        let ln0 = self.grid.local_niche(4, 4);
        founder.recompute_phenome(
            &self.controls,
            le0,
            lb0,
            ln0,
            0.0,
            BiomeZone::Alpha,
            1.0,
            count,
            self.env_phase,
        );
        self.population.push(founder.clone());

        let multi_trophic_stage = regime >= RegimeId::R30FullyDarwinian
            && (self.controls.trophic_coupling >= 0.05
                || self.controls.cross_feeding_rate >= 0.05);

        for k in 1..count {
            let cid = self.alloc_id();
            let (gx, gy) = if regime >= RegimeId::R30FullyDarwinian && k % 2 == 1 {
                (11 + (k % 3), 11 + ((k / 3) % 3))
            } else {
                (3 + (k % 3), 3 + ((k / 3) % 3))
            };
            let biome = self.grid.biome_at(gx, gy);
            let le = self.grid.local_exergy(gx, gy);
            let lb = self.grid.local_byproduct(gx, gy);
            let ln = self.grid.local_niche(gx, gy);

            if multi_trophic_stage && (k == 4 || k == 9 || k == 13) {
                let syn = ModeledUnit::new_founder_with_role(
                    cid,
                    self.organism,
                    regime,
                    TrophicGuild::Syntroph,
                    gx,
                    gy,
                    biome,
                    self.env_phase,
                    rng,
                );
                self.population.push(syn);
            } else if multi_trophic_stage && (k == 7 || k == 14) {
                let pred = ModeledUnit::new_founder_with_role(
                    cid,
                    self.organism,
                    regime,
                    TrophicGuild::PredatorParasite,
                    gx,
                    gy,
                    biome,
                    self.env_phase,
                    rng,
                );
                self.population.push(pred);
            } else if regime >= RegimeId::R20ProtoDarwinian {
                let (child, _, _, _) = founder.replicate_with_partner(
                    None,
                    cid,
                    gx,
                    gy,
                    biome,
                    le,
                    lb,
                    ln,
                    &self.controls,
                    0.88,
                    count,
                    self.env_phase,
                    rng,
                );
                self.population.push(child);
            } else {
                let mut copy = founder.clone();
                copy.id = cid;
                copy.grid_x = gx;
                copy.grid_y = gy;
                copy.recompute_phenome(
                    &self.controls,
                    le,
                    lb,
                    ln,
                    0.0,
                    biome,
                    1.0,
                    count,
                    self.env_phase,
                );
                self.population.push(copy);
            }
        }
    }

    pub fn step(&mut self) {
        if self.mode == SimulationMode::GuidedCurriculum
            && self.stage_step >= self.steps_per_stage
        {
            let next_idx = (self.target_regime.index() + 1) % RegimeId::ALL.len();
            let preserve = next_idx > 2;
            self.apply_regime_preset(RegimeId::ALL[next_idx], preserve);
        }

        self.step_count += 1;
        self.stage_step += 1;
        self.phase_time += 0.16;
        self.env_phase += self.controls.env_drift_rate;

        let mut rng = self.rng.clone();

        self.grid.step_reaction_diffusion(
            self.controls.bath_exergy_flux,
            self.controls.exergy_diffusion,
            self.controls.niche_decay_rate,
            0.035,
        );

        let grid_cells = self.grid.width * self.grid.height;
        let mut host_counts = vec![0u16; grid_cells];
        let mut predator_counts = vec![0u16; grid_cells];
        let mut syntroph_counts = vec![0u16; grid_cells];
        for u in &self.population {
            let idx = self.grid.idx(u.grid_x, u.grid_y);
            match u.guild {
                TrophicGuild::Producer | TrophicGuild::Holobiont => {
                    host_counts[idx] = host_counts[idx].saturating_add(1);
                }
                TrophicGuild::Syntroph => {
                    host_counts[idx] = host_counts[idx].saturating_add(1);
                    syntroph_counts[idx] = syntroph_counts[idx].saturating_add(1);
                }
                TrophicGuild::PredatorParasite => {
                    predator_counts[idx] = predator_counts[idx].saturating_add(1);
                }
            }
        }

        let adami_rho = self.last_adami_rho;
        let rep_pop_n = if self.controls.replication_rate >= 0.06 {
            self.population.len()
        } else {
            1
        };

        let mut plastic_flux_sum = 0.0;
        let mut step_cross_feed_flux = 0.0;

        for unit in &mut self.population {
            unit.step_spatial_movement(
                &self.grid,
                &self.controls,
                &host_counts,
                &predator_counts,
                &mut rng,
            );
            let gx = unit.grid_x;
            let gy = unit.grid_y;
            let loc_e = self.grid.local_exergy(gx, gy);
            let loc_b = self.grid.local_byproduct(gx, gy);
            let loc_n = self.grid.local_niche(gx, gy);
            let biome = self.grid.biome_at(gx, gy);

            let mut local_syn_neighbors = 0u16;
            for ny in gy.saturating_sub(1)..=(gy + 1).min(self.grid.height - 1) {
                for nx in gx.saturating_sub(1)..=(gx + 1).min(self.grid.width - 1) {
                    local_syn_neighbors = local_syn_neighbors
                        .saturating_add(syntroph_counts[self.grid.idx(nx, ny)]);
                }
            }
            let syn_relief = (local_syn_neighbors as f64) * 0.28 * self.controls.cross_feeding_rate;

            let culls_before = unit.cumulative_apoptotic_culls;
            let pol_before = unit.cumulative_policed_eggs;
            let sh_before = unit.cumulative_antimicrobial_shields;
            let aud_before = unit.cumulative_peer_verifications;
            let prf_before = unit.cumulative_formal_proofs;
            let integ_before = unit.integrity;
            let m_before = if unit.lifetime_memory.is_empty() {
                0.0
            } else {
                unit.lifetime_memory.iter().sum::<f64>() / (unit.lifetime_memory.len() as f64)
            };

            unit.step_lifetime(
                &self.controls,
                loc_e,
                loc_b,
                loc_n,
                syn_relief,
                biome,
                adami_rho,
                rep_pop_n,
                self.env_phase,
                &mut rng,
            );

            if unit.cumulative_apoptotic_culls > culls_before {
                self.total_apoptotic_culls +=
                    (unit.cumulative_apoptotic_culls - culls_before) as usize;
            }
            if integ_before >= 0.30 && unit.integrity < 0.25 && unit.somatic_cheater_load >= 0.70 {
                self.total_neoplastic_collapses += 1;
            }
            if unit.cumulative_policed_eggs > pol_before {
                self.total_policed_eggs += (unit.cumulative_policed_eggs - pol_before) as usize;
            }
            if integ_before >= 0.30 && unit.integrity < 0.25 && unit.rebel_worker_load >= 0.70 {
                self.total_colony_anarchy_collapses += 1;
            }
            if unit.cumulative_antimicrobial_shields > sh_before {
                self.total_antimicrobial_shields +=
                    (unit.cumulative_antimicrobial_shields - sh_before) as usize;
            }
            if integ_before >= 0.30
                && unit.integrity < 0.25
                && unit.escovopsis_pathogen_load >= 0.70
            {
                self.total_garden_blight_collapses += 1;
            }
            if integ_before >= 0.30
                && unit.integrity < 0.25
                && unit.cultural_drift_erosion >= 0.70
            {
                self.total_cultural_drift_collapses += 1;
            }
            if unit.cumulative_peer_verifications > aud_before {
                self.total_peer_retraction_audits +=
                    (unit.cumulative_peer_verifications - aud_before) as usize;
            }
            if integ_before >= 0.30
                && unit.integrity < 0.25
                && unit.archive_false_discovery_load >= 0.70
            {
                self.total_archive_corruption_collapses += 1;
            }
            if unit.cumulative_formal_proofs > prf_before {
                self.total_formal_proofs_verified +=
                    (unit.cumulative_formal_proofs - prf_before) as usize;
            }
            if integ_before >= 0.30
                && unit.integrity < 0.25
                && unit.recursive_autophagy_load >= 0.70
            {
                self.total_recursive_model_collapses += 1;
            }
            step_cross_feed_flux += unit.syntrophic_work;
            let external_consumed_b = if unit.guild == TrophicGuild::Holobiont {
                unit.syntrophic_work * 0.25
            } else if unit.guild == TrophicGuild::PredatorParasite && unit.has_active_kleptoplast()
            {
                unit.syntrophic_work * 0.65
            } else {
                unit.syntrophic_work
            };
            self.grid.apply_trophic_flux(
                gx,
                gy,
                unit.primary_work,
                unit.byproduct_output,
                external_consumed_b,
                unit.niche_output,
            );

            let m_after = if unit.lifetime_memory.is_empty() {
                0.0
            } else {
                unit.lifetime_memory.iter().sum::<f64>() / (unit.lifetime_memory.len() as f64)
            };
            let dm = (m_after - m_before).max(0.0);
            let c_ratio = ((unit.complexity_c - unit.base_complexity_c).max(0.0))
                / unit.base_complexity_c.max(0.1);
            plastic_flux_sum += c_ratio + 2.8 * dm * self.controls.memory_plasticity;
        }
        let step_plastic = plastic_flux_sum / (self.population.len().max(1) as f64);
        self.smoothed_plastic_signal = 0.68 * self.smoothed_plastic_signal + 0.32 * step_plastic;

        // 3a. Spatial Multi-Trophic Encounters + Mixotrophic Kleptoplasty + Apicoplast Invasion!
        let mut step_trophic_flux = 0.0;
        let mut klepto_events_this_step = 0usize;
        if self.controls.bath_exergy_flux > 0.03
            && self.controls.trophic_coupling >= 0.05
            && self.population.len() >= 2
        {
            let pop_len = self.population.len();
            let mut work_deltas = vec![0.0_f64; pop_len];
            let mut best_prey_for_pred: Vec<Option<usize>> = vec![None; pop_len];
            let a_pred = self.controls.trophic_coupling.clamp(0.0, 1.5);
            let eta_trophic = self.controls.trophic_efficiency.clamp(0.08, 0.80);
            let handling_time = 0.55;

            for p_idx in 0..pop_len {
                if self.population[p_idx].guild != TrophicGuild::PredatorParasite {
                    continue;
                }
                let px = self.population[p_idx].grid_x;
                let py = self.population[p_idx].grid_y;

                let mut prey_list: Vec<(usize, f64)> = Vec::new();
                let mut available_prey_work = 0.0;
                let mut top_prey: Option<(usize, f64)> = None;
                for h_idx in 0..pop_len {
                    if p_idx == h_idx
                        || self.population[h_idx].guild == TrophicGuild::PredatorParasite
                    {
                        continue;
                    }
                    let hx = self.population[h_idx].grid_x;
                    let hy = self.population[h_idx].grid_y;
                    if hx.abs_diff(px) <= 1 && hy.abs_diff(py) <= 1 {
                        let susc =
                            self.population[h_idx].susceptibility_to(&self.population[p_idx]);
                        let avail = (self.population[h_idx].harvested_work
                            - self.population[h_idx].landauer_cost * 1.05)
                            .max(0.0);
                        if avail > 0.02 {
                            let eff_prey = avail * susc;
                            available_prey_work += eff_prey;
                            prey_list.push((h_idx, eff_prey));
                            if top_prey.map_or(true, |(_, w)| eff_prey > w) {
                                top_prey = Some((h_idx, eff_prey));
                            }
                        }
                    }
                }

                if available_prey_work > 1e-4 && !prey_list.is_empty() {
                    best_prey_for_pred[p_idx] = top_prey.map(|(i, _)| i);
                    let total_extracted = (a_pred * available_prey_work)
                        / (1.0 + a_pred * handling_time * available_prey_work);
                    for (h_idx, eff_w) in prey_list {
                        let share = eff_w / available_prey_work;
                        let loss = (total_extracted * share)
                            .min(0.45 * self.population[h_idx].harvested_work);
                        work_deltas[h_idx] -= loss;
                        let assimilated_gain = loss * eta_trophic;
                        work_deltas[p_idx] += assimilated_gain;
                        step_trophic_flux += assimilated_gain;
                    }
                }
            }

            for (idx, delta) in work_deltas.into_iter().enumerate() {
                if delta.abs() > 1e-6 {
                    let capped = if delta < 0.0 {
                        delta.max(-0.55 * self.population[idx].harvested_work)
                    } else {
                        delta
                    };
                    self.population[idx].apply_trophic_encounter_delta(capped, &self.controls);
                }
            }

            if self.controls.kleptoplasty_rate >= 0.03 || self.controls.klepto_hgt_rate >= 0.03 {
                let mut tipping_transitions: Vec<(
                    u64,
                    u64,
                    Option<u64>,
                    usize,
                    BiomeZone,
                    u8,
                    f64,
                    usize,
                    usize,
                    usize,
                    f64,
                )> = Vec::new();

                for (p_idx, maybe_prey_idx) in best_prey_for_pred.into_iter().enumerate() {
                    let Some(h_idx) = maybe_prey_idx else {
                        continue;
                    };
                    if self.population[p_idx].has_apicoplast()
                        || self.population[p_idx].apicoplast_eliminated
                    {
                        continue;
                    }
                    let do_steal = rng.bernoulli(self.controls.kleptoplasty_rate * 0.72);
                    let do_hgt = rng.bernoulli(self.controls.klepto_hgt_rate * 0.42);
                    if !do_steal && !do_hgt {
                        continue;
                    }
                    let prey_snap = self.population[h_idx].clone();
                    let old_pred_id = self.population[p_idx].id;
                    let (_stole, did_hgt, crossed) = self.population[p_idx]
                        .steal_kleptoplast_from_prey(&prey_snap, do_steal, do_hgt, true, &mut rng);
                    if did_hgt {
                        klepto_events_this_step += 1;
                    }
                    if crossed {
                        let new_holo_id = self.next_unit_id;
                        self.next_unit_id += 1;
                        let gx = self.population[p_idx].grid_x;
                        let gy = self.population[p_idx].grid_y;
                        let biome = self.grid.biome_at(gx, gy);
                        let le = self.grid.local_exergy(gx, gy);
                        let lb = self.grid.local_byproduct(gx, gy);
                        let ln = self.grid.local_niche(gx, gy);
                        self.population[p_idx].id = new_holo_id;
                        self.population[p_idx].recompute_phenome(
                            &self.controls,
                            le,
                            lb,
                            ln,
                            0.30,
                            biome,
                            adami_rho,
                            pop_len,
                            self.env_phase,
                        );
                        let u_ref = &self.population[p_idx];
                        let gamma = (u_ref.endosymbiont_sequence.len() as f64)
                            / ((u_ref.sequence.len() + u_ref.endosymbiont_sequence.len()).max(1)
                                as f64);
                        let lev = u_ref.complexity_c / u_ref.structure_s.max(2.0);
                        tipping_transitions.push((
                            new_holo_id,
                            old_pred_id,
                            u_ref.endosymbiont_origin_id,
                            u_ref.species_id,
                            biome,
                            u_ref.endo_tier,
                            gamma,
                            u_ref.klepto_hgt_loci,
                            u_ref.endosymbiont_sequence.len(),
                            u_ref.nucleomorph_sequence.len(),
                            lev,
                        ));
                        klepto_events_this_step += 1;
                    }
                }

                for (hid, pid, did, sp, bm, tier, gm, hgt_l, org_rem, nm_l, lv) in
                    tipping_transitions
                {
                    self.phylogeny.record_klepto_tipping_fixation(
                        hid,
                        pid,
                        did,
                        self.step_count,
                        sp,
                        bm,
                        tier,
                        gm,
                        hgt_l,
                        org_rem,
                        nm_l,
                        lv,
                    );
                }
            }
        }
        self.last_trophic_flux = 0.60 * self.last_trophic_flux + 0.40 * step_trophic_flux;
        self.last_cross_feed_flux = 0.60 * self.last_cross_feed_flux + 0.40 * step_cross_feed_flux;

        // 3b. Serial `1°/2°/3°` Endosymbiotic Fusions + Apicoplast & Anaerobic MRO Transitions!
        let mut fusion_events_this_step = klepto_events_this_step;
        if self.controls.bath_exergy_flux > 0.03
            && self.controls.endosymbiosis_rate >= 0.04
            && self.population.len() >= 4
        {
            let cap = self.controls.carrying_capacity.max(1);
            let current_holos = self
                .population
                .iter()
                .filter(|u| u.guild == TrophicGuild::Holobiont)
                .count();
            let current_prods = self
                .population
                .iter()
                .filter(|u| u.guild == TrophicGuild::Producer)
                .count();
            let current_syns = self
                .population
                .iter()
                .filter(|u| u.guild == TrophicGuild::Syntroph)
                .count();
            let current_apicos = self
                .population
                .iter()
                .filter(|u| u.has_apicoplast() || u.apicoplast_eliminated)
                .count();

            if self.controls.apicoplast_reduction_rate >= 0.04
                && current_holos >= 3
                && current_apicos < (cap / 5).max(3)
                && rng.bernoulli((self.controls.apicoplast_reduction_rate * 0.24).clamp(0.0, 0.35))
            {
                if let Some((hi, _)) = self
                    .population
                    .iter()
                    .enumerate()
                    .find(|(_, u)| u.guild == TrophicGuild::Holobiont)
                {
                    let old_holo_id = self.population[hi].id;
                    if self.population[hi].convert_holobiont_to_apicoplast_parasite(&mut rng) {
                        let new_apico_id = self.next_unit_id;
                        self.next_unit_id += 1;
                        let gx = self.population[hi].grid_x;
                        let gy = self.population[hi].grid_y;
                        let biome = self.grid.biome_at(gx, gy);
                        let le = self.grid.local_exergy(gx, gy);
                        let lb = self.grid.local_byproduct(gx, gy);
                        let ln = self.grid.local_niche(gx, gy);
                        let pop_n = self.population.len();
                        self.population[hi].id = new_apico_id;
                        self.population[hi].recompute_phenome(
                            &self.controls,
                            le,
                            lb,
                            ln,
                            0.28,
                            biome,
                            adami_rho,
                            pop_n,
                            self.env_phase,
                        );
                        let u_ref = &self.population[hi];
                        let lev = u_ref.complexity_c / u_ref.structure_s.max(2.0);
                        self.phylogeny.record_apicoplast_transition(
                            new_apico_id,
                            old_holo_id,
                            self.step_count,
                            u_ref.species_id,
                            biome,
                            u_ref.endo_tier,
                            u_ref.egt_transferred_loci,
                            u_ref.endosymbiont_sequence.len(),
                            u_ref.apicoplast_pathway_count(),
                            lev,
                        );
                        fusion_events_this_step += 1;
                    }
                }
            }

            let mut did_serial = false;
            if self.controls.serial_endosymbiosis_rate >= 0.04 && current_holos >= 2 {
                let mut serial_pair: Option<(usize, usize)> = None;
                for (hi, host_u) in self.population.iter().enumerate() {
                    if !matches!(
                        host_u.guild,
                        TrophicGuild::Holobiont | TrophicGuild::PredatorParasite
                    ) || host_u.endo_tier >= 3
                        || host_u.has_apicoplast()
                    {
                        continue;
                    }
                    if let Some((ei, _)) = self.population.iter().enumerate().find(|(j, cand)| {
                        *j != hi
                            && cand.guild == TrophicGuild::Holobiont
                            && cand.endo_tier <= 2
                            && cand.grid_x.abs_diff(host_u.grid_x) <= 2
                            && cand.grid_y.abs_diff(host_u.grid_y) <= 2
                    }) {
                        let p_ser =
                            (self.controls.serial_endosymbiosis_rate * 0.28).clamp(0.0, 0.42);
                        if rng.bernoulli(p_ser) {
                            serial_pair = Some((hi, ei));
                            break;
                        }
                    }
                }

                if let Some((hi, ei)) = serial_pair {
                    let host = self.population[hi].clone();
                    let engulfed = self.population[ei].clone();
                    let cid = self.next_unit_id;
                    self.next_unit_id += 1;
                    let biome = self.grid.biome_at(host.grid_x, host.grid_y);
                    let le = self.grid.local_exergy(host.grid_x, host.grid_y);
                    let lb = self.grid.local_byproduct(host.grid_x, host.grid_y);
                    let ln = self.grid.local_niche(host.grid_x, host.grid_y);
                    let pop_n = self.population.len();
                    let (chimera, gamma) = ModeledUnit::fuse_serial_holobiont(
                        &host,
                        &engulfed,
                        cid,
                        biome,
                        le,
                        lb,
                        ln,
                        &self.controls,
                        adami_rho,
                        pop_n,
                        self.env_phase,
                    );
                    let lev = chimera.complexity_c / chimera.structure_s.max(2.0);
                    self.phylogeny.record_serial_endosymbiotic_fusion(
                        chimera.id,
                        host.id,
                        engulfed.id,
                        self.step_count,
                        chimera.species_id,
                        biome,
                        chimera.endo_tier,
                        gamma,
                        chimera.endosymbiont_sequence.len(),
                        chimera.nucleomorph_sequence.len(),
                        chimera.egt_transferred_loci,
                        lev,
                    );
                    let rem_hi = hi.max(ei);
                    let rem_lo = hi.min(ei);
                    self.population.remove(rem_hi);
                    self.population.remove(rem_lo);
                    self.population.push(chimera);
                    fusion_events_this_step += 1;
                    did_serial = true;
                }
            }

            if !did_serial
                && current_holos < (cap / 3).max(4)
                && current_prods >= 3
                && current_syns >= 2
            {
                let mut candidate_pair: Option<(usize, usize)> = None;
                for (pi, u) in self.population.iter().enumerate() {
                    if u.guild != TrophicGuild::Producer
                        || u.harvested_work <= u.landauer_cost * 1.15
                    {
                        continue;
                    }
                    if let Some((si, _)) =
                        self.population.iter().enumerate().find(|(j, cand)| {
                            *j != pi
                                && cand.guild == TrophicGuild::Syntroph
                                && cand.grid_x.abs_diff(u.grid_x) <= 2
                                && cand.grid_y.abs_diff(u.grid_y) <= 2
                        })
                    {
                        let loc_n = self.grid.local_niche(u.grid_x, u.grid_y);
                        let loc_b = self.grid.local_byproduct(u.grid_x, u.grid_y);
                        let p_fuse = (self.controls.endosymbiosis_rate
                            * (0.24 + 0.35 * loc_n + 0.25 * loc_b))
                            .clamp(0.0, 0.55);
                        if rng.bernoulli(p_fuse) {
                            candidate_pair = Some((pi, si));
                            break;
                        }
                    }
                }

                if let Some((pi, si)) = candidate_pair {
                    let host = self.population[pi].clone();
                    let sym = self.population[si].clone();
                    let hid = self.next_unit_id;
                    self.next_unit_id += 1;
                    let biome = self.grid.biome_at(host.grid_x, host.grid_y);
                    let le = self.grid.local_exergy(host.grid_x, host.grid_y);
                    let lb = self.grid.local_byproduct(host.grid_x, host.grid_y);
                    let ln = self.grid.local_niche(host.grid_x, host.grid_y);
                    let pop_n = self.population.len();
                    let (holo, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
                        &host,
                        &sym,
                        hid,
                        biome,
                        le,
                        lb,
                        ln,
                        &self.controls,
                        adami_rho,
                        pop_n,
                        self.env_phase,
                    );
                    let lev = holo.complexity_c / holo.structure_s.max(2.0);
                    self.phylogeny.record_endosymbiotic_fusion(
                        holo.id,
                        host.id,
                        sym.id,
                        self.step_count,
                        holo.species_id,
                        biome,
                        gamma,
                        holo.endosymbiont_sequence.len(),
                        lev,
                    );
                    let hi = pi.max(si);
                    let lo = pi.min(si);
                    self.population.remove(hi);
                    self.population.remove(lo);
                    self.population.push(holo);
                    fusion_events_this_step += 1;
                }
            }
        }

        // 4. Spatially Local Replication + Germline EGT/NM/Apico/MRO Reduction + Price Decomposition
        let start_pop_len = self.population.len().max(1);
        let mut births_this_step = 0usize;
        let mut assim_events = 0usize;
        let mut recomb_events = fusion_events_this_step;
        let mut price_records: Vec<(f64, f64, f64)> = Vec::new();

        let mean_margin = self
            .population
            .iter()
            .map(|u| u.harvested_work - u.landauer_cost)
            .sum::<f64>()
            / (start_pop_len as f64);

        if self.controls.bath_exergy_flux > 0.03
            && mean_margin > 0.01
            && self.controls.replication_rate >= 0.06
        {
            let cap = self.controls.carrying_capacity.max(1);
            let crowding = if self.controls.interaction_strength >= 0.20 {
                ((start_pop_len as f64) / (cap as f64)).clamp(0.0, 1.5)
            } else {
                0.0
            };

            let mean_raw_s = self.mean_raw_structure_s();
            let theta_f =
                compute_critical_fidelity(mean_raw_s, self.controls.selective_superiority);
            let f_eff = self.controls.effective_fidelity();
            let darwinian_ratchet = (self.controls.heritable_variation >= 0.02
                || self.controls.soft_inheritance_coupling >= 0.10
                || self.controls.recombination_rate >= 0.10
                || self.controls.endosymbiosis_rate >= 0.05
                || self.controls.klepto_hgt_rate >= 0.05)
                && self.controls.interaction_strength >= 0.20
                && f_eff >= theta_f;

            let m_fit = self.mean_fitness().max(0.05);
            let snapshot_parents = self.population.clone();
            let mut newborns = Vec::new();
            let mut f_obs_sum = 0.0;

            for (idx, parent) in snapshot_parents.iter().enumerate() {
                if start_pop_len + newborns.len() >= cap * 2 {
                    break;
                }
                let unit_margin = parent.harvested_work - parent.landauer_cost;
                let z_p = parent.heritable_leverage(&self.controls);
                if unit_margin <= 0.0 {
                    price_records.push((0.05, z_p, z_p));
                    continue;
                }

                let streamline_speed =
                    (24.0 / (parent.sequence.len().max(10) as f64)).clamp(0.85, 1.30);
                let fit_factor = if darwinian_ratchet {
                    ((parent.fitness / m_fit) * streamline_speed).clamp(0.20, 2.40)
                } else {
                    1.0
                };
                let effective_birth = self.controls.replication_rate
                    * fit_factor
                    * (1.0 - 0.55 * crowding * self.controls.interaction_strength)
                        .clamp(0.12, 1.0);

                if rng.bernoulli(effective_birth) {
                    let cid = self.next_unit_id;
                    self.next_unit_id += 1;

                    let (cx, cy) = if darwinian_ratchet && rng.bernoulli(0.06) {
                        (
                            rng.gen_index(self.grid.width),
                            rng.gen_index(self.grid.height),
                        )
                    } else {
                        let dx = (rng.gen_index(3) as isize) - 1;
                        let dy = (rng.gen_index(3) as isize) - 1;
                        (
                            (parent.grid_x as isize + dx)
                                .clamp(0, (self.grid.width - 1) as isize)
                                as usize,
                            (parent.grid_y as isize + dy)
                                .clamp(0, (self.grid.height - 1) as isize)
                                as usize,
                        )
                    };

                    let donor = if snapshot_parents.len() >= 2
                        && self.controls.recombination_rate >= 0.05
                    {
                        snapshot_parents
                            .iter()
                            .enumerate()
                            .find(|(j, cand)| {
                                *j != idx
                                    && cand.grid_x.abs_diff(parent.grid_x)
                                        + cand.grid_y.abs_diff(parent.grid_y)
                                        <= 5
                                    && normalized_sequence_divergence(
                                        &parent.sequence,
                                        &cand.sequence,
                                    ) <= self.controls.speciation_threshold
                            })
                            .map(|(_, cand)| cand)
                    } else {
                        None
                    };

                    let c_biome = self.grid.biome_at(cx, cy);
                    let c_le = self.grid.local_exergy(cx, cy);
                    let c_lb = self.grid.local_byproduct(cx, cy);
                    let c_ln = self.grid.local_niche(cx, cy);

                    if self.population[idx].has_active_kleptoplast() {
                        let shield =
                            ((self.population[idx].klepto_hgt_loci as f64) / 4.0).clamp(0.0, 1.0);
                        self.population[idx].kleptoplast_integrity *= 0.56 + 0.12 * shield;
                        if self.population[idx].kleptoplast_integrity < 0.06 {
                            self.population[idx].kleptoplast_sequence.clear();
                            self.population[idx].kleptoplast_integrity = 0.0;
                            self.population[idx].kleptoplast_donor_id = None;
                        }
                    }

                    let (child, f_obs, assim, recomb) = parent.replicate_with_partner(
                        donor,
                        cid,
                        cx,
                        cy,
                        c_biome,
                        c_le,
                        c_lb,
                        c_ln,
                        &self.controls,
                        adami_rho,
                        start_pop_len,
                        self.env_phase,
                        &mut rng,
                    );
                    let z_c = child.heritable_leverage(&self.controls);
                    price_records.push((1.0 + fit_factor, z_p, z_c));

                    let seq_div =
                        normalized_sequence_divergence(&parent.sequence, &child.sequence);
                    let donor_uid = if recomb { donor.map(|d| d.id) } else { None };
                    self.phylogeny.record_birth_full(
                        child.id,
                        parent.id,
                        donor_uid,
                        self.step_count,
                        child.species_id,
                        child.guild,
                        c_biome,
                        parent.species_id,
                        parent.guild,
                        seq_div,
                        self.controls.speciation_threshold,
                        z_c,
                        child.egt_transferred_loci,
                        child.endosymbiont_sequence.len(),
                        child.endo_tier,
                        child.nucleomorph_sequence.len(),
                        child.nucleomorph_eliminated && !parent.nucleomorph_eliminated,
                    );

                    births_this_step += 1;
                    if assim {
                        assim_events += 1;
                    }
                    if recomb
                        || child.egt_transferred_loci > parent.egt_transferred_loci
                        || child.nucleomorph_transferred_loci > parent.nucleomorph_transferred_loci
                        || (child.has_apicoplast() && !parent.has_apicoplast())
                        || (child.mro_class != parent.mro_class)
                    {
                        recomb_events += 1;
                    }
                    f_obs_sum += f_obs;
                    newborns.push(child);
                } else {
                    price_records.push((0.5 * fit_factor, z_p, z_p));
                }
            }

            if births_this_step > 0 {
                let batch_f = f_obs_sum / (births_this_step as f64);
                self.last_observed_fidelity = Some(match self.last_observed_fidelity {
                    Some(prev) => 0.65 * prev + 0.35 * batch_f,
                    None => batch_f,
                });
            }

            if price_records.len() >= 2 {
                let step_price = PriceDecomposition::from_parent_offspring_records(&price_records);
                self.last_price = PriceDecomposition {
                    selection_covariance: 0.65 * self.last_price.selection_covariance
                        + 0.35 * step_price.selection_covariance,
                    transmission_bias: 0.65 * self.last_price.transmission_bias
                        + 0.35 * step_price.transmission_bias,
                    selection_share: 0.65 * self.last_price.selection_share
                        + 0.35 * step_price.selection_share,
                };
            }

            self.population.extend(newborns);

            if self.population.len() > cap {
                if darwinian_ratchet {
                    let max_pred = if self.controls.trophic_coupling >= 0.05 {
                        (cap / 4).max(2)
                    } else {
                        0
                    };
                    let max_syn = if self.controls.cross_feeding_rate >= 0.05 {
                        (cap / 4).max(2)
                    } else {
                        0
                    };
                    let max_holo = (cap / 3).max(4);

                    let mut prod_alpha = Vec::new();
                    let mut prod_beta = Vec::new();
                    let mut syn_pool = Vec::new();
                    let mut holo_pool = Vec::new();
                    let mut pred_pool = Vec::new();

                    for u in self.population.drain(..) {
                        let score = (if self.controls.artificial_selection {
                            u.fitness * 1.35
                        } else {
                            u.fitness
                        }) + rng.range_f64(-0.02, 0.02);

                        match u.guild {
                            TrophicGuild::Producer => match self.grid.biome_at(u.grid_x, u.grid_y) {
                                BiomeZone::Alpha => prod_alpha.push((score, u)),
                                BiomeZone::Beta => prod_beta.push((score, u)),
                            },
                            TrophicGuild::Syntroph => syn_pool.push((score, u)),
                            TrophicGuild::Holobiont => holo_pool.push((score, u)),
                            TrophicGuild::PredatorParasite => {
                                if u.harvested_work > u.landauer_cost * 0.85 {
                                    pred_pool.push((score, u));
                                }
                            }
                        }
                    }

                    prod_alpha.sort_by(|a, b| b.0.total_cmp(&a.0));
                    prod_beta.sort_by(|a, b| b.0.total_cmp(&a.0));
                    syn_pool.sort_by(|a, b| b.0.total_cmp(&a.0));
                    holo_pool.sort_by(|a, b| b.0.total_cmp(&a.0));
                    pred_pool.sort_by(|a, b| b.0.total_cmp(&a.0));

                    let take_pred = pred_pool.len().min(max_pred);
                    let take_syn = syn_pool.len().min(max_syn);
                    let take_holo = holo_pool.len().min(max_holo);
                    let prod_slots = cap
                        .saturating_sub(take_pred + take_syn + take_holo)
                        .max(4);
                    let half_prod = (prod_slots / 2).max(1);

                    let mut survivors = Vec::with_capacity(cap);
                    let mut overflow = Vec::new();

                    for (i, item) in pred_pool.into_iter().enumerate() {
                        if i < take_pred {
                            survivors.push(item.1);
                        } else {
                            overflow.push(item);
                        }
                    }
                    for (i, item) in syn_pool.into_iter().enumerate() {
                        if i < take_syn {
                            survivors.push(item.1);
                        } else {
                            overflow.push(item);
                        }
                    }
                    for (i, item) in holo_pool.into_iter().enumerate() {
                        if i < take_holo {
                            survivors.push(item.1);
                        } else {
                            overflow.push(item);
                        }
                    }
                    let take_pa = prod_alpha.len().min(half_prod);
                    let take_pb = prod_beta.len().min(prod_slots.saturating_sub(take_pa));
                    for (i, item) in prod_alpha.into_iter().enumerate() {
                        if i < take_pa {
                            survivors.push(item.1);
                        } else {
                            overflow.push(item);
                        }
                    }
                    for (i, item) in prod_beta.into_iter().enumerate() {
                        if i < take_pb {
                            survivors.push(item.1);
                        } else {
                            overflow.push(item);
                        }
                    }
                    if survivors.len() < cap && !overflow.is_empty() {
                        overflow.sort_by(|a, b| b.0.total_cmp(&a.0));
                        let need = cap - survivors.len();
                        survivors.extend(overflow.into_iter().take(need).map(|(_, u)| u));
                    }
                    if survivors.len() > cap {
                        survivors.truncate(cap);
                    }
                    self.population = survivors;
                } else {
                    while self.population.len() > cap {
                        let victim = rng.gen_index(self.population.len());
                        self.population.swap_remove(victim);
                    }
                }
            }
        } else {
            self.smoothed_birth_rate *= 0.60;
        }

        let inst_birth = (births_this_step as f64) / (start_pop_len as f64);
        if self.controls.replication_rate >= 0.06 && self.controls.bath_exergy_flux > 0.03 {
            self.smoothed_birth_rate = 0.72 * self.smoothed_birth_rate + 0.28 * inst_birth.max(0.05);
            let b_denom = (births_this_step + fusion_events_this_step).max(1) as f64;
            let inst_assim = if births_this_step > 0 {
                (assim_events as f64) / b_denom
            } else {
                self.smoothed_assimilation
            };
            let inst_recomb = if births_this_step + fusion_events_this_step > 0 {
                (recomb_events as f64) / b_denom
            } else {
                self.smoothed_recombination
            };
            self.smoothed_assimilation = 0.72 * self.smoothed_assimilation + 0.28 * inst_assim;
            self.smoothed_recombination = 0.72 * self.smoothed_recombination + 0.28 * inst_recomb;
        } else {
            self.smoothed_birth_rate = 0.0;
            self.smoothed_assimilation = 0.0;
            self.smoothed_recombination = 0.0;
        }

        if self.population.is_empty() {
            let fid = self.alloc_id();
            let fallback = ModeledUnit::new_founder(
                fid,
                self.organism,
                self.target_regime,
                4,
                4,
                BiomeZone::Alpha,
                self.env_phase,
                &mut rng,
            );
            self.population.push(fallback);
            self.reset_phylogeny_from_population();
        }

        self.rng = rng;

        if self.coupled_biosphere_mode
            && self.controls.biogeochemical_feedback_gain > 0.001
            && !self.population.is_empty()
        {
            let pop_f = self.population.len().max(1) as f64;
            let (mut s_o2, mut s_rec, mut s_weath, mut s_tech, mut s_stew) =
                (0.0, 0.0, 0.0, 0.0, 0.0);
            for u in &mut self.population {
                let (ao, sr, sw, tb, vs) =
                    u.compute_planetary_guild_fluxes(&self.controls, &self.biosphere);
                s_o2 += ao;
                s_rec += sr;
                s_weath += sw;
                s_tech += tb;
                s_stew += vs;
            }
            let (_, burnout) = self.biosphere.step_biogeochemical_feedbacks(
                s_o2 / pop_f,
                s_rec / pop_f,
                s_weath / pop_f,
                s_tech / pop_f,
                s_stew / pop_f,
                &self.controls,
            );
            self.controls.ambient_oxygen_level = self.biosphere.atmospheric_o2;
            if burnout {
                for u in &mut self.population {
                    if u.apply_planetary_boundary_stress(&self.biosphere) {
                        self.total_planetary_burnout_collapses += 1;
                    }
                }
            }
        }

        self.record_metrics();
    }

    pub fn trigger_network_transition(&mut self) -> RegimeId {
        let edges = self.classified_regime.outgoing_transitions();
        if edges.is_empty() {
            return self.classified_regime;
        }
        let total_w: f64 = edges.iter().map(|(_, _, w)| *w).sum();
        let mut pick = self.rng.next_f64() * total_w;
        let mut chosen = edges[0].0;
        for &(target, _, w) in edges {
            if pick <= w {
                chosen = target;
                break;
            }
            pick -= w;
        }
        self.apply_regime_preset(chosen, chosen.index() >= 2);
        self.record_metrics();
        chosen
    }

    pub fn best_unit(&self) -> Option<&ModeledUnit> {
        self.population
            .iter()
            .max_by(|a, b| a.fitness.total_cmp(&b.fitness))
    }

    pub fn phylo_summary(&self) -> PhyloSummary {
        self.phylogeny.compute_summary(self.step_count)
    }

    pub fn export_newick(&self, include_nhx: bool) -> String {
        self.phylogeny.to_newick(self.step_count, include_nhx)
    }

    pub fn export_phylogram_svg(&self) -> String {
        let title = format!(
            "ITME v12.0 Anaerobic MROs (Hydrogenosomes & Mitosomes), Apicoplasts & Symbiogenesis — {} @ {} (Step {})",
            self.organism.title(),
            self.classified_regime.full_title(),
            self.step_count
        );
        self.phylogeny.to_svg_phylogram(&title, self.step_count)
    }

    pub fn current_phylogram_layout(&self) -> PhylogramLayout {
        self.phylogeny.extract_layout(self.step_count)
    }

    fn mean_raw_structure_s(&self) -> f64 {
        if self.population.is_empty() {
            return 8.0;
        }
        self.population
            .iter()
            .map(|u| u.structure_raw_s)
            .sum::<f64>()
            / (self.population.len() as f64)
    }

    fn mean_fitness(&self) -> f64 {
        if self.population.is_empty() {
            return 0.1;
        }
        self.population.iter().map(|u| u.fitness).sum::<f64>() / (self.population.len() as f64)
    }

    fn record_metrics(&mut self) {
        let sigma_size = self.organism.alphabet_symbols().len();

        let (species_count, fst_index, adami_rho) =
            if self.controls.replication_rate >= 0.06 && self.population.len() >= 2 {
                let refs: Vec<&[u16]> = self
                    .population
                    .iter()
                    .map(|u| u.sequence.as_slice())
                    .collect();
                let unit_biomes: Vec<BiomeZone> = self
                    .population
                    .iter()
                    .map(|u| self.grid.biome_at(u.grid_x, u.grid_y))
                    .collect();
                let spec = analyze_spatial_speciation(
                    &refs,
                    &unit_biomes,
                    sigma_size,
                    self.controls.speciation_threshold,
                );
                for (u, sp) in self.population.iter_mut().zip(spec.unit_species.into_iter()) {
                    u.species_id = sp;
                }
                (
                    spec.species_count,
                    spec.fst_index,
                    spec.within_species_adami_rho,
                )
            } else {
                for u in &mut self.population {
                    u.species_id = 0;
                }
                (1, 0.0, 1.0)
            };

        self.last_adami_rho = adami_rho;
        self.last_species_count = species_count;
        self.last_fst = fst_index;

        let extant_states: Vec<ExtantPhyloState> = self
            .population
            .iter()
            .map(|u| {
                let b = self.grid.biome_at(u.grid_x, u.grid_y);
                let lev = u.complexity_c / u.structure_s.max(2.0);
                ExtantPhyloState {
                    taxon_id: u.id,
                    species_id: u.species_id,
                    guild: u.guild,
                    biome: b,
                    leverage_li: lev,
                    egt_loci: u.egt_transferred_loci,
                    organelle_rem: u.endosymbiont_sequence.len(),
                    endo_tier: u.endo_tier,
                    nucleomorph_loci: u.nucleomorph_sequence.len(),
                    nucleomorph_eliminated: u.nucleomorph_eliminated,
                    klepto_hgt_loci: u.klepto_hgt_loci,
                    has_kleptoplast: u.has_active_kleptoplast(),
                    has_apicoplast: u.has_apicoplast(),
                    apicoplast_pathways: u.apicoplast_pathway_count(),
                    apicoplast_eliminated: u.apicoplast_eliminated,
                    mro_class: u.mro_class,
                    mro_pathways: u.mro_pathway_count(),
                    mro_mtdna_loci: u.mro_mtdna_loci,
                    multicell_state: u.multicell_state,
                    cell_count: u.cell_count,
                    soma_fraction: u.soma_fraction,
                    eusocial_state: u.eusocial_state,
                    colony_size: u.colony_size,
                    sterile_caste_fraction: u.total_sterile_caste_fraction(),
                    agriculture_state: u.agriculture_state,
                    cultivar_polyploidy: u.cultivar_polyploidy,
                    gongylidia_yield: u.gongylidia_yield,
                    cultural_state: u.cultural_state,
                    cultural_ratchet_tier: u.cultural_ratchet_tier,
                    cultural_repertoire_bits: u.cultural_repertoire_bits(),
                    exosomatic_state: u.exosomatic_state,
                    exosomatic_archive_bits: u.exosomatic_archive_bits(),
                    machine_synthesis_state: u.machine_synthesis_state,
                    verified_theorem_bits: u.verified_theorem_bits(),
                    distribution_tail_diversity: u.distribution_tail_diversity,
                    planetary_epoch: u.planetary_epoch,
                    chemical_disequilibrium_bits: self.biosphere.chemical_disequilibrium_bits,
                    abiogenesis_stage: u.abiogenesis_stage,
                    raf_closure_score: u.raf_closure_score,
                    chirality_code_stage: u.chirality_code_stage,
                    genetic_code_mp: u.genetic_code_mp,
                    genome_substrate_stage: u.genome_substrate_stage,
                    chromosome_capacity_kb: u.chromosome_capacity_kb,
                }
            })
            .collect();
        self.phylogeny
            .sync_extant_population_full(self.step_count, &extant_states);
        let phylo_sum = self.phylogeny.compute_summary(self.step_count);

        let mut w_sum = 0.0;
        let mut s_raw_sum = 0.0;
        let mut s_eff_sum = 0.0;
        let mut c_sum = 0.0;
        let mut li_sum = 0.0;
        let mut margin_sum = 0.0;
        let mut landauer_sum = 0.0;
        let mut entropy_prod_sum = 0.0;
        let mut producer_count = 0usize;
        let mut syntroph_count = 0usize;
        let mut predator_count = 0usize;
        let mut aerobic_mito_count = 0usize;
        let mut hydrogenosome_count = 0usize;
        let mut mitosome_count = 0usize;
        let mut amiteukaryote_count = 0usize;
        let mut clonal_multicell_count = 0usize;
        let mut germ_soma_count = 0usize;
        let mut cell_count_sum = 0.0;
        let mut soma_frac_sum = 0.0;
        let mut cheater_load_sum = 0.0;
        let mut multicell_bearers = 0usize;
        let mut eusocial_colony_count = 0usize;
        let mut superorganism_count = 0usize;
        let mut colony_size_sum = 0.0;
        let mut caste_frac_sum = 0.0;
        let mut rebel_load_sum = 0.0;
        let mut social_bearers = 0usize;
        let mut farming_colony_count = 0usize;
        let mut obligate_agro_count = 0usize;
        let mut gongylidia_sum = 0.0;
        let mut escovopsis_sum = 0.0;
        let mut cultural_lineage_count = 0usize;
        let mut symbolic_culture_count = 0usize;
        let mut cultural_bits_sum = 0.0;
        let mut neural_cap_sum = 0.0;
        let mut archival_lineage_count = 0usize;
        let mut verified_institution_count = 0usize;
        let mut archive_bits_sum = 0.0;
        let mut archive_fdr_sum = 0.0;
        let mut computational_symbiont_count = 0usize;
        let mut formal_symbio_count = 0usize;
        let mut theorem_bits_sum = 0.0;
        let mut tail_div_sum = 0.0;
        let mut mro_pathways_sum = 0.0;
        let mut mro_bearers = 0usize;
        let mut apicoplast_parasite_count = 0usize;
        let mut apico_pathways_sum = 0.0;
        let mut kleptoplastic_host_count = 0usize;
        let mut kp_integrity_sum = 0.0;
        let mut klepto_hgt_sum = 0.0;
        let mut klepto_hgt_bearers = 0usize;
        let mut mosaic_donor_sum = 0.0;
        let mut mosaic_bearers = 0usize;
        let mut holobiont_count = 0usize;
        let mut enslaved_holobiont_count = 0usize;
        let mut secondary_holobiont_count = 0usize;
        let mut tertiary_holobiont_count = 0usize;
        let mut active_nucleomorph_count = 0usize;
        let mut nm_len_sum = 0.0;
        let mut egt_frac_sum = 0.0;
        let mut org_len_sum = 0.0;
        let mut tl_sum = 0.0;

        for u in &self.population {
            let w = u.maturation_weight();
            w_sum += w;
            s_raw_sum += w * u.structure_raw_s;
            s_eff_sum += w * u.structure_s;
            c_sum += w * u.complexity_c;
            li_sum += w * (u.complexity_c / u.structure_s.max(2.0));

            let margin = u.harvested_work - u.landauer_cost;
            margin_sum += margin;
            landauer_sum += u.landauer_cost;
            entropy_prod_sum += u.landauer_cost + u.harvested_work.abs() * 0.42;

            colony_size_sum += u.colony_size.max(1) as f64;
            if u.eusocial_state.is_social() {
                social_bearers += 1;
                caste_frac_sum += u.total_sterile_caste_fraction();
                rebel_load_sum += u.rebel_worker_load;
                if u.eusocial_state.is_eusocial() {
                    eusocial_colony_count += 1;
                }
                if u.eusocial_state == EusocialState::ObligateSuperorganism {
                    superorganism_count += 1;
                }
            }
            if u.agriculture_state.is_farming() {
                farming_colony_count += 1;
                gongylidia_sum += u.gongylidia_yield;
                escovopsis_sum += u.escovopsis_pathogen_load;
                if u.agriculture_state.is_obligate_monoculture() {
                    obligate_agro_count += 1;
                }
            }
            if u.cultural_state.is_cultural() {
                cultural_lineage_count += 1;
                cultural_bits_sum += u.cultural_repertoire_bits();
                neural_cap_sum += u.neural_learning_capacity;
                if u.cultural_state == CulturalState::SymbolicCumulativeCulture {
                    symbolic_culture_count += 1;
                }
            }
            if u.exosomatic_state.is_archival() {
                archival_lineage_count += 1;
                archive_bits_sum += u.exosomatic_archive_bits();
                archive_fdr_sum += u.archive_false_discovery_load;
                if u.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                    verified_institution_count += 1;
                }
            }
            if u.machine_synthesis_state.has_computational_symbiont() {
                computational_symbiont_count += 1;
                theorem_bits_sum += u.verified_theorem_bits();
                tail_div_sum += u.distribution_tail_diversity;
                if u.machine_synthesis_state.is_formally_verified() {
                    formal_symbio_count += 1;
                }
            }
            cell_count_sum += u.cell_count.max(1) as f64;
            if u.multicell_state.is_multicellular() {
                multicell_bearers += 1;
                soma_frac_sum += u.soma_fraction;
                cheater_load_sum += u.somatic_cheater_load;
                if u.multicell_state.is_clonal() {
                    clonal_multicell_count += 1;
                }
                if u.multicell_state == MulticellularState::GermSomaMulticellular {
                    germ_soma_count += 1;
                }
            }

            match u.mro_class {
                MroClass::None => {}
                MroClass::AerobicMitochondrion | MroClass::AnaerobicMitochondrion => {
                    aerobic_mito_count += 1;
                    mro_pathways_sum += u.mro_pathway_count() as f64;
                    mro_bearers += 1;
                }
                MroClass::Hydrogenosome => {
                    hydrogenosome_count += 1;
                    mro_pathways_sum += u.mro_pathway_count() as f64;
                    mro_bearers += 1;
                }
                MroClass::Mitosome => {
                    mitosome_count += 1;
                    mro_pathways_sum += u.mro_pathway_count() as f64;
                    mro_bearers += 1;
                }
                MroClass::Amiteukaryote => {
                    amiteukaryote_count += 1;
                    mro_bearers += 1;
                }
            }

            if (u.klepto_hgt_loci > 0 || u.guild == TrophicGuild::PredatorParasite)
                && !u.has_apicoplast()
            {
                klepto_hgt_sum += u.klepto_hgt_loci as f64;
                klepto_hgt_bearers += 1;
            }
            if u.mosaic_donor_count() > 0 {
                mosaic_donor_sum += u.mosaic_donor_count() as f64;
                mosaic_bearers += 1;
            }

            match u.guild {
                TrophicGuild::Producer => {
                    producer_count += 1;
                    tl_sum += u.guild.trophic_level();
                }
                TrophicGuild::Syntroph => {
                    syntroph_count += 1;
                    tl_sum += u.guild.trophic_level();
                }
                TrophicGuild::PredatorParasite => {
                    predator_count += 1;
                    if u.has_apicoplast() {
                        apicoplast_parasite_count += 1;
                        apico_pathways_sum += u.apicoplast_pathway_count() as f64;
                        tl_sum += 2.08;
                    } else if u.has_active_kleptoplast() {
                        kleptoplastic_host_count += 1;
                        kp_integrity_sum += u.kleptoplast_integrity;
                        tl_sum += 1.65;
                    } else {
                        tl_sum += u.guild.trophic_level();
                    }
                }
                TrophicGuild::Holobiont => {
                    holobiont_count += 1;
                    if u.is_organelle_enslaved() {
                        enslaved_holobiont_count += 1;
                    }
                    if u.endo_tier == 2 {
                        secondary_holobiont_count += 1;
                    } else if u.endo_tier >= 3 {
                        tertiary_holobiont_count += 1;
                    }
                    if !u.nucleomorph_sequence.is_empty() {
                        active_nucleomorph_count += 1;
                    }
                    if u.endo_tier >= 2 {
                        nm_len_sum += u.nucleomorph_sequence.len() as f64;
                    }
                    egt_frac_sum += u.egt_nuclear_fraction();
                    org_len_sum += u.endosymbiont_sequence.len() as f64;
                    tl_sum += u.guild.trophic_level();
                }
            }
        }

        let denom = w_sum.max(1e-6);
        let pop_n = self.population.len();
        let pop_f = pop_n.max(1) as f64;
        let mean_s_raw = s_raw_sum / denom;
        let mean_s_eff = s_eff_sum / denom;
        let mean_c = c_sum / denom;
        let leverage_li = li_sum / denom;
        let net_exergy_margin = margin_sum / pop_f;
        let mean_landauer = landauer_sum / pop_f;
        let mean_grid_e = self.grid.mean_exergy();
        let mean_niche_n = self.grid.mean_niche();
        let mean_byproduct_b = self.grid.mean_byproduct();
        let mean_trophic_level = if pop_n == 0 { 1.0 } else { tl_sum / pop_f };
        let mean_mro_pathways = if mro_bearers == 0 {
            0.0
        } else {
            mro_pathways_sum / (mro_bearers as f64)
        };
        let mean_cell_count = cell_count_sum / pop_f;
        let mean_soma_fraction = if multicell_bearers == 0 {
            0.0
        } else {
            soma_frac_sum / (multicell_bearers as f64)
        };
        let mean_somatic_cheater_load = if multicell_bearers == 0 {
            0.0
        } else {
            cheater_load_sum / (multicell_bearers as f64)
        };
        let mean_colony_size = colony_size_sum / pop_f;
        let mean_sterile_caste_fraction = if social_bearers == 0 {
            0.0
        } else {
            caste_frac_sum / (social_bearers as f64)
        };
        let mean_rebel_worker_load = if social_bearers == 0 {
            0.0
        } else {
            rebel_load_sum / (social_bearers as f64)
        };
        let mean_gongylidia_yield = if farming_colony_count == 0 {
            0.0
        } else {
            gongylidia_sum / (farming_colony_count as f64)
        };
        let mean_escovopsis_load = if farming_colony_count == 0 {
            0.0
        } else {
            escovopsis_sum / (farming_colony_count as f64)
        };
        let mean_cultural_repertoire_bits = if cultural_lineage_count == 0 {
            0.0
        } else {
            cultural_bits_sum / (cultural_lineage_count as f64)
        };
        let mean_neural_teaching_capacity = if cultural_lineage_count == 0 {
            0.0
        } else {
            neural_cap_sum / (cultural_lineage_count as f64)
        };
        let mean_exosomatic_archive_bits = if archival_lineage_count == 0 {
            0.0
        } else {
            archive_bits_sum / (archival_lineage_count as f64)
        };
        let mean_archive_fdr_load = if archival_lineage_count == 0 {
            0.0
        } else {
            archive_fdr_sum / (archival_lineage_count as f64)
        };
        let mean_verified_theorem_bits = if computational_symbiont_count == 0 {
            0.0
        } else {
            theorem_bits_sum / (computational_symbiont_count as f64)
        };
        let mean_tail_diversity = if computational_symbiont_count == 0 {
            1.0
        } else {
            tail_div_sum / (computational_symbiont_count as f64)
        };
        let mean_apicoplast_pathways = if apicoplast_parasite_count == 0 {
            0.0
        } else {
            apico_pathways_sum / (apicoplast_parasite_count as f64)
        };
        let mean_kleptoplast_integrity = if kleptoplastic_host_count == 0 {
            0.0
        } else {
            kp_integrity_sum / (kleptoplastic_host_count as f64)
        };
        let mean_klepto_hgt_loci = if klepto_hgt_bearers == 0 {
            0.0
        } else {
            klepto_hgt_sum / (klepto_hgt_bearers as f64)
        };
        let mosaic_donor_diversity = if mosaic_bearers == 0 {
            0.0
        } else {
            mosaic_donor_sum / (mosaic_bearers as f64)
        };
        let mean_egt_fraction = if holobiont_count == 0 {
            0.0
        } else {
            egt_frac_sum / (holobiont_count as f64)
        };
        let mean_organelle_len = if holobiont_count == 0 {
            0.0
        } else {
            org_len_sum / (holobiont_count as f64)
        };
        let serial_holo_n = secondary_holobiont_count + tertiary_holobiont_count;
        let mean_nucleomorph_len = if serial_holo_n == 0 {
            0.0
        } else {
            nm_len_sum / (serial_holo_n as f64)
        };
        let serial_endosymbiotic_fusions =
            phylo_sum.secondary_fusions + phylo_sum.tertiary_fusions;

        let food_web_entropy_bits = compute_food_web_entropy(
            producer_count,
            syntroph_count,
            predator_count,
            holobiont_count,
        );

        let theta_f = compute_critical_fidelity(mean_s_raw, self.controls.selective_superiority);
        let f_eff = self.controls.effective_fidelity();
        let f_obs = self.last_observed_fidelity.unwrap_or(f_eff);
        let blended_f = if self.smoothed_birth_rate >= 0.04 {
            0.65 * f_eff + 0.35 * f_obs
        } else {
            0.0
        };

        let cultural_teach_f = self.controls.symbolic_teaching_fidelity.clamp(0.0, 0.999);
        let cultural_heritable_mi = if mean_cultural_repertoire_bits > 0.05 && cultural_teach_f >= 0.50 {
            compute_heritable_mutual_info(cultural_teach_f, sigma_size, mean_cultural_repertoire_bits)
        } else {
            0.0
        };
        let archival_f = (0.96 + 0.038 * self.controls.peer_verification_rigor.clamp(0.0, 1.0)
            - 0.45 * mean_archive_fdr_load)
            .clamp(0.40, 0.999);
        let exosomatic_heritable_mi = if mean_exosomatic_archive_bits > 0.05 {
            compute_heritable_mutual_info(archival_f, sigma_size, mean_exosomatic_archive_bits)
        } else {
            0.0
        };
        let formal_kernel_f = (0.985 * mean_tail_diversity
            * self.controls.formal_verifier_strictness.clamp(0.45, 1.0))
        .clamp(0.40, 0.9995);
        let machine_theorem_mi = if mean_verified_theorem_bits > 0.05 {
            compute_heritable_mutual_info(formal_kernel_f, sigma_size, mean_verified_theorem_bits)
        } else {
            0.0
        };
        let heritable_mi = if self.smoothed_birth_rate >= 0.04 {
            compute_heritable_mutual_info(blended_f, sigma_size, mean_s_raw)
                + cultural_heritable_mi
                + exosomatic_heritable_mi
                + machine_theorem_mi
        } else {
            cultural_heritable_mi + exosomatic_heritable_mi + machine_theorem_mi
        };

        let thermo_hazard = compute_thermodynamic_hazard(
            self.controls.bath_exergy_flux,
            net_exergy_margin,
            mean_niche_n,
        );
        let meltdown_hazard = if self.smoothed_birth_rate >= 0.04
            && (1.0 - adami_rho) > 0.05
            && blended_f < theta_f
        {
            ((theta_f - blended_f) * 2.4).clamp(0.18, 0.58)
        } else {
            0.0
        };

        let (mut stability_st, mut volatility_nu, mut order_o, mut v_o) =
            self.telemetry_window.push_and_evaluate(
                leverage_li,
                self.controls.stability_epsilon,
                thermo_hazard,
                meltdown_hazard,
            );

        if self.controls.bath_exergy_flux.abs() <= 0.03 && net_exergy_margin.abs() <= 0.03 {
            stability_st = 1.0;
            volatility_nu = 0.0;
            order_o = leverage_li;
            v_o = 0.0;
        }

        let species_diversity_bonus = 1.0 + 0.16 * ((species_count.saturating_sub(1)) as f64);
        let trophic_web_bonus = 1.0 + 0.18 * food_web_entropy_bits;
        let holo_individuality_bonus = 1.0
            + 0.03 * (holobiont_count as f64) * (1.0 + 0.35 * mean_egt_fraction)
            + 0.04 * (serial_holo_n as f64)
            + 0.02 * (apicoplast_parasite_count as f64)
            + 0.02 * ((hydrogenosome_count + mitosome_count) as f64)
            + 0.045 * (clonal_multicell_count as f64)
            + 0.065 * (germ_soma_count as f64)
            + 0.055 * (eusocial_colony_count as f64)
            + 0.085 * (superorganism_count as f64)
            + 0.065 * (farming_colony_count as f64)
            + 0.095 * (obligate_agro_count as f64)
            + 0.075 * (cultural_lineage_count as f64)
            + 0.110 * (symbolic_culture_count as f64)
            + 0.090 * (archival_lineage_count as f64)
            + 0.135 * (verified_institution_count as f64)
            + 0.105 * (computational_symbiont_count as f64)
            + 0.155 * (formal_symbio_count as f64)
            + 0.015 * (kleptoplastic_host_count as f64);
        let ecosystem_order_pop = order_o
            * ((1 + pop_n) as f64).ln()
            * species_diversity_bonus
            * trophic_web_bonus
            * holo_individuality_bonus;
        let crowding_pressure = self.controls.interaction_strength
            * ((pop_n as f64) / (self.controls.carrying_capacity.max(1) as f64));

        let obs_var = (1.0 - adami_rho)
            .max(fst_index * 0.5)
            .max(food_web_entropy_bits * 0.15);
        self.smoothed_mut_diversity = 0.65 * self.smoothed_mut_diversity + 0.35 * obs_var;
        let effective_adami_for_var =
            (1.0 - self.smoothed_mut_diversity.max(obs_var)).clamp(0.0, 1.0);

        let (mut flags, coords, diagnostic, classified) = classify_emergent_state(
            self.controls.bath_exergy_flux,
            net_exergy_margin,
            v_o,
            self.smoothed_plastic_signal,
            self.smoothed_birth_rate,
            crowding_pressure,
            effective_adami_for_var,
            blended_f,
            theta_f,
            self.smoothed_assimilation,
            self.smoothed_recombination,
            &self.last_price,
            self.controls.artificial_selection,
        );
        flags.kleptoplasty =
            kleptoplastic_host_count > 0 || phylo_sum.klepto_fixation_events > 0;
        flags.klepto_tipping_fixation = phylo_sum.klepto_fixation_events > 0;
        flags.endosymbiotic_fusion =
            holobiont_count > 0 || phylo_sum.endosymbiotic_fusions > 0;
        flags.nuclear_egt = phylo_sum.egt_transfer_events > 0 || mean_egt_fraction >= 0.05;
        flags.serial_endosymbiosis = serial_holo_n > 0 || serial_endosymbiotic_fusions > 0;
        flags.nucleomorph_reduction = phylo_sum.nucleomorph_loss_events > 0
            || self
                .population
                .iter()
                .any(|u| u.nucleomorph_transferred_loci > 0);
        flags.apicoplast_reduction =
            apicoplast_parasite_count > 0 || phylo_sum.apicoplast_transition_events > 0;
        flags.apicoplast_loss = phylo_sum.apicoplast_loss_events > 0;
        flags.hydrogenosome_fermentation =
            hydrogenosome_count > 0 || phylo_sum.hydrogenosome_transition_events > 0;
        flags.mitosome_isc_retention =
            mitosome_count > 0 || phylo_sum.mitosome_transition_events > 0;
        flags.amiteukaryote_mro_loss =
            amiteukaryote_count > 0 || phylo_sum.mro_complete_loss_events > 0;
        flags.clonal_multicellularity =
            clonal_multicell_count > 0 || phylo_sum.clonal_multicell_transitions > 0;
        flags.germ_soma_division =
            germ_soma_count > 0 || phylo_sum.germ_soma_transitions > 0;
        flags.apoptotic_cheater_policing = self.total_apoptotic_culls > 0;
        flags.eusocial_colony =
            eusocial_colony_count > 0 || phylo_sum.eusocial_transitions > 0;
        flags.sterile_worker_castes =
            superorganism_count > 0 || phylo_sum.superorganism_transitions > 0;
        flags.worker_kin_policing = self.total_policed_eggs > 0;
        flags.agricultural_domestication =
            farming_colony_count > 0 || phylo_sum.agriculture_transitions > 0;
        flags.clonal_crop_monoculture =
            obligate_agro_count > 0 || phylo_sum.obligate_crop_transitions > 0;
        flags.defensive_microbiome_symbiosis = self.total_antimicrobial_shields > 0;
        flags.dual_gene_culture_inheritance = cultural_lineage_count > 0;
        flags.cumulative_cultural_ratchet = phylo_sum.cultural_ratchet_transitions > 0
            || symbolic_culture_count > 0;
        flags.symbolic_pedagogical_teaching = symbolic_culture_count > 0
            || phylo_sum.symbolic_teaching_transitions > 0;
        flags.exosomatic_memory_offloading = archival_lineage_count > 0;
        flags.written_archival_storage =
            phylo_sum.written_archive_transitions > 0 || verified_institution_count > 0;
        flags.institutional_peer_verification = verified_institution_count > 0
            || phylo_sum.peer_institution_transitions > 0;
        flags.computational_symbiont = computational_symbiont_count > 0;
        flags.formal_proof_verification =
            formal_symbio_count > 0 || phylo_sum.formal_verifier_transitions > 0;
        flags.model_collapse_safeguard = formal_symbio_count > 0
            && mean_tail_diversity >= 0.90
            && self.controls.empirical_grounding_floor >= 0.12;
        flags.oxygenic_goe_regulation =
            self.biosphere.goe_crossed || phylo_sum.goe_transitions > 0;
        flags.syntrophic_nutrient_closure = self.biosphere.nutrient_pool_np >= 0.65;
        flags.planetary_homeostatic_awakening = self.biosphere.homeostatic_awakening_crossed
            || phylo_sum.homeostatic_awakening_transitions > 0;
        let mean_oligomer_nt = self
            .population
            .iter()
            .map(|u| u.rna_oligomer_length_nt as f64)
            .sum::<f64>()
            / pop_f;
        let mean_raf_closure = self
            .population
            .iter()
            .map(|u| u.raf_closure_score)
            .sum::<f64>()
            / pop_f;
        flags.prebiotic_raf_closure = mean_raf_closure >= 0.50;
        flags.ribozyme_hypercycle_active = self
            .population
            .iter()
            .any(|u| u.abiogenesis_stage.has_ribozyme_heredity());
        flags.darwinian_threshold_luca = self
            .population
            .iter()
            .any(|u| u.abiogenesis_stage.crossed_darwinian_threshold());
        flags.homochiral_d_rna_l_aa =
            self.chiral_reactor.ribose_ee_d >= 0.85 && self.chiral_reactor.amino_ee_l >= 0.80;
        flags.stereoselective_aminoacylation = self
            .population
            .iter()
            .any(|u| u.chirality_code_stage.has_stereoselective_peptides());
        flags.error_minimized_universal_sgc = self
            .population
            .iter()
            .any(|u| u.chirality_code_stage.is_crystallized_sgc());
        flags.rnr_deoxyribose_synthesis = self
            .population
            .iter()
            .any(|u| u.genome_substrate_stage.has_deoxyribose_rnr());
        flags.thymine_udg_deamination_repair = self
            .population
            .iter()
            .any(|u| u.genome_substrate_stage.has_thymine_udg_repair());
        flags.virosphere_dna_coevolution = self
            .population
            .iter()
            .any(|u| u.genome_substrate_stage.is_coevolved_megabase_dna());
        flags.niche_construction = mean_niche_n >= 0.05;
        flags.cross_feeding =
            (syntroph_count + holobiont_count) > 0 && self.last_cross_feed_flux >= 0.02;
        flags.trophic_predation = predator_count > 0 && self.last_trophic_flux >= 0.015;
        self.classified_regime = classified;

        let phase = OrderPhase::from_velocity(v_o, 0.006);
        let metrics = ItmeMetrics {
            step: self.step_count,
            alphabet_size: sigma_size,
            structure_raw_bits: mean_s_raw,
            structure_s: mean_s_eff,
            adami_conservation_rho: adami_rho,
            complexity_c: mean_c,
            leverage_li,
            stability_st,
            volatility_nu,
            order_o,
            ecosystem_order_pop,
            order_velocity_vo: v_o,
            net_exergy_phi: net_exergy_margin,
            landauer_cost_q: mean_landauer,
            entropy_production_rate: entropy_prod_sum,
            mean_grid_exergy: mean_grid_e,
            mean_niche_capital: mean_niche_n,
            mean_byproduct_pool: mean_byproduct_b,
            species_count,
            fst_fixation_index: fst_index,
            producer_count,
            syntroph_count,
            predator_count,
            aerobic_mito_count,
            hydrogenosome_count,
            mitosome_count,
            amiteukaryote_count,
            mean_mro_pathways,
            hydrogenosome_transition_events: phylo_sum.hydrogenosome_transition_events,
            mitosome_transition_events: phylo_sum.mitosome_transition_events,
            mro_complete_loss_events: phylo_sum.mro_complete_loss_events,
            clonal_multicell_count,
            germ_soma_count,
            mean_cell_count,
            mean_soma_fraction,
            mean_somatic_cheater_load,
            clonal_multicell_transitions: phylo_sum.clonal_multicell_transitions,
            germ_soma_transitions: phylo_sum.germ_soma_transitions,
            apoptotic_cull_events: self.total_apoptotic_culls,
            neoplastic_collapse_events: self.total_neoplastic_collapses,
            eusocial_colony_count,
            superorganism_count,
            mean_colony_size,
            mean_sterile_caste_fraction,
            mean_rebel_worker_load,
            eusocial_transitions: phylo_sum.eusocial_transitions,
            superorganism_transitions: phylo_sum.superorganism_transitions,
            worker_policing_events: self.total_policed_eggs,
            colony_anarchy_collapses: self.total_colony_anarchy_collapses,
            farming_colony_count,
            obligate_agro_count,
            mean_gongylidia_yield,
            mean_escovopsis_load,
            agriculture_transitions: phylo_sum.agriculture_transitions,
            obligate_crop_transitions: phylo_sum.obligate_crop_transitions,
            antimicrobial_suppression_events: self.total_antimicrobial_shields,
            garden_blight_collapses: self.total_garden_blight_collapses,
            cultural_lineage_count,
            symbolic_culture_count,
            mean_cultural_repertoire_bits,
            mean_neural_teaching_capacity,
            cultural_ratchet_transitions: phylo_sum.cultural_ratchet_transitions,
            symbolic_teaching_transitions: phylo_sum.symbolic_teaching_transitions,
            cross_deme_recombinations: self.total_cross_deme_recombinations
                + phylo_sum.cross_deme_cultural_transfers,
            cultural_drift_collapses: self.total_cultural_drift_collapses,
            archival_lineage_count,
            verified_institution_count,
            mean_exosomatic_archive_bits,
            mean_archive_fdr_load,
            written_archive_transitions: phylo_sum.written_archive_transitions,
            peer_institution_transitions: phylo_sum.peer_institution_transitions,
            peer_retraction_audits: self.total_peer_retraction_audits,
            archive_corruption_collapses: self.total_archive_corruption_collapses,
            computational_symbiont_count,
            formal_symbio_count,
            mean_verified_theorem_bits,
            mean_tail_diversity,
            machine_symbiosis_transitions: phylo_sum.machine_symbiosis_transitions,
            formal_verifier_transitions: phylo_sum.formal_verifier_transitions,
            formal_proofs_verified: self.total_formal_proofs_verified,
            recursive_model_collapses: self.total_recursive_model_collapses,
            planetary_o2_level: self.biosphere.atmospheric_o2,
            planetary_co2_ppm: self.biosphere.atmospheric_co2_ppm,
            planetary_nutrient_pool: self.biosphere.nutrient_pool_np,
            planetary_exergy_budget: self.biosphere.planetary_exergy_budget,
            chemical_disequilibrium_bits: self.biosphere.chemical_disequilibrium_bits,
            boundary_transgression_load: self.biosphere.boundary_transgression_load,
            goe_transitions: phylo_sum.goe_transitions,
            homeostatic_awakening_transitions: phylo_sum.homeostatic_awakening_transitions,
            planetary_burnout_collapses: self.total_planetary_burnout_collapses,
            mean_rna_oligomer_nt: mean_oligomer_nt,
            raf_catalytic_closure: mean_raf_closure,
            eigen_max_length_nt: self.hadean_bath.eigen_max_length_nt,
            raf_protocell_transitions: phylo_sum.raf_protocell_transitions,
            darwinian_luca_transitions: phylo_sum.darwinian_threshold_transitions,
            eigen_error_catastrophes: self.total_eigen_error_catastrophes,
            mg_vesicle_collapses: self.total_mg_vesicle_collapses,
            ribose_ee_d: self.chiral_reactor.ribose_ee_d,
            amino_ee_l: self.chiral_reactor.amino_ee_l,
            sgc_error_minimization_mp: self.chiral_reactor.sgc_error_minimization_mp,
            encoded_amino_acids: self.chiral_reactor.encoded_amino_acids,
            homochiral_transitions: phylo_sum.homochiral_symmetry_transitions,
            sgc_crystallization_transitions: phylo_sum.sgc_crystallization_transitions,
            racemic_cross_inhibitions: self.total_racemic_cross_inhibitions,
            mistranslation_babel_collapses: self.total_mistranslation_babel_collapses,
            genome_capacity_kb: self.rnr_virosphere.stable_genome_capacity_kb,
            deamination_lesion_load: self.rnr_virosphere.deamination_lesion_load,
            dntp_pool_balance: self.rnr_virosphere.dntp_allosteric_balance,
            rnr_dna_transitions: phylo_sum.rnr_dna_transitions,
            udg_thymine_transitions: phylo_sum.thymine_udg_transitions,
            cytosine_deamination_meltdowns: self.total_cytosine_deamination_meltdowns,
            lytic_viral_collapses: self.total_lytic_viral_collapses,
            apicoplast_parasite_count,
            mean_apicoplast_pathways,
            apicoplast_transition_events: phylo_sum.apicoplast_transition_events,
            apicoplast_loss_events: phylo_sum.apicoplast_loss_events,
            kleptoplastic_host_count,
            mean_kleptoplast_integrity,
            mean_klepto_hgt_loci,
            klepto_hgt_events: phylo_sum.klepto_hgt_events,
            klepto_fixation_events: phylo_sum.klepto_fixation_events,
            mosaic_donor_diversity,
            holobiont_count,
            enslaved_holobiont_count,
            secondary_holobiont_count,
            tertiary_holobiont_count,
            active_nucleomorph_count,
            mean_nucleomorph_len,
            nucleomorph_loss_events: phylo_sum.nucleomorph_loss_events,
            serial_endosymbiotic_fusions,
            mean_egt_fraction,
            mean_organelle_len,
            egt_transfer_events: phylo_sum.egt_transfer_events,
            trophic_transfer_flux: self.last_trophic_flux,
            cross_feeding_flux: self.last_cross_feed_flux,
            mean_trophic_level,
            food_web_entropy_bits,
            tmrca_steps: phylo_sum.tmrca_steps,
            sackin_imbalance: phylo_sum.sackin_imbalance,
            speciation_events: phylo_sum.total_speciation_splits,
            guild_transition_events: phylo_sum.total_guild_transitions,
            endosymbiotic_fusions: phylo_sum.endosymbiotic_fusions,
            active_reticulations: phylo_sum.active_reticulations,
            fidelity_f: f_eff,
            observed_fidelity: f_obs,
            heritable_mutual_info_bits: heritable_mi,
            critical_fidelity_theta: theta_f,
            price: self.last_price,
            coords,
            population_size: pop_n,
            mean_fitness: self.mean_fitness(),
            flags,
            diagnostic,
            phase,
        };

        self.history.push(metrics);
        self.regime_history.push(classified);
        if self.history.len() > 780 {
            self.history.remove(0);
            self.regime_history.remove(0);
        }
    }

    pub fn latest_metrics(&self) -> &ItmeMetrics {
        self.history
            .last()
            .expect("metrics initialized on world creation")
    }

    pub fn current_scene(&self) -> StructureScene {
        let best = self
            .best_unit()
            .expect("population always maintains >= 1 unit");
        build_structure_scene(
            best,
            self.classified_regime,
            self.population.len(),
            self.phase_time,
        )
    }

    pub const TELEMETRY_CSV_HEADER: &'static str = "organism,regime_code,regime_name,step,S_raw_bits,S_eff_bits,Adami_rho_within,C_emergent,L_I,S_T_minus,Order_O,O_pop,v_O,Phi_net,Q_landauer,Grid_Exergy,Byproduct_Pool,Niche_Capital,Producers,Syntrophs,Predators,Aero_Mitos,Hydrogenosomes,Mitosomes,Amiteukaryotes,Mean_MRO_Paths,Hydro_Transitions,MitoS_Transitions,MRO_Losses,Apico_Parasites,Mean_Apico_Paths,Apico_Transitions,Apico_Losses,Klepto_Hosts,Mean_KP_Integ,Mean_Klepto_HGT,Klepto_HGT_Events,Klepto_Fixations,Mosaic_Donors,Holobionts,Enslaved_Holobionts,Secondary_Holobionts,Tertiary_Holobionts,Active_Nucleomorphs,Mean_NM_Len,NM_Loss_Events,Serial_Fusions,Mean_EGT_Frac,Mean_Organelle_Len,EGT_Events,Mean_TL,FoodWeb_H_bits,Trophic_Flux,CrossFeed_Flux,Species_Count,F_ST,TMRCA_steps,Sackin_Imbalance,Speciation_Splits,Guild_Transitions,Endosymbiotic_Fusions,Active_Reticulations,MI_bits,Price_Cov,Price_Trans,Pop,Emergent_Dc\n";

    pub fn csv_escape(raw: &str) -> String {
        let escaped = raw.replace('"', "\"\"");
        format!("\"{escaped}\"")
    }

    pub fn format_telemetry_csv_row(
        org: OrganismKind,
        regime: RegimeId,
        m: &ItmeMetrics,
    ) -> String {
        format!(
            "{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.5},{:.4},{:.4},{:.4},{:.4},{:.4},{},{},{},{},{},{},{},{:.2},{},{},{},{},{:.2},{},{},{},{:.3},{:.2},{},{},{:.2},{},{},{},{},{},{:.2},{},{},{:.4},{:.2},{},{:.3},{:.3},{:.4},{:.4},{},{:.4},{},{:.4},{},{},{},{},{:.3},{:.5},{:.5},{},{}\n",
            org.slug(),
            regime.code(),
            Self::csv_escape(regime.name()),
            m.step,
            m.structure_raw_bits,
            m.structure_s,
            m.adami_conservation_rho,
            m.complexity_c,
            m.leverage_li,
            m.stability_st,
            m.order_o,
            m.ecosystem_order_pop,
            m.order_velocity_vo,
            m.net_exergy_phi,
            m.landauer_cost_q,
            m.mean_grid_exergy,
            m.mean_byproduct_pool,
            m.mean_niche_capital,
            m.producer_count,
            m.syntroph_count,
            m.predator_count,
            m.aerobic_mito_count,
            m.hydrogenosome_count,
            m.mitosome_count,
            m.amiteukaryote_count,
            m.mean_mro_pathways,
            m.hydrogenosome_transition_events,
            m.mitosome_transition_events,
            m.mro_complete_loss_events,
            m.apicoplast_parasite_count,
            m.mean_apicoplast_pathways,
            m.apicoplast_transition_events,
            m.apicoplast_loss_events,
            m.kleptoplastic_host_count,
            m.mean_kleptoplast_integrity,
            m.mean_klepto_hgt_loci,
            m.klepto_hgt_events,
            m.klepto_fixation_events,
            m.mosaic_donor_diversity,
            m.holobiont_count,
            m.enslaved_holobiont_count,
            m.secondary_holobiont_count,
            m.tertiary_holobiont_count,
            m.active_nucleomorph_count,
            m.mean_nucleomorph_len,
            m.nucleomorph_loss_events,
            m.serial_endosymbiotic_fusions,
            m.mean_egt_fraction,
            m.mean_organelle_len,
            m.egt_transfer_events,
            m.mean_trophic_level,
            m.food_web_entropy_bits,
            m.trophic_transfer_flux,
            m.cross_feeding_flux,
            m.species_count,
            m.fst_fixation_index,
            m.tmrca_steps,
            m.sackin_imbalance,
            m.speciation_events,
            m.guild_transition_events,
            m.endosymbiotic_fusions,
            m.active_reticulations,
            m.heritable_mutual_info_bits,
            m.price.selection_covariance,
            m.price.transmission_bias,
            m.population_size,
            Self::csv_escape(&m.diagnostic.to_string())
        )
    }

    pub fn export_history_csv(&self) -> String {
        let mut out = String::from(Self::TELEMETRY_CSV_HEADER);
        for (idx, m) in self.history.iter().enumerate() {
            let reg = self
                .regime_history
                .get(idx)
                .copied()
                .unwrap_or(self.classified_regime);
            out.push_str(&Self::format_telemetry_csv_row(self.organism, reg, m));
        }
        out
    }

    pub fn export_grand_trajectory_csv(epochs: &[GrandTrajectoryEpochSummary]) -> String {
        let mut out = String::from(
            "epoch,age_ga,tick_end,population,channels,mutual_info_bits,leverage_li,o2,co2_ppm,exergy,landmark\n",
        );
        for s in epochs {
            out.push_str(&format!(
                "{},{:.3},{},{},{},{:.3},{:.3},{:.3},{:.1},{:.3},{}\n",
                s.epoch.code(),
                s.geologic_age_ga,
                s.tick_end,
                s.population_size,
                s.heritable_channels_active,
                s.heritable_mutual_info_bits,
                s.informational_leverage_li,
                s.planetary_o2,
                s.planetary_co2_ppm,
                s.planetary_exergy_budget,
                Self::csv_escape(&s.landmark_summary)
            ));
        }
        out
    }

    pub fn build_stage_atlas_svg(&self) -> String {
        let mut probe = Self::new(
            self.organism,
            RegimeId::ALL[0],
            SimulationMode::OrganicSandbox,
            12,
            20250813,
        );
        probe.coupled_biosphere_mode = self.coupled_biosphere_mode;
        let mut panels = Vec::with_capacity(RegimeId::ALL.len());
        for (r_idx, &regime) in RegimeId::ALL.iter().enumerate() {
            let preserve = r_idx > 2;
            probe.apply_regime_preset(regime, preserve);
            if !self.coupled_biosphere_mode {
                probe.controls = regime.domain_preset_controls(self.organism);
            }
            for _ in 0..10 {
                probe.step();
            }
            let m = probe.latest_metrics();
            let panel_title = format!("{} | {}", self.organism.short_title(), regime.code());
            let metric_line = if self.coupled_biosphere_mode {
                format!(
                    "O={:.2} O2={:.2} CO2={:.0}ppm",
                    m.order_o, m.planetary_o2_level, m.planetary_co2_ppm
                )
            } else {
                format!(
                    "O={:.2} Hydro={} MitoS={} Apico={}",
                    m.order_o,
                    m.hydrogenosome_count,
                    m.mitosome_count,
                    m.apicoplast_parasite_count
                )
            };
            panels.push((panel_title, metric_line, probe.current_scene()));
        }
        let title = if self.coupled_biosphere_mode {
            "ITME v23.0 Coupled 4.25 Ga Planetary Biosphere Stage Atlas".to_string()
        } else {
            format!("ITME v23.0 {} — 13-Regime Stage Atlas", self.organism.title())
        };
        render_stage_atlas_svg(&title, &panels)
    }
}
