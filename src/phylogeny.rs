//! ITME v12.0 (`v1.2.0`) Reticulate Phylogenetic Network: Anaerobic Mitochondrial Reduction
//! (`Mito->Hydro`, `Hydro->MitoS`, and *Monocercomonoides* `MRO->0`), Reductive Apicoplast Evolution
//! (`Holo->Apico` & *Cryptosporidium* `Apico->0`), Mixotrophic Kleptoplasty -> HGT Tipping-Point
//! Fixation (`KP->FIX`), Serial `1°/2°/3°` Endosymbiosis, Periplastidial Nucleomorph Loss (`NM->0`),
//! and Nuclear `EGT` (`src/phylogeny.rs`).
//!
//! Tracks explicit **Rooted Directed Acyclic Phylogenetic Networks** (`eNewick`, Cardona, Rosselló,
//! & Valiente, 2008) including:
//! - **Anaerobic MRO Reductive Transitions (`Mito->Hydro`, `Hydro->MitoS`, `MRO->0`)**: tracks
//!   stepwise conversion of `AerobicMitochondrion` holobionts under hypoxia into genome-less `H2`-producing
//!   `Hydrogenosome`s (*Trichomonas*), `ISC`-only `Mitosome`s (*Giardia* / microsporidia), and
//!   cytosolic-`SUF`-rescued `Amiteukaryote`s (*Monocercomonoides*; Stairs et al., 2015; Allen, 2015;
//!   Tsaousis, 2019; Karnkowska et al., 2016) while preserving ancestral `#Hk` endosymbiont donor arcs!
//! - **Reductive Apicoplast Transitions (`Holo->Apico` & `Apico->0`)**: secondary loss of photosynthesis
//!   retaining `MEP/IPP + SUF/Fe-S` pathways unless rescued by exogenous `IPP`.

use crate::model::{AbiogenesisStage, AgricultureState, ChiralityCodeStage, CulturalState, EusocialState, ExosomaticState, GenomeSubstrateStage, MachineSynthesisState, MroClass, MulticellularState, PlanetaryState, TrophicGuild};
use crate::spatial::BiomeZone;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ReticulationKind {
    HorizontalTransfer,
    EndosymbioticFusion,
    SecondaryEndosymbiosis,
    TertiaryEndosymbiosis,
    KleptoTippingFixation,
}

impl ReticulationKind {
    pub fn is_any_endosymbiosis(self) -> bool {
        matches!(
            self,
            Self::EndosymbioticFusion
                | Self::SecondaryEndosymbiosis
                | Self::TertiaryEndosymbiosis
                | Self::KleptoTippingFixation
        )
    }

    pub fn endo_tier(self) -> u8 {
        match self {
            Self::HorizontalTransfer => 0,
            Self::EndosymbioticFusion | Self::KleptoTippingFixation => 1,
            Self::SecondaryEndosymbiosis => 2,
            Self::TertiaryEndosymbiosis => 3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PhyloNode {
    pub taxon_id: u64,
    pub parent_id: Option<u64>,
    pub donor_id: Option<u64>,
    pub children: Vec<u64>,
    pub symbio_children: Vec<u64>,
    pub birth_step: usize,
    pub death_step: Option<usize>,
    pub extant: bool,
    pub species_id: usize,
    pub guild: TrophicGuild,
    pub biome: BiomeZone,
    pub is_speciation_split: bool,
    pub guild_transition: Option<(TrophicGuild, TrophicGuild)>,
    pub reticulation_kind: Option<ReticulationKind>,
    pub hybrid_id: Option<usize>,
    pub genome_share_gamma: f64,
    pub endo_tier: u8,
    pub nucleomorph_loci: usize,
    pub is_nucleomorph_loss: bool,
    pub klepto_hgt_loci: usize,
    pub has_kleptoplast: bool,
    pub is_klepto_tipping_point: bool,
    pub has_apicoplast: bool,
    pub apicoplast_pathways: u8,
    pub is_apicoplast_transition: bool,
    pub is_apicoplast_complete_loss: bool,
    pub mro_class: MroClass,
    pub mro_pathways: u8,
    pub mro_mtdna_loci: usize,
    /// Landmark flag marking transition to an anaerobic `H2`-producing `Hydrogenosome` (`mtDNA = 0`)!
    pub is_hydrogenosome_transition: bool,
    /// Landmark flag marking reduction to an `ISC`-only `Mitosome` (`ATP = 0, mtDNA = 0`)!
    pub is_mitosome_transition: bool,
    /// Landmark flag marking cytosolic-`SUF`-rescued complete loss of the MRO (`MRO -> 0`, *Monocercomonoides*)!
    pub is_mro_complete_loss: bool,
    pub multicell_state: MulticellularState,
    pub cell_count: u16,
    pub soma_fraction: f64,
    pub is_clonal_multicell_transition: bool,
    pub is_germ_soma_transition: bool,
    pub eusocial_state: EusocialState,
    pub colony_size: u16,
    pub sterile_caste_fraction: f64,
    pub is_eusocial_transition: bool,
    pub is_superorganism_transition: bool,
    pub agriculture_state: AgricultureState,
    pub cultivar_polyploidy: u8,
    pub gongylidia_yield: f64,
    pub is_agriculture_transition: bool,
    pub is_obligate_crop_transition: bool,
    pub cultural_state: CulturalState,
    pub cultural_ratchet_tier: u8,
    pub cultural_repertoire_bits: f64,
    pub is_cultural_ratchet_transition: bool,
    pub is_symbolic_teaching_transition: bool,
    pub exosomatic_state: ExosomaticState,
    pub exosomatic_archive_bits: f64,
    pub is_written_archive_transition: bool,
    pub is_peer_institution_transition: bool,
    pub machine_synthesis_state: MachineSynthesisState,
    pub verified_theorem_bits: f64,
    pub distribution_tail_diversity: f64,
    pub is_machine_symbiosis_transition: bool,
    pub is_formal_verifier_transition: bool,
    pub planetary_epoch: PlanetaryState,
    pub chemical_disequilibrium_bits: f64,
    pub is_goe_transition: bool,
    pub is_homeostatic_awakening_transition: bool,
    pub abiogenesis_stage: AbiogenesisStage,
    pub raf_closure_score: f64,
    pub is_raf_protocell_transition: bool,
    pub is_darwinian_threshold_transition: bool,
    pub chirality_code_stage: ChiralityCodeStage,
    pub genetic_code_mp: f64,
    pub is_homochiral_symmetry_transition: bool,
    pub is_sgc_crystallization_transition: bool,
    pub genome_substrate_stage: GenomeSubstrateStage,
    pub chromosome_capacity_kb: f64,
    pub is_rnr_dna_transition: bool,
    pub is_thymine_udg_transition: bool,
    pub egt_transferred_loci: usize,
    pub organelle_loci_remaining: usize,
    pub egt_nuclear_fraction: f64,
    pub is_egt_lockin: bool,
    pub seq_divergence: f64,
    pub leverage_li: f64,
}

#[derive(Copy, Clone, Debug, Default)]
pub struct PhyloSummary {
    pub mrca_id: u64,
    pub mrca_birth_step: usize,
    pub tmrca_steps: usize,
    pub active_node_count: usize,
    pub sackin_imbalance: f64,
    pub total_speciation_splits: usize,
    pub total_guild_transitions: usize,
    pub prod_to_syntroph_count: usize,
    pub to_predator_count: usize,
    pub endosymbiotic_fusions: usize,
    pub secondary_fusions: usize,
    pub tertiary_fusions: usize,
    pub klepto_fixation_events: usize,
    pub klepto_hgt_events: usize,
    pub apicoplast_transition_events: usize,
    pub apicoplast_loss_events: usize,
    pub hydrogenosome_transition_events: usize,
    pub mitosome_transition_events: usize,
    pub mro_complete_loss_events: usize,
    pub clonal_multicell_transitions: usize,
    pub germ_soma_transitions: usize,
    pub eusocial_transitions: usize,
    pub superorganism_transitions: usize,
    pub agriculture_transitions: usize,
    pub obligate_crop_transitions: usize,
    pub cultural_ratchet_transitions: usize,
    pub symbolic_teaching_transitions: usize,
    pub cross_deme_cultural_transfers: usize,
    pub written_archive_transitions: usize,
    pub peer_institution_transitions: usize,
    pub machine_symbiosis_transitions: usize,
    pub formal_verifier_transitions: usize,
    pub goe_transitions: usize,
    pub homeostatic_awakening_transitions: usize,
    pub raf_protocell_transitions: usize,
    pub darwinian_threshold_transitions: usize,
    pub homochiral_symmetry_transitions: usize,
    pub sgc_crystallization_transitions: usize,
    pub rnr_dna_transitions: usize,
    pub thymine_udg_transitions: usize,
    pub active_reticulations: usize,
    pub horizontal_transfers: usize,
    pub egt_transfer_events: usize,
    pub egt_lockin_events: usize,
    pub nucleomorph_loss_events: usize,
    pub mean_extant_egt_fraction: f64,
    pub mean_extant_nucleomorph_len: f64,
}

#[derive(Clone, Debug)]
pub struct PhylogramNodeVis {
    pub taxon_id: u64,
    pub x_norm: f64,
    pub y_norm: f64,
    pub extant: bool,
    pub species_id: usize,
    pub guild: TrophicGuild,
    pub biome: BiomeZone,
    pub is_speciation_split: bool,
    pub guild_transition: Option<(TrophicGuild, TrophicGuild)>,
    pub hybrid_id: Option<usize>,
    pub endo_tier: u8,
    pub nucleomorph_loci: usize,
    pub is_nucleomorph_loss: bool,
    pub has_kleptoplast: bool,
    pub is_klepto_tipping_point: bool,
    pub has_apicoplast: bool,
    pub apicoplast_pathways: u8,
    pub is_apicoplast_transition: bool,
    pub is_apicoplast_complete_loss: bool,
    pub mro_class: MroClass,
    pub mro_pathways: u8,
    pub is_hydrogenosome_transition: bool,
    pub is_mitosome_transition: bool,
    pub is_mro_complete_loss: bool,
    pub multicell_state: MulticellularState,
    pub cell_count: u16,
    pub is_clonal_multicell_transition: bool,
    pub is_germ_soma_transition: bool,
    pub eusocial_state: EusocialState,
    pub colony_size: u16,
    pub is_eusocial_transition: bool,
    pub is_superorganism_transition: bool,
    pub agriculture_state: AgricultureState,
    pub cultivar_polyploidy: u8,
    pub is_agriculture_transition: bool,
    pub is_obligate_crop_transition: bool,
    pub cultural_state: CulturalState,
    pub cultural_ratchet_tier: u8,
    pub is_cultural_ratchet_transition: bool,
    pub is_symbolic_teaching_transition: bool,
    pub exosomatic_state: ExosomaticState,
    pub is_written_archive_transition: bool,
    pub is_peer_institution_transition: bool,
    pub machine_synthesis_state: MachineSynthesisState,
    pub is_machine_symbiosis_transition: bool,
    pub is_formal_verifier_transition: bool,
    pub planetary_epoch: PlanetaryState,
    pub is_goe_transition: bool,
    pub is_homeostatic_awakening_transition: bool,
    pub abiogenesis_stage: AbiogenesisStage,
    pub is_raf_protocell_transition: bool,
    pub is_darwinian_threshold_transition: bool,
    pub chirality_code_stage: ChiralityCodeStage,
    pub is_homochiral_symmetry_transition: bool,
    pub is_sgc_crystallization_transition: bool,
    pub genome_substrate_stage: GenomeSubstrateStage,
    pub is_rnr_dna_transition: bool,
    pub is_thymine_udg_transition: bool,
    pub egt_nuclear_fraction: f64,
    pub is_egt_lockin: bool,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct PhylogramBranchVis {
    pub parent_x: f64,
    pub parent_y: f64,
    pub child_x: f64,
    pub child_y: f64,
    pub child_guild: TrophicGuild,
    pub child_species: usize,
    pub child_endo_tier: u8,
    pub child_has_apicoplast: bool,
    pub child_mro_class: MroClass,
    pub is_speciation_split: bool,
    pub guild_transition: Option<(TrophicGuild, TrophicGuild)>,
    pub child_egt_fraction: f64,
    pub is_egt_lockin: bool,
    pub is_nucleomorph_loss: bool,
    pub is_klepto_tipping_point: bool,
    pub is_apicoplast_transition: bool,
    pub is_apicoplast_complete_loss: bool,
    pub is_hydrogenosome_transition: bool,
    pub is_mitosome_transition: bool,
    pub is_mro_complete_loss: bool,
    pub child_multicell_state: MulticellularState,
    pub is_clonal_multicell_transition: bool,
    pub is_germ_soma_transition: bool,
    pub child_eusocial_state: EusocialState,
    pub is_eusocial_transition: bool,
    pub is_superorganism_transition: bool,
    pub child_agriculture_state: AgricultureState,
    pub is_agriculture_transition: bool,
    pub is_obligate_crop_transition: bool,
    pub child_cultural_state: CulturalState,
    pub is_cultural_ratchet_transition: bool,
    pub is_symbolic_teaching_transition: bool,
    pub child_exosomatic_state: ExosomaticState,
    pub is_written_archive_transition: bool,
    pub is_peer_institution_transition: bool,
    pub child_machine_synthesis_state: MachineSynthesisState,
    pub is_machine_symbiosis_transition: bool,
    pub is_formal_verifier_transition: bool,
    pub child_planetary_epoch: PlanetaryState,
    pub is_goe_transition: bool,
    pub is_homeostatic_awakening_transition: bool,
    pub child_abiogenesis_stage: AbiogenesisStage,
    pub is_raf_protocell_transition: bool,
    pub is_darwinian_threshold_transition: bool,
    pub child_chirality_code_stage: ChiralityCodeStage,
    pub is_homochiral_symmetry_transition: bool,
    pub is_sgc_crystallization_transition: bool,
    pub child_genome_substrate_stage: GenomeSubstrateStage,
    pub is_rnr_dna_transition: bool,
    pub is_thymine_udg_transition: bool,
}

#[derive(Clone, Debug)]
pub struct PhylogramReticulationVis {
    pub from_x: f64,
    pub from_y: f64,
    pub to_x: f64,
    pub to_y: f64,
    pub is_endosymbiosis: bool,
    pub is_klepto_tipping: bool,
    pub endo_tier: u8,
    pub hybrid_id: Option<usize>,
    pub gamma: f64,
}

#[derive(Clone, Debug, Default)]
pub struct PhylogramLayout {
    pub min_step: usize,
    pub max_step: usize,
    pub mrca_id: u64,
    pub nodes: Vec<PhylogramNodeVis>,
    pub branches: Vec<PhylogramBranchVis>,
    pub reticulations: Vec<PhylogramReticulationVis>,
}

#[derive(Copy, Clone, Debug)]
pub struct ExtantPhyloState {
    pub taxon_id: u64,
    pub species_id: usize,
    pub guild: TrophicGuild,
    pub biome: BiomeZone,
    pub leverage_li: f64,
    pub egt_loci: usize,
    pub organelle_rem: usize,
    pub endo_tier: u8,
    pub nucleomorph_loci: usize,
    pub nucleomorph_eliminated: bool,
    pub klepto_hgt_loci: usize,
    pub has_kleptoplast: bool,
    pub has_apicoplast: bool,
    pub apicoplast_pathways: u8,
    pub apicoplast_eliminated: bool,
    pub mro_class: MroClass,
    pub mro_pathways: u8,
    pub mro_mtdna_loci: usize,
    pub multicell_state: MulticellularState,
    pub cell_count: u16,
    pub soma_fraction: f64,
    pub eusocial_state: EusocialState,
    pub colony_size: u16,
    pub sterile_caste_fraction: f64,
    pub agriculture_state: AgricultureState,
    pub cultivar_polyploidy: u8,
    pub gongylidia_yield: f64,
    pub cultural_state: CulturalState,
    pub cultural_ratchet_tier: u8,
    pub cultural_repertoire_bits: f64,
    pub exosomatic_state: ExosomaticState,
    pub exosomatic_archive_bits: f64,
    pub machine_synthesis_state: MachineSynthesisState,
    pub verified_theorem_bits: f64,
    pub distribution_tail_diversity: f64,
    pub planetary_epoch: PlanetaryState,
    pub chemical_disequilibrium_bits: f64,
    pub abiogenesis_stage: AbiogenesisStage,
    pub raf_closure_score: f64,
    pub chirality_code_stage: ChiralityCodeStage,
    pub genetic_code_mp: f64,
    pub genome_substrate_stage: GenomeSubstrateStage,
    pub chromosome_capacity_kb: f64,
}

#[derive(Clone, Debug)]
pub struct PhylogenyTracker {
    pub nodes: BTreeMap<u64, PhyloNode>,
    pub root_id: u64,
    pub total_speciation_splits: usize,
    pub total_guild_transitions: usize,
    pub prod_to_syntroph_count: usize,
    pub to_predator_count: usize,
    pub endosymbiotic_fusions: usize,
    pub secondary_fusions: usize,
    pub tertiary_fusions: usize,
    pub klepto_fixation_events: usize,
    pub klepto_hgt_events: usize,
    pub apicoplast_transition_events: usize,
    pub apicoplast_loss_events: usize,
    pub hydrogenosome_transition_events: usize,
    pub mitosome_transition_events: usize,
    pub mro_complete_loss_events: usize,
    pub clonal_multicell_transitions: usize,
    pub germ_soma_transitions: usize,
    pub eusocial_transitions: usize,
    pub superorganism_transitions: usize,
    pub agriculture_transitions: usize,
    pub obligate_crop_transitions: usize,
    pub cultural_ratchet_transitions: usize,
    pub symbolic_teaching_transitions: usize,
    pub cross_deme_cultural_transfers: usize,
    pub written_archive_transitions: usize,
    pub peer_institution_transitions: usize,
    pub machine_symbiosis_transitions: usize,
    pub formal_verifier_transitions: usize,
    pub goe_transitions: usize,
    pub homeostatic_awakening_transitions: usize,
    pub raf_protocell_transitions: usize,
    pub darwinian_threshold_transitions: usize,
    pub homochiral_symmetry_transitions: usize,
    pub sgc_crystallization_transitions: usize,
    pub rnr_dna_transitions: usize,
    pub thymine_udg_transitions: usize,
    pub horizontal_transfers: usize,
    pub egt_transfer_events: usize,
    pub egt_lockin_events: usize,
    pub nucleomorph_loss_events: usize,
    pub next_hybrid_id: usize,
}

impl Default for PhylogenyTracker {
    fn default() -> Self {
        Self::new()
    }
}

impl PhylogenyTracker {
    pub fn new() -> Self {
        let mut nodes = BTreeMap::new();
        nodes.insert(
            0,
            PhyloNode {
                taxon_id: 0,
                parent_id: None,
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: 0,
                death_step: Some(0),
                extant: false,
                species_id: 0,
                guild: TrophicGuild::Producer,
                biome: BiomeZone::Alpha,
                is_speciation_split: false,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: 0,
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::None,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: 0,
                organelle_loci_remaining: 0,
                egt_nuclear_fraction: 0.0,
                is_egt_lockin: false,
                seq_divergence: 0.0,
                leverage_li: 1.0,
            },
        );
        Self {
            nodes,
            root_id: 0,
            total_speciation_splits: 0,
            total_guild_transitions: 0,
            prod_to_syntroph_count: 0,
            to_predator_count: 0,
            endosymbiotic_fusions: 0,
            secondary_fusions: 0,
            tertiary_fusions: 0,
            klepto_fixation_events: 0,
            klepto_hgt_events: 0,
            apicoplast_transition_events: 0,
            apicoplast_loss_events: 0,
            hydrogenosome_transition_events: 0,
            mitosome_transition_events: 0,
            mro_complete_loss_events: 0,
            clonal_multicell_transitions: 0,
            germ_soma_transitions: 0,
            eusocial_transitions: 0,
            superorganism_transitions: 0,
            agriculture_transitions: 0,
            obligate_crop_transitions: 0,
            cultural_ratchet_transitions: 0,
            symbolic_teaching_transitions: 0,
            cross_deme_cultural_transfers: 0,
            written_archive_transitions: 0,
            peer_institution_transitions: 0,
            machine_symbiosis_transitions: 0,
            formal_verifier_transitions: 0,
            goe_transitions: 0,
            homeostatic_awakening_transitions: 0,
            raf_protocell_transitions: 0,
            darwinian_threshold_transitions: 0,
            homochiral_symmetry_transitions: 0,
            sgc_crystallization_transitions: 0,
            rnr_dna_transitions: 0,
            thymine_udg_transitions: 0,
            horizontal_transfers: 0,
            egt_transfer_events: 0,
            egt_lockin_events: 0,
            nucleomorph_loss_events: 0,
            next_hybrid_id: 1,
        }
    }

    pub fn reset_with_founders(
        &mut self,
        step: usize,
        founders: &[(u64, usize, TrophicGuild, BiomeZone, f64)],
    ) {
        self.nodes.clear();
        self.root_id = 0;
        self.total_speciation_splits = 0;
        self.total_guild_transitions = 0;
        self.prod_to_syntroph_count = 0;
        self.to_predator_count = 0;
        self.endosymbiotic_fusions = 0;
        self.secondary_fusions = 0;
        self.tertiary_fusions = 0;
        self.klepto_fixation_events = 0;
        self.klepto_hgt_events = 0;
        self.apicoplast_transition_events = 0;
        self.apicoplast_loss_events = 0;
        self.hydrogenosome_transition_events = 0;
        self.mitosome_transition_events = 0;
        self.mro_complete_loss_events = 0;
        self.clonal_multicell_transitions = 0;
        self.germ_soma_transitions = 0;
        self.eusocial_transitions = 0;
        self.superorganism_transitions = 0;
        self.agriculture_transitions = 0;
        self.obligate_crop_transitions = 0;
        self.cultural_ratchet_transitions = 0;
        self.symbolic_teaching_transitions = 0;
        self.cross_deme_cultural_transfers = 0;
        self.written_archive_transitions = 0;
        self.peer_institution_transitions = 0;
        self.machine_symbiosis_transitions = 0;
        self.formal_verifier_transitions = 0;
        self.goe_transitions = 0;
        self.homeostatic_awakening_transitions = 0;
        self.raf_protocell_transitions = 0;
        self.darwinian_threshold_transitions = 0;
        self.homochiral_symmetry_transitions = 0;
        self.sgc_crystallization_transitions = 0;
        self.rnr_dna_transitions = 0;
        self.thymine_udg_transitions = 0;
        self.horizontal_transfers = 0;
        self.egt_transfer_events = 0;
        self.egt_lockin_events = 0;
        self.nucleomorph_loss_events = 0;
        self.next_hybrid_id = 1;

        let mut root_children = Vec::with_capacity(founders.len());
        for &(uid, _sp, _g, _b, _l) in founders {
            root_children.push(uid);
        }

        self.nodes.insert(
            0,
            PhyloNode {
                taxon_id: 0,
                parent_id: None,
                donor_id: None,
                children: root_children,
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: Some(step),
                extant: false,
                species_id: 0,
                guild: TrophicGuild::Producer,
                biome: BiomeZone::Alpha,
                is_speciation_split: false,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: 0,
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::None,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: 0,
                organelle_loci_remaining: 0,
                egt_nuclear_fraction: 0.0,
                is_egt_lockin: false,
                seq_divergence: 0.0,
                leverage_li: 1.0,
            },
        );

        for &(uid, sp_id, guild, biome, lev) in founders {
            let shift = if guild != TrophicGuild::Producer {
                self.record_guild_shift_counters(TrophicGuild::Producer, guild);
                Some((TrophicGuild::Producer, guild))
            } else {
                None
            };
            let is_holo = guild == TrophicGuild::Holobiont;
            self.nodes.insert(
                uid,
                PhyloNode {
                    taxon_id: uid,
                    parent_id: Some(0),
                    donor_id: None,
                    children: Vec::new(),
                    symbio_children: Vec::new(),
                    birth_step: step,
                    death_step: None,
                    extant: true,
                    species_id: sp_id,
                    guild,
                    biome,
                    is_speciation_split: sp_id > 0,
                    guild_transition: shift,
                    reticulation_kind: None,
                    hybrid_id: None,
                    genome_share_gamma: 0.0,
                    endo_tier: u8::from(is_holo),
                    nucleomorph_loci: 0,
                    is_nucleomorph_loss: false,
                    klepto_hgt_loci: 0,
                    has_kleptoplast: false,
                    is_klepto_tipping_point: false,
                    has_apicoplast: false,
                    apicoplast_pathways: if is_holo { 4 } else { 0 },
                    is_apicoplast_transition: false,
                    is_apicoplast_complete_loss: false,
                    mro_class: if is_holo {
                        MroClass::AerobicMitochondrion
                    } else {
                        MroClass::None
                    },
                    mro_pathways: if is_holo { 4 } else { 0 },
                    mro_mtdna_loci: if is_holo { 3 } else { 0 },
                    is_hydrogenosome_transition: false,
                    is_mitosome_transition: false,
                    is_mro_complete_loss: false,
                    multicell_state: MulticellularState::Unicellular,
                    cell_count: 1,
                    soma_fraction: 0.0,
                    is_clonal_multicell_transition: false,
                    is_germ_soma_transition: false,
                    eusocial_state: EusocialState::Solitary,
                    colony_size: 1,
                    sterile_caste_fraction: 0.0,
                    is_eusocial_transition: false,
                    is_superorganism_transition: false,
                    agriculture_state: AgricultureState::HunterForager,
                    cultivar_polyploidy: 1,
                    gongylidia_yield: 0.0,
                    is_agriculture_transition: false,
                    is_obligate_crop_transition: false,
                    cultural_state: CulturalState::AsocialTrialError,
                    cultural_ratchet_tier: 0,
                    cultural_repertoire_bits: 0.0,
                    is_cultural_ratchet_transition: false,
                    is_symbolic_teaching_transition: false,
                    exosomatic_state: ExosomaticState::EndosomaticBound,
                    exosomatic_archive_bits: 0.0,
                    is_written_archive_transition: false,
                    is_peer_institution_transition: false,
                    machine_synthesis_state: MachineSynthesisState::ManualArchival,
                    verified_theorem_bits: 0.0,
                    distribution_tail_diversity: 1.0,
                    is_machine_symbiosis_transition: false,
                    is_formal_verifier_transition: false,
                    planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                    chemical_disequilibrium_bits: 14.5,
                    is_goe_transition: false,
                    is_homeostatic_awakening_transition: false,
                    abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                    raf_closure_score: 0.96,
                    is_raf_protocell_transition: false,
                    is_darwinian_threshold_transition: false,
                    chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                    genetic_code_mp: 0.79,
                    is_homochiral_symmetry_transition: false,
                    is_sgc_crystallization_transition: false,
                    genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                    chromosome_capacity_kb: 1850.0,
                    is_rnr_dna_transition: false,
                    is_thymine_udg_transition: false,
                    egt_transferred_loci: 0,
                    organelle_loci_remaining: 0,
                    egt_nuclear_fraction: 0.0,
                    is_egt_lockin: false,
                    seq_divergence: 0.0,
                    leverage_li: lev,
                },
            );
        }
    }

    fn record_guild_shift_counters(&mut self, from: TrophicGuild, to: TrophicGuild) {
        if from == to {
            return;
        }
        self.total_guild_transitions += 1;
        if from == TrophicGuild::Producer && to == TrophicGuild::Syntroph {
            self.prod_to_syntroph_count += 1;
        }
        if to == TrophicGuild::PredatorParasite {
            self.to_predator_count += 1;
        }
    }

    /// Registers an **Anaerobic Mitochondrial Reduction Transition** (`Mito -> Hydrogenosome`,
    /// `Hydrogenosome -> Mitosome`, or cytosolic-`SUF`-rescued `Mitosome -> Amiteukaryote MRO->0`)
    /// while preserving ancestral `#Hk` endosymbiont donor reticulations in the phylogenetic DAG!
    #[allow(clippy::too_many_arguments)]
    pub fn record_mro_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let (prev_guild, prev_mro, prev_apico, prev_apico_paths) = self
            .nodes
            .get(&host_id)
            .map(|n| (n.guild, n.mro_class, n.has_apicoplast, n.apicoplast_pathways))
            .unwrap_or((guild, MroClass::AerobicMitochondrion, false, 0));

        let is_hydro_tr =
            mro_class == MroClass::Hydrogenosome && prev_mro != MroClass::Hydrogenosome;
        let is_mitos_tr = mro_class == MroClass::Mitosome && prev_mro != MroClass::Mitosome;
        let is_mro_zero =
            mro_class == MroClass::Amiteukaryote && prev_mro != MroClass::Amiteukaryote;

        if is_hydro_tr {
            self.hydrogenosome_transition_events += 1;
        }
        if is_mitos_tr {
            self.mitosome_transition_events += 1;
        }
        if is_mro_zero {
            self.mro_complete_loss_events += 1;
        }
        self.total_speciation_splits += 1;
        let shift = if prev_guild != guild {
            self.record_guild_shift_counters(prev_guild, guild);
            Some((prev_guild, guild))
        } else {
            None
        };

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: shift,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: prev_apico,
                apicoplast_pathways: prev_apico_paths,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: is_hydro_tr,
                is_mitosome_transition: is_mitos_tr,
                is_mro_complete_loss: is_mro_zero,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50 || mro_mtdna_loci == 0,
                seq_divergence: 0.18,
                leverage_li,
            },
        );
    }

    /// Registers a **Clonal Multicellularity (`Holo -> Clonal`)** or **Soma-Germline Division (`Clonal -> GermSoma`)**
    /// evolutionary transition in individuality while preserving `#Hk` endosymbiotic donor reticulations!
    #[allow(clippy::too_many_arguments)]
    pub fn record_multicellular_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_mc = self
            .nodes
            .get(&host_id)
            .map(|n| n.multicell_state)
            .unwrap_or(MulticellularState::Unicellular);

        let is_clonal_tr = !prev_mc.is_multicellular() && multicell_state.is_multicellular();
        let is_gs_tr = prev_mc != MulticellularState::GermSomaMulticellular
            && multicell_state == MulticellularState::GermSomaMulticellular;
        if is_clonal_tr {
            self.clonal_multicell_transitions += 1;
        }
        if is_gs_tr {
            self.germ_soma_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: is_clonal_tr,
                is_germ_soma_transition: is_gs_tr,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers a **Major Evolutionary Transition to Eusociality (`GermSoma -> Eusocial`)** or
    /// **Point-of-No-Return Superorganismality (`Eusocial -> Superorg`)** while preserving `#Hk`
    /// endosymbiotic donor arcs and multicellular transition landmarks!
    #[allow(clippy::too_many_arguments)]
    pub fn record_superorganism_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_eu = self
            .nodes
            .get(&host_id)
            .map(|n| n.eusocial_state)
            .unwrap_or(EusocialState::Solitary);

        let is_eusoc_tr = !prev_eu.is_eusocial() && eusocial_state.is_eusocial();
        let is_super_tr = prev_eu != EusocialState::ObligateSuperorganism
            && eusocial_state == EusocialState::ObligateSuperorganism;
        if is_eusoc_tr {
            self.eusocial_transitions += 1;
        }
        if is_super_tr {
            self.superorganism_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: is_eusoc_tr,
                is_superorganism_transition: is_super_tr,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers an **Agricultural Ectosymbiosis Origin (`Superorg -> LowerAg`)** or **Obligate Clonal
    /// Crop / Leaf-Cutter Domestication (`LowerAg -> GongylCrop / LeafCutter`)** while preserving all
    /// upstream `#Hk` endosymbiotic, multicellular, and superorganismal landmarks!
    #[allow(clippy::too_many_arguments)]
    pub fn record_agriculture_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        agriculture_state: AgricultureState,
        cultivar_polyploidy: u8,
        gongylidia_yield: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_ag = self
            .nodes
            .get(&host_id)
            .map(|n| n.agriculture_state)
            .unwrap_or(AgricultureState::HunterForager);

        let is_ag_tr = !prev_ag.is_farming() && agriculture_state.is_farming();
        let is_crop_tr = !prev_ag.is_obligate_monoculture()
            && agriculture_state.is_obligate_monoculture();
        if is_ag_tr {
            self.agriculture_transitions += 1;
        }
        if is_crop_tr {
            self.obligate_crop_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state,
                cultivar_polyploidy: cultivar_polyploidy.max(1),
                gongylidia_yield,
                is_agriculture_transition: is_ag_tr,
                is_obligate_crop_transition: is_crop_tr,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_cultural_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        teacher_donor_id: Option<u64>,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        agriculture_state: AgricultureState,
        cultivar_polyploidy: u8,
        gongylidia_yield: f64,
        cultural_state: CulturalState,
        cultural_ratchet_tier: u8,
        cultural_repertoire_bits: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_cult = self
            .nodes
            .get(&host_id)
            .map(|n| n.cultural_state)
            .unwrap_or(CulturalState::AsocialTrialError);

        let is_ratchet_tr =
            !prev_cult.is_cumulative_ratchet() && cultural_state.is_cumulative_ratchet();
        let is_sym_tr = prev_cult != CulturalState::SymbolicCumulativeCulture
            && cultural_state == CulturalState::SymbolicCumulativeCulture;
        if is_ratchet_tr {
            self.cultural_ratchet_transitions += 1;
        }
        if is_sym_tr {
            self.symbolic_teaching_transitions += 1;
        }
        self.total_speciation_splits += 1;

        let valid_teacher = teacher_donor_id.filter(|&tid| {
            tid != host_id && tid != new_taxon_id && self.nodes.contains_key(&tid)
        });
        if let Some(tid) = valid_teacher {
            self.cross_deme_cultural_transfers += 1;
            self.horizontal_transfers += 1;
            if let Some(d_node) = self.nodes.get_mut(&tid) {
                if !d_node.symbio_children.contains(&new_taxon_id) {
                    d_node.symbio_children.push(new_taxon_id);
                }
            }
        }

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: valid_teacher,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: valid_teacher.map(|_| ReticulationKind::HorizontalTransfer),
                hybrid_id: None,
                genome_share_gamma: if valid_teacher.is_some() { 0.25 } else { 0.0 },
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state,
                cultivar_polyploidy: cultivar_polyploidy.max(1),
                gongylidia_yield,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state,
                cultural_ratchet_tier,
                cultural_repertoire_bits,
                is_cultural_ratchet_transition: is_ratchet_tr,
                is_symbolic_teaching_transition: is_sym_tr,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers a **Written Archival Transition (`Symbolic -> WrittenArch`)** or **Peer-Verified
    /// Epistemic Institution Transition (`WrittenArch -> VerifiedInst`)**, plus optional asynchronous
    /// archival citation transfer (`archive_donor_id`) across generations/demes!
    #[allow(clippy::too_many_arguments)]
    pub fn record_exosomatic_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        archive_donor_id: Option<u64>,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        agriculture_state: AgricultureState,
        cultivar_polyploidy: u8,
        gongylidia_yield: f64,
        cultural_state: CulturalState,
        cultural_ratchet_tier: u8,
        cultural_repertoire_bits: f64,
        exosomatic_state: ExosomaticState,
        exosomatic_archive_bits: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_exo = self
            .nodes
            .get(&host_id)
            .map(|n| n.exosomatic_state)
            .unwrap_or(ExosomaticState::EndosomaticBound);

        let is_arch_tr = !prev_exo.has_written_archive() && exosomatic_state.has_written_archive();
        let is_inst_tr = prev_exo != ExosomaticState::PeerVerifiedInstitution
            && exosomatic_state == ExosomaticState::PeerVerifiedInstitution;
        if is_arch_tr {
            self.written_archive_transitions += 1;
        }
        if is_inst_tr {
            self.peer_institution_transitions += 1;
        }
        self.total_speciation_splits += 1;

        let valid_archive_donor = archive_donor_id.filter(|&aid| {
            aid != host_id && aid != new_taxon_id && self.nodes.contains_key(&aid)
        });
        if let Some(aid) = valid_archive_donor {
            self.horizontal_transfers += 1;
            if let Some(d_node) = self.nodes.get_mut(&aid) {
                if !d_node.symbio_children.contains(&new_taxon_id) {
                    d_node.symbio_children.push(new_taxon_id);
                }
            }
        }

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: valid_archive_donor,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: valid_archive_donor.map(|_| ReticulationKind::HorizontalTransfer),
                hybrid_id: None,
                genome_share_gamma: if valid_archive_donor.is_some() { 0.30 } else { 0.0 },
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state,
                cultivar_polyploidy: cultivar_polyploidy.max(1),
                gongylidia_yield,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state,
                cultural_ratchet_tier,
                cultural_repertoire_bits,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state,
                exosomatic_archive_bits,
                is_written_archive_transition: is_arch_tr,
                is_peer_institution_transition: is_inst_tr,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers a **Computational Symbiont Coevolution Transition (`Arch -> AlgoMed`)** or
    /// **Machine-Verified Neuro-Symbolic Synthesis Transition (`AlgoMed -> FormalSymbio`)**,
    /// optionally recording a `#LGT` computational-symbiont fusion arc from `symbiont_ai_donor_id`!
    #[allow(clippy::too_many_arguments)]
    pub fn record_machine_synthesis_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        symbiont_ai_donor_id: Option<u64>,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        agriculture_state: AgricultureState,
        cultivar_polyploidy: u8,
        gongylidia_yield: f64,
        cultural_state: CulturalState,
        cultural_ratchet_tier: u8,
        cultural_repertoire_bits: f64,
        exosomatic_state: ExosomaticState,
        exosomatic_archive_bits: f64,
        machine_synthesis_state: MachineSynthesisState,
        verified_theorem_bits: f64,
        distribution_tail_diversity: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_ms = self
            .nodes
            .get(&host_id)
            .map(|n| n.machine_synthesis_state)
            .unwrap_or(MachineSynthesisState::ManualArchival);

        let is_ms_tr = !prev_ms.has_computational_symbiont()
            && machine_synthesis_state.has_computational_symbiont();
        let is_formal_tr =
            !prev_ms.is_formally_verified() && machine_synthesis_state.is_formally_verified();
        if is_ms_tr {
            self.machine_symbiosis_transitions += 1;
        }
        if is_formal_tr {
            self.formal_verifier_transitions += 1;
        }
        self.total_speciation_splits += 1;

        let valid_ai_donor = symbiont_ai_donor_id.filter(|&aid| {
            aid != host_id && aid != new_taxon_id && self.nodes.contains_key(&aid)
        });
        if let Some(aid) = valid_ai_donor {
            self.horizontal_transfers += 1;
            if let Some(d_node) = self.nodes.get_mut(&aid) {
                if !d_node.symbio_children.contains(&new_taxon_id) {
                    d_node.symbio_children.push(new_taxon_id);
                }
            }
        }

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: valid_ai_donor,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: valid_ai_donor.map(|_| ReticulationKind::HorizontalTransfer),
                hybrid_id: None,
                genome_share_gamma: if valid_ai_donor.is_some() { 0.34 } else { 0.0 },
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state,
                cultivar_polyploidy: cultivar_polyploidy.max(1),
                gongylidia_yield,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state,
                cultural_ratchet_tier,
                cultural_repertoire_bits,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state,
                exosomatic_archive_bits,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state,
                verified_theorem_bits,
                distribution_tail_diversity,
                is_machine_symbiosis_transition: is_ms_tr,
                is_formal_verifier_transition: is_formal_tr,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers a **Great Oxidation Event (`Archean -> GOE/Phanerozoic`)** or **Gaia-2.0
    /// Homeostatic Awakening (`Techno -> Gaia2Homeostasis`)** landmark in the reticulate phylogeny!
    #[allow(clippy::too_many_arguments)]
    pub fn record_planetary_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        endo_tier: u8,
        mro_class: MroClass,
        mro_pathways: u8,
        mro_mtdna_loci: usize,
        multicell_state: MulticellularState,
        cell_count: u16,
        soma_fraction: f64,
        eusocial_state: EusocialState,
        colony_size: u16,
        sterile_caste_fraction: f64,
        agriculture_state: AgricultureState,
        cultivar_polyploidy: u8,
        gongylidia_yield: f64,
        cultural_state: CulturalState,
        cultural_ratchet_tier: u8,
        cultural_repertoire_bits: f64,
        exosomatic_state: ExosomaticState,
        exosomatic_archive_bits: f64,
        machine_synthesis_state: MachineSynthesisState,
        verified_theorem_bits: f64,
        distribution_tail_diversity: f64,
        planetary_epoch: PlanetaryState,
        chemical_disequilibrium_bits: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_pl = self
            .nodes
            .get(&host_id)
            .map(|n| n.planetary_epoch)
            .unwrap_or(PlanetaryState::AnoxicArchean);

        let is_goe_tr = !prev_pl.is_oxygenated() && planetary_epoch.is_oxygenated();
        let is_gaia2_tr =
            !prev_pl.is_homeostatic_awakened() && planetary_epoch.is_homeostatic_awakened();
        if is_goe_tr {
            self.goe_transitions += 1;
        }
        if is_gaia2_tr {
            self.homeostatic_awakening_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class,
                mro_pathways,
                mro_mtdna_loci,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state,
                cell_count: cell_count.max(1),
                soma_fraction,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state,
                colony_size: colony_size.max(1),
                sterile_caste_fraction,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state,
                cultivar_polyploidy: cultivar_polyploidy.max(1),
                gongylidia_yield,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state,
                cultural_ratchet_tier,
                cultural_repertoire_bits,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state,
                exosomatic_archive_bits,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state,
                verified_theorem_bits,
                distribution_tail_diversity,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch,
                chemical_disequilibrium_bits,
                is_goe_transition: is_goe_tr,
                is_homeostatic_awakening_transition: is_gaia2_tr,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.16,
                leverage_li,
            },
        );
    }

    /// Registers a **Hadean RAF Protocell (`Soup -> RAFProtocell`)** or **Woese Darwinian Threshold
    /// (`Progenote -> LUCA`)** landmark in the reticulate phylogeny!
    #[allow(clippy::too_many_arguments)]
    pub fn record_abiogenesis_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        abiogenesis_stage: AbiogenesisStage,
        raf_closure_score: f64,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_ab = self
            .nodes
            .get(&host_id)
            .map(|n| n.abiogenesis_stage)
            .unwrap_or(AbiogenesisStage::PrebioticMonomerSoup);

        let is_raf_tr = !prev_ab.has_raf_closure() && abiogenesis_stage.has_raf_closure();
        let is_luca_tr = !prev_ab.crossed_darwinian_threshold()
            && abiogenesis_stage.crossed_darwinian_threshold();
        if is_raf_tr {
            self.raf_protocell_transitions += 1;
        }
        if is_luca_tr {
            self.darwinian_threshold_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: 1,
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::None,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::AnoxicArchean,
                chemical_disequilibrium_bits: 2.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage,
                raf_closure_score,
                is_raf_protocell_transition: is_raf_tr,
                is_darwinian_threshold_transition: is_luca_tr,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: 0,
                organelle_loci_remaining: 0,
                egt_nuclear_fraction: 0.0,
                is_egt_lockin: false,
                seq_divergence: 0.14,
                leverage_li,
            },
        );
    }

    /// Registers **D/L Homochiral Symmetry Breaking (`Racemic -> Homochiral`)** or
    /// **Collective Error-Minimizing `SGC` Crystallization (`StatCode -> UniversalSGC`)**!
    #[allow(clippy::too_many_arguments)]
    pub fn record_chirality_code_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        chirality_code_stage: ChiralityCodeStage,
        genetic_code_mp: f64,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_ch = self
            .nodes
            .get(&host_id)
            .map(|n| n.chirality_code_stage)
            .unwrap_or(ChiralityCodeStage::RacemicCrossInhibited);

        let is_homo_tr = !prev_ch.is_homochiral_rna() && chirality_code_stage.is_homochiral_rna();
        let is_sgc_tr =
            !prev_ch.is_crystallized_sgc() && chirality_code_stage.is_crystallized_sgc();
        if is_homo_tr {
            self.homochiral_symmetry_transitions += 1;
        }
        if is_sgc_tr {
            self.sgc_crystallization_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: 1,
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::None,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::AnoxicArchean,
                chemical_disequilibrium_bits: 2.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: if chirality_code_stage.is_crystallized_sgc() {
                    AbiogenesisStage::LucaCellularLineage
                } else {
                    AbiogenesisStage::ProgenoteCommunity
                },
                raf_closure_score: 0.92,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage,
                genetic_code_mp,
                is_homochiral_symmetry_transition: is_homo_tr,
                is_sgc_crystallization_transition: is_sgc_tr,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: 0,
                organelle_loci_remaining: 0,
                egt_nuclear_fraction: 0.0,
                is_egt_lockin: false,
                seq_divergence: 0.14,
                leverage_li,
            },
        );
    }

    /// Registers **Anaerobic `RNR` `RNA -> U-DNA`** or **`ThyX + UDG` `U-DNA -> Megabase T-DNA`**
    /// milestones in the reticulate phylogeny!
    #[allow(clippy::too_many_arguments)]
    pub fn record_rna_to_dna_transition(
        &mut self,
        new_taxon_id: u64,
        parent_taxon_id: u64,
        step: usize,
        species_id: usize,
        guild: TrophicGuild,
        biome: BiomeZone,
        genome_substrate_stage: GenomeSubstrateStage,
        chromosome_capacity_kb: f64,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_taxon_id) {
            parent_taxon_id
        } else {
            self.root_id
        };
        let prev_dn = self
            .nodes
            .get(&host_id)
            .map(|n| n.genome_substrate_stage)
            .unwrap_or(GenomeSubstrateStage::LabileRnaChromosome);

        let is_rnr_tr =
            !prev_dn.has_deoxyribose_rnr() && genome_substrate_stage.has_deoxyribose_rnr();
        let is_udg_tr =
            !prev_dn.has_thymine_udg_repair() && genome_substrate_stage.has_thymine_udg_repair();
        if is_rnr_tr {
            self.rnr_dna_transitions += 1;
        }
        if is_udg_tr {
            self.thymine_udg_transitions += 1;
        }
        self.total_speciation_splits += 1;

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&new_taxon_id) {
                p_node.children.push(new_taxon_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        self.nodes.insert(
            new_taxon_id,
            PhyloNode {
                taxon_id: new_taxon_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild,
                biome,
                is_speciation_split: true,
                guild_transition: None,
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: 1,
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::None,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::AnoxicArchean,
                chemical_disequilibrium_bits: 2.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage,
                chromosome_capacity_kb,
                is_rnr_dna_transition: is_rnr_tr,
                is_thymine_udg_transition: is_udg_tr,
                egt_transferred_loci: 0,
                organelle_loci_remaining: 0,
                egt_nuclear_fraction: 0.0,
                is_egt_lockin: false,
                seq_divergence: 0.14,
                leverage_li,
            },
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_apicoplast_transition(
        &mut self,
        apico_parasite_id: u64,
        parent_holobiont_id: u64,
        step: usize,
        species_id: usize,
        biome: BiomeZone,
        endo_tier: u8,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        apicoplast_pathways: u8,
        leverage_li: f64,
    ) {
        let host_id = if self.nodes.contains_key(&parent_holobiont_id) {
            parent_holobiont_id
        } else {
            self.root_id
        };
        let (prev_guild, prev_mro, prev_mro_paths, prev_mtdna) = self
            .nodes
            .get(&host_id)
            .map(|n| (n.guild, n.mro_class, n.mro_pathways, n.mro_mtdna_loci))
            .unwrap_or((TrophicGuild::Holobiont, MroClass::AerobicMitochondrion, 4, 3));

        self.apicoplast_transition_events += 1;
        self.total_speciation_splits += 1;
        self.record_guild_shift_counters(prev_guild, TrophicGuild::PredatorParasite);

        if let Some(p_node) = self.nodes.get_mut(&host_id) {
            if !p_node.children.contains(&apico_parasite_id) {
                p_node.children.push(apico_parasite_id);
            }
            if p_node.extant {
                p_node.extant = false;
                p_node.death_step = Some(step);
            }
        }

        let tot = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (egt_transferred_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            apico_parasite_id,
            PhyloNode {
                taxon_id: apico_parasite_id,
                parent_id: Some(host_id),
                donor_id: None,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild: TrophicGuild::PredatorParasite,
                biome,
                is_speciation_split: true,
                guild_transition: Some((prev_guild, TrophicGuild::PredatorParasite)),
                reticulation_kind: None,
                hybrid_id: None,
                genome_share_gamma: 0.0,
                endo_tier: endo_tier.max(1),
                nucleomorph_loci: 0,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: true,
                apicoplast_pathways,
                is_apicoplast_transition: true,
                is_apicoplast_complete_loss: false,
                mro_class: prev_mro,
                mro_pathways: prev_mro_paths,
                mro_mtdna_loci: prev_mtdna,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: 0.22,
                leverage_li,
            },
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_birth_full(
        &mut self,
        child_id: u64,
        parent_id: u64,
        donor_id: Option<u64>,
        step: usize,
        child_species: usize,
        child_guild: TrophicGuild,
        child_biome: BiomeZone,
        parent_species: usize,
        parent_guild: TrophicGuild,
        seq_divergence: f64,
        speciation_threshold: f64,
        leverage_li: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
        endo_tier: u8,
        nucleomorph_loci: usize,
        nucleomorph_eliminated: bool,
    ) {
        let attach_parent = if self.nodes.contains_key(&parent_id) {
            parent_id
        } else {
            self.root_id
        };
        let (
            parent_egt_frac,
            parent_egt_loci,
            parent_nm_loci,
            parent_hgt,
            parent_kp,
            parent_apico,
            parent_apico_paths,
            parent_mro,
            parent_mro_paths,
            parent_mtdna,
        ) = self
            .nodes
            .get(&attach_parent)
            .map(|n| {
                (
                    n.egt_nuclear_fraction,
                    n.egt_transferred_loci,
                    n.nucleomorph_loci,
                    n.klepto_hgt_loci,
                    n.has_kleptoplast,
                    n.has_apicoplast,
                    n.apicoplast_pathways,
                    n.mro_class,
                    n.mro_pathways,
                    n.mro_mtdna_loci,
                )
            })
            .unwrap_or((0.0, 0, 0, 0, false, false, 0, MroClass::None, 0, 0));

        if egt_transferred_loci > parent_egt_loci {
            self.egt_transfer_events += egt_transferred_loci - parent_egt_loci;
        }

        let tot_endo = egt_transferred_loci + organelle_loci_remaining;
        let egt_frac = if (child_guild == TrophicGuild::Holobiont
            || parent_apico
            || parent_mro != MroClass::None)
            && tot_endo > 0
        {
            (egt_transferred_loci as f64) / (tot_endo as f64)
        } else {
            0.0
        };
        let crossed_lockin =
            child_guild == TrophicGuild::Holobiont && egt_frac >= 0.50 && parent_egt_frac < 0.50;
        if crossed_lockin {
            self.egt_lockin_events += 1;
        }

        let crossed_nm_loss = child_guild == TrophicGuild::Holobiont
            && endo_tier >= 2
            && nucleomorph_loci == 0
            && (parent_nm_loci > 0 || nucleomorph_eliminated);
        if crossed_nm_loss {
            self.nucleomorph_loss_events += 1;
        }

        let is_holo_to_apico = parent_guild == TrophicGuild::Holobiont
            && child_guild == TrophicGuild::PredatorParasite
            && organelle_loci_remaining >= 1
            && parent_mro != MroClass::Hydrogenosome;
        if is_holo_to_apico {
            self.apicoplast_transition_events += 1;
        }
        let child_has_apico = is_holo_to_apico
            || (parent_apico
                && child_guild == TrophicGuild::PredatorParasite
                && organelle_loci_remaining >= 1);
        let crossed_apico_loss =
            parent_apico && !child_has_apico && child_guild == TrophicGuild::PredatorParasite;
        if crossed_apico_loss {
            self.apicoplast_loss_events += 1;
        }

        let valid_donor =
            donor_id.filter(|did| self.nodes.contains_key(did) && *did != attach_parent);

        let is_spec_split =
            child_species != parent_species || seq_divergence >= speciation_threshold * 0.65;
        if is_spec_split {
            self.total_speciation_splits += 1;
        }

        let guild_transition = if child_guild != parent_guild {
            self.record_guild_shift_counters(parent_guild, child_guild);
            Some((parent_guild, child_guild))
        } else {
            None
        };

        if let Some(p_node) = self.nodes.get_mut(&attach_parent) {
            if !p_node.children.contains(&child_id) {
                p_node.children.push(child_id);
            }
        }
        if let Some(did) = valid_donor {
            self.horizontal_transfers += 1;
            if let Some(d_node) = self.nodes.get_mut(&did) {
                if !d_node.symbio_children.contains(&child_id) {
                    d_node.symbio_children.push(child_id);
                }
            }
        }

        let effective_tier = if child_guild == TrophicGuild::Holobiont
            || child_has_apico
            || parent_mro != MroClass::None
        {
            endo_tier.max(1)
        } else {
            0
        };

        let child_mro = if parent_mro != MroClass::None {
            parent_mro
        } else if child_guild == TrophicGuild::Holobiont {
            MroClass::AerobicMitochondrion
        } else {
            MroClass::None
        };
        let child_mro_paths = if parent_mro != MroClass::None {
            parent_mro_paths
        } else if child_guild == TrophicGuild::Holobiont {
            4
        } else {
            0
        };
        let child_mtdna = if parent_mro != MroClass::None {
            parent_mtdna
        } else if child_guild == TrophicGuild::Holobiont {
            organelle_loci_remaining.min(3)
        } else {
            0
        };

        self.nodes.insert(
            child_id,
            PhyloNode {
                taxon_id: child_id,
                parent_id: Some(attach_parent),
                donor_id: valid_donor,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id: child_species,
                guild: child_guild,
                biome: child_biome,
                is_speciation_split: is_spec_split,
                guild_transition,
                reticulation_kind: valid_donor.map(|_| ReticulationKind::HorizontalTransfer),
                hybrid_id: None,
                genome_share_gamma: if valid_donor.is_some() { 0.14 } else { 0.0 },
                endo_tier: effective_tier,
                nucleomorph_loci,
                is_nucleomorph_loss: crossed_nm_loss,
                klepto_hgt_loci: parent_hgt,
                has_kleptoplast: parent_kp
                    && child_guild == TrophicGuild::PredatorParasite
                    && !child_has_apico,
                is_klepto_tipping_point: false,
                has_apicoplast: child_has_apico,
                apicoplast_pathways: if child_has_apico {
                    parent_apico_paths.max(2)
                } else if child_guild == TrophicGuild::Holobiont {
                    4
                } else {
                    0
                },
                is_apicoplast_transition: is_holo_to_apico,
                is_apicoplast_complete_loss: crossed_apico_loss,
                mro_class: child_mro,
                mro_pathways: child_mro_paths,
                mro_mtdna_loci: child_mtdna,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: crossed_lockin,
                seq_divergence,
                leverage_li,
            },
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_birth_with_egt(
        &mut self,
        child_id: u64,
        parent_id: u64,
        donor_id: Option<u64>,
        step: usize,
        child_species: usize,
        child_guild: TrophicGuild,
        child_biome: BiomeZone,
        parent_species: usize,
        parent_guild: TrophicGuild,
        seq_divergence: f64,
        speciation_threshold: f64,
        leverage_li: f64,
        egt_transferred_loci: usize,
        organelle_loci_remaining: usize,
    ) {
        let inherited_tier = if child_guild == TrophicGuild::Holobiont {
            self.nodes
                .get(&parent_id)
                .map(|n| n.endo_tier.max(1))
                .unwrap_or(1)
        } else {
            0
        };
        let inherited_nm = self
            .nodes
            .get(&parent_id)
            .map(|n| n.nucleomorph_loci)
            .unwrap_or(0);
        self.record_birth_full(
            child_id,
            parent_id,
            donor_id,
            step,
            child_species,
            child_guild,
            child_biome,
            parent_species,
            parent_guild,
            seq_divergence,
            speciation_threshold,
            leverage_li,
            egt_transferred_loci,
            organelle_loci_remaining,
            inherited_tier,
            inherited_nm,
            false,
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_birth(
        &mut self,
        child_id: u64,
        parent_id: u64,
        donor_id: Option<u64>,
        step: usize,
        child_species: usize,
        child_guild: TrophicGuild,
        child_biome: BiomeZone,
        parent_species: usize,
        parent_guild: TrophicGuild,
        seq_divergence: f64,
        speciation_threshold: f64,
        leverage_li: f64,
    ) {
        self.record_birth_full(
            child_id,
            parent_id,
            donor_id,
            step,
            child_species,
            child_guild,
            child_biome,
            parent_species,
            parent_guild,
            seq_divergence,
            speciation_threshold,
            leverage_li,
            0,
            0,
            u8::from(child_guild == TrophicGuild::Holobiont),
            0,
            false,
        );
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_klepto_tipping_fixation(
        &mut self,
        new_holobiont_id: u64,
        klepto_host_parent_id: u64,
        plastid_donor_id: Option<u64>,
        step: usize,
        species_id: usize,
        biome: BiomeZone,
        endo_tier: u8,
        genome_share_gamma: f64,
        klepto_hgt_loci: usize,
        organelle_loci_remaining: usize,
        nucleomorph_loci: usize,
        leverage_li: f64,
    ) -> usize {
        let host_id = if self.nodes.contains_key(&klepto_host_parent_id) {
            klepto_host_parent_id
        } else {
            self.root_id
        };
        let prev_guild = self
            .nodes
            .get(&host_id)
            .map(|n| n.guild)
            .unwrap_or(TrophicGuild::PredatorParasite);

        let valid_donor =
            plastid_donor_id.filter(|did| self.nodes.contains_key(did) && *did != host_id);

        let h_idx = self.next_hybrid_id;
        self.next_hybrid_id += 1;
        self.endosymbiotic_fusions += 1;
        self.klepto_fixation_events += 1;
        let tier = endo_tier.clamp(1, 3);
        if tier == 2 {
            self.secondary_fusions += 1;
        } else if tier >= 3 {
            self.tertiary_fusions += 1;
        }
        self.record_guild_shift_counters(prev_guild, TrophicGuild::Holobiont);

        if let Some(h_node) = self.nodes.get_mut(&host_id) {
            if !h_node.children.contains(&new_holobiont_id) {
                h_node.children.push(new_holobiont_id);
            }
            if h_node.extant {
                h_node.extant = false;
                h_node.death_step = Some(step);
            }
        }
        if let Some(did) = valid_donor {
            if let Some(d_node) = self.nodes.get_mut(&did) {
                if !d_node.symbio_children.contains(&new_holobiont_id) {
                    d_node.symbio_children.push(new_holobiont_id);
                }
            }
        }

        let tot = klepto_hgt_loci + organelle_loci_remaining;
        let egt_frac = if tot > 0 {
            (klepto_hgt_loci as f64) / (tot as f64)
        } else {
            0.0
        };

        self.nodes.insert(
            new_holobiont_id,
            PhyloNode {
                taxon_id: new_holobiont_id,
                parent_id: Some(host_id),
                donor_id: valid_donor,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id,
                guild: TrophicGuild::Holobiont,
                biome,
                is_speciation_split: true,
                guild_transition: Some((prev_guild, TrophicGuild::Holobiont)),
                reticulation_kind: Some(ReticulationKind::KleptoTippingFixation),
                hybrid_id: Some(h_idx),
                genome_share_gamma: genome_share_gamma.clamp(0.10, 0.90),
                endo_tier: tier,
                nucleomorph_loci,
                is_nucleomorph_loss: false,
                klepto_hgt_loci,
                has_kleptoplast: false,
                is_klepto_tipping_point: true,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::AerobicMitochondrion,
                mro_pathways: 4,
                mro_mtdna_loci: 3,
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: klepto_hgt_loci,
                organelle_loci_remaining,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: genome_share_gamma,
                leverage_li,
            },
        );
        self.total_speciation_splits += 1;
        h_idx
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_serial_endosymbiotic_fusion(
        &mut self,
        holobiont_id: u64,
        host_parent_id: u64,
        symbiont_donor_id: u64,
        step: usize,
        holobiont_species: usize,
        biome: BiomeZone,
        endo_tier: u8,
        genome_share_gamma: f64,
        initial_organelle_len: usize,
        initial_nucleomorph_len: usize,
        initial_egt_loci: usize,
        leverage_li: f64,
    ) -> usize {
        let host_id = if self.nodes.contains_key(&host_parent_id) {
            host_parent_id
        } else {
            self.root_id
        };
        let host_prev_guild = self
            .nodes
            .get(&host_id)
            .map(|n| n.guild)
            .unwrap_or(TrophicGuild::Producer);

        let sym_id = if self.nodes.contains_key(&symbiont_donor_id) && symbiont_donor_id != host_id
        {
            Some(symbiont_donor_id)
        } else {
            None
        };

        let tier = endo_tier.clamp(1, 3);
        let h_idx = self.next_hybrid_id;
        self.next_hybrid_id += 1;
        self.endosymbiotic_fusions += 1;
        if tier == 2 {
            self.secondary_fusions += 1;
        } else if tier >= 3 {
            self.tertiary_fusions += 1;
        }
        self.record_guild_shift_counters(host_prev_guild, TrophicGuild::Holobiont);

        if let Some(h_node) = self.nodes.get_mut(&host_id) {
            if !h_node.children.contains(&holobiont_id) {
                h_node.children.push(holobiont_id);
            }
            if h_node.extant {
                h_node.extant = false;
                h_node.death_step = Some(step);
            }
        }

        if let Some(sid) = sym_id {
            if let Some(s_node) = self.nodes.get_mut(&sid) {
                if !s_node.symbio_children.contains(&holobiont_id) {
                    s_node.symbio_children.push(holobiont_id);
                }
                if s_node.extant {
                    s_node.extant = false;
                    s_node.death_step = Some(step);
                }
            }
        }

        let ret_kind = match tier {
            1 => ReticulationKind::EndosymbioticFusion,
            2 => ReticulationKind::SecondaryEndosymbiosis,
            _ => ReticulationKind::TertiaryEndosymbiosis,
        };
        let tot_endo = initial_egt_loci + initial_organelle_len;
        let egt_frac = if tot_endo > 0 {
            (initial_egt_loci as f64) / (tot_endo as f64)
        } else {
            0.0
        };
        let g_shift = if host_prev_guild != TrophicGuild::Holobiont {
            Some((host_prev_guild, TrophicGuild::Holobiont))
        } else {
            None
        };

        self.nodes.insert(
            holobiont_id,
            PhyloNode {
                taxon_id: holobiont_id,
                parent_id: Some(host_id),
                donor_id: sym_id,
                children: Vec::new(),
                symbio_children: Vec::new(),
                birth_step: step,
                death_step: None,
                extant: true,
                species_id: holobiont_species,
                guild: TrophicGuild::Holobiont,
                biome,
                is_speciation_split: true,
                guild_transition: g_shift,
                reticulation_kind: Some(ret_kind),
                hybrid_id: Some(h_idx),
                genome_share_gamma: genome_share_gamma.clamp(0.10, 0.90),
                endo_tier: tier,
                nucleomorph_loci: initial_nucleomorph_len,
                is_nucleomorph_loss: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                is_klepto_tipping_point: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                is_apicoplast_transition: false,
                is_apicoplast_complete_loss: false,
                mro_class: MroClass::AerobicMitochondrion,
                mro_pathways: 4,
                mro_mtdna_loci: initial_organelle_len.min(3).max(1),
                is_hydrogenosome_transition: false,
                is_mitosome_transition: false,
                is_mro_complete_loss: false,
                multicell_state: MulticellularState::Unicellular,
                cell_count: 1,
                soma_fraction: 0.0,
                is_clonal_multicell_transition: false,
                is_germ_soma_transition: false,
                eusocial_state: EusocialState::Solitary,
                colony_size: 1,
                sterile_caste_fraction: 0.0,
                is_eusocial_transition: false,
                is_superorganism_transition: false,
                agriculture_state: AgricultureState::HunterForager,
                cultivar_polyploidy: 1,
                gongylidia_yield: 0.0,
                is_agriculture_transition: false,
                is_obligate_crop_transition: false,
                cultural_state: CulturalState::AsocialTrialError,
                cultural_ratchet_tier: 0,
                cultural_repertoire_bits: 0.0,
                is_cultural_ratchet_transition: false,
                is_symbolic_teaching_transition: false,
                exosomatic_state: ExosomaticState::EndosomaticBound,
                exosomatic_archive_bits: 0.0,
                is_written_archive_transition: false,
                is_peer_institution_transition: false,
                machine_synthesis_state: MachineSynthesisState::ManualArchival,
                verified_theorem_bits: 0.0,
                distribution_tail_diversity: 1.0,
                is_machine_symbiosis_transition: false,
                is_formal_verifier_transition: false,
                planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                chemical_disequilibrium_bits: 14.5,
                is_goe_transition: false,
                is_homeostatic_awakening_transition: false,
                abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                raf_closure_score: 0.96,
                is_raf_protocell_transition: false,
                is_darwinian_threshold_transition: false,
                chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                genetic_code_mp: 0.79,
                is_homochiral_symmetry_transition: false,
                is_sgc_crystallization_transition: false,
                genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                chromosome_capacity_kb: 1850.0,
                is_rnr_dna_transition: false,
                is_thymine_udg_transition: false,
                egt_transferred_loci: initial_egt_loci,
                organelle_loci_remaining: initial_organelle_len,
                egt_nuclear_fraction: egt_frac,
                is_egt_lockin: egt_frac >= 0.50,
                seq_divergence: genome_share_gamma,
                leverage_li,
            },
        );
        self.total_speciation_splits += 1;
        h_idx
    }

    #[allow(clippy::too_many_arguments)]
    pub fn record_endosymbiotic_fusion(
        &mut self,
        holobiont_id: u64,
        host_parent_id: u64,
        symbiont_donor_id: u64,
        step: usize,
        holobiont_species: usize,
        biome: BiomeZone,
        genome_share_gamma: f64,
        initial_organelle_len: usize,
        leverage_li: f64,
    ) -> usize {
        self.record_serial_endosymbiotic_fusion(
            holobiont_id,
            host_parent_id,
            symbiont_donor_id,
            step,
            holobiont_species,
            biome,
            1,
            genome_share_gamma,
            initial_organelle_len,
            0,
            0,
            leverage_li,
        )
    }

    pub fn sync_extant_population_full(
        &mut self,
        current_step: usize,
        extant_units: &[ExtantPhyloState],
    ) {
        let mut extant_map = BTreeMap::new();
        for &st in extant_units {
            extant_map.insert(st.taxon_id, st);
            if !self.nodes.contains_key(&st.taxon_id) {
                if let Some(r_node) = self.nodes.get_mut(&self.root_id) {
                    r_node.children.push(st.taxon_id);
                }
                let tot = st.egt_loci + st.organelle_rem;
                let frac = if (st.guild == TrophicGuild::Holobiont
                    || st.has_apicoplast
                    || st.mro_class != MroClass::None)
                    && tot > 0
                {
                    (st.egt_loci as f64) / (tot as f64)
                } else {
                    0.0
                };
                self.nodes.insert(
                    st.taxon_id,
                    PhyloNode {
                        taxon_id: st.taxon_id,
                        parent_id: Some(self.root_id),
                        donor_id: None,
                        children: Vec::new(),
                        symbio_children: Vec::new(),
                        birth_step: current_step,
                        death_step: None,
                        extant: true,
                        species_id: st.species_id,
                        guild: st.guild,
                        biome: st.biome,
                        is_speciation_split: false,
                        guild_transition: None,
                        reticulation_kind: None,
                        hybrid_id: None,
                        genome_share_gamma: 0.0,
                        endo_tier: st.endo_tier,
                        nucleomorph_loci: st.nucleomorph_loci,
                        is_nucleomorph_loss: st.nucleomorph_eliminated,
                        klepto_hgt_loci: st.klepto_hgt_loci,
                        has_kleptoplast: st.has_kleptoplast,
                        is_klepto_tipping_point: false,
                        has_apicoplast: st.has_apicoplast,
                        apicoplast_pathways: st.apicoplast_pathways,
                        is_apicoplast_transition: false,
                        is_apicoplast_complete_loss: st.apicoplast_eliminated,
                        mro_class: st.mro_class,
                        mro_pathways: st.mro_pathways,
                        mro_mtdna_loci: st.mro_mtdna_loci,
                        is_hydrogenosome_transition: false,
                        is_mitosome_transition: false,
                        is_mro_complete_loss: st.mro_class == MroClass::Amiteukaryote,
                        multicell_state: st.multicell_state,
                        cell_count: st.cell_count.max(1),
                        soma_fraction: st.soma_fraction,
                        is_clonal_multicell_transition: false,
                        is_germ_soma_transition: false,
                        eusocial_state: st.eusocial_state,
                        colony_size: st.colony_size.max(1),
                        sterile_caste_fraction: st.sterile_caste_fraction,
                        is_eusocial_transition: false,
                        is_superorganism_transition: false,
                        agriculture_state: st.agriculture_state,
                        cultivar_polyploidy: st.cultivar_polyploidy.max(1),
                        gongylidia_yield: st.gongylidia_yield,
                        is_agriculture_transition: false,
                        is_obligate_crop_transition: false,
                        cultural_state: st.cultural_state,
                        cultural_ratchet_tier: st.cultural_ratchet_tier,
                        cultural_repertoire_bits: st.cultural_repertoire_bits,
                        is_cultural_ratchet_transition: false,
                        is_symbolic_teaching_transition: false,
                        exosomatic_state: st.exosomatic_state,
                        exosomatic_archive_bits: st.exosomatic_archive_bits,
                        is_written_archive_transition: false,
                        is_peer_institution_transition: false,
                        machine_synthesis_state: st.machine_synthesis_state,
                        verified_theorem_bits: st.verified_theorem_bits,
                        distribution_tail_diversity: st.distribution_tail_diversity,
                        is_machine_symbiosis_transition: false,
                        is_formal_verifier_transition: false,
                        planetary_epoch: st.planetary_epoch,
                        chemical_disequilibrium_bits: st.chemical_disequilibrium_bits,
                        is_goe_transition: false,
                        is_homeostatic_awakening_transition: false,
                        abiogenesis_stage: st.abiogenesis_stage,
                        raf_closure_score: st.raf_closure_score,
                        is_raf_protocell_transition: false,
                        is_darwinian_threshold_transition: false,
                        chirality_code_stage: st.chirality_code_stage,
                        genetic_code_mp: st.genetic_code_mp,
                        is_homochiral_symmetry_transition: false,
                        is_sgc_crystallization_transition: false,
                        genome_substrate_stage: st.genome_substrate_stage,
                        chromosome_capacity_kb: st.chromosome_capacity_kb,
                        is_rnr_dna_transition: false,
                        is_thymine_udg_transition: false,
                        egt_transferred_loci: st.egt_loci,
                        organelle_loci_remaining: st.organelle_rem,
                        egt_nuclear_fraction: frac,
                        is_egt_lockin: frac >= 0.50,
                        seq_divergence: 0.0,
                        leverage_li: st.leverage_li,
                    },
                );
            }
        }

        let all_ids: Vec<u64> = self.nodes.keys().copied().collect();
        for id in all_ids {
            if let Some(&st) = extant_map.get(&id) {
                let mut shift_detected = None;
                let mut split_detected = false;
                let mut egt_delta = 0usize;
                let mut hgt_delta = 0usize;
                let mut lockin_crossed = false;
                let mut nm_loss_crossed = false;
                let mut apico_loss_crossed = false;
                let mut hydro_crossed = false;
                let mut mitos_crossed = false;
                let mut mro_zero_crossed = false;
                let mut clonal_mc_crossed = false;
                let mut gs_mc_crossed = false;
                let mut eusoc_crossed = false;
                let mut superorg_crossed = false;
                let mut ag_crossed = false;
                let mut crop_crossed = false;
                let mut ratchet_crossed = false;
                let mut sym_crossed = false;
                let mut arch_crossed = false;
                let mut inst_crossed = false;
                let mut ms_crossed = false;
                let mut formal_crossed = false;
                let mut goe_crossed = false;
                let mut gaia2_crossed = false;
                let mut raf_crossed = false;
                let mut luca_crossed = false;
                let mut homo_crossed = false;
                let mut sgc_crossed = false;
                let mut rnr_crossed = false;
                let mut udg_crossed = false;

                if let Some(node) = self.nodes.get_mut(&id) {
                    node.extant = true;
                    node.death_step = None;
                    if node.guild != st.guild && node.guild_transition.is_none() {
                        shift_detected = Some((node.guild, st.guild));
                        node.guild_transition = shift_detected;
                    }
                    if node.species_id != st.species_id
                        && !node.is_speciation_split
                        && node.birth_step > 0
                    {
                        node.is_speciation_split = true;
                        split_detected = true;
                    }
                    if st.egt_loci > node.egt_transferred_loci {
                        egt_delta = st.egt_loci - node.egt_transferred_loci;
                    }
                    if st.klepto_hgt_loci > node.klepto_hgt_loci {
                        hgt_delta = st.klepto_hgt_loci - node.klepto_hgt_loci;
                    }
                    let tot = st.egt_loci + st.organelle_rem;
                    let new_frac = if (st.guild == TrophicGuild::Holobiont
                        || st.has_apicoplast
                        || st.mro_class != MroClass::None)
                        && tot > 0
                    {
                        (st.egt_loci as f64) / (tot as f64)
                    } else {
                        0.0
                    };
                    if st.guild == TrophicGuild::Holobiont
                        && new_frac >= 0.50
                        && node.egt_nuclear_fraction < 0.50
                        && !node.is_egt_lockin
                    {
                        node.is_egt_lockin = true;
                        lockin_crossed = true;
                    }
                    if st.guild == TrophicGuild::Holobiont
                        && st.endo_tier >= 2
                        && st.nucleomorph_loci == 0
                        && (node.nucleomorph_loci > 0 || st.nucleomorph_eliminated)
                        && !node.is_nucleomorph_loss
                    {
                        node.is_nucleomorph_loss = true;
                        nm_loss_crossed = true;
                    }
                    if !st.has_apicoplast
                        && (node.has_apicoplast || st.apicoplast_eliminated)
                        && !node.is_apicoplast_complete_loss
                    {
                        node.is_apicoplast_complete_loss = true;
                        apico_loss_crossed = true;
                    }
                    if st.mro_class == MroClass::Hydrogenosome
                        && node.mro_class != MroClass::Hydrogenosome
                        && !node.is_hydrogenosome_transition
                    {
                        node.is_hydrogenosome_transition = true;
                        hydro_crossed = true;
                    }
                    if st.mro_class == MroClass::Mitosome
                        && node.mro_class != MroClass::Mitosome
                        && !node.is_mitosome_transition
                    {
                        node.is_mitosome_transition = true;
                        mitos_crossed = true;
                    }
                    if st.mro_class == MroClass::Amiteukaryote
                        && node.mro_class != MroClass::Amiteukaryote
                        && !node.is_mro_complete_loss
                    {
                        node.is_mro_complete_loss = true;
                        mro_zero_crossed = true;
                    }
                    node.endo_tier = st.endo_tier;
                    node.nucleomorph_loci = st.nucleomorph_loci;
                    node.klepto_hgt_loci = st.klepto_hgt_loci;
                    node.has_kleptoplast = st.has_kleptoplast;
                    node.has_apicoplast = st.has_apicoplast;
                    node.apicoplast_pathways = st.apicoplast_pathways;
                    if !node.multicell_state.is_multicellular()
                        && st.multicell_state.is_multicellular()
                        && !node.is_clonal_multicell_transition
                    {
                        node.is_clonal_multicell_transition = true;
                        clonal_mc_crossed = true;
                    }
                    if node.multicell_state != MulticellularState::GermSomaMulticellular
                        && st.multicell_state == MulticellularState::GermSomaMulticellular
                        && !node.is_germ_soma_transition
                    {
                        node.is_germ_soma_transition = true;
                        gs_mc_crossed = true;
                    }
                    if !node.eusocial_state.is_eusocial()
                        && st.eusocial_state.is_eusocial()
                        && !node.is_eusocial_transition
                    {
                        node.is_eusocial_transition = true;
                        eusoc_crossed = true;
                    }
                    if node.eusocial_state != EusocialState::ObligateSuperorganism
                        && st.eusocial_state == EusocialState::ObligateSuperorganism
                        && !node.is_superorganism_transition
                    {
                        node.is_superorganism_transition = true;
                        superorg_crossed = true;
                    }
                    if !node.agriculture_state.is_farming()
                        && st.agriculture_state.is_farming()
                        && !node.is_agriculture_transition
                    {
                        node.is_agriculture_transition = true;
                        ag_crossed = true;
                    }
                    if !node.agriculture_state.is_obligate_monoculture()
                        && st.agriculture_state.is_obligate_monoculture()
                        && !node.is_obligate_crop_transition
                    {
                        node.is_obligate_crop_transition = true;
                        crop_crossed = true;
                    }
                    if !node.cultural_state.is_cumulative_ratchet()
                        && st.cultural_state.is_cumulative_ratchet()
                        && !node.is_cultural_ratchet_transition
                    {
                        node.is_cultural_ratchet_transition = true;
                        ratchet_crossed = true;
                    }
                    if node.cultural_state != CulturalState::SymbolicCumulativeCulture
                        && st.cultural_state == CulturalState::SymbolicCumulativeCulture
                        && !node.is_symbolic_teaching_transition
                    {
                        node.is_symbolic_teaching_transition = true;
                        sym_crossed = true;
                    }
                    if !node.exosomatic_state.has_written_archive()
                        && st.exosomatic_state.has_written_archive()
                        && !node.is_written_archive_transition
                    {
                        node.is_written_archive_transition = true;
                        arch_crossed = true;
                    }
                    if node.exosomatic_state != ExosomaticState::PeerVerifiedInstitution
                        && st.exosomatic_state == ExosomaticState::PeerVerifiedInstitution
                        && !node.is_peer_institution_transition
                    {
                        node.is_peer_institution_transition = true;
                        inst_crossed = true;
                    }
                    if !node.machine_synthesis_state.has_computational_symbiont()
                        && st.machine_synthesis_state.has_computational_symbiont()
                        && !node.is_machine_symbiosis_transition
                    {
                        node.is_machine_symbiosis_transition = true;
                        ms_crossed = true;
                    }
                    if !node.machine_synthesis_state.is_formally_verified()
                        && st.machine_synthesis_state.is_formally_verified()
                        && !node.is_formal_verifier_transition
                    {
                        node.is_formal_verifier_transition = true;
                        formal_crossed = true;
                    }
                    if !node.planetary_epoch.is_oxygenated()
                        && st.planetary_epoch.is_oxygenated()
                        && !node.is_goe_transition
                    {
                        node.is_goe_transition = true;
                        goe_crossed = true;
                    }
                    if !node.planetary_epoch.is_homeostatic_awakened()
                        && st.planetary_epoch.is_homeostatic_awakened()
                        && !node.is_homeostatic_awakening_transition
                    {
                        node.is_homeostatic_awakening_transition = true;
                        gaia2_crossed = true;
                    }
                    if !node.abiogenesis_stage.has_raf_closure()
                        && st.abiogenesis_stage.has_raf_closure()
                        && !node.is_raf_protocell_transition
                    {
                        node.is_raf_protocell_transition = true;
                        raf_crossed = true;
                    }
                    if !node.abiogenesis_stage.crossed_darwinian_threshold()
                        && st.abiogenesis_stage.crossed_darwinian_threshold()
                        && !node.is_darwinian_threshold_transition
                    {
                        node.is_darwinian_threshold_transition = true;
                        luca_crossed = true;
                    }
                    if !node.chirality_code_stage.is_homochiral_rna()
                        && st.chirality_code_stage.is_homochiral_rna()
                        && !node.is_homochiral_symmetry_transition
                    {
                        node.is_homochiral_symmetry_transition = true;
                        homo_crossed = true;
                    }
                    if !node.chirality_code_stage.is_crystallized_sgc()
                        && st.chirality_code_stage.is_crystallized_sgc()
                        && !node.is_sgc_crystallization_transition
                    {
                        node.is_sgc_crystallization_transition = true;
                        sgc_crossed = true;
                    }
                    if !node.genome_substrate_stage.has_deoxyribose_rnr()
                        && st.genome_substrate_stage.has_deoxyribose_rnr()
                        && !node.is_rnr_dna_transition
                    {
                        node.is_rnr_dna_transition = true;
                        rnr_crossed = true;
                    }
                    if !node.genome_substrate_stage.has_thymine_udg_repair()
                        && st.genome_substrate_stage.has_thymine_udg_repair()
                        && !node.is_thymine_udg_transition
                    {
                        node.is_thymine_udg_transition = true;
                        udg_crossed = true;
                    }
                    node.genome_substrate_stage = st.genome_substrate_stage;
                    node.chromosome_capacity_kb = st.chromosome_capacity_kb;
                    node.chirality_code_stage = st.chirality_code_stage;
                    node.genetic_code_mp = st.genetic_code_mp;
                    node.abiogenesis_stage = st.abiogenesis_stage;
                    node.raf_closure_score = st.raf_closure_score;
                    node.planetary_epoch = st.planetary_epoch;
                    node.chemical_disequilibrium_bits = st.chemical_disequilibrium_bits;
                    node.machine_synthesis_state = st.machine_synthesis_state;
                    node.verified_theorem_bits = st.verified_theorem_bits;
                    node.distribution_tail_diversity = st.distribution_tail_diversity;
                    node.exosomatic_state = st.exosomatic_state;
                    node.exosomatic_archive_bits = st.exosomatic_archive_bits;
                    node.cultural_state = st.cultural_state;
                    node.cultural_ratchet_tier = st.cultural_ratchet_tier;
                    node.cultural_repertoire_bits = st.cultural_repertoire_bits;
                    node.agriculture_state = st.agriculture_state;
                    node.cultivar_polyploidy = st.cultivar_polyploidy.max(1);
                    node.gongylidia_yield = st.gongylidia_yield;
                    node.eusocial_state = st.eusocial_state;
                    node.colony_size = st.colony_size.max(1);
                    node.sterile_caste_fraction = st.sterile_caste_fraction;
                    node.multicell_state = st.multicell_state;
                    node.cell_count = st.cell_count.max(1);
                    node.soma_fraction = st.soma_fraction;
                    node.mro_class = st.mro_class;
                    node.mro_pathways = st.mro_pathways;
                    node.mro_mtdna_loci = st.mro_mtdna_loci;
                    node.egt_transferred_loci = st.egt_loci;
                    node.organelle_loci_remaining = st.organelle_rem;
                    node.egt_nuclear_fraction = new_frac;
                    node.species_id = st.species_id;
                    node.guild = st.guild;
                    node.biome = st.biome;
                    node.leverage_li = st.leverage_li;
                }
                if let Some((from_g, to_g)) = shift_detected {
                    self.record_guild_shift_counters(from_g, to_g);
                }
                if split_detected {
                    self.total_speciation_splits += 1;
                }
                self.egt_transfer_events += egt_delta;
                self.klepto_hgt_events += hgt_delta;
                if lockin_crossed {
                    self.egt_lockin_events += 1;
                }
                if nm_loss_crossed {
                    self.nucleomorph_loss_events += 1;
                }
                if apico_loss_crossed {
                    self.apicoplast_loss_events += 1;
                }
                if hydro_crossed {
                    self.hydrogenosome_transition_events += 1;
                }
                if mitos_crossed {
                    self.mitosome_transition_events += 1;
                }
                if mro_zero_crossed {
                    self.mro_complete_loss_events += 1;
                }
                if clonal_mc_crossed {
                    self.clonal_multicell_transitions += 1;
                }
                if gs_mc_crossed {
                    self.germ_soma_transitions += 1;
                }
                if eusoc_crossed {
                    self.eusocial_transitions += 1;
                }
                if superorg_crossed {
                    self.superorganism_transitions += 1;
                }
                if ag_crossed {
                    self.agriculture_transitions += 1;
                }
                if crop_crossed {
                    self.obligate_crop_transitions += 1;
                }
                if ratchet_crossed {
                    self.cultural_ratchet_transitions += 1;
                }
                if sym_crossed {
                    self.symbolic_teaching_transitions += 1;
                }
                if arch_crossed {
                    self.written_archive_transitions += 1;
                }
                if inst_crossed {
                    self.peer_institution_transitions += 1;
                }
                if ms_crossed {
                    self.machine_symbiosis_transitions += 1;
                }
                if formal_crossed {
                    self.formal_verifier_transitions += 1;
                }
                if goe_crossed {
                    self.goe_transitions += 1;
                }
                if gaia2_crossed {
                    self.homeostatic_awakening_transitions += 1;
                }
                if raf_crossed {
                    self.raf_protocell_transitions += 1;
                }
                if luca_crossed {
                    self.darwinian_threshold_transitions += 1;
                }
                if homo_crossed {
                    self.homochiral_symmetry_transitions += 1;
                }
                if sgc_crossed {
                    self.sgc_crystallization_transitions += 1;
                }
                if rnr_crossed {
                    self.rnr_dna_transitions += 1;
                }
                if udg_crossed {
                    self.thymine_udg_transitions += 1;
                }
            } else if let Some(node) = self.nodes.get_mut(&id) {
                if node.extant {
                    node.extant = false;
                    node.death_step = Some(current_step);
                }
            }
        }

        self.prune_extinct_and_collapse_chains();
    }

    pub fn sync_extant_population_with_egt(
        &mut self,
        current_step: usize,
        extant_units: &[(u64, usize, TrophicGuild, BiomeZone, f64, usize, usize)],
    ) {
        let states: Vec<ExtantPhyloState> = extant_units
            .iter()
            .map(|&(id, sp, g, b, lev, egt_loci, org_rem)| {
                let (
                    tier,
                    nm_loci,
                    nm_lost,
                    hgt_loci,
                    has_kp,
                    has_ap,
                    ap_paths,
                    ap_lost,
                    mro_c,
                    mro_p,
                    mtdna,
                ) = self
                    .nodes
                    .get(&id)
                    .map(|n| {
                        (
                            n.endo_tier,
                            n.nucleomorph_loci,
                            n.is_nucleomorph_loss,
                            n.klepto_hgt_loci,
                            n.has_kleptoplast,
                            n.has_apicoplast,
                            n.apicoplast_pathways,
                            n.is_apicoplast_complete_loss,
                            n.mro_class,
                            n.mro_pathways,
                            n.mro_mtdna_loci,
                        )
                    })
                    .unwrap_or((
                        u8::from(g == TrophicGuild::Holobiont),
                        0,
                        false,
                        0,
                        false,
                        false,
                        if g == TrophicGuild::Holobiont { 4 } else { 0 },
                        false,
                        if g == TrophicGuild::Holobiont {
                            MroClass::AerobicMitochondrion
                        } else {
                            MroClass::None
                        },
                        if g == TrophicGuild::Holobiont { 4 } else { 0 },
                        if g == TrophicGuild::Holobiont { 3 } else { 0 },
                    ));
                ExtantPhyloState {
                    taxon_id: id,
                    species_id: sp,
                    guild: g,
                    biome: b,
                    leverage_li: lev,
                    egt_loci,
                    organelle_rem: org_rem,
                    endo_tier: tier,
                    nucleomorph_loci: nm_loci,
                    nucleomorph_eliminated: nm_lost,
                    klepto_hgt_loci: hgt_loci,
                    has_kleptoplast: has_kp,
                    has_apicoplast: has_ap,
                    apicoplast_pathways: ap_paths,
                    apicoplast_eliminated: ap_lost,
                    mro_class: mro_c,
                    mro_pathways: mro_p,
                    mro_mtdna_loci: mtdna,
                    multicell_state: MulticellularState::Unicellular,
                    cell_count: 1,
                    soma_fraction: 0.0,
                    eusocial_state: EusocialState::Solitary,
                    colony_size: 1,
                    sterile_caste_fraction: 0.0,
                    agriculture_state: AgricultureState::HunterForager,
                    cultivar_polyploidy: 1,
                    gongylidia_yield: 0.0,
                    cultural_state: CulturalState::AsocialTrialError,
                    cultural_ratchet_tier: 0,
                    cultural_repertoire_bits: 0.0,
                    exosomatic_state: ExosomaticState::EndosomaticBound,
                    exosomatic_archive_bits: 0.0,
                    machine_synthesis_state: MachineSynthesisState::ManualArchival,
                    verified_theorem_bits: 0.0,
                    distribution_tail_diversity: 1.0,
                    planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
                    chemical_disequilibrium_bits: 14.5,
                    abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
                    raf_closure_score: 0.96,
                    chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
                    genetic_code_mp: 0.79,
                    genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
                    chromosome_capacity_kb: 1850.0,
                }
            })
            .collect();
        self.sync_extant_population_full(current_step, &states);
    }

    pub fn sync_extant_population(
        &mut self,
        current_step: usize,
        extant_units: &[(u64, usize, TrophicGuild, BiomeZone, f64)],
    ) {
        let expanded: Vec<(u64, usize, TrophicGuild, BiomeZone, f64, usize, usize)> = extant_units
            .iter()
            .map(|&(id, sp, g, b, lev)| {
                let (e_loci, o_rem) = self
                    .nodes
                    .get(&id)
                    .map(|n| (n.egt_transferred_loci, n.organelle_loci_remaining))
                    .unwrap_or((0, 0));
                (id, sp, g, b, lev, e_loci, o_rem)
            })
            .collect();
        self.sync_extant_population_with_egt(current_step, &expanded);
    }

    fn prune_extinct_and_collapse_chains(&mut self) {
        loop {
            let dead_leaves: Vec<u64> = self
                .nodes
                .values()
                .filter(|n| {
                    !n.extant
                        && n.children.is_empty()
                        && n.symbio_children.is_empty()
                        && n.taxon_id != self.root_id
                })
                .map(|n| n.taxon_id)
                .collect();
            if dead_leaves.is_empty() {
                break;
            }
            for leaf_id in dead_leaves {
                if let Some(removed) = self.nodes.remove(&leaf_id) {
                    if let Some(pid) = removed.parent_id {
                        if let Some(p_node) = self.nodes.get_mut(&pid) {
                            p_node.children.retain(|&c| c != leaf_id);
                        }
                    }
                    if let Some(did) = removed.donor_id {
                        if let Some(d_node) = self.nodes.get_mut(&did) {
                            d_node.symbio_children.retain(|&c| c != leaf_id);
                        }
                    }
                }
            }
        }

        let candidates: Vec<u64> = self
            .nodes
            .values()
            .filter(|n| {
                !n.extant
                    && n.taxon_id != self.root_id
                    && n.children.len() == 1
                    && n.symbio_children.is_empty()
                    && n.reticulation_kind.is_none()
                    && !n.is_speciation_split
                    && !n.is_egt_lockin
                    && !n.is_nucleomorph_loss
                    && !n.is_klepto_tipping_point
                    && !n.is_apicoplast_transition
                    && !n.is_apicoplast_complete_loss
                    && !n.is_hydrogenosome_transition
                    && !n.is_mitosome_transition
                    && !n.is_mro_complete_loss
                    && !n.is_clonal_multicell_transition
                    && !n.is_germ_soma_transition
                    && !n.is_eusocial_transition
                    && !n.is_superorganism_transition
                    && !n.is_agriculture_transition
                    && !n.is_obligate_crop_transition
                    && !n.is_cultural_ratchet_transition
                    && !n.is_symbolic_teaching_transition
                    && !n.is_written_archive_transition
                    && !n.is_peer_institution_transition
                    && !n.is_machine_symbiosis_transition
                    && !n.is_formal_verifier_transition
                    && !n.is_goe_transition
                    && !n.is_homeostatic_awakening_transition
                    && !n.is_raf_protocell_transition
                    && !n.is_darwinian_threshold_transition
                    && !n.is_homochiral_symmetry_transition
                    && !n.is_sgc_crystallization_transition
                    && !n.is_rnr_dna_transition
                    && !n.is_thymine_udg_transition
                    && n.guild_transition.is_none()
            })
            .map(|n| n.taxon_id)
            .collect();

        for mid_id in candidates {
            let Some(mid_node) = self.nodes.get(&mid_id).cloned() else {
                continue;
            };
            if mid_node.children.len() != 1 || !mid_node.symbio_children.is_empty() {
                continue;
            }
            let Some(pid) = mid_node.parent_id else {
                continue;
            };
            let child_id = mid_node.children[0];
            if let Some(p_node) = self.nodes.get_mut(&pid) {
                for slot in &mut p_node.children {
                    if *slot == mid_id {
                        *slot = child_id;
                    }
                }
            }
            if let Some(c_node) = self.nodes.get_mut(&child_id) {
                c_node.parent_id = Some(pid);
            }
            self.nodes.remove(&mid_id);
        }
    }

    pub fn find_extant_mrca(&self) -> u64 {
        let mut terminal_ids: Vec<u64> = self
            .nodes
            .values()
            .filter(|n| {
                n.extant
                    || (n.children.is_empty()
                        && !n.symbio_children.is_empty()
                        && n.taxon_id != self.root_id)
            })
            .map(|n| n.taxon_id)
            .collect();
        terminal_ids.sort_unstable();
        terminal_ids.dedup();

        if terminal_ids.is_empty() {
            return self.root_id;
        }
        if terminal_ids.len() == 1 {
            return terminal_ids[0];
        }

        let mut first_chain = Vec::new();
        let mut cur = Some(terminal_ids[0]);
        let mut guard = 0usize;
        while let Some(id) = cur {
            first_chain.push(id);
            cur = self.nodes.get(&id).and_then(|n| n.parent_id);
            guard += 1;
            if guard > self.nodes.len() + 2 {
                break;
            }
        }

        let mut common: BTreeSet<u64> = first_chain.iter().copied().collect();
        for &ext_id in &terminal_ids[1..] {
            let mut anc = BTreeSet::new();
            let mut c = Some(ext_id);
            let mut g = 0usize;
            while let Some(id) = c {
                anc.insert(id);
                c = self.nodes.get(&id).and_then(|n| n.parent_id);
                g += 1;
                if g > self.nodes.len() + 2 {
                    break;
                }
            }
            common.retain(|id| anc.contains(id));
        }

        first_chain
            .into_iter()
            .find(|id| common.contains(id))
            .unwrap_or(self.root_id)
    }

    pub fn compute_summary(&self, current_step: usize) -> PhyloSummary {
        let mrca_id = self.find_extant_mrca();
        let mrca_birth = self
            .nodes
            .get(&mrca_id)
            .map(|n| n.birth_step)
            .unwrap_or(current_step);
        let tmrca_steps = current_step.saturating_sub(mrca_birth);

        let mut leaf_depths = Vec::new();
        self.collect_leaf_depths(mrca_id, 0, &mut leaf_depths);
        let n_leaves = leaf_depths.len();
        let sackin_imbalance = if n_leaves <= 2 {
            0.0
        } else {
            let raw_sackin: usize = leaf_depths.iter().sum();
            let nf = n_leaves as f64;
            let min_bal = nf * nf.log2().max(1.0);
            let max_cat = (nf * (nf + 1.0)) * 0.5 - 1.0;
            if max_cat > min_bal + 1e-6 {
                (((raw_sackin as f64) - min_bal) / (max_cat - min_bal)).clamp(0.0, 1.0)
            } else {
                0.0
            }
        };

        let active_reticulations = self
            .nodes
            .values()
            .filter(|n| n.reticulation_kind.is_some() && n.donor_id.is_some())
            .count();

        let holo_extant: Vec<&PhyloNode> = self
            .nodes
            .values()
            .filter(|n| n.extant && n.guild == TrophicGuild::Holobiont)
            .collect();
        let mean_extant_egt_fraction = if holo_extant.is_empty() {
            0.0
        } else {
            holo_extant.iter().map(|n| n.egt_nuclear_fraction).sum::<f64>()
                / (holo_extant.len() as f64)
        };
        let serial_extant: Vec<usize> = holo_extant
            .iter()
            .filter(|n| n.endo_tier >= 2)
            .map(|n| n.nucleomorph_loci)
            .collect();
        let mean_extant_nucleomorph_len = if serial_extant.is_empty() {
            0.0
        } else {
            (serial_extant.iter().sum::<usize>() as f64) / (serial_extant.len() as f64)
        };

        PhyloSummary {
            mrca_id,
            mrca_birth_step: mrca_birth,
            tmrca_steps,
            active_node_count: self.nodes.len(),
            sackin_imbalance,
            total_speciation_splits: self.total_speciation_splits,
            total_guild_transitions: self.total_guild_transitions,
            prod_to_syntroph_count: self.prod_to_syntroph_count,
            to_predator_count: self.to_predator_count,
            endosymbiotic_fusions: self.endosymbiotic_fusions,
            secondary_fusions: self.secondary_fusions,
            tertiary_fusions: self.tertiary_fusions,
            klepto_fixation_events: self.klepto_fixation_events,
            klepto_hgt_events: self.klepto_hgt_events,
            apicoplast_transition_events: self.apicoplast_transition_events,
            apicoplast_loss_events: self.apicoplast_loss_events,
            hydrogenosome_transition_events: self.hydrogenosome_transition_events,
            mitosome_transition_events: self.mitosome_transition_events,
            mro_complete_loss_events: self.mro_complete_loss_events,
            clonal_multicell_transitions: self.clonal_multicell_transitions,
            germ_soma_transitions: self.germ_soma_transitions,
            eusocial_transitions: self.eusocial_transitions,
            superorganism_transitions: self.superorganism_transitions,
            agriculture_transitions: self.agriculture_transitions,
            obligate_crop_transitions: self.obligate_crop_transitions,
            cultural_ratchet_transitions: self.cultural_ratchet_transitions,
            symbolic_teaching_transitions: self.symbolic_teaching_transitions,
            cross_deme_cultural_transfers: self.cross_deme_cultural_transfers,
            written_archive_transitions: self.written_archive_transitions,
            peer_institution_transitions: self.peer_institution_transitions,
            machine_symbiosis_transitions: self.machine_symbiosis_transitions,
            formal_verifier_transitions: self.formal_verifier_transitions,
            goe_transitions: self.goe_transitions,
            homeostatic_awakening_transitions: self.homeostatic_awakening_transitions,
            raf_protocell_transitions: self.raf_protocell_transitions,
            darwinian_threshold_transitions: self.darwinian_threshold_transitions,
            homochiral_symmetry_transitions: self.homochiral_symmetry_transitions,
            sgc_crystallization_transitions: self.sgc_crystallization_transitions,
            rnr_dna_transitions: self.rnr_dna_transitions,
            thymine_udg_transitions: self.thymine_udg_transitions,
            active_reticulations,
            horizontal_transfers: self.horizontal_transfers,
            egt_transfer_events: self.egt_transfer_events,
            egt_lockin_events: self.egt_lockin_events,
            nucleomorph_loss_events: self.nucleomorph_loss_events,
            mean_extant_egt_fraction,
            mean_extant_nucleomorph_len,
        }
    }

    fn collect_leaf_depths(&self, node_id: u64, depth: usize, out: &mut Vec<usize>) {
        let Some(node) = self.nodes.get(&node_id) else {
            return;
        };
        if node.children.is_empty() {
            if node.extant {
                out.push(depth);
            }
            return;
        }
        let next_depth = if node.children.len() >= 2 {
            depth + 1
        } else {
            depth
        };
        if node.extant {
            out.push(next_depth);
        }
        for &cid in &node.children {
            self.collect_leaf_depths(cid, next_depth, out);
        }
    }

    pub fn to_newick(&self, current_step: usize, include_nhx: bool) -> String {
        let mrca = self.find_extant_mrca();
        let mut buf = String::with_capacity(3600);
        let parent_birth = self
            .nodes
            .get(&mrca)
            .map(|n| n.birth_step)
            .unwrap_or(current_step);
        self.write_newick_subtree(mrca, parent_birth, current_step, include_nhx, &mut buf);
        buf.push(';');
        buf
    }

    fn write_newick_subtree(
        &self,
        node_id: u64,
        parent_step: usize,
        current_step: usize,
        include_nhx: bool,
        out: &mut String,
    ) {
        let Some(node) = self.nodes.get(&node_id) else {
            return;
        };
        let node_time = self.effective_node_time(node, current_step);
        let branch_len = (node_time.saturating_sub(parent_step)) as f64;

        let active_symbio: Vec<u64> = node
            .symbio_children
            .iter()
            .copied()
            .filter(|cid| self.nodes.contains_key(cid))
            .collect();

        let has_group = !node.children.is_empty() || !active_symbio.is_empty();
        let mut sub_items_written = 0usize;

        let hybrid_suffix = match (node.hybrid_id, node.reticulation_kind) {
            (Some(hid), Some(rk)) if rk.is_any_endosymbiosis() => format!("#H{hid}"),
            (_, Some(ReticulationKind::HorizontalTransfer)) => {
                format!("#LGT{}", node.taxon_id)
            }
            _ => String::new(),
        };

        if has_group {
            out.push('(');
            if node.extant {
                let self_len = (current_step.saturating_sub(node_time)) as f64;
                let _ = write!(
                    out,
                    "U{}_{}_Sp{}:{:.1}",
                    node.taxon_id,
                    node.guild.short_code(),
                    node.species_id,
                    self_len
                );
                if include_nhx {
                    self.append_nhx_tag(node, out);
                }
                sub_items_written += 1;
            }
            for &cid in &node.children {
                if sub_items_written > 0 {
                    out.push(',');
                }
                self.write_newick_subtree(cid, node_time, current_step, include_nhx, out);
                sub_items_written += 1;
            }
            for &hid_target in &active_symbio {
                if let Some(h_node) = self.nodes.get(&hid_target) {
                    if sub_items_written > 0 {
                        out.push(',');
                    }
                    let sym_len = (h_node.birth_step.saturating_sub(node_time)) as f64;
                    let h_tag = match (h_node.hybrid_id, h_node.reticulation_kind) {
                        (Some(k), Some(rk)) if rk.is_any_endosymbiosis() => format!("#H{k}"),
                        _ => format!("#LGT{}", h_node.taxon_id),
                    };
                    let _ = write!(out, "{h_tag}:{sym_len:.1}");
                    if include_nhx {
                        let donor_role = match h_node.reticulation_kind {
                            Some(ReticulationKind::KleptoTippingFixation) => {
                                "Kleptoplast_Fixed_Donor"
                            }
                            Some(ReticulationKind::SecondaryEndosymbiosis) => {
                                "Secondary_Eukaryote_Donor"
                            }
                            Some(ReticulationKind::TertiaryEndosymbiosis) => {
                                "Tertiary_Complex_Donor"
                            }
                            _ => "Endosymbiont_Donor",
                        };
                        let _ = write!(
                            out,
                            "[&&NHX:ROLE={}:TIER={}:HOST=N{}:GAMMA={:.2}]",
                            donor_role,
                            h_node.endo_tier.max(1),
                            h_node.parent_id.unwrap_or(0),
                            h_node.genome_share_gamma
                        );
                    }
                    sub_items_written += 1;
                }
            }
            out.push(')');
            let _ = write!(out, "N{}{}:{:.1}", node.taxon_id, hybrid_suffix, branch_len);
        } else {
            let end_step = if node.extant {
                current_step
            } else {
                node.death_step.unwrap_or(current_step)
            };
            let tip_len = (end_step.saturating_sub(parent_step)) as f64;
            let _ = write!(
                out,
                "U{}_{}_Sp{}{}:{:.1}",
                node.taxon_id,
                node.guild.short_code(),
                node.species_id,
                hybrid_suffix,
                tip_len
            );
        }

        if include_nhx {
            self.append_nhx_tag(node, out);
        }
    }

    fn effective_node_time(&self, node: &PhyloNode, current_step: usize) -> usize {
        if node.extant && node.children.is_empty() && node.symbio_children.is_empty() {
            current_step
        } else {
            let earliest_vert = node
                .children
                .iter()
                .filter_map(|c| self.nodes.get(c).map(|cn| cn.birth_step))
                .min()
                .unwrap_or(current_step);
            let earliest_sym = node
                .symbio_children
                .iter()
                .filter_map(|c| self.nodes.get(c).map(|cn| cn.birth_step))
                .min()
                .unwrap_or(current_step);
            node.death_step
                .unwrap_or(current_step)
                .min(earliest_vert)
                .min(earliest_sym)
                .max(node.birth_step)
        }
    }

    fn append_nhx_tag(&self, node: &PhyloNode, out: &mut String) {
        let biome_str = match node.biome {
            BiomeZone::Alpha => "Alpha",
            BiomeZone::Beta => "Beta",
        };
        let spec_flag = u8::from(node.is_speciation_split);
        let shift_str = match node.guild_transition {
            Some((f, t)) => format!("{}_to_{}", f.short_code(), t.short_code()),
            None => "none".to_string(),
        };
        let mut extra = match (node.hybrid_id, node.reticulation_kind) {
            (Some(k), Some(ReticulationKind::EndosymbioticFusion)) => format!(
                ":HYBRID=H{}:FUSION=Prod+Synt->Holo:TIER=1:GAMMA={:.2}",
                k, node.genome_share_gamma
            ),
            (Some(k), Some(ReticulationKind::SecondaryEndosymbiosis)) => format!(
                ":HYBRID=H{}:FUSION=Secondary_Euk+Holo->2degHolo:TIER=2:GAMMA={:.2}",
                k, node.genome_share_gamma
            ),
            (Some(k), Some(ReticulationKind::TertiaryEndosymbiosis)) => format!(
                ":HYBRID=H{}:FUSION=Tertiary_Holo+2degHolo->3degHolo:TIER=3:GAMMA={:.2}",
                k, node.genome_share_gamma
            ),
            (Some(k), Some(ReticulationKind::KleptoTippingFixation)) => format!(
                ":HYBRID=H{}:FUSION=Klepto_Tipping_Fixation:TIER={}:KLEPTO_TIP=1:HGT={}:GAMMA={:.2}",
                k,
                node.endo_tier.max(1),
                node.klepto_hgt_loci,
                node.genome_share_gamma
            ),
            (_, Some(ReticulationKind::HorizontalTransfer)) => {
                format!(":LGT=1:GAMMA={:.2}", node.genome_share_gamma)
            }
            _ => String::new(),
        };
        if node.has_kleptoplast
            || (node.klepto_hgt_loci > 0 && node.guild == TrophicGuild::PredatorParasite)
        {
            let _ = write!(
                extra,
                ":KLEPTO={}:HGT={}",
                u8::from(node.has_kleptoplast),
                node.klepto_hgt_loci
            );
        }
        if node.has_apicoplast
            || node.is_apicoplast_transition
            || node.is_apicoplast_complete_loss
        {
            let _ = write!(
                extra,
                ":APICO={}:APICO_PATHS={}:APICO_TRANS={}:APICO_LOST={}",
                u8::from(node.has_apicoplast),
                node.apicoplast_pathways,
                u8::from(node.is_apicoplast_transition),
                u8::from(node.is_apicoplast_complete_loss)
            );
        }
        if node.mro_class != MroClass::None
            || node.is_hydrogenosome_transition
            || node.is_mitosome_transition
            || node.is_mro_complete_loss
        {
            let _ = write!(
                extra,
                ":MRO={}:MRO_PATHS={}:MTDNA={}:HYDRO_TR={}:MITOS_TR={}:MRO_LOST={}",
                node.mro_class.short_code(),
                node.mro_pathways,
                node.mro_mtdna_loci,
                u8::from(node.is_hydrogenosome_transition),
                u8::from(node.is_mitosome_transition),
                u8::from(node.is_mro_complete_loss)
            );
        }
        if node.genome_substrate_stage != GenomeSubstrateStage::CoevolvedCellularVirosphere
            || node.is_rnr_dna_transition
            || node.is_thymine_udg_transition
        {
            let _ = write!(
                extra,
                ":DNA={}:KB={:.0}:RNR_TR={}:UDG_TR={}",
                node.genome_substrate_stage.short_code(),
                node.chromosome_capacity_kb,
                u8::from(node.is_rnr_dna_transition),
                u8::from(node.is_thymine_udg_transition)
            );
        }
        if node.chirality_code_stage != ChiralityCodeStage::UniversalErrorMinimizedSgc
            || node.is_homochiral_symmetry_transition
            || node.is_sgc_crystallization_transition
        {
            let _ = write!(
                extra,
                ":CHIRAL={}:MP={:.2}:HOMO_TR={}:SGC_TR={}",
                node.chirality_code_stage.short_code(),
                node.genetic_code_mp,
                u8::from(node.is_homochiral_symmetry_transition),
                u8::from(node.is_sgc_crystallization_transition)
            );
        }
        if node.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage
            || node.is_raf_protocell_transition
            || node.is_darwinian_threshold_transition
        {
            let _ = write!(
                extra,
                ":ABIO={}:RAF={:.2}:RAF_TR={}:DARWIN_TR={}",
                node.abiogenesis_stage.short_code(),
                node.raf_closure_score,
                u8::from(node.is_raf_protocell_transition),
                u8::from(node.is_darwinian_threshold_transition)
            );
        }
        if node.planetary_epoch != PlanetaryState::OxygenatedPhanerozoic
            || node.is_goe_transition
            || node.is_homeostatic_awakening_transition
        {
            let _ = write!(
                extra,
                ":GAIA={}:DISEQ={:.1}:GOE_TR={}:GAIA2_TR={}",
                node.planetary_epoch.short_code(),
                node.chemical_disequilibrium_bits,
                u8::from(node.is_goe_transition),
                u8::from(node.is_homeostatic_awakening_transition)
            );
        }
        if node.machine_synthesis_state.has_computational_symbiont()
            || node.is_machine_symbiosis_transition
            || node.is_formal_verifier_transition
        {
            let _ = write!(
                extra,
                ":MSYN={}:THMBITS={:.1}:TAILS={:.2}:MSYN_TR={}:FORMAL_TR={}",
                node.machine_synthesis_state.short_code(),
                node.verified_theorem_bits,
                node.distribution_tail_diversity,
                u8::from(node.is_machine_symbiosis_transition),
                u8::from(node.is_formal_verifier_transition)
            );
        }
        if node.exosomatic_state.is_archival()
            || node.is_written_archive_transition
            || node.is_peer_institution_transition
        {
            let _ = write!(
                extra,
                ":EXO={}:ABITS={:.1}:ARCH_TR={}:INST_TR={}",
                node.exosomatic_state.short_code(),
                node.exosomatic_archive_bits,
                u8::from(node.is_written_archive_transition),
                u8::from(node.is_peer_institution_transition)
            );
        }
        if node.cultural_state.is_cultural()
            || node.is_cultural_ratchet_transition
            || node.is_symbolic_teaching_transition
        {
            let _ = write!(
                extra,
                ":CULT={}:RATCHET=T{}:CBITS={:.1}:RATCHET_TR={}:SYM_TR={}",
                node.cultural_state.short_code(),
                node.cultural_ratchet_tier,
                node.cultural_repertoire_bits,
                u8::from(node.is_cultural_ratchet_transition),
                u8::from(node.is_symbolic_teaching_transition)
            );
        }
        if node.agriculture_state.is_farming()
            || node.is_agriculture_transition
            || node.is_obligate_crop_transition
        {
            let _ = write!(
                extra,
                ":AGRI={}:PLOIDY={}n:GONGYL={:.2}:AGRI_TR={}:CROP_TR={}",
                node.agriculture_state.short_code(),
                node.cultivar_polyploidy.max(1),
                node.gongylidia_yield,
                u8::from(node.is_agriculture_transition),
                u8::from(node.is_obligate_crop_transition)
            );
        }
        if node.eusocial_state.is_social()
            || node.is_eusocial_transition
            || node.is_superorganism_transition
        {
            let _ = write!(
                extra,
                ":EUSOC={}:COLONY={}:CASTES={:.2}:EUSOC_TR={}:SUPER_TR={}",
                node.eusocial_state.short_code(),
                node.colony_size,
                node.sterile_caste_fraction,
                u8::from(node.is_eusocial_transition),
                u8::from(node.is_superorganism_transition)
            );
        }
        if node.multicell_state.is_multicellular()
            || node.is_clonal_multicell_transition
            || node.is_germ_soma_transition
        {
            let _ = write!(
                extra,
                ":MC={}:CELLS={}:SOMA={:.2}:CLONAL_TR={}:GS_TR={}",
                node.multicell_state.short_code(),
                node.cell_count,
                node.soma_fraction,
                u8::from(node.is_clonal_multicell_transition),
                u8::from(node.is_germ_soma_transition)
            );
        }
        if node.guild == TrophicGuild::Holobiont
            || node.egt_transferred_loci > 0
            || node.organelle_loci_remaining > 0
        {
            let _ = write!(
                extra,
                ":TIER_CUR={}:EGT={}:ORG_LEN={}:NUC_INT={:.2}:LOCK={}:NM_LEN={}:NM_LOST={}",
                node.endo_tier.max(1),
                node.egt_transferred_loci,
                node.organelle_loci_remaining,
                node.egt_nuclear_fraction,
                u8::from(node.is_egt_lockin),
                node.nucleomorph_loci,
                u8::from(node.is_nucleomorph_loss)
            );
        }
        let _ = write!(
            out,
            "[&&NHX:G={}:S={}:B={}:SPEC={}:SHIFT={}{}:L={:.2}]",
            node.guild.short_code(),
            node.species_id,
            biome_str,
            spec_flag,
            shift_str,
            extra,
            node.leverage_li
        );
    }

    pub fn extract_layout(&self, current_step: usize) -> PhylogramLayout {
        let mrca = self.find_extant_mrca();
        let min_step = self
            .nodes
            .get(&mrca)
            .map(|n| n.birth_step)
            .unwrap_or(0);
        let max_step = current_step.max(min_step + 1);
        let span_t = (max_step - min_step).max(1) as f64;

        let mut leaf_counter = 0usize;
        let total_leaves = self.count_visual_leaves(mrca).max(1);
        let mut y_coords = BTreeMap::new();
        let mut x_coords = BTreeMap::new();
        let mut vis_nodes = Vec::new();
        let mut vis_branches = Vec::new();

        self.layout_dfs(
            mrca,
            min_step,
            max_step,
            span_t,
            total_leaves,
            &mut leaf_counter,
            &mut x_coords,
            &mut y_coords,
            &mut vis_nodes,
            &mut vis_branches,
        );

        let mut reticulations = Vec::new();
        for n in self.nodes.values() {
            if let Some(did) = n.donor_id {
                if let (Some(&fx), Some(&fy), Some(&tx), Some(&ty)) = (
                    x_coords.get(&did),
                    y_coords.get(&did),
                    x_coords.get(&n.taxon_id),
                    y_coords.get(&n.taxon_id),
                ) {
                    let rk = n.reticulation_kind;
                    reticulations.push(PhylogramReticulationVis {
                        from_x: fx,
                        from_y: fy,
                        to_x: tx,
                        to_y: ty,
                        is_endosymbiosis: rk.map(|k| k.is_any_endosymbiosis()).unwrap_or(false),
                        is_klepto_tipping: rk == Some(ReticulationKind::KleptoTippingFixation),
                        endo_tier: n.endo_tier.max(rk.map(|k| k.endo_tier()).unwrap_or(0)),
                        hybrid_id: n.hybrid_id,
                        gamma: n.genome_share_gamma,
                    });
                }
            }
        }

        PhylogramLayout {
            min_step,
            max_step,
            mrca_id: mrca,
            nodes: vis_nodes,
            branches: vis_branches,
            reticulations,
        }
    }

    fn is_symbio_donor_terminal(&self, node: &PhyloNode) -> bool {
        !node.extant
            && node.children.is_empty()
            && node
                .symbio_children
                .iter()
                .any(|c| self.nodes.contains_key(c))
    }

    fn count_visual_leaves(&self, node_id: u64) -> usize {
        let Some(node) = self.nodes.get(&node_id) else {
            return 0;
        };
        if node.children.is_empty() {
            return usize::from(node.extant || self.is_symbio_donor_terminal(node));
        }
        let mut sum = usize::from(node.extant);
        for &cid in &node.children {
            sum += self.count_visual_leaves(cid);
        }
        sum
    }

    #[allow(clippy::too_many_arguments)]
    fn layout_dfs(
        &self,
        node_id: u64,
        min_step: usize,
        max_step: usize,
        span_t: f64,
        total_leaves: usize,
        leaf_counter: &mut usize,
        x_coords: &mut BTreeMap<u64, f64>,
        y_coords: &mut BTreeMap<u64, f64>,
        vis_nodes: &mut Vec<PhylogramNodeVis>,
        vis_branches: &mut Vec<PhylogramBranchVis>,
    ) -> Option<f64> {
        let node = self.nodes.get(&node_id)?;
        let donor_terminal = self.is_symbio_donor_terminal(node);
        let node_t = if node.children.is_empty() && node.extant {
            max_step
        } else if donor_terminal {
            node.symbio_children
                .iter()
                .filter_map(|c| self.nodes.get(c).map(|cn| cn.birth_step))
                .max()
                .unwrap_or(node.death_step.unwrap_or(max_step))
                .max(node.birth_step)
                .min(max_step)
        } else {
            self.effective_node_time(node, max_step).min(max_step)
        };
        let x_norm = ((node_t.saturating_sub(min_step)) as f64 / span_t).clamp(0.0, 1.0);

        let format_tip = |nd: &PhyloNode| -> String {
            if nd.planetary_epoch == PlanetaryState::HomeostaticAwakening {
                return format!(
                    "U{} Gaia2Homeostasis({:.1}bDiseq)",
                    nd.taxon_id, nd.chemical_disequilibrium_bits
                );
            }
            if nd.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis {
                return format!(
                    "U{} FormalSymbio({:.0}b,T{:.0}%)",
                    nd.taxon_id,
                    nd.verified_theorem_bits,
                    nd.distribution_tail_diversity * 100.0
                );
            } else if nd.machine_synthesis_state == MachineSynthesisState::AlgorithmicMediation {
                return format!("U{} AlgoMed({:.0}b)", nd.taxon_id, nd.verified_theorem_bits);
            }
            if nd.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                return format!(
                    "U{} VerifiedInst({:.0}b)",
                    nd.taxon_id, nd.exosomatic_archive_bits
                );
            } else if nd.exosomatic_state == ExosomaticState::WrittenArchive {
                return format!(
                    "U{} WrittenArch({:.0}b)",
                    nd.taxon_id, nd.exosomatic_archive_bits
                );
            }
            if nd.cultural_state == CulturalState::SymbolicCumulativeCulture {
                return format!(
                    "U{} SymbolicCCE(T{},{:.0}b)",
                    nd.taxon_id,
                    nd.cultural_ratchet_tier,
                    nd.cultural_repertoire_bits
                );
            } else if nd.cultural_state == CulturalState::ImitativeRatchet {
                return format!(
                    "U{} CultRatchet(T{},{:.0}b)",
                    nd.taxon_id,
                    nd.cultural_ratchet_tier,
                    nd.cultural_repertoire_bits
                );
            }
            if nd.agriculture_state == AgricultureState::LeafCutterAgroEcosystem {
                return format!(
                    "U{} LeafCutter({}n,{:.0}%G)",
                    nd.taxon_id,
                    nd.cultivar_polyploidy,
                    nd.gongylidia_yield * 100.0
                );
            } else if nd.agriculture_state == AgricultureState::HigherDomesticatedCrop {
                return format!(
                    "U{} GongylCrop({}n,{:.0}%G)",
                    nd.taxon_id,
                    nd.cultivar_polyploidy,
                    nd.gongylidia_yield * 100.0
                );
            } else if nd.agriculture_state == AgricultureState::LowerFacultativeFarming {
                return format!("U{} LowerAg({}n)", nd.taxon_id, nd.cultivar_polyploidy);
            }
            if nd.eusocial_state == EusocialState::ObligateSuperorganism {
                return format!(
                    "U{} Superorg({}col,{:.0}%C)",
                    nd.taxon_id,
                    nd.colony_size,
                    nd.sterile_caste_fraction * 100.0
                );
            } else if nd.eusocial_state == EusocialState::FacultativeEusocial {
                return format!("U{} FacEusoc({}col)", nd.taxon_id, nd.colony_size);
            }
            if nd.multicell_state == MulticellularState::GermSomaMulticellular {
                return format!(
                    "U{} GermSoma({}c,{:.0}%S)",
                    nd.taxon_id,
                    nd.cell_count,
                    nd.soma_fraction * 100.0
                );
            } else if nd.multicell_state == MulticellularState::UndifferentiatedClonal {
                return format!("U{} Clonal({}c)", nd.taxon_id, nd.cell_count);
            }
            match nd.mro_class {
                MroClass::Hydrogenosome => {
                    return format!("U{} Hydro(H2+ISC)", nd.taxon_id);
                }
                MroClass::Mitosome => {
                    return format!("U{} MitoS(ISC=1)", nd.taxon_id);
                }
                MroClass::Amiteukaryote => {
                    return format!("U{} Amiteuk(SUF)", nd.taxon_id);
                }
                _ => {}
            }
            if nd.guild == TrophicGuild::Holobiont {
                if nd.endo_tier >= 2 {
                    format!("U{} {}°Holo NM:{}", nd.taxon_id, nd.endo_tier, nd.nucleomorph_loci)
                } else {
                    format!(
                        "U{} 1°Holo EGT:{:.0}%",
                        nd.taxon_id,
                        nd.egt_nuclear_fraction * 100.0
                    )
                }
            } else if nd.has_apicoplast {
                format!("U{} Apico(P:{})", nd.taxon_id, nd.apicoplast_pathways)
            } else if nd.is_apicoplast_complete_loss {
                format!("U{} Crypto(Ap0)", nd.taxon_id)
            } else if nd.has_kleptoplast {
                format!("U{} Klepto HGT:{}", nd.taxon_id, nd.klepto_hgt_loci)
            } else {
                format!("U{} {} Sp{}", nd.taxon_id, nd.guild.short_code(), nd.species_id)
            }
        };

        let y_norm = if node.children.is_empty() {
            if !node.extant && !donor_terminal {
                return None;
            }
            let yn = if total_leaves <= 1 {
                0.5
            } else {
                (*leaf_counter as f64) / ((total_leaves - 1) as f64)
            };
            *leaf_counter += 1;
            yn
        } else {
            let mut child_ys = Vec::new();
            let mut child_records = Vec::new();

            if node.extant {
                let bud_y = if total_leaves <= 1 {
                    0.5
                } else {
                    (*leaf_counter as f64) / ((total_leaves - 1) as f64)
                };
                *leaf_counter += 1;
                child_ys.push(bud_y);
                vis_nodes.push(PhylogramNodeVis {
                    taxon_id: node.taxon_id,
                    x_norm: 1.0,
                    y_norm: bud_y,
                    extant: true,
                    species_id: node.species_id,
                    guild: node.guild,
                    biome: node.biome,
                    is_speciation_split: false,
                    guild_transition: None,
                    hybrid_id: node.hybrid_id,
                    endo_tier: node.endo_tier,
                    nucleomorph_loci: node.nucleomorph_loci,
                    is_nucleomorph_loss: false,
                    has_kleptoplast: node.has_kleptoplast,
                    is_klepto_tipping_point: false,
                    has_apicoplast: node.has_apicoplast,
                    apicoplast_pathways: node.apicoplast_pathways,
                    is_apicoplast_transition: false,
                    is_apicoplast_complete_loss: false,
                    mro_class: node.mro_class,
                    mro_pathways: node.mro_pathways,
                    is_hydrogenosome_transition: false,
                    is_mitosome_transition: false,
                    is_mro_complete_loss: false,
                    multicell_state: node.multicell_state,
                    cell_count: node.cell_count,
                    is_clonal_multicell_transition: false,
                    is_germ_soma_transition: false,
                    eusocial_state: node.eusocial_state,
                    colony_size: node.colony_size,
                    is_eusocial_transition: false,
                    is_superorganism_transition: false,
                    agriculture_state: node.agriculture_state,
                    cultivar_polyploidy: node.cultivar_polyploidy,
                    is_agriculture_transition: false,
                    is_obligate_crop_transition: false,
                    cultural_state: node.cultural_state,
                    cultural_ratchet_tier: node.cultural_ratchet_tier,
                    is_cultural_ratchet_transition: false,
                    is_symbolic_teaching_transition: false,
                    exosomatic_state: node.exosomatic_state,
                    is_written_archive_transition: false,
                    is_peer_institution_transition: false,
                    machine_synthesis_state: node.machine_synthesis_state,
                    is_machine_symbiosis_transition: false,
                    is_formal_verifier_transition: false,
                    planetary_epoch: node.planetary_epoch,
                    is_goe_transition: false,
                    is_homeostatic_awakening_transition: false,
                    abiogenesis_stage: node.abiogenesis_stage,
                    is_raf_protocell_transition: false,
                    is_darwinian_threshold_transition: false,
                    chirality_code_stage: node.chirality_code_stage,
                    is_homochiral_symmetry_transition: false,
                    is_sgc_crystallization_transition: false,
                    genome_substrate_stage: node.genome_substrate_stage,
                    is_rnr_dna_transition: false,
                    is_thymine_udg_transition: false,
                    egt_nuclear_fraction: node.egt_nuclear_fraction,
                    is_egt_lockin: false,
                    label: format_tip(node),
                });
                child_records.push(PhylogramBranchVis {
                    parent_x: x_norm,
                    parent_y: 0.0,
                    child_x: 1.0,
                    child_y: bud_y,
                    child_guild: node.guild,
                    child_species: node.species_id,
                    child_endo_tier: node.endo_tier,
                    child_has_apicoplast: node.has_apicoplast,
                    child_mro_class: node.mro_class,
                    is_speciation_split: false,
                    guild_transition: None,
                    child_egt_fraction: node.egt_nuclear_fraction,
                    is_egt_lockin: false,
                    is_nucleomorph_loss: false,
                    is_klepto_tipping_point: false,
                    is_apicoplast_transition: false,
                    is_apicoplast_complete_loss: false,
                    is_hydrogenosome_transition: false,
                    is_mitosome_transition: false,
                    is_mro_complete_loss: false,
                    child_multicell_state: node.multicell_state,
                    is_clonal_multicell_transition: false,
                    is_germ_soma_transition: false,
                    child_eusocial_state: node.eusocial_state,
                    is_eusocial_transition: false,
                    is_superorganism_transition: false,
                    child_agriculture_state: node.agriculture_state,
                    is_agriculture_transition: false,
                    is_obligate_crop_transition: false,
                    child_cultural_state: node.cultural_state,
                    is_cultural_ratchet_transition: false,
                    is_symbolic_teaching_transition: false,
                    child_exosomatic_state: node.exosomatic_state,
                    is_written_archive_transition: false,
                    is_peer_institution_transition: false,
                    child_machine_synthesis_state: node.machine_synthesis_state,
                    is_machine_symbiosis_transition: false,
                    is_formal_verifier_transition: false,
                    child_planetary_epoch: node.planetary_epoch,
                    is_goe_transition: false,
                    is_homeostatic_awakening_transition: false,
                    child_abiogenesis_stage: node.abiogenesis_stage,
                    is_raf_protocell_transition: false,
                    is_darwinian_threshold_transition: false,
                    child_chirality_code_stage: node.chirality_code_stage,
                    is_homochiral_symmetry_transition: false,
                    is_sgc_crystallization_transition: false,
                    child_genome_substrate_stage: node.genome_substrate_stage,
                    is_rnr_dna_transition: false,
                    is_thymine_udg_transition: false,
                });
            }

            for &cid in &node.children {
                if let Some(cy) = self.layout_dfs(
                    cid,
                    min_step,
                    max_step,
                    span_t,
                    total_leaves,
                    leaf_counter,
                    x_coords,
                    y_coords,
                    vis_nodes,
                    vis_branches,
                ) {
                    if let Some(cn) = self.nodes.get(&cid) {
                        let cx = x_coords.get(&cid).copied().unwrap_or(x_norm);
                        child_ys.push(cy);
                        child_records.push(PhylogramBranchVis {
                            parent_x: x_norm,
                            parent_y: 0.0,
                            child_x: cx,
                            child_y: cy,
                            child_guild: cn.guild,
                            child_species: cn.species_id,
                            child_endo_tier: cn.endo_tier,
                            child_has_apicoplast: cn.has_apicoplast,
                            child_mro_class: cn.mro_class,
                            is_speciation_split: cn.is_speciation_split,
                            guild_transition: cn.guild_transition,
                            child_egt_fraction: cn.egt_nuclear_fraction,
                            is_egt_lockin: cn.is_egt_lockin,
                            is_nucleomorph_loss: cn.is_nucleomorph_loss,
                            is_klepto_tipping_point: cn.is_klepto_tipping_point,
                            is_apicoplast_transition: cn.is_apicoplast_transition,
                            is_apicoplast_complete_loss: cn.is_apicoplast_complete_loss,
                            is_hydrogenosome_transition: cn.is_hydrogenosome_transition,
                            is_mitosome_transition: cn.is_mitosome_transition,
                            is_mro_complete_loss: cn.is_mro_complete_loss,
                            child_multicell_state: cn.multicell_state,
                            is_clonal_multicell_transition: cn.is_clonal_multicell_transition,
                            is_germ_soma_transition: cn.is_germ_soma_transition,
                            child_eusocial_state: cn.eusocial_state,
                            is_eusocial_transition: cn.is_eusocial_transition,
                            is_superorganism_transition: cn.is_superorganism_transition,
                            child_agriculture_state: cn.agriculture_state,
                            is_agriculture_transition: cn.is_agriculture_transition,
                            is_obligate_crop_transition: cn.is_obligate_crop_transition,
                            child_cultural_state: cn.cultural_state,
                            is_cultural_ratchet_transition: cn.is_cultural_ratchet_transition,
                            is_symbolic_teaching_transition: cn.is_symbolic_teaching_transition,
                            child_exosomatic_state: cn.exosomatic_state,
                            is_written_archive_transition: cn.is_written_archive_transition,
                            is_peer_institution_transition: cn.is_peer_institution_transition,
                            child_machine_synthesis_state: cn.machine_synthesis_state,
                            is_machine_symbiosis_transition: cn.is_machine_symbiosis_transition,
                            is_formal_verifier_transition: cn.is_formal_verifier_transition,
                            child_planetary_epoch: cn.planetary_epoch,
                            is_goe_transition: cn.is_goe_transition,
                            is_homeostatic_awakening_transition: cn.is_homeostatic_awakening_transition,
                            child_abiogenesis_stage: cn.abiogenesis_stage,
                            is_raf_protocell_transition: cn.is_raf_protocell_transition,
                            is_darwinian_threshold_transition: cn.is_darwinian_threshold_transition,
                            child_chirality_code_stage: cn.chirality_code_stage,
                            is_homochiral_symmetry_transition: cn.is_homochiral_symmetry_transition,
                            is_sgc_crystallization_transition: cn.is_sgc_crystallization_transition,
                            child_genome_substrate_stage: cn.genome_substrate_stage,
                            is_rnr_dna_transition: cn.is_rnr_dna_transition,
                            is_thymine_udg_transition: cn.is_thymine_udg_transition,
                        });
                    }
                }
            }

            if child_ys.is_empty() {
                return None;
            }
            let mean_y = child_ys.iter().sum::<f64>() / (child_ys.len() as f64);
            for mut br in child_records {
                br.parent_y = mean_y;
                vis_branches.push(br);
            }
            mean_y
        };

        x_coords.insert(node.taxon_id, x_norm);
        y_coords.insert(node.taxon_id, y_norm);

        vis_nodes.push(PhylogramNodeVis {
            taxon_id: node.taxon_id,
            x_norm,
            y_norm,
            extant: node.extant && node.children.is_empty(),
            species_id: node.species_id,
            guild: node.guild,
            biome: node.biome,
            is_speciation_split: node.is_speciation_split,
            guild_transition: node.guild_transition,
            hybrid_id: node.hybrid_id,
            endo_tier: node.endo_tier,
            nucleomorph_loci: node.nucleomorph_loci,
            is_nucleomorph_loss: node.is_nucleomorph_loss,
            has_kleptoplast: node.has_kleptoplast,
            is_klepto_tipping_point: node.is_klepto_tipping_point,
            has_apicoplast: node.has_apicoplast,
            apicoplast_pathways: node.apicoplast_pathways,
            is_apicoplast_transition: node.is_apicoplast_transition,
            is_apicoplast_complete_loss: node.is_apicoplast_complete_loss,
            mro_class: node.mro_class,
            mro_pathways: node.mro_pathways,
            is_hydrogenosome_transition: node.is_hydrogenosome_transition,
            is_mitosome_transition: node.is_mitosome_transition,
            is_mro_complete_loss: node.is_mro_complete_loss,
            multicell_state: node.multicell_state,
            cell_count: node.cell_count,
            is_clonal_multicell_transition: node.is_clonal_multicell_transition,
            is_germ_soma_transition: node.is_germ_soma_transition,
            eusocial_state: node.eusocial_state,
            colony_size: node.colony_size,
            is_eusocial_transition: node.is_eusocial_transition,
            is_superorganism_transition: node.is_superorganism_transition,
            agriculture_state: node.agriculture_state,
            cultivar_polyploidy: node.cultivar_polyploidy,
            is_agriculture_transition: node.is_agriculture_transition,
            is_obligate_crop_transition: node.is_obligate_crop_transition,
            cultural_state: node.cultural_state,
            cultural_ratchet_tier: node.cultural_ratchet_tier,
            is_cultural_ratchet_transition: node.is_cultural_ratchet_transition,
            is_symbolic_teaching_transition: node.is_symbolic_teaching_transition,
            exosomatic_state: node.exosomatic_state,
            is_written_archive_transition: node.is_written_archive_transition,
            is_peer_institution_transition: node.is_peer_institution_transition,
            machine_synthesis_state: node.machine_synthesis_state,
            is_machine_symbiosis_transition: node.is_machine_symbiosis_transition,
            is_formal_verifier_transition: node.is_formal_verifier_transition,
            planetary_epoch: node.planetary_epoch,
            is_goe_transition: node.is_goe_transition,
            is_homeostatic_awakening_transition: node.is_homeostatic_awakening_transition,
            abiogenesis_stage: node.abiogenesis_stage,
            is_raf_protocell_transition: node.is_raf_protocell_transition,
            is_darwinian_threshold_transition: node.is_darwinian_threshold_transition,
            chirality_code_stage: node.chirality_code_stage,
            is_homochiral_symmetry_transition: node.is_homochiral_symmetry_transition,
            is_sgc_crystallization_transition: node.is_sgc_crystallization_transition,
            genome_substrate_stage: node.genome_substrate_stage,
            is_rnr_dna_transition: node.is_rnr_dna_transition,
            is_thymine_udg_transition: node.is_thymine_udg_transition,
            egt_nuclear_fraction: node.egt_nuclear_fraction,
            is_egt_lockin: node.is_egt_lockin,
            label: if node.children.is_empty() && node.extant {
                format_tip(node)
            } else if let Some(k) = node.hybrid_id {
                if node.is_klepto_tipping_point {
                    format!("N{}#H{}(KP→FIX)", node.taxon_id, k)
                } else {
                    format!("N{}#H{}({}°)", node.taxon_id, k, node.endo_tier.max(1))
                }
            } else {
                format!("N{}", node.taxon_id)
            },
        });

        Some(y_norm)
    }

    pub fn to_svg_phylogram(&self, title: &str, current_step: usize) -> String {
        let layout = self.extract_layout(current_step);
        let summary = self.compute_summary(current_step);
        let w = 1240.0_f64;
        let leaf_n = layout.nodes.iter().filter(|n| n.extant).count().max(8);
        let h = (210.0 + (leaf_n as f64) * 23.0).clamp(515.0, 1480.0);
        let pad_l = 74.0;
        let pad_r = 264.0;
        let pad_t = 98.0;
        let pad_b = 58.0;
        let plot_w = w - pad_l - pad_r;
        let plot_h = h - pad_t - pad_b;

        let guild_hex = |g: TrophicGuild| -> &'static str {
            match g {
                TrophicGuild::Producer => "#4cc9f0",
                TrophicGuild::Syntroph => "#50e389",
                TrophicGuild::PredatorParasite => "#ff5c8a",
                TrophicGuild::Holobiont => "#c77dff",
            }
        };

        let mut svg = String::with_capacity(48_000);
        let _ = writeln!(
            svg,
            r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:.0} {h:.0}" width="{w:.0}" height="{h:.0}">"#
        );
        let _ = writeln!(svg, r##"<rect width="100%" height="100%" fill="#0b1020"/>"##);
        let _ = writeln!(
            svg,
            r##"<text x="36" y="36" fill="#f0f4ff" font-family="monospace, sans-serif" font-size="18" font-weight="bold">{}</text>"##,
            escape_svg_xml(title)
        );
        let _ = writeln!(
            svg,
            r##"<text x="36" y="60" fill="#9fb6e0" font-family="monospace, sans-serif" font-size="12">MRCA: N{} | #H: {} | StatCode→UniversalSGC: {} | RNA→U-DNA(RNR): {} | U-DNA→T-DNA(UDG): {} | Techno→Gaia2Homeostasis: {}</text>"##,
            summary.mrca_id,
            summary.endosymbiotic_fusions,
            summary.sgc_crystallization_transitions,
            summary.rnr_dna_transitions,
            summary.thymine_udg_transitions,
            summary.homeostatic_awakening_transitions
        );

        let axis_y = h - 28.0;
        let _ = writeln!(
            svg,
            r##"<line x1="{pad_l:.1}" y1="{axis_y:.1}" x2="{:.1}" y2="{axis_y:.1}" stroke="#3a4f7a" stroke-width="1.2"/>"##,
            pad_l + plot_w
        );
        for k in 0..=4 {
            let f = (k as f64) / 4.0;
            let tx = pad_l + f * plot_w;
            let step_lbl =
                layout.min_step + (((layout.max_step - layout.min_step) as f64) * f).round() as usize;
            let _ = writeln!(
                svg,
                r##"<line x1="{tx:.1}" y1="{pad_t:.1}" x2="{tx:.1}" y2="{axis_y:.1}" stroke="#1c2744" stroke-dasharray="3,4" stroke-width="1"/>"##
            );
            let _ = writeln!(
                svg,
                r##"<text x="{tx:.1}" y="{:.1}" text-anchor="middle" fill="#8ca4d0" font-family="monospace" font-size="11">t={step_lbl}</text>"##,
                axis_y + 18.0
            );
        }

        for br in &layout.branches {
            let x1 = pad_l + br.parent_x * plot_w;
            let y1 = pad_t + br.parent_y * plot_h;
            let x2 = pad_l + br.child_x * plot_w;
            let y2 = pad_t + br.child_y * plot_h;
            let col = match br.child_mro_class {
                MroClass::Hydrogenosome => "#06d6a0",
                MroClass::Mitosome => "#ffd166",
                MroClass::Amiteukaryote => "#38bdf8",
                _ if br.child_guild == TrophicGuild::Holobiont => match br.child_endo_tier {
                    3..=u8::MAX => "#f72585",
                    2 => "#ff9f1c",
                    _ if br.child_egt_fraction >= 0.50 => "#9d8df1",
                    _ => "#c77dff",
                },
                _ if br.child_has_apicoplast => "#e066ff",
                _ => guild_hex(br.child_guild),
            };
            let sw = if br.child_mro_class.is_reduced_anaerobic_mro()
                || br.child_endo_tier >= 2
                || br.child_has_apicoplast
            {
                2.9
            } else if br.child_guild == TrophicGuild::Holobiont {
                2.6
            } else {
                2.1
            };
            let _ = writeln!(
                svg,
                r#"<path d="M {x1:.1} {y1:.1} L {x1:.1} {y2:.1} L {x2:.1} {y2:.1}" fill="none" stroke="{col}" stroke-width="{sw:.1}"/>"#
            );

            let mx = (x1 + x2) * 0.5;
            if br.is_speciation_split {
                let _ = writeln!(
                    svg,
                    r##"<polygon points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}" fill="#ffd166" stroke="#0b1020" stroke-width="1"/>"##,
                    x1,
                    y2 - 5.5,
                    x1 + 5.5,
                    y2,
                    x1,
                    y2 + 5.5,
                    x1 - 5.5,
                    y2
                );
            }
            if br.is_thymine_udg_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.3" fill="#10b981" stroke="#fef08a" stroke-width="2.4"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a7f3d0" font-family="monospace" font-size="9.5">U-DNA→T-DNA(UDG)</text>"##,
                    y2 - 8.5
                );
            } else if br.is_rnr_dna_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.8" fill="#06b6d4" stroke="#cffafe" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a5f3fc" font-family="monospace" font-size="9.5">RNA→U-DNA(RNR)</text>"##,
                    y2 - 8.5
                );
            } else if br.is_sgc_crystallization_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.2" fill="#38bdf8" stroke="#fef08a" stroke-width="2.4"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#bae6fd" font-family="monospace" font-size="9.5">StatCode→UniversalSGC</text>"##,
                    y2 - 8.5
                );
            } else if br.is_homochiral_symmetry_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.7" fill="#a855f7" stroke="#f3e8ff" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#e9d5ff" font-family="monospace" font-size="9.5">Racemic→Homochiral</text>"##,
                    y2 - 8.5
                );
            } else if br.is_darwinian_threshold_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.2" fill="#f59e0b" stroke="#fef08a" stroke-width="2.4"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#fde68a" font-family="monospace" font-size="9.5">Progenote→LUCA</text>"##,
                    y2 - 8.5
                );
            } else if br.is_raf_protocell_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.6" fill="#ec4899" stroke="#fbcfe8" stroke-width="2.1"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#f9a8d4" font-family="monospace" font-size="9.5">Soup→RAFProtocell</text>"##,
                    y2 - 8.5
                );
            } else if br.is_homeostatic_awakening_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.6" fill="#10b981" stroke="#fef08a" stroke-width="2.5"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a7f3d0" font-family="monospace" font-size="9.5">Techno→Gaia2Homeostasis</text>"##,
                    y2 - 8.5
                );
            } else if br.is_goe_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.9" fill="#38bdf8" stroke="#bae6fd" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#7dd3fc" font-family="monospace" font-size="9.5">Archean→GOE</text>"##,
                    y2 - 8.5
                );
            } else if br.is_formal_verifier_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.5" fill="#8b5cf6" stroke="#ede9fe" stroke-width="2.4"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#ddd6fe" font-family="monospace" font-size="9.5">AlgoMed→FormalSymbio</text>"##,
                    y2 - 8.5
                );
            } else if br.is_machine_symbiosis_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.8" fill="#60a5fa" stroke="#dbeafe" stroke-width="2.1"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#bfdbfe" font-family="monospace" font-size="9.5">Inst→AlgoMed</text>"##,
                    y2 - 8.5
                );
            } else if br.is_peer_institution_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.4" fill="#14b8a6" stroke="#ccfbf1" stroke-width="2.4"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#99f6e4" font-family="monospace" font-size="9.5">WrittenArch→VerifiedInst</text>"##,
                    y2 - 8.5
                );
            } else if br.is_written_archive_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.8" fill="#f8fafc" stroke="#93c5fd" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#e2e8f0" font-family="monospace" font-size="9.5">Symbolic→WrittenArch</text>"##,
                    y2 - 8.5
                );
            } else if br.is_symbolic_teaching_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.2" fill="#facc15" stroke="#fef9c3" stroke-width="2.3"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#fef08a" font-family="monospace" font-size="9.5">Ratchet→SymbolicCCE</text>"##,
                    y2 - 8.5
                );
            } else if br.is_cultural_ratchet_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.6" fill="#38bdf8" stroke="#e0f2fe" stroke-width="2.1"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#bae6fd" font-family="monospace" font-size="9.5">Tradition→Ratchet</text>"##,
                    y2 - 8.5
                );
            } else if br.is_obligate_crop_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="7.0" fill="#84cc16" stroke="#ecfccb" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#d9f99d" font-family="monospace" font-size="9.5">LowerAg→GongylCrop</text>"##,
                    y2 - 8.5
                );
            } else if br.is_agriculture_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.5" fill="#22d3ee" stroke="#cffafe" stroke-width="2.0"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a5f3fc" font-family="monospace" font-size="9.5">Superorg→LowerAg</text>"##,
                    y2 - 8.5
                );
            } else if br.is_superorganism_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.8" fill="#f59e0b" stroke="#fef08a" stroke-width="2.2"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#fde68a" font-family="monospace" font-size="9.5">Eusocial→Superorg</text>"##,
                    y2 - 8.5
                );
            } else if br.is_eusocial_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.4" fill="#10b981" stroke="#a7f3d0" stroke-width="2.0"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a7f3d0" font-family="monospace" font-size="9.5">GermSoma→Eusocial</text>"##,
                    y2 - 8.5
                );
            } else if br.is_germ_soma_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.6" fill="#ffd166" stroke="#38bdf8" stroke-width="2.1"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#fef08a" font-family="monospace" font-size="9.5">Clonal→GermSoma</text>"##,
                    y2 - 8.5
                );
            } else if br.is_clonal_multicell_transition {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.3" fill="#38bdf8" stroke="#c77dff" stroke-width="2.0"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#bae6fd" font-family="monospace" font-size="9.5">Holo→Clonal</text>"##,
                    y2 - 8.5
                );
            } else {
                let mut has_mro_badge = false;
                if br.is_hydrogenosome_transition {
                    has_mro_badge = true;
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.4" fill="#06d6a0" stroke="#ffd166" stroke-width="2.0"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#a7f3d0" font-family="monospace" font-size="9.5">Mito→Hydro</text>"##,
                        y2 - 8.5
                    );
                } else if br.is_mitosome_transition {
                    has_mro_badge = true;
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.4" fill="#ffd166" stroke="#ff5c8a" stroke-width="2.0"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#fde68a" font-family="monospace" font-size="9.5">Hydro→MitoS</text>"##,
                        y2 - 8.5
                    );
                } else if br.is_mro_complete_loss {
                    has_mro_badge = true;
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.2" fill="#38bdf8" stroke="#ffd166" stroke-width="2.0"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#bae6fd" font-family="monospace" font-size="9.5">MRO→0</text>"##,
                        y2 - 8.5
                    );
                }

                let has_apico_badge =
                    br.is_apicoplast_transition || br.is_apicoplast_complete_loss;
                let apico_y = if has_mro_badge { y2 + 12.5 } else { y2 - 8.5 };
                if br.is_apicoplast_transition {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.4" fill="#ff5c8a" stroke="#c77dff" stroke-width="2.1"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<text x="{mx:.1}" y="{apico_y:.1}" text-anchor="middle" fill="#f3c4fb" font-family="monospace" font-size="9.5">Holo→Apico</text>"##
                    );
                } else if br.is_apicoplast_complete_loss {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.2" fill="#00f5d4" stroke="#ffd166" stroke-width="2.0"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<text x="{mx:.1}" y="{apico_y:.1}" text-anchor="middle" fill="#9bf6ff" font-family="monospace" font-size="9.5">Apico→0</text>"##
                    );
                }

                if !has_mro_badge && !has_apico_badge {
                    if br.is_klepto_tipping_point {
                        let _ = writeln!(
                            svg,
                            r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.4" fill="#80ed99" stroke="#00f5d4" stroke-width="2.0"/>"##
                        );
                        let _ = writeln!(
                            svg,
                            r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#b7efc5" font-family="monospace" font-size="9.5">KP→FIX</text>"##,
                            y2 - 8.5
                        );
                    } else if br.is_nucleomorph_loss {
                        let _ = writeln!(
                            svg,
                            r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="6.2" fill="#ffd166" stroke="#f72585" stroke-width="1.9"/>"##
                        );
                        let _ = writeln!(
                            svg,
                            r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#ffe66d" font-family="monospace" font-size="9.5">NM→0</text>"##,
                            y2 - 8.5
                        );
                    } else if br.is_egt_lockin {
                        let _ = writeln!(
                            svg,
                            r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="5.8" fill="#00f5d4" stroke="#c77dff" stroke-width="1.8"/>"##
                        );
                        let _ = writeln!(
                            svg,
                            r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#80ffdb" font-family="monospace" font-size="9.5">EGT→Nuc</text>"##,
                            y2 - 8.5
                        );
                    } else if let Some((from_g, to_g)) = br.guild_transition {
                        let badge_col = guild_hex(to_g);
                        let _ = writeln!(
                            svg,
                            r##"<circle cx="{mx:.1}" cy="{y2:.1}" r="5.2" fill="{badge_col}" stroke="#ffffff" stroke-width="1.3"/>"##
                        );
                        let _ = writeln!(
                            svg,
                            r##"<text x="{mx:.1}" y="{:.1}" text-anchor="middle" fill="#e6eeff" font-family="monospace" font-size="9.5">{}→{}</text>"##,
                            y2 - 8.0,
                            from_g.short_code(),
                            to_g.short_code()
                        );
                    }
                }
            }
        }

        for ret in &layout.reticulations {
            let x1 = pad_l + ret.from_x * plot_w;
            let y1 = pad_t + ret.from_y * plot_h;
            let x2 = pad_l + ret.to_x * plot_w;
            let y2 = pad_t + ret.to_y * plot_h;
            if ret.is_endosymbiosis {
                let c1x = x1 + (x2 - x1).abs().max(18.0) * 0.55;
                let c2x = x2 - (x2 - x1).abs().max(18.0) * 0.35;
                let arc_col = if ret.is_klepto_tipping {
                    "#80ed99"
                } else {
                    match ret.endo_tier {
                        3..=u8::MAX => "#f72585",
                        2 => "#ff9f1c",
                        _ => "#50e389",
                    }
                };
                let _ = writeln!(
                    svg,
                    r##"<path d="M {x1:.1} {y1:.1} C {c1x:.1} {y1:.1}, {c2x:.1} {y2:.1}, {x2:.1} {y2:.1}" fill="none" stroke="{arc_col}" stroke-width="2.4" stroke-dasharray="6,3"/>"##
                );
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{x2:.1}" cy="{y2:.1}" r="7.0" fill="{arc_col}" stroke="#ffd166" stroke-width="1.8"/>"##
                );
                if let Some(hid) = ret.hybrid_id {
                    let tag_lbl = if ret.is_klepto_tipping {
                        format!("#H{hid} (KP→FIX)")
                    } else {
                        format!("#H{hid} ({}° γ={:.2})", ret.endo_tier.max(1), ret.gamma)
                    };
                    let _ = writeln!(
                        svg,
                        r##"<text x="{:.1}" y="{:.1}" fill="#ffd166" font-family="monospace" font-size="10.5" font-weight="bold">{}</text>"##,
                        x2 + 8.5,
                        y2 - 6.5,
                        escape_svg_xml(&tag_lbl)
                    );
                }
            } else {
                let cx = (x1 + x2) * 0.5;
                let cy = y1.min(y2) - 18.0;
                let _ = writeln!(
                    svg,
                    r##"<path d="M {x1:.1} {y1:.1} Q {cx:.1} {cy:.1} {x2:.1} {y2:.1}" fill="none" stroke="#80ffdb" stroke-width="1.4" stroke-dasharray="4,3"/>"##
                );
            }
        }

        for nd in &layout.nodes {
            let nx = pad_l + nd.x_norm * plot_w;
            let ny = pad_t + nd.y_norm * plot_h;
            let col = guild_hex(nd.guild);
            if nd.extant {
                if nd.mro_class == MroClass::Hydrogenosome {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="5.8" fill="#06d6a0" stroke="#ffd166" stroke-width="2.0"/>"##
                    );
                } else if nd.mro_class == MroClass::Mitosome {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="5.8" fill="#ffd166" stroke="#ff5c8a" stroke-width="2.0"/>"##
                    );
                } else if nd.mro_class == MroClass::Amiteukaryote {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="5.6" fill="#38bdf8" stroke="#ffd166" stroke-width="1.8"/>"##
                    );
                } else if nd.guild == TrophicGuild::Holobiont {
                    let outer_col = match nd.endo_tier {
                        3..=u8::MAX => "#f72585",
                        2 => "#ff9f1c",
                        _ => "#c77dff",
                    };
                    let ring_col = if nd.nucleomorph_loci > 0 {
                        "#ffe66d"
                    } else if nd.egt_nuclear_fraction >= 0.50 {
                        "#00f5d4"
                    } else {
                        "#50e389"
                    };
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="6.0" fill="{outer_col}" stroke="{ring_col}" stroke-width="1.8"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="2.3" fill="{ring_col}"/>"##
                    );
                } else if nd.has_apicoplast {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="5.8" fill="#ff5c8a" stroke="#c77dff" stroke-width="2.0"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="2.2" fill="#ffd166"/>"##
                    );
                } else if nd.has_kleptoplast {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="5.6" fill="#ff5c8a" stroke="#80ed99" stroke-width="1.8"/>"##
                    );
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="2.1" fill="#80ed99"/>"##
                    );
                } else {
                    let _ = writeln!(
                        svg,
                        r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="4.4" fill="{col}" stroke="#0b1020" stroke-width="1"/>"##
                    );
                }
                let lbl_col = match nd.mro_class {
                    MroClass::Hydrogenosome => "#a7f3d0",
                    MroClass::Mitosome => "#fde68a",
                    MroClass::Amiteukaryote => "#bae6fd",
                    _ if nd.has_apicoplast => "#f3c4fb",
                    _ => col,
                };
                let _ = writeln!(
                    svg,
                    r##"<text x="{:.1}" y="{:.1}" fill="{lbl_col}" font-family="monospace" font-size="11">{}</text>"##,
                    nx + 8.5,
                    ny + 3.8,
                    escape_svg_xml(&nd.label)
                );
            } else {
                let _ = writeln!(
                    svg,
                    r##"<circle cx="{nx:.1}" cy="{ny:.1}" r="2.8" fill="#9fb6e0"/>"##
                );
            }
        }

        let _ = writeln!(svg, "</svg>");
        svg
    }
}

fn escape_svg_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
