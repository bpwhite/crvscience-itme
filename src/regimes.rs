//! ITME v23.0 (`v2.3.0`) Emergent Regime Taxonomy & Purely Observable Phase Classifier.
//!
//! Classifies evolutionary regimes **emergently** from live measured invariants across
//! anaerobic MRO reduction (`AerobicMito -> Hydrogenosome -> Mitosome -> MRO->0`), reductive
//! apicoplast evolution (`Chromera -> Plasmodium -> Cryptosporidium`), mixotrophic kleptoplasty,
//! shopping-bag HGT tipping points, `1°/2°/3°` serial endosymbiosis, nucleomorph compaction,
//! nuclear `EGT`, multi-trophic spatial webs, and speciating demes.

use crate::model::{
    CapabilityFlags, ContinuousPhaseCoords, DiagnosticVector, EnvironmentControls,
    PriceDecomposition,
};
use std::fmt;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(usize)]
pub enum RegimeId {
    R00Static = 0,
    R05Decaying = 1,
    R10AdaptiveDissipation = 2,
    R15Autoadaptive = 3,
    R15Replicative = 4,
    R20InteractiveReplicative = 5,
    R20AutoadaptiveReplicative = 6,
    R20ProtoDarwinian = 7,
    R30FullyDarwinian = 8,
    R30SoftInheritance = 9,
    R35Symbiogenetic = 10,
    R30PseudoDarwinian = 11,
    R30PseudoSoftInheritance = 12,
}

impl RegimeId {
    pub const ALL: [RegimeId; 13] = [
        RegimeId::R00Static,
        RegimeId::R05Decaying,
        RegimeId::R10AdaptiveDissipation,
        RegimeId::R15Autoadaptive,
        RegimeId::R15Replicative,
        RegimeId::R20InteractiveReplicative,
        RegimeId::R20AutoadaptiveReplicative,
        RegimeId::R20ProtoDarwinian,
        RegimeId::R30FullyDarwinian,
        RegimeId::R30SoftInheritance,
        RegimeId::R35Symbiogenetic,
        RegimeId::R30PseudoDarwinian,
        RegimeId::R30PseudoSoftInheritance,
    ];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn code(self) -> &'static str {
        match self {
            Self::R00Static => "R(0.0)",
            Self::R05Decaying => "R(0.5)",
            Self::R10AdaptiveDissipation => "R(1.0)",
            Self::R15Autoadaptive => "R(1.5-AA)",
            Self::R15Replicative => "R(1.5-RP)",
            Self::R20InteractiveReplicative => "R(2.0-EC)",
            Self::R20AutoadaptiveReplicative => "R(2.0-BAL)",
            Self::R20ProtoDarwinian => "R(2.0-QS)",
            Self::R30FullyDarwinian => "R(3.0-DAR)",
            Self::R30SoftInheritance => "R(3.0-CUL)",
            Self::R35Symbiogenetic => "R(3.5-SYM)",
            Self::R30PseudoDarwinian => "R(4.0-ENG)",
            Self::R30PseudoSoftInheritance => "R(4.0-AUT)",
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::R00Static => "Static Equilibrium",
            Self::R05Decaying => "Entropic Decay",
            Self::R10AdaptiveDissipation => "Dissipative Self-Organization",
            Self::R15Autoadaptive => "Autoadaptive Learning",
            Self::R15Replicative => "Unbounded Clonal Replication",
            Self::R20InteractiveReplicative => "Ecological Clonal Competition",
            Self::R20AutoadaptiveReplicative => "Baldwinian Autoadaptive-Replicator",
            Self::R20ProtoDarwinian => "Quasispecies Error Catastrophe",
            Self::R30FullyDarwinian => "Weismannian Darwinian Selection",
            Self::R30SoftInheritance => "Cumulative Soft / Cultural Evolution",
            Self::R35Symbiogenetic => "Symbiogenesis -> Hydrogenosomes, Mitosomes & Apicoplasts",
            Self::R30PseudoDarwinian => "Directed / Meta-Engineered Darwinian",
            Self::R30PseudoSoftInheritance => "Recursive Self-Distilling & Organelle-Streamlined",
        }
    }

    pub fn full_title(self) -> String {
        format!("{} {}", self.code(), self.name())
    }

    pub fn core_concept(self) -> &'static str {
        match self {
            Self::R00Static => "Zero net exergy flux and zero entropy production at equilibrium",
            Self::R05Decaying => "Landauer maintenance cost exceeds harvested work (Phi_net < 0)",
            Self::R10AdaptiveDissipation => "Driven dissipative structure maintained by positive exergy margin",
            Self::R15Autoadaptive => "Single entity increases functional work extraction via lifetime plasticity",
            Self::R15Replicative => "High-MI template copying in an uncrowded spatial resource grid",
            Self::R20InteractiveReplicative => "Crowded zero-variance clones competing on a finite exergy grid",
            Self::R20AutoadaptiveReplicative => "Plastic foraging & niche use improve survival without assimilation",
            Self::R20ProtoDarwinian => "Negative Price transmission bias overwhelms selection below Eigen's bound",
            Self::R30FullyDarwinian => "Multi-trophic host-parasite & syntrophic co-evolution + speciation",
            Self::R30SoftInheritance => "Multi-mentor pedagogy + niche & cross-feeding inheritance",
            Self::R35Symbiogenetic => "1°/2°/3° endosymbiosis, anaerobic hydrogenosomes/mitosomes (ISC) & apicoplasts",
            Self::R30PseudoDarwinian => "External meta-objective sculpts multi-trophic spatial variation",
            Self::R30PseudoSoftInheritance => "Recursive self-distillation + SUF/IPP-rescued organelle streamlining",
        }
    }

    pub fn transition_threshold_note(self) -> &'static str {
        match self {
            Self::R00Static => "Supply bath exergy flux J_in > Q_landauer -> R(1.0)",
            Self::R05Decaying => "Restore harvested work W_harvest > Q_landauer -> R(1.0)",
            Self::R10AdaptiveDissipation => "Gain plastic state M -> R(1.5-AA) or template channel R -> R(1.5-RP)",
            Self::R15Autoadaptive => "Couple plastic entity to template replication -> R(2.0-BAL)",
            Self::R15Replicative => "Saturate spatial carrying capacity -> R(2.0-EC)",
            Self::R20InteractiveReplicative => "Introduce heritable variation mu > 0 -> R(2.0-QS / R3.0-DAR)",
            Self::R20AutoadaptiveReplicative => "Cross Eigen bound with variation -> R(3.0-DAR) or pedagogy -> R(3.0-CUL)",
            Self::R20ProtoDarwinian => "Lift channel capacity I(S;S') above theta_F -> R(3.0-DAR)",
            Self::R30FullyDarwinian => "Open oblique assimilation -> R(3.0-CUL) or symbiogenesis -> R(3.5-SYM)",
            Self::R30SoftInheritance => "Fix endosymbionts & evolve MRO / apicoplast lineages -> R(3.5-SYM)",
            Self::R35Symbiogenetic => "Couple MRO/plastid lineages to recursive meta-objective -> R(4.0-AUT)",
            Self::R30PseudoDarwinian => "Enable teacher-student distillation + organelle streamlining -> R(4.0-AUT)",
            Self::R30PseudoSoftInheritance => "Autonomous open-ended symbiogenetic & reductive-organelle evolution",
        }
    }

    pub const fn canonical_dc(self) -> DiagnosticVector {
        match self {
            Self::R00Static => DiagnosticVector::new(0, 0, 0, 0, 0),
            Self::R05Decaying => DiagnosticVector::new(1, 0, 0, 0, 0),
            Self::R10AdaptiveDissipation => DiagnosticVector::new(2, 0, 0, 0, 1),
            Self::R15Autoadaptive => DiagnosticVector::new(2, 2, 0, 0, 1),
            Self::R15Replicative => DiagnosticVector::new(2, 1, 3, 0, 1),
            Self::R20InteractiveReplicative => DiagnosticVector::new(2, 1, 3, 0, 2),
            Self::R20AutoadaptiveReplicative => DiagnosticVector::new(2, 2, 3, 0, 2),
            Self::R20ProtoDarwinian => DiagnosticVector::new(2, 1, 2, 1, 3),
            Self::R30FullyDarwinian => DiagnosticVector::new(2, 2, 4, 2, 3),
            Self::R30SoftInheritance => DiagnosticVector::new(2, 2, 4, 3, 3),
            Self::R35Symbiogenetic => DiagnosticVector::new(2, 2, 5, 3, 3),
            Self::R30PseudoDarwinian => DiagnosticVector::new(2, 2, 4, 2, 4),
            Self::R30PseudoSoftInheritance => DiagnosticVector::new(2, 2, 5, 3, 4),
        }
    }

    pub fn from_slug(input: &str) -> Option<Self> {
        let norm = input.trim().to_ascii_uppercase();
        match norm.as_str() {
            "0" | "0.0" | "R(0.0)" | "STATIC" => Some(Self::R00Static),
            "1" | "0.5" | "R(0.5)" | "DECAY" | "DECAYING" => Some(Self::R05Decaying),
            "2" | "1.0" | "R(1.0)" | "DISSIPATION" => Some(Self::R10AdaptiveDissipation),
            "3" | "1.5-AA" | "R(1.5-AA)" | "AA" => Some(Self::R15Autoadaptive),
            "4" | "1.5-RP" | "R(1.5-RP)" | "RP" => Some(Self::R15Replicative),
            "5" | "2.0-EC" | "2.0-IF" | "2.0-HF" | "R(2.0-EC)" => {
                Some(Self::R20InteractiveReplicative)
            }
            "6" | "2.0-BAL" | "2.0-AR" | "R(2.0-BAL)" => Some(Self::R20AutoadaptiveReplicative),
            "7" | "2.0-QS" | "2.0-PD" | "R(2.0-QS)" => Some(Self::R20ProtoDarwinian),
            "8" | "3.0-DAR" | "3.0-FD" | "R(3.0-DAR)" => Some(Self::R30FullyDarwinian),
            "9" | "3.0-CUL" | "3.0-SI" | "R(3.0-CUL)" => Some(Self::R30SoftInheritance),
            "10" | "3.5-SYM" | "R(3.5-SYM)" | "SYM" => Some(Self::R35Symbiogenetic),
            "11" | "4.0-ENG" | "3.0-PSD" | "R(4.0-ENG)" => Some(Self::R30PseudoDarwinian),
            "12" | "4.0-AUT" | "3.0-PSI" | "R(4.0-AUT)" => Some(Self::R30PseudoSoftInheritance),
            _ => None,
        }
    }

    pub fn preset_controls(self) -> EnvironmentControls {
        match self {
            Self::R00Static => EnvironmentControls {
                bath_exergy_flux: 0.0,
                landauer_bit_cost: 0.0,
                niche_construction_rate: 0.0,
                carrying_capacity: 1,
                ..EnvironmentControls::default()
            },
            Self::R05Decaying => EnvironmentControls {
                bath_exergy_flux: -0.35,
                landauer_bit_cost: 0.0035,
                niche_construction_rate: 0.0,
                carrying_capacity: 4,
                ..EnvironmentControls::default()
            },
            Self::R10AdaptiveDissipation => EnvironmentControls {
                bath_exergy_flux: 0.72,
                niche_construction_rate: 0.04,
                carrying_capacity: 6,
                ..EnvironmentControls::default()
            },
            Self::R15Autoadaptive => EnvironmentControls {
                bath_exergy_flux: 0.80,
                memory_plasticity: 0.80,
                niche_construction_rate: 0.12,
                carrying_capacity: 1,
                ..EnvironmentControls::default()
            },
            Self::R15Replicative => EnvironmentControls {
                bath_exergy_flux: 0.82,
                memory_plasticity: 0.22,
                replication_rate: 0.40,
                replication_fidelity: 0.995,
                niche_construction_rate: 0.08,
                carrying_capacity: 26,
                ..EnvironmentControls::default()
            },
            Self::R20InteractiveReplicative => EnvironmentControls {
                bath_exergy_flux: 0.84,
                memory_plasticity: 0.22,
                replication_rate: 0.44,
                replication_fidelity: 0.995,
                interaction_strength: 0.80,
                niche_construction_rate: 0.10,
                carrying_capacity: 24,
                ..EnvironmentControls::default()
            },
            Self::R20AutoadaptiveReplicative => EnvironmentControls {
                bath_exergy_flux: 0.88,
                memory_plasticity: 0.82,
                replication_rate: 0.42,
                replication_fidelity: 0.995,
                interaction_strength: 0.80,
                niche_construction_rate: 0.16,
                carrying_capacity: 24,
                ..EnvironmentControls::default()
            },
            Self::R20ProtoDarwinian => EnvironmentControls {
                bath_exergy_flux: 0.88,
                memory_plasticity: 0.25,
                replication_rate: 0.50,
                replication_fidelity: 0.65,
                heritable_variation: 0.34,
                interaction_strength: 0.84,
                niche_construction_rate: 0.06,
                carrying_capacity: 28,
                selective_superiority: 5.5,
                ..EnvironmentControls::default()
            },
            Self::R30FullyDarwinian => EnvironmentControls {
                bath_exergy_flux: 0.98,
                memory_plasticity: 0.82,
                replication_rate: 0.54,
                replication_fidelity: 0.965,
                heritable_variation: 0.08,
                interaction_strength: 0.88,
                niche_construction_rate: 0.18,
                trophic_coupling: 0.36,
                trophic_efficiency: 0.38,
                cross_feeding_rate: 0.28,
                speciation_threshold: 0.24,
                carrying_capacity: 32,
                selective_superiority: 10.5,
                env_drift_rate: 0.015,
                ..EnvironmentControls::default()
            },
            Self::R30SoftInheritance => EnvironmentControls {
                bath_exergy_flux: 1.04,
                memory_plasticity: 0.88,
                replication_rate: 0.56,
                replication_fidelity: 0.74,
                heritable_variation: 0.09,
                interaction_strength: 0.90,
                niche_construction_rate: 0.24,
                trophic_coupling: 0.34,
                trophic_efficiency: 0.40,
                cross_feeding_rate: 0.32,
                speciation_threshold: 0.24,
                carrying_capacity: 32,
                soft_inheritance_coupling: 0.48,
                horizontal_pedagogy: 0.86,
                teacher_count: 5,
                selective_superiority: 11.5,
                env_drift_rate: 0.025,
                ..EnvironmentControls::default()
            },
            Self::R35Symbiogenetic => EnvironmentControls {
                bath_exergy_flux: 1.15,
                memory_plasticity: 0.90,
                replication_rate: 0.58,
                replication_fidelity: 0.78,
                heritable_variation: 0.09,
                interaction_strength: 0.92,
                niche_construction_rate: 0.28,
                trophic_coupling: 0.36,
                trophic_efficiency: 0.42,
                cross_feeding_rate: 0.40,
                kleptoplasty_rate: 0.48,
                klepto_hgt_rate: 0.46,
                endosymbiosis_rate: 0.46,
                egt_rate: 0.44,
                serial_endosymbiosis_rate: 0.38,
                nucleomorph_reduction_rate: 0.46,
                apicoplast_reduction_rate: 0.34,
                host_metabolite_rescue: 0.15,
                ambient_oxygen_level: 0.58,
                anaerobic_mro_reduction_rate: 0.32,
                cytosolic_suf_hgt_rescue: 0.15,
                multicellular_adhesion_rate: 0.44,
                germ_soma_differentiation_rate: 0.40,
                apoptosis_policing_rate: 0.65,
                single_cell_bottleneck: true,
                queen_monogamy_index: 0.85,
                eusocial_caste_canalization_rate: 0.42,
                worker_policing_rate: 0.72,
                crop_domestication_rate: 0.44,
                agro_pathogen_pressure: 0.38,
                defensive_symbiont_potency: 0.76,
                symbolic_teaching_fidelity: 0.89,
                inter_deme_connectivity: 0.36,
                cultural_innovation_rate: 0.46,
                archival_externalization_rate: 0.48,
                peer_verification_rigor: 0.78,
                unverified_publish_pressure: 0.32,
                synthetic_generation_rate: 0.50,
                formal_verifier_strictness: 0.82,
                empirical_grounding_floor: 0.28,
                biogeochemical_feedback_gain: 0.52,
                planetary_stewardship_coupling: 0.78,
                extractive_burn_intensity: 0.20,
                speciation_threshold: 0.26,
                carrying_capacity: 34,
                soft_inheritance_coupling: 0.52,
                horizontal_pedagogy: 0.88,
                teacher_count: 5,
                recombination_rate: 0.45,
                selective_superiority: 12.5,
                env_drift_rate: 0.035,
                ..EnvironmentControls::default()
            },
            Self::R30PseudoDarwinian => EnvironmentControls {
                bath_exergy_flux: 1.12,
                memory_plasticity: 0.85,
                replication_rate: 0.58,
                replication_fidelity: 0.975,
                heritable_variation: 0.09,
                interaction_strength: 0.94,
                niche_construction_rate: 0.22,
                trophic_coupling: 0.34,
                trophic_efficiency: 0.40,
                cross_feeding_rate: 0.30,
                speciation_threshold: 0.26,
                carrying_capacity: 32,
                selective_superiority: 13.0,
                artificial_selection: true,
                ..EnvironmentControls::default()
            },
            Self::R30PseudoSoftInheritance => EnvironmentControls {
                bath_exergy_flux: 1.26,
                memory_plasticity: 0.94,
                replication_rate: 0.62,
                replication_fidelity: 0.84,
                heritable_variation: 0.10,
                interaction_strength: 0.96,
                niche_construction_rate: 0.32,
                trophic_coupling: 0.34,
                trophic_efficiency: 0.44,
                cross_feeding_rate: 0.44,
                kleptoplasty_rate: 0.52,
                klepto_hgt_rate: 0.50,
                endosymbiosis_rate: 0.50,
                egt_rate: 0.50,
                serial_endosymbiosis_rate: 0.44,
                nucleomorph_reduction_rate: 0.52,
                apicoplast_reduction_rate: 0.40,
                host_metabolite_rescue: 0.72,
                ambient_oxygen_level: 0.48,
                anaerobic_mro_reduction_rate: 0.38,
                cytosolic_suf_hgt_rescue: 0.72,
                multicellular_adhesion_rate: 0.52,
                germ_soma_differentiation_rate: 0.48,
                apoptosis_policing_rate: 0.76,
                single_cell_bottleneck: true,
                queen_monogamy_index: 0.90,
                eusocial_caste_canalization_rate: 0.52,
                worker_policing_rate: 0.85,
                crop_domestication_rate: 0.54,
                agro_pathogen_pressure: 0.42,
                defensive_symbiont_potency: 0.86,
                symbolic_teaching_fidelity: 0.94,
                inter_deme_connectivity: 0.38,
                cultural_innovation_rate: 0.62,
                archival_externalization_rate: 0.66,
                peer_verification_rigor: 0.88,
                unverified_publish_pressure: 0.36,
                synthetic_generation_rate: 0.72,
                formal_verifier_strictness: 0.92,
                empirical_grounding_floor: 0.35,
                biogeochemical_feedback_gain: 0.68,
                planetary_stewardship_coupling: 0.88,
                extractive_burn_intensity: 0.22,
                speciation_threshold: 0.26,
                carrying_capacity: 34,
                soft_inheritance_coupling: 0.68,
                horizontal_pedagogy: 0.94,
                teacher_count: 6,
                recombination_rate: 0.52,
                selective_superiority: 15.0,
                env_drift_rate: 0.04,
                artificial_selection: true,
                ..EnvironmentControls::default()
            },
        }
    }

    pub fn outgoing_transitions(self) -> &'static [(RegimeId, &'static str, f64)] {
        match self {
            Self::R00Static => &[
                (Self::R10AdaptiveDissipation, "Bath exergy J_in > Q_landauer", 0.85),
                (Self::R30PseudoDarwinian, "External fabrication", 0.15),
            ],
            Self::R05Decaying => &[
                (Self::R00Static, "Thermal exhaustion", 0.55),
                (Self::R10AdaptiveDissipation, "Exergy restoration", 0.45),
            ],
            Self::R10AdaptiveDissipation => &[
                (Self::R05Decaying, "Landauer starvation", 0.20),
                (Self::R15Autoadaptive, "Conformational/neural plasticity", 0.40),
                (Self::R15Replicative, "Autocatalytic template emergence", 0.40),
            ],
            Self::R15Autoadaptive => &[
                (Self::R20AutoadaptiveReplicative, "Couple plasticity to copying", 0.85),
                (Self::R05Decaying, "Exergy collapse", 0.15),
            ],
            Self::R15Replicative => &[
                (Self::R20InteractiveReplicative, "Carrying-capacity crowding", 0.65),
                (Self::R20ProtoDarwinian, "Unproofread mutation burst", 0.35),
            ],
            Self::R20InteractiveReplicative => &[
                (Self::R20AutoadaptiveReplicative, "Plasticity acquisition", 0.40),
                (Self::R20ProtoDarwinian, "Mutational exploration", 0.60),
            ],
            Self::R20AutoadaptiveReplicative => &[
                (Self::R30FullyDarwinian, "Proofread germline variation", 0.50),
                (Self::R30SoftInheritance, "Pedagogical skill assimilation", 0.50),
            ],
            Self::R20ProtoDarwinian => &[
                (Self::R05Decaying, "Runaway Landauer meltdown", 0.25),
                (Self::R30FullyDarwinian, "Exonuclease proofreading", 0.40),
                (Self::R30SoftInheritance, "Multi-mentor social correction", 0.35),
            ],
            Self::R30FullyDarwinian => &[
                (Self::R30SoftInheritance, "Evolve oblique pedagogy", 0.45),
                (Self::R35Symbiogenetic, "Endosymbiosis + anaerobic MRO evolution", 0.30),
                (Self::R30PseudoDarwinian, "Directed selection", 0.25),
            ],
            Self::R30SoftInheritance => &[
                (Self::R35Symbiogenetic, "Organelle fixation & MRO/apicoplast reduction", 0.55),
                (Self::R30PseudoSoftInheritance, "Recursive meta-distillation", 0.45),
            ],
            Self::R35Symbiogenetic => &[
                (Self::R30PseudoSoftInheritance, "Recursive meta-architecture + SUF rescue", 0.65),
                (Self::R30SoftInheritance, "Lineage isolation", 0.35),
            ],
            Self::R30PseudoDarwinian => &[
                (Self::R30PseudoSoftInheritance, "Add distillation + MRO/plastid reduction", 0.70),
                (Self::R30FullyDarwinian, "Autonomous ecological release", 0.30),
            ],
            Self::R30PseudoSoftInheritance => &[
                (Self::R35Symbiogenetic, "Open ecosystem release", 0.70),
                (Self::R00Static, "Power shutdown", 0.30),
            ],
        }
    }
}

impl fmt::Display for RegimeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.code(), self.name())
    }
}

/// Purely observable-driven ITME v12.0 (`v1.2.0`) classifier!
#[allow(clippy::too_many_arguments)]
pub fn classify_emergent_state(
    bath_flux: f64,
    net_exergy_margin: f64,
    order_velocity_vo: f64,
    plastic_gain_ratio: f64,
    empirical_birth_activity: f64,
    crowding_pressure: f64,
    adami_rho: f64,
    blended_fidelity: f64,
    critical_theta_f: f64,
    assimilation_activity: f64,
    recombination_activity: f64,
    price: &PriceDecomposition,
    artificial_meta_context: bool,
) -> (CapabilityFlags, ContinuousPhaseCoords, DiagnosticVector, RegimeId) {
    let is_static = bath_flux.abs() <= 0.03 && net_exergy_margin.abs() <= 0.03;
    let is_decaying = !is_static && (net_exergy_margin < -0.02 || order_velocity_vo < -0.22);
    let exergy_driven = !is_static && !is_decaying;

    let v_oc: u8 = if is_static {
        0
    } else if is_decaying {
        1
    } else {
        2
    };

    let m_c: u8 = if !exergy_driven || plastic_gain_ratio < 0.03 {
        0
    } else if plastic_gain_ratio < 0.16 {
        1
    } else {
        2
    };

    let has_rep = exergy_driven && empirical_birth_activity >= 0.04;
    let has_var = has_rep
        && ((1.0 - adami_rho) > 0.025
            || (price.selection_covariance.abs() > 0.002 && (1.0 - adami_rho) > 0.004)
            || assimilation_activity > 0.05);
    let above_eigen = has_rep && blended_fidelity >= critical_theta_f;
    let has_soft = has_rep && m_c == 2 && assimilation_activity >= 0.08 && above_eigen;
    let has_symbio = has_rep && recombination_activity >= 0.10 && above_eigen;

    let r_c: u8 = if !has_rep {
        0
    } else if !has_var {
        if blended_fidelity >= 0.85 {
            3
        } else {
            1
        }
    } else if !above_eigen {
        2
    } else if has_symbio {
        5
    } else {
        4
    };

    let v_c: u8 = if !has_var {
        0
    } else if !above_eigen {
        1
    } else if has_soft {
        3
    } else {
        2
    };

    let has_interaction = exergy_driven && crowding_pressure >= 0.18;
    let p_c: u8 = if !exergy_driven {
        0
    } else if artificial_meta_context && has_rep && has_var && above_eigen {
        4
    } else if has_rep && has_var && has_interaction {
        3
    } else if has_rep && has_interaction {
        2
    } else {
        1
    };

    let flags = CapabilityFlags {
        exergy_driven,
        memory: m_c > 0,
        replication: has_rep,
        interaction: has_interaction,
        variation: has_var,
        selection: has_rep && has_var && has_interaction && above_eigen,
        soft_inheritance: has_soft,
        symbiogenesis: has_symbio,
        kleptoplasty: false,
        klepto_tipping_fixation: false,
        endosymbiotic_fusion: false,
        nuclear_egt: false,
        serial_endosymbiosis: false,
        nucleomorph_reduction: false,
        apicoplast_reduction: false,
        apicoplast_loss: false,
        hydrogenosome_fermentation: false,
        mitosome_isc_retention: false,
        amiteukaryote_mro_loss: false,
        clonal_multicellularity: false,
        germ_soma_division: false,
        apoptotic_cheater_policing: false,
        eusocial_colony: false,
        sterile_worker_castes: false,
        worker_kin_policing: false,
        agricultural_domestication: false,
        clonal_crop_monoculture: false,
        defensive_microbiome_symbiosis: false,
        dual_gene_culture_inheritance: false,
        cumulative_cultural_ratchet: false,
        symbolic_pedagogical_teaching: false,
        exosomatic_memory_offloading: false,
        written_archival_storage: false,
        institutional_peer_verification: false,
        computational_symbiont: false,
        formal_proof_verification: false,
        model_collapse_safeguard: false,
        oxygenic_goe_regulation: false,
        syntrophic_nutrient_closure: false,
        planetary_homeostatic_awakening: false,
        prebiotic_raf_closure: false,
        ribozyme_hypercycle_active: false,
        darwinian_threshold_luca: false,
        homochiral_d_rna_l_aa: false,
        stereoselective_aminoacylation: false,
        error_minimized_universal_sgc: false,
        rnr_deoxyribose_synthesis: false,
        thymine_udg_deamination_repair: false,
        virosphere_dna_coevolution: false,
        niche_construction: false,
        cross_feeding: false,
        trophic_predation: false,
    };

    let dc = DiagnosticVector::new(v_oc, m_c, r_c, v_c, p_c);

    let regime = if v_oc == 0 {
        RegimeId::R00Static
    } else if v_oc == 1 {
        RegimeId::R05Decaying
    } else if !has_rep {
        if m_c >= 2 {
            RegimeId::R15Autoadaptive
        } else {
            RegimeId::R10AdaptiveDissipation
        }
    } else if !has_var {
        if !has_interaction {
            RegimeId::R15Replicative
        } else if m_c >= 2 {
            RegimeId::R20AutoadaptiveReplicative
        } else {
            RegimeId::R20InteractiveReplicative
        }
    } else if !above_eigen {
        RegimeId::R20ProtoDarwinian
    } else if artificial_meta_context {
        if has_soft || has_symbio {
            RegimeId::R30PseudoSoftInheritance
        } else {
            RegimeId::R30PseudoDarwinian
        }
    } else if has_symbio {
        RegimeId::R35Symbiogenetic
    } else if has_soft {
        RegimeId::R30SoftInheritance
    } else {
        RegimeId::R30FullyDarwinian
    };

    let signed_price = if above_eigen {
        price.selection_share
    } else {
        -price.transmission_bias.abs().clamp(0.05, 1.0)
    };

    let coords = ContinuousPhaseCoords {
        net_exergy_margin,
        order_velocity: order_velocity_vo,
        plasticity_fraction: plastic_gain_ratio,
        heritable_capacity: if has_rep { blended_fidelity } else { 0.0 },
        price_selection_index: signed_price,
    };

    (flags, coords, dc, regime)
}
