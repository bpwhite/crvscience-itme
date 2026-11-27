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

//! ITME v23.0 (`v2.3.0`) Core Invariants: Continuous 4.25 Ga 8-Epoch Hadean-to-Gaia Grand Trajectory Engine,
//! RNR RNA-to-DNA Transition, Uracil-to-Thymine UDG Repair & Virosphere Coevolution,
//! D/L Homochiral Symmetry Breaking, Stereoselective Aminoacylation & SGC Crystallization,
//! Hadean Early-Earth Abiogenesis, Prebiotic Geochemistry, Eigen Threshold & LUCA,
//! Planetary Biogeochemical Feedbacks, O2/CO2 Co-Regulation & Gaia-2.0 Homeostasis,
//! Autonomous Machine-Verified Exosomatic Synthesis & Model-Collapse Safeguards,
//! Exosomatic Memory, Written Archives & Peer-Verification Ratchets,
//! Cumulative Cultural Evolution, Dual Inheritance & Symbolic Teaching,
//! Obligate Agricultural Ectosymbioses, Clonal Crops & Defensive Microbiomes,
//! Eusocial Superorganisms, Sterile Castes & Worker Kin-Policing,
//! Clonal Multicellularity, Soma-Germline Division of Labor &
//! Mitochondrial Apoptosis alongside Anaerobic Reductive Evolution of Mitochondria into
//! Hydrogenosomes, Mitosomes (`ISC` Fe-S core) & Cytosolic-`SUF` Amiteukaryotes (`MRO -> 0`),
//! Reductive Apicoplasts (`MEP/IPP + SUF`), Mixotrophic Kleptoplasty & Shopping-Bag `HGT` Tipping
//! Points, Serial `1°/2°/3°` Endosymbiosis, Periplastidial Nucleomorph Elimination, Nuclear `EGT`
//! & CoRR Core Retention, Reticulate Phylogenetic Networks, Spatial Reaction-Diffusion Exergy,
//! Landauer Thermodynamics, Price Equation, and `F_ST`.

use std::collections::VecDeque;
use std::fmt;

/// Deterministic Xoshiro256** PRNG.
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn from_seed(seed: u64) -> Self {
        let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut s = [0u64; 4];
        for slot in &mut s {
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            *slot = z ^ (z >> 31);
            if *slot == 0 {
                *slot = 0xCAFE_BABE_DEAD_BEEF;
            }
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
        }
        Self { s }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = (self.s[1].wrapping_mul(5)).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    #[inline]
    pub fn range_f64(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.next_f64()
    }

    #[inline]
    pub fn gen_index(&mut self, len: usize) -> usize {
        if len <= 1 {
            0
        } else {
            (self.next_u64() as usize) % len
        }
    }

    #[inline]
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.next_f64() < p.clamp(0.0, 1.0)
    }
}

/// Evolutionary continuum of **Mitochondrion-Related Organelles (`MRO`s)** across oxygen gradients
/// (Müller et al., 2012; Stairs, Leger, & Roger, 2015; Karnkowska et al., 2016; Tsaousis, 2019):
/// - `None`: Pre-endosymbiotic / basal non-mitochondriate lineage.
/// - `AerobicMitochondrion` (Class 1): Cristate mitochondrion with proton-pumping `OXPHOS` ETC + `mtDNA` (`CoRR>=3`).
/// - `AnaerobicMitochondrion` (Class 2/3, *Blastocystis*): Facultative fumarate/rhodoquinone MRO retaining vestigial `mtDNA`.
/// - `Hydrogenosome` (Class 4, *Trichomonas* / chytrids): Genome-less (`mtDNA=0`) anaerobic organelle generating `ATP`
///   via substrate-level phosphorylation and evolving `H2` + acetate via `PFO` + `[FeFe]`-hydrogenase!
/// - `Mitosome` (Class 5, *Giardia* / *Entamoeba* / microsporidia): Genome-less (`mtDNA=0`), non-ATP-generating
///   minimalist double-membrane vesicle retained solely for essential `ISC` (`[2Fe-2S]`) iron-sulfur cluster assembly!
/// - `Amiteukaryote` (*Monocercomonoides exilis*): Complete secondary loss of the mitochondrial organelle (`MRO -> 0`)
///   unlocked only when bacterial cytosolic `SUF` (`SufSU/SufB/SufC`) lateral gene transfer replaces `ISC`!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum MroClass {
    #[default]
    None = 0,
    AerobicMitochondrion = 1,
    AnaerobicMitochondrion = 2,
    Hydrogenosome = 3,
    Mitosome = 4,
    Amiteukaryote = 5,
}

impl MroClass {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::None => "None",
            Self::AerobicMitochondrion => "AeroMito",
            Self::AnaerobicMitochondrion => "AnaeroMito",
            Self::Hydrogenosome => "Hydro",
            Self::Mitosome => "MitoS",
            Self::Amiteukaryote => "Amiteuk",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::None => "No MRO",
            Self::AerobicMitochondrion => "Class 1 Aerobic Mitochondrion (OXPHOS + mtDNA)",
            Self::AnaerobicMitochondrion => "Class 2/3 Anaerobic Mitochondrion (RQ + vestigial mtDNA)",
            Self::Hydrogenosome => "Class 4 Hydrogenosome (PFO/HydA H2 + ATP + ISC, mtDNA=0)",
            Self::Mitosome => "Class 5 Mitosome (ISC Fe-S Core Only, ATP=0, mtDNA=0)",
            Self::Amiteukaryote => "Cytosolic-SUF Amiteukaryote (MRO->0, Monocercomonoides)",
        }
    }

    pub fn is_reduced_anaerobic_mro(self) -> bool {
        matches!(
            self,
            Self::AnaerobicMitochondrion | Self::Hydrogenosome | Self::Mitosome
        )
    }
}

/// Evolutionary transition in individuality from unicellular chimeric holobionts to differentiated
/// clonal multicellular organisms (Ratcliff et al., 2012; Radzvilavicius et al., 2016; Aktipis et al., 2015;
/// Arias Del Angel et al., 2020; König & Nedelcu, 2023):
/// - `Unicellular`: Solitary cell (`cell_count = 1`).
/// - `AggregativeChimerism`: Non-clonal "coming together" aggregate (`r < 1.0`, vulnerable to horizontal cheaters).
/// - `UndifferentiatedClonal`: Clonal "staying together" cluster via incomplete cytokinesis (`r = 1.0`, baseline scission apoptosis).
/// - `GermSomaMulticellular`: Differentiated clonal organism (`Volvox` / metazoan grade) with sequestered quiescent
///   germline, disposable foraging/defensive soma (`soma_fraction`), zygotic `N=1` bottleneck, and mitochondrial apoptosis!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum MulticellularState {
    #[default]
    Unicellular = 0,
    AggregativeChimerism = 1,
    UndifferentiatedClonal = 2,
    GermSomaMulticellular = 3,
}

impl MulticellularState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::Unicellular => "Uni",
            Self::AggregativeChimerism => "Aggr",
            Self::UndifferentiatedClonal => "Clonal",
            Self::GermSomaMulticellular => "GermSoma",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Unicellular => "Unicellular Holobiont / Protist (N=1)",
            Self::AggregativeChimerism => "Aggregative Chimeric Mound (Coming Together, r<1)",
            Self::UndifferentiatedClonal => "Undifferentiated Clonal Cluster (Staying Together, r=1)",
            Self::GermSomaMulticellular => "Clonal Germ-Soma Organism (regA + Quiescent Germline + Apoptosis)",
        }
    }

    pub fn is_multicellular(self) -> bool {
        !matches!(self, Self::Unicellular)
    }

    pub fn is_clonal(self) -> bool {
        matches!(self, Self::UndifferentiatedClonal | Self::GermSomaMulticellular)
    }
}

/// Evolutionary transition in individuality from differentiated germ-soma multicellular organisms
/// to **Eusocial Superorganisms** (Boomsma, 2009; West, Fisher, Gardner, & Kiers, 2015; Qiu et al., 2018;
/// Dai et al., 2024; Nonacs, 2006):
/// - `Solitary`: Independent germ-soma multicellular (or unicellular) individual (`colony_size = 1`).
/// - `SubsocialFamily`: Ancestral closed maternal-care brood nest (`colony_size = 6`).
/// - `FacultativeEusocial`: Behaviorally sterile helper caste (*Polistes* / halictid grade, `colony_size = 18`);
///   reversible unless strict lifetime monogamy (`r_sib ~ 0.50`) unlocks the morphological point of no return!
/// - `ObligateSuperorganism`: Irreversibly canalized **queen colony germline** + morphologically **sterile foraging
///   worker & defensive soldier castes** (`colony_size = 64`), sustained against rebel worker egg-laying by **worker kin-policing**!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum EusocialState {
    #[default]
    Solitary = 0,
    SubsocialFamily = 1,
    FacultativeEusocial = 2,
    ObligateSuperorganism = 3,
}

impl EusocialState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::Solitary => "Solitary",
            Self::SubsocialFamily => "Subsoc",
            Self::FacultativeEusocial => "FacEusoc",
            Self::ObligateSuperorganism => "Superorg",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Solitary => "Solitary Organism (Colony=1)",
            Self::SubsocialFamily => "Subsocial Maternal-Care Family (Colony=6)",
            Self::FacultativeEusocial => "Facultative Eusocial Colony (Reversible Helpers, Colony=18)",
            Self::ObligateSuperorganism => {
                "Obligate Superorganism (Queen Germline + Sterile Worker/Soldier Castes, Colony=64)"
            }
        }
    }

    pub fn is_social(self) -> bool {
        !matches!(self, Self::Solitary)
    }

    pub fn is_eusocial(self) -> bool {
        matches!(self, Self::FacultativeEusocial | Self::ObligateSuperorganism)
    }
}

/// Multi-kingdom **Agricultural Ectosymbiosis & Crop Domestication** grade in eusocial colonies
/// (Schultz & Brady, 2008; Branstetter et al., 2017; Kooij et al., 2015; Leal-Dutra et al., 2023;
/// Li et al., 2018; Andersen et al., 2015):
/// - `HunterForager`: Non-farming eusocial or solitary lineage (`cultivar_polyploidy = 1`).
/// - `LowerFacultativeFarming`: Lower-attine detritus cultivation of facultatively free-living diploid (`2n`) fungi.
/// - `HigherDomesticatedCrop`: Higher-attine obligate clonal crop monoculture (`4n`) evolving autophagic
///   **gongylidia** nutritional bulbs, vertical gyne-pellet inheritance, and host metabolic streamlining.
/// - `LeafCutterAgroEcosystem`: *Atta*/*Acromyrmex*—*Leucoagaricus gongylophorus* (`8n`)—*Pseudonocardia*
///   herbivorous agro-ecosystem converting fresh plant cellulose into gongylidia while cuticular actinobacteria
///   and sterile-worker weeding suppress *Escovopsis* crop mycoparasites!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AgricultureState {
    #[default]
    HunterForager = 0,
    LowerFacultativeFarming = 1,
    HigherDomesticatedCrop = 2,
    LeafCutterAgroEcosystem = 3,
}

impl AgricultureState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::HunterForager => "Forager",
            Self::LowerFacultativeFarming => "LowerAg",
            Self::HigherDomesticatedCrop => "GongylCrop",
            Self::LeafCutterAgroEcosystem => "LeafCutter",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::HunterForager => "Hunter-Forager (No Cultivated Crop)",
            Self::LowerFacultativeFarming => {
                "Lower Facultative Farming (2n Detritus Cultivar Pool)"
            }
            Self::HigherDomesticatedCrop => {
                "Higher Domesticated Clonal Crop (4n Gongylidia Monoculture)"
            }
            Self::LeafCutterAgroEcosystem => {
                "Leaf-Cutter Obligate Agro-Ecosystem (8n Polyploid Crop + Pseudonocardia Shield)"
            }
        }
    }

    pub fn is_farming(self) -> bool {
        !matches!(self, Self::HunterForager)
    }

    pub fn is_obligate_monoculture(self) -> bool {
        matches!(
            self,
            Self::HigherDomesticatedCrop | Self::LeafCutterAgroEcosystem
        )
    }
}

/// **Cumulative Cultural Evolution & Dual Genetic–Cultural Inheritance** grade across structured populations
/// (Lewis & Laland, 2012; Derex & Boyd, 2016; Legare, 2017; Muthukrishna et al., 2018; Whitehead et al., 2019):
/// - `AsocialTrialError`: Purely genetic inheritance + non-transmitted individual learning (`ratchet_tier = 0`).
/// - `LowFidelityTradition`: Stimulus enhancement / low-fidelity emulation (`fidelity < 0.82`, `ratchet_tier = 1`);
///   generates local animal traditions that turn over without multi-generational ratcheting.
/// - `ImitativeRatchet`: High-resolution procedural imitation (`fidelity >= 0.82`, `ratchet_tier = 2`) crossing
///   Lewis & Laland's (2012) fidelity threshold and initiating gene-culture neural expansion.
/// - `SymbolicCumulativeCulture`: Compositional **symbolic language + active pedagogical teaching**
///   (`fidelity >= 0.88`, `ratchet_tier = 3..=6`) across partially connected structured demes, unlocking
///   cross-deme cultural trait recombination and dual-channel genetic-cultural mutual information!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CulturalState {
    #[default]
    AsocialTrialError = 0,
    LowFidelityTradition = 1,
    ImitativeRatchet = 2,
    SymbolicCumulativeCulture = 3,
}

impl CulturalState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::AsocialTrialError => "Asocial",
            Self::LowFidelityTradition => "Tradition",
            Self::ImitativeRatchet => "Ratchet",
            Self::SymbolicCumulativeCulture => "SymbolicCCE",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::AsocialTrialError => "Asocial Trial-and-Error (Genetic Channel Only)",
            Self::LowFidelityTradition => "Low-Fidelity Social Tradition (Non-Cumulative, Tier 1)",
            Self::ImitativeRatchet => "High-Fidelity Imitative Ratchet (Cumulative Culture, Tier 2)",
            Self::SymbolicCumulativeCulture => {
                "Symbolic Pedagogical Cumulative Culture (Dual Inheritance + Cross-Deme Recomb, Tier 4-6)"
            }
        }
    }

    pub fn is_cultural(self) -> bool {
        !matches!(self, Self::AsocialTrialError)
    }

    pub fn is_cumulative_ratchet(self) -> bool {
        matches!(self, Self::ImitativeRatchet | Self::SymbolicCumulativeCulture)
    }
}

/// **Externalized Exosomatic Memory, Written Archives & Institutional Peer-Verification** grade
/// (Pegado, 2022; Grinschgl et al., 2021; González-Forero, 2025; Smaldino & McElreath, 2016; Tiokhin et al., 2019):
/// - `EndosomaticBound`: Cultural repertoire stored exclusively in living neural tissue (`archive = 0`).
/// - `MaterialNotation`: External tallies, tokens & physical mnemonics (`8` external loci, `12%` offloading).
/// - `WrittenArchive`: Durable orthographic script & manuscript libraries (`22` external loci, `28%` cognitive
///   offloading); decouples storage from individual brain size and rescues isolated demes from Tasmanian drift,
///   but remains vulnerable to unverified false-discovery pollution!
/// - `PeerVerifiedInstitution`: Institutionalized replication, open audit & peer-verification ratchets
///   (`38..=48` verified external loci, `40%` cognitive offloading); purges false discoveries while streamlining
///   individual metabolic brain cost (`neural_learning_capacity: 0.90 -> 0.68`)!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ExosomaticState {
    #[default]
    EndosomaticBound = 0,
    MaterialNotation = 1,
    WrittenArchive = 2,
    PeerVerifiedInstitution = 3,
}

impl ExosomaticState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::EndosomaticBound => "Endosomatic",
            Self::MaterialNotation => "Notation",
            Self::WrittenArchive => "WrittenArch",
            Self::PeerVerifiedInstitution => "VerifiedInst",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::EndosomaticBound => "Endosomatic Neural Storage Only (Brain-Bound)",
            Self::MaterialNotation => "Material Tallies & Mnemonic Notation (8 Exosomatic Loci)",
            Self::WrittenArchive => "Durable Written Archive (22 Loci, Brain-Decoupled)",
            Self::PeerVerifiedInstitution => {
                "Peer-Verified Epistemic Institution (38-48 Verified Loci + Replication Ratchet)"
            }
        }
    }

    pub fn is_archival(self) -> bool {
        !matches!(self, Self::EndosomaticBound)
    }

    pub fn has_written_archive(self) -> bool {
        matches!(self, Self::WrittenArchive | Self::PeerVerifiedInstitution)
    }
}

/// **Autonomous Machine-Verified Exosomatic Synthesis & Computational Symbiont** grade
/// (Shumailov et al., 2024; Romera-Paredes et al., 2024; Trinh et al., 2024; PMC13487184):
/// - `ManualArchival`: Purely biological curation of written/peer-reviewed archives (`0` machine loci).
/// - `AlgorithmicMediation`: Search & recommender indexing accelerating cross-archive retrieval (`10` loci).
/// - `UnverifiedRecursiveGen`: High-throughput synthetic generation without formal proof-checking or
///   empirical provenance; vulnerable to **Shumailov (2024) Recursive Model Collapse** (tail loss -> delta mode)!
/// - `MachineVerifiedSymbiosis`: Distributed **FunSearch / AlphaGeometry** loop coupling broad-horizon
///   generative samplers with deterministic formal evaluators (`formal_verifier_strictness >= 0.72`) and
///   empirical grounding anchors (`empirical_grounding_floor >= 0.12`), accumulating `36..=64` formally
///   verified theorems/programs with preserved distribution tails!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum MachineSynthesisState {
    #[default]
    ManualArchival = 0,
    AlgorithmicMediation = 1,
    UnverifiedRecursiveGen = 2,
    MachineVerifiedSymbiosis = 3,
}

impl MachineSynthesisState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::ManualArchival => "ManualArch",
            Self::AlgorithmicMediation => "AlgoMed",
            Self::UnverifiedRecursiveGen => "UnverifiedGen",
            Self::MachineVerifiedSymbiosis => "FormalSymbio",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::ManualArchival => "Manual Biological Archival Curation Only",
            Self::AlgorithmicMediation => "Algorithmic Search & Recommender Mediation (10 Loci)",
            Self::UnverifiedRecursiveGen => {
                "Unverified Recursive Synthetic Generation (Vulnerable to Model Collapse!)"
            }
            Self::MachineVerifiedSymbiosis => {
                "Machine-Verified Neuro-Symbolic Symbiosis (Formal Proof Kernel + Empirical Anchor)"
            }
        }
    }

    pub fn has_computational_symbiont(self) -> bool {
        !matches!(self, Self::ManualArchival)
    }

    pub fn is_formally_verified(self) -> bool {
        matches!(self, Self::MachineVerifiedSymbiosis)
    }
}

/// **Planetary-Scale Biogeochemical & Technosphere Homeostasis** epoch
/// (Kleidon, 2012; Krissansen-Totton et al., 2018; Wong & Bartlett, 2022; Richardson et al., 2023; PMC11848429):
/// - `AnoxicArchean`: Reducing/hypoxic atmosphere (`O2 < 0.25`, high `CO2/CH4`), dominated by anaerobic syntrophs & MROs.
/// - `OxygenatedPhanerozoic`: Autotroph + syntroph + superorganism co-regulated `O2 >= 0.55`, closed `N/P` cycles,
///   and stable Holocene-like carbon/exergy balance.
/// - `AnthropoceneDisequilibrium`: Superlinear unguided technospheric extraction breaches safe planetary boundaries
///   (`CO2 > 700 ppm`, `boundary_transgression_load >= 0.45`), risking **Asymptotic Planetary Burnout**!
/// - `HomeostaticAwakening`: Machine-verified planetary stewardship (`Gaia 2.0`) couples formal Earth-system models
///   with closed-loop carbon drawdown and renewable solar/geothermal exergy transduction!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum PlanetaryState {
    AnoxicArchean = 0,
    #[default]
    OxygenatedPhanerozoic = 1,
    AnthropoceneDisequilibrium = 2,
    HomeostaticAwakening = 3,
}

impl PlanetaryState {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::AnoxicArchean => "Archean",
            Self::OxygenatedPhanerozoic => "Phanerozoic",
            Self::AnthropoceneDisequilibrium => "BoundaryBreach",
            Self::HomeostaticAwakening => "Gaia2Homeostasis",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::AnoxicArchean => "Anoxic Archean Biosphere (Hypoxic O2 < 0.25, High CO2/CH4)",
            Self::OxygenatedPhanerozoic => {
                "Oxygenated Phanerozoic Biosphere (Autotroph-Syntroph-Superorg Regulated)"
            }
            Self::AnthropoceneDisequilibrium => {
                "Planetary Boundary Disequilibrium (Carbon/Nutrient Overshoot — Burnout Risk!)"
            }
            Self::HomeostaticAwakening => {
                "Gaia-2.0 Homeostatic Awakening (Machine-Verified Carbon/Exergy Co-Regulation)"
            }
        }
    }

    pub fn is_oxygenated(self) -> bool {
        !matches!(self, Self::AnoxicArchean)
    }

    pub fn is_homeostatic_awakened(self) -> bool {
        matches!(self, Self::HomeostaticAwakening)
    }
}

/// Prognostic **Planetary Atmosphere, Nutrient & Free-Energy Reservoir** co-regulated by biological
/// and machine-verified guilds (Kleidon, 2012; Krissansen-Totton et al., 2018; Wong & Bartlett, 2022).
/// Note: `atmospheric_o2` tracks a normalized `O2 (PAL Index)` (Present Atmospheric Level index in `[0.02, 0.94]`).
#[derive(Clone, Debug)]
pub struct PlanetaryBiosphere {
    pub state: PlanetaryState,
    /// Normalized `O2 (PAL Index)` (Present Atmospheric Level index; `[0.02, 0.94]` dynamic range).
    pub atmospheric_o2: f64,
    pub atmospheric_co2_ppm: f64,
    pub nutrient_pool_np: f64,
    pub planetary_exergy_budget: f64,
    pub chemical_disequilibrium_bits: f64,
    pub boundary_transgression_load: f64,
    pub goe_crossed: bool,
    pub homeostatic_awakening_crossed: bool,
}

impl PlanetaryBiosphere {
    pub fn new_for_controls(controls: &EnvironmentControls) -> Self {
        let o2 = controls.ambient_oxygen_level.clamp(0.02, 0.98);
        let init_state = if o2 < 0.25 {
            PlanetaryState::AnoxicArchean
        } else {
            PlanetaryState::OxygenatedPhanerozoic
        };
        let co2 = if o2 < 0.25 { 1180.0 } else { 320.0 };
        let diseq = (2.4 + 14.5 * o2).clamp(1.0, 22.0);
        Self {
            state: init_state,
            atmospheric_o2: o2,
            atmospheric_co2_ppm: co2,
            nutrient_pool_np: 0.82,
            planetary_exergy_budget: controls.bath_exergy_flux.max(0.45),
            chemical_disequilibrium_bits: diseq,
            boundary_transgression_load: 0.0,
            goe_crossed: o2 >= 0.25,
            homeostatic_awakening_crossed: false,
        }
    }

    /// Advances planetary `O2`, `CO2`, `N/P` nutrients, chemical disequilibrium, and planetary boundary
    /// load by 1 tick from aggregated guild fluxes; returns `(new_state, burnout_collapse_triggered)`.
    pub fn step_biogeochemical_feedbacks(
        &mut self,
        autotroph_flux: f64,
        syntroph_recycling_flux: f64,
        superorg_weathering_flux: f64,
        unguided_techno_flux: f64,
        verified_stewardship_flux: f64,
        controls: &EnvironmentControls,
    ) -> (PlanetaryState, bool) {
        let gain = controls.biogeochemical_feedback_gain.clamp(0.0, 1.0);
        if gain <= 0.001 {
            self.atmospheric_o2 = controls.ambient_oxygen_level.clamp(0.02, 0.98);
            self.planetary_exergy_budget = controls.bath_exergy_flux.max(0.25);
            return (self.state, false);
        }

        let steward_c = controls.planetary_stewardship_coupling.clamp(0.0, 1.0);
        let burn_i = controls.extractive_burn_intensity.clamp(0.0, 1.0);

        // 1. Syntrophic + Superorganismal + Stewardship N/P Nutrient Loop Closure vs Autotrophic Uptake
        let nutrient_uptake = 0.055 * autotroph_flux.clamp(0.0, 2.5);
        let nutrient_return = 0.115 * syntroph_recycling_flux.clamp(0.0, 2.5)
            + 0.065 * superorg_weathering_flux.clamp(0.0, 2.5)
            + 0.080 * steward_c * verified_stewardship_flux.clamp(0.0, 2.5);
        self.nutrient_pool_np =
            (self.nutrient_pool_np + gain * (nutrient_return - nutrient_uptake)).clamp(0.08, 1.0);

        // 2. Atmospheric O2 Production (Nutrient-Gated Autotrophy) vs Respiration & Extractive Combustion
        let eff_autotrophy = autotroph_flux.clamp(0.0, 3.0) * (0.25 + 0.75 * self.nutrient_pool_np);
        let o2_source = 0.085 * eff_autotrophy + 0.035 * steward_c * verified_stewardship_flux;
        let o2_sink = 0.028
            + 0.045 * self.atmospheric_o2
            + 0.095 * burn_i * unguided_techno_flux * (1.0 - 0.75 * steward_c);
        self.atmospheric_o2 =
            (self.atmospheric_o2 + gain * (o2_source - o2_sink)).clamp(0.02, 0.94);

        if self.atmospheric_o2 >= 0.25 {
            self.goe_crossed = true;
        }

        // 3. Atmospheric CO2 Budget: Volcanic/Extractive Emissions vs Photosynthetic + Silicate Weathering + DAC Drawdown
        let emissions_ppm = 14.0
            + 135.0 * burn_i * (0.35 + unguided_techno_flux) * (1.0 - 0.82 * steward_c);
        let drawdown_ppm = 38.0 * eff_autotrophy
            + 26.0 * superorg_weathering_flux
            + 95.0 * steward_c * (0.25 + verified_stewardship_flux);
        self.atmospheric_co2_ppm =
            (self.atmospheric_co2_ppm + gain * (emissions_ppm - drawdown_ppm)).clamp(220.0, 1800.0);

        // 4. Planetary Boundaries Load (Richardson et al., 2023) & Wong-Bartlett Burnout vs Awakening
        let co2_excess = ((self.atmospheric_co2_ppm - 420.0).max(0.0) / 650.0).clamp(0.0, 1.2);
        let nutrient_disrupt = (0.35 - self.nutrient_pool_np).max(0.0) * 0.85;
        let boundary_influx = 0.24 * co2_excess + nutrient_disrupt;
        let stewardship_repair = steward_c * (0.14 + 0.22 * verified_stewardship_flux);
        self.boundary_transgression_load = (self.boundary_transgression_load
            + gain * (boundary_influx - stewardship_repair))
            .clamp(0.0, 1.0);

        // 5. Gibbs Chemical Disequilibrium Biosignature (Krissansen-Totton et al., 2018) & Planetary Exergy Budget (Kleidon, 2012)
        self.chemical_disequilibrium_bits = (2.2
            + 15.6 * self.atmospheric_o2 * self.nutrient_pool_np
            + 4.8 * steward_c * verified_stewardship_flux.min(1.5))
        .clamp(1.2, 26.0);

        let thermal_penalty = (1.0 - 0.72 * self.boundary_transgression_load).clamp(0.22, 1.0);
        let solar_biotic_transduction = 0.72
            + 0.48 * self.nutrient_pool_np * (0.4 + 0.6 * self.atmospheric_o2)
            + 0.52 * steward_c * verified_stewardship_flux.min(1.4);
        self.planetary_exergy_budget = (controls.bath_exergy_flux.max(0.60)
            * solar_biotic_transduction
            * thermal_penalty)
            .clamp(0.22, 2.35);

        self.state = if self.boundary_transgression_load >= 0.42 {
            PlanetaryState::AnthropoceneDisequilibrium
        } else if verified_stewardship_flux >= 0.25
            && steward_c >= 0.65
            && self.atmospheric_o2 >= 0.55
            && self.atmospheric_co2_ppm <= 380.0
        {
            self.homeostatic_awakening_crossed = true;
            PlanetaryState::HomeostaticAwakening
        } else if self.atmospheric_o2 < 0.25 {
            PlanetaryState::AnoxicArchean
        } else {
            PlanetaryState::OxygenatedPhanerozoic
        };

        let burnout = self.boundary_transgression_load >= 0.70;
        (self.state, burnout)
    }
}

/// **Early-Earth Prebiotic Geochemical Cradle** (Damer & Deamer, 2020; Mulkidjanian et al., 2012;
/// Herschy et al., 2014; Lane, 2014; Pearce et al., 2020):
/// - `HotSpringWetDry`: Subaerial volcanic geothermal pools cycling through dry multilamellar lipid films,
///   hydrated protocell budding, and moist hydrogel phases with high `K+/Na+ > 1.0`, `Zn2+`, and phosphate.
/// - `AlkalineHydrothermalVent`: Submarine serpentinization microporous chimneys with `Fe(Ni)S` (greigite/mackinawite)
///   catalytic walls and natural `Delta-pH ~ 3.5` gradients driving `H2 + CO2 -> Acetyl-CoA` fixation.
/// - `HybridHydrothermalContinent`: Coupled volcanic island hydrology uniting wet-dry condensation polymers
///   and chelated `Mg2+` lipid vesicles with hydrothermal `FeS` redox cofactors!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AbiogenesisCradle {
    HotSpringWetDry = 0,
    AlkalineHydrothermalVent = 1,
    #[default]
    HybridHydrothermalContinent = 2,
}

impl AbiogenesisCradle {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::HotSpringWetDry => "HotSpringWetDry",
            Self::AlkalineHydrothermalVent => "AlkalineVentFeS",
            Self::HybridHydrothermalContinent => "CoupledHotSpringVent",
        }
    }
}

/// **Stepwise Hadean Abiogenesis Phase Transition** (`~4.2 - 4.0 Ga`):
/// - `PrebioticMonomerSoup`: Abiotic `CO2/N2/H2/H2S/NH3/PO4` monomer pool (`L <= 10 nt`).
/// - `AutocatalyticRafProtocell`: Amphiphilic vesicle encapsulating a Reflexively Autocatalytic
///   Food-Generated (`RAF`) reaction network stabilized by citrate/peptide `Mg2+` chelation
///   (Hordijk et al., 2018; Xavier et al., 2020; Adamala & Szostak, 2013).
/// - `RnaRibozymeHypercycle`: Cooperative polymerase + ligase ribozyme hypercycle overcoming
///   Eigen's Error Catastrophe (`L_oligomer < ln(sigma) / mu`, PMC4824418; Papastavrou et al., 2024).
/// - `ProgenoteCommunity`: Woese's communal proto-ribosomal RNA-peptide network evolving via
///   pervasive horizontal gene/innovation sharing (`HGT`) before vertical lineage crystallization.
/// - `LucaCellularLineage`: Crosses **Woese's Darwinian Threshold** into vertical RNP/DNA heredity,
///   universal genetic code, `K+/Na+` ion-tight membranes, and rotor ATPase (`LUCA`)!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AbiogenesisStage {
    PrebioticMonomerSoup = 0,
    AutocatalyticRafProtocell = 1,
    RnaRibozymeHypercycle = 2,
    ProgenoteCommunity = 3,
    #[default]
    LucaCellularLineage = 4,
}

impl AbiogenesisStage {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::PrebioticMonomerSoup => "MonomerSoup",
            Self::AutocatalyticRafProtocell => "RAFProtocell",
            Self::RnaRibozymeHypercycle => "RibozymeHypercycle",
            Self::ProgenoteCommunity => "ProgenotePool",
            Self::LucaCellularLineage => "LUCA_Cellular",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::PrebioticMonomerSoup => "Hadean Prebiotic Monomer Soup (Abiotic C/H/O/N/P/S)",
            Self::AutocatalyticRafProtocell => {
                "Encapsulated RAF Protocell (Citrate-Mg2+ Lipid Vesicle + FeS Cofactors)"
            }
            Self::RnaRibozymeHypercycle => {
                "Cooperative RNA Ribozyme Hypercycle (Above Eigen Error Cliff)"
            }
            Self::ProgenoteCommunity => {
                "Woese Communal Progenote Network (Horizontal Proto-Ribosomal Pool)"
            }
            Self::LucaCellularLineage => {
                "LUCA Vertical Cellular Lineage (Darwinian Threshold Crossed)"
            }
        }
    }

    pub fn has_raf_closure(self) -> bool {
        (self as u8) >= (Self::AutocatalyticRafProtocell as u8)
    }

    pub fn has_ribozyme_heredity(self) -> bool {
        (self as u8) >= (Self::RnaRibozymeHypercycle as u8)
    }

    pub fn crossed_darwinian_threshold(self) -> bool {
        matches!(self, Self::LucaCellularLineage)
    }
}

/// Prognostic **Hadean Geochemical & Prebiotic Reactor State** (~4.2 Ga Early Earth).
#[derive(Clone, Debug)]
pub struct HadeanGeochemicalBath {
    pub cradle: AbiogenesisCradle,
    /// Hydration-dehydration-gel phase in `[0.0, 1.0]` (`<0.33` Dry Multilamellar, `0.33..0.67` Hydrated Budding, `>=0.67` Moist Gel).
    pub wet_dry_phase: f64,
    pub wet_dry_amplitude: f64,
    /// Submarine/geothermal electrochemical proton gradient `Delta-pH` (`0.0..4.5`, ~210 mV at 3.5 pH units).
    pub vent_delta_ph: f64,
    /// Geothermal `K+/Na+` molar ratio (`~0.03` bulk Hadean ocean vs `1.2..2.4` anoxic geothermal vapor condensate; Mulkidjanian et al., 2012).
    pub k_to_na_ratio: f64,
    /// Catalytic divalent `Mg2+` concentration (`mM`, typically `50 mM` for template-directed RNA extension).
    pub mg2_molar_mm: f64,
    /// Citrate / short acidic peptide `Mg2+` chelation fraction `[0.0, 1.0]` (Adamala & Szostak, 2013).
    pub citrate_chelation_fraction: f64,
    /// Per-nucleotide replication error probability `mu` (`0.04` non-enzymatic -> `0.006` ribozyme -> `0.0006` RNP/DNA).
    pub per_base_error_mu: f64,
    /// Reflexively Autocatalytic Food-Generated (`RAF`) catalytic closure `[0.0, 1.0]`.
    pub raf_catalytic_closure: f64,
    /// Eigen error-threshold ceiling `L_max = floor(ln(sigma) / mu)` in nucleotides (`nt`).
    pub eigen_max_length_nt: usize,
    pub darwinian_threshold_crossed: bool,
}

impl HadeanGeochemicalBath {
    pub fn new_for_controls(controls: &EnvironmentControls) -> Self {
        let wd = controls.wet_dry_cycle_amplitude.clamp(0.0, 1.0);
        let dph = controls.hydrothermal_delta_ph.clamp(0.0, 4.5);
        let chelat = controls.citrate_mg_chelation.clamp(0.0, 1.0);
        let cradle = if wd >= 0.45 && dph >= 1.8 {
            AbiogenesisCradle::HybridHydrothermalContinent
        } else if wd >= dph / 4.5 {
            AbiogenesisCradle::HotSpringWetDry
        } else {
            AbiogenesisCradle::AlkalineHydrothermalVent
        };
        let k_na = (0.08 + 1.85 * wd + 0.15 * dph).clamp(0.03, 2.50);
        let mu = (0.055 * (1.0 - 0.96 * controls.ribozyme_polymerase_fidelity.clamp(0.0, 1.0)))
            .clamp(0.0006, 0.065);
        let sigma = (1.65 + 1.85 * wd + 0.42 * dph).clamp(1.25, 4.8);
        let eigen_l = ((sigma.ln() / mu).floor() as usize).clamp(12, 4800);
        let raf = (0.28 * wd + 0.14 * dph + 0.25 * chelat).clamp(0.05, 0.99);
        Self {
            cradle,
            wet_dry_phase: 0.0,
            wet_dry_amplitude: wd,
            vent_delta_ph: dph,
            k_to_na_ratio: k_na,
            mg2_molar_mm: 50.0,
            citrate_chelation_fraction: chelat,
            per_base_error_mu: mu,
            raf_catalytic_closure: raf,
            eigen_max_length_nt: eigen_l,
            darwinian_threshold_crossed: !controls.hadean_abiogenesis_mode,
        }
    }

    /// Advances the Hadean geochemical reactor by 1 hydration-dehydration-gel & hydrothermal redox cycle!
    pub fn step_geochemical_cycle(
        &mut self,
        controls: &EnvironmentControls,
        mean_raf_feedback: f64,
        mean_ribozyme_fraction: f64,
    ) {
        self.wet_dry_amplitude = controls.wet_dry_cycle_amplitude.clamp(0.0, 1.0);
        self.vent_delta_ph = controls.hydrothermal_delta_ph.clamp(0.0, 4.5);
        self.citrate_chelation_fraction = controls.citrate_mg_chelation.clamp(0.0, 1.0);
        self.wet_dry_phase = (self.wet_dry_phase + 0.34) % 1.0;
        self.k_to_na_ratio =
            (0.08 + 1.85 * self.wet_dry_amplitude + 0.15 * self.vent_delta_ph).clamp(0.03, 2.50);

        let geochem_drive = 0.42 * self.wet_dry_amplitude + 0.13 * self.vent_delta_ph;
        let target_raf = (geochem_drive
            + 0.26 * self.citrate_chelation_fraction
            + 0.28 * mean_raf_feedback.clamp(0.0, 1.0))
        .clamp(0.05, 0.99);
        self.raf_catalytic_closure =
            (0.65 * self.raf_catalytic_closure + 0.35 * target_raf).clamp(0.05, 0.99);

        let fid = (controls.ribozyme_polymerase_fidelity.clamp(0.0, 1.0) * 0.75
            + 0.25 * mean_ribozyme_fraction.clamp(0.0, 1.0))
        .clamp(0.0, 0.995);
        self.per_base_error_mu = (0.055 * (1.0 - 0.985 * fid)).clamp(0.0006, 0.065);
        let sigma = (1.60
            + 1.75 * self.wet_dry_amplitude
            + 0.40 * self.vent_delta_ph
            + 1.40 * self.raf_catalytic_closure)
        .clamp(1.25, 5.50);
        self.eigen_max_length_nt =
            ((sigma.ln() / self.per_base_error_mu).floor() as usize).clamp(12, 6000);
    }

    /// Free-energy condensation & redox power multiplier relative to hydrolytic bulk water!
    pub fn condensation_exergy_drive(&self) -> f64 {
        let phase_pulse = if self.wet_dry_phase < 0.34 {
            1.28 // Dry multilamellar ester/phosphodiester condensation film
        } else if self.wet_dry_phase < 0.67 {
            1.08 // Hydrated amphiphilic vesicle budding
        } else {
            1.16 // Moist hydrogel combinatorial exchange
        };
        let base = 0.58 * self.wet_dry_amplitude * phase_pulse + 0.16 * self.vent_delta_ph;
        base.clamp(0.04, 1.85)
    }
}

/// **Prebiotic Homochirality, Stereoselective Aminoacylation & Genetic Code Crystallization** stage:
/// - `RacemicCrossInhibited`: Symmetric `50:50 D/L` monomer bath (`|ee_D| < 0.25`); opposite-handed
///   nucleotides cap `3'` primers via Joyce-Orgel enantiomeric cross-inhibition (PMC6120706; PMC6974302).
/// - `HomochiralRnaSymmetryBroken`: Magnetite `Fe3O4` spin-polarized `RAO` crystallization (`CISS`,
///   Ozturk et al., 2023, PMC10246896) + Frank autocatalysis amplifies `ee_D_ribose >= 0.85`.
/// - `StereoselectiveMinihelixAcylated`: Tamura-Schimmel `3'-NCCA` `D-RNA` minihelices (`105 deg`
///   Burgi-Dunitz steric clash, ~4:1 `L:D` preference, PMC1564265; PMC10057131) propagate chirality to
///   `L-aminoacyl` peptides (`ee_L_amino >= 0.80`), enabling folded Class I/II urzyme synthetases!
/// - `AmbiguousStatisticalCode`: Woese statistical proto-ribosome (`10` early prebiotic amino acids,
///   `MP ~ 0.42`), vulnerable to fragmented "Tower of Babel" suboptimal code freezing without `HGT`.
/// - `UniversalErrorMinimizedSgc`: Collective progenote `HGT` optimization (Vetsigian, Woese & Goldenfeld,
///   2006, PMC1502294; Koonin & Novozhilov, 2009, PMC3293468) crystallizes the universal `64-codon / 20-AA`
///   error-minimizing Standard Genetic Code (`MP >= 0.76`, `I(Codon;AA) >= 4.18 bits`)!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ChiralityCodeStage {
    RacemicCrossInhibited = 0,
    HomochiralRnaSymmetryBroken = 1,
    StereoselectiveMinihelixAcylated = 2,
    AmbiguousStatisticalCode = 3,
    #[default]
    UniversalErrorMinimizedSgc = 4,
}

impl ChiralityCodeStage {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::RacemicCrossInhibited => "RacemicInhibited",
            Self::HomochiralRnaSymmetryBroken => "Homochiral_D_RNA",
            Self::StereoselectiveMinihelixAcylated => "Minihelix_L_AA",
            Self::AmbiguousStatisticalCode => "StatisticalCode",
            Self::UniversalErrorMinimizedSgc => "Universal_SGC",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::RacemicCrossInhibited => {
                "Racemic 50:50 D/L Pool (Enantiomeric Cross-Inhibition Stall)"
            }
            Self::HomochiralRnaSymmetryBroken => {
                "Homochiral D-Ribose RNA (Magnetite CISS + Frank Symmetry Breaking)"
            }
            Self::StereoselectiveMinihelixAcylated => {
                "Stereoselective D-RNA Minihelix -> L-Aminoacyl Peptides (4:1 Steric Bias)"
            }
            Self::AmbiguousStatisticalCode => {
                "Woese Statistical Proto-Code (10-AA Ambiguous Translation)"
            }
            Self::UniversalErrorMinimizedSgc => {
                "Crystallized Universal 64-Codon / 20-AA Error-Minimized SGC (MP >= 78%)"
            }
        }
    }

    pub fn is_homochiral_rna(self) -> bool {
        (self as u8) >= (Self::HomochiralRnaSymmetryBroken as u8)
    }

    pub fn has_stereoselective_peptides(self) -> bool {
        (self as u8) >= (Self::StereoselectiveMinihelixAcylated as u8)
    }

    pub fn is_crystallized_sgc(self) -> bool {
        matches!(self, Self::UniversalErrorMinimizedSgc)
    }
}

/// Prognostic **Prebiotic Chiral Symmetry-Breaking, Minihelix Aminoacylation & SGC Optimizer** state.
#[derive(Clone, Debug)]
pub struct ChiralCodeReactor {
    /// Magnetite (`Fe3O4`) spin-polarized `CISS` + Frank mutual antagonism bias `[0.0, 1.0]`.
    pub magnetite_ciss_bias: f64,
    /// Tamura-Schimmel `3'-NCCA` RNA minihelix steric `L-amino acid` selectivity `[0.0, 1.0]`.
    pub minihelix_steric_selectivity: f64,
    /// Vetsigian-Woese-Goldenfeld communal progenote `HGT` code-sharing coupling `[0.0, 1.0]`.
    pub communal_code_hgt_coupling: f64,
    /// Enantiomeric excess of `D-ribose` nucleotides in `[0.0, 1.0]` (`0.0` racemic, `1.0` homochiral `D`).
    pub ribose_ee_d: f64,
    /// Enantiomeric excess of `L-aminoacyl` peptides in `[0.0, 1.0]` (`0.0` racemic, `1.0` homochiral `L`).
    pub amino_ee_l: f64,
    /// Probability of `3'` RNA primer chain-termination by opposite-handed `L-ribonucleotides`!
    pub cross_inhibition_prob: f64,
    /// Active encoded canonical `L-amino acid` alphabet size (`0 -> 4 -> 10 -> 20`).
    pub encoded_amino_acids: u8,
    /// Polar-requirement **Error-Minimization Percentage (`MP`)** in `[0.0, 0.88]` (`SGC ~ 0.78`).
    pub sgc_error_minimization_mp: f64,
    /// Fraction of population sharing a single mutually translatable codon table (`< 0.35` = Babel Fragmentation).
    pub code_universality_consensus: f64,
    /// Shannon mutual information `I(Codon ; AminoAcid)` in bits (`log2(20) = 4.322` max).
    pub translational_mutual_info_bits: f64,
    pub sgc_crystallized: bool,
}

impl ChiralCodeReactor {
    pub fn new_for_controls(controls: &EnvironmentControls) -> Self {
        if !controls.chiral_symmetry_breaking_mode {
            return Self {
                magnetite_ciss_bias: controls.magnetite_ciss_spin_bias.clamp(0.0, 1.0),
                minihelix_steric_selectivity: controls.minihelix_chiral_selectivity.clamp(0.0, 1.0),
                communal_code_hgt_coupling: controls.communal_code_hgt_rate.clamp(0.0, 1.0),
                ribose_ee_d: 0.99,
                amino_ee_l: 0.99,
                cross_inhibition_prob: 0.005,
                encoded_amino_acids: 20,
                sgc_error_minimization_mp: 0.79,
                code_universality_consensus: 0.98,
                translational_mutual_info_bits: 4.26,
                sgc_crystallized: true,
            };
        }
        let ciss = controls.magnetite_ciss_spin_bias.clamp(0.0, 1.0);
        let ster = controls.minihelix_chiral_selectivity.clamp(0.0, 1.0);
        let hgt = controls.communal_code_hgt_rate.clamp(0.0, 1.0);
        let init_ee_d = (0.04 + 0.58 * ciss).clamp(0.02, 0.99);
        let init_ee_l = (init_ee_d * 0.62 * ster).clamp(0.01, 0.98);
        let cross_inh = (0.78 * (1.0 - init_ee_d).powi(2)).clamp(0.005, 0.80);
        Self {
            magnetite_ciss_bias: ciss,
            minihelix_steric_selectivity: ster,
            communal_code_hgt_coupling: hgt,
            ribose_ee_d: init_ee_d,
            amino_ee_l: init_ee_l,
            cross_inhibition_prob: cross_inh,
            encoded_amino_acids: if init_ee_l >= 0.45 { 4 } else { 0 },
            sgc_error_minimization_mp: 0.18,
            code_universality_consensus: (0.20 + 0.76 * hgt).clamp(0.12, 0.99),
            translational_mutual_info_bits: 0.35,
            sgc_crystallized: false,
        }
    }

    /// Advances chiral amplification (`Frank + Magnetite CISS`), Tamura-Schimmel `4:1` minihelix
    /// `D-RNA -> L-AA` chiral transfer, and Vetsigian-Woese communal `SGC` error minimization!
    pub fn step_chiral_and_code_dynamics(&mut self, controls: &EnvironmentControls) {
        self.magnetite_ciss_bias = controls.magnetite_ciss_spin_bias.clamp(0.0, 1.0);
        self.minihelix_steric_selectivity = controls.minihelix_chiral_selectivity.clamp(0.0, 1.0);
        self.communal_code_hgt_coupling = controls.communal_code_hgt_rate.clamp(0.0, 1.0);

        // 1. Frank mutual antagonism + Ozturk-Sasselov Fe3O4 RAO spin-selective crystallization:
        if self.magnetite_ciss_bias < 0.12 {
            // Racemization erosion outweighs tiny chiral fluctuations -> trapped near 50:50 racemic!
            self.ribose_ee_d = (self.ribose_ee_d * 0.84 + 0.03 * self.magnetite_ciss_bias).clamp(0.01, 0.22);
        } else {
            // Bifurcation above Frank critical threshold -> rapid homochiral D-ribose amplification!
            let frank_gain = 0.26 * self.magnetite_ciss_bias + 0.22 * self.ribose_ee_d * (1.0 - self.ribose_ee_d);
            self.ribose_ee_d = (self.ribose_ee_d + frank_gain - 0.012 * (1.0 - self.magnetite_ciss_bias))
                .clamp(0.02, 0.995);
        }
        self.cross_inhibition_prob = (0.82 * (1.0 - self.ribose_ee_d).powi(2)).clamp(0.004, 0.82);

        // 2. Central Dogma of Homochirality (Tamura & Schimmel, 2006: D-RNA minihelix -> 4:1 L-AA steric bias):
        if self.ribose_ee_d >= 0.65 && self.minihelix_steric_selectivity >= 0.25 {
            let target_l = (self.ribose_ee_d * (0.35 + 0.65 * self.minihelix_steric_selectivity)).clamp(0.05, 0.99);
            self.amino_ee_l = (0.58 * self.amino_ee_l + 0.42 * target_l).clamp(0.02, 0.99);
        } else {
            self.amino_ee_l = (self.amino_ee_l * 0.82 + 0.08 * self.minihelix_steric_selectivity * self.ribose_ee_d)
                .clamp(0.01, 0.24);
        }

        // 3. Coevolutionary amino acid recruitment (0 -> 4 GADV -> 10 early abiotic -> 20 SGC) & Error Minimization MP:
        self.encoded_amino_acids = if self.amino_ee_l < 0.40 {
            0
        } else if self.amino_ee_l < 0.72 {
            4
        } else if self.sgc_error_minimization_mp < 0.68 {
            10
        } else {
            20
        };

        self.code_universality_consensus =
            (0.16 + 0.83 * self.communal_code_hgt_coupling).clamp(0.12, 0.99);

        if self.amino_ee_l >= 0.72 {
            // Without communal HGT (`hgt < 0.20`), isolated codes freeze at a rugged local optimum (`MP <= 0.44`)!
            let mp_ceiling = if self.communal_code_hgt_coupling < 0.20 {
                0.34 + 0.45 * self.communal_code_hgt_coupling
            } else {
                0.48 + 0.36 * self.communal_code_hgt_coupling
            };
            self.sgc_error_minimization_mp = (0.64 * self.sgc_error_minimization_mp + 0.36 * mp_ceiling)
                .clamp(0.12, 0.84);
        } else {
            self.sgc_error_minimization_mp = (self.sgc_error_minimization_mp * 0.90).max(0.12);
        }

        let aa_f = (self.encoded_amino_acids.max(1)) as f64;
        let max_bits = if self.encoded_amino_acids == 0 {
            0.25
        } else {
            aa_f.log2()
        };
        self.translational_mutual_info_bits =
            (max_bits * self.amino_ee_l * (0.35 + 0.65 * self.sgc_error_minimization_mp / 0.80))
                .clamp(0.10, 4.32);

        self.sgc_crystallized = self.ribose_ee_d >= 0.90
            && self.amino_ee_l >= 0.85
            && self.encoded_amino_acids == 20
            && self.sgc_error_minimization_mp >= 0.74
            && self.code_universality_consensus >= 0.65;
    }
}

/// **Ribonucleotide Reductase (`RNR`) RNA-to-DNA Transition, `U -> T` Deamination Repair & Virosphere Stage**:
/// - `LabileRnaChromosome`: Ribose `2'-OH` chromosome capped below `~28 kb` by intramolecular
///   transesterification hydrolysis and vulnerable to host/viral RNases.
/// - `ViralModifiedUDna`: Anaerobic Class III `4Fe-4S` thiyl-radical `urRNR` (`2'-OH -> 2'-H`,
///   Lundin et al., 2015, PMC4390871) yields chemically stable `Uracil-DNA (U-DNA)` first selected in
///   viruses to evade host RNA-nuclease defenses (Forterre, 2006, PMC1450140), **but** vulnerable to
///   unrecognizable `C -> U` cytosine deamination transitions!
/// - `DeaminationRepairableTDna`: Thymidylate synthase (`ThyA/ThyX`, `dUMP -> dTMP`) + `dUTPase` +
///   **Uracil-DNA Glycosylase (`UDG`)** base-excision repair (Sousa et al., 2009, PMC2660913;
///   Vertessy & Toth, 2009, PMC2732909) distinguishes deaminated `U` lesions from native `T`!
/// - `CoevolvedCellularVirosphere`: Temperate viral/plasmid transduction of replicative DNA polymerases
///   (`PolC/DnaE` & `PolB/PolD`) + allosteric `dNTP` balance (Torrents, 2014, PMC4009431) + adaptive
///   antiviral immunity (`CRISPR/restriction`, Koonin et al., 2020, PMC7062200) stabilizes Mb chromosomes!
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum GenomeSubstrateStage {
    LabileRnaChromosome = 0,
    ViralModifiedUDna = 1,
    DeaminationRepairableTDna = 2,
    #[default]
    CoevolvedCellularVirosphere = 3,
}

impl GenomeSubstrateStage {
    pub fn short_code(self) -> &'static str {
        match self {
            Self::LabileRnaChromosome => "RNA_2PrimeOH",
            Self::ViralModifiedUDna => "Viral_U_DNA",
            Self::DeaminationRepairableTDna => "Repair_T_DNA",
            Self::CoevolvedCellularVirosphere => "Cellular_TDNA_Virosphere",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::LabileRnaChromosome => {
                "Labile 2'-OH RNA Chromosome (Transesterification Cleavage Cap < 28 kb)"
            }
            Self::ViralModifiedUDna => {
                "RNR Deoxyribose Uracil-DNA (RNase-Resistant, Vulnerable to C->U Deamination)"
            }
            Self::DeaminationRepairableTDna => {
                "Thymine-DNA + UDG Base-Excision Repair (Cytosine Deamination Protected)"
            }
            Self::CoevolvedCellularVirosphere => {
                "Coevolved Megabase T-DNA Chromosome + Regulated dNTPs & Antiviral Immunity"
            }
        }
    }

    pub fn has_deoxyribose_rnr(self) -> bool {
        (self as u8) >= (Self::ViralModifiedUDna as u8)
    }

    pub fn has_thymine_udg_repair(self) -> bool {
        (self as u8) >= (Self::DeaminationRepairableTDna as u8)
    }

    pub fn is_coevolved_megabase_dna(self) -> bool {
        matches!(self, Self::CoevolvedCellularVirosphere)
    }
}

/// Prognostic **Ribonucleotide Reductase (`RNR`), `ThyX/UDG` Deamination Repair & Virosphere Arms-Race** state.
#[derive(Clone, Debug)]
pub struct RnrVirosphereReactor {
    /// Anaerobic `4Fe-4S` radical-SAM `urRNR` thiyl-radical `2'-OH -> 2'-H` catalysis `[0.0, 1.0]`.
    pub rnr_thiyl_radical_activity: f64,
    /// Coupled `ThyA/ThyX` (`dUMP -> dTMP`) + `dUTPase` + `UDG` excision efficiency `[0.0, 1.0]`.
    pub thymidylate_udg_repair: f64,
    /// Temperate virocell / plasmid gene-transfer & adaptive antiviral defense coupling `[0.0, 1.0]`.
    pub virosphere_coevolution_coupling: f64,
    /// Lytic viral burst & parasitic replication pressure `[0.0, 1.0]`.
    pub lytic_viral_pressure: f64,
    /// Allosteric `dNTP` pool stoichiometry balance `[0.0, 1.0]` (Torrents, 2014).
    pub dntp_allosteric_balance: f64,
    /// Unrepaired `C -> U` cytosine deamination mutational load `[0.0, 1.0]`.
    pub deamination_lesion_load: f64,
    /// Maximum stable non-melting chromosome capacity in kilobases (`kb`, `12 kb` RNA -> `2400 kb` T-DNA).
    pub stable_genome_capacity_kb: f64,
    pub megabase_tdna_established: bool,
}

impl RnrVirosphereReactor {
    pub fn new_for_controls(controls: &EnvironmentControls) -> Self {
        if !controls.rna_to_dna_virosphere_mode {
            return Self {
                rnr_thiyl_radical_activity: controls.rnr_thiyl_radical_activity.clamp(0.0, 1.0),
                thymidylate_udg_repair: controls.thymidylate_udg_repair_efficiency.clamp(0.0, 1.0),
                virosphere_coevolution_coupling: controls
                    .virosphere_arms_race_coupling
                    .clamp(0.0, 1.0),
                lytic_viral_pressure: controls.lytic_viral_burst_pressure.clamp(0.0, 1.0),
                dntp_allosteric_balance: 0.96,
                deamination_lesion_load: 0.01,
                stable_genome_capacity_kb: 1850.0,
                megabase_tdna_established: true,
            };
        }
        let rnr = controls.rnr_thiyl_radical_activity.clamp(0.0, 1.0);
        let udg = controls.thymidylate_udg_repair_efficiency.clamp(0.0, 1.0);
        let vir_c = controls.virosphere_arms_race_coupling.clamp(0.0, 1.0);
        let lyt = controls.lytic_viral_burst_pressure.clamp(0.0, 1.0);
        Self {
            rnr_thiyl_radical_activity: rnr,
            thymidylate_udg_repair: udg,
            virosphere_coevolution_coupling: vir_c,
            lytic_viral_pressure: lyt,
            dntp_allosteric_balance: (0.35 + 0.62 * rnr).clamp(0.20, 0.99),
            deamination_lesion_load: (0.48 * (1.0 - udg)).clamp(0.01, 0.85),
            stable_genome_capacity_kb: if rnr < 0.25 { 18.0 } else { 85.0 },
            megabase_tdna_established: false,
        }
    }

    /// Advances `RNR` deoxyribonucleotide synthesis, `C -> U` hydrolytic deamination vs. `ThyX + UDG`
    /// base-excision repair, and viral-cellular polymerase transduction vs. lytic epidemic load!
    pub fn step_rnr_and_virosphere_cycle(&mut self, controls: &EnvironmentControls) {
        self.rnr_thiyl_radical_activity = controls.rnr_thiyl_radical_activity.clamp(0.0, 1.0);
        self.thymidylate_udg_repair = controls.thymidylate_udg_repair_efficiency.clamp(0.0, 1.0);
        self.virosphere_coevolution_coupling =
            controls.virosphere_arms_race_coupling.clamp(0.0, 1.0);
        self.lytic_viral_pressure = controls.lytic_viral_burst_pressure.clamp(0.0, 1.0);

        self.dntp_allosteric_balance = (0.32
            + 0.44 * self.rnr_thiyl_radical_activity
            + 0.24 * self.virosphere_coevolution_coupling)
            .clamp(0.18, 0.99);

        // Spontaneous hydrolytic C -> U deamination grows with single-stranded replication forks unless
        // Thymine + dUTPase + UDG base-excision repair excises uracil!
        let raw_deam_influx = if self.rnr_thiyl_radical_activity >= 0.35 {
            0.34
        } else {
            0.18
        };
        let udg_clearance = 0.56 * self.thymidylate_udg_repair;
        self.deamination_lesion_load = (0.62 * self.deamination_lesion_load
            + raw_deam_influx
            - udg_clearance)
            .clamp(0.008, 0.95);

        let target_kb = if self.rnr_thiyl_radical_activity < 0.30 {
            // 2'-OH transesterification cleavage caps RNA chromosomes below ~28 kb!
            14.0 + 42.0 * self.rnr_thiyl_radical_activity
        } else if self.thymidylate_udg_repair < 0.35 {
            // U-DNA resists 2'-OH hydrolysis (~110 kb), but C->U deamination blocks megabase scaling!
            65.0 + 130.0 * self.thymidylate_udg_repair
        } else {
            // T-DNA + UDG + viral PolB/PolC/PolD proofreading unlocks megabase cellular chromosomes!
            240.0
                + 1950.0
                    * self.rnr_thiyl_radical_activity
                    * self.thymidylate_udg_repair
                    * (0.35 + 0.65 * self.virosphere_coevolution_coupling)
        };
        self.stable_genome_capacity_kb =
            (0.55 * self.stable_genome_capacity_kb + 0.45 * target_kb).clamp(10.0, 2400.0);

        self.megabase_tdna_established = self.rnr_thiyl_radical_activity >= 0.70
            && self.thymidylate_udg_repair >= 0.72
            && self.virosphere_coevolution_coupling >= 0.45
            && self.deamination_lesion_load <= 0.12
            && self.stable_genome_capacity_kb >= 1000.0;
    }
}

/// **8-Epoch Continuous 4.25-Billion-Year ITME Grand Trajectory** stages (Szathmary, 2015, PMC4547294;
/// Kempes et al., 2017, PMC5686401; Wong & Bartlett, 2022, PMC9065981):
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum GrandEpochStage {
    #[default]
    Epoch1HadeanHomochirality = 0,
    Epoch2SgcAndRafProtocells = 1,
    Epoch3RnaToDnaAndLuca = 2,
    Epoch4GoeAndSymbiogenesis = 3,
    Epoch5ClonalMulticellularity = 4,
    Epoch6SuperorganismAgriculture = 5,
    Epoch7SymbolicExosomaticArchives = 6,
    Epoch8Gaia2PlanetaryHomeostasis = 7,
}

impl GrandEpochStage {
    pub const ALL: [GrandEpochStage; 8] = [
        Self::Epoch1HadeanHomochirality,
        Self::Epoch2SgcAndRafProtocells,
        Self::Epoch3RnaToDnaAndLuca,
        Self::Epoch4GoeAndSymbiogenesis,
        Self::Epoch5ClonalMulticellularity,
        Self::Epoch6SuperorganismAgriculture,
        Self::Epoch7SymbolicExosomaticArchives,
        Self::Epoch8Gaia2PlanetaryHomeostasis,
    ];

    pub fn code(self) -> &'static str {
        match self {
            Self::Epoch1HadeanHomochirality => "E1_Homochirality",
            Self::Epoch2SgcAndRafProtocells => "E2_SGC_RAF",
            Self::Epoch3RnaToDnaAndLuca => "E3_TDNA_LUCA",
            Self::Epoch4GoeAndSymbiogenesis => "E4_GOE_LECA",
            Self::Epoch5ClonalMulticellularity => "E5_GermSoma",
            Self::Epoch6SuperorganismAgriculture => "E6_Superorg_Agri",
            Self::Epoch7SymbolicExosomaticArchives => "E7_Archives_FormalAI",
            Self::Epoch8Gaia2PlanetaryHomeostasis => "E8_Gaia2_Homeostasis",
        }
    }

    pub fn geologic_age_ga(self) -> f64 {
        match self {
            Self::Epoch1HadeanHomochirality => 4.25,
            Self::Epoch2SgcAndRafProtocells => 4.18,
            Self::Epoch3RnaToDnaAndLuca => 4.05,
            Self::Epoch4GoeAndSymbiogenesis => 1.85,
            Self::Epoch5ClonalMulticellularity => 0.80,
            Self::Epoch6SuperorganismAgriculture => 0.06,
            Self::Epoch7SymbolicExosomaticArchives => 0.005,
            Self::Epoch8Gaia2PlanetaryHomeostasis => 0.000,
        }
    }

    pub fn active_heritable_channels(self) -> u8 {
        match self {
            Self::Epoch1HadeanHomochirality => 1,
            Self::Epoch2SgcAndRafProtocells => 1,
            Self::Epoch3RnaToDnaAndLuca => 2,
            Self::Epoch4GoeAndSymbiogenesis => 2,
            Self::Epoch5ClonalMulticellularity => 3,
            Self::Epoch6SuperorganismAgriculture => 3,
            Self::Epoch7SymbolicExosomaticArchives => 4,
            Self::Epoch8Gaia2PlanetaryHomeostasis => 5,
        }
    }
}

/// Telemetry checkpoint emitted at the end of each epoch during `run_continuous_grand_trajectory`.
#[derive(Clone, Debug)]
pub struct GrandTrajectoryEpochSummary {
    pub epoch: GrandEpochStage,
    pub tick_end: usize,
    pub geologic_age_ga: f64,
    pub population_size: usize,
    pub heritable_channels_active: u8,
    pub heritable_mutual_info_bits: f64,
    pub informational_leverage_li: f64,
    pub order_velocity_vo: f64,
    pub planetary_o2: f64,
    pub planetary_co2_ppm: f64,
    pub planetary_exergy_budget: f64,
    pub landmark_summary: String,
}

/// Per-epoch comparative planetary habitability & biodiversity turnover summary.
#[derive(Clone, Debug)]
pub struct EpochHabitabilityTurnoverSummary {
    pub epoch_code: String,
    pub geologic_age_ga: f64,
    pub tick_end: usize,
    pub population_size: usize,
    pub heritable_channels_active: u8,
    pub heritable_mutual_info_bits: f64,
    pub informational_leverage_li: f64,
    pub order_velocity_vo: f64,
    pub planetary_o2_pal: f64,
    pub planetary_co2_ppm: f64,
    pub planetary_exergy_budget: f64,
    pub habitability_index: f64,
    pub biodiversity_turnover_index: f64,
    pub delta_o2_vs_prev: f64,
    pub delta_mi_vs_prev: f64,
    pub landmark_summary: String,
}

/// Multi-epoch comparative habitability & biodiversity turnover digest across the 4.25 Ga trajectory.
#[derive(Clone, Debug)]
pub struct MultiEpochComparativeDigest {
    pub seed: u64,
    pub ticks_per_epoch: usize,
    pub epoch_count: usize,
    pub mean_habitability_index: f64,
    pub peak_habitability_epoch: Option<String>,
    pub cumulative_o2_rise_pal: f64,
    pub mean_turnover_index: f64,
    pub total_mi_gain_bits: f64,
    pub epochs: Vec<EpochHabitabilityTurnoverSummary>,
}

/// Emergent ecological & individuality guild determined by sequence motifs, mixotrophic kleptoplasty,
/// primary/secondary/tertiary endosymbiotic organelle fixation, anaerobic MRO reduction, and secondary
/// reductive evolution of photosynthetic holobionts into parasitic apicoplast-bearing lineages.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum TrophicGuild {
    Producer = 0,
    Syntroph = 1,
    PredatorParasite = 2,
    Holobiont = 3,
}

impl TrophicGuild {
    pub fn label(self) -> &'static str {
        match self {
            Self::Producer => "Producer (TL 1.0)",
            Self::Syntroph => "Syntroph Cross-Feeder (TL 1.5)",
            Self::PredatorParasite => "Predator / Klepto / Apico / Mitosome-Parasite (TL 1.65-2.1)",
            Self::Holobiont => "Endosymbiotic / MRO Holobiont (TL 1.25-1.55)",
        }
    }

    pub fn short_code(self) -> &'static str {
        match self {
            Self::Producer => "Prod",
            Self::Syntroph => "Synt",
            Self::PredatorParasite => "Pred",
            Self::Holobiont => "Holo",
        }
    }

    pub fn trophic_level(self) -> f64 {
        match self {
            Self::Producer => 1.0,
            Self::Holobiont => 1.25,
            Self::Syntroph => 1.5,
            Self::PredatorParasite => 2.0,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum OrderPhase {
    Adaptive,
    Stasis,
    Regressive,
}

impl OrderPhase {
    pub fn from_velocity(v_o: f64, eps: f64) -> Self {
        if v_o > eps {
            Self::Adaptive
        } else if v_o < -eps {
            Self::Regressive
        } else {
            Self::Stasis
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Adaptive => "Adaptive (v_O > 0)",
            Self::Stasis => "Dynamic Stasis (v_O ~ 0)",
            Self::Regressive => "Regressive (v_O < 0)",
        }
    }
}

#[derive(Copy, Clone, Debug, Default)]
pub struct PriceDecomposition {
    pub selection_covariance: f64,
    pub transmission_bias: f64,
    pub selection_share: f64,
}

impl PriceDecomposition {
    pub fn from_parent_offspring_records(records: &[(f64, f64, f64)]) -> Self {
        if records.len() < 2 {
            return Self::default();
        }
        let n = records.len() as f64;
        let mean_w = (records.iter().map(|r| r.0).sum::<f64>() / n).max(1e-6);
        let mean_z = records.iter().map(|r| r.1).sum::<f64>() / n;

        let mut cov_wz = 0.0;
        let mut trans_wz = 0.0;
        for &(w, z_p, z_c) in records {
            cov_wz += (w - mean_w) * (z_p - mean_z);
            trans_wz += w * (z_c - z_p);
        }
        cov_wz /= n * mean_w;
        trans_wz /= n * mean_w;

        let denom = cov_wz.abs() + trans_wz.abs() + 1e-4;
        let selection_share = if cov_wz > 0.0 {
            (cov_wz / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };

        Self {
            selection_covariance: cov_wz,
            transmission_bias: trans_wz,
            selection_share,
        }
    }
}

#[derive(Copy, Clone, Debug, Default)]
pub struct ContinuousPhaseCoords {
    pub net_exergy_margin: f64,
    pub order_velocity: f64,
    pub plasticity_fraction: f64,
    pub heritable_capacity: f64,
    pub price_selection_index: f64,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct DiagnosticVector {
    pub v_oc: u8,
    pub m_c: u8,
    pub r_c: u8,
    pub v_c: u8,
    pub p_c: u8,
}

impl DiagnosticVector {
    pub const fn new(v_oc: u8, m_c: u8, r_c: u8, v_c: u8, p_c: u8) -> Self {
        Self {
            v_oc,
            m_c,
            r_c,
            v_c,
            p_c,
        }
    }
}

impl fmt::Display for DiagnosticVector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Dc = {{{},{},{},{},{}}}",
            self.v_oc, self.m_c, self.r_c, self.v_c, self.p_c
        )
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct CapabilityFlags {
    pub exergy_driven: bool,
    pub memory: bool,
    pub replication: bool,
    pub interaction: bool,
    pub variation: bool,
    pub selection: bool,
    pub soft_inheritance: bool,
    pub symbiogenesis: bool,
    pub kleptoplasty: bool,
    pub klepto_tipping_fixation: bool,
    pub endosymbiotic_fusion: bool,
    pub nuclear_egt: bool,
    pub serial_endosymbiosis: bool,
    pub nucleomorph_reduction: bool,
    pub apicoplast_reduction: bool,
    pub apicoplast_loss: bool,
    pub hydrogenosome_fermentation: bool,
    pub mitosome_isc_retention: bool,
    pub amiteukaryote_mro_loss: bool,
    pub clonal_multicellularity: bool,
    pub germ_soma_division: bool,
    pub apoptotic_cheater_policing: bool,
    pub eusocial_colony: bool,
    pub sterile_worker_castes: bool,
    pub worker_kin_policing: bool,
    pub agricultural_domestication: bool,
    pub clonal_crop_monoculture: bool,
    pub defensive_microbiome_symbiosis: bool,
    pub dual_gene_culture_inheritance: bool,
    pub cumulative_cultural_ratchet: bool,
    pub symbolic_pedagogical_teaching: bool,
    pub exosomatic_memory_offloading: bool,
    pub written_archival_storage: bool,
    pub institutional_peer_verification: bool,
    pub computational_symbiont: bool,
    pub formal_proof_verification: bool,
    pub model_collapse_safeguard: bool,
    pub oxygenic_goe_regulation: bool,
    pub syntrophic_nutrient_closure: bool,
    pub planetary_homeostatic_awakening: bool,
    pub prebiotic_raf_closure: bool,
    pub ribozyme_hypercycle_active: bool,
    pub darwinian_threshold_luca: bool,
    pub homochiral_d_rna_l_aa: bool,
    pub stereoselective_aminoacylation: bool,
    pub error_minimized_universal_sgc: bool,
    pub rnr_deoxyribose_synthesis: bool,
    pub thymine_udg_deamination_repair: bool,
    pub virosphere_dna_coevolution: bool,
    pub niche_construction: bool,
    pub cross_feeding: bool,
    pub trophic_predation: bool,
}

#[derive(Clone, Debug)]
pub struct EnvironmentControls {
    pub bath_exergy_flux: f64,
    pub landauer_bit_cost: f64,
    pub exergy_diffusion: f64,
    pub niche_construction_rate: f64,
    pub niche_decay_rate: f64,
    pub speciation_threshold: f64,
    pub trophic_coupling: f64,
    pub trophic_efficiency: f64,
    pub cross_feeding_rate: f64,
    pub kleptoplasty_rate: f64,
    pub klepto_hgt_rate: f64,
    pub endosymbiosis_rate: f64,
    pub egt_rate: f64,
    pub serial_endosymbiosis_rate: f64,
    pub nucleomorph_reduction_rate: f64,
    pub apicoplast_reduction_rate: f64,
    pub host_metabolite_rescue: f64,
    /// Ambient dissolved oxygen / redox acceptor availability `[0.0 = anoxic, 1.0 = normoxic]`;
    /// hypoxia (`< 0.45`) favors shedding `OXPHOS` ETC & `mtDNA` into `Hydrogenosome`s & `Mitosome`s!
    pub ambient_oxygen_level: f64,
    /// Rate of anaerobic MRO reductive evolution (`AerobicMito -> AnaerobicMito -> Hydrogenosome -> Mitosome`).
    pub anaerobic_mro_reduction_rate: f64,
    /// Bacterial cytosolic `SUF` (`SufSU/SufB/SufC`) lateral-gene-transfer replacement `[0.0, 1.0]`;
    /// when `>= 0.65`, bypasses mitochondrial `ISC` essentiality to unlock *Monocercomonoides* `MRO -> 0`!
    pub cytosolic_suf_hgt_rescue: f64,
    /// Rate of incomplete cytokinesis / cadherin clonal staying-together (`Holo -> UndifferentiatedClonal`).
    pub multicellular_adhesion_rate: f64,
    /// Rate of `regA`-like soma-germline differentiation & mitochondrial dirty-work germline sequestration.
    pub germ_soma_differentiation_rate: f64,
    /// Mitochondrial `cytochrome c` / `MOMP` / metacaspase apoptosis rate policing somatic cheaters & cleaving propagules.
    pub apoptosis_policing_rate: f64,
    /// Enforces a single-cell (`N = 1`) zygotic bottleneck at reproduction (`r -> 1.0`, purging somatic cheaters).
    pub single_cell_bottleneck: bool,
    /// Queen lifetime monogamy commitment `[0.0 = promiscuous/polyandrous, 1.0 = strict monandry r_sib = 0.50]`;
    /// `>= 0.75` opens Boomsma's (2009) Monogamy Window to cross into `ObligateSuperorganism`!
    pub queen_monogamy_index: f64,
    /// Juvenile-hormone / blastogenic caste canalization drive (`GermSoma -> Subsocial -> Eusocial -> Superorganism`).
    pub eusocial_caste_canalization_rate: f64,
    /// Mutual worker kin-policing rate `[0.0, 1.0]` destroying selfish rebel worker-laid eggs!
    pub worker_policing_rate: f64,
    /// Ectosymbiotic cultivar domestication & polyploidization drive (`Forager -> LowerAg -> GongylCrop -> LeafCutter`).
    pub crop_domestication_rate: f64,
    /// Specialized *Escovopsis*-like agro-pathogen pressure `[0.0, 1.0]` targeting clonal crop monocultures!
    pub agro_pathogen_pressure: f64,
    /// Cuticular exocrine-crypt *Pseudonocardia* / *Streptomyces* antimicrobial potency `[0.0, 1.0]`.
    pub defensive_symbiont_potency: f64,
    /// Social learning & active symbolic teaching fidelity `[0.0, 1.0]` (`>= 0.82` cumulative ratchet; `>= 0.88` symbolic teaching).
    pub symbolic_teaching_fidelity: f64,
    /// Inter-deme migration & oblique contact rate `[0.0, 1.0]` across structured sub-populations (`0.15..=0.70` Derex-Boyd optimum).
    pub inter_deme_connectivity: f64,
    /// Cultural exploration, modification & cross-deme recombination rate `[0.0, 1.0]`.
    pub cultural_innovation_rate: f64,
    /// Drive to externalize symbolic repertoires into physical/digital written archives `[0.0, 1.0]`.
    pub archival_externalization_rate: f64,
    /// Institutional replication & peer-verification rigor `[0.0, 1.0]` purging false-positive claims!
    pub peer_verification_rigor: f64,
    /// Smaldino-McElreath volume-incentive unverified publication/deposit pressure `[0.0, 1.0]`.
    pub unverified_publish_pressure: f64,
    /// Autonomous computational sampler / program-synthesis generation rate `[0.0, 1.0]`.
    pub synthetic_generation_rate: f64,
    /// FunSearch / AlphaGeometry formal proof-checker & deterministic evaluator strictness `[0.0, 1.0]`.
    pub formal_verifier_strictness: f64,
    /// Shumailov (2024) fresh/preserved empirical data anchor fraction `gamma_real in [0.0, 1.0]` (`>= 0.12` prevents tail collapse!).
    pub empirical_grounding_floor: f64,
    /// Strength of biotic & technospheric planetary O2/CO2/N-P feedback coupling `[0.0, 1.0]`.
    pub biogeochemical_feedback_gain: f64,
    /// Gaia-2.0 machine-verified planetary stewardship & renewable exergy coupling `[0.0, 1.0]`.
    pub planetary_stewardship_coupling: f64,
    /// Unguided superlinear fossil/carbon combustion intensity `[0.0, 1.0]` (Wong & Bartlett, 2022).
    pub extractive_burn_intensity: f64,
    /// Enables de-novo Hadean prebiotic monomer-to-LUCA abiogenesis dynamics when `true`.
    pub hadean_abiogenesis_mode: bool,
    /// Hot-spring hydration-dehydration-gel condensation amplitude `[0.0, 1.0]` (Damer & Deamer, 2020).
    pub wet_dry_cycle_amplitude: f64,
    /// Submarine/geothermal `Fe(Ni)S` natural proton gradient `Delta-pH` `[0.0, 4.5]` (Lane & Martin, 2012).
    pub hydrothermal_delta_ph: f64,
    /// Citrate / acidic-peptide `Mg2+` chelation `[0.0, 1.0]` protecting fatty-acid vesicles (Adamala & Szostak, 2013).
    pub citrate_mg_chelation: f64,
    /// Ribozyme polymerase + RNP error-correction fidelity `[0.0, 1.0]` governing Eigen's Error Cliff!
    pub ribozyme_polymerase_fidelity: f64,
    /// Enables explicit racemic-to-homochiral `D-RNA / L-AA` & `SGC` crystallization dynamics when `true`.
    pub chiral_symmetry_breaking_mode: bool,
    /// Ozturk-Sasselov spin-polarized magnetite (`Fe3O4`) `CISS` + Frank chiral bias `[0.0, 1.0]`.
    pub magnetite_ciss_spin_bias: f64,
    /// Tamura-Schimmel `3'-NCCA` RNA minihelix steric `L-aminoacyl` selectivity `[0.0, 1.0]`.
    pub minihelix_chiral_selectivity: f64,
    /// Vetsigian-Woese-Goldenfeld communal progenote `HGT` code-sharing rate `[0.0, 1.0]`.
    pub communal_code_hgt_rate: f64,
    /// Enables explicit `RNR` `RNA -> U-DNA -> T-DNA` + `UDG` repair + virosphere dynamics when `true`.
    pub rna_to_dna_virosphere_mode: bool,
    /// Anaerobic Class III `4Fe-4S` radical-SAM `urRNR` `2'-OH -> 2'-H` activity `[0.0, 1.0]`.
    pub rnr_thiyl_radical_activity: f64,
    /// Thymidylate synthase (`ThyA/ThyX`) + `dUTPase` + `UDG` `C->U` repair efficiency `[0.0, 1.0]`.
    pub thymidylate_udg_repair_efficiency: f64,
    /// Viral-to-cellular polymerase/RNR symbiogenesis & antiviral defense coupling `[0.0, 1.0]`.
    pub virosphere_arms_race_coupling: f64,
    /// Lytic viral burst & predation pressure `[0.0, 1.0]`.
    pub lytic_viral_burst_pressure: f64,
    pub stability_epsilon: f64,
    pub memory_plasticity: f64,
    pub replication_rate: f64,
    pub replication_fidelity: f64,
    pub heritable_variation: f64,
    pub interaction_strength: f64,
    pub carrying_capacity: usize,
    pub soft_inheritance_coupling: f64,
    pub horizontal_pedagogy: f64,
    pub teacher_count: usize,
    pub recombination_rate: f64,
    pub selective_superiority: f64,
    pub env_drift_rate: f64,
    pub artificial_selection: bool,
}

impl Default for EnvironmentControls {
    fn default() -> Self {
        Self {
            bath_exergy_flux: 0.0,
            landauer_bit_cost: 0.0022,
            exergy_diffusion: 0.10,
            niche_construction_rate: 0.14,
            niche_decay_rate: 0.04,
            speciation_threshold: 0.28,
            trophic_coupling: 0.0,
            trophic_efficiency: 0.36,
            cross_feeding_rate: 0.0,
            kleptoplasty_rate: 0.0,
            klepto_hgt_rate: 0.0,
            endosymbiosis_rate: 0.0,
            egt_rate: 0.0,
            serial_endosymbiosis_rate: 0.0,
            nucleomorph_reduction_rate: 0.0,
            apicoplast_reduction_rate: 0.0,
            host_metabolite_rescue: 0.0,
            ambient_oxygen_level: 0.85,
            anaerobic_mro_reduction_rate: 0.0,
            cytosolic_suf_hgt_rescue: 0.0,
            multicellular_adhesion_rate: 0.0,
            germ_soma_differentiation_rate: 0.0,
            apoptosis_policing_rate: 0.0,
            single_cell_bottleneck: true,
            queen_monogamy_index: 0.0,
            eusocial_caste_canalization_rate: 0.0,
            worker_policing_rate: 0.0,
            crop_domestication_rate: 0.0,
            agro_pathogen_pressure: 0.25,
            defensive_symbiont_potency: 0.0,
            symbolic_teaching_fidelity: 0.0,
            inter_deme_connectivity: 0.35,
            cultural_innovation_rate: 0.0,
            archival_externalization_rate: 0.0,
            peer_verification_rigor: 0.0,
            unverified_publish_pressure: 0.30,
            synthetic_generation_rate: 0.0,
            formal_verifier_strictness: 0.0,
            empirical_grounding_floor: 0.20,
            biogeochemical_feedback_gain: 0.0,
            planetary_stewardship_coupling: 0.0,
            extractive_burn_intensity: 0.15,
            hadean_abiogenesis_mode: false,
            wet_dry_cycle_amplitude: 0.82,
            hydrothermal_delta_ph: 3.20,
            citrate_mg_chelation: 0.84,
            ribozyme_polymerase_fidelity: 0.88,
            chiral_symmetry_breaking_mode: false,
            magnetite_ciss_spin_bias: 0.88,
            minihelix_chiral_selectivity: 0.86,
            communal_code_hgt_rate: 0.85,
            rna_to_dna_virosphere_mode: false,
            rnr_thiyl_radical_activity: 0.90,
            thymidylate_udg_repair_efficiency: 0.92,
            virosphere_arms_race_coupling: 0.84,
            lytic_viral_burst_pressure: 0.24,
            stability_epsilon: 0.08,
            memory_plasticity: 0.0,
            replication_rate: 0.0,
            replication_fidelity: 0.0,
            heritable_variation: 0.0,
            interaction_strength: 0.0,
            carrying_capacity: 32,
            soft_inheritance_coupling: 0.0,
            horizontal_pedagogy: 0.0,
            teacher_count: 1,
            recombination_rate: 0.0,
            selective_superiority: 9.0,
            env_drift_rate: 0.0,
            artificial_selection: false,
        }
    }
}

impl EnvironmentControls {
    pub fn effective_fidelity(&self) -> f64 {
        let f0 = self.replication_fidelity.clamp(0.0, 1.0);
        if self.horizontal_pedagogy <= 0.001 || self.teacher_count <= 1 {
            return f0;
        }
        let exponent =
            1.0 + 0.85 * self.horizontal_pedagogy * ((1 + self.teacher_count) as f64).ln();
        (1.0 - (1.0 - f0).powf(exponent)).clamp(0.0, 0.9995)
    }
}

/// Complete ITME v12.0 (`v1.2.0`) telemetry snapshot including Anaerobic MRO Reductive Evolution
/// (`Mito -> Hydrogenosome -> Mitosome -> MRO->0`), Apicoplast Reductive Evolution, Mixotrophic
/// Kleptoplasty, Shopping-Bag `HGT` Tipping Points, `1°/2°/3°` Serial Endosymbiosis & `eNewick` DAGs.
#[derive(Clone, Debug)]
pub struct ItmeMetrics {
    pub step: usize,
    pub alphabet_size: usize,
    pub structure_raw_bits: f64,
    pub structure_s: f64,
    pub adami_conservation_rho: f64,
    pub complexity_c: f64,
    pub leverage_li: f64,
    pub stability_st: f64,
    pub volatility_nu: f64,
    pub order_o: f64,
    pub ecosystem_order_pop: f64,
    pub order_velocity_vo: f64,
    pub net_exergy_phi: f64,
    pub landauer_cost_q: f64,
    pub entropy_production_rate: f64,
    pub mean_grid_exergy: f64,
    pub mean_niche_capital: f64,
    pub mean_byproduct_pool: f64,
    pub species_count: usize,
    pub fst_fixation_index: f64,
    pub producer_count: usize,
    pub syntroph_count: usize,
    pub predator_count: usize,
    pub aerobic_mito_count: usize,
    pub hydrogenosome_count: usize,
    pub mitosome_count: usize,
    pub amiteukaryote_count: usize,
    pub mean_mro_pathways: f64,
    pub hydrogenosome_transition_events: usize,
    pub mitosome_transition_events: usize,
    pub mro_complete_loss_events: usize,
    pub clonal_multicell_count: usize,
    pub germ_soma_count: usize,
    pub mean_cell_count: f64,
    pub mean_soma_fraction: f64,
    pub mean_somatic_cheater_load: f64,
    pub clonal_multicell_transitions: usize,
    pub germ_soma_transitions: usize,
    pub apoptotic_cull_events: usize,
    pub neoplastic_collapse_events: usize,
    pub eusocial_colony_count: usize,
    pub superorganism_count: usize,
    pub mean_colony_size: f64,
    pub mean_sterile_caste_fraction: f64,
    pub mean_rebel_worker_load: f64,
    pub eusocial_transitions: usize,
    pub superorganism_transitions: usize,
    pub worker_policing_events: usize,
    pub colony_anarchy_collapses: usize,
    pub farming_colony_count: usize,
    pub obligate_agro_count: usize,
    pub mean_gongylidia_yield: f64,
    pub mean_escovopsis_load: f64,
    pub agriculture_transitions: usize,
    pub obligate_crop_transitions: usize,
    pub antimicrobial_suppression_events: usize,
    pub garden_blight_collapses: usize,
    pub cultural_lineage_count: usize,
    pub symbolic_culture_count: usize,
    pub mean_cultural_repertoire_bits: f64,
    pub mean_neural_teaching_capacity: f64,
    pub cultural_ratchet_transitions: usize,
    pub symbolic_teaching_transitions: usize,
    pub cross_deme_recombinations: usize,
    pub cultural_drift_collapses: usize,
    pub archival_lineage_count: usize,
    pub verified_institution_count: usize,
    pub mean_exosomatic_archive_bits: f64,
    pub mean_archive_fdr_load: f64,
    pub written_archive_transitions: usize,
    pub peer_institution_transitions: usize,
    pub peer_retraction_audits: usize,
    pub archive_corruption_collapses: usize,
    pub computational_symbiont_count: usize,
    pub formal_symbio_count: usize,
    pub mean_verified_theorem_bits: f64,
    pub mean_tail_diversity: f64,
    pub machine_symbiosis_transitions: usize,
    pub formal_verifier_transitions: usize,
    pub formal_proofs_verified: usize,
    pub recursive_model_collapses: usize,
    pub planetary_o2_level: f64,
    pub planetary_co2_ppm: f64,
    pub planetary_nutrient_pool: f64,
    pub planetary_exergy_budget: f64,
    pub chemical_disequilibrium_bits: f64,
    pub boundary_transgression_load: f64,
    pub goe_transitions: usize,
    pub homeostatic_awakening_transitions: usize,
    pub planetary_burnout_collapses: usize,
    pub mean_rna_oligomer_nt: f64,
    pub raf_catalytic_closure: f64,
    pub eigen_max_length_nt: usize,
    pub raf_protocell_transitions: usize,
    pub darwinian_luca_transitions: usize,
    pub eigen_error_catastrophes: usize,
    pub mg_vesicle_collapses: usize,
    pub ribose_ee_d: f64,
    pub amino_ee_l: f64,
    pub sgc_error_minimization_mp: f64,
    pub encoded_amino_acids: u8,
    pub homochiral_transitions: usize,
    pub sgc_crystallization_transitions: usize,
    pub racemic_cross_inhibitions: usize,
    pub mistranslation_babel_collapses: usize,
    pub genome_capacity_kb: f64,
    pub deamination_lesion_load: f64,
    pub dntp_pool_balance: f64,
    pub rnr_dna_transitions: usize,
    pub udg_thymine_transitions: usize,
    pub cytosine_deamination_meltdowns: usize,
    pub lytic_viral_collapses: usize,
    pub apicoplast_parasite_count: usize,
    pub mean_apicoplast_pathways: f64,
    pub apicoplast_transition_events: usize,
    pub apicoplast_loss_events: usize,
    pub kleptoplastic_host_count: usize,
    pub mean_kleptoplast_integrity: f64,
    pub mean_klepto_hgt_loci: f64,
    pub klepto_hgt_events: usize,
    pub klepto_fixation_events: usize,
    pub mosaic_donor_diversity: f64,
    pub holobiont_count: usize,
    pub enslaved_holobiont_count: usize,
    pub secondary_holobiont_count: usize,
    pub tertiary_holobiont_count: usize,
    pub active_nucleomorph_count: usize,
    pub mean_nucleomorph_len: f64,
    pub nucleomorph_loss_events: usize,
    pub serial_endosymbiotic_fusions: usize,
    pub mean_egt_fraction: f64,
    pub mean_organelle_len: f64,
    pub egt_transfer_events: usize,
    pub trophic_transfer_flux: f64,
    pub cross_feeding_flux: f64,
    pub mean_trophic_level: f64,
    pub food_web_entropy_bits: f64,
    pub tmrca_steps: usize,
    pub sackin_imbalance: f64,
    pub speciation_events: usize,
    pub guild_transition_events: usize,
    pub endosymbiotic_fusions: usize,
    pub active_reticulations: usize,
    pub fidelity_f: f64,
    pub observed_fidelity: f64,
    pub heritable_mutual_info_bits: f64,
    pub critical_fidelity_theta: f64,
    pub price: PriceDecomposition,
    pub coords: ContinuousPhaseCoords,
    pub population_size: usize,
    pub mean_fitness: f64,
    pub flags: CapabilityFlags,
    pub diagnostic: DiagnosticVector,
    pub phase: OrderPhase,
}

pub fn compute_food_web_entropy(
    producers: usize,
    syntrophs: usize,
    predators: usize,
    holobionts: usize,
) -> f64 {
    let total = producers + syntrophs + predators + holobionts;
    if total == 0 {
        return 0.0;
    }
    let nf = total as f64;
    let mut h = 0.0;
    for c in [producers, syntrophs, predators, holobionts] {
        if c > 0 {
            let p = (c as f64) / nf;
            h -= p * p.log2();
        }
    }
    h
}

pub fn lz76_complexity(seq: &[u16]) -> usize {
    let n = seq.len();
    if n == 0 {
        return 0;
    }
    let mut c = 1usize;
    let mut l = 1usize;
    let mut i = 0usize;
    let mut k = 1usize;
    let mut k_max = 1usize;

    while l + k <= n {
        if seq[i + k - 1] == seq[l + k - 1] {
            k += 1;
            if l + k > n {
                c += 1;
                break;
            }
        } else {
            if k > k_max {
                k_max = k;
            }
            i += 1;
            if i == l {
                c += 1;
                l += k_max;
                if l >= n {
                    break;
                }
                i = 0;
                k = 1;
                k_max = 1;
            } else {
                k = 1;
            }
        }
    }
    c
}

pub fn compute_raw_structure_bits(
    seq: &[u16],
    alphabet_size: usize,
    essential_edges: usize,
) -> f64 {
    if seq.is_empty() {
        return 2.0;
    }
    let sigma_bits = (alphabet_size.max(2) as f64).log2();
    let phrases = lz76_complexity(seq) as f64;
    let len_log = (seq.len().max(2) as f64).log2();
    let seq_bits = phrases * (sigma_bits + len_log * 0.35);
    let edge_bits = (essential_edges as f64) * (len_log * 0.25);
    (seq_bits + edge_bits).max(2.0)
}

pub fn compute_adami_conservation(sequences: &[&[u16]], alphabet_size: usize) -> f64 {
    if sequences.len() < 2 || alphabet_size < 2 {
        return 1.0;
    }
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    if max_len == 0 {
        return 1.0;
    }
    let h_max = (alphabet_size as f64).log2();
    if !h_max.is_finite() || h_max <= 1e-9 {
        return 1.0;
    }
    let mut entropy_sum = 0.0;
    let mut evaluated_loci = 0usize;
    let mut counts = vec![0usize; alphabet_size];

    for pos in 0..max_len {
        counts.fill(0);
        let mut present = 0usize;
        for seq in sequences {
            if let Some(&sym) = seq.get(pos) {
                counts[(sym as usize) % alphabet_size] += 1;
                present += 1;
            }
        }
        if present >= 2 {
            evaluated_loci += 1;
            let mut h_pos = 0.0;
            let n_f = present as f64;
            for &c in &counts {
                if c > 0 {
                    let p = (c as f64) / n_f;
                    h_pos -= p * p.log2();
                }
            }
            // Adami et al. (2000) physical complexity normalizes site entropy by log2(A).
            let norm_h = (h_pos / h_max).clamp(0.0, 1.0);
            entropy_sum += norm_h;
        }
    }

    if evaluated_loci == 0 {
        return 1.0;
    }
    (1.0 - entropy_sum / (evaluated_loci as f64)).clamp(0.08, 1.0)
}

pub fn compute_effective_structure_bits(
    s_raw: f64,
    rho_phys: f64,
    rho_ens: f64,
    replicating_pop_size: usize,
) -> f64 {
    let s_min = 2.0;
    let lambda_p = if replicating_pop_size <= 1 {
        0.0
    } else {
        (replicating_pop_size as f64) / ((replicating_pop_size as f64) + 8.0)
    };
    let blended_rho =
        ((1.0 - lambda_p) * rho_phys.clamp(0.05, 1.0) + lambda_p * rho_ens.clamp(0.08, 1.0))
            .clamp(0.06, 1.0);
    s_min + (s_raw - s_min).max(0.0) * blended_rho
}

pub fn compute_critical_fidelity(structure_raw_s: f64, selective_superiority: f64) -> f64 {
    let ln_sigma = selective_superiority.max(1.25).ln();
    let effective_digits = (structure_raw_s * 0.22).max(4.5);
    (1.0 - ln_sigma / effective_digits).clamp(0.72, 0.999)
}

pub fn compute_heritable_mutual_info(f_obs: f64, alphabet_size: usize, s_raw: f64) -> f64 {
    let q = alphabet_size.max(2) as f64;
    let f = f_obs.clamp(1.0 / q, 0.999_999);
    let err = 1.0 - f;
    let h_cond = if err <= 1e-9 {
        0.0
    } else {
        -(f * f.log2() + err * (err / (q - 1.0)).log2())
    };
    let cap_frac = (1.0 - h_cond / q.log2()).clamp(0.0, 1.0);
    s_raw * cap_frac
}

pub fn compute_thermodynamic_hazard(
    bath_flux: f64,
    net_exergy_margin: f64,
    mean_niche: f64,
) -> f64 {
    if bath_flux.abs() <= 0.03 && net_exergy_margin.abs() <= 0.03 {
        0.0
    } else if net_exergy_margin < -0.02 {
        (0.42 + 0.55 * (-net_exergy_margin)).clamp(0.35, 0.92)
    } else {
        let niche_shield = 1.0 / (1.0 + 0.65 * mean_niche.max(0.0));
        ((0.042 / (0.55 + net_exergy_margin)) * niche_shield).clamp(0.008, 0.065)
    }
}

#[derive(Clone, Debug)]
pub struct TelemetryWindow {
    window_dt: usize,
    li_history: VecDeque<f64>,
    order_history: VecDeque<f64>,
}

impl TelemetryWindow {
    pub fn new(window_dt: usize) -> Self {
        Self {
            window_dt: window_dt.max(4),
            li_history: VecDeque::with_capacity(window_dt + 2),
            order_history: VecDeque::with_capacity(window_dt + 2),
        }
    }

    pub fn reset(&mut self) {
        self.li_history.clear();
        self.order_history.clear();
    }

    pub fn push_and_evaluate(
        &mut self,
        leverage_li: f64,
        epsilon: f64,
        thermo_hazard: f64,
        meltdown_hazard: f64,
    ) -> (f64, f64, f64, f64) {
        self.li_history.push_back(leverage_li);
        if self.li_history.len() > self.window_dt {
            self.li_history.pop_front();
        }

        let eps = epsilon.max(1e-4);
        let mut downside_losses = 0usize;
        let mut total_jumps = 0usize;
        let mut pairs = 0usize;

        if self.li_history.len() >= 2 {
            let slice: Vec<f64> = self.li_history.iter().copied().collect();
            for w in slice.windows(2) {
                pairs += 1;
                let delta = w[1] - w[0];
                if delta < -eps {
                    downside_losses += 1;
                }
                if delta.abs() > eps {
                    total_jumps += 1;
                }
            }
        }

        let p_downside = if pairs == 0 {
            0.0
        } else {
            ((0.25 + downside_losses as f64) / (4.25 + pairs as f64)).clamp(0.0, 0.96)
        };

        let volatility_nu = if pairs == 0 {
            0.0
        } else {
            ((total_jumps as f64) / (pairs as f64)).clamp(0.0, 1.0)
        };

        let combined_hazard = (thermo_hazard + meltdown_hazard).clamp(0.0, 0.94);
        let stability_st = ((1.0 - combined_hazard) * (1.0 - p_downside)).clamp(0.02, 1.0);
        let order_o = stability_st * leverage_li;

        self.order_history.push_back(order_o);
        if self.order_history.len() > self.window_dt {
            self.order_history.pop_front();
        }

        let v_o = if self.order_history.len() >= 2 {
            let n = self.order_history.len() as f64;
            let mean_x = (n - 1.0) * 0.5;
            let mean_y = self.order_history.iter().sum::<f64>() / n;
            let mut num = 0.0;
            let mut den = 0.0;
            for (idx, &y) in self.order_history.iter().enumerate() {
                let dx = (idx as f64) - mean_x;
                num += dx * (y - mean_y);
                den += dx * dx;
            }
            if den > 1e-9 {
                num / den
            } else {
                0.0
            }
        } else {
            0.0
        };

        (stability_st, volatility_nu, order_o, v_o)
    }
}
