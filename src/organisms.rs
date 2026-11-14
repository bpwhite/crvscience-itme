//! ITME v12.0 (`v1.2.0`) Anaerobic Reductive Evolution of Mitochondria into Hydrogenosomes,
//! Mitosomes (`ISC` Fe-S Core) & Cytosolic-`SUF` Amiteukaryotes (`MRO -> 0`), alongside Reductive
//! Apicoplasts, Mixotrophic Kleptoplasty & Serial `1°/2°/3°` Endosymbiosis (`src/organisms.rs`):
//! - Models the **5-Class Mitochondrion-Related Organelle (`MRO`) Continuum** across oxygen gradients
//!   (`step_anaerobic_mro_reduction`; Stairs, Leger, & Roger, 2015; Shiflett & Johnson, 2010; Hjort et al., 2010):
//!   `AerobicMitochondrion` (`OXPHOS` + `3`-locus `CoRR` `mtDNA`) -> `AnaerobicMitochondrion` (*Blastocystis*)
//!   -> genome-less `H2`-producing `Hydrogenosome` (*Trichomonas*: `PFO/HydA` + substrate `ATP` + `ISC`, `mtDNA=0`
//!   confirming Allen 2015 `CoRR`) -> genome-less non-ATP `Mitosome` (*Giardia* / *Entamoeba* / microsporidia:
//!   retaining only `MRO_PATH_ISC` `[2Fe-2S]` assembly for cytosolic `CIA` & nuclear replicative enzymes;
//!   Tsaousis, 2019; Braymer & Lill, 2017; Motyčková et al., 2023).
//! - Models **Cytosolic `SUF` HGT Rescue & Complete Mitochondrial Loss (`MRO -> 0`)** (*Monocercomonoides exilis*;
//!   Karnkowska et al., 2016; Hampl et al., 2016; Treitli et al., 2019): disrupting `ISC` without replacement
//!   triggers **CIA `[Fe-S]` Replicative Collapse**, whereas `cytosolic_suf_hgt_rescue >= 0.65` unlocks viable
//!   complete mitosome elimination (`MroClass::Amiteukaryote`)!

use crate::model::{
    compute_effective_structure_bits, compute_raw_structure_bits, AbiogenesisStage,
    AgricultureState, ChiralCodeReactor, ChiralityCodeStage, CulturalState, EnvironmentControls,
    EusocialState, ExosomaticState, GenomeSubstrateStage, HadeanGeochemicalBath,
    MachineSynthesisState, MroClass, MulticellularState, PlanetaryBiosphere, PlanetaryState,
    RnrVirosphereReactor, Rng, TrophicGuild,
};
use crate::regimes::RegimeId;
use crate::spatial::{BiomeZone, SpatialGrid};
use std::fmt;

/// Minimal bioenergetic organelle core retained in situ for rapid redox regulation (Allen, 2015 CoRR).
pub const CORR_CORE_LOCI: usize = 3;

/// Periplastidial eukaryotic housekeeping core retained in a nucleomorph while >=1 plastid gene remains.
pub const NUCLEOMORPH_HOUSEKEEPING_LOCI: usize = 4;

/// Nuclear HGT targeting/repair/division locus threshold ("targeting-early" tipping point, Hehenberger et al., 2019).
pub const KLEPTO_FIXATION_HGT_THRESHOLD: usize = 4;

/// Minimal non-photosynthetic apicoplast housekeeping genome (`sufB` + `clpC` + rRNAs, ~35 kb circular genome).
pub const APICOPLAST_MIN_LOCI: usize = 1;

/// Essential MEP / DOXP Isoprenoid precursor pathway (`IPP` / `DMAPP`; Yeh & DeRisi, 2011).
pub const APICO_PATH_MEP_IPP: u8 = 0b0001;
/// Essential `SUF` Iron-Sulfur cluster assembly pathway (`sufB` plastid-encoded cofactor for MEP; Okada & Sigala, 2023).
pub const APICO_PATH_SUF_FES: u8 = 0b0010;
/// Type II Fatty Acid Synthesis (`FASII`, dispensable in blood stages when scavenging host lipids; McFadden & Yeh, 2017).
pub const APICO_PATH_FASII: u8 = 0b0100;
/// Co-localized Plastid-Mitochondrial Heme Biosynthesis (`Heme`, dispensable when scavenging host heme).
pub const APICO_PATH_HEME: u8 = 0b1000;
/// Indispensable blood-stage apicoplast core (`MEP/IPP + SUF/Fe-S`).
pub const APICO_PATH_ESSENTIAL_CORE: u8 = APICO_PATH_MEP_IPP | APICO_PATH_SUF_FES;
/// Complete 4-pathway anabolic plastid complement.
pub const APICO_PATH_ALL: u8 = 0b1111;

/// Indispensable mitochondrial Iron-Sulfur Cluster (`ISC`: `IscS/Nfs1`, `IscU/Isu1`) `[2Fe-2S]` assembly
/// feeding cytosolic `CIA` maturation of nuclear DNA polymerases & helicases (Tsaousis, 2019; Braymer & Lill, 2017).
pub const MRO_PATH_ISC: u8 = 0b0001;
/// Organelle substrate-level phosphorylation (`ASCT` / succinyl-CoA synthetase `ATP` generation).
pub const MRO_PATH_SUBSTRATE_ATP: u8 = 0b0010;
/// Anaerobic pyruvate:ferredoxin oxidoreductase (`PFO`) + `[FeFe]`-hydrogenase (`HydA`) evolving `H2` + acetate!
pub const MRO_PATH_PFO_HYD: u8 = 0b0100;
/// Aerobic proton-pumping Electron Transport Chain (`Complexes I/III/IV` + `O2` oxidase) requiring `CoRR` `mtDNA`!
pub const MRO_PATH_OXPHOS_ETC: u8 = 0b1000;
/// Canonical 3-pathway genome-less `Hydrogenosome` complement (`ISC | SUBSTRATE_ATP | PFO_HYD`).
pub const MRO_PATH_HYDROGENOSOME: u8 = MRO_PATH_ISC | MRO_PATH_SUBSTRATE_ATP | MRO_PATH_PFO_HYD;
/// Minimalist 1-pathway genome-less `Mitosome` core (`ISC` Fe-S cluster assembly only; *Giardia* / microsporidia).
pub const MRO_PATH_MITOSOME_CORE: u8 = MRO_PATH_ISC;
/// Complete 4-pathway `MRO` bioenergetic + biosynthetic complement.
pub const MRO_PATH_ALL: u8 = 0b1111;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum OrganismKind {
    RnaDnaProtein,
    HominidPopulation,
    ArtificialIntelligence,
}

impl OrganismKind {
    pub const ALL: [OrganismKind; 3] = [
        OrganismKind::RnaDnaProtein,
        OrganismKind::HominidPopulation,
        OrganismKind::ArtificialIntelligence,
    ];

    pub fn slug(self) -> &'static str {
        match self {
            Self::RnaDnaProtein => "rna",
            Self::HominidPopulation => "hominid",
            Self::ArtificialIntelligence => "ai",
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::RnaDnaProtein => "Small RNA/DNA/Protein Complexes",
            Self::HominidPopulation => "Hominid Populations",
            Self::ArtificialIntelligence => "Artificial Intelligences",
        }
    }

    pub fn short_title(self) -> &'static str {
        match self {
            Self::RnaDnaProtein => "RNA/DNA/Protein",
            Self::HominidPopulation => "Hominids",
            Self::ArtificialIntelligence => "AI Systems",
        }
    }

    pub fn from_slug(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "1" | "rna" | "dna" | "protein" | "rnp" | "molecular" => Some(Self::RnaDnaProtein),
            "2" | "hominid" | "hominids" | "human" | "tribe" | "cultural" => {
                Some(Self::HominidPopulation)
            }
            "3" | "ai" | "artificial" | "neural" | "agent" => {
                Some(Self::ArtificialIntelligence)
            }
            _ => None,
        }
    }

    pub fn alphabet_symbols(self) -> &'static [&'static str] {
        match self {
            Self::RnaDnaProtein => &[
                "A", "U/T", "G", "C", "Ala", "Cys", "Asp", "Glu", "His", "Lys", "Arg", "Trp",
            ],
            Self::HominidPopulation => &[
                "Forage",
                "Strike",
                "Rotate",
                "Thin",
                "Symmetrize",
                "Haft",
                "Hearth",
                "Share",
                "Imitate",
                "Norm",
            ],
            Self::ArtificialIntelligence => &[
                "START", "Linear", "ReLU", "Conv", "Attn", "MemRead", "MemWrite", "Residual",
                "Verify", "END",
            ],
        }
    }

    pub fn alphabet_size(self) -> usize {
        self.alphabet_symbols().len()
    }

    pub fn alphabet_description(self) -> &'static str {
        match self {
            Self::RnaDnaProtein => {
                "Sigma = {A, U/T, G, C} + catalytic/structural residues {Ala, Cys, Asp, Glu, His, Lys, Arg, Trp}"
            }
            Self::HominidPopulation => {
                "Sigma = {Forage, Strike, Rotate, Thin, Symmetrize, Haft, Hearth, Share, Imitate, Norm}"
            }
            Self::ArtificialIntelligence => {
                "Sigma = {START, Linear, ReLU, Conv, Attn, MemRead, MemWrite, Residual, Verify, END}"
            }
        }
    }

    pub fn niche_artifact_label(self) -> &'static str {
        match self {
            Self::RnaDnaProtein => "Catalytic Lipid-Mineral Biofilm",
            Self::HominidPopulation => "Hearths, Lithic Caches & Shelters",
            Self::ArtificialIntelligence => "Shared External Vector Memory & Tool Index",
        }
    }

    pub fn biome_target_symbols(self, biome: BiomeZone, env_phase: f64) -> [u16; 5] {
        let shift = ((env_phase.sin().abs() * 2.0) as u16) % 2;
        match (self, biome) {
            (Self::RnaDnaProtein, BiomeZone::Alpha) => [2, 3, 8, 6, 5 + shift],
            (Self::RnaDnaProtein, BiomeZone::Beta) => [0, 1, 7, 9, 10 + shift],

            (Self::HominidPopulation, BiomeZone::Alpha) => [1, 2, 3, 4, 6 + shift],
            (Self::HominidPopulation, BiomeZone::Beta) => [2, 3, 5, 7, 8 + shift],

            (Self::ArtificialIntelligence, BiomeZone::Alpha) => [3, 4, 6, 5, 8],
            (Self::ArtificialIntelligence, BiomeZone::Beta) => [1, 2, 7, 4 + shift, 8],
        }
    }

    pub fn regime_narrative(self, regime: RegimeId) -> &'static str {
        match (self, regime) {
            (Self::RnaDnaProtein, RegimeId::R00Static) => {
                "Frozen nucleotide-peptide crystal on an equilibrium lattice with zero exergy diffusion."
            }
            (Self::RnaDnaProtein, RegimeId::R05Decaying) => {
                "Landauer maintenance deficit hydrolyzing oligomers as biofilm niches weather away."
            }
            (Self::RnaDnaProtein, RegimeId::R10AdaptiveDissipation) => {
                "Geothermal reaction-diffusion rings sustained over hydrothermal exergy hotspots."
            }
            (Self::RnaDnaProtein, RegimeId::R15Autoadaptive) => {
                "Chemotactic allosteric ribozyme climbing spatial exergy gradients within its lifetime."
            }
            (Self::RnaDnaProtein, RegimeId::R15Replicative) => {
                "Clonal duplex wave expanding across an uncrowded reaction-diffusion substrate grid."
            }
            (Self::RnaDnaProtein, RegimeId::R20InteractiveReplicative) => {
                "Local resource depletion halos create spatial competition among identical RNA clones."
            }
            (Self::RnaDnaProtein, RegimeId::R20AutoadaptiveReplicative) => {
                "Chemotactic ribozymes construct local chelator niches while daughter strands reset."
            }
            (Self::RnaDnaProtein, RegimeId::R20ProtoDarwinian) => {
                "Mutational meltdown prevents stable specialization across Alpha/Beta vent biomes."
            }
            (Self::RnaDnaProtein, RegimeId::R30FullyDarwinian) => {
                "Takeuchi-Hogeweg replicase-parasite waves + Black-Queen syntrophic cross-feeders!"
            }
            (Self::RnaDnaProtein, RegimeId::R30SoftInheritance) => {
                "Biofilm niche inheritance + epigenetic defense stabilize multi-trophic RNP consortia."
            }
            (Self::RnaDnaProtein, RegimeId::R35Symbiogenetic) => {
                "Aerobic mitochondria reducing under hypoxia to H2 hydrogenosomes & ISC-only mitosomes!"
            }
            (Self::RnaDnaProtein, RegimeId::R30PseudoDarwinian) => {
                "Spatially patterned microfluidic bioreactors balancing synthetic producer-syntroph webs."
            }
            (Self::RnaDnaProtein, RegimeId::R30PseudoSoftInheritance) => {
                "Cytosolic SUF & IPP rescue unlocking Monocercomonoides/Cryptosporidium organelle loss!"
            }

            (Self::HominidPopulation, RegimeId::R00Static) => {
                "Inert river cobble landscape at zero caloric flux."
            }
            (Self::HominidPopulation, RegimeId::R05Decaying) => {
                "Drought & abandonment eroding hearths, lithic caches, and toolmaking chains."
            }
            (Self::HominidPopulation, RegimeId::R10AdaptiveDissipation) => {
                "Seasonal foraging loop tracking reaction-diffusion riparian resource oases."
            }
            (Self::HominidPopulation, RegimeId::R15Autoadaptive) => {
                "Individual landscape learning & hearth tending along resource gradients."
            }
            (Self::HominidPopulation, RegimeId::R15Replicative) => {
                "Band fission wave colonizing open territory from an ancestral valley oasis."
            }
            (Self::HominidPopulation, RegimeId::R20InteractiveReplicative) => {
                "Local game/forage depletion around band camps creates territorial competition."
            }
            (Self::HominidPopulation, RegimeId::R20AutoadaptiveReplicative) => {
                "Bands build local hearth/shelter niches while unpedagogical offspring reset skills."
            }
            (Self::HominidPopulation, RegimeId::R20ProtoDarwinian) => {
                "Irreversible reductive flaking errors prevent regional lithic tradition divergence."
            }
            (Self::HominidPopulation, RegimeId::R30FullyDarwinian) => {
                "Multi-guild specialization: primary quarry-foragers, exchange specialists & raiders!"
            }
            (Self::HominidPopulation, RegimeId::R30SoftInheritance) => {
                "Cooperative norm defense + cross-feeding exchange stabilize regional cultural guilds."
            }
            (Self::HominidPopulation, RegimeId::R35Symbiogenetic) => {
                "Smelting hubs streamlining in scarce-fuel zones into minimalist iron-cofactor guilds!"
            }
            (Self::HominidPopulation, RegimeId::R30PseudoDarwinian) => {
                "Planned regional settlement networks balancing trophic production and specialization."
            }
            (Self::HominidPopulation, RegimeId::R30PseudoSoftInheritance) => {
                "Lateral craft adoption allowing full decentralization of legacy smelting guilds."
            }

            (Self::ArtificialIntelligence, RegimeId::R00Static) => {
                "Unpowered compute grid with zero power flux and zero cached memory state."
            }
            (Self::ArtificialIntelligence, RegimeId::R05Decaying) => {
                "Grid brownout causing DRAM decoherence and decay of shared vector indices."
            }
            (Self::ArtificialIntelligence, RegimeId::R10AdaptiveDissipation) => {
                "Oscillatory compute load balancing across datacenter power hotspots."
            }
            (Self::ArtificialIntelligence, RegimeId::R15Autoadaptive) => {
                "Single agent migrating toward high-compute nodes while tuning internal weights."
            }
            (Self::ArtificialIntelligence, RegimeId::R15Replicative) => {
                "Containerized worm propagating across adjacent cluster nodes."
            }
            (Self::ArtificialIntelligence, RegimeId::R20InteractiveReplicative) => {
                "Local GPU/power throttling around dense clone clusters."
            }
            (Self::ArtificialIntelligence, RegimeId::R20AutoadaptiveReplicative) => {
                "Agents write local external cache indices (niche) while spawns boot cold weights."
            }
            (Self::ArtificialIntelligence, RegimeId::R20ProtoDarwinian) => {
                "High architecture mutation shatters operator pipelines across compute zones."
            }
            (Self::ArtificialIntelligence, RegimeId::R30FullyDarwinian) => {
                "Multi-trophic AI ecology: foundation generators, distilled edge syntrophs & adversaries!"
            }
            (Self::ArtificialIntelligence, RegimeId::R30SoftInheritance) => {
                "Cryptographic Verify defense + shared embedding streams sustain multi-tier AI guilds."
            }
            (Self::ArtificialIntelligence, RegimeId::R35Symbiogenetic) => {
                "Edge AI stacks shedding GPU training blocks into minimalist hardware root-of-trust enclaves!"
            }
            (Self::ArtificialIntelligence, RegimeId::R30PseudoDarwinian) => {
                "Adversarial red-team / blue-team NAS arenas co-evolving robust AI architectures."
            }
            (Self::ArtificialIntelligence, RegimeId::R30PseudoSoftInheritance) => {
                "Pure-software cryptographic SUF verification replacing hardware enclave co-processors."
            }
        }
    }
}

impl fmt::Display for OrganismKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.title())
    }
}

#[derive(Clone, Debug)]
pub struct ModeledUnit {
    pub id: u64,
    pub kind: OrganismKind,
    pub generation: u32,
    pub age: u32,
    pub grid_x: usize,
    pub grid_y: usize,
    pub species_id: usize,
    pub guild: TrophicGuild,
    pub endo_tier: u8,
    pub organelle_membranes: u8,
    pub sequence: Vec<u16>,
    pub kleptoplast_sequence: Vec<u16>,
    pub kleptoplast_integrity: f64,
    pub kleptoplast_donor_id: Option<u64>,
    pub kleptoplast_donor_tier: u8,
    pub klepto_hgt_loci: usize,
    pub klepto_donor_mask: u32,
    pub has_apicoplast_flag: bool,
    pub apicoplast_pathways_mask: u8,
    pub apicoplast_eliminated: bool,
    /// Mitochondrion-Related Organelle (`MRO`) evolutionary state (`AerobicMito -> Hydrogenosome -> Mitosome -> Amiteukaryote`).
    pub mro_class: MroClass,
    /// Bitmask of retained MRO pathways (`MRO_PATH_ISC`, `MRO_PATH_SUBSTRATE_ATP`, `MRO_PATH_PFO_HYD`, `MRO_PATH_OXPHOS_ETC`).
    pub mro_pathways_mask: u8,
    /// Retained mitochondrial `mtDNA` `CoRR` loci (`>=3` in aerobic mitochondria, `0` in hydrogenosomes & mitosomes!).
    pub mro_mtdna_loci: usize,
    /// Individuality grade (`Unicellular`, `AggregativeChimerism`, `UndifferentiatedClonal`, `GermSomaMulticellular`).
    pub multicell_state: MulticellularState,
    /// Number of constituent cells (`1` in unicellular protists/holobionts, `8..=64` in clonal multicellular bodies).
    pub cell_count: u16,
    /// Fraction of cells differentiated into sterile foraging/motile/defensive **soma** (`[0.0, 0.88]`).
    pub soma_fraction: f64,
    /// Within-organism genetic relatedness `r in [0.0, 1.0]` (`1.0` with clonal `N=1` zygotic bottleneck).
    pub kin_relatedness_r: f64,
    /// Selfish proliferating **somatic cheater / neoplasm** fraction (`[0.0, 1.0]`); policed by mitochondrial apoptosis!
    pub somatic_cheater_load: f64,
    /// Cumulative mitochondrial `cytochrome c` / metacaspase apoptotic cheater-purge & propagule-scission events.
    pub cumulative_apoptotic_culls: u32,
    /// Eusociality / superorganismality grade (`Solitary`, `SubsocialFamily`, `FacultativeEusocial`, `ObligateSuperorganism`).
    pub eusocial_state: EusocialState,
    /// Number of multicellular organisms integrated into the colony (`1` solitary, `6..=64` eusocial colony).
    pub colony_size: u16,
    /// Morphologically / developmentally sterile **foraging & brood-care worker caste** fraction (`[0.0, 0.76]`).
    pub sterile_worker_fraction: f64,
    /// Specialized sterile **defensive soldier caste** fraction (`[0.0, 0.22]`).
    pub soldier_fraction: f64,
    /// Full/half-sibling brood relatedness `r_sib in [0.25, 0.50]` (`0.50` under strict lifetime monogamy).
    pub colony_relatedness_r: f64,
    /// Selfish **rebel worker reproduction / unpoliced male-egg** defection load (`[0.0, 1.0]`).
    pub rebel_worker_load: f64,
    /// Cumulative rebel worker-laid eggs destroyed by **mutual worker kin-policing**!
    pub cumulative_policed_eggs: u32,
    /// Ectosymbiotic agriculture grade (`HunterForager`, `LowerFacultativeFarming`, `HigherDomesticatedCrop`, `LeafCutterAgroEcosystem`).
    pub agriculture_state: AgricultureState,
    /// Domesticated fungal/crop ploidy level (`1` none, `2` diploid lower-attine, `4` tetraploid higher-attine, `8` octoploid leaf-cutter; Kooij et al., 2015).
    pub cultivar_polyploidy: u8,
    /// Genetic clonal uniformity `[0.0, 1.0]` driven by vertical maternal infrabuccal-pellet propagation.
    pub crop_clonal_uniformity: f64,
    /// Autophagic **gongylidia** nutritional reward yield `[0.0, 1.0]` converting cellulose into host energy (Leal-Dutra et al., 2023).
    pub gongylidia_yield: f64,
    /// Specialized *Escovopsis*-like crop mycoparasite burden `[0.0, 1.0]` exploiting clonal monocultures.
    pub escovopsis_pathogen_load: f64,
    /// Cuticular exocrine-crypt *Pseudonocardia* / *Streptomyces* biofilm coverage `[0.0, 1.0]` (Li et al., 2018; Andersen et al., 2015).
    pub defensive_actinobacteria_cover: f64,
    /// Cumulative *Escovopsis* infection waves neutralized by antimicrobial actinobacteria + worker weeding!
    pub cumulative_antimicrobial_shields: u32,
    /// Cumulative cultural evolution grade (`AsocialTrialError`, `LowFidelityTradition`, `ImitativeRatchet`, `SymbolicCumulativeCulture`).
    pub cultural_state: CulturalState,
    /// Structured population deme index (`0..4`) partitioning spatial/social learning sub-networks (Derex & Boyd, 2016).
    pub deme_id: u8,
    /// Second heritable informational channel: cumulative cultural repertoire / symbolic corpus loci!
    pub cultural_sequence: Vec<u16>,
    /// Cumulative technological / symbolic ratchet tier (`0..=6`).
    pub cultural_ratchet_tier: u8,
    /// Genetically co-evolved neural plasticity & active pedagogical teaching capacity `[0.0, 1.0]` (Muthukrishna et al., 2018).
    pub neural_learning_capacity: f64,
    /// Transmission-error & demographic-isolation skill loss `[0.0, 1.0]` (Tasmanian ratchet erosion).
    pub cultural_drift_erosion: f64,
    /// Cumulative cross-deme oblique/horizontal cultural recombination events!
    pub cumulative_cross_deme_recombinations: u32,
    /// Externalized exosomatic archive grade (`EndosomaticBound`, `MaterialNotation`, `WrittenArchive`, `PeerVerifiedInstitution`).
    pub exosomatic_state: ExosomaticState,
    /// Durable written/digital exosomatic archive loci decoupled from individual cranial volume!
    pub archive_sequence: Vec<u16>,
    /// Fraction of internal neural working-memory load offloaded onto external archives `[0.0, 0.45]` (Holocene encephalization relief).
    pub cognitive_offloading_fraction: f64,
    /// False-discovery / unverified noise fraction `[0.0, 1.0]` in the externalized archive (Smaldino & McElreath, 2016).
    pub archive_false_discovery_load: f64,
    /// Cumulative replication audits & false-positive retractions executed by epistemic institutions!
    pub cumulative_peer_verifications: u32,
    /// Neuro-symbolic / machine-verified synthesis grade (`ManualArchival`, `AlgorithmicMediation`, `UnverifiedRecursiveGen`, `MachineVerifiedSymbiosis`).
    pub machine_synthesis_state: MachineSynthesisState,
    /// Fourth heritable informational channel: machine-synthesized executable programs & formally checked proofs!
    pub computational_symbiont_sequence: Vec<u16>,
    /// Empirical tail coverage `[0.0, 1.0]` (`1.0` full support; `< 0.25` Shumailov early tail loss; `< 0.10` delta collapse!).
    pub distribution_tail_diversity: f64,
    /// Recursive synthetic self-training variance degradation `[0.0, 1.0]` (Shumailov et al., 2024).
    pub recursive_autophagy_load: f64,
    /// Cumulative candidate programs/theorems formally checked by deterministic evaluators!
    pub cumulative_formal_proofs: u32,
    /// Current planetary biogeochemical epoch experienced and co-regulated by this guild lineage!
    pub planetary_epoch: PlanetaryState,
    /// Net oxygenic production minus respiration flux contributed by this unit.
    pub guild_o2_flux: f64,
    /// Net atmospheric carbon drawdown + silicate weathering minus extractive combustion.
    pub guild_carbon_drawdown: f64,
    /// Syntrophic & superorganismal N/P detritus remineralization flux.
    pub guild_nutrient_recycling: f64,
    /// Cumulative Gaia-2.0 machine-verified decarbonization & solar-exergy stewardship cycles!
    pub cumulative_stewardship_actions: u32,
    /// Prebiotic-to-cellular abiogenesis stage (`MonomerSoup -> RAFProtocell -> RibozymeHypercycle -> ProgenotePool -> LUCA_Cellular`).
    pub abiogenesis_stage: AbiogenesisStage,
    /// Longest functional RNA/RNP/DNA informational polymer length (`nt`).
    pub rna_oligomer_length_nt: u16,
    /// Local Reflexively Autocatalytic Food-Generated (`RAF`) metabolic closure `[0.0, 1.0]`.
    pub raf_closure_score: f64,
    /// Amphiphilic / isoprenoid / phospholipid protocell membrane stability `[0.0, 1.0]`.
    pub vesicle_membrane_stability: f64,
    /// Horizontal communal innovation transfer fraction (`Woese Progenote` `~0.72` -> `LUCA` `<= 0.14`).
    pub progenote_hgt_flux: f64,
    /// Chiral symmetry-breaking, minihelix stereoselection & genetic code crystallization stage!
    pub chirality_code_stage: ChiralityCodeStage,
    /// Enantiomeric excess of `D-ribose` nucleotides `[0.0, 1.0]`.
    pub ribose_ee: f64,
    /// Enantiomeric excess of `L-aminoacyl` peptides `[0.0, 1.0]`.
    pub peptide_l_ee: f64,
    /// Polar-requirement Error-Minimization Percentage `MP` of this lineage's genetic code `[0.0, 0.88]`.
    pub genetic_code_mp: f64,
    /// Number of canonically encoded `L-amino acids` (`0 -> 4 -> 10 -> 20`).
    pub amino_acid_alphabet_size: u8,
    /// Nucleic-acid chromosome substrate (`RNA_2PrimeOH -> Viral_U_DNA -> Repair_T_DNA -> Cellular_TDNA_Virosphere`).
    pub genome_substrate_stage: GenomeSubstrateStage,
    /// Stable chromosome coding capacity in kilobases (`kb`, `12 kb` RNA -> `1850+ kb` T-DNA).
    pub chromosome_capacity_kb: f64,
    /// Accumulated unrepaired `C -> U` cytosine deamination burden `[0.0, 1.0]`.
    pub deamination_uracil_load: f64,
    /// Adaptive cellular antiviral defense (`CRISPR/restriction-modification` vs lytic virions) `[0.0, 1.0]`.
    pub antiviral_defense_titer: f64,
    pub nucleomorph_sequence: Vec<u16>,
    pub nucleomorph_transferred_loci: usize,
    pub nucleomorph_eliminated: bool,
    pub endosymbiont_sequence: Vec<u16>,
    pub endosymbiont_origin_id: Option<u64>,
    pub egt_transferred_loci: usize,
    pub wiring: Vec<(usize, usize)>,
    pub lifetime_memory: Vec<f64>,
    pub integrity: f64,
    pub structure_raw_s: f64,
    pub structure_s: f64,
    pub base_complexity_c: f64,
    pub complexity_c: f64,
    pub primary_work: f64,
    pub syntrophic_work: f64,
    pub net_trophic_flux: f64,
    pub harvested_work: f64,
    pub landauer_cost: f64,
    pub byproduct_output: f64,
    pub niche_output: f64,
    pub offense_trait: f64,
    pub defense_trait: f64,
    pub fitness: f64,
}

impl ModeledUnit {
    #[allow(clippy::too_many_arguments)]
    pub fn new_founder_with_role(
        id: u64,
        kind: OrganismKind,
        regime: RegimeId,
        role_hint: TrophicGuild,
        grid_x: usize,
        grid_y: usize,
        biome: BiomeZone,
        env_phase: f64,
        rng: &mut Rng,
    ) -> Self {
        let sigma = kind.alphabet_symbols().len() as u16;
        let base_len: usize = match regime {
            RegimeId::R00Static => 16,
            RegimeId::R05Decaying => 20,
            RegimeId::R10AdaptiveDissipation => 18,
            RegimeId::R15Autoadaptive => 20,
            RegimeId::R15Replicative | RegimeId::R20InteractiveReplicative => 18,
            RegimeId::R20AutoadaptiveReplicative => 20,
            RegimeId::R20ProtoDarwinian => 24,
            RegimeId::R30FullyDarwinian => 24,
            RegimeId::R30SoftInheritance => 26,
            RegimeId::R35Symbiogenetic => 28,
            RegimeId::R30PseudoDarwinian => 26,
            RegimeId::R30PseudoSoftInheritance => 28,
        };
        let target_len = match role_hint {
            TrophicGuild::Producer | TrophicGuild::Holobiont => base_len,
            TrophicGuild::Syntroph => (base_len.saturating_sub(5)).max(12),
            TrophicGuild::PredatorParasite => (base_len.saturating_sub(4)).max(13),
        };

        let biome_targets = kind.biome_target_symbols(biome, env_phase);
        let syntroph_motifs: &[u16] = match kind {
            OrganismKind::RnaDnaProtein => &[1, 4, 7],
            OrganismKind::HominidPopulation => &[0, 7, 8],
            OrganismKind::ArtificialIntelligence => &[1, 5, 7],
        };
        let predator_motifs: &[u16] = match kind {
            OrganismKind::RnaDnaProtein => &[2, 3, 9, 10],
            OrganismKind::HominidPopulation => &[1, 3, 5],
            OrganismKind::ArtificialIntelligence => &[3, 4],
        };

        let mut sequence = Vec::with_capacity(target_len);
        match (kind, regime) {
            (_, RegimeId::R00Static) => {
                for i in 0..target_len {
                    if kind == OrganismKind::ArtificialIntelligence && i == 0 {
                        sequence.push(0);
                    } else if kind == OrganismKind::ArtificialIntelligence && i + 1 == target_len {
                        sequence.push(sigma - 1);
                    } else {
                        sequence.push(((i % 2) + 1) as u16 % sigma);
                    }
                }
            }
            (OrganismKind::ArtificialIntelligence, _) => {
                sequence.push(0);
                for i in 1..(target_len - 1) {
                    let sym = match role_hint {
                        TrophicGuild::Syntroph if i % 2 == 0 => {
                            syntroph_motifs[i % syntroph_motifs.len()]
                        }
                        TrophicGuild::PredatorParasite if i % 2 == 0 => {
                            predator_motifs[i % predator_motifs.len()]
                        }
                        _ if regime >= RegimeId::R30FullyDarwinian && i % 2 == 0 => {
                            biome_targets[(i / 2) % biome_targets.len()]
                        }
                        _ => 1 + (rng.gen_index((sigma - 2) as usize) as u16),
                    };
                    sequence.push(sym);
                }
                sequence.push(sigma - 1);
            }
            _ => {
                for i in 0..target_len {
                    let sym = match role_hint {
                        TrophicGuild::Syntroph if i % 2 == 0 => {
                            syntroph_motifs[i % syntroph_motifs.len()]
                        }
                        TrophicGuild::PredatorParasite if i % 2 == 0 => {
                            predator_motifs[i % predator_motifs.len()]
                        }
                        _ if regime >= RegimeId::R30FullyDarwinian && i % 2 == 0 => {
                            biome_targets[(i / 2) % biome_targets.len()]
                        }
                        _ => rng.gen_index(sigma as usize) as u16,
                    };
                    sequence.push(sym);
                }
            }
        }

        let is_holo = role_hint == TrophicGuild::Holobiont;
        let endosymbiont_sequence = if is_holo {
            syntroph_motifs.iter().cycle().take(11).copied().collect()
        } else {
            Vec::new()
        };

        let wiring = Self::build_wiring(kind, &sequence, regime);
        let lifetime_memory = vec![0.0; sequence.len()];
        let ctrl = regime.preset_controls();

        let mut unit = Self {
            id,
            kind,
            generation: 0,
            age: 0,
            grid_x,
            grid_y,
            species_id: 0,
            guild: TrophicGuild::Producer,
            endo_tier: u8::from(is_holo),
            organelle_membranes: if is_holo { 2 } else { 0 },
            sequence,
            kleptoplast_sequence: Vec::new(),
            kleptoplast_integrity: 0.0,
            kleptoplast_donor_id: None,
            kleptoplast_donor_tier: 0,
            klepto_hgt_loci: 0,
            klepto_donor_mask: 0,
            has_apicoplast_flag: false,
            apicoplast_pathways_mask: if is_holo { APICO_PATH_ALL } else { 0 },
            apicoplast_eliminated: false,
            mro_class: if is_holo {
                MroClass::AerobicMitochondrion
            } else {
                MroClass::None
            },
            mro_pathways_mask: if is_holo { MRO_PATH_ALL } else { 0 },
            mro_mtdna_loci: if is_holo { CORR_CORE_LOCI } else { 0 },
            multicell_state: MulticellularState::Unicellular,
            cell_count: 1,
            soma_fraction: 0.0,
            kin_relatedness_r: 1.0,
            somatic_cheater_load: 0.0,
            cumulative_apoptotic_culls: 0,
            eusocial_state: EusocialState::Solitary,
            colony_size: 1,
            sterile_worker_fraction: 0.0,
            soldier_fraction: 0.0,
            colony_relatedness_r: 0.50,
            rebel_worker_load: 0.0,
            cumulative_policed_eggs: 0,
            agriculture_state: AgricultureState::HunterForager,
            cultivar_polyploidy: 1,
            crop_clonal_uniformity: 0.0,
            gongylidia_yield: 0.0,
            escovopsis_pathogen_load: 0.0,
            defensive_actinobacteria_cover: 0.0,
            cumulative_antimicrobial_shields: 0,
            cultural_state: CulturalState::AsocialTrialError,
            deme_id: (((grid_x / 8) * 2 + (grid_y / 8)) % 4) as u8,
            cultural_sequence: Vec::new(),
            cultural_ratchet_tier: 0,
            neural_learning_capacity: 0.0,
            cultural_drift_erosion: 0.0,
            cumulative_cross_deme_recombinations: 0,
            exosomatic_state: ExosomaticState::EndosomaticBound,
            archive_sequence: Vec::new(),
            cognitive_offloading_fraction: 0.0,
            archive_false_discovery_load: 0.0,
            cumulative_peer_verifications: 0,
            machine_synthesis_state: MachineSynthesisState::ManualArchival,
            computational_symbiont_sequence: Vec::new(),
            distribution_tail_diversity: 1.0,
            recursive_autophagy_load: 0.0,
            cumulative_formal_proofs: 0,
            planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
            guild_o2_flux: 0.0,
            guild_carbon_drawdown: 0.0,
            guild_nutrient_recycling: 0.0,
            cumulative_stewardship_actions: 0,
            abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
            rna_oligomer_length_nt: 1280,
            raf_closure_score: 0.96,
            vesicle_membrane_stability: 0.95,
            progenote_hgt_flux: 0.08,
            chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
            ribose_ee: 0.99,
            peptide_l_ee: 0.99,
            genetic_code_mp: 0.79,
            amino_acid_alphabet_size: 20,
            genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
            chromosome_capacity_kb: 1850.0,
            deamination_uracil_load: 0.01,
            antiviral_defense_titer: 0.88,
            nucleomorph_sequence: Vec::new(),
            nucleomorph_transferred_loci: 0,
            nucleomorph_eliminated: false,
            endosymbiont_sequence,
            endosymbiont_origin_id: if is_holo { Some(id) } else { None },
            egt_transferred_loci: 0,
            wiring,
            lifetime_memory,
            integrity: 1.0,
            structure_raw_s: 8.0,
            structure_s: 8.0,
            base_complexity_c: 8.0,
            complexity_c: 8.0,
            primary_work: 0.0,
            syntrophic_work: 0.0,
            net_trophic_flux: 0.0,
            harvested_work: 0.0,
            landauer_cost: 0.0,
            byproduct_output: 0.0,
            niche_output: 0.0,
            offense_trait: 0.0,
            defense_trait: 0.0,
            fitness: 0.5,
        };
        unit.recompute_phenome(
            &ctrl,
            ctrl.bath_exergy_flux,
            0.0,
            0.0,
            0.0,
            biome,
            1.0,
            1,
            env_phase,
        );
        unit
    }

    pub fn new_founder(
        id: u64,
        kind: OrganismKind,
        regime: RegimeId,
        grid_x: usize,
        grid_y: usize,
        biome: BiomeZone,
        env_phase: f64,
        rng: &mut Rng,
    ) -> Self {
        Self::new_founder_with_role(
            id,
            kind,
            regime,
            TrophicGuild::Producer,
            grid_x,
            grid_y,
            biome,
            env_phase,
            rng,
        )
    }

    pub fn has_active_kleptoplast(&self) -> bool {
        !self.kleptoplast_sequence.is_empty() && self.kleptoplast_integrity > 0.05
    }

    pub fn mosaic_donor_count(&self) -> usize {
        self.klepto_donor_mask.count_ones() as usize
    }

    pub fn has_apicoplast(&self) -> bool {
        self.has_apicoplast_flag && !self.apicoplast_eliminated
    }

    pub fn apicoplast_pathway_count(&self) -> u8 {
        self.apicoplast_pathways_mask.count_ones() as u8
    }

    pub fn has_intact_apicoplast_essential_core(&self) -> bool {
        (self.apicoplast_pathways_mask & APICO_PATH_ESSENTIAL_CORE) == APICO_PATH_ESSENTIAL_CORE
            && self.endosymbiont_sequence.len() >= APICOPLAST_MIN_LOCI
    }

    /// Number of retained MRO pathways (`0..=4`: `ISC`, `SUBSTRATE_ATP`, `PFO_HYD`, `OXPHOS_ETC`).
    pub fn mro_pathway_count(&self) -> u8 {
        self.mro_pathways_mask.count_ones() as u8
    }

    /// True when mitochondrial `ISC` (`[2Fe-2S]`) is intact OR replaced by bacterial cytosolic `SUF` HGT!
    /// Required for cytosolic `CIA` assembly of `[4Fe-4S]` nuclear DNA polymerases & repair helicases!
    pub fn has_intact_fes_cluster_biogenesis(&self, cytosolic_suf_hgt_rescue: f64) -> bool {
        match self.mro_class {
            MroClass::None => true,
            MroClass::Amiteukaryote => cytosolic_suf_hgt_rescue >= 0.65,
            _ => {
                (self.mro_pathways_mask & MRO_PATH_ISC) != 0 || cytosolic_suf_hgt_rescue >= 0.65
            }
        }
    }

    /// Computes local dissolved oxygen availability for this organism's spatial patch (`Alpha` oxic vs `Beta` microoxic seep).
    pub fn local_oxygen_at(ambient_oxygen_level: f64, biome: BiomeZone) -> f64 {
        let base = ambient_oxygen_level.clamp(0.0, 1.0);
        match biome {
            BiomeZone::Alpha => base,
            BiomeZone::Beta => (base * 0.72).clamp(0.0, 1.0),
        }
    }

    /// **Stepwise Anaerobic Reductive Evolution of Mitochondria into `Hydrogenosome`s, `Mitosome`s & `Amiteukaryote`s**
    /// (Stairs, Leger, & Roger, 2015; Allen, 2015; Tsaousis, 2019; Karnkowska et al., 2016):
    /// 1. `AerobicMitochondrion -> AnaerobicMitochondrion`: sheds `O2` terminal oxidase (`MRO_PATH_OXPHOS_ETC = 0`),
    ///    trims `mro_mtdna_loci` from `3 -> 2` (*Blastocystis* facultative rhodoquinone intermediate).
    /// 2. `AnaerobicMitochondrion -> Hydrogenosome`: relinquishes membrane electron transport completely ->
    ///    **sweeps `mro_mtdna_loci` to `0`** (confirming Allen's `CoRR` prediction!), retaining `3` nuclear-encoded
    ///    pathways (`MRO_PATH_ISC | MRO_PATH_SUBSTRATE_ATP | MRO_PATH_PFO_HYD = 0b0111`) that generate anaerobic `ATP`
    ///    and evolve `H2` + acetate syntrophic cross-feeding byproducts (*Trichomonas* / anaerobic chytrids)!
    /// 3. `Hydrogenosome -> Mitosome`: sheds `PFO/HydA` (`H2` evolution) and organelle `ATP` generation,
    ///    **strictly arresting at the 1-pathway essential core `MRO_PATH_ISC` (`0b0001`, `mtDNA = 0`)** (*Giardia*,
    ///    *Entamoeba*, *Cryptosporidium*, microsporidia)!
    /// 4. At `Mitosome`:
    ///    - If `cytosolic_suf_hgt_rescue >= 0.65` (*Monocercomonoides exilis* lateral acquisition of bacterial
    ///      cytosolic `SufSU/SufB/SufC`), safely eliminates the mitosome (`MroClass::Amiteukaryote`, `MRO -> 0`)!
    ///    - Else if `force_disrupt_isc` is true, strips `MRO_PATH_ISC` without cytosolic `SUF`, triggering lethal
    ///      **CIA `[Fe-S]` Nuclear Replicative Collapse**!
    ///    - Otherwise, purifying selection preserves `MRO_PATH_ISC` intact!
    /// Returns `(transitioned, new_mro_class, triggered_cia_fes_collapse)`.
    pub fn step_anaerobic_mro_reduction(
        &mut self,
        cytosolic_suf_hgt_rescue: f64,
        force_disrupt_isc: bool,
    ) -> (bool, MroClass, bool) {
        match self.mro_class {
            MroClass::None => (false, MroClass::None, false),
            MroClass::AerobicMitochondrion => {
                self.mro_class = MroClass::AnaerobicMitochondrion;
                self.mro_pathways_mask = MRO_PATH_HYDROGENOSOME;
                self.mro_mtdna_loci = 2;
                (true, self.mro_class, false)
            }
            MroClass::AnaerobicMitochondrion => {
                // Relinquishing proton-pumping ETC eliminates selective pressure for mtDNA (Allen, 2015 CoRR)!
                self.mro_class = MroClass::Hydrogenosome;
                self.mro_pathways_mask = MRO_PATH_HYDROGENOSOME;
                self.mro_mtdna_loci = 0;
                self.egt_transferred_loci = self.egt_transferred_loci.max(3);
                (true, self.mro_class, false)
            }
            MroClass::Hydrogenosome => {
                // Shed fermentative PFO/[FeFe]-hydrogenase & organelle substrate-level ATP -> arrest at ISC-only Mitosome!
                self.mro_class = MroClass::Mitosome;
                self.mro_pathways_mask = MRO_PATH_MITOSOME_CORE;
                self.mro_mtdna_loci = 0;
                (true, self.mro_class, false)
            }
            MroClass::Mitosome => {
                if cytosolic_suf_hgt_rescue >= 0.65 {
                    // Bacterial cytosolic SUF lateral transfer replaces mitochondrial ISC -> Monocercomonoides MRO->0!
                    self.mro_class = MroClass::Amiteukaryote;
                    self.mro_pathways_mask = 0;
                    self.mro_mtdna_loci = 0;
                    return (true, self.mro_class, false);
                }
                if force_disrupt_isc {
                    // Ablating mitosomal ISC without cytosolic SUF starves CIA nuclear DNA polymerases/helicases!
                    self.mro_pathways_mask = 0;
                    self.integrity = (self.integrity * 0.20).max(0.02);
                    return (true, self.mro_class, true);
                }
                (false, self.mro_class, false)
            }
            MroClass::Amiteukaryote => {
                if cytosolic_suf_hgt_rescue < 0.65 && force_disrupt_isc {
                    self.integrity = (self.integrity * 0.20).max(0.02);
                    return (false, self.mro_class, true);
                }
                (false, self.mro_class, false)
            }
        }
    }

    /// Bioenergetic & signaling efficiency of **mitochondrial intrinsic apoptosis** (`cytochrome c` / `MOMP` / metacaspases;
    /// Reece et al., 2011; Hernández-Aguilera et al., 2013): cristate respiratory mitochondria provide full apoptotic
    /// cheater-surveillance & scission fidelity (`1.0`), hydrogenosomes retain partial metacaspase execution (`0.55`),
    /// whereas non-respiratory mitosomes/amiteukaryotes lack `cytochrome c` bioenergetic sensing (`0.15`).
    pub fn mitochondrial_apoptosis_efficiency(&self) -> f64 {
        match self.mro_class {
            MroClass::AerobicMitochondrion | MroClass::AnaerobicMitochondrion => 1.0,
            MroClass::Hydrogenosome => 0.55,
            MroClass::Mitosome | MroClass::Amiteukaryote | MroClass::None => 0.15,
        }
    }

    /// Advances a chimeric eukaryotic `Holobiont` along the **Evolutionary Transition in Individuality (`ETI`)**
    /// to **Clonal Multicellularity** (`UndifferentiatedClonal`, *Saccharomyces* snowflake / choanoflagellate rosette)
    /// and **Soma-Germline Division of Labor** (`GermSomaMulticellular`, *Volvox* `regA` / metazoan quiescent germline;
    /// Ratcliff et al., 2012; Radzvilavicius et al., 2016; Arias Del Angel et al., 2020; König & Nedelcu, 2023)!
    /// When `single_cell_bottleneck` is `false`, unicellular cells instead form an `AggregativeChimerism` mound (`r < 1`).
    pub fn step_clonal_multicellularity(
        &mut self,
        single_cell_bottleneck: bool,
        rng: &mut Rng,
    ) -> (bool, MulticellularState) {
        if self.guild != TrophicGuild::Holobiont {
            return (false, self.multicell_state);
        }
        let has_bioenergetic_organelle = matches!(
            self.mro_class,
            MroClass::AerobicMitochondrion
                | MroClass::AnaerobicMitochondrion
                | MroClass::Hydrogenosome
        ) || !self.endosymbiont_sequence.is_empty();
        if !has_bioenergetic_organelle {
            return (false, self.multicell_state);
        }

        match self.multicell_state {
            MulticellularState::Unicellular => {
                if single_cell_bottleneck {
                    self.multicell_state = MulticellularState::UndifferentiatedClonal;
                    self.cell_count = 12;
                    self.soma_fraction = 0.0;
                    self.kin_relatedness_r = 1.0;
                    self.somatic_cheater_load = 0.02;
                } else {
                    self.multicell_state = MulticellularState::AggregativeChimerism;
                    self.cell_count = 10;
                    self.soma_fraction = 0.15;
                    self.kin_relatedness_r = 0.48;
                    self.somatic_cheater_load = 0.22;
                }
                (true, self.multicell_state)
            }
            MulticellularState::UndifferentiatedClonal | MulticellularState::AggregativeChimerism => {
                // Co-opt regA-like stress suppression + sequester low-ROS quiescent germline!
                self.multicell_state = MulticellularState::GermSomaMulticellular;
                self.cell_count = 32;
                self.soma_fraction = 0.75;
                self.kin_relatedness_r = if single_cell_bottleneck { 1.0 } else { 0.54 };
                if single_cell_bottleneck {
                    self.somatic_cheater_load = self.somatic_cheater_load.min(0.04);
                }
                let adh_sym: u16 = match self.kind {
                    OrganismKind::RnaDnaProtein => 5,
                    OrganismKind::HominidPopulation => 9,
                    OrganismKind::ArtificialIntelligence => 8,
                };
                self.insert_locus_into_host_nucleus(adh_sym, rng);
                (true, self.multicell_state)
            }
            MulticellularState::GermSomaMulticellular => {
                if self.cell_count < 48 {
                    self.cell_count = (self.cell_count + 8).min(48);
                }
                (false, self.multicell_state)
            }
        }
    }

    /// Executes 1 tick of **Multicellular Somatic Turnover, Cheater Mutation & Mitochondrial Apoptotic Policing**
    /// (Aktipis et al., 2015; Pentz, Taylor, & Ratcliff, 2016; Boddy et al., 2017):
    /// - Somatic cell divisions and non-clonal chimerism (`1 - r`) generate selfish proliferating cheaters.
    /// - Mitochondrial `cytochrome c` / metacaspase apoptosis (`apoptosis_policing_rate * mito_eff`) culls rogue
    ///   somatic cheaters and cleaves senescent branch junctions for propagule scission!
    /// - If apoptosis is disabled or polyclonal/chimeric transmission lacks a single-cell bottleneck,
    ///   `somatic_cheater_load` crosses `0.70`, triggering **Somatic Neoplastic Collapse**!
    pub fn step_multicellular_homeostasis(
        &mut self,
        apoptosis_policing_rate: f64,
        single_cell_bottleneck: bool,
    ) -> (u32, bool) {
        if !self.multicell_state.is_multicellular() {
            return (0, false);
        }
        if !single_cell_bottleneck {
            self.kin_relatedness_r = (self.kin_relatedness_r * 0.94).clamp(0.32, 0.78);
        }
        let mito_pcd = self.mitochondrial_apoptosis_efficiency();
        let chimera_conflict = (1.0 - self.kin_relatedness_r).max(0.0) * 0.22;
        let bottleneck_leak = if single_cell_bottleneck { 0.0 } else { 0.085 };
        let de_novo_cheat =
            0.034 + 0.0011 * (self.cell_count as f64) + chimera_conflict + bottleneck_leak;

        let purge_strength = apoptosis_policing_rate.clamp(0.0, 1.0)
            * mito_pcd
            * (0.16 + 0.12 * self.soma_fraction);
        let culls = if purge_strength >= 0.035
            && (self.somatic_cheater_load > 0.015 || self.cell_count >= 12)
        {
            1u32 + u32::from(purge_strength >= 0.11 && self.somatic_cheater_load > 0.06)
        } else {
            0
        };
        self.cumulative_apoptotic_culls = self.cumulative_apoptotic_culls.saturating_add(culls);
        self.somatic_cheater_load =
            (self.somatic_cheater_load + de_novo_cheat - purge_strength).clamp(0.0, 1.0);

        if self.somatic_cheater_load >= 0.70 {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (culls, true);
        }
        (culls, false)
    }

    /// Total nonreproductive colony soma fraction (`sterile_worker_fraction + soldier_fraction`).
    pub fn total_sterile_caste_fraction(&self) -> f64 {
        (self.sterile_worker_fraction + self.soldier_fraction).clamp(0.0, 0.95)
    }

    /// Advances a differentiated `GermSomaMulticellular` organism along the **Major Evolutionary Transition
    /// to Eusocial Superorganismality** (`Solitary -> SubsocialFamily -> FacultativeEusocial -> ObligateSuperorganism`;
    /// Boomsma, 2009; West, Fisher, Gardner, & Kiers, 2015; Qiu et al., 2018; Dai et al., 2024):
    /// - Requires `MulticellularState::GermSomaMulticellular` (building hierarchically on soma-germline animals!).
    /// - Crossing the **Point of No Return** from `FacultativeEusocial` (*Polistes* wasp grade) into
    ///   `ObligateSuperorganism` (*Monomorium* / *Atta* / *Apis* grade with morphologically irreversible worker
    ///   & soldier castes) strictly requires **Boomsma's Lifetime Monogamy Window (`queen_monogamy_index >= 0.75`,
    ///   `r_sib >= 0.4375 ~ 0.50`)** so Hamilton's rule `r * B > C` locks in sterile castes!
    pub fn step_eusocial_superorganism(
        &mut self,
        queen_monogamy_index: f64,
        rng: &mut Rng,
    ) -> (bool, EusocialState) {
        if self.multicell_state != MulticellularState::GermSomaMulticellular {
            return (false, self.eusocial_state);
        }
        let r_sib = (0.25 + 0.25 * queen_monogamy_index.clamp(0.0, 1.0)).clamp(0.25, 0.50);
        self.colony_relatedness_r = r_sib;

        match self.eusocial_state {
            EusocialState::Solitary => {
                self.eusocial_state = EusocialState::SubsocialFamily;
                self.colony_size = 6;
                self.sterile_worker_fraction = 0.25;
                self.soldier_fraction = 0.0;
                self.rebel_worker_load = if queen_monogamy_index >= 0.75 { 0.02 } else { 0.10 };
                (true, self.eusocial_state)
            }
            EusocialState::SubsocialFamily => {
                if queen_monogamy_index < 0.50 {
                    // Promiscuity erodes sibling relatedness below Hamilton's threshold for helper altruism!
                    return (false, self.eusocial_state);
                }
                self.eusocial_state = EusocialState::FacultativeEusocial;
                self.colony_size = 18;
                self.sterile_worker_fraction = 0.60;
                self.soldier_fraction = 0.08;
                self.rebel_worker_load = if queen_monogamy_index >= 0.75 { 0.03 } else { 0.14 };
                (true, self.eusocial_state)
            }
            EusocialState::FacultativeEusocial => {
                // Boomsma's (2009) Lifetime Monogamy Window: irreversible morphological sterile castes
                // evolve ONLY within broods of single lifetime-monogamous mothers (`queen_monogamy >= 0.75`)!
                if queen_monogamy_index < 0.75 {
                    return (false, self.eusocial_state);
                }
                self.eusocial_state = EusocialState::ObligateSuperorganism;
                self.colony_size = 64;
                self.sterile_worker_fraction = 0.72;
                self.soldier_fraction = 0.18;
                self.rebel_worker_load = 0.01;
                let caste_grn_sym: u16 = match self.kind {
                    OrganismKind::RnaDnaProtein => 8,
                    OrganismKind::HominidPopulation => 7,
                    OrganismKind::ArtificialIntelligence => 6,
                };
                self.insert_locus_into_host_nucleus(caste_grn_sym, rng);
                (true, self.eusocial_state)
            }
            EusocialState::ObligateSuperorganism => (false, self.eusocial_state),
        }
    }

    /// Executes 1 tick of **Colony Kin-Conflict & Mutual Worker Policing** (Nonacs, 2006; Olejarz et al., 2016;
    /// Wang et al., 2020; Bourke, 2023):
    /// - Unmated / selfish workers attempt to lay haploid male eggs (`rebel_pressure`), especially when
    ///   polyandry lowers sister-sister relatedness (`0.50 - colony_relatedness_r`).
    /// - **Mutual worker kin-policing** (`worker_policing_rate`) destroys rebel worker-laid eggs in favor of
    ///   the queen's brood, enabling crown superorganisms (like *Apis mellifera*) to remain harmonious even
    ///   under secondary polyandry!
    /// - Conversely, if worker policing is ablated (`worker_policing_rate < 0.15`) in a non-monogamous colony,
    ///   `rebel_worker_load` crosses `0.70`, triggering **Colony Reproductive Anarchy Collapse**!
    pub fn step_colony_kin_policing(
        &mut self,
        queen_monogamy_index: f64,
        worker_policing_rate: f64,
    ) -> (u32, bool) {
        if !self.eusocial_state.is_social() {
            return (0, false);
        }
        let r_sib = (0.25 + 0.25 * queen_monogamy_index.clamp(0.0, 1.0)).clamp(0.25, 0.50);
        self.colony_relatedness_r = r_sib;

        let polyandry_conflict = (0.50 - r_sib).max(0.0) * 0.34;
        let pre_ponr_instability = if self.eusocial_state == EusocialState::FacultativeEusocial
            && queen_monogamy_index < 0.55
        {
            0.055
        } else {
            0.0
        };
        let rebel_pressure = 0.028 + polyandry_conflict + pre_ponr_instability;

        let police_cap = worker_policing_rate.clamp(0.0, 1.0)
            * (0.15 + 0.16 * self.sterile_worker_fraction);
        let policed = if police_cap >= 0.035
            && (self.rebel_worker_load > 0.012 || self.colony_size >= 18)
        {
            1u32 + u32::from(police_cap >= 0.12 && self.rebel_worker_load > 0.05)
        } else {
            0
        };
        self.cumulative_policed_eggs = self.cumulative_policed_eggs.saturating_add(policed);
        self.rebel_worker_load =
            (self.rebel_worker_load + rebel_pressure - police_cap).clamp(0.0, 1.0);

        if self.rebel_worker_load >= 0.70 {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (policed, true);
        }
        (policed, false)
    }

    /// Advances a eusocial colony along the **Multi-Kingdom Agricultural Ectosymbiosis Continuum**
    /// (`HunterForager -> LowerFacultativeFarming -> HigherDomesticatedCrop -> LeafCutterAgroEcosystem`;
    /// Schultz & Brady, 2008; Branstetter et al., 2017; Kooij et al., 2015; Leal-Dutra et al., 2023):
    /// - Requires `EusocialState::is_eusocial()` to initiate `LowerFacultativeFarming` (`2n` detritus fungus).
    /// - Crossing into **obligate clonal crop monoculture** (`HigherDomesticatedCrop` `4n` gongylidia and
    ///   `LeafCutterAgroEcosystem` `8n` *Leucoagaricus gongylophorus* cellulose-bioreactor) strictly requires
    ///   `EusocialState::ObligateSuperorganism` division of labor!
    /// - Triggers reciprocal **host metabolic streamlining** (losing redundant host amino-acid biosynthesis
    ///   loci once fungal gongylidia reliably provision essential nutrients; Shik et al., 2024).
    pub fn step_agricultural_domestication(
        &mut self,
        defensive_symbiont_potency: f64,
        rng: &mut Rng,
    ) -> (bool, AgricultureState) {
        if !self.eusocial_state.is_eusocial() {
            return (false, self.agriculture_state);
        }
        let actino = defensive_symbiont_potency.clamp(0.0, 1.0);

        match self.agriculture_state {
            AgricultureState::HunterForager => {
                self.agriculture_state = AgricultureState::LowerFacultativeFarming;
                self.cultivar_polyploidy = 2;
                self.crop_clonal_uniformity = 0.45;
                self.gongylidia_yield = 0.28;
                self.escovopsis_pathogen_load = 0.04;
                self.defensive_actinobacteria_cover = (actino * 0.75).clamp(0.0, 1.0);
                (true, self.agriculture_state)
            }
            AgricultureState::LowerFacultativeFarming => {
                if self.eusocial_state != EusocialState::ObligateSuperorganism {
                    return (false, self.agriculture_state);
                }
                self.agriculture_state = AgricultureState::HigherDomesticatedCrop;
                self.cultivar_polyploidy = 4;
                self.crop_clonal_uniformity = 0.92;
                self.gongylidia_yield = 0.74;
                self.escovopsis_pathogen_load = 0.03;
                self.defensive_actinobacteria_cover = actino;
                // Reciprocal host metabolic streamlining (loss of redundant arginine biosynthesis locus)!
                if self.sequence.len() > 8 {
                    self.sequence.pop();
                }
                (true, self.agriculture_state)
            }
            AgricultureState::HigherDomesticatedCrop => {
                if self.eusocial_state != EusocialState::ObligateSuperorganism {
                    return (false, self.agriculture_state);
                }
                self.agriculture_state = AgricultureState::LeafCutterAgroEcosystem;
                self.cultivar_polyploidy = 8;
                self.crop_clonal_uniformity = 0.98;
                self.gongylidia_yield = 0.96;
                self.escovopsis_pathogen_load = 0.02;
                self.defensive_actinobacteria_cover = actino;
                let agro_sym: u16 = match self.kind {
                    OrganismKind::RnaDnaProtein => 9,
                    OrganismKind::HominidPopulation => 8,
                    OrganismKind::ArtificialIntelligence => 7,
                };
                self.insert_locus_into_host_nucleus(agro_sym, rng);
                (true, self.agriculture_state)
            }
            AgricultureState::LeafCutterAgroEcosystem => (false, self.agriculture_state),
        }
    }

    /// Executes 1 tick of **Clonal Monoculture Red-Queen Vulnerability vs. Cuticular Antimicrobial Defense**
    /// (Li et al., 2018; Andersen et al., 2015; Goldstein & Klassen, 2020; Ortega et al., 2019):
    /// - Vertical clonal propagation (`crop_clonal_uniformity -> 0.98`) prevents cultivar fruiting conflict
    ///   and locks in polyploid `gongylidia` expression, but leaves the genetically uniform garden acutely
    ///   susceptible to specialized *Escovopsis* mycoparasites (`monoculture_amplification`).
    /// - Cuticular exocrine-crypt actinobacteria (*Pseudonocardia* / *Streptomyces* secreting targeted
    ///   antifungal depsipeptides) plus sterile-worker prophylactic weeding suppress `escovopsis_pathogen_load`.
    /// - Ablating the defensive microbiome (`defensive_symbiont_potency < 0.15`) in a clonal monoculture
    ///   causes `escovopsis_pathogen_load >= 0.70` -> catastrophic **Escovopsis Garden Blight Collapse**!
    pub fn step_agro_pathogen_and_defensive_microbiome(
        &mut self,
        agro_pathogen_pressure: f64,
        defensive_symbiont_potency: f64,
    ) -> (u32, bool) {
        if !self.agriculture_state.is_farming() {
            return (0, false);
        }
        let actino = defensive_symbiont_potency.clamp(0.0, 1.0);
        self.defensive_actinobacteria_cover = actino;

        let monoculture_amplification =
            0.35 + 1.15 * self.crop_clonal_uniformity * self.crop_clonal_uniformity;
        let pathogen_influx = 0.018
            + agro_pathogen_pressure.clamp(0.0, 1.0) * (0.09 + 0.16 * monoculture_amplification);

        let worker_weeding = 0.035 * self.sterile_worker_fraction;
        let antimicrobial_suppression =
            actino * (0.18 + 0.16 * self.sterile_worker_fraction) + worker_weeding;

        let shields = if actino >= 0.20
            && (self.escovopsis_pathogen_load > 0.01 || agro_pathogen_pressure > 0.10)
        {
            1u32 + u32::from(
                actino >= 0.65 && self.agriculture_state.is_obligate_monoculture(),
            )
        } else {
            0
        };
        self.cumulative_antimicrobial_shields =
            self.cumulative_antimicrobial_shields.saturating_add(shields);

        self.escovopsis_pathogen_load =
            (self.escovopsis_pathogen_load + pathogen_influx - antimicrobial_suppression)
                .clamp(0.0, 1.0);

        if self.escovopsis_pathogen_load >= 0.70 {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (shields, true);
        }
        (shields, false)
    }

    /// Effective non-redundant informational content (in bits) of the heritable cultural repertoire.
    pub fn cultural_repertoire_bits(&self) -> f64 {
        if self.cultural_sequence.is_empty() {
            return 0.0;
        }
        let raw = compute_raw_structure_bits(
            &self.cultural_sequence,
            self.kind.alphabet_size().max(8),
            0,
        );
        let health = (1.0 - 0.85 * self.cultural_drift_erosion).clamp(0.10, 1.0);
        (raw * (1.0 + 0.16 * (self.cultural_ratchet_tier as f64)) * health).max(0.0)
    }

    /// Advances a lineage along **Cumulative Cultural Evolution & Dual Gene–Culture Inheritance**
    /// (`AsocialTrialError -> LowFidelityTradition -> ImitativeRatchet -> SymbolicCumulativeCulture`;
    /// Lewis & Laland, 2012; Derex & Boyd, 2016; Legare, 2017; Muthukrishna et al., 2018).
    pub fn step_cumulative_culture_transition(
        &mut self,
        symbolic_teaching_fidelity: f64,
        inter_deme_connectivity: f64,
        rng: &mut Rng,
    ) -> (bool, CulturalState) {
        let alpha = self.kind.alphabet_size().max(8) as u16;
        let deme_offset = (self.deme_id as u16) * 3;

        match self.cultural_state {
            CulturalState::AsocialTrialError => {
                self.cultural_state = CulturalState::LowFidelityTradition;
                self.cultural_ratchet_tier = 1;
                self.neural_learning_capacity = 0.28;
                self.cultural_drift_erosion = 0.04;
                self.cultural_sequence = (0..4u16)
                    .map(|k| (100 + deme_offset * 4 + k) % alpha.max(64))
                    .collect();
                (true, self.cultural_state)
            }
            CulturalState::LowFidelityTradition => {
                if symbolic_teaching_fidelity < 0.82 {
                    return (false, self.cultural_state);
                }
                self.cultural_state = CulturalState::ImitativeRatchet;
                self.cultural_ratchet_tier = 2;
                self.neural_learning_capacity = 0.62;
                self.cultural_drift_erosion = 0.02;
                for k in 4..8u16 {
                    self.cultural_sequence
                        .push((100 + deme_offset * 4 + k) % alpha.max(64));
                }
                (true, self.cultural_state)
            }
            CulturalState::ImitativeRatchet => {
                if symbolic_teaching_fidelity < 0.88 || inter_deme_connectivity < 0.10 {
                    return (false, self.cultural_state);
                }
                self.cultural_state = CulturalState::SymbolicCumulativeCulture;
                self.cultural_ratchet_tier = 4;
                self.neural_learning_capacity = 0.90;
                self.cultural_drift_erosion = 0.01;
                for k in 8..14u16 {
                    self.cultural_sequence
                        .push((100 + deme_offset * 4 + k) % alpha.max(64));
                }
                if self.sequence.len() > 8 {
                    self.sequence.pop();
                }
                let sym_pedagogy_gene: u16 = match self.kind {
                    OrganismKind::RnaDnaProtein => 10,
                    OrganismKind::HominidPopulation => 9,
                    OrganismKind::ArtificialIntelligence => 8,
                };
                self.insert_locus_into_host_nucleus(sym_pedagogy_gene, rng);
                (true, self.cultural_state)
            }
            CulturalState::SymbolicCumulativeCulture => (false, self.cultural_state),
        }
    }

    /// Executes **Cross-Deme Oblique/Horizontal Cultural Recombination** between structured sub-populations
    /// (Derex & Boyd, 2016; Derex, Perreault, & Boyd, 2018; Muthukrishna et al., 2021).
    pub fn fuse_cross_deme_cultural_recombination(
        &mut self,
        teacher: &Self,
        inter_deme_connectivity: f64,
        symbolic_teaching_fidelity: f64,
    ) -> usize {
        if !self.cultural_state.is_cumulative_ratchet()
            || !teacher.cultural_state.is_cultural()
            || symbolic_teaching_fidelity < 0.82
            || inter_deme_connectivity < 0.08
        {
            return 0;
        }
        let max_transfer = if inter_deme_connectivity > 0.82 {
            0usize
        } else if (0.15..=0.72).contains(&inter_deme_connectivity) {
            4usize
        } else {
            2usize
        };
        if max_transfer == 0 {
            return 0;
        }

        let mut added = 0usize;
        for &sym in &teacher.cultural_sequence {
            if added >= max_transfer || self.cultural_sequence.len() >= 22 {
                break;
            }
            if !self.cultural_sequence.contains(&sym) {
                self.cultural_sequence.push(sym);
                added += 1;
            }
        }
        if added > 0 {
            self.cultural_ratchet_tier = (self.cultural_ratchet_tier + 1).clamp(3, 6);
            self.cumulative_cross_deme_recombinations =
                self.cumulative_cross_deme_recombinations.saturating_add(1);
        }
        added
    }

    /// Executes 1 tick of **Dual-Inheritance Teaching Maintenance vs. Demographic/Fidelity Ratchet Drift**
    /// (Lewis & Laland, 2012; Muthukrishna et al., 2018; Whitehead et al., 2019).
    pub fn step_dual_inheritance_and_cultural_drift(
        &mut self,
        symbolic_teaching_fidelity: f64,
        inter_deme_connectivity: f64,
        effective_deme_peers: usize,
    ) -> (u32, bool) {
        if !self.cultural_state.is_cultural() {
            return (0, false);
        }
        let fid = symbolic_teaching_fidelity.clamp(0.0, 1.0);
        let conn = inter_deme_connectivity.clamp(0.0, 1.0);
        // Autonomous written archives act as asynchronous cross-generational teachers (Pegado, 2022; PMC10603192),
        // rescuing even tiny isolated demes (`peers = 1, conn = 0.0`) from Tasmanian demographic skill loss!
        let archival_virtual_teachers = if self.exosomatic_state.has_written_archive()
            && self.archive_false_discovery_load < 0.45
        {
            8.5 * (1.0 - self.archive_false_discovery_load)
        } else {
            0.0
        };
        let n_social = (effective_deme_peers as f64) + 16.0 * conn + archival_virtual_teachers;

        let complexity_demand = 0.020 * (self.cultural_ratchet_tier as f64)
            + 0.004 * (self.cultural_sequence.len() as f64);
        let fidelity_slippage = (0.84 - fid).max(0.0) * 0.42;
        let isolation_loss = if self.cultural_ratchet_tier >= 2 && n_social < 4.5 {
            ((4.5 - n_social) / 4.5) * 0.18
        } else {
            0.0
        };
        let drift_influx = complexity_demand + fidelity_slippage + isolation_loss;

        let network_scale = if n_social < 2.2 {
            0.06
        } else {
            (n_social / 6.0).clamp(0.30, 1.28)
        };
        let teaching_repair =
            fid * (0.15 + 0.19 * self.neural_learning_capacity) * network_scale;

        let taught_events = if fid >= 0.80 && n_social >= 2.5 {
            1u32 + u32::from(
                fid >= 0.88 && self.cultural_state == CulturalState::SymbolicCumulativeCulture,
            )
        } else {
            0
        };

        self.cultural_drift_erosion =
            (self.cultural_drift_erosion + drift_influx - teaching_repair).clamp(0.0, 1.0);

        if self.cultural_drift_erosion > 0.35 && self.cultural_sequence.len() > 3 {
            self.cultural_sequence.pop();
        }
        if self.cultural_drift_erosion >= 0.70 && self.cultural_state.is_cumulative_ratchet() {
            self.integrity = (self.integrity * 0.24).max(0.02);
            return (taught_events, true);
        }
        (taught_events, false)
    }

    /// Verified informational content (in bits) of the externalized exosomatic archive,
    /// discounted by unverified false-discovery pollution (`archive_false_discovery_load`).
    pub fn exosomatic_archive_bits(&self) -> f64 {
        if self.archive_sequence.is_empty() {
            return 0.0;
        }
        let raw = compute_raw_structure_bits(
            &self.archive_sequence,
            self.kind.alphabet_size().max(16),
            0,
        );
        let epistemic_purity = (1.0 - 0.90 * self.archive_false_discovery_load).clamp(0.05, 1.0);
        let inst_mult = match self.exosomatic_state {
            ExosomaticState::PeerVerifiedInstitution => 2.15,
            ExosomaticState::WrittenArchive => 1.65,
            ExosomaticState::MaterialNotation => 1.20,
            ExosomaticState::EndosomaticBound => 0.0,
        };
        (raw * inst_mult * epistemic_purity).max(0.0)
    }

    /// Advances a cumulative cultural lineage along **Externalized Exosomatic Memory & Institutional
    /// Peer-Verification** (`EndosomaticBound -> MaterialNotation -> WrittenArchive -> PeerVerifiedInstitution`;
    /// Pegado, 2022; Grinschgl et al., 2021; González-Forero, 2025; Smaldino & McElreath, 2016):
    /// - Decouples total cultural storage from individual brain size: `archive_sequence` scales up to `38..=48`
    ///   loci while individual `neural_learning_capacity` streamlines (`0.90 -> 0.68`, Holocene encephalization
    ///   relief via cognitive offloading!), cutting per-individual metabolic Landauer brain maintenance cost!
    pub fn step_exosomatic_archive_transition(
        &mut self,
        peer_verification_rigor: f64,
        _rng: &mut Rng,
    ) -> (bool, ExosomaticState) {
        if !self.cultural_state.is_cumulative_ratchet() {
            return (false, self.exosomatic_state);
        }
        let alpha = self.kind.alphabet_size().max(64) as u16;
        let deme_seed = (self.deme_id as u16) * 7;

        match self.exosomatic_state {
            ExosomaticState::EndosomaticBound => {
                self.exosomatic_state = ExosomaticState::MaterialNotation;
                self.cognitive_offloading_fraction = 0.12;
                self.neural_learning_capacity = (self.neural_learning_capacity * 0.94).max(0.55);
                self.archive_false_discovery_load = 0.04;
                self.archive_sequence = (0..8u16)
                    .map(|k| (31 + deme_seed + k * 3) % alpha)
                    .collect();
                (true, self.exosomatic_state)
            }
            ExosomaticState::MaterialNotation => {
                if self.cultural_state != CulturalState::SymbolicCumulativeCulture {
                    return (false, self.exosomatic_state);
                }
                self.exosomatic_state = ExosomaticState::WrittenArchive;
                self.cognitive_offloading_fraction = 0.28;
                // Holocene cognitive-offloading encephalization relief (DeSilva et al., 2021; PMC12198904)!
                self.neural_learning_capacity = 0.74;
                self.archive_false_discovery_load = 0.05;
                self.archive_sequence = (0..22u16)
                    .map(|k| (17 + deme_seed + k * 5) % alpha)
                    .collect();
                (true, self.exosomatic_state)
            }
            ExosomaticState::WrittenArchive => {
                if self.cultural_state != CulturalState::SymbolicCumulativeCulture
                    || peer_verification_rigor < 0.65
                {
                    return (false, self.exosomatic_state);
                }
                self.exosomatic_state = ExosomaticState::PeerVerifiedInstitution;
                self.cognitive_offloading_fraction = 0.40;
                self.neural_learning_capacity = 0.68;
                self.archive_false_discovery_load = 0.01;
                self.archive_sequence = (0..38u16)
                    .map(|k| (11 + deme_seed + k * 7) % alpha)
                    .collect();
                (true, self.exosomatic_state)
            }
            ExosomaticState::PeerVerifiedInstitution => (false, self.exosomatic_state),
        }
    }

    /// Executes 1 tick of **Institutional Peer-Verification & Replication vs. Smaldino–McElreath
    /// False-Discovery Archive Pollution** (Smaldino & McElreath, 2016; Tiokhin et al., 2019; Kohrt et al., 2023):
    /// - High publication/deposit volume pressure (`unverified_publish_pressure`) injects unverified claims
    ///   into the shared external archive (`fdr_influx`).
    /// - Institutional peer audit & replication (`peer_verification_rigor`) retracts false positives and
    ///   accrues verified archival theorems (up to `48` loci).
    /// - If `peer_verification_rigor` is ablated (`< 0.20`) under high publish pressure, `archive_false_discovery_load`
    ///   crosses `0.70`, triggering **Epistemic Archive Corruption Collapse** in cognitive-offloaded populations!
    pub fn step_archive_verification_and_fdr_dynamics(
        &mut self,
        unverified_publish_pressure: f64,
        peer_verification_rigor: f64,
    ) -> (u32, bool) {
        if !self.exosomatic_state.is_archival() {
            return (0, false);
        }
        let pub_p = unverified_publish_pressure.clamp(0.0, 1.0);
        let rigor = peer_verification_rigor.clamp(0.0, 1.0);

        let fdr_influx = 0.022 + pub_p * (0.12 + 0.18 * self.cognitive_offloading_fraction);
        let inst_coeff = if self.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
            0.35
        } else {
            0.22
        };
        let audit_suppression = rigor * inst_coeff;

        let audits = if rigor >= 0.25 && (self.archive_false_discovery_load > 0.01 || pub_p > 0.12) {
            1u32 + u32::from(
                rigor >= 0.70 && self.exosomatic_state == ExosomaticState::PeerVerifiedInstitution,
            )
        } else {
            0
        };
        self.cumulative_peer_verifications =
            self.cumulative_peer_verifications.saturating_add(audits);

        self.archive_false_discovery_load =
            (self.archive_false_discovery_load + fdr_influx - audit_suppression).clamp(0.0, 1.0);

        if rigor >= 0.70
            && self.archive_false_discovery_load < 0.12
            && self.exosomatic_state == ExosomaticState::PeerVerifiedInstitution
            && self.archive_sequence.len() < 48
        {
            let alpha = self.kind.alphabet_size().max(64) as u16;
            let next_sym = ((self.archive_sequence.len() as u16) * 13 + (self.deme_id as u16) * 5 + 3) % alpha;
            self.archive_sequence.push(next_sym);
        }

        if self.archive_false_discovery_load >= 0.70 && self.exosomatic_state.has_written_archive() {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (audits, true);
        }
        (audits, false)
    }

    /// Formally verified program & theorem information (in bits) maintained by computational symbionts,
    /// scaled by preserved tail diversity (`distribution_tail_diversity`) and penalized by autophagy!
    pub fn verified_theorem_bits(&self) -> f64 {
        if self.computational_symbiont_sequence.is_empty() {
            return 0.0;
        }
        let raw = compute_raw_structure_bits(
            &self.computational_symbiont_sequence,
            self.kind.alphabet_size().max(16),
            0,
        );
        let tail_health = self.distribution_tail_diversity.clamp(0.05, 1.0)
            * (1.0 - 0.90 * self.recursive_autophagy_load).clamp(0.05, 1.0);
        let mult = match self.machine_synthesis_state {
            MachineSynthesisState::MachineVerifiedSymbiosis => 2.45,
            MachineSynthesisState::UnverifiedRecursiveGen => 1.15,
            MachineSynthesisState::AlgorithmicMediation => 1.25,
            MachineSynthesisState::ManualArchival => 0.0,
        };
        (raw * mult * tail_health).max(0.0)
    }

    /// Advances an archival lineage along **Autonomous Machine-Verified Exosomatic Synthesis**
    /// (`ManualArchival -> AlgorithmicMediation -> UnverifiedRecursiveGen / MachineVerifiedSymbiosis`;
    /// Romera-Paredes et al., 2024; Trinh et al., 2024; Shumailov et al., 2024; PMC13487184):
    /// - Requires `ExosomaticState::has_written_archive()` to host computational symbionts.
    /// - Reaching `MachineVerifiedSymbiosis` strictly requires `ExosomaticState::PeerVerifiedInstitution`
    ///   plus **formal evaluator strictness `>= 0.72`** AND **empirical grounding anchor `>= 0.12`**;
    ///   otherwise generative scaling lands in `UnverifiedRecursiveGen` (exposed to model collapse!).
    pub fn step_machine_synthesis_transition(
        &mut self,
        formal_verifier_strictness: f64,
        empirical_grounding_floor: f64,
        _rng: &mut Rng,
    ) -> (bool, MachineSynthesisState) {
        if !self.exosomatic_state.has_written_archive() {
            return (false, self.machine_synthesis_state);
        }
        let alpha = self.kind.alphabet_size().max(16) as u16;
        let d_seed = (self.deme_id as u16) * 11;
        let safeguarded = formal_verifier_strictness >= 0.72
            && empirical_grounding_floor >= 0.12
            && self.exosomatic_state == ExosomaticState::PeerVerifiedInstitution;

        match self.machine_synthesis_state {
            MachineSynthesisState::ManualArchival => {
                self.machine_synthesis_state = MachineSynthesisState::AlgorithmicMediation;
                self.distribution_tail_diversity = 0.95;
                self.recursive_autophagy_load = 0.02;
                self.computational_symbiont_sequence = (0..10u16)
                    .map(|k| (7 + d_seed + k * 5) % alpha)
                    .collect();
                (true, self.machine_synthesis_state)
            }
            MachineSynthesisState::AlgorithmicMediation => {
                if safeguarded {
                    self.machine_synthesis_state = MachineSynthesisState::MachineVerifiedSymbiosis;
                    self.distribution_tail_diversity = 0.98;
                    self.recursive_autophagy_load = 0.01;
                    self.computational_symbiont_sequence = (0..38u16)
                        .map(|k| (3 + d_seed + k * 9) % alpha)
                        .collect();
                } else {
                    self.machine_synthesis_state = MachineSynthesisState::UnverifiedRecursiveGen;
                    self.distribution_tail_diversity = 0.68;
                    self.recursive_autophagy_load = 0.18;
                    self.computational_symbiont_sequence = (0..22u16)
                        .map(|k| (5 + d_seed + k * 3) % alpha)
                        .collect();
                }
                (true, self.machine_synthesis_state)
            }
            MachineSynthesisState::UnverifiedRecursiveGen => {
                if !safeguarded {
                    return (false, self.machine_synthesis_state);
                }
                self.machine_synthesis_state = MachineSynthesisState::MachineVerifiedSymbiosis;
                self.distribution_tail_diversity = 0.97;
                self.recursive_autophagy_load = 0.01;
                self.computational_symbiont_sequence = (0..40u16)
                    .map(|k| (3 + d_seed + k * 9) % alpha)
                    .collect();
                (true, self.machine_synthesis_state)
            }
            MachineSynthesisState::MachineVerifiedSymbiosis => (false, self.machine_synthesis_state),
        }
    }

    /// Executes 1 tick of **Shumailov (2024) Recursive Model Collapse vs. FunSearch / AlphaGeometry
    /// Formal Proof & Empirical Grounding Safeguards**:
    /// - High synthetic generation (`synthetic_generation_rate`) without formal verification
    ///   (`1.0 - formal_verifier_strictness`) or without a preserved real-world empirical anchor
    ///   (`empirical_grounding_floor < 0.12`, Shumailov `\gamma_i -> 0`) drives **early tail loss**
    ///   (`distribution_tail_diversity -> 0.0`) followed by **late variance collapse**
    ///   (`recursive_autophagy_load >= 0.70`), corrupting the coupled archive and collapsing integrity!
    /// - Conversely, pairing deterministic formal evaluators (`>= 0.72`) with empirical grounding (`>= 0.12`)
    ///   suppresses autophagy to `0.0`, preserves distribution tails (`>= 0.96`), and discovers novel
    ///   machine-verified theorems up to `64` loci!
    pub fn step_recursive_synthesis_and_collapse_safeguards(
        &mut self,
        synthetic_generation_rate: f64,
        formal_verifier_strictness: f64,
        empirical_grounding_floor: f64,
    ) -> (u32, bool) {
        if !self.machine_synthesis_state.has_computational_symbiont() {
            return (0, false);
        }
        let gen_r = synthetic_generation_rate.clamp(0.0, 1.0);
        let ver_s = formal_verifier_strictness.clamp(0.0, 1.0);
        let emp_g = empirical_grounding_floor.clamp(0.0, 1.0);

        let unverified_leak = 0.30 * (1.0 - ver_s);
        let ungrounded_leak = if emp_g < 0.14 {
            0.26 * ((0.14 - emp_g) / 0.14)
        } else {
            0.0
        };
        let autophagy_influx = 0.020 + gen_r * (unverified_leak + ungrounded_leak);

        let safeguard_suppression = if ver_s >= 0.72 && emp_g >= 0.12 {
            0.28 * ver_s + 0.24 * emp_g.min(0.50)
        } else {
            0.07 * ver_s * (emp_g / 0.12).clamp(0.0, 1.0)
        };

        let proofs = if ver_s >= 0.70 && gen_r >= 0.10 {
            1u32 + u32::from(
                ver_s >= 0.85
                    && emp_g >= 0.15
                    && self.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis,
            )
        } else {
            0
        };
        self.cumulative_formal_proofs = self.cumulative_formal_proofs.saturating_add(proofs);

        self.recursive_autophagy_load =
            (self.recursive_autophagy_load + autophagy_influx - safeguard_suppression)
                .clamp(0.0, 1.0);
        self.distribution_tail_diversity =
            (1.0 - 0.94 * self.recursive_autophagy_load).clamp(0.02, 0.99);

        if ver_s >= 0.75
            && emp_g >= 0.15
            && self.recursive_autophagy_load < 0.08
            && self.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis
            && self.computational_symbiont_sequence.len() < 64
        {
            let alpha = self.kind.alphabet_size().max(16) as u16;
            let sym = ((self.computational_symbiont_sequence.len() as u16) * 17
                + (self.deme_id as u16) * 7
                + 5)
                % alpha;
            self.computational_symbiont_sequence.push(sym);
        }

        if self.recursive_autophagy_load >= 0.70
            && matches!(
                self.machine_synthesis_state,
                MachineSynthesisState::UnverifiedRecursiveGen
                    | MachineSynthesisState::MachineVerifiedSymbiosis
            )
        {
            self.archive_false_discovery_load =
                (self.archive_false_discovery_load + 0.45).clamp(0.0, 1.0);
            self.integrity = (self.integrity * 0.20).max(0.02);
            return (proofs, true);
        }
        (proofs, false)
    }

    /// Computes this unit's **Guild-Specific Planetary Biogeochemical Fluxes**
    /// `(autotroph_o2_flux, syntroph_recycle_flux, superorg_weather_flux, unguided_techno_burn, verified_stewardship_flux)`
    /// across autotrophic, syntrophic, superorganismal, and machine-verified guilds!
    pub fn compute_planetary_guild_fluxes(
        &mut self,
        controls: &EnvironmentControls,
        biosphere: &PlanetaryBiosphere,
    ) -> (f64, f64, f64, f64, f64) {
        self.planetary_epoch = biosphere.state;
        let auto_o2 = match self.guild {
            TrophicGuild::Producer => 1.15,
            TrophicGuild::Holobiont if !self.has_apicoplast_flag => {
                0.95 + 0.22 * (self.endo_tier as f64)
            }
            TrophicGuild::PredatorParasite if self.has_active_kleptoplast() => 0.45,
            _ => 0.0,
        };
        let syntroph_rec = match self.guild {
            TrophicGuild::Syntroph => 1.35,
            TrophicGuild::Holobiont => {
                0.35 + if self.mro_class == MroClass::Hydrogenosome {
                    0.45
                } else {
                    0.15
                }
            }
            _ => 0.08,
        };
        let superorg_weath = if self.agriculture_state == AgricultureState::LeafCutterAgroEcosystem {
            1.45
        } else if self.eusocial_state == EusocialState::ObligateSuperorganism {
            0.95
        } else if self.eusocial_state.is_social() {
            0.35
        } else {
            0.0
        };
        let techno_scale = if self.machine_synthesis_state.has_computational_symbiont() {
            1.25
        } else if self.exosomatic_state.has_written_archive() {
            0.65
        } else {
            0.0
        };
        let verified_steward =
            if self.machine_synthesis_state == MachineSynthesisState::MachineVerifiedSymbiosis {
                1.40 * self.distribution_tail_diversity
            } else if self.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                0.65
            } else {
                0.0
            };
        if verified_steward > 0.25 && controls.planetary_stewardship_coupling >= 0.60 {
            self.cumulative_stewardship_actions =
                self.cumulative_stewardship_actions.saturating_add(1);
        }

        self.guild_o2_flux = auto_o2;
        self.guild_nutrient_recycling = syntroph_rec + 0.45 * superorg_weath;
        self.guild_carbon_drawdown = 0.60 * auto_o2
            + 0.45 * superorg_weath
            + controls.planetary_stewardship_coupling * verified_steward
            - controls.extractive_burn_intensity * techno_scale;

        (
            auto_o2,
            syntroph_rec,
            superorg_weath,
            techno_scale,
            verified_steward,
        )
    }

    /// Applies **Planetary Boundary Breach / Asymptotic Burnout Stress** (Wong & Bartlett, 2022;
    /// Richardson et al., 2023) when `boundary_transgression_load >= 0.70`!
    pub fn apply_planetary_boundary_stress(&mut self, biosphere: &PlanetaryBiosphere) -> bool {
        self.planetary_epoch = biosphere.state;
        if biosphere.boundary_transgression_load >= 0.70 {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return true;
        }
        false
    }

    /// Initializes this unit as an abiotic Hadean prebiotic monomer-vesicle candidate (`~4.2 Ga`).
    pub fn seed_hadean_prebiotic_pool(&mut self, bath: &HadeanGeochemicalBath) {
        self.abiogenesis_stage = AbiogenesisStage::PrebioticMonomerSoup;
        self.planetary_epoch = PlanetaryState::AnoxicArchean;
        self.rna_oligomer_length_nt = 8;
        self.raf_closure_score = bath.raf_catalytic_closure.min(0.32);
        self.vesicle_membrane_stability = (0.35 + 0.55 * bath.citrate_chelation_fraction).clamp(0.15, 0.95);
        self.progenote_hgt_flux = 0.85;
    }

    /// Advances **Early-Earth Hadean Abiogenesis** across 4 biophysical & information-theoretic gates:
    /// 1. **Water-Paradox Condensation & Citrate-`Mg2+` Protocell Encapsulation** (Damer & Deamer, 2020; Adamala & Szostak, 2013)
    ///    -> `AutocatalyticRafProtocell` (`RAF >= 0.52`, `K+/Na+ >= 0.50`, chelated `Mg2+`).
    /// 2. **Eigen's Error Catastrophe & Ribozyme Hypercycle** (Eigen, 1971; Papastavrou et al., 2024)
    ///    -> `RnaRibozymeHypercycle` (`L >= 48 nt` AND `L <= L_eigen_max`).
    /// 3. **Woese Communal Progenote** (`L >= 165 nt`, proto-ribosomal peptide cofactors, high horizontal `HGT`).
    /// 4. **Darwinian Threshold -> LUCA** (`L >= 420 nt`, vertical RNP/DNA fidelity `mu <= 0.004`, `HGT` crystallization `<= 0.18`).
    /// Returns `(new_stage, eigen_error_catastrophe, mg_vesicle_precipitated)`.
    pub fn step_hadean_abiogenesis(
        &mut self,
        bath: &HadeanGeochemicalBath,
        rng: &mut Rng,
    ) -> (AbiogenesisStage, bool, bool) {
        // 1. Check unchelated Mg2+ fatty-acid vesicle precipitation & RNA strand cleavage (Adamala & Szostak, 2013):
        if self.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage
            && bath.mg2_molar_mm >= 25.0
            && bath.citrate_chelation_fraction < 0.35
        {
            self.vesicle_membrane_stability = (self.vesicle_membrane_stability * 0.35).max(0.04);
            self.integrity = (self.integrity * 0.24).max(0.02);
            return (self.abiogenesis_stage, false, true);
        }

        self.vesicle_membrane_stability = (0.30
            + 0.48 * bath.citrate_chelation_fraction
            + 0.22 * (bath.k_to_na_ratio / 2.0).clamp(0.0, 1.0))
        .clamp(0.15, 0.99);

        // 2. Water-paradox check: dilute bulk water without wet-dry films or vent Delta-pH stalls at <= 10 nt!
        let cond_drive = bath.condensation_exergy_drive();
        if cond_drive < 0.24 && self.abiogenesis_stage == AbiogenesisStage::PrebioticMonomerSoup {
            self.rna_oligomer_length_nt = self.rna_oligomer_length_nt.min(10);
            self.raf_closure_score = (self.raf_closure_score * 0.92).max(0.08);
            return (self.abiogenesis_stage, false, false);
        }

        // 3. Polymer elongation & RAF catalytic closure amplification:
        self.raf_closure_score = (0.55 * self.raf_closure_score
            + 0.45 * bath.raf_catalytic_closure
            + 0.06 * cond_drive)
            .clamp(0.08, 0.99);

        let step_elong = match self.abiogenesis_stage {
            AbiogenesisStage::PrebioticMonomerSoup => (14.0 * cond_drive) as u16 + 4,
            AbiogenesisStage::AutocatalyticRafProtocell => (28.0 * cond_drive) as u16 + 10,
            AbiogenesisStage::RnaRibozymeHypercycle => (68.0 * cond_drive) as u16 + 24,
            AbiogenesisStage::ProgenoteCommunity => (145.0 * cond_drive) as u16 + 55,
            AbiogenesisStage::LucaCellularLineage => 18,
        };
        let jitter = (rng.next_f64() * 4.0) as u16;
        self.rna_oligomer_length_nt = self
            .rna_oligomer_length_nt
            .saturating_add(step_elong + jitter)
            .min(2400);

        // 4. Eigen's Error Catastrophe (Eigen Cliff: L > ln(sigma) / mu melts quasispecies sequence info!):
        if (self.rna_oligomer_length_nt as usize) > bath.eigen_max_length_nt
            && self.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage
        {
            self.rna_oligomer_length_nt = (bath.eigen_max_length_nt as u16).max(12) / 2;
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (self.abiogenesis_stage, true, false);
        }

        // 5. Stepwise Abiogenesis Phase Transitions:
        match self.abiogenesis_stage {
            AbiogenesisStage::PrebioticMonomerSoup => {
                if self.raf_closure_score >= 0.50
                    && self.vesicle_membrane_stability >= 0.58
                    && self.rna_oligomer_length_nt >= 18
                {
                    self.abiogenesis_stage = AbiogenesisStage::AutocatalyticRafProtocell;
                    self.progenote_hgt_flux = 0.78;
                }
            }
            AbiogenesisStage::AutocatalyticRafProtocell => {
                if self.rna_oligomer_length_nt >= 48
                    && bath.eigen_max_length_nt >= 64
                    && self.raf_closure_score >= 0.62
                {
                    self.abiogenesis_stage = AbiogenesisStage::RnaRibozymeHypercycle;
                    self.progenote_hgt_flux = 0.68;
                }
            }
            AbiogenesisStage::RnaRibozymeHypercycle => {
                if self.rna_oligomer_length_nt >= 150
                    && bath.k_to_na_ratio >= 0.45
                    && bath.eigen_max_length_nt >= 200
                {
                    self.abiogenesis_stage = AbiogenesisStage::ProgenoteCommunity;
                    self.progenote_hgt_flux = 0.58;
                }
            }
            AbiogenesisStage::ProgenoteCommunity => {
                if self.rna_oligomer_length_nt >= 380
                    && bath.per_base_error_mu <= 0.0045
                    && bath.eigen_max_length_nt >= 500
                {
                    self.abiogenesis_stage = AbiogenesisStage::LucaCellularLineage;
                    self.progenote_hgt_flux = 0.12;
                    for slot in self.sequence.iter_mut().take(12) {
                        *slot = (*slot + 1) % 8;
                    }
                }
            }
            AbiogenesisStage::LucaCellularLineage => {}
        }

        (self.abiogenesis_stage, false, false)
    }

    /// Initializes this candidate in a **racemic 50:50 D/L prebiotic pool** prior to chiral symmetry breaking.
    pub fn seed_racemic_prebiotic_pool(&mut self, reactor: &ChiralCodeReactor) {
        self.chirality_code_stage = ChiralityCodeStage::RacemicCrossInhibited;
        self.ribose_ee = reactor.ribose_ee_d.min(0.18);
        self.peptide_l_ee = reactor.amino_ee_l.min(0.10);
        self.genetic_code_mp = 0.16;
        self.amino_acid_alphabet_size = 0;
    }

    /// Advances **Prebiotic Homochiral Symmetry Breaking, Stereoselective Minihelix Aminoacylation,
    /// and Error-Minimizing Genetic Code (`SGC`) Crystallization**:
    /// 1. **Racemic Enantiomeric Cross-Inhibition** (`ribose_ee < 0.35`): `L-ribonucleotides` cap `3'`
    ///    ends and terminate RNA template extension (`racemic_terminated = true`).
    /// 2. **Magnetite `CISS` + Frank Amplification** (`ribose_ee >= 0.80` -> `HomochiralRnaSymmetryBroken`).
    /// 3. **Tamura-Schimmel `3'-NCCA` Minihelix Steric Filter** (`peptide_l_ee >= 0.72` -> `StereoselectiveMinihelixAcylated`).
    /// 4. **Woese Statistical Code vs. Babel Trap**: Isolated vertical lineages (`code_universality < 0.35`)
    ///    freeze into fragmented suboptimal codes (`MP <= 0.44`), suffering mistranslation load (`babel_mistranslated = true`).
    /// 5. **Collective `HGT` SGC Crystallization** (`MP >= 0.74`, `20 L-AA`, `consensus >= 0.65` -> `UniversalErrorMinimizedSgc`)!
    pub fn step_homochiral_code_crystallization(
        &mut self,
        reactor: &ChiralCodeReactor,
        rng: &mut Rng,
    ) -> (ChiralityCodeStage, bool, bool) {
        self.ribose_ee = reactor.ribose_ee_d;
        self.peptide_l_ee = reactor.amino_ee_l;
        self.genetic_code_mp = reactor.sgc_error_minimization_mp;
        self.amino_acid_alphabet_size = reactor.encoded_amino_acids;

        // 1. Joyce-Orgel Enantiomeric Cross-Inhibition in racemic baths (`ee_D < 0.35`):
        if self.chirality_code_stage != ChiralityCodeStage::UniversalErrorMinimizedSgc
            && reactor.cross_inhibition_prob >= 0.38
        {
            self.rna_oligomer_length_nt = self.rna_oligomer_length_nt.min(12);
            self.integrity = (self.integrity * 0.24).max(0.02);
            return (self.chirality_code_stage, true, false);
        }

        // 2. Stepwise Chiral & Genetic Code Transitions:
        match self.chirality_code_stage {
            ChiralityCodeStage::RacemicCrossInhibited => {
                if self.ribose_ee >= 0.78 {
                    self.chirality_code_stage = ChiralityCodeStage::HomochiralRnaSymmetryBroken;
                }
            }
            ChiralityCodeStage::HomochiralRnaSymmetryBroken => {
                if self.peptide_l_ee >= 0.70 && self.amino_acid_alphabet_size >= 4 {
                    self.chirality_code_stage =
                        ChiralityCodeStage::StereoselectiveMinihelixAcylated;
                }
            }
            ChiralityCodeStage::StereoselectiveMinihelixAcylated => {
                if self.amino_acid_alphabet_size >= 10 && self.genetic_code_mp >= 0.32 {
                    self.chirality_code_stage = ChiralityCodeStage::AmbiguousStatisticalCode;
                }
            }
            ChiralityCodeStage::AmbiguousStatisticalCode => {
                // Check Vetsigian-Woese-Goldenfeld Babel Fragmentation Trap when communal HGT is absent:
                if reactor.code_universality_consensus < 0.35 && self.genetic_code_mp < 0.46 {
                    let jitter = rng.next_f64();
                    if jitter < 0.85 {
                        self.integrity = (self.integrity * 0.26).max(0.02);
                        return (self.chirality_code_stage, false, true);
                    }
                }
                if reactor.sgc_crystallized {
                    self.chirality_code_stage = ChiralityCodeStage::UniversalErrorMinimizedSgc;
                    for slot in self.sequence.iter_mut().take(10) {
                        *slot = (*slot + 1) % 8;
                    }
                }
            }
            ChiralityCodeStage::UniversalErrorMinimizedSgc => {}
        }

        (self.chirality_code_stage, false, false)
    }

    /// Initializes this lineage as an **ancestral ribocell (`2'-OH` labile RNA chromosome, ~16 kb)**
    /// coevolving with the primordial virosphere prior to `RNR` deoxyribonucleotide synthesis!
    pub fn seed_labile_rna_cell(&mut self, reactor: &RnrVirosphereReactor) {
        self.genome_substrate_stage = GenomeSubstrateStage::LabileRnaChromosome;
        self.chromosome_capacity_kb = 16.0;
        self.deamination_uracil_load = reactor.deamination_lesion_load.max(0.18);
        self.antiviral_defense_titer =
            (0.22 + 0.68 * reactor.virosphere_coevolution_coupling).clamp(0.10, 0.96);
    }

    /// Advances **Ribonucleotide Reductase (`RNR`) `RNA -> U-DNA -> T-DNA` Transition,
    /// `ThyA/ThyX + UDG` Cytosine Deamination Repair, and Virosphere Arms-Race Coevolution**:
    /// 1. **`2'-OH` Hydrolysis Ceiling** (`rnr < 0.30`): RNA chromosomes stay capped `< 28 kb`.
    /// 2. **Unrepaired `C -> U` Deamination Meltdown in `U-DNA`** (`deamination_lesion_load >= 0.52`
    ///    once genomes expand into `U-DNA` without `ThyX + UDG` excision -> `deamination_melted = true`).
    /// 3. **Lytic Virosphere Epidemic Collapse** (`lytic_viral_pressure > antiviral_defense_titer + 0.40`
    ///    -> `lytic_lysed = true`).
    /// 4. **Viral Polymerase Domestication & Megabase `T-DNA` Crown** (`CoevolvedCellularVirosphere`)!
    pub fn step_rnr_dna_and_virosphere_coevolution(
        &mut self,
        reactor: &RnrVirosphereReactor,
        _rng: &mut Rng,
    ) -> (GenomeSubstrateStage, bool, bool) {
        self.deamination_uracil_load = reactor.deamination_lesion_load;
        self.antiviral_defense_titer =
            (0.20 + 0.75 * reactor.virosphere_coevolution_coupling).clamp(0.10, 0.98);

        // 1. Check runaway lytic phage epidemic when host defense / arms-race coupling is overwhelmed:
        if self.genome_substrate_stage != GenomeSubstrateStage::CoevolvedCellularVirosphere
            && reactor.lytic_viral_pressure > self.antiviral_defense_titer + 0.40
        {
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (self.genome_substrate_stage, false, true);
        }

        self.chromosome_capacity_kb = reactor.stable_genome_capacity_kb;

        // 2. Check Cytosine Deamination (`C -> U`) Mutational Meltdown in expanded U-DNA lacking Thymine + UDG!
        if self.genome_substrate_stage == GenomeSubstrateStage::ViralModifiedUDna
            && reactor.deamination_lesion_load >= 0.50
        {
            self.chromosome_capacity_kb = self.chromosome_capacity_kb.min(42.0);
            self.integrity = (self.integrity * 0.24).max(0.02);
            return (self.genome_substrate_stage, true, false);
        }

        // 3. Stepwise Genome Substrate & Virosphere Transitions:
        match self.genome_substrate_stage {
            GenomeSubstrateStage::LabileRnaChromosome => {
                if reactor.rnr_thiyl_radical_activity >= 0.45
                    && reactor.stable_genome_capacity_kb >= 45.0
                {
                    self.genome_substrate_stage = GenomeSubstrateStage::ViralModifiedUDna;
                }
            }
            GenomeSubstrateStage::ViralModifiedUDna => {
                if reactor.thymidylate_udg_repair >= 0.55
                    && reactor.deamination_lesion_load <= 0.22
                    && reactor.stable_genome_capacity_kb >= 220.0
                {
                    self.genome_substrate_stage = GenomeSubstrateStage::DeaminationRepairableTDna;
                }
            }
            GenomeSubstrateStage::DeaminationRepairableTDna => {
                if reactor.megabase_tdna_established {
                    self.genome_substrate_stage =
                        GenomeSubstrateStage::CoevolvedCellularVirosphere;
                    for slot in self.sequence.iter_mut().take(12) {
                        *slot = (*slot + 1) % 8;
                    }
                }
            }
            GenomeSubstrateStage::CoevolvedCellularVirosphere => {}
        }

        (self.genome_substrate_stage, false, false)
    }

    pub fn egt_nuclear_fraction(&self) -> f64 {
        let tot = self.egt_transferred_loci + self.endosymbiont_sequence.len();
        if tot == 0 || (self.endosymbiont_sequence.is_empty() && self.egt_transferred_loci == 0) {
            0.0
        } else {
            ((self.egt_transferred_loci as f64) / (tot as f64)).clamp(0.0, 1.0)
        }
    }

    pub fn is_organelle_enslaved(&self) -> bool {
        self.egt_transferred_loci >= 2
            || self.egt_nuclear_fraction() >= 0.50
            || self.endo_tier >= 2
            || self.mro_class.is_reduced_anaerobic_mro()
    }

    pub fn nucleomorph_plastid_loci_remaining(&self) -> usize {
        self.nucleomorph_sequence
            .len()
            .saturating_sub(NUCLEOMORPH_HOUSEKEEPING_LOCI)
    }

    fn insert_locus_into_host_nucleus(&mut self, sym: u16, rng: &mut Rng) {
        let ins_min = usize::from(self.kind == OrganismKind::ArtificialIntelligence);
        let ins_max = if self.kind == OrganismKind::ArtificialIntelligence {
            self.sequence.len().saturating_sub(1).max(ins_min)
        } else {
            self.sequence.len()
        };

        if self.sequence.len() < 48 {
            let span = (ins_max - ins_min).max(1);
            let dst = ins_min + rng.gen_index(span + 1);
            let clamped_dst = dst.min(ins_max);
            self.sequence.insert(clamped_dst, sym);
            let mem_dst = clamped_dst.min(self.lifetime_memory.len());
            self.lifetime_memory.insert(mem_dst, 0.44);
        } else if ins_max > ins_min {
            let dst = ins_min + rng.gen_index(ins_max - ins_min);
            self.sequence[dst] = sym;
        }
        self.wiring = Self::build_wiring(self.kind, &self.sequence, RegimeId::R35Symbiogenetic);
    }

    pub fn convert_holobiont_to_apicoplast_parasite(&mut self, rng: &mut Rng) -> bool {
        if self.guild != TrophicGuild::Holobiont || self.has_apicoplast_flag {
            return false;
        }

        let suf_b_locus = self.endosymbiont_sequence.first().copied().unwrap_or(1);
        self.endosymbiont_sequence.clear();
        self.endosymbiont_sequence.push(suf_b_locus);

        if !self.nucleomorph_sequence.is_empty() {
            self.nucleomorph_sequence.clear();
            self.nucleomorph_eliminated = true;
        }

        self.endo_tier = self.endo_tier.max(2);
        self.organelle_membranes = 4;
        self.has_apicoplast_flag = true;
        self.apicoplast_eliminated = false;
        self.apicoplast_pathways_mask = APICO_PATH_ALL;
        if self.mro_class == MroClass::None {
            self.mro_class = MroClass::AerobicMitochondrion;
            self.mro_pathways_mask = MRO_PATH_ALL;
            self.mro_mtdna_loci = CORR_CORE_LOCI;
        }
        self.egt_transferred_loci = self.egt_transferred_loci.max(4);
        self.kleptoplast_sequence.clear();
        self.kleptoplast_integrity = 0.0;

        let inv_motifs: &[u16] = match self.kind {
            OrganismKind::RnaDnaProtein => &[2, 3, 9, 10],
            OrganismKind::HominidPopulation => &[1, 3, 5],
            OrganismKind::ArtificialIntelligence => &[3, 4],
        };
        let start_i = usize::from(self.kind == OrganismKind::ArtificialIntelligence);
        let end_i = if self.kind == OrganismKind::ArtificialIntelligence {
            self.sequence.len().saturating_sub(1)
        } else {
            self.sequence.len()
        };
        if end_i > start_i + 2 {
            for k in 0..2 {
                let pos = start_i + rng.gen_index(end_i - start_i);
                self.sequence[pos] = inv_motifs[k % inv_motifs.len()];
            }
        }
        self.wiring = Self::build_wiring(self.kind, &self.sequence, RegimeId::R35Symbiogenetic);
        self.guild = TrophicGuild::PredatorParasite;
        true
    }

    pub fn streamline_apicoplast_pathway(
        &mut self,
        host_metabolite_rescue: f64,
        force_disrupt_essential: bool,
    ) -> (bool, bool, bool) {
        if !self.has_apicoplast() {
            return (false, false, false);
        }

        if (self.apicoplast_pathways_mask & APICO_PATH_HEME) != 0 {
            self.apicoplast_pathways_mask &= !APICO_PATH_HEME;
            return (true, false, false);
        }

        if (self.apicoplast_pathways_mask & APICO_PATH_FASII) != 0 {
            self.apicoplast_pathways_mask &= !APICO_PATH_FASII;
            return (true, false, false);
        }

        if host_metabolite_rescue >= 0.65 {
            self.apicoplast_pathways_mask = 0;
            self.endosymbiont_sequence.clear();
            self.has_apicoplast_flag = false;
            self.apicoplast_eliminated = true;
            self.organelle_membranes = 0;
            return (true, true, false);
        }

        if force_disrupt_essential {
            self.apicoplast_pathways_mask = 0;
            self.endosymbiont_sequence.clear();
            self.integrity = (self.integrity * 0.22).max(0.02);
            return (true, false, true);
        }

        (false, false, false)
    }

    pub fn steal_kleptoplast_from_prey(
        &mut self,
        prey: &ModeledUnit,
        steal_plastid: bool,
        recruit_hgt: bool,
        allow_fixation: bool,
        rng: &mut Rng,
    ) -> (bool, bool, bool) {
        if self.guild != TrophicGuild::PredatorParasite
            || self.has_apicoplast()
            || self.apicoplast_eliminated
        {
            return (false, false, false);
        }

        let donor_bit = 1u32 << ((prey.species_id as u32 + (prey.guild as u32) * 5) % 24);
        let mut did_steal = false;
        let mut did_hgt = false;

        if steal_plastid {
            let src_seq = if !prey.endosymbiont_sequence.is_empty() && !prey.has_apicoplast() {
                &prey.endosymbiont_sequence
            } else {
                &prey.sequence
            };
            if !src_seq.is_empty() {
                let kp_len = src_seq.len().clamp(6, 9).min(src_seq.len());
                self.kleptoplast_sequence = src_seq[..kp_len].to_vec();
                self.kleptoplast_integrity = 1.0;
                self.kleptoplast_donor_id = Some(prey.id);
                self.kleptoplast_donor_tier = if prey.guild == TrophicGuild::Holobiont {
                    prey.endo_tier
                } else {
                    0
                };
                self.klepto_donor_mask |= donor_bit;
                did_steal = true;
            }
        }

        if recruit_hgt && !prey.sequence.is_empty() && self.klepto_hgt_loci < 8 {
            let syn_motifs: &[u16] = match self.kind {
                OrganismKind::RnaDnaProtein => &[1, 4, 7, 8],
                OrganismKind::HominidPopulation => &[0, 6, 7, 8],
                OrganismKind::ArtificialIntelligence => &[1, 5, 7, 8],
            };
            let picked_sym = prey
                .sequence
                .iter()
                .copied()
                .find(|s| syn_motifs.contains(s))
                .unwrap_or_else(|| prey.sequence[rng.gen_index(prey.sequence.len())]);
            self.insert_locus_into_host_nucleus(picked_sym, rng);
            self.klepto_hgt_loci = self.klepto_hgt_loci.saturating_add(1);
            self.klepto_donor_mask |= donor_bit;
            if self.has_active_kleptoplast() {
                self.kleptoplast_integrity = (self.kleptoplast_integrity + 0.18).clamp(0.0, 1.0);
            }
            did_hgt = true;
        }

        let crossed = if allow_fixation {
            self.try_cross_kleptoplasty_tipping_point()
        } else {
            false
        };
        (did_steal, did_hgt, crossed)
    }

    pub fn try_cross_kleptoplasty_tipping_point(&mut self) -> bool {
        if self.guild != TrophicGuild::PredatorParasite
            || self.has_apicoplast()
            || self.apicoplast_eliminated
            || self.klepto_hgt_loci < KLEPTO_FIXATION_HGT_THRESHOLD
            || !self.has_active_kleptoplast()
            || self.kleptoplast_integrity < 0.25
        {
            return false;
        }

        let new_tier = (self.kleptoplast_donor_tier + 1).clamp(1, 3);
        let mut fixed_org = std::mem::take(&mut self.kleptoplast_sequence);
        while fixed_org.len() < CORR_CORE_LOCI {
            fixed_org.push(1);
        }

        self.endosymbiont_sequence = fixed_org;
        self.endosymbiont_origin_id = self.kleptoplast_donor_id.take();
        self.endo_tier = new_tier;
        self.organelle_membranes = if new_tier >= 2 { 4 } else { 2 };
        if new_tier >= 2 {
            self.nucleomorph_sequence.clear();
            self.nucleomorph_eliminated = true;
        }
        self.apicoplast_pathways_mask = APICO_PATH_ALL;
        if self.mro_class == MroClass::None {
            self.mro_class = MroClass::AerobicMitochondrion;
            self.mro_pathways_mask = MRO_PATH_ALL;
            self.mro_mtdna_loci = CORR_CORE_LOCI;
        }
        self.egt_transferred_loci = self.egt_transferred_loci.saturating_add(self.klepto_hgt_loci);
        self.kleptoplast_integrity = 0.0;
        self.kleptoplast_donor_tier = 0;
        self.guild = TrophicGuild::Holobiont;
        true
    }

    pub fn migrate_endosymbiont_locus_to_host(
        &mut self,
        allow_below_corr: bool,
        rng: &mut Rng,
    ) -> bool {
        if self.has_apicoplast() {
            return false;
        }
        let floor = if allow_below_corr { 0 } else { CORR_CORE_LOCI };
        if self.endosymbiont_sequence.len() <= floor {
            return false;
        }

        let pick_idx = self.endosymbiont_sequence.len() - 1;
        let migrated_sym = self.endosymbiont_sequence.remove(pick_idx);
        self.insert_locus_into_host_nucleus(migrated_sym, rng);
        self.egt_transferred_loci = self.egt_transferred_loci.saturating_add(1);
        true
    }

    pub fn migrate_nucleomorph_locus_to_host(
        &mut self,
        allow_final_sweep: bool,
        rng: &mut Rng,
    ) -> (usize, bool) {
        if self.nucleomorph_sequence.is_empty() {
            return (0, false);
        }

        let bonsai_floor = NUCLEOMORPH_HOUSEKEEPING_LOCI + 1;
        if self.nucleomorph_sequence.len() > bonsai_floor {
            let migrated_sym = self.nucleomorph_sequence.pop().unwrap_or(1);
            self.insert_locus_into_host_nucleus(migrated_sym, rng);
            self.nucleomorph_transferred_loci = self.nucleomorph_transferred_loci.saturating_add(1);
            return (1, false);
        }

        if !allow_final_sweep {
            return (0, false);
        }

        let final_plastid_sym = self.nucleomorph_sequence.pop().unwrap_or(1);
        self.insert_locus_into_host_nucleus(final_plastid_sym, rng);
        self.nucleomorph_transferred_loci = self.nucleomorph_transferred_loci.saturating_add(1);
        self.nucleomorph_sequence.clear();
        self.nucleomorph_eliminated = true;
        (1, true)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn fuse_endosymbiotic_holobiont(
        host: &ModeledUnit,
        symbiont: &ModeledUnit,
        holobiont_id: u64,
        biome: BiomeZone,
        local_exergy: f64,
        local_byproduct: f64,
        local_niche: f64,
        controls: &EnvironmentControls,
        adami_rho: f64,
        pop_size: usize,
        env_phase: f64,
    ) -> (ModeledUnit, f64) {
        let host_seq = host.sequence.clone();
        let target_endo_len = ((symbiont.sequence.len() * 3) / 4).clamp(10, 14);
        let mut endo_seq = Vec::with_capacity(target_endo_len);
        for (idx, &sym) in symbiont.sequence.iter().enumerate() {
            if endo_seq.len() >= target_endo_len {
                break;
            }
            if idx % 4 != 3 || endo_seq.len() + (symbiont.sequence.len() - idx) <= target_endo_len {
                endo_seq.push(sym);
            }
        }
        while endo_seq.len() < 10 {
            let fallback = symbiont
                .sequence
                .get(endo_seq.len() % symbiont.sequence.len().max(1))
                .copied()
                .unwrap_or(1);
            endo_seq.push(fallback);
        }

        let gamma =
            (endo_seq.len() as f64) / ((host_seq.len() + endo_seq.len()).max(1) as f64);

        let wiring = Self::build_wiring(host.kind, &host_seq, RegimeId::R35Symbiogenetic);
        let mut lifetime_memory = vec![0.0; host_seq.len()];
        for (i, slot) in lifetime_memory.iter_mut().enumerate() {
            let hm = host.lifetime_memory.get(i).copied().unwrap_or(0.35);
            let sm = symbiont
                .lifetime_memory
                .get(i % symbiont.lifetime_memory.len().max(1))
                .copied()
                .unwrap_or(0.35);
            *slot = (0.65 * hm + 0.35 * sm).clamp(0.0, 1.35);
        }

        let mut holo = ModeledUnit {
            id: holobiont_id,
            kind: host.kind,
            generation: host.generation.max(symbiont.generation).saturating_add(1),
            age: 0,
            grid_x: host.grid_x,
            grid_y: host.grid_y,
            species_id: host.species_id,
            guild: TrophicGuild::Holobiont,
            endo_tier: 1,
            organelle_membranes: 2,
            sequence: host_seq,
            kleptoplast_sequence: Vec::new(),
            kleptoplast_integrity: 0.0,
            kleptoplast_donor_id: None,
            kleptoplast_donor_tier: 0,
            klepto_hgt_loci: host.klepto_hgt_loci,
            klepto_donor_mask: host.klepto_donor_mask,
            has_apicoplast_flag: false,
            apicoplast_pathways_mask: APICO_PATH_ALL,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways_mask: MRO_PATH_ALL,
            mro_mtdna_loci: CORR_CORE_LOCI,
            multicell_state: MulticellularState::Unicellular,
            cell_count: 1,
            soma_fraction: 0.0,
            kin_relatedness_r: 1.0,
            somatic_cheater_load: 0.0,
            cumulative_apoptotic_culls: 0,
            eusocial_state: EusocialState::Solitary,
            colony_size: 1,
            sterile_worker_fraction: 0.0,
            soldier_fraction: 0.0,
            colony_relatedness_r: 0.50,
            rebel_worker_load: 0.0,
            cumulative_policed_eggs: 0,
            agriculture_state: AgricultureState::HunterForager,
            cultivar_polyploidy: 1,
            crop_clonal_uniformity: 0.0,
            gongylidia_yield: 0.0,
            escovopsis_pathogen_load: 0.0,
            defensive_actinobacteria_cover: 0.0,
            cumulative_antimicrobial_shields: 0,
            cultural_state: host.cultural_state,
            deme_id: host.deme_id,
            cultural_sequence: host.cultural_sequence.clone(),
            cultural_ratchet_tier: host.cultural_ratchet_tier,
            neural_learning_capacity: host.neural_learning_capacity,
            cultural_drift_erosion: host.cultural_drift_erosion,
            cumulative_cross_deme_recombinations: host.cumulative_cross_deme_recombinations,
            exosomatic_state: host.exosomatic_state,
            archive_sequence: host.archive_sequence.clone(),
            cognitive_offloading_fraction: host.cognitive_offloading_fraction,
            archive_false_discovery_load: host.archive_false_discovery_load,
            cumulative_peer_verifications: host.cumulative_peer_verifications,
            machine_synthesis_state: host.machine_synthesis_state,
            computational_symbiont_sequence: host.computational_symbiont_sequence.clone(),
            distribution_tail_diversity: host.distribution_tail_diversity,
            recursive_autophagy_load: host.recursive_autophagy_load,
            cumulative_formal_proofs: host.cumulative_formal_proofs,
            planetary_epoch: host.planetary_epoch,
            guild_o2_flux: host.guild_o2_flux,
            guild_carbon_drawdown: host.guild_carbon_drawdown,
            guild_nutrient_recycling: host.guild_nutrient_recycling,
            cumulative_stewardship_actions: host.cumulative_stewardship_actions,
            abiogenesis_stage: host.abiogenesis_stage,
            rna_oligomer_length_nt: host.rna_oligomer_length_nt,
            raf_closure_score: host.raf_closure_score,
            vesicle_membrane_stability: host.vesicle_membrane_stability,
            progenote_hgt_flux: host.progenote_hgt_flux,
            chirality_code_stage: host.chirality_code_stage,
            ribose_ee: host.ribose_ee,
            peptide_l_ee: host.peptide_l_ee,
            genetic_code_mp: host.genetic_code_mp,
            amino_acid_alphabet_size: host.amino_acid_alphabet_size,
            genome_substrate_stage: host.genome_substrate_stage,
            chromosome_capacity_kb: host.chromosome_capacity_kb,
            deamination_uracil_load: host.deamination_uracil_load,
            antiviral_defense_titer: host.antiviral_defense_titer,
            nucleomorph_sequence: Vec::new(),
            nucleomorph_transferred_loci: 0,
            nucleomorph_eliminated: false,
            endosymbiont_sequence: endo_seq,
            endosymbiont_origin_id: Some(symbiont.id),
            egt_transferred_loci: 0,
            wiring,
            lifetime_memory,
            integrity: 1.0,
            structure_raw_s: 12.0,
            structure_s: 12.0,
            base_complexity_c: 12.0,
            complexity_c: 12.0,
            primary_work: 0.0,
            syntrophic_work: 0.0,
            net_trophic_flux: 0.0,
            harvested_work: 0.0,
            landauer_cost: 0.0,
            byproduct_output: 0.0,
            niche_output: 0.0,
            offense_trait: 0.0,
            defense_trait: 0.0,
            fitness: 0.8,
        };
        holo.recompute_phenome(
            controls,
            local_exergy,
            local_byproduct,
            local_niche,
            0.35,
            biome,
            adami_rho,
            pop_size,
            env_phase,
        );
        (holo, gamma)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn fuse_serial_holobiont(
        outer_host: &ModeledUnit,
        engulfed_holobiont: &ModeledUnit,
        chimera_id: u64,
        biome: BiomeZone,
        local_exergy: f64,
        local_byproduct: f64,
        local_niche: f64,
        controls: &EnvironmentControls,
        adami_rho: f64,
        pop_size: usize,
        env_phase: f64,
    ) -> (ModeledUnit, f64) {
        let new_tier = (engulfed_holobiont.endo_tier.max(1) + 1).min(3);
        let host_seq = outer_host.sequence.clone();

        let mut inner_org = engulfed_holobiont.endosymbiont_sequence.clone();
        while inner_org.len() < CORR_CORE_LOCI {
            let sym = engulfed_holobiont
                .sequence
                .get(inner_org.len() % engulfed_holobiont.sequence.len().max(1))
                .copied()
                .unwrap_or(1);
            inner_org.push(sym);
        }

        let target_nm_len = NUCLEOMORPH_HOUSEKEEPING_LOCI + 4;
        let mut nm_seq = Vec::with_capacity(target_nm_len);
        let src_len = engulfed_holobiont.sequence.len().max(1);
        for i in 0..target_nm_len {
            let pick = engulfed_holobiont.sequence[(i * 2 + 1) % src_len];
            nm_seq.push(pick);
        }

        let endo_bits = inner_org.len() + nm_seq.len();
        let gamma = (endo_bits as f64) / ((host_seq.len() + endo_bits).max(1) as f64);

        let wiring = Self::build_wiring(outer_host.kind, &host_seq, RegimeId::R35Symbiogenetic);
        let mut lifetime_memory = vec![0.0; host_seq.len()];
        for (i, slot) in lifetime_memory.iter_mut().enumerate() {
            let hm = outer_host.lifetime_memory.get(i).copied().unwrap_or(0.40);
            let em = engulfed_holobiont
                .lifetime_memory
                .get(i % engulfed_holobiont.lifetime_memory.len().max(1))
                .copied()
                .unwrap_or(0.40);
            *slot = (0.60 * hm + 0.40 * em).clamp(0.0, 1.40);
        }

        let inherited_egt = outer_host
            .egt_transferred_loci
            .max(engulfed_holobiont.egt_transferred_loci.min(6));
        let host_mro = if outer_host.mro_class != MroClass::None {
            outer_host.mro_class
        } else {
            MroClass::AerobicMitochondrion
        };
        let host_mro_paths = if outer_host.mro_class != MroClass::None {
            outer_host.mro_pathways_mask
        } else {
            MRO_PATH_ALL
        };
        let host_mtdna = if outer_host.mro_class != MroClass::None {
            outer_host.mro_mtdna_loci
        } else {
            CORR_CORE_LOCI
        };

        let mut chimera = ModeledUnit {
            id: chimera_id,
            kind: outer_host.kind,
            generation: outer_host
                .generation
                .max(engulfed_holobiont.generation)
                .saturating_add(1),
            age: 0,
            grid_x: outer_host.grid_x,
            grid_y: outer_host.grid_y,
            species_id: outer_host.species_id,
            guild: TrophicGuild::Holobiont,
            endo_tier: new_tier,
            organelle_membranes: 4,
            sequence: host_seq,
            kleptoplast_sequence: Vec::new(),
            kleptoplast_integrity: 0.0,
            kleptoplast_donor_id: None,
            kleptoplast_donor_tier: 0,
            klepto_hgt_loci: outer_host.klepto_hgt_loci,
            klepto_donor_mask: outer_host.klepto_donor_mask | engulfed_holobiont.klepto_donor_mask,
            has_apicoplast_flag: false,
            apicoplast_pathways_mask: APICO_PATH_ALL,
            apicoplast_eliminated: false,
            mro_class: host_mro,
            mro_pathways_mask: host_mro_paths,
            mro_mtdna_loci: host_mtdna,
            multicell_state: outer_host.multicell_state,
            cell_count: outer_host.cell_count.max(1),
            soma_fraction: outer_host.soma_fraction,
            kin_relatedness_r: outer_host.kin_relatedness_r,
            somatic_cheater_load: outer_host.somatic_cheater_load,
            cumulative_apoptotic_culls: outer_host.cumulative_apoptotic_culls,
            eusocial_state: outer_host.eusocial_state,
            colony_size: outer_host.colony_size.max(1),
            sterile_worker_fraction: outer_host.sterile_worker_fraction,
            soldier_fraction: outer_host.soldier_fraction,
            colony_relatedness_r: outer_host.colony_relatedness_r,
            rebel_worker_load: outer_host.rebel_worker_load,
            cumulative_policed_eggs: outer_host.cumulative_policed_eggs,
            agriculture_state: outer_host.agriculture_state,
            cultivar_polyploidy: outer_host.cultivar_polyploidy.max(1),
            crop_clonal_uniformity: outer_host.crop_clonal_uniformity,
            gongylidia_yield: outer_host.gongylidia_yield,
            escovopsis_pathogen_load: outer_host.escovopsis_pathogen_load,
            defensive_actinobacteria_cover: outer_host.defensive_actinobacteria_cover,
            cumulative_antimicrobial_shields: outer_host.cumulative_antimicrobial_shields,
            cultural_state: outer_host.cultural_state,
            deme_id: outer_host.deme_id,
            cultural_sequence: outer_host.cultural_sequence.clone(),
            cultural_ratchet_tier: outer_host.cultural_ratchet_tier,
            neural_learning_capacity: outer_host.neural_learning_capacity,
            cultural_drift_erosion: outer_host.cultural_drift_erosion,
            cumulative_cross_deme_recombinations: outer_host.cumulative_cross_deme_recombinations,
            exosomatic_state: outer_host.exosomatic_state,
            archive_sequence: outer_host.archive_sequence.clone(),
            cognitive_offloading_fraction: outer_host.cognitive_offloading_fraction,
            archive_false_discovery_load: outer_host.archive_false_discovery_load,
            cumulative_peer_verifications: outer_host.cumulative_peer_verifications,
            machine_synthesis_state: outer_host.machine_synthesis_state,
            computational_symbiont_sequence: outer_host.computational_symbiont_sequence.clone(),
            distribution_tail_diversity: outer_host.distribution_tail_diversity,
            recursive_autophagy_load: outer_host.recursive_autophagy_load,
            cumulative_formal_proofs: outer_host.cumulative_formal_proofs,
            planetary_epoch: outer_host.planetary_epoch,
            guild_o2_flux: outer_host.guild_o2_flux,
            guild_carbon_drawdown: outer_host.guild_carbon_drawdown,
            guild_nutrient_recycling: outer_host.guild_nutrient_recycling,
            cumulative_stewardship_actions: outer_host.cumulative_stewardship_actions,
            abiogenesis_stage: outer_host.abiogenesis_stage,
            rna_oligomer_length_nt: outer_host.rna_oligomer_length_nt,
            raf_closure_score: outer_host.raf_closure_score,
            vesicle_membrane_stability: outer_host.vesicle_membrane_stability,
            progenote_hgt_flux: outer_host.progenote_hgt_flux,
            chirality_code_stage: outer_host.chirality_code_stage,
            ribose_ee: outer_host.ribose_ee,
            peptide_l_ee: outer_host.peptide_l_ee,
            genetic_code_mp: outer_host.genetic_code_mp,
            amino_acid_alphabet_size: outer_host.amino_acid_alphabet_size,
            genome_substrate_stage: outer_host.genome_substrate_stage,
            chromosome_capacity_kb: outer_host.chromosome_capacity_kb,
            deamination_uracil_load: outer_host.deamination_uracil_load,
            antiviral_defense_titer: outer_host.antiviral_defense_titer,
            nucleomorph_sequence: nm_seq,
            nucleomorph_transferred_loci: 0,
            nucleomorph_eliminated: false,
            endosymbiont_sequence: inner_org,
            endosymbiont_origin_id: Some(engulfed_holobiont.id),
            egt_transferred_loci: inherited_egt,
            wiring,
            lifetime_memory,
            integrity: 1.0,
            structure_raw_s: 14.0,
            structure_s: 14.0,
            base_complexity_c: 15.0,
            complexity_c: 15.0,
            primary_work: 0.0,
            syntrophic_work: 0.0,
            net_trophic_flux: 0.0,
            harvested_work: 0.0,
            landauer_cost: 0.0,
            byproduct_output: 0.0,
            niche_output: 0.0,
            offense_trait: 0.0,
            defense_trait: 0.0,
            fitness: 0.9,
        };
        chimera.recompute_phenome(
            controls,
            local_exergy,
            local_byproduct,
            local_niche,
            0.38,
            biome,
            adami_rho,
            pop_size,
            env_phase,
        );
        (chimera, gamma)
    }

    pub fn evaluate_trophic_traits(
        &self,
        controls: &EnvironmentControls,
    ) -> (TrophicGuild, f64, f64, f64) {
        let n = self.sequence.len().max(1) as f64;
        let (off_syms, def_syms, syn_syms): (&[u16], &[u16], &[u16]) = match self.kind {
            OrganismKind::RnaDnaProtein => (&[2, 3, 9, 10], &[5, 8, 11], &[1, 4, 7]),
            OrganismKind::HominidPopulation => (&[1, 3, 5], &[6, 9], &[0, 7, 8]),
            OrganismKind::ArtificialIntelligence => (&[3, 4], &[7, 8], &[1, 5]),
        };

        let off_frac = self.sequence.iter().filter(|s| off_syms.contains(s)).count() as f64 / n;
        let def_frac = self.sequence.iter().filter(|s| def_syms.contains(s)).count() as f64 / n;
        let syn_frac = self.sequence.iter().filter(|s| syn_syms.contains(s)).count() as f64 / n;

        if self.has_apicoplast_flag || self.apicoplast_eliminated {
            let ipp_invasion_boost = if self.has_intact_apicoplast_essential_core()
                || (self.apicoplast_eliminated && controls.host_metabolite_rescue >= 0.65)
            {
                0.32
            } else {
                -0.25
            };
            let offense = (off_frac * 1.45 + 0.34 + ipp_invasion_boost).clamp(0.12, 1.0);
            let defense = (def_frac * 1.35 + 0.24).clamp(0.15, 0.95);
            let syntroph_aff = (syn_frac * 0.65).clamp(0.05, 0.65);
            return (TrophicGuild::PredatorParasite, offense, defense, syntroph_aff);
        }

        if !self.endosymbiont_sequence.is_empty()
            || self.egt_transferred_loci >= 2
            || self.endo_tier >= 2
            || self.mro_class != MroClass::None
        {
            let en = self.endosymbiont_sequence.len().max(1) as f64;
            let endo_syn = if self.endosymbiont_sequence.is_empty() {
                syn_frac
            } else {
                self.endosymbiont_sequence
                    .iter()
                    .filter(|s| syn_syms.contains(s))
                    .count() as f64
                    / en
            };
            let egt_f = self.egt_nuclear_fraction();
            let tier_shield = 0.08 * ((self.endo_tier.saturating_sub(1)) as f64);
            let offense = (off_frac * 1.15).clamp(0.05, 0.88);
            let defense = (def_frac * 1.50 + 0.25 * egt_f + tier_shield + 0.22).clamp(0.20, 1.0);
            let syntroph_aff =
                (0.45 * syn_frac + 1.15 * endo_syn + 0.28 * egt_f + 0.20).clamp(0.25, 1.0);
            return (TrophicGuild::Holobiont, offense, defense, syntroph_aff);
        }

        let streamline = ((24.0 - n) / 12.0).clamp(0.0, 1.0);
        let offense = (off_frac * 1.45 + 0.25 * streamline).clamp(0.05, 1.0);
        let defense = (def_frac * 1.65 + 0.15 * (1.0 - streamline)).clamp(0.05, 1.0);
        let syntroph_aff = (syn_frac * 1.55 + 0.30 * streamline).clamp(0.05, 1.0);

        if self.has_active_kleptoplast()
            || (self.klepto_hgt_loci > 0 && self.guild == TrophicGuild::PredatorParasite)
        {
            return (
                TrophicGuild::PredatorParasite,
                offense.max(0.55),
                defense,
                syntroph_aff,
            );
        }

        let allow_pred = controls.trophic_coupling >= 0.05;
        let allow_syn = controls.cross_feeding_rate >= 0.05;

        let guild =
            if allow_pred && offense >= 0.56 && offense > syntroph_aff + 0.04 && def_frac < 0.26 {
                TrophicGuild::PredatorParasite
            } else if allow_syn && syntroph_aff >= 0.52 && self.sequence.len() <= 23 {
                TrophicGuild::Syntroph
            } else {
                TrophicGuild::Producer
            };

        (guild, offense, defense, syntroph_aff)
    }

    pub fn susceptibility_to(&self, predator: &ModeledUnit) -> f64 {
        let k = self.sequence.len().min(predator.sequence.len()).min(6);
        let mut tag_matches = 0usize;
        for i in 0..k {
            if (self.sequence[i] + predator.sequence[i]) % 3 == 0 {
                tag_matches += 1;
            }
        }
        let lock_key = if k == 0 {
            0.35
        } else {
            0.30 + 0.70 * ((tag_matches as f64) / (k as f64))
        };
        let multicell_size_refuge = match self.multicell_state {
            MulticellularState::GermSomaMulticellular => {
                (0.56 - 0.14 * self.soma_fraction + 0.35 * self.somatic_cheater_load).clamp(0.42, 0.95)
            }
            MulticellularState::UndifferentiatedClonal => 0.72,
            MulticellularState::AggregativeChimerism => 0.82,
            MulticellularState::Unicellular => 1.0,
        };
        let soldier_fortress_shield = match self.eusocial_state {
            EusocialState::ObligateSuperorganism => {
                (0.44 - 0.55 * self.soldier_fraction + 0.38 * self.rebel_worker_load).clamp(0.28, 0.90)
            }
            EusocialState::FacultativeEusocial => 0.68,
            EusocialState::SubsocialFamily => 0.84,
            EusocialState::Solitary => 1.0,
        };
        let agro_hygiene_shield = if self.agriculture_state.is_farming() {
            (0.88 - 0.22 * self.defensive_actinobacteria_cover + 0.38 * self.escovopsis_pathogen_load)
                .clamp(0.64, 1.25)
        } else {
            1.0
        };
        let cultural_tool_shield = if self.cultural_state.is_cumulative_ratchet() {
            (0.84 - 0.04 * (self.cultural_ratchet_tier as f64) + 0.32 * self.cultural_drift_erosion)
                .clamp(0.58, 1.15)
        } else {
            1.0
        };
        let compartment_shield = (if self.guild == TrophicGuild::Holobiont {
            (0.78 - 0.08 * self.egt_nuclear_fraction()
                - 0.06 * ((self.endo_tier.saturating_sub(1)) as f64))
                .max(0.54)
        } else {
            1.0
        }) * multicell_size_refuge
            * soldier_fortress_shield
            * agro_hygiene_shield
            * cultural_tool_shield;
        let apical_invasion_mult =
            if predator.has_apicoplast() && predator.has_intact_apicoplast_essential_core() {
                1.28
            } else if predator.apicoplast_eliminated {
                1.22
            } else {
                1.0
            };
        (lock_key
            * compartment_shield
            * apical_invasion_mult
            * (0.45 + 0.85 * predator.offense_trait)
            * (1.15 - 0.72 * self.defense_trait))
            .clamp(0.08, 1.45)
    }

    pub fn build_wiring(
        kind: OrganismKind,
        sequence: &[u16],
        regime: RegimeId,
    ) -> Vec<(usize, usize)> {
        let n = sequence.len();
        if n < 4 || regime == RegimeId::R00Static {
            return Vec::new();
        }
        let mut edges = Vec::new();
        match kind {
            OrganismKind::RnaDnaProtein => {
                for i in 0..(n / 2) {
                    let j = n - 1 - i;
                    if j > i + 2 {
                        let a = sequence[i];
                        let b = sequence[j];
                        let paired = (a == 0 && b == 1)
                            || (a == 1 && b == 0)
                            || (a == 2 && b == 3)
                            || (a == 3 && b == 2)
                            || (a >= 4 && b >= 4 && (a + b) % 3 == 0);
                        if paired {
                            edges.push((i, j));
                        }
                    }
                }
            }
            OrganismKind::HominidPopulation => {
                for i in 0..n.saturating_sub(3) {
                    if (sequence[i] + sequence[i + 3]) % 2 == 1 {
                        edges.push((i, i + 3));
                    }
                }
                if n >= 8 && regime >= RegimeId::R20InteractiveReplicative {
                    edges.push((1, n - 2));
                }
            }
            OrganismKind::ArtificialIntelligence => {
                for i in 1..n.saturating_sub(3) {
                    let sym = sequence[i];
                    if (4..=8).contains(&sym) {
                        let target = (i + 2 + (sym as usize % 3)).min(n - 1);
                        if target > i + 1 {
                            edges.push((i, target));
                        }
                    }
                }
            }
        }
        edges
    }

    pub fn formatted_blueprint(&self, max_tokens: usize) -> String {
        let symbols = self.kind.alphabet_symbols();
        let shown: Vec<&str> = self
            .sequence
            .iter()
            .take(max_tokens)
            .map(|&idx| symbols.get(idx as usize).copied().unwrap_or("?"))
            .collect();
        let sep = match self.kind {
            OrganismKind::RnaDnaProtein => "-",
            OrganismKind::HominidPopulation => " -> ",
            OrganismKind::ArtificialIntelligence => " :: ",
        };
        let mut out = shown.join(sep);
        if self.sequence.len() > max_tokens {
            out.push_str(&format!(" ... (+{})", self.sequence.len() - max_tokens));
        }
        if self.mro_class != MroClass::None && self.mro_class != MroClass::AerobicMitochondrion {
            out.push_str(&format!(
                " || [MRO:{} | Paths:{}/4 (ISC={}) | mtDNA={}]",
                self.mro_class.short_code(),
                self.mro_pathway_count(),
                u8::from((self.mro_pathways_mask & MRO_PATH_ISC) != 0),
                self.mro_mtdna_loci
            ));
        }
        if self.has_apicoplast() {
            out.push_str(&format!(
                " || [Apicoplast 4m | sufB={} | Paths:{}/4 (MEP={},SUF={},FASII={},Heme={})]",
                self.endosymbiont_sequence.len(),
                self.apicoplast_pathway_count(),
                u8::from((self.apicoplast_pathways_mask & APICO_PATH_MEP_IPP) != 0),
                u8::from((self.apicoplast_pathways_mask & APICO_PATH_SUF_FES) != 0),
                u8::from((self.apicoplast_pathways_mask & APICO_PATH_FASII) != 0),
                u8::from((self.apicoplast_pathways_mask & APICO_PATH_HEME) != 0),
            ));
        } else if self.apicoplast_eliminated {
            out.push_str(" || [Cryptosporidium-mode: Apicoplast=lost (IPP-rescued)]");
        } else if self.has_active_kleptoplast()
            || (self.klepto_hgt_loci > 0 && self.guild == TrophicGuild::PredatorParasite)
        {
            out.push_str(&format!(
                " || [Kleptoplast int={:.0}% | HGT:{}/{} | Donors:{}]",
                self.kleptoplast_integrity * 100.0,
                self.klepto_hgt_loci,
                KLEPTO_FIXATION_HGT_THRESHOLD,
                self.mosaic_donor_count()
            ));
        } else if !self.endosymbiont_sequence.is_empty()
            || self.egt_transferred_loci > 0
            || self.endo_tier >= 2
        {
            let endo_shown: Vec<&str> = self
                .endosymbiont_sequence
                .iter()
                .take(3)
                .map(|&idx| symbols.get(idx as usize).copied().unwrap_or("?"))
                .collect();
            let nm_note = if self.endo_tier >= 2 {
                if self.nucleomorph_sequence.is_empty() {
                    format!(" | {}°({}m) NM:lost", self.endo_tier, self.organelle_membranes)
                } else {
                    format!(
                        " | {}°({}m) NM:{}loci",
                        self.endo_tier,
                        self.organelle_membranes,
                        self.nucleomorph_sequence.len()
                    )
                }
            } else {
                String::from(" | 1°(2m)")
            };
            out.push_str(&format!(
                " || [CoRR-Org:{} len={} | Nuc-EGT:{} ({:.0}%){}]",
                endo_shown.join(sep),
                self.endosymbiont_sequence.len(),
                self.egt_transferred_loci,
                self.egt_nuclear_fraction() * 100.0,
                nm_note
            ));
        }
        out
    }

    pub fn heritable_leverage(&self, controls: &EnvironmentControls) -> f64 {
        let inherited_ready = if controls.soft_inheritance_coupling >= 0.10 {
            1.0 + 0.32 * controls.soft_inheritance_coupling * controls.memory_plasticity
        } else {
            1.0
        };
        (self.base_complexity_c * inherited_ready) / self.structure_s.max(2.0)
    }

    pub fn maturation_weight(&self) -> f64 {
        1.0 - (-((self.age as f64 + 1.0) / 3.0)).exp()
    }

    pub fn physical_constraint_fraction(&self, bath_flux: f64) -> f64 {
        if bath_flux.abs() <= 0.03 {
            0.92
        } else {
            let wire_bonus = ((self.wiring.len() as f64) / (self.sequence.len().max(1) as f64))
                .clamp(0.0, 0.35);
            let holo_bonus = if self.has_apicoplast() {
                0.08
            } else if !self.endosymbiont_sequence.is_empty()
                || self.egt_transferred_loci >= 2
                || self.endo_tier >= 2
                || self.mro_class != MroClass::None
            {
                0.06
                    + 0.05 * self.egt_nuclear_fraction()
                    + 0.03 * ((self.endo_tier.saturating_sub(1)) as f64)
            } else if self.has_active_kleptoplast() {
                0.04 * self.kleptoplast_integrity
            } else {
                0.0
            };
            (self.integrity * (0.72 + 0.80 * wire_bonus + holo_bonus)).clamp(0.04, 1.0)
        }
    }

    pub fn step_spatial_movement(
        &mut self,
        grid: &SpatialGrid,
        controls: &EnvironmentControls,
        host_counts: &[u16],
        predator_counts: &[u16],
        rng: &mut Rng,
    ) {
        if controls.bath_exergy_flux <= 0.03 {
            return;
        }
        let pursuit_prob = if self.guild == TrophicGuild::PredatorParasite {
            if self.has_apicoplast() {
                0.92
            } else {
                0.85
            }
        } else if self.guild == TrophicGuild::Holobiont {
            0.82
        } else {
            0.75
        };
        if (controls.memory_plasticity >= 0.45 || self.guild != TrophicGuild::Producer)
            && rng.bernoulli(pursuit_prob)
        {
            let (bx, by) = grid.best_guild_neighbor(
                self.grid_x,
                self.grid_y,
                self.guild,
                host_counts,
                predator_counts,
            );
            self.grid_x = bx;
            self.grid_y = by;
        } else if rng.bernoulli(0.35) {
            let dx = (rng.gen_index(3) as isize) - 1;
            let dy = (rng.gen_index(3) as isize) - 1;
            self.grid_x = (self.grid_x as isize + dx).clamp(0, (grid.width - 1) as isize) as usize;
            self.grid_y = (self.grid_y as isize + dy).clamp(0, (grid.height - 1) as isize) as usize;
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_lifetime(
        &mut self,
        controls: &EnvironmentControls,
        local_exergy: f64,
        local_byproduct: f64,
        local_niche: f64,
        syntroph_relief_bonus: f64,
        biome: BiomeZone,
        adami_rho: f64,
        pop_size: usize,
        env_phase: f64,
        rng: &mut Rng,
    ) {
        self.age = self.age.saturating_add(1);
        let net_margin = self.harvested_work - self.landauer_cost;

        if self.has_active_kleptoplast() {
            let hgt_shield = ((self.klepto_hgt_loci as f64)
                / (KLEPTO_FIXATION_HGT_THRESHOLD as f64))
                .clamp(0.0, 0.85);
            let kp_decay = 0.22 * (1.0 - 0.72 * hgt_shield);
            self.kleptoplast_integrity = (self.kleptoplast_integrity - kp_decay).max(0.0);
            if self.kleptoplast_integrity <= 0.05 {
                self.kleptoplast_sequence.clear();
                self.kleptoplast_integrity = 0.0;
                self.kleptoplast_donor_id = None;
                self.kleptoplast_donor_tier = 0;
            }
        }

        if self.has_apicoplast()
            && !self.has_intact_apicoplast_essential_core()
            && controls.host_metabolite_rescue < 0.65
        {
            self.integrity = (self.integrity * 0.45).max(0.02);
        }

        // CIA `[Fe-S]` essentiality check: losing mitochondrial `ISC` without bacterial cytosolic `SUF`
        // replacement causes immediate nuclear DNA polymerase & Rli1 collapse!
        if !self.has_intact_fes_cluster_biogenesis(controls.cytosolic_suf_hgt_rescue) {
            self.integrity = (self.integrity * 0.42).max(0.02);
        }

        // Multicellular somatic mutation turnover & mitochondrial cytochrome-c apoptotic policing!
        if self.multicell_state.is_multicellular() {
            self.step_multicellular_homeostasis(
                controls.apoptosis_policing_rate,
                controls.single_cell_bottleneck,
            );
        }
        // Eusocial superorganism kin-conflict & mutual worker egg-policing!
        if self.eusocial_state.is_social() {
            self.step_colony_kin_policing(
                controls.queen_monogamy_index,
                controls.worker_policing_rate,
            );
        }
        // Obligate agricultural monoculture *Escovopsis* pressure & *Pseudonocardia* antimicrobial shield!
        if self.agriculture_state.is_farming() {
            self.step_agro_pathogen_and_defensive_microbiome(
                controls.agro_pathogen_pressure,
                controls.defensive_symbiont_potency,
            );
        }
        // Dual genetic-cultural pedagogical maintenance vs. isolation/slippage drift!
        if self.cultural_state.is_cultural() {
            self.step_dual_inheritance_and_cultural_drift(
                controls.symbolic_teaching_fidelity,
                controls.inter_deme_connectivity,
                5,
            );
        }
        // Exosomatic archive replication audit vs. Smaldino-McElreath false-discovery pollution!
        if self.exosomatic_state.is_archival() {
            self.step_archive_verification_and_fdr_dynamics(
                controls.unverified_publish_pressure,
                controls.peer_verification_rigor,
            );
        }
        // Autonomous machine-verified synthesis vs. Shumailov recursive model collapse!
        if self.machine_synthesis_state.has_computational_symbiont() {
            self.step_recursive_synthesis_and_collapse_safeguards(
                controls.synthetic_generation_rate,
                controls.formal_verifier_strictness,
                controls.empirical_grounding_floor,
            );
        }

        if controls.bath_exergy_flux < -0.03
            || (controls.bath_exergy_flux > 0.03 && net_margin < -0.05)
        {
            let deficit = (-controls.bath_exergy_flux).max(-net_margin).max(0.08);
            let decay_rate = (deficit * 0.06).clamp(0.008, 0.13);
            self.integrity = (self.integrity * (1.0 - decay_rate)).max(0.015);
            if !self.is_organelle_enslaved()
                && !self.has_apicoplast()
                && !self.endosymbiont_sequence.is_empty()
                && rng.bernoulli(decay_rate * 0.85)
            {
                self.endosymbiont_sequence.pop();
                if self.endosymbiont_sequence.len() < 4 {
                    self.endosymbiont_sequence.clear();
                    self.endosymbiont_origin_id = None;
                    self.endo_tier = 0;
                    self.organelle_membranes = 0;
                    self.mro_class = MroClass::None;
                    self.mro_pathways_mask = 0;
                    self.mro_mtdna_loci = 0;
                }
            }
            if self.sequence.len() > 4 && rng.bernoulli(decay_rate * 0.65) {
                let drop_idx = if self.kind == OrganismKind::ArtificialIntelligence {
                    1 + rng.gen_index(self.sequence.len().saturating_sub(2).max(1))
                } else {
                    rng.gen_index(self.sequence.len())
                };
                if drop_idx < self.sequence.len()
                    && !(self.kind == OrganismKind::ArtificialIntelligence
                        && (drop_idx == 0 || drop_idx + 1 == self.sequence.len()))
                {
                    self.sequence.remove(drop_idx);
                    if drop_idx < self.lifetime_memory.len() {
                        self.lifetime_memory.remove(drop_idx);
                    }
                }
            }
            for m in &mut self.lifetime_memory {
                *m *= 0.80;
            }
        } else if controls.bath_exergy_flux > 0.03 {
            self.integrity = (self.integrity + 0.08).clamp(0.15, 1.0);
            if controls.egt_rate >= 0.03
                && !self.has_apicoplast()
                && self.endosymbiont_sequence.len() > CORR_CORE_LOCI
                && rng.bernoulli((controls.egt_rate * 0.22).clamp(0.0, 0.35))
            {
                self.migrate_endosymbiont_locus_to_host(false, rng);
            }
            if controls.nucleomorph_reduction_rate >= 0.03
                && !self.nucleomorph_sequence.is_empty()
                && rng.bernoulli((controls.nucleomorph_reduction_rate * 0.28).clamp(0.0, 0.40))
            {
                let allow_sweep = controls.nucleomorph_reduction_rate >= 0.32;
                self.migrate_nucleomorph_locus_to_host(allow_sweep, rng);
            }
            if controls.apicoplast_reduction_rate >= 0.03
                && self.has_apicoplast()
                && rng.bernoulli((controls.apicoplast_reduction_rate * 0.28).clamp(0.0, 0.42))
            {
                self.streamline_apicoplast_pathway(controls.host_metabolite_rescue, false);
            }
            // Hypoxia-driven anaerobic MRO reduction (`AerobicMito -> Hydrogenosome -> Mitosome -> MRO->0`)!
            let loc_o2 = Self::local_oxygen_at(controls.ambient_oxygen_level, biome);
            if controls.anaerobic_mro_reduction_rate >= 0.03
                && loc_o2 < 0.45
                && self.mro_class != MroClass::None
                && self.mro_class != MroClass::Amiteukaryote
                && !self.multicell_state.is_multicellular()
                && rng.bernoulli(
                    (controls.anaerobic_mro_reduction_rate * (0.55 - loc_o2)).clamp(0.0, 0.38),
                )
            {
                self.step_anaerobic_mro_reduction(controls.cytosolic_suf_hgt_rescue, false);
            }
            // Normoxic chimeric Holobionts evolve clonal multicellularity & germ-soma division of labor!
            if self.guild == TrophicGuild::Holobiont && loc_o2 >= 0.35 {
                if self.multicell_state == MulticellularState::Unicellular
                    && controls.multicellular_adhesion_rate >= 0.04
                    && rng.bernoulli((controls.multicellular_adhesion_rate * 0.22).clamp(0.0, 0.35))
                {
                    self.step_clonal_multicellularity(controls.single_cell_bottleneck, rng);
                } else if matches!(
                    self.multicell_state,
                    MulticellularState::UndifferentiatedClonal
                        | MulticellularState::AggregativeChimerism
                ) && controls.germ_soma_differentiation_rate >= 0.04
                    && rng.bernoulli(
                        (controls.germ_soma_differentiation_rate * 0.25).clamp(0.0, 0.38),
                    )
                {
                    self.step_clonal_multicellularity(controls.single_cell_bottleneck, rng);
                } else if self.multicell_state == MulticellularState::GermSomaMulticellular
                    && self.eusocial_state != EusocialState::ObligateSuperorganism
                    && controls.eusocial_caste_canalization_rate >= 0.04
                    && rng.bernoulli(
                        (controls.eusocial_caste_canalization_rate * 0.26).clamp(0.0, 0.40),
                    )
                {
                    self.step_eusocial_superorganism(controls.queen_monogamy_index, rng);
                } else if self.eusocial_state == EusocialState::ObligateSuperorganism
                    && self.agriculture_state != AgricultureState::LeafCutterAgroEcosystem
                    && controls.crop_domestication_rate >= 0.04
                    && rng.bernoulli(
                        (controls.crop_domestication_rate * 0.28).clamp(0.0, 0.44),
                    )
                {
                    self.step_agricultural_domestication(
                        controls.defensive_symbiont_potency,
                        rng,
                    );
                } else if self.cultural_state != CulturalState::SymbolicCumulativeCulture
                    && controls.cultural_innovation_rate >= 0.04
                    && rng.bernoulli(
                        (controls.cultural_innovation_rate * 0.28).clamp(0.0, 0.45),
                    )
                {
                    self.step_cumulative_culture_transition(
                        controls.symbolic_teaching_fidelity,
                        controls.inter_deme_connectivity,
                        rng,
                    );
                } else if self.cultural_state == CulturalState::SymbolicCumulativeCulture
                    && self.exosomatic_state != ExosomaticState::PeerVerifiedInstitution
                    && controls.archival_externalization_rate >= 0.04
                    && rng.bernoulli(
                        (controls.archival_externalization_rate * 0.30).clamp(0.0, 0.48),
                    )
                {
                    self.step_exosomatic_archive_transition(
                        controls.peer_verification_rigor,
                        rng,
                    );
                } else if self.exosomatic_state.has_written_archive()
                    && self.machine_synthesis_state != MachineSynthesisState::MachineVerifiedSymbiosis
                    && controls.synthetic_generation_rate >= 0.04
                    && rng.bernoulli(
                        (controls.synthetic_generation_rate * 0.30).clamp(0.0, 0.48),
                    )
                {
                    self.step_machine_synthesis_transition(
                        controls.formal_verifier_strictness,
                        controls.empirical_grounding_floor,
                        rng,
                    );
                }
            }
            if controls.memory_plasticity >= 0.10 {
                let resource_signal =
                    (local_exergy + 0.55 * local_byproduct + 0.35 * local_niche).min(1.5);
                let learn_gain = controls.memory_plasticity * resource_signal * 0.16;
                for (i, mem) in self.lifetime_memory.iter_mut().enumerate() {
                    let target = 0.72 + 0.28 * (((i as f64) * 0.7 + env_phase).sin().abs());
                    let noise = rng.range_f64(-0.015, 0.025);
                    *mem = (*mem + learn_gain * (target - *mem) + noise).clamp(0.0, 1.40);
                }
            } else {
                for mem in &mut self.lifetime_memory {
                    *mem *= 0.86;
                }
            }
        }

        let approx_regime = if controls.bath_exergy_flux.abs() <= 0.03 {
            RegimeId::R00Static
        } else {
            RegimeId::R20InteractiveReplicative
        };
        self.wiring = Self::build_wiring(self.kind, &self.sequence, approx_regime);
        self.net_trophic_flux = 0.0;
        self.recompute_phenome(
            controls,
            local_exergy,
            local_byproduct,
            local_niche,
            syntroph_relief_bonus,
            biome,
            adami_rho,
            pop_size,
            env_phase,
        );
    }

    pub fn apply_trophic_encounter_delta(
        &mut self,
        delta_work: f64,
        controls: &EnvironmentControls,
    ) {
        let scaled_delta = if delta_work > 0.0 {
            if self.has_apicoplast() && self.has_intact_apicoplast_essential_core() {
                delta_work * 1.24
            } else if self.apicoplast_eliminated && controls.host_metabolite_rescue >= 0.65 {
                delta_work * 1.22
            } else if self.has_apicoplast() && !self.has_intact_apicoplast_essential_core() {
                delta_work * 0.25
            } else {
                delta_work
            }
        } else {
            delta_work
        };

        self.net_trophic_flux += scaled_delta;
        self.harvested_work =
            (self.primary_work + self.syntrophic_work + self.net_trophic_flux).max(0.02);

        let net_exergy_margin = self.harvested_work - self.landauer_cost;
        let exergy_factor = if controls.bath_exergy_flux.abs() <= 0.03 {
            0.30
        } else if net_exergy_margin < -0.02 {
            (0.28 * self.integrity).clamp(0.008, 0.32)
        } else {
            (0.76 + 0.48 * net_exergy_margin.min(1.8)).clamp(0.65, 1.75)
        };

        let mean_mem = if self.lifetime_memory.is_empty() {
            0.0
        } else {
            self.lifetime_memory.iter().sum::<f64>() / (self.lifetime_memory.len() as f64)
        };
        let plastic_amplification = 1.0 + 0.72 * mean_mem * controls.memory_plasticity.max(0.0);
        self.complexity_c = (self.base_complexity_c * (exergy_factor / 1.15) * plastic_amplification)
            .max(0.08);
        let leverage = self.complexity_c / self.structure_s.max(2.0);
        let energy_viability = (0.55 + 0.68 * net_exergy_margin).clamp(0.04, 2.05);
        self.fitness = ((0.48 * (leverage / 2.3) + 0.52 * energy_viability) * self.integrity)
            .clamp(0.02, 4.5);
    }

    #[allow(clippy::too_many_arguments)]
    pub fn replicate_with_partner(
        &self,
        donor: Option<&ModeledUnit>,
        child_id: u64,
        child_x: usize,
        child_y: usize,
        child_biome: BiomeZone,
        local_exergy: f64,
        local_byproduct: f64,
        local_niche: f64,
        controls: &EnvironmentControls,
        adami_rho: f64,
        pop_size: usize,
        env_phase: f64,
        rng: &mut Rng,
    ) -> (ModeledUnit, f64, bool, bool) {
        let sigma = self.kind.alphabet_symbols().len() as u16;
        let mut child_seq = self.sequence.clone();
        let mut child_endo = self.endosymbiont_sequence.clone();
        let mut child_nm = self.nucleomorph_sequence.clone();
        let child_egt_loci = self.egt_transferred_loci;
        let child_nm_loci_xfers = self.nucleomorph_transferred_loci;
        let child_nm_eliminated = self.nucleomorph_eliminated;
        let child_klepto_hgt = self.klepto_hgt_loci;
        let child_donor_mask = self.klepto_donor_mask;
        let hgt_div_shield =
            ((self.klepto_hgt_loci as f64) / (KLEPTO_FIXATION_HGT_THRESHOLD as f64)).clamp(0.0, 0.85);
        let diluted_kp_int = if self.has_active_kleptoplast() {
            self.kleptoplast_integrity * (0.48 + 0.14 * hgt_div_shield)
        } else {
            0.0
        };
        let (child_kp_seq, child_kp_int, child_kp_donor, child_kp_tier) = if diluted_kp_int > 0.10 {
            (
                self.kleptoplast_sequence.clone(),
                diluted_kp_int,
                self.kleptoplast_donor_id,
                self.kleptoplast_donor_tier,
            )
        } else {
            (Vec::new(), 0.0, None, 0)
        };

        // Nuclear DNA polymerases depend on ISC -> CIA [Fe-S] cofactors, AND germline sequestration
        // shields template copying from OXPHOS `ROS` dirty-work mutations (Radzvilavicius et al., 2016)!
        let fes_ok = self.has_intact_fes_cluster_biogenesis(controls.cytosolic_suf_hgt_rescue);
        let f_base = if fes_ok {
            controls.effective_fidelity()
        } else {
            (controls.effective_fidelity() * 0.58).clamp(0.25, 0.68)
        };
        let queen_germline_shield = if self.eusocial_state == EusocialState::ObligateSuperorganism {
            0.35 * (1.0 - self.rebel_worker_load).clamp(0.0, 1.0)
        } else if self.eusocial_state == EusocialState::FacultativeEusocial {
            0.15
        } else {
            0.0
        };
        let f_eff = match self.multicell_state {
            MulticellularState::GermSomaMulticellular => {
                // Metabolically quiescent germline + sheltered eusocial queen isolation!
                let gs_f = 1.0 - (1.0 - f_base) * (1.0 - 0.62 * self.soma_fraction);
                (1.0 - (1.0 - gs_f) * (1.0 - queen_germline_shield)).clamp(0.25, 0.9997)
            }
            MulticellularState::UndifferentiatedClonal | MulticellularState::AggregativeChimerism => {
                // Totipotent dividing cells perform active OXPHOS dirty work simultaneously -> higher ROS load!
                (f_base - 0.026).clamp(0.25, 0.995)
            }
            MulticellularState::Unicellular => f_base,
        };

        let start_idx = usize::from(self.kind == OrganismKind::ArtificialIntelligence);
        let end_idx = if self.kind == OrganismKind::ArtificialIntelligence {
            child_seq.len().saturating_sub(1)
        } else {
            child_seq.len()
        };

        let mut recombined = false;
        if let Some(partner) = donor {
            if controls.recombination_rate >= 0.05
                && partner.id != self.id
                && partner.sequence.len() >= 8
                && rng.bernoulli(controls.recombination_rate * 0.55)
            {
                let p_start = usize::from(self.kind == OrganismKind::ArtificialIntelligence);
                let p_end = if self.kind == OrganismKind::ArtificialIntelligence {
                    partner.sequence.len().saturating_sub(1)
                } else {
                    partner.sequence.len()
                };
                if p_end > p_start + 3 && end_idx > start_idx + 3 {
                    let mod_len = 3.min(p_end - p_start).min(end_idx - start_idx);
                    let src_pos = p_start + rng.gen_index(p_end - p_start - mod_len + 1);
                    let dst_pos = start_idx + rng.gen_index(end_idx - start_idx - mod_len + 1);
                    child_seq[dst_pos..(dst_pos + mod_len)]
                        .copy_from_slice(&partner.sequence[src_pos..(src_pos + mod_len)]);
                    recombined = true;
                }
            }
        }

        let mutable_span =
            (end_idx.saturating_sub(start_idx) + child_endo.len() + child_nm.len()).max(1);
        let mut unchanged_loci = 0usize;
        let mut assimilated = false;
        let biome_targets = self.kind.biome_target_symbols(child_biome, env_phase);

        if controls.heritable_variation >= 0.01 || f_eff < 0.99 || recombined {
            let infid = (1.0 - f_eff).clamp(0.0, 0.65);
            let per_locus_mut =
                (infid * 0.65 + controls.heritable_variation * 0.20).clamp(0.0, 0.65);

            for idx in start_idx..end_idx {
                let old_sym = self.sequence[idx];
                if rng.bernoulli(per_locus_mut) {
                    child_seq[idx] = if f_eff >= 0.88 && rng.bernoulli(0.35) {
                        biome_targets[rng.gen_index(biome_targets.len())]
                    } else if self.kind == OrganismKind::ArtificialIntelligence {
                        1 + (rng.gen_index((sigma - 2).max(1) as usize) as u16)
                    } else {
                        rng.gen_index(sigma as usize) as u16
                    };
                }
                if child_seq[idx] == old_sym {
                    unchanged_loci += 1;
                }
            }

            let syn_motifs: &[u16] = match self.kind {
                OrganismKind::RnaDnaProtein => &[1, 4, 7, 8],
                OrganismKind::HominidPopulation => &[0, 6, 7, 8],
                OrganismKind::ArtificialIntelligence => &[1, 5, 7, 8],
            };
            let org_mut = (per_locus_mut * 1.25).clamp(0.0, 0.75);
            for (e_idx, slot) in child_endo.iter_mut().enumerate() {
                let old_e = self.endosymbiont_sequence[e_idx];
                if rng.bernoulli(org_mut) {
                    *slot = if f_eff >= 0.88 && rng.bernoulli(0.38) {
                        syn_motifs[rng.gen_index(syn_motifs.len())]
                    } else {
                        rng.gen_index(sigma as usize) as u16
                    };
                }
                if *slot == old_e {
                    unchanged_loci += 1;
                }
            }

            let nm_mut = (per_locus_mut * 1.32).clamp(0.0, 0.78);
            for (n_idx, slot) in child_nm.iter_mut().enumerate() {
                let old_n = self.nucleomorph_sequence[n_idx];
                if rng.bernoulli(nm_mut) {
                    *slot = rng.gen_index(sigma as usize) as u16;
                }
                if *slot == old_n {
                    unchanged_loci += 1;
                }
            }

            let mean_parent_mem = if self.lifetime_memory.is_empty() {
                0.0
            } else {
                self.lifetime_memory.iter().sum::<f64>() / (self.lifetime_memory.len() as f64)
            };
            if controls.soft_inheritance_coupling >= 0.08 && mean_parent_mem > 0.12 {
                let assimilate_prob =
                    (controls.soft_inheritance_coupling * mean_parent_mem * 0.55).clamp(0.0, 0.65);
                if rng.bernoulli(assimilate_prob) && end_idx > start_idx {
                    let locus = start_idx + rng.gen_index(end_idx - start_idx);
                    child_seq[locus] = biome_targets[rng.gen_index(biome_targets.len())];
                    assimilated = true;
                }
            }

            let del_bias = if controls.heritable_variation < 0.01 && f_eff >= 0.99 {
                0.0
            } else if controls.cross_feeding_rate >= 0.10 || controls.trophic_coupling >= 0.10 {
                0.09
            } else {
                0.03
            };
            if child_seq.len() < 44 && rng.bernoulli(controls.heritable_variation * 0.16) {
                let insert_pos = if self.kind == OrganismKind::ArtificialIntelligence {
                    end_idx.max(1)
                } else {
                    rng.gen_index(child_seq.len().max(1))
                };
                let new_sym = if controls.soft_inheritance_coupling >= 0.10 {
                    biome_targets[rng.gen_index(biome_targets.len())]
                } else if self.kind == OrganismKind::ArtificialIntelligence {
                    1 + (rng.gen_index((sigma - 2).max(1) as usize) as u16)
                } else {
                    rng.gen_index(sigma as usize) as u16
                };
                child_seq.insert(insert_pos.min(child_seq.len()), new_sym);
            } else if child_endo.is_empty()
                && child_egt_loci == 0
                && self.endo_tier == 0
                && child_klepto_hgt == 0
                && !self.has_apicoplast_flag
                && self.mro_class == MroClass::None
                && child_seq.len() > 11
                && rng.bernoulli(infid * 0.28 + del_bias)
            {
                let del_pos = if self.kind == OrganismKind::ArtificialIntelligence {
                    start_idx + rng.gen_index((end_idx - start_idx).max(1))
                } else {
                    rng.gen_index(child_seq.len())
                };
                if del_pos < end_idx && child_seq.len() > 9 {
                    child_seq.remove(del_pos);
                }
            }
        } else {
            unchanged_loci = mutable_span;
        }

        let f_obs = ((unchanged_loci as f64) / (mutable_span as f64)).clamp(0.0, 1.0);
        let wiring = Self::build_wiring(self.kind, &child_seq, RegimeId::R30FullyDarwinian);

        let inherited_mem_scale = if controls.soft_inheritance_coupling >= 0.10 {
            (controls.soft_inheritance_coupling * 0.55 + 0.08 * local_niche.min(1.5))
                .clamp(0.0, 0.65)
        } else {
            0.0
        };
        let mut lifetime_memory = vec![0.0; child_seq.len()];
        if inherited_mem_scale > 0.0 {
            assimilated = true;
            for (i, slot) in lifetime_memory.iter_mut().enumerate() {
                let parent_m = self.lifetime_memory.get(i).copied().unwrap_or(0.35);
                *slot = parent_m * inherited_mem_scale;
            }
        }

        let mut child = ModeledUnit {
            id: child_id,
            kind: self.kind,
            generation: self.generation.saturating_add(1),
            age: 0,
            grid_x: child_x,
            grid_y: child_y,
            species_id: self.species_id,
            guild: self.guild,
            endo_tier: self.endo_tier,
            organelle_membranes: self.organelle_membranes,
            sequence: child_seq,
            kleptoplast_sequence: child_kp_seq,
            kleptoplast_integrity: child_kp_int,
            kleptoplast_donor_id: child_kp_donor,
            kleptoplast_donor_tier: child_kp_tier,
            klepto_hgt_loci: child_klepto_hgt,
            klepto_donor_mask: child_donor_mask,
            has_apicoplast_flag: self.has_apicoplast_flag,
            apicoplast_pathways_mask: self.apicoplast_pathways_mask,
            apicoplast_eliminated: self.apicoplast_eliminated,
            mro_class: self.mro_class,
            mro_pathways_mask: self.mro_pathways_mask,
            mro_mtdna_loci: self.mro_mtdna_loci,
            multicell_state: self.multicell_state,
            cell_count: self.cell_count.max(1),
            soma_fraction: self.soma_fraction,
            kin_relatedness_r: if !self.multicell_state.is_multicellular()
                || controls.single_cell_bottleneck
            {
                1.0
            } else {
                (self.kin_relatedness_r * 0.82).clamp(0.28, 0.72)
            },
            somatic_cheater_load: if !self.multicell_state.is_multicellular()
                || controls.single_cell_bottleneck
            {
                0.0
            } else {
                (self.somatic_cheater_load * 0.85 + 0.12).clamp(0.0, 1.0)
            },
            cumulative_apoptotic_culls: 0,
            eusocial_state: self.eusocial_state,
            colony_size: self.colony_size.max(1),
            sterile_worker_fraction: self.sterile_worker_fraction,
            soldier_fraction: self.soldier_fraction,
            colony_relatedness_r: (0.25 + 0.25 * controls.queen_monogamy_index.clamp(0.0, 1.0)).clamp(0.25, 0.50),
            rebel_worker_load: if self.eusocial_state == EusocialState::ObligateSuperorganism
                && controls.worker_policing_rate >= 0.45
            {
                0.0
            } else {
                (self.rebel_worker_load * 0.55).clamp(0.0, 1.0)
            },
            cumulative_policed_eggs: 0,
            agriculture_state: self.agriculture_state,
            cultivar_polyploidy: self.cultivar_polyploidy.max(1),
            crop_clonal_uniformity: self.crop_clonal_uniformity,
            gongylidia_yield: self.gongylidia_yield,
            escovopsis_pathogen_load: if self.agriculture_state.is_farming()
                && controls.defensive_symbiont_potency >= 0.45
            {
                (self.escovopsis_pathogen_load * 0.20).clamp(0.0, 1.0)
            } else {
                (self.escovopsis_pathogen_load * 0.65).clamp(0.0, 1.0)
            },
            defensive_actinobacteria_cover: if self.agriculture_state.is_farming() {
                controls.defensive_symbiont_potency.clamp(0.0, 1.0)
            } else {
                0.0
            },
            cumulative_antimicrobial_shields: 0,
            cultural_state: self.cultural_state,
            deme_id: if rng.bernoulli((controls.inter_deme_connectivity * 0.22).clamp(0.0, 0.25)) {
                (self.deme_id + 1) % 4
            } else {
                self.deme_id
            },
            cultural_sequence: {
                let mut c_seq = self.cultural_sequence.clone();
                if controls.symbolic_teaching_fidelity < 0.78 && c_seq.len() > 4 {
                    c_seq.pop();
                } else if let Some(p) = donor {
                    if controls.symbolic_teaching_fidelity >= 0.84
                        && controls.inter_deme_connectivity >= 0.12
                        && controls.inter_deme_connectivity <= 0.78
                        && c_seq.len() < 22
                    {
                        for &sym in &p.cultural_sequence {
                            if !c_seq.contains(&sym) && c_seq.len() < 22 {
                                c_seq.push(sym);
                                break;
                            }
                        }
                    }
                }
                c_seq
            },
            cultural_ratchet_tier: self.cultural_ratchet_tier,
            neural_learning_capacity: self.neural_learning_capacity,
            cultural_drift_erosion: if controls.symbolic_teaching_fidelity >= 0.84 {
                (self.cultural_drift_erosion * 0.35).clamp(0.0, 1.0)
            } else {
                (self.cultural_drift_erosion * 0.85).clamp(0.0, 1.0)
            },
            cumulative_cross_deme_recombinations: 0,
            exosomatic_state: self.exosomatic_state,
            archive_sequence: {
                let mut a_seq = self.archive_sequence.clone();
                if let Some(p) = donor {
                    if self.exosomatic_state.has_written_archive()
                        && controls.peer_verification_rigor >= 0.60
                        && a_seq.len() < 48
                    {
                        for &sym in &p.archive_sequence {
                            if !a_seq.contains(&sym) && a_seq.len() < 48 {
                                a_seq.push(sym);
                                break;
                            }
                        }
                    }
                }
                a_seq
            },
            cognitive_offloading_fraction: self.cognitive_offloading_fraction,
            archive_false_discovery_load: if controls.peer_verification_rigor >= 0.65 {
                (self.archive_false_discovery_load * 0.30).clamp(0.0, 1.0)
            } else {
                (self.archive_false_discovery_load * 0.85).clamp(0.0, 1.0)
            },
            cumulative_peer_verifications: 0,
            machine_synthesis_state: self.machine_synthesis_state,
            computational_symbiont_sequence: self.computational_symbiont_sequence.clone(),
            distribution_tail_diversity: self.distribution_tail_diversity,
            recursive_autophagy_load: if controls.formal_verifier_strictness >= 0.72
                && controls.empirical_grounding_floor >= 0.12
            {
                (self.recursive_autophagy_load * 0.25).clamp(0.0, 1.0)
            } else {
                (self.recursive_autophagy_load * 0.90).clamp(0.0, 1.0)
            },
            cumulative_formal_proofs: 0,
            planetary_epoch: self.planetary_epoch,
            guild_o2_flux: self.guild_o2_flux,
            guild_carbon_drawdown: self.guild_carbon_drawdown,
            guild_nutrient_recycling: self.guild_nutrient_recycling,
            cumulative_stewardship_actions: 0,
            abiogenesis_stage: self.abiogenesis_stage,
            rna_oligomer_length_nt: self.rna_oligomer_length_nt,
            raf_closure_score: self.raf_closure_score,
            vesicle_membrane_stability: self.vesicle_membrane_stability,
            progenote_hgt_flux: self.progenote_hgt_flux,
            chirality_code_stage: self.chirality_code_stage,
            ribose_ee: self.ribose_ee,
            peptide_l_ee: self.peptide_l_ee,
            genetic_code_mp: self.genetic_code_mp,
            amino_acid_alphabet_size: self.amino_acid_alphabet_size,
            genome_substrate_stage: self.genome_substrate_stage,
            chromosome_capacity_kb: self.chromosome_capacity_kb,
            deamination_uracil_load: self.deamination_uracil_load,
            antiviral_defense_titer: self.antiviral_defense_titer,
            nucleomorph_sequence: child_nm,
            nucleomorph_transferred_loci: child_nm_loci_xfers,
            nucleomorph_eliminated: child_nm_eliminated,
            endosymbiont_sequence: child_endo,
            endosymbiont_origin_id: self.endosymbiont_origin_id,
            egt_transferred_loci: child_egt_loci,
            wiring,
            lifetime_memory,
            integrity: 1.0,
            structure_raw_s: 8.0,
            structure_s: 8.0,
            base_complexity_c: 8.0,
            complexity_c: 8.0,
            primary_work: 0.0,
            syntrophic_work: 0.0,
            net_trophic_flux: 0.0,
            harvested_work: 0.0,
            landauer_cost: 0.0,
            byproduct_output: 0.0,
            niche_output: 0.0,
            offense_trait: 0.0,
            defense_trait: 0.0,
            fitness: 0.5,
        };

        if controls.egt_rate >= 0.03
            && !child.has_apicoplast()
            && child.endosymbiont_sequence.len() > CORR_CORE_LOCI
            && rng.bernoulli((controls.egt_rate * 0.55).clamp(0.0, 0.70))
        {
            child.migrate_endosymbiont_locus_to_host(false, rng);
        }

        if controls.nucleomorph_reduction_rate >= 0.03
            && !child.nucleomorph_sequence.is_empty()
            && rng.bernoulli((controls.nucleomorph_reduction_rate * 0.58).clamp(0.0, 0.75))
        {
            let allow_sweep = controls.nucleomorph_reduction_rate >= 0.30;
            child.migrate_nucleomorph_locus_to_host(allow_sweep, rng);
        }

        if controls.apicoplast_reduction_rate >= 0.04 {
            if child.guild == TrophicGuild::Holobiont
                && rng.bernoulli((controls.apicoplast_reduction_rate * 0.24).clamp(0.0, 0.36))
            {
                child.convert_holobiont_to_apicoplast_parasite(rng);
            } else if child.has_apicoplast()
                && rng.bernoulli((controls.apicoplast_reduction_rate * 0.45).clamp(0.0, 0.65))
            {
                child.streamline_apicoplast_pathway(controls.host_metabolite_rescue, false);
            }
        }

        // Germline anaerobic MRO reduction in hypoxic patches (`AerobicMito -> Hydrogenosome -> Mitosome -> MRO->0`)!
        let c_loc_o2 = Self::local_oxygen_at(controls.ambient_oxygen_level, child_biome);
        if controls.anaerobic_mro_reduction_rate >= 0.04
            && c_loc_o2 < 0.45
            && child.mro_class != MroClass::None
            && child.mro_class != MroClass::Amiteukaryote
            && !child.multicell_state.is_multicellular()
            && rng.bernoulli((controls.anaerobic_mro_reduction_rate * 0.38).clamp(0.0, 0.55))
        {
            child.step_anaerobic_mro_reduction(controls.cytosolic_suf_hgt_rescue, false);
        }

        // Developmental transition to clonal multicellularity & germ-soma differentiation!
        if child.guild == TrophicGuild::Holobiont && c_loc_o2 >= 0.35 {
            if child.multicell_state == MulticellularState::Unicellular
                && controls.multicellular_adhesion_rate >= 0.04
                && rng.bernoulli((controls.multicellular_adhesion_rate * 0.28).clamp(0.0, 0.45))
            {
                child.step_clonal_multicellularity(controls.single_cell_bottleneck, rng);
            } else if matches!(
                child.multicell_state,
                MulticellularState::UndifferentiatedClonal
                    | MulticellularState::AggregativeChimerism
            ) && controls.germ_soma_differentiation_rate >= 0.04
                && rng.bernoulli(
                    (controls.germ_soma_differentiation_rate * 0.32).clamp(0.0, 0.50),
                )
            {
                child.step_clonal_multicellularity(controls.single_cell_bottleneck, rng);
            } else if child.multicell_state == MulticellularState::GermSomaMulticellular
                && child.eusocial_state != EusocialState::ObligateSuperorganism
                && controls.eusocial_caste_canalization_rate >= 0.04
                && rng.bernoulli(
                    (controls.eusocial_caste_canalization_rate * 0.30).clamp(0.0, 0.48),
                )
            {
                child.step_eusocial_superorganism(controls.queen_monogamy_index, rng);
            } else if child.eusocial_state == EusocialState::ObligateSuperorganism
                && child.agriculture_state != AgricultureState::LeafCutterAgroEcosystem
                && controls.crop_domestication_rate >= 0.04
                && rng.bernoulli(
                    (controls.crop_domestication_rate * 0.32).clamp(0.0, 0.50),
                )
            {
                child.step_agricultural_domestication(
                    controls.defensive_symbiont_potency,
                    rng,
                );
            } else if child.cultural_state != CulturalState::SymbolicCumulativeCulture
                && controls.cultural_innovation_rate >= 0.04
                && rng.bernoulli(
                    (controls.cultural_innovation_rate * 0.32).clamp(0.0, 0.50),
                )
            {
                child.step_cumulative_culture_transition(
                    controls.symbolic_teaching_fidelity,
                    controls.inter_deme_connectivity,
                    rng,
                );
            } else if child.cultural_state == CulturalState::SymbolicCumulativeCulture
                && child.exosomatic_state != ExosomaticState::PeerVerifiedInstitution
                && controls.archival_externalization_rate >= 0.04
                && rng.bernoulli(
                    (controls.archival_externalization_rate * 0.34).clamp(0.0, 0.52),
                )
            {
                child.step_exosomatic_archive_transition(
                    controls.peer_verification_rigor,
                    rng,
                );
            } else if child.exosomatic_state.has_written_archive()
                && child.machine_synthesis_state != MachineSynthesisState::MachineVerifiedSymbiosis
                && controls.synthetic_generation_rate >= 0.04
                && rng.bernoulli(
                    (controls.synthetic_generation_rate * 0.34).clamp(0.0, 0.52),
                )
            {
                child.step_machine_synthesis_transition(
                    controls.formal_verifier_strictness,
                    controls.empirical_grounding_floor,
                    rng,
                );
            }
        }

        child.recompute_phenome(
            controls,
            local_exergy,
            local_byproduct,
            local_niche,
            0.0,
            child_biome,
            adami_rho,
            pop_size,
            env_phase,
        );
        (child, f_obs, assimilated, recombined)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn recompute_phenome(
        &mut self,
        controls: &EnvironmentControls,
        local_exergy: f64,
        local_byproduct: f64,
        local_niche: f64,
        syntroph_relief_bonus: f64,
        biome: BiomeZone,
        adami_rho: f64,
        replicating_pop_size: usize,
        env_phase: f64,
    ) {
        let sigma_len = self.kind.alphabet_symbols().len();
        let host_raw_s = compute_raw_structure_bits(&self.sequence, sigma_len, self.wiring.len());
        let endo_raw_s = if self.endosymbiont_sequence.is_empty() {
            0.0
        } else if self.has_apicoplast() {
            compute_raw_structure_bits(&self.endosymbiont_sequence, sigma_len, 1) * 0.55
        } else {
            compute_raw_structure_bits(&self.endosymbiont_sequence, sigma_len, 2) * 1.18
        };
        let nm_raw_s = if self.nucleomorph_sequence.is_empty() {
            0.0
        } else {
            compute_raw_structure_bits(&self.nucleomorph_sequence, sigma_len, 2) * 1.42
        };
        let kp_somatic_s = if self.has_active_kleptoplast() {
            compute_raw_structure_bits(&self.kleptoplast_sequence, sigma_len, 1) * 0.28
        } else {
            0.0
        };
        // Landauer maintenance bits for MRO `mtDNA` + `TOM/TIM` protein-import machinery:
        // shedding `OXPHOS` + `mtDNA` into a `Hydrogenosome` (3 pathways, mtDNA=0), `Mitosome` (1 ISC pathway, mtDNA=0),
        // or `Amiteukaryote` (0 pathways, MRO->0) progressively drops organelle Landauer heat `Q`!
        let mro_maintenance_bits = match self.mro_class {
            MroClass::None | MroClass::Amiteukaryote => 0.0,
            MroClass::AerobicMitochondrion => {
                (self.mro_mtdna_loci as f64) * 1.15 + (self.mro_pathway_count() as f64) * 0.45
            }
            MroClass::AnaerobicMitochondrion => {
                (self.mro_mtdna_loci as f64) * 0.95 + (self.mro_pathway_count() as f64) * 0.40
            }
            MroClass::Hydrogenosome => (self.mro_pathway_count() as f64) * 0.36,
            MroClass::Mitosome => (self.mro_pathway_count() as f64) * 0.28,
        };
        let per_nm_import_cost = if self.endo_tier >= 3 { 0.64 } else { 0.52 };
        let apico_pathway_import_bits = if self.has_apicoplast() {
            (self.apicoplast_pathway_count() as f64) * 0.42
        } else {
            0.0
        };
        let eusocial_comm_bits = match self.eusocial_state {
            EusocialState::Solitary => 0.0,
            EusocialState::SubsocialFamily => 0.65,
            EusocialState::FacultativeEusocial => 1.25,
            EusocialState::ObligateSuperorganism => {
                1.85 + 0.90 * self.total_sterile_caste_fraction()
            }
        };
        // Cuticular exocrine crypt maintenance partially offset by reciprocal host metabolic streamlining!
        let agro_crypt_net_bits = match self.agriculture_state {
            AgricultureState::HunterForager => 0.0,
            AgricultureState::LowerFacultativeFarming => 0.35,
            AgricultureState::HigherDomesticatedCrop => {
                0.42 + 0.25 * self.defensive_actinobacteria_cover - 0.30
            }
            AgricultureState::LeafCutterAgroEcosystem => {
                0.55 + 0.30 * self.defensive_actinobacteria_cover - 0.38
            }
        };
        // Cognitive offloading onto non-metabolic written/digital archives slashes per-individual
        // Landauer neural storage overhead (DeSilva et al., 2021; Grinschgl et al., 2021)!
        let offloaded_neural_load = self.neural_learning_capacity
            * (1.0 - 0.58 * self.cognitive_offloading_fraction).clamp(0.42, 1.0);
        let archival_curation_bits = match self.exosomatic_state {
            ExosomaticState::EndosomaticBound => 0.0,
            ExosomaticState::MaterialNotation => 0.06,
            ExosomaticState::WrittenArchive => 0.12,
            ExosomaticState::PeerVerifiedInstitution => 0.16,
        };
        // FunSearch distributed CPU evaluators keep Landauer verification cost low (Romera-Paredes et al., 2024)!
        let verifier_compute_bits = match self.machine_synthesis_state {
            MachineSynthesisState::ManualArchival => 0.0,
            MachineSynthesisState::AlgorithmicMediation => 0.08,
            MachineSynthesisState::UnverifiedRecursiveGen => 0.24,
            MachineSynthesisState::MachineVerifiedSymbiosis => 0.18,
        };
        let neural_brain_bits =
            1.45 * offloaded_neural_load + archival_curation_bits + verifier_compute_bits;
        let multicell_coord_bits = (match self.multicell_state {
            MulticellularState::Unicellular => 0.0,
            MulticellularState::AggregativeChimerism => 1.45,
            MulticellularState::UndifferentiatedClonal => 1.15,
            MulticellularState::GermSomaMulticellular => 1.65 + 0.85 * self.soma_fraction,
        }) + eusocial_comm_bits
            + agro_crypt_net_bits
            + neural_brain_bits;
        let import_overhead_bits = (if self.apicoplast_eliminated {
            0.0
        } else {
            (self.egt_transferred_loci as f64) * 0.35
                + (self.nucleomorph_transferred_loci as f64) * per_nm_import_cost
                + (self.klepto_hgt_loci as f64) * 0.28
                + apico_pathway_import_bits
        }) + mro_maintenance_bits
            + multicell_coord_bits;

        let encephalization_s = if self.cultural_state.is_cultural()
            || self.exosomatic_state.is_archival()
            || self.machine_synthesis_state.has_computational_symbiont()
        {
            neural_brain_bits
        } else {
            0.0
        };
        self.structure_raw_s =
            host_raw_s + endo_raw_s * 0.65 + nm_raw_s * 0.55 + encephalization_s;

        let rho_phys = self.physical_constraint_fraction(controls.bath_exergy_flux);
        self.structure_s = compute_effective_structure_bits(
            self.structure_raw_s,
            rho_phys,
            adami_rho,
            replicating_pop_size,
        );

        let (guild, offense, defense, syntroph_aff) = self.evaluate_trophic_traits(controls);
        self.guild = guild;
        if self.guild == TrophicGuild::Holobiont && self.endo_tier == 0 {
            self.endo_tier = 1;
            self.organelle_membranes = 2;
        }
        self.offense_trait = offense;
        self.defense_trait = defense;

        let epistatic_score =
            self.compute_causal_epistatic_score(controls.artificial_selection, biome, env_phase);
        let mean_mem = if self.lifetime_memory.is_empty() {
            0.0
        } else {
            self.lifetime_memory.iter().sum::<f64>() / (self.lifetime_memory.len() as f64)
        };
        let mem_l1: f64 = self.lifetime_memory.iter().map(|m| m.abs()).sum();

        self.landauer_cost = if controls.bath_exergy_flux.abs() <= 0.03 {
            0.0
        } else {
            controls.landauer_bit_cost
                * (host_raw_s
                    + endo_raw_s
                    + nm_raw_s
                    + kp_somatic_s
                    + import_overhead_bits
                    + 2.8 * mem_l1)
        };

        let loc_o2 = Self::local_oxygen_at(controls.ambient_oxygen_level, biome);
        let mro_bioenergetic_mult = match self.mro_class {
            MroClass::None => 1.0,
            MroClass::AerobicMitochondrion => {
                if loc_o2 >= 0.45 {
                    1.0 + 0.14 * (loc_o2 - 0.45)
                } else {
                    // Severe hypoxia stalls aerobic Complex IV and leaks ROS!
                    (0.58 + 0.92 * loc_o2).clamp(0.58, 1.0)
                }
            }
            MroClass::AnaerobicMitochondrion => (1.04 - 0.12 * loc_o2).clamp(0.90, 1.06),
            MroClass::Hydrogenosome => {
                // Substrate-level ATP + PFO/[FeFe]-hydrogenase thrives in hypoxia/anoxia!
                (1.16 - 0.24 * loc_o2).clamp(0.90, 1.16)
            }
            MroClass::Mitosome => {
                if self.has_intact_fes_cluster_biogenesis(controls.cytosolic_suf_hgt_rescue) {
                    (1.06 - 0.10 * loc_o2).clamp(0.94, 1.06)
                } else {
                    0.22
                }
            }
            MroClass::Amiteukaryote => {
                if controls.cytosolic_suf_hgt_rescue >= 0.65 {
                    (1.10 - 0.10 * loc_o2).clamp(0.98, 1.10)
                } else {
                    0.20
                }
            }
        };

        let niche_harvest_boost = 1.0 + 0.36 * local_niche.clamp(0.0, 2.0);
        let base_eff =
            ((0.48 + 0.55 * epistatic_score + 0.28 * mean_mem * controls.memory_plasticity)
                * niche_harvest_boost
                * mro_bioenergetic_mult)
                .clamp(0.15, 2.55);

        if controls.bath_exergy_flux <= 0.03 {
            self.primary_work = controls.bath_exergy_flux;
            self.syntrophic_work = 0.0;
            self.byproduct_output = 0.0;
        } else {
            match self.guild {
                TrophicGuild::Producer => {
                    let prod_inhib = 1.0
                        / (1.0 + 0.20 * (local_byproduct - syntroph_relief_bonus).max(0.0));
                    let detox_boost = 1.0 + 0.22 * syntroph_relief_bonus.clamp(0.0, 1.2);
                    self.primary_work = local_exergy.max(0.12)
                        * base_eff
                        * prod_inhib
                        * detox_boost
                        * self.integrity;
                    self.syntrophic_work = 0.0;
                    self.byproduct_output =
                        controls.cross_feeding_rate * 0.42 * self.primary_work;
                }
                TrophicGuild::Syntroph => {
                    self.primary_work =
                        0.34 * local_exergy.max(0.10) * base_eff * self.integrity;
                    self.syntrophic_work = (1.35 * local_byproduct + 0.16 * local_niche)
                        * (0.65 + 0.85 * syntroph_aff)
                        * (1.0 + 0.45 * controls.cross_feeding_rate)
                        * self.integrity;
                    self.byproduct_output = 0.0;
                }
                TrophicGuild::Holobiont => {
                    let egt_f = self.egt_nuclear_fraction();
                    let corr_redox_coupling = if self.endosymbiont_sequence.len() >= CORR_CORE_LOCI {
                        1.0 + 0.22 * egt_f
                    } else {
                        0.66
                    };
                    let nm_streamline_bonus = if self.nucleomorph_eliminated {
                        0.08
                    } else {
                        0.0
                    };
                    let mosaic_bonus = 0.03 * (self.mosaic_donor_count().min(4) as f64);
                    let serial_niche_breadth = match self.endo_tier {
                        3..=u8::MAX => 1.24 + nm_streamline_bonus + mosaic_bonus,
                        2 => 1.15 + nm_streamline_bonus + mosaic_bonus,
                        _ => 1.0 + mosaic_bonus,
                    };
                    let coop_health =
                        (1.0 - 0.82 * self.somatic_cheater_load).clamp(0.12, 1.0);
                    let colony_harmony =
                        (1.0 - 0.85 * self.rebel_worker_load).clamp(0.14, 1.0);
                    let superorg_foraging_mult = match self.eusocial_state {
                        EusocialState::ObligateSuperorganism => {
                            // Polyandrous genetic worker diversity mildly boosts task foraging WHEN policed!
                            let gen_div_bonus = 0.08 * (0.50 - self.colony_relatedness_r).max(0.0) * 4.0;
                            (1.0 + 0.44 * self.sterile_worker_fraction + gen_div_bonus)
                                * colony_harmony
                        }
                        EusocialState::FacultativeEusocial => 1.16 * colony_harmony,
                        EusocialState::SubsocialFamily => 1.06 * colony_harmony,
                        EusocialState::Solitary => 1.0,
                    };
                    let garden_health =
                        (1.0 - 0.88 * self.escovopsis_pathogen_load).clamp(0.10, 1.0);
                    let agro_cellulose_mult = match self.agriculture_state {
                        AgricultureState::LeafCutterAgroEcosystem => {
                            (1.0 + 0.52 * self.gongylidia_yield) * garden_health
                        }
                        AgricultureState::HigherDomesticatedCrop => {
                            (1.0 + 0.36 * self.gongylidia_yield) * garden_health
                        }
                        AgricultureState::LowerFacultativeFarming => {
                            (1.0 + 0.16 * self.gongylidia_yield) * garden_health
                        }
                        AgricultureState::HunterForager => 1.0,
                    };
                    let culture_health =
                        (1.0 - 0.86 * self.cultural_drift_erosion).clamp(0.12, 1.0);
                    let conformity_factor = if controls.inter_deme_connectivity > 0.80 {
                        (1.0 - 0.34 * ((controls.inter_deme_connectivity - 0.80) / 0.20))
                            .clamp(0.64, 1.0)
                    } else if controls.inter_deme_connectivity < 0.08 {
                        0.82
                    } else {
                        1.0
                    };
                    let recomb_bonus = if (0.12..=0.75).contains(&controls.inter_deme_connectivity)
                        && self.cumulative_cross_deme_recombinations > 0
                    {
                        0.22
                    } else {
                        0.0
                    };
                    let cultural_work_mult = match self.cultural_state {
                        CulturalState::SymbolicCumulativeCulture => {
                            (1.0 + 0.46 * (self.cultural_ratchet_tier as f64 / 4.0) + recomb_bonus)
                                * conformity_factor
                                * culture_health
                        }
                        CulturalState::ImitativeRatchet => {
                            1.24 * conformity_factor * culture_health
                        }
                        CulturalState::LowFidelityTradition => 1.07 * culture_health,
                        CulturalState::AsocialTrialError => 1.0,
                    };
                    let multicell_foraging_mult = (match self.multicell_state {
                        MulticellularState::GermSomaMulticellular => {
                            (1.0 + 0.36 * self.soma_fraction) * coop_health
                        }
                        MulticellularState::UndifferentiatedClonal => 1.10 * coop_health,
                        MulticellularState::AggregativeChimerism => 0.96 * coop_health,
                        MulticellularState::Unicellular => 1.0,
                    }) * superorg_foraging_mult
                        * agro_cellulose_mult
                        * cultural_work_mult
                        * {
                            let arch_purity =
                                (1.0 - 0.90 * self.archive_false_discovery_load).clamp(0.10, 1.0);
                            let exo_m = match self.exosomatic_state {
                                ExosomaticState::PeerVerifiedInstitution => 1.52 * arch_purity,
                                ExosomaticState::WrittenArchive => 1.28 * arch_purity,
                                ExosomaticState::MaterialNotation => 1.09 * arch_purity,
                                ExosomaticState::EndosomaticBound => 1.0,
                            };
                            let tail_ok = self.distribution_tail_diversity.clamp(0.08, 1.0)
                                * (1.0 - 0.90 * self.recursive_autophagy_load).clamp(0.08, 1.0);
                            let symbio_ai_m = match self.machine_synthesis_state {
                                MachineSynthesisState::MachineVerifiedSymbiosis => 1.58 * tail_ok,
                                MachineSynthesisState::UnverifiedRecursiveGen => 1.12 * tail_ok,
                                MachineSynthesisState::AlgorithmicMediation => 1.16 * tail_ok,
                                MachineSynthesisState::ManualArchival => 1.0,
                            };
                            let gaia_m = match self.planetary_epoch {
                                PlanetaryState::HomeostaticAwakening => 1.38,
                                PlanetaryState::OxygenatedPhanerozoic => 1.0,
                                PlanetaryState::AnoxicArchean => 0.88,
                                PlanetaryState::AnthropoceneDisequilibrium => 0.72,
                            };
                            let abio_m = match self.abiogenesis_stage {
                                AbiogenesisStage::PrebioticMonomerSoup => 0.28,
                                AbiogenesisStage::AutocatalyticRafProtocell => 0.52,
                                AbiogenesisStage::RnaRibozymeHypercycle => 0.74,
                                AbiogenesisStage::ProgenoteCommunity => 0.88,
                                AbiogenesisStage::LucaCellularLineage => 1.0,
                            };
                            let chiral_m = match self.chirality_code_stage {
                                ChiralityCodeStage::RacemicCrossInhibited => 0.22,
                                ChiralityCodeStage::HomochiralRnaSymmetryBroken => 0.54,
                                ChiralityCodeStage::StereoselectiveMinihelixAcylated => 0.76,
                                ChiralityCodeStage::AmbiguousStatisticalCode => 0.86,
                                ChiralityCodeStage::UniversalErrorMinimizedSgc => 1.0,
                            };
                            let dna_m = match self.genome_substrate_stage {
                                GenomeSubstrateStage::LabileRnaChromosome => 0.46,
                                GenomeSubstrateStage::ViralModifiedUDna => 0.68,
                                GenomeSubstrateStage::DeaminationRepairableTDna => 0.88,
                                GenomeSubstrateStage::CoevolvedCellularVirosphere => 1.0,
                            };
                            exo_m * symbio_ai_m * gaia_m * abio_m * chiral_m * dna_m
                        };
                    self.primary_work = 1.06
                        * local_exergy.max(0.12)
                        * base_eff
                        * corr_redox_coupling
                        * serial_niche_breadth
                        * multicell_foraging_mult
                        * self.integrity;
                    let internal_pool =
                        0.42 * self.primary_work * (0.65 + 0.75 * controls.cross_feeding_rate);
                    let ambient_uptake =
                        0.55 * local_byproduct * syntroph_aff * serial_niche_breadth;
                    self.syntrophic_work = (0.88 * internal_pool + ambient_uptake)
                        * (0.85 + 0.55 * syntroph_aff)
                        * corr_redox_coupling
                        * self.integrity;
                    // Hydrogenosomes evolve abundant molecular H2 + acetate via PFO + [FeFe]-hydrogenase,
                    // feeding neighboring syntrophic methanogens/cross-feeders!
                    let h2_efflux_mult = if self.mro_class == MroClass::Hydrogenosome {
                        2.15
                    } else {
                        1.0
                    };
                    self.byproduct_output = 0.10 * internal_pool * h2_efflux_mult;
                }
                TrophicGuild::PredatorParasite => {
                    if self.has_apicoplast() || self.apicoplast_eliminated {
                        let anabolic_viability = if self.has_apicoplast() {
                            if self.has_intact_apicoplast_essential_core() {
                                1.0
                            } else if controls.host_metabolite_rescue >= 0.65 {
                                0.96
                            } else {
                                0.18
                            }
                        } else if controls.host_metabolite_rescue >= 0.65 {
                            1.02
                        } else {
                            0.20
                        };
                        self.primary_work = 0.12
                            * local_exergy.max(0.08)
                            * base_eff
                            * anabolic_viability
                            * self.integrity;
                        self.syntrophic_work = (0.38 * local_byproduct
                            + 0.24 * controls.host_metabolite_rescue)
                            * anabolic_viability
                            * self.integrity;
                        self.byproduct_output = if self.mro_class == MroClass::Hydrogenosome {
                            0.14 * self.syntrophic_work
                        } else {
                            0.0
                        };
                    } else {
                        let kp_transporter_boost = if self.has_active_kleptoplast() {
                            self.kleptoplast_integrity
                                * (0.62 + 0.15 * (self.klepto_hgt_loci.min(5) as f64))
                        } else {
                            0.0
                        };
                        self.primary_work = (0.22 + 0.68 * kp_transporter_boost)
                            * local_exergy.max(0.08)
                            * base_eff
                            * self.integrity;
                        self.syntrophic_work = (0.18 + 0.44 * kp_transporter_boost)
                            * local_byproduct
                            * self.integrity;
                        self.byproduct_output =
                            controls.cross_feeding_rate * 0.18 * self.primary_work;
                    }
                }
            }
        }

        self.harvested_work = self.primary_work + self.syntrophic_work + self.net_trophic_flux;

        let nc_motif_fraction = self.niche_construction_trait_fraction();
        let guild_nc_mult = match self.guild {
            TrophicGuild::Producer => 1.0,
            TrophicGuild::Syntroph => 1.15,
            TrophicGuild::Holobiont => {
                let ecm_soma_boost = if self.multicell_state == MulticellularState::GermSomaMulticellular {
                    0.42 * self.soma_fraction * (1.0 - 0.75 * self.somatic_cheater_load).clamp(0.1, 1.0)
                } else if self.multicell_state == MulticellularState::UndifferentiatedClonal {
                    0.16
                } else {
                    0.0
                };
                let nest_fortress_boost = match self.eusocial_state {
                    EusocialState::ObligateSuperorganism => {
                        0.58 * self.total_sterile_caste_fraction()
                            * (1.0 - 0.80 * self.rebel_worker_load).clamp(0.1, 1.0)
                    }
                    EusocialState::FacultativeEusocial => 0.24,
                    EusocialState::SubsocialFamily => 0.10,
                    EusocialState::Solitary => 0.0,
                };
                let fungus_garden_boost = match self.agriculture_state {
                    AgricultureState::LeafCutterAgroEcosystem => {
                        0.64 * self.gongylidia_yield
                            * (1.0 - 0.85 * self.escovopsis_pathogen_load).clamp(0.1, 1.0)
                    }
                    AgricultureState::HigherDomesticatedCrop => {
                        0.40 * self.gongylidia_yield
                            * (1.0 - 0.85 * self.escovopsis_pathogen_load).clamp(0.1, 1.0)
                    }
                    AgricultureState::LowerFacultativeFarming => {
                        0.18 * (1.0 - 0.80 * self.escovopsis_pathogen_load).clamp(0.1, 1.0)
                    }
                    AgricultureState::HunterForager => 0.0,
                };
                let epistemic_niche_boost = match self.cultural_state {
                    CulturalState::SymbolicCumulativeCulture => {
                        (0.46 + 0.08 * (self.cultural_ratchet_tier as f64))
                            * (1.0 - 0.85 * self.cultural_drift_erosion).clamp(0.1, 1.0)
                    }
                    CulturalState::ImitativeRatchet => {
                        0.22 * (1.0 - 0.85 * self.cultural_drift_erosion).clamp(0.1, 1.0)
                    }
                    CulturalState::LowFidelityTradition => 0.08,
                    CulturalState::AsocialTrialError => 0.0,
                };
                let archival_niche_boost = match self.exosomatic_state {
                    ExosomaticState::PeerVerifiedInstitution => {
                        0.72 * (1.0 - 0.88 * self.archive_false_discovery_load).clamp(0.1, 1.0)
                    }
                    ExosomaticState::WrittenArchive => {
                        0.44 * (1.0 - 0.88 * self.archive_false_discovery_load).clamp(0.1, 1.0)
                    }
                    ExosomaticState::MaterialNotation => 0.16,
                    ExosomaticState::EndosomaticBound => 0.0,
                };
                let formal_symbio_niche_boost = match self.machine_synthesis_state {
                    MachineSynthesisState::MachineVerifiedSymbiosis => {
                        0.84 * self.distribution_tail_diversity.clamp(0.1, 1.0)
                            * (1.0 - 0.88 * self.recursive_autophagy_load).clamp(0.1, 1.0)
                    }
                    MachineSynthesisState::AlgorithmicMediation => 0.24,
                    MachineSynthesisState::UnverifiedRecursiveGen => {
                        0.18 * (1.0 - 0.90 * self.recursive_autophagy_load).clamp(0.05, 1.0)
                    }
                    MachineSynthesisState::ManualArchival => 0.0,
                };
                1.45
                    + 0.18 * self.egt_nuclear_fraction()
                    + 0.10 * ((self.endo_tier.saturating_sub(1)) as f64)
                    + ecm_soma_boost
                    + nest_fortress_boost
                    + fungus_garden_boost
                    + epistemic_niche_boost
                    + archival_niche_boost
                    + formal_symbio_niche_boost
            }
            TrophicGuild::PredatorParasite => {
                if self.has_active_kleptoplast() {
                    0.65
                } else {
                    0.25
                }
            }
        };
        self.niche_output = if controls.bath_exergy_flux > 0.03
            && self.harvested_work > self.landauer_cost
        {
            controls.niche_construction_rate
                * guild_nc_mult
                * nc_motif_fraction
                * (0.5 + 0.5 * mean_mem + 0.4 * epistatic_score)
                * self.integrity
        } else {
            0.0
        };

        let net_exergy_margin = self.harvested_work - self.landauer_cost;
        let exergy_factor = if controls.bath_exergy_flux.abs() <= 0.03 {
            0.30
        } else if net_exergy_margin < -0.02 {
            (0.28 * self.integrity).clamp(0.008, 0.32)
        } else {
            (0.76 + 0.48 * net_exergy_margin.min(1.6)).clamp(0.65, 1.72)
        };

        let engineered_bonus = if controls.artificial_selection {
            1.16
        } else {
            1.0
        };
        let mro_streamline_synergy = match self.mro_class {
            MroClass::None | MroClass::AerobicMitochondrion => 1.0,
            MroClass::AnaerobicMitochondrion => 1.04,
            MroClass::Hydrogenosome => {
                if loc_o2 < 0.45 {
                    1.12
                } else {
                    0.96
                }
            }
            MroClass::Mitosome => {
                if self.has_intact_fes_cluster_biogenesis(controls.cytosolic_suf_hgt_rescue) {
                    1.16
                } else {
                    0.35
                }
            }
            MroClass::Amiteukaryote => {
                if controls.cytosolic_suf_hgt_rescue >= 0.65 {
                    1.24
                } else {
                    0.32
                }
            }
        };
        let planetary_gaia_synergy = match self.planetary_epoch {
            PlanetaryState::HomeostaticAwakening => 1.36,
            PlanetaryState::OxygenatedPhanerozoic => 1.0,
            PlanetaryState::AnoxicArchean => 0.92,
            PlanetaryState::AnthropoceneDisequilibrium => 0.74,
        };
        let machine_symbio_synergy = planetary_gaia_synergy * match self.machine_synthesis_state {
            MachineSynthesisState::MachineVerifiedSymbiosis => {
                1.54 * self.distribution_tail_diversity.clamp(0.12, 1.0)
                    * (1.0 - 0.90 * self.recursive_autophagy_load).clamp(0.10, 1.0)
            }
            MachineSynthesisState::AlgorithmicMediation => 1.15,
            MachineSynthesisState::UnverifiedRecursiveGen => {
                1.10 * (1.0 - 0.90 * self.recursive_autophagy_load).clamp(0.08, 1.0)
            }
            MachineSynthesisState::ManualArchival => 1.0,
        };
        let exosomatic_eti_synergy = machine_symbio_synergy * match self.exosomatic_state {
            ExosomaticState::PeerVerifiedInstitution => {
                1.48 * (1.0 - 0.88 * self.archive_false_discovery_load).clamp(0.12, 1.0)
            }
            ExosomaticState::WrittenArchive => {
                1.24 * (1.0 - 0.88 * self.archive_false_discovery_load).clamp(0.14, 1.0)
            }
            ExosomaticState::MaterialNotation => 1.08,
            ExosomaticState::EndosomaticBound => 1.0,
        };
        let cultural_eti_synergy = exosomatic_eti_synergy * match self.cultural_state {
            CulturalState::SymbolicCumulativeCulture => {
                let ch = (1.0 - 0.86 * self.cultural_drift_erosion).clamp(0.12, 1.0);
                let cf = if controls.inter_deme_connectivity > 0.80 {
                    0.74
                } else if controls.inter_deme_connectivity < 0.08 {
                    0.82
                } else {
                    1.0 + 0.14 * f64::from(self.cumulative_cross_deme_recombinations > 0)
                };
                (1.0 + 0.36 * (self.cultural_ratchet_tier as f64 / 4.0)) * cf * ch
            }
            CulturalState::ImitativeRatchet => {
                1.18 * (1.0 - 0.85 * self.cultural_drift_erosion).clamp(0.16, 1.0)
            }
            CulturalState::LowFidelityTradition => 1.06,
            CulturalState::AsocialTrialError => 1.0,
        };
        let agro_eti_synergy = cultural_eti_synergy * match self.agriculture_state {
            AgricultureState::LeafCutterAgroEcosystem => {
                let gh = (1.0 - 0.88 * self.escovopsis_pathogen_load).clamp(0.12, 1.0);
                (1.0 + 0.38 * self.gongylidia_yield) * gh
            }
            AgricultureState::HigherDomesticatedCrop => {
                let gh = (1.0 - 0.88 * self.escovopsis_pathogen_load).clamp(0.14, 1.0);
                (1.0 + 0.24 * self.gongylidia_yield) * gh
            }
            AgricultureState::LowerFacultativeFarming => {
                1.08 * (1.0 - 0.80 * self.escovopsis_pathogen_load).clamp(0.18, 1.0)
            }
            AgricultureState::HunterForager => 1.0,
        };
        let superorg_eti_synergy = agro_eti_synergy * match self.eusocial_state {
            EusocialState::ObligateSuperorganism => {
                let harmony = (1.0 - 0.86 * self.rebel_worker_load).clamp(0.16, 1.0);
                (1.0 + 0.42 * self.total_sterile_caste_fraction()) * harmony
            }
            EusocialState::FacultativeEusocial => {
                1.14 * (1.0 - 0.80 * self.rebel_worker_load).clamp(0.20, 1.0)
            }
            EusocialState::SubsocialFamily => 1.05,
            EusocialState::Solitary => 1.0,
        };
        let multicell_eti_synergy = superorg_eti_synergy * match self.multicell_state {
            MulticellularState::GermSomaMulticellular => {
                let coop = (1.0 - 0.85 * self.somatic_cheater_load).clamp(0.18, 1.0);
                (1.0 + 0.34 * self.soma_fraction) * coop
            }
            MulticellularState::UndifferentiatedClonal => {
                1.10 * (1.0 - 0.75 * self.somatic_cheater_load).clamp(0.22, 1.0)
            }
            MulticellularState::AggregativeChimerism => {
                0.95 * (1.0 - 0.85 * self.somatic_cheater_load).clamp(0.18, 1.0)
            }
            MulticellularState::Unicellular => 1.0,
        };
        let holobiont_synergy = (if self.guild == TrophicGuild::Holobiont {
            if self.endosymbiont_sequence.len() >= CORR_CORE_LOCI {
                let tier_gain = 0.12 * ((self.endo_tier.saturating_sub(1)) as f64);
                let nm_clean_gain = if self.nucleomorph_eliminated {
                    0.10
                } else {
                    0.0
                };
                (1.26 + 0.22 * self.egt_nuclear_fraction() + tier_gain + nm_clean_gain)
                    * multicell_eti_synergy
            } else {
                0.98 * multicell_eti_synergy
            }
        } else if self.has_apicoplast() {
            if self.has_intact_apicoplast_essential_core() {
                let streamline_gain =
                    0.06 * ((4u8.saturating_sub(self.apicoplast_pathway_count())) as f64);
                1.22 + streamline_gain
            } else if controls.host_metabolite_rescue >= 0.65 {
                1.20
            } else {
                0.42
            }
        } else if self.apicoplast_eliminated {
            if controls.host_metabolite_rescue >= 0.65 {
                1.36
            } else {
                0.40
            }
        } else if self.has_active_kleptoplast() {
            1.0 + 0.16 * self.kleptoplast_integrity
        } else {
            1.0
        }) * mro_streamline_synergy;

        let sublinear_scaffold = 4.25 * self.structure_s.powf(0.62);
        let niche_scaffold = 1.0 + 0.22 * local_niche.clamp(0.0, 2.0);
        self.base_complexity_c = (sublinear_scaffold
            * niche_scaffold
            * (0.38 + 1.45 * epistatic_score)
            * exergy_factor
            * engineered_bonus
            * holobiont_synergy
            * self.integrity)
            .max(0.08);

        let plastic_amplification = 1.0 + 0.72 * mean_mem * controls.memory_plasticity.max(0.0);
        self.complexity_c = (self.base_complexity_c * plastic_amplification).max(0.08);

        let leverage = self.complexity_c / self.structure_s.max(2.0);
        let energy_viability = (0.55 + 0.65 * net_exergy_margin).clamp(0.05, 1.95);
        self.fitness =
            ((0.34 * epistatic_score + 0.46 * (leverage / 2.4) + 0.20 * energy_viability)
                * self.integrity)
                .clamp(0.02, 4.6);
    }

    fn niche_construction_trait_fraction(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let builders: &[u16] = match self.kind {
            OrganismKind::RnaDnaProtein => &[4, 5, 8, 11],
            OrganismKind::HominidPopulation => &[6, 7, 9],
            OrganismKind::ArtificialIntelligence => &[6, 7, 8],
        };
        let hits = self
            .sequence
            .iter()
            .chain(self.endosymbiont_sequence.iter())
            .chain(self.nucleomorph_sequence.iter())
            .filter(|s| builders.contains(s))
            .count() as f64;
        let denom = (self.sequence.len()
            + self.endosymbiont_sequence.len()
            + self.nucleomorph_sequence.len())
        .max(1) as f64;
        (hits / (denom * 0.35)).clamp(0.05, 1.35)
    }

    fn compute_causal_epistatic_score(
        &self,
        artificial_target: bool,
        biome: BiomeZone,
        env_phase: f64,
    ) -> f64 {
        let n = self.sequence.len();
        if n == 0 {
            return 0.05;
        }
        let mut unique_mask = 0u32;
        for &sym in self
            .sequence
            .iter()
            .chain(self.endosymbiont_sequence.iter())
            .chain(self.nucleomorph_sequence.iter())
            .chain(self.kleptoplast_sequence.iter())
        {
            if (sym as u32) < 31 {
                unique_mask |= 1u32 << (sym as u32);
            }
        }
        let diversity = (unique_mask.count_ones() as f64)
            / (self.kind.alphabet_symbols().len().max(1) as f64);
        let wiring_density = (self.wiring.len() as f64) / (n.max(1) as f64);

        let biome_targets = self.kind.biome_target_symbols(biome, env_phase);
        let raw_biome_match = self
            .sequence
            .iter()
            .filter(|&&s| biome_targets.contains(&s))
            .count() as f64
            / (n as f64);
        let biome_match = if self.endo_tier >= 2 {
            raw_biome_match.max(0.36 + 0.06 * ((self.endo_tier - 1) as f64))
        } else {
            raw_biome_match
        };
        let division_of_labor = if self.has_apicoplast() {
            0.22
        } else if !self.endosymbiont_sequence.is_empty()
            || self.egt_transferred_loci >= 2
            || self.endo_tier >= 2
            || self.mro_class != MroClass::None
        {
            0.16
                + 0.14 * self.egt_nuclear_fraction()
                + 0.07 * ((self.endo_tier.saturating_sub(1)) as f64)
        } else if self.has_active_kleptoplast() {
            0.09 * self.kleptoplast_integrity
        } else {
            0.0
        };

        match self.kind {
            OrganismKind::RnaDnaProtein => {
                let stem_energy = (wiring_density * 1.35).min(0.45);
                let has_nucleic = self.sequence.iter().any(|&x| x <= 3);
                let has_his_or_glu = self.sequence.contains(&8) || self.sequence.contains(&7);
                let has_acid = self.sequence.contains(&6) || self.sequence.contains(&9);
                let has_nuc_res = self.sequence.iter().any(|&x| x == 5 || x == 10);
                let triad_epistasis = match (has_his_or_glu, has_acid, has_nuc_res) {
                    (true, true, true) => 0.40,
                    (true, true, false) | (true, false, true) => 0.18,
                    _ => 0.04,
                };
                let rnp_coop = if has_nucleic && triad_epistasis > 0.2 {
                    0.22
                } else {
                    0.0
                };
                let selex = if artificial_target && biome_match > 0.45 {
                    0.20
                } else {
                    0.0
                };
                (0.16
                    + 0.26 * diversity
                    + stem_energy
                    + triad_epistasis
                    + rnp_coop
                    + 0.42 * biome_match
                    + division_of_labor
                    + selex)
                    .clamp(0.08, 2.08)
            }
            OrganismKind::HominidPopulation => {
                let causal_chain = [1u16, 2, 3, 4, 5];
                let mut stage_cursor = 0usize;
                let mut reductive_fractures = 0usize;
                for &op in &self.sequence {
                    if stage_cursor < causal_chain.len() && op == causal_chain[stage_cursor] {
                        stage_cursor += 1;
                    } else if op == 5 && stage_cursor < 3 {
                        reductive_fractures += 1;
                    }
                }
                let chain_progress = (stage_cursor as f64) / (causal_chain.len() as f64);
                let fracture_penalty = (reductive_fractures as f64 * 0.06).min(0.24);
                let social_institutions = [6u16, 7, 8, 9]
                    .iter()
                    .filter(|&&s| {
                        self.sequence.contains(&s)
                            || self.endosymbiont_sequence.contains(&s)
                            || self.nucleomorph_sequence.contains(&s)
                    })
                    .count() as f64
                    / 4.0;
                let norm_synergy = if stage_cursor >= 4 && self.sequence.contains(&8) {
                    0.30
                } else {
                    0.0
                };
                let inst_bonus = if artificial_target && stage_cursor == 5 {
                    0.18
                } else {
                    0.0
                };
                (0.16
                    + 0.46 * chain_progress
                    + 0.32 * social_institutions
                    + 0.32 * wiring_density
                    + norm_synergy
                    + 0.38 * biome_match
                    + division_of_labor
                    + inst_bonus
                    - fracture_penalty)
                    .clamp(0.08, 2.08)
            }
            OrganismKind::ArtificialIntelligence => {
                let valid_framing =
                    self.sequence.first() == Some(&0) && self.sequence.last() == Some(&9);
                let pos_attn = self.sequence.iter().position(|&x| x == 4);
                let pos_mem_w = self.sequence.iter().position(|&x| x == 6);
                let pos_mem_r = self.sequence.iter().rposition(|&x| x == 5);
                let has_verify = self.sequence.contains(&8);

                let causal_mem_bonus = match (pos_mem_w, pos_mem_r) {
                    (Some(w), Some(r)) if w < r => 0.32,
                    (Some(_), Some(_)) => 0.15,
                    _ => 0.0,
                };
                let pipeline_epistasis = (if valid_framing { 0.20 } else { 0.0 })
                    + (if pos_attn.is_some() { 0.20 } else { 0.0 })
                    + causal_mem_bonus
                    + (if has_verify { 0.16 } else { 0.0 });
                let meta_bonus = if artificial_target && has_verify {
                    0.20
                } else {
                    0.0
                };
                (0.14
                    + 0.26 * diversity
                    + 0.42 * wiring_density
                    + pipeline_epistasis
                    + 0.38 * biome_match
                    + division_of_labor
                    + meta_bonus)
                    .clamp(0.08, 2.10)
            }
        }
    }
}
