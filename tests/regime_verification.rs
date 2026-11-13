//! Automated verification suite for ITME v12.0 (`v1.2.0`):
//! 1. Price Equation exact algebraic identity `Delta bar_z = Cov(w,z)/bar_w + E[w dz]/bar_w`
//! 2. Landauer thermodynamic starvation when maintenance cost exceeds harvested work
//! 3. Purely observable emergent classification across all 13 regimes and 3 organisms
//! 4. Spatial niche construction accumulation & stability buffering
//! 5. Dual-biome parapatric/allopatric speciation (`species_count >= 2`, `F_ST > 0.15`)
//! 6. Multi-trophic Holling Type II predation/parasitism & Black-Queen syntrophic cross-feeding
//! 7. Real-time pruned phylogenetic lineage tracking, `T_MRCA`, Sackin balance, Newick/NHX & SVG export
//! 8. Endosymbiotic `Producer + Syntroph -> 1° Holobiont` fusion & Cardona eNewick `#H1` reticulate network
//! 9. Gradual Endosymbiotic Gene Transfer (`EGT`), John F. Allen's 3-locus CoRR redox core retention,
//!    Landauer copy-economy payoff, and irreversible organellogenesis lock-in (`EGT->Nuc`)
//! 10. Serial Secondary (`2°`) & Tertiary (`3°`) Endosymbiosis, 4-membrane `SELMA` complex plastids,
//!     periplastidial Nucleomorph compaction (`5`-locus bonsai arrest), and `NM -> 0` elimination
//! 11. Mixotrophic Kleptoplasty turnover, Shopping-Bag `HGT` mosaicism, and `KP->FIX` tipping point
//! 12. Reductive Evolution of Photosynthetic Holobionts into Parasitic `Apicoplast`-Bearing Lineages
//!     (`Chromera -> Toxoplasma / Plasmodium`), `MEP/IPP + SUF/Fe-S` delayed-death essentiality vs.
//!     `IPP`-rescued *Cryptosporidium* loss (`Apico -> 0`), and algal `#H` reticulation preservation
//! 13. Anaerobic Mitochondrial Reduction (`AerobicMito -> Hydrogenosome -> Mitosome -> Amiteukaryote`),
//!     `CoRR` `mtDNA` sweep upon respiratory ETC loss, `H2` syntrophic cross-feeding, indispensable
//!     `ISC` `[2Fe-2S]` -> cytosolic `CIA` nuclear replicative fidelity floor vs. bacterial cytosolic
//!     `SUF` HGT rescue (*Monocercomonoides* `MRO -> 0`), and `#H` reticulation preservation!

use itme_evolution_sim::engine::{SimulationMode, SimulationWorld};
use itme_evolution_sim::model::{MroClass, PriceDecomposition, Rng, TrophicGuild, MulticellularState, EusocialState, AgricultureState, CulturalState, ExosomaticState, MachineSynthesisState, PlanetaryBiosphere, PlanetaryState, AbiogenesisStage, HadeanGeochemicalBath, ChiralCodeReactor, ChiralityCodeStage, GenomeSubstrateStage, RnrVirosphereReactor};
use itme_evolution_sim::organisms::{
    ModeledUnit, OrganismKind, APICOPLAST_MIN_LOCI, APICO_PATH_ALL, APICO_PATH_ESSENTIAL_CORE,
    CORR_CORE_LOCI, KLEPTO_FIXATION_HGT_THRESHOLD, MRO_PATH_ALL, MRO_PATH_HYDROGENOSOME,
    MRO_PATH_MITOSOME_CORE, NUCLEOMORPH_HOUSEKEEPING_LOCI,
};
use itme_evolution_sim::phylogeny::{ExtantPhyloState, PhylogenyTracker};
use itme_evolution_sim::regimes::RegimeId;
use itme_evolution_sim::spatial::{analyze_spatial_speciation, BiomeZone};

#[test]
fn price_equation_satisfies_exact_evolutionary_identity() {
    let records = [
        (0.6, 1.1, 1.05),
        (1.2, 1.8, 1.85),
        (2.1, 2.6, 2.75),
        (0.8, 1.3, 1.25),
    ];
    let p = PriceDecomposition::from_parent_offspring_records(&records);
    let n = records.len() as f64;
    let w_bar = records.iter().map(|r| r.0).sum::<f64>() / n;
    let z_bar = records.iter().map(|r| r.1).sum::<f64>() / n;
    let z_next = records.iter().map(|r| (r.0 / w_bar) * r.2).sum::<f64>() / n;
    let lhs = z_next - z_bar;
    let rhs = p.selection_covariance + p.transmission_bias;
    assert!(
        (lhs - rhs).abs() < 1e-10,
        "Price identity violated: lhs={lhs}, rhs={rhs}"
    );
}

#[test]
fn landauer_maintenance_overload_drives_entropic_decay() {
    let mut world = SimulationWorld::new(
        OrganismKind::ArtificialIntelligence,
        RegimeId::R10AdaptiveDissipation,
        SimulationMode::OrganicSandbox,
        20,
        77,
    );
    world.controls.landauer_bit_cost = 0.045;
    for _ in 0..12 {
        world.step();
    }
    assert_eq!(
        world.classified_regime,
        RegimeId::R05Decaying,
        "Excessive Landauer bit cost must trigger emergent R(0.5) Entropic Decay"
    );
}

#[test]
fn all_thirteen_emergent_regimes_classify_cleanly_across_all_three_organisms() {
    for &org in &OrganismKind::ALL {
        let mut world = SimulationWorld::new(
            org,
            RegimeId::R00Static,
            SimulationMode::OrganicSandbox,
            20,
            42,
        );
        for &regime in &RegimeId::ALL {
            world.apply_regime_preset(regime, false);
            for _ in 0..16 {
                world.step();
            }
            assert_eq!(
                world.classified_regime,
                regime,
                "Emergent classifier mismatch for {:?} in {}: got {}",
                org,
                regime.code(),
                world.classified_regime.code()
            );
        }
    }
}

#[test]
fn spatial_niche_construction_accumulates_and_lifts_ecosystem_order() {
    let mut world_no_niche = SimulationWorld::new(
        OrganismKind::HominidPopulation,
        RegimeId::R30SoftInheritance,
        SimulationMode::OrganicSandbox,
        40,
        2026,
    );
    world_no_niche.controls.niche_construction_rate = 0.0;

    let mut world_niche = SimulationWorld::new(
        OrganismKind::HominidPopulation,
        RegimeId::R30SoftInheritance,
        SimulationMode::OrganicSandbox,
        40,
        2026,
    );
    world_niche.controls.niche_construction_rate = 0.34;

    for _ in 0..36 {
        world_no_niche.step();
        world_niche.step();
    }

    let m_off = world_no_niche.latest_metrics();
    let m_on = world_niche.latest_metrics();
    assert!(
        m_on.mean_niche_capital > 0.08 && m_on.flags.niche_construction,
        "Expected active niche construction capital > 0.08, got {:.4}",
        m_on.mean_niche_capital
    );
    assert!(
        m_on.order_o > m_off.order_o,
        "Constructed niche should lift Order O ({:.3} vs {:.3})",
        m_on.order_o,
        m_off.order_o
    );
}

#[test]
fn divergent_biomes_yield_multiple_species_and_high_fst_while_preserving_within_species_adami() {
    let alpha_ecotype_1: [u16; 12] = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1];
    let alpha_ecotype_2: [u16; 12] = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 2];
    let beta_ecotype_1: [u16; 12] = [6, 7, 8, 9, 5, 6, 7, 8, 9, 5, 6, 7];
    let beta_ecotype_2: [u16; 12] = [6, 7, 8, 9, 5, 6, 7, 8, 9, 5, 6, 8];

    let seqs: Vec<&[u16]> = vec![
        &alpha_ecotype_1,
        &alpha_ecotype_2,
        &beta_ecotype_1,
        &beta_ecotype_2,
    ];
    let biomes = [
        BiomeZone::Alpha,
        BiomeZone::Alpha,
        BiomeZone::Beta,
        BiomeZone::Beta,
    ];

    let spec = analyze_spatial_speciation(&seqs, &biomes, 10, 0.25);
    assert_eq!(spec.species_count, 2, "Expected 2 reproductively isolated species");
    assert!(spec.fst_index > 0.65);
    assert!(spec.within_species_adami_rho > 0.85);
}

#[test]
fn multi_trophic_predation_and_cross_feeding_establish_active_food_web() {
    let mut world = SimulationWorld::new(
        OrganismKind::RnaDnaProtein,
        RegimeId::R30FullyDarwinian,
        SimulationMode::OrganicSandbox,
        40,
        202605,
    );
    world.controls.trophic_coupling = 0.42;
    world.controls.cross_feeding_rate = 0.36;

    for _ in 0..32 {
        world.step();
    }

    let m = world.latest_metrics();
    assert!(m.producer_count > 0 && (m.syntroph_count > 0 || m.predator_count > 0));
    assert!(m.mean_trophic_level > 1.05 && m.food_web_entropy_bits > 0.25);
}

#[test]
fn phylogeny_tracker_records_speciation_guild_shifts_prunes_extinct_and_exports_nhx_svg() {
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(0, &[(1, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.4)]);

    tracker.record_birth(
        2,
        1,
        None,
        5,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        0,
        TrophicGuild::Producer,
        0.08,
        0.25,
        1.5,
    );
    tracker.record_birth(
        3,
        1,
        None,
        8,
        1,
        TrophicGuild::Syntroph,
        BiomeZone::Beta,
        0,
        TrophicGuild::Producer,
        0.34,
        0.25,
        1.8,
    );
    tracker.record_birth(
        4,
        2,
        None,
        10,
        0,
        TrophicGuild::PredatorParasite,
        BiomeZone::Alpha,
        0,
        TrophicGuild::Producer,
        0.30,
        0.25,
        1.6,
    );
    tracker.record_birth(
        5,
        2,
        Some(3),
        14,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        0,
        TrophicGuild::Producer,
        0.05,
        0.25,
        1.9,
    );

    tracker.sync_extant_population(
        20,
        &[
            (3, 1, TrophicGuild::Syntroph, BiomeZone::Beta, 1.8),
            (5, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.9),
        ],
    );

    assert!(
        !tracker.nodes.contains_key(&4),
        "Extinct dead-end taxon #4 should be pruned"
    );
    let sum = tracker.compute_summary(20);
    assert_eq!(sum.mrca_id, 1, "MRCA of #3 and #5 must be founder #1");
    assert_eq!(sum.tmrca_steps, 20);
    assert!(sum.prod_to_syntroph_count >= 1 && sum.to_predator_count >= 1);
    assert!(sum.total_speciation_splits >= 1);

    let nhx = tracker.to_newick(20, true);
    assert!(nhx.ends_with(';') && nhx.contains("&&NHX:") && nhx.contains("SHIFT=Prod_to_Synt"));
    let open_p = nhx.chars().filter(|&c| c == '(').count();
    let close_p = nhx.chars().filter(|&c| c == ')').count();
    assert_eq!(open_p, close_p, "Unbalanced parentheses in NHX output: {nhx}");

    let svg = tracker.to_svg_phylogram("Unit Test Phylogram", 20);
    assert!(svg.starts_with("<svg") && svg.contains("</svg>") && svg.contains("Prod→Synt"));
}

#[test]
fn endosymbiotic_fusion_merges_producer_and_syntroph_into_holobiont_and_preserves_reticulate_enewick() {
    let mut rng = Rng::from_seed(20261112);
    let ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    let host = ModeledUnit::new_founder_with_role(
        10,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let symbiont = ModeledUnit::new_founder_with_role(
        11,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );

    let (holo, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &symbiont,
        12,
        BiomeZone::Alpha,
        1.1,
        0.35,
        0.45,
        &ctrl,
        0.9,
        12,
        0.4,
    );
    assert_eq!(holo.guild, TrophicGuild::Holobiont);
    assert_eq!(holo.endo_tier, 1);
    assert_eq!(holo.organelle_membranes, 2);
    assert!(!holo.endosymbiont_sequence.is_empty());
    assert!(gamma > 0.15 && gamma < 0.65);
    assert!(holo.harvested_work > host.harvested_work);

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(0, &[(1, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.4)]);
    tracker.record_birth(
        2,
        1,
        None,
        5,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        0,
        TrophicGuild::Producer,
        0.06,
        0.25,
        1.5,
    );
    tracker.record_birth(
        3,
        1,
        None,
        6,
        1,
        TrophicGuild::Syntroph,
        BiomeZone::Beta,
        0,
        TrophicGuild::Producer,
        0.32,
        0.25,
        1.6,
    );
    let hid = tracker.record_endosymbiotic_fusion(
        6,
        2,
        3,
        14,
        0,
        BiomeZone::Alpha,
        gamma,
        holo.endosymbiont_sequence.len(),
        2.3,
    );
    assert_eq!(hid, 1);

    tracker.sync_extant_population_with_egt(
        24,
        &[(
            6,
            0,
            TrophicGuild::Holobiont,
            BiomeZone::Alpha,
            2.3,
            0,
            holo.endosymbiont_sequence.len(),
        )],
    );
    assert!(tracker.nodes.contains_key(&3));

    let enewick = tracker.to_newick(24, true);
    assert!(enewick.matches("#H1").count() == 2);
    assert!(enewick.contains("FUSION=Prod+Synt->Holo") && enewick.contains("ROLE=Endosymbiont_Donor"));
}

#[test]
fn endosymbiotic_gene_transfer_migrates_loci_to_host_respects_corr_core_and_locks_organellogenesis() {
    let mut rng = Rng::from_seed(20261115);
    let ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    let host = ModeledUnit::new_founder_with_role(
        20,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        5,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let sym = ModeledUnit::new_founder_with_role(
        21,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        5,
        6,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );

    let (mut holo, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &sym,
        22,
        BiomeZone::Alpha,
        1.1,
        0.35,
        0.45,
        &ctrl,
        0.92,
        14,
        0.4,
    );
    let init_host_len = holo.sequence.len();
    let init_org_len = holo.endosymbiont_sequence.len();
    let init_work = holo.harvested_work;

    let mut migrations = 0usize;
    while holo.migrate_endosymbiont_locus_to_host(false, &mut rng) {
        migrations += 1;
    }
    assert_eq!(holo.endosymbiont_sequence.len(), CORR_CORE_LOCI);
    assert_eq!(migrations, init_org_len - CORR_CORE_LOCI);
    assert_eq!(holo.sequence.len(), init_host_len + migrations);
    assert!(holo.is_organelle_enslaved());

    holo.recompute_phenome(&ctrl, 1.1, 0.35, 0.45, 0.35, BiomeZone::Alpha, 0.92, 14, 0.4);
    let corr_work = holo.harvested_work;
    assert!(corr_work > init_work);

    while holo.migrate_endosymbiont_locus_to_host(true, &mut rng) {}
    assert_eq!(holo.endosymbiont_sequence.len(), 0);
    holo.recompute_phenome(&ctrl, 1.1, 0.35, 0.45, 0.35, BiomeZone::Alpha, 0.92, 14, 0.4);
    assert!(holo.harvested_work < corr_work);

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(0, &[(1, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5)]);
    tracker.record_birth(
        2,
        1,
        None,
        4,
        0,
        TrophicGuild::Syntroph,
        BiomeZone::Beta,
        0,
        TrophicGuild::Producer,
        0.30,
        0.25,
        1.6,
    );
    tracker.record_endosymbiotic_fusion(3, 1, 2, 10, 0, BiomeZone::Alpha, gamma, 10, 2.1);
    tracker.record_birth_with_egt(
        4,
        3,
        None,
        18,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        0,
        TrophicGuild::Holobiont,
        0.08,
        0.25,
        2.5,
        7,
        3,
    );
    tracker.sync_extant_population_with_egt(
        24,
        &[(4, 0, TrophicGuild::Holobiont, BiomeZone::Alpha, 2.5, 7, 3)],
    );

    let nhx = tracker.to_newick(24, true);
    assert!(nhx.contains("EGT=7:ORG_LEN=3:NUC_INT=0.70:LOCK=1"));
    let svg = tracker.to_svg_phylogram("EGT Lock-In Test", 24);
    assert!(svg.contains("EGT→Nuc"));
}

#[test]
fn serial_secondary_and_tertiary_endosymbiosis_compacts_and_eliminates_nucleomorph_at_housekeeping_threshold() {
    let mut rng = Rng::from_seed(20261118);
    let ctrl = RegimeId::R35Symbiogenetic.preset_controls();

    let prod = ModeledUnit::new_founder_with_role(
        30,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let synt = ModeledUnit::new_founder_with_role(
        31,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let (mut holo_1deg, g1) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &prod,
        &synt,
        32,
        BiomeZone::Alpha,
        1.1,
        0.35,
        0.45,
        &ctrl,
        0.92,
        16,
        0.4,
    );
    while holo_1deg.migrate_endosymbiont_locus_to_host(false, &mut rng) {}

    let outer_euk = ModeledUnit::new_founder_with_role(
        33,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::PredatorParasite,
        5,
        5,
        BiomeZone::Beta,
        0.4,
        &mut rng,
    );
    let (mut holo_2deg, g2) = ModeledUnit::fuse_serial_holobiont(
        &outer_euk,
        &holo_1deg,
        34,
        BiomeZone::Beta,
        1.15,
        0.38,
        0.48,
        &ctrl,
        0.92,
        16,
        0.4,
    );
    assert_eq!(holo_2deg.endo_tier, 2);
    assert_eq!(holo_2deg.organelle_membranes, 4);
    let q_with_full_nm = holo_2deg.landauer_cost;

    loop {
        let (moved, swept) = holo_2deg.migrate_nucleomorph_locus_to_host(false, &mut rng);
        assert!(!swept);
        if moved == 0 {
            break;
        }
    }
    assert_eq!(
        holo_2deg.nucleomorph_sequence.len(),
        NUCLEOMORPH_HOUSEKEEPING_LOCI + 1
    );

    let (moved_last, swept_now) = holo_2deg.migrate_nucleomorph_locus_to_host(true, &mut rng);
    assert_eq!(moved_last, 1);
    assert!(swept_now);
    assert_eq!(holo_2deg.nucleomorph_sequence.len(), 0);
    holo_2deg.recompute_phenome(&ctrl, 1.15, 0.38, 0.48, 0.35, BiomeZone::Beta, 0.92, 16, 0.4);
    assert!(holo_2deg.landauer_cost < q_with_full_nm);

    let (holo_3deg, g3) = ModeledUnit::fuse_serial_holobiont(
        &holo_1deg,
        &holo_2deg,
        35,
        BiomeZone::Alpha,
        1.18,
        0.40,
        0.50,
        &ctrl,
        0.92,
        16,
        0.4,
    );
    assert_eq!(holo_3deg.endo_tier, 3);

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (1, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (2, 1, TrophicGuild::Syntroph, BiomeZone::Beta, 1.6),
            (3, 1, TrophicGuild::PredatorParasite, BiomeZone::Beta, 1.4),
        ],
    );
    tracker.record_endosymbiotic_fusion(4, 1, 2, 8, 0, BiomeZone::Alpha, g1, 3, 2.2);
    tracker.record_serial_endosymbiotic_fusion(5, 3, 4, 16, 1, BiomeZone::Beta, 2, g2, 3, 8, 7, 2.6);
    tracker.record_birth_full(
        6,
        5,
        None,
        24,
        1,
        TrophicGuild::Holobiont,
        BiomeZone::Beta,
        1,
        TrophicGuild::Holobiont,
        0.09,
        0.25,
        2.9,
        7,
        3,
        2,
        0,
        true,
    );
    tracker.record_serial_endosymbiotic_fusion(7, 4, 6, 30, 0, BiomeZone::Alpha, 3, g3, 3, 8, 7, 3.1);
    tracker.sync_extant_population_full(
        34,
        &[
            ExtantPhyloState {
                taxon_id: 6,
                species_id: 1,
                guild: TrophicGuild::Holobiont,
                biome: BiomeZone::Beta,
                leverage_li: 2.9,
                egt_loci: 7,
                organelle_rem: 3,
                endo_tier: 2,
                nucleomorph_loci: 0,
                nucleomorph_eliminated: true,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                apicoplast_eliminated: false,
                mro_class: MroClass::AerobicMitochondrion,
                mro_pathways: 4,
                mro_mtdna_loci: 3,
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
            },
            ExtantPhyloState {
                taxon_id: 7,
                species_id: 0,
                guild: TrophicGuild::Holobiont,
                biome: BiomeZone::Alpha,
                leverage_li: 3.1,
                egt_loci: 7,
                organelle_rem: 3,
                endo_tier: 3,
                nucleomorph_loci: 8,
                nucleomorph_eliminated: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: false,
                apicoplast_pathways: 4,
                apicoplast_eliminated: false,
                mro_class: MroClass::AerobicMitochondrion,
                mro_pathways: 4,
                mro_mtdna_loci: 3,
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
            },
        ],
    );

    let nhx = tracker.to_newick(34, true);
    assert!(nhx.contains("TIER=2") && nhx.contains("TIER=3") && nhx.contains("NM_LOST=1"));
    let svg = tracker.to_svg_phylogram("Serial Endosymbiosis Test", 34);
    assert!(svg.contains("NM→0") && svg.contains("#H3"));
}

#[test]
fn mixotrophic_kleptoplasty_decays_without_hgt_and_crosses_tipping_point_at_threshold() {
    let mut rng = Rng::from_seed(20261124);
    let ctrl = RegimeId::R35Symbiogenetic.preset_controls();

    let mut prey_alpha = ModeledUnit::new_founder_with_role(
        40,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    prey_alpha.species_id = 0;

    let mut prey_beta = ModeledUnit::new_founder_with_role(
        41,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        5,
        5,
        BiomeZone::Beta,
        0.4,
        &mut rng,
    );
    prey_beta.species_id = 1;

    let mut predator = ModeledUnit::new_founder_with_role(
        42,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::PredatorParasite,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    predator.recompute_phenome(&ctrl, 1.0, 0.35, 0.45, 0.25, BiomeZone::Alpha, 0.92, 16, 0.4);
    let pure_pred_work = predator.harvested_work;

    let (stole, hgt0, crossed0) =
        predator.steal_kleptoplast_from_prey(&prey_alpha, true, false, true, &mut rng);
    assert!(stole && !hgt0 && !crossed0);
    predator.recompute_phenome(&ctrl, 1.0, 0.35, 0.45, 0.25, BiomeZone::Alpha, 0.92, 16, 0.4);
    assert!(predator.harvested_work > pure_pred_work);

    for _ in 0..7 {
        predator.step_lifetime(
            &ctrl,
            1.0,
            0.35,
            0.45,
            0.25,
            BiomeZone::Alpha,
            0.92,
            16,
            0.4,
            &mut rng,
        );
    }
    assert!(!predator.has_active_kleptoplast());

    for step_i in 0..KLEPTO_FIXATION_HGT_THRESHOLD {
        let donor = if step_i % 2 == 0 {
            &prey_alpha
        } else {
            &prey_beta
        };
        let (_s, did_hgt, crossed) =
            predator.steal_kleptoplast_from_prey(donor, true, true, true, &mut rng);
        assert!(did_hgt);
        if step_i + 1 < KLEPTO_FIXATION_HGT_THRESHOLD {
            assert!(!crossed);
        } else {
            assert!(crossed);
        }
    }
    assert_eq!(predator.guild, TrophicGuild::Holobiont);
    assert!(predator.mosaic_donor_count() >= 2);

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (40, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (41, 1, TrophicGuild::Syntroph, BiomeZone::Beta, 1.6),
            (42, 0, TrophicGuild::PredatorParasite, BiomeZone::Alpha, 1.4),
        ],
    );
    tracker.record_klepto_tipping_fixation(
        43,
        42,
        Some(41),
        18,
        0,
        BiomeZone::Alpha,
        1,
        0.28,
        KLEPTO_FIXATION_HGT_THRESHOLD,
        predator.endosymbiont_sequence.len(),
        0,
        2.7,
    );
    tracker.sync_extant_population_full(
        26,
        &[ExtantPhyloState {
            taxon_id: 43,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 2.7,
            egt_loci: KLEPTO_FIXATION_HGT_THRESHOLD,
            organelle_rem: predator.endosymbiont_sequence.len(),
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: KLEPTO_FIXATION_HGT_THRESHOLD,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
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
        }],
    );

    let nhx = tracker.to_newick(26, true);
    assert!(nhx.contains("FUSION=Klepto_Tipping_Fixation") && nhx.contains("KLEPTO_TIP=1"));
    let svg = tracker.to_svg_phylogram("Kleptoplasty Tipping Point Test", 26);
    assert!(svg.contains("KP→FIX") && svg.contains("#H1"));
}

#[test]
fn reductive_apicoplast_evolution_retains_essential_mep_suf_core_unless_ipp_rescued_and_preserves_algal_reticulation() {
    let mut rng = Rng::from_seed(20261202);
    let mut ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    ctrl.host_metabolite_rescue = 0.15;

    let prod = ModeledUnit::new_founder_with_role(
        50,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let synt = ModeledUnit::new_founder_with_role(
        51,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let (mut red_alga_1deg, g1) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &prod,
        &synt,
        52,
        BiomeZone::Alpha,
        1.1,
        0.35,
        0.45,
        &ctrl,
        0.92,
        16,
        0.4,
    );
    while red_alga_1deg.migrate_endosymbiont_locus_to_host(false, &mut rng) {}

    let alveolate_host = ModeledUnit::new_founder_with_role(
        53,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::PredatorParasite,
        5,
        5,
        BiomeZone::Beta,
        0.4,
        &mut rng,
    );
    let (mut chromera_2deg, g2) = ModeledUnit::fuse_serial_holobiont(
        &alveolate_host,
        &red_alga_1deg,
        54,
        BiomeZone::Beta,
        1.15,
        0.38,
        0.48,
        &ctrl,
        0.92,
        16,
        0.4,
    );
    let photo_q = chromera_2deg.landauer_cost;
    let photo_primary_work = chromera_2deg.primary_work;
    assert_eq!(chromera_2deg.guild, TrophicGuild::Holobiont);
    assert_eq!(chromera_2deg.apicoplast_pathways_mask, APICO_PATH_ALL);

    assert!(chromera_2deg.convert_holobiont_to_apicoplast_parasite(&mut rng));
    chromera_2deg.recompute_phenome(&ctrl, 1.15, 0.38, 0.48, 0.30, BiomeZone::Beta, 0.92, 16, 0.4);
    assert_eq!(chromera_2deg.guild, TrophicGuild::PredatorParasite);
    assert!(chromera_2deg.has_apicoplast());
    assert_eq!(chromera_2deg.organelle_membranes, 4);
    assert_eq!(chromera_2deg.endosymbiont_sequence.len(), APICOPLAST_MIN_LOCI);
    assert_eq!(chromera_2deg.nucleomorph_sequence.len(), 0);
    assert!(
        chromera_2deg.landauer_cost < photo_q && chromera_2deg.primary_work < photo_primary_work * 0.35
    );

    let (shed_heme, lost0, dd0) = chromera_2deg.streamline_apicoplast_pathway(0.15, false);
    assert!(shed_heme && !lost0 && !dd0 && chromera_2deg.apicoplast_pathway_count() == 3);
    let (shed_fas2, lost1, dd1) = chromera_2deg.streamline_apicoplast_pathway(0.15, false);
    assert!(shed_fas2 && !lost1 && !dd1 && chromera_2deg.apicoplast_pathway_count() == 2);
    assert_eq!(chromera_2deg.apicoplast_pathways_mask, APICO_PATH_ESSENTIAL_CORE);

    let (shed_core, _, _) = chromera_2deg.streamline_apicoplast_pathway(0.15, false);
    assert!(!shed_core && chromera_2deg.has_apicoplast());

    let mut unrescued_clone = chromera_2deg.clone();
    let (_, _, delayed_death) = unrescued_clone.streamline_apicoplast_pathway(0.15, true);
    assert!(delayed_death && unrescued_clone.integrity < 0.30);

    let q_with_apico = chromera_2deg.landauer_cost;
    let mut ctrl_rescued = ctrl.clone();
    ctrl_rescued.host_metabolite_rescue = 0.80;
    let (shed_ipp_rescued, crypto_lost, dd_rescued) =
        chromera_2deg.streamline_apicoplast_pathway(0.80, false);
    assert!(shed_ipp_rescued && crypto_lost && !dd_rescued);
    assert!(!chromera_2deg.has_apicoplast() && chromera_2deg.apicoplast_eliminated);
    chromera_2deg.recompute_phenome(
        &ctrl_rescued,
        1.15,
        0.38,
        0.48,
        0.30,
        BiomeZone::Beta,
        0.92,
        16,
        0.4,
    );
    assert!(chromera_2deg.landauer_cost < q_with_apico && chromera_2deg.integrity > 0.90);

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (50, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (51, 1, TrophicGuild::Syntroph, BiomeZone::Beta, 1.6),
            (53, 1, TrophicGuild::PredatorParasite, BiomeZone::Beta, 1.4),
        ],
    );
    tracker.record_endosymbiotic_fusion(52, 50, 51, 8, 0, BiomeZone::Alpha, g1, 3, 2.2);
    tracker.record_serial_endosymbiotic_fusion(54, 53, 52, 16, 1, BiomeZone::Beta, 2, g2, 3, 8, 7, 2.7);
    tracker.record_apicoplast_transition(55, 54, 24, 1, BiomeZone::Beta, 2, 7, 1, 2, 2.9);
    tracker.record_birth_full(
        56,
        55,
        None,
        32,
        1,
        TrophicGuild::PredatorParasite,
        BiomeZone::Beta,
        1,
        TrophicGuild::PredatorParasite,
        0.08,
        0.25,
        3.1,
        7,
        0,
        2,
        0,
        true,
    );
    tracker.sync_extant_population_full(
        38,
        &[
            ExtantPhyloState {
                taxon_id: 55,
                species_id: 1,
                guild: TrophicGuild::PredatorParasite,
                biome: BiomeZone::Beta,
                leverage_li: 2.9,
                egt_loci: 7,
                organelle_rem: 1,
                endo_tier: 2,
                nucleomorph_loci: 0,
                nucleomorph_eliminated: true,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: true,
                apicoplast_pathways: 2,
                apicoplast_eliminated: false,
                mro_class: MroClass::AerobicMitochondrion,
                mro_pathways: 4,
                mro_mtdna_loci: 3,
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
            },
            ExtantPhyloState {
                taxon_id: 56,
                species_id: 1,
                guild: TrophicGuild::PredatorParasite,
                biome: BiomeZone::Beta,
                leverage_li: 3.1,
                egt_loci: 7,
                organelle_rem: 0,
                endo_tier: 2,
                nucleomorph_loci: 0,
                nucleomorph_eliminated: true,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                apicoplast_eliminated: true,
                mro_class: MroClass::Mitosome,
                mro_pathways: 1,
                mro_mtdna_loci: 0,
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
            },
        ],
    );

    assert!(tracker.nodes.contains_key(&52) && tracker.nodes.contains_key(&51));
    let sum = tracker.compute_summary(38);
    assert_eq!(sum.apicoplast_transition_events, 1);
    assert!(sum.apicoplast_loss_events >= 1);

    let nhx = tracker.to_newick(38, true);
    assert!(nhx.contains("APICO=1:APICO_PATHS=2") && nhx.contains("APICO_LOST=1") && nhx.contains("#H2"));
    let svg = tracker.to_svg_phylogram("Apicoplast Reductive Evolution Test", 38);
    assert!(svg.contains("Holo→Apico") && svg.contains("Apico→0") && svg.contains("#H2"));
}

#[test]
fn anaerobic_mro_reduction_streamlines_mitochondria_to_hydrogenosomes_and_mitosomes_retaining_isc_core_unless_cytosolic_suf_rescued() {
    let mut rng = Rng::from_seed(20261212);
    let mut ctrl_anoxic = RegimeId::R35Symbiogenetic.preset_controls();
    ctrl_anoxic.ambient_oxygen_level = 0.08;
    ctrl_anoxic.cytosolic_suf_hgt_rescue = 0.15;

    // 1. Fuse a primary aerobic mitochondriate Holobiont (`AerobicMitochondrion`, `OXPHOS + mtDNA=3`)
    let host = ModeledUnit::new_founder_with_role(
        60,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Beta,
        0.4,
        &mut rng,
    );
    let alpha_proteo = ModeledUnit::new_founder_with_role(
        61,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Beta,
        0.4,
        &mut rng,
    );
    let (mut euk, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &alpha_proteo,
        62,
        BiomeZone::Beta,
        1.1,
        0.35,
        0.45,
        &ctrl_anoxic,
        0.92,
        16,
        0.4,
    );
    assert_eq!(euk.mro_class, MroClass::AerobicMitochondrion);
    assert_eq!(euk.mro_pathways_mask, MRO_PATH_ALL);
    assert_eq!(euk.mro_mtdna_loci, CORR_CORE_LOCI);
    let aero_anoxic_work = euk.harvested_work;
    let aero_q = euk.landauer_cost;
    let aero_byp = euk.byproduct_output;

    // 2. Step 1 -> Class 2/3 Anaerobic Mitochondrion (*Blastocystis* intermediate)
    let (tr1, cls1, col1) = euk.step_anaerobic_mro_reduction(0.15, false);
    assert!(tr1 && cls1 == MroClass::AnaerobicMitochondrion && !col1);
    assert_eq!(euk.mro_mtdna_loci, 2);

    // 3. Step 2 -> Class 4 Hydrogenosome (*Trichomonas*): relinquishes ETC -> sweeps mtDNA to 0 (CoRR!),
    //    boosts H2/acetate syntrophic byproduct output, and outperforms aerobic OXPHOS under anoxia!
    let (tr2, cls2, col2) = euk.step_anaerobic_mro_reduction(0.15, false);
    assert!(tr2 && cls2 == MroClass::Hydrogenosome && !col2);
    assert_eq!(euk.mro_mtdna_loci, 0, "Hydrogenosomes must relinquish mtDNA once membrane ETC is lost (CoRR)!");
    assert_eq!(euk.mro_pathways_mask, MRO_PATH_HYDROGENOSOME);
    euk.recompute_phenome(&ctrl_anoxic, 1.1, 0.35, 0.45, 0.30, BiomeZone::Beta, 0.92, 16, 0.4);
    let hydro_q = euk.landauer_cost;
    assert!(euk.harvested_work > aero_anoxic_work, "Hydrogenosome must outperform stalled aerobic OXPHOS in anoxia!");
    assert!(euk.byproduct_output > aero_byp * 1.8, "Hydrogenosome PFO/HydA must boost H2 syntrophic cross-feeding efflux!");
    assert!(hydro_q < aero_q, "Eliminating mtDNA + OXPHOS must reduce Landauer maintenance Q!");

    // 4. Step 3 -> Class 5 Mitosome (*Giardia* / microsporidia): sheds PFO/HydA & organelle ATP,
    //    strictly arresting at `MRO_PATH_MITOSOME_CORE` (`ISC` Fe-S cluster assembly only)!
    let (tr3, cls3, col3) = euk.step_anaerobic_mro_reduction(0.15, false);
    assert!(tr3 && cls3 == MroClass::Mitosome && !col3);
    assert_eq!(euk.mro_mtdna_loci, 0);
    assert_eq!(euk.mro_pathways_mask, MRO_PATH_MITOSOME_CORE);
    euk.recompute_phenome(&ctrl_anoxic, 1.1, 0.35, 0.45, 0.30, BiomeZone::Beta, 0.92, 16, 0.4);
    let mitos_q = euk.landauer_cost;
    assert!(mitos_q < hydro_q, "Mitosome must further cut Landauer maintenance Q below Hydrogenosome!");

    // Verify purifying selection strictly arrests reduction at Mitosome `ISC` when cytosolic SUF < 0.65!
    let (tr_blocked, cls_stay, col_stay) = euk.step_anaerobic_mro_reduction(0.15, false);
    assert!(!tr_blocked && cls_stay == MroClass::Mitosome && !col_stay);

    // 5. Verify CIA `[Fe-S]` Nuclear Replicative Collapse when `ISC` is disrupted WITHOUT cytosolic SUF
    //    vs. viable *Monocercomonoides* `MRO -> 0` complete loss WITH cytosolic `SUF` (`>= 0.65`)!
    let mut isc_ablated = euk.clone();
    let (_, _, cia_collapse) = isc_ablated.step_anaerobic_mro_reduction(0.15, true);
    assert!(cia_collapse && isc_ablated.integrity < 0.25);

    let mut ctrl_suf = ctrl_anoxic.clone();
    ctrl_suf.cytosolic_suf_hgt_rescue = 0.80;
    let (tr_suf, cls_amito, col_suf) = euk.step_anaerobic_mro_reduction(0.80, false);
    assert!(tr_suf && cls_amito == MroClass::Amiteukaryote && !col_suf);
    assert_eq!(euk.mro_pathways_mask, 0);
    euk.recompute_phenome(&ctrl_suf, 1.1, 0.35, 0.45, 0.30, BiomeZone::Beta, 0.92, 16, 0.4);
    assert!(
        euk.landauer_cost < mitos_q && euk.integrity > 0.90,
        "Cytosolic SUF HGT replacement must unlock viable Monocercomonoides MRO->0 complete loss!"
    );

    // 6. Verify `PhylogenyTracker` preserves ancestral `#H1` alpha-proteobacterial donor reticulation
    //    across `Mito -> Hydro`, `Hydro -> MitoS`, and `MRO -> 0` transitions!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (60, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (61, 1, TrophicGuild::Syntroph, BiomeZone::Beta, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(62, 60, 61, 8, 0, BiomeZone::Beta, gamma, 3, 2.2);
    tracker.record_mro_transition(
        63,
        62,
        16,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Beta,
        1,
        MroClass::Hydrogenosome,
        3,
        0,
        7,
        0,
        2.6,
    );
    tracker.record_mro_transition(
        64,
        63,
        24,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Beta,
        1,
        MroClass::Mitosome,
        1,
        0,
        7,
        0,
        2.9,
    );
    tracker.record_mro_transition(
        65,
        64,
        32,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Beta,
        1,
        MroClass::Amiteukaryote,
        0,
        0,
        7,
        0,
        3.2,
    );
    tracker.sync_extant_population_full(
        38,
        &[
            ExtantPhyloState {
                taxon_id: 64,
                species_id: 0,
                guild: TrophicGuild::Holobiont,
                biome: BiomeZone::Beta,
                leverage_li: 2.9,
                egt_loci: 7,
                organelle_rem: 0,
                endo_tier: 1,
                nucleomorph_loci: 0,
                nucleomorph_eliminated: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                apicoplast_eliminated: false,
                mro_class: MroClass::Mitosome,
                mro_pathways: 1,
                mro_mtdna_loci: 0,
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
            },
            ExtantPhyloState {
                taxon_id: 65,
                species_id: 0,
                guild: TrophicGuild::Holobiont,
                biome: BiomeZone::Beta,
                leverage_li: 3.2,
                egt_loci: 7,
                organelle_rem: 0,
                endo_tier: 1,
                nucleomorph_loci: 0,
                nucleomorph_eliminated: false,
                klepto_hgt_loci: 0,
                has_kleptoplast: false,
                has_apicoplast: false,
                apicoplast_pathways: 0,
                apicoplast_eliminated: false,
                mro_class: MroClass::Amiteukaryote,
                mro_pathways: 0,
                mro_mtdna_loci: 0,
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
            },
        ],
    );

    assert!(
        tracker.nodes.contains_key(&61),
        "Ancestral alpha-proteobacterial donor (#61) must be preserved in DAG even after mtDNA=0 and MRO->0!"
    );
    let sum = tracker.compute_summary(38);
    assert_eq!(sum.hydrogenosome_transition_events, 1);
    assert_eq!(sum.mitosome_transition_events, 1);
    assert_eq!(sum.mro_complete_loss_events, 1);

    let nhx = tracker.to_newick(38, true);
    assert!(nhx.contains("MRO=MitoS:MRO_PATHS=1:MTDNA=0") && nhx.contains("MRO=Amiteuk") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Anaerobic MRO Reductive Evolution Test", 38);
    assert!(
        svg.contains("Mito→Hydro")
            && svg.contains("Hydro→MitoS")
            && svg.contains("MRO→0")
            && svg.contains("#H1")
    );
}

#[test]
fn clonal_multicellularity_evolves_germ_soma_division_shields_germline_fidelity_and_polices_cheaters_via_mitochondrial_apoptosis() {
    let mut rng = Rng::from_seed(20261218);
    let mut ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.65;
    ctrl.germ_soma_differentiation_rate = 0.65;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;

    let host = ModeledUnit::new_founder_with_role(
        70,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        5,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let alpha_sym = ModeledUnit::new_founder_with_role(
        71,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        5,
        6,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let predator = ModeledUnit::new_founder_with_role(
        72,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::PredatorParasite,
        6,
        6,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );

    let (mut holo, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &alpha_sym,
        73,
        BiomeZone::Alpha,
        1.15,
        0.35,
        0.45,
        &ctrl,
        0.92,
        18,
        0.4,
    );
    assert_eq!(holo.multicell_state, MulticellularState::Unicellular);
    let uni_work = holo.harvested_work;
    let uni_susc = holo.susceptibility_to(&predator);

    // 1. Step 1 -> Undifferentiated Clonal Cluster ("staying together", N=12, r=1.0)
    let (tr1, mc1) = holo.step_clonal_multicellularity(true, &mut rng);
    assert!(tr1 && mc1 == MulticellularState::UndifferentiatedClonal);
    assert_eq!(holo.kin_relatedness_r, 1.0);
    holo.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 18, 0.4);
    let clonal_undiff = holo.clone();

    // 2. Step 2 -> Clonal Germ-Soma Organism (Volvox / metazoan grade: soma_fraction = 0.75, N=32)
    let (tr2, mc2) = holo.step_clonal_multicellularity(true, &mut rng);
    assert!(tr2 && mc2 == MulticellularState::GermSomaMulticellular);
    assert!(holo.soma_fraction >= 0.70);
    holo.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 18, 0.4);
    assert!(
        holo.harvested_work > uni_work && holo.susceptibility_to(&predator) < uni_susc * 0.70,
        "Germ-soma division must boost cooperative foraging work and grant a size-refuge predation shield!"
    );

    // 3. Verify dirty-work germline fidelity shield (Radzvilavicius et al., 2016) & N=1 zygotic bottleneck!
    let mut f_undiff_sum = 0.0;
    let mut f_gs_sum = 0.0;
    for k in 0..24 {
        let (_, f_u, _, _) = clonal_undiff.replicate_with_partner(
            None,
            200 + k,
            5,
            5,
            BiomeZone::Alpha,
            1.15,
            0.35,
            0.45,
            &ctrl,
            0.92,
            18,
            0.4,
            &mut rng,
        );
        let (gs_child, f_g, _, _) = holo.replicate_with_partner(
            None,
            300 + k,
            5,
            5,
            BiomeZone::Alpha,
            1.15,
            0.35,
            0.45,
            &ctrl,
            0.92,
            18,
            0.4,
            &mut rng,
        );
        assert_eq!(gs_child.somatic_cheater_load, 0.0, "Zygotic N=1 bottleneck must reset somatic cheaters!");
        f_undiff_sum += f_u;
        f_gs_sum += f_g;
    }
    assert!(
        f_gs_sum > f_undiff_sum,
        "Sequestered quiescent germline must out-proofread totipotent OXPHOS dirty-work clusters!"
    );

    // 4. Verify mitochondrial apoptosis suppresses somatic cheaters vs. neoplastic collapse when apoptosis & bottleneck fail!
    let mut policed = holo.clone();
    let mut unpoliced_chimera = holo.clone();
    unpoliced_chimera.kin_relatedness_r = 0.48;
    let mut collapsed = false;
    for _ in 0..16 {
        let (_, col_ok) = policed.step_multicellular_homeostasis(0.85, true);
        assert!(!col_ok);
        let (_, col_bad) = unpoliced_chimera.step_multicellular_homeostasis(0.0, false);
        if col_bad {
            collapsed = true;
        }
    }
    assert!(policed.cumulative_apoptotic_culls > 0 && policed.somatic_cheater_load < 0.15);
    assert!(
        collapsed && unpoliced_chimera.integrity < 0.25,
        "Ablating mitochondrial apoptosis without a single-cell bottleneck must trigger somatic neoplastic collapse!"
    );

    // 5. Verify PhylogenyTracker preserves #H1 mitochondrial donor arc across Holo->Clonal and Clonal->GermSoma!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (70, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (71, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(73, 70, 71, 6, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_multicellular_transition(
        74,
        73,
        14,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::UndifferentiatedClonal,
        12,
        0.0,
        5,
        3,
        2.7,
    );
    tracker.record_multicellular_transition(
        75,
        74,
        22,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        6,
        3,
        3.5,
    );
    tracker.sync_extant_population_full(
        30,
        &[ExtantPhyloState {
            taxon_id: 75,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 3.5,
            egt_loci: 6,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
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
        }],
    );

    assert!(tracker.nodes.contains_key(&71), "Endosymbiotic mitochondrial donor #71 must remain preserved!");
    let sum = tracker.compute_summary(30);
    assert_eq!(sum.clonal_multicell_transitions, 1);
    assert_eq!(sum.germ_soma_transitions, 1);

    let nhx = tracker.to_newick(30, true);
    assert!(nhx.contains("MC=GermSoma:CELLS=32:SOMA=0.75") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Clonal Multicellularity & Germ-Soma Test", 30);
    assert!(svg.contains("Holo→Clonal") && svg.contains("Clonal→GermSoma") && svg.contains("#H1"));
}

#[test]
fn eusocial_superorganisms_require_monogamy_window_evolve_sterile_worker_soldier_castes_and_enforce_kin_policing() {
    let mut rng = Rng::from_seed(20261228);
    let mut ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.70;
    ctrl.germ_soma_differentiation_rate = 0.70;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;
    ctrl.queen_monogamy_index = 0.90;
    ctrl.eusocial_caste_canalization_rate = 0.75;
    ctrl.worker_policing_rate = 0.85;

    let host = ModeledUnit::new_founder_with_role(
        80,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        5,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let alpha_sym = ModeledUnit::new_founder_with_role(
        81,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        5,
        6,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let predator = ModeledUnit::new_founder_with_role(
        82,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::PredatorParasite,
        6,
        6,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );

    let (mut org, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &alpha_sym,
        83,
        BiomeZone::Alpha,
        1.15,
        0.35,
        0.45,
        &ctrl,
        0.92,
        20,
        0.4,
    );
    org.step_clonal_multicellularity(true, &mut rng);
    org.step_clonal_multicellularity(true, &mut rng);
    assert_eq!(org.multicell_state, MulticellularState::GermSomaMulticellular);
    org.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 20, 0.4);
    let solitary_gs_work = org.harvested_work;
    let solitary_gs_niche = org.niche_output;
    let solitary_gs_susc = org.susceptibility_to(&predator);

    // 1. Advance to SubsocialFamily -> FacultativeEusocial (*Polistes* wasp grade)
    let (e1, s1) = org.step_eusocial_superorganism(0.60, &mut rng);
    assert!(e1 && s1 == EusocialState::SubsocialFamily);
    let (e2, s2) = org.step_eusocial_superorganism(0.60, &mut rng);
    assert!(e2 && s2 == EusocialState::FacultativeEusocial);

    // 2. Verify Boomsma's (2009) Lifetime Monogamy Window: promiscuity (`queen_monogamy = 0.60 < 0.75`)
    //    strictly blocks crossing the morphological Point of No Return into `ObligateSuperorganism`!
    let (blocked, s_stay) = org.step_eusocial_superorganism(0.60, &mut rng);
    assert!(!blocked && s_stay == EusocialState::FacultativeEusocial);

    // 3. Open the Lifetime Monogamy Window (`queen_monogamy = 0.92 >= 0.75`, r_sib ~ 0.48-0.50):
    //    crosses Point of No Return into `ObligateSuperorganism` with sterile worker + soldier castes!
    let (e3, s3) = org.step_eusocial_superorganism(0.92, &mut rng);
    assert!(e3 && s3 == EusocialState::ObligateSuperorganism);
    assert_eq!(org.colony_size, 64);
    assert!(org.sterile_worker_fraction >= 0.70 && org.soldier_fraction >= 0.15);
    org.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 20, 0.4);
    assert!(
        org.harvested_work > solitary_gs_work
            && org.niche_output > solitary_gs_niche
            && org.susceptibility_to(&predator) < solitary_gs_susc * 0.55,
        "Obligate superorganism sterile worker/soldier castes must boost foraging, nest-fortress niche construction, and predation defense!"
    );

    // 4. Verify secondary polyandry (*Apis mellifera* mode, `queen_monogamy = 0.30`) remains harmonious
    //    when mutual worker kin-policing is active (`worker_policing_rate = 0.85`), whereas ablating
    //    worker policing (`0.0`) unleashes Rebel Worker Reproductive Anarchy Collapse!
    let mut policed_polyandrous = org.clone();
    let mut unpoliced_polyandrous = org.clone();
    let mut anarchy_collapsed = false;
    for _ in 0..16 {
        let (_, col_ok) = policed_polyandrous.step_colony_kin_policing(0.30, 0.85);
        assert!(!col_ok);
        let (_, col_bad) = unpoliced_polyandrous.step_colony_kin_policing(0.30, 0.0);
        if col_bad {
            anarchy_collapsed = true;
        }
    }
    assert!(
        policed_polyandrous.cumulative_policed_eggs > 0
            && policed_polyandrous.rebel_worker_load < 0.14
    );
    assert!(
        anarchy_collapsed && unpoliced_polyandrous.integrity < 0.25,
        "Ablating worker kin-policing under polyandry must trigger Colony Reproductive Anarchy Collapse!"
    );

    // 5. Verify PhylogenyTracker preserves #H1 mitochondrial donor arc across the entire 4-level ETI chain:
    //    `#H1 -> Holo->Clonal -> Clonal->GermSoma -> GermSoma->Eusocial -> Eusocial->Superorg`!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (80, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (81, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(83, 80, 81, 5, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_multicellular_transition(
        84,
        83,
        10,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::UndifferentiatedClonal,
        12,
        0.0,
        5,
        3,
        2.7,
    );
    tracker.record_multicellular_transition(
        85,
        84,
        16,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        6,
        3,
        3.4,
    );
    tracker.record_superorganism_transition(
        86,
        85,
        22,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::FacultativeEusocial,
        18,
        0.68,
        6,
        3,
        3.9,
    );
    tracker.record_superorganism_transition(
        87,
        86,
        28,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::ObligateSuperorganism,
        64,
        0.90,
        7,
        3,
        4.6,
    );
    tracker.sync_extant_population_full(
        36,
        &[ExtantPhyloState {
            taxon_id: 87,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 4.6,
            egt_loci: 7,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::ObligateSuperorganism,
            colony_size: 64,
            sterile_caste_fraction: 0.90,
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
        }],
    );

    assert!(tracker.nodes.contains_key(&81), "Ancestral #H1 mitochondrial donor must remain preserved!");
    let sum = tracker.compute_summary(36);
    assert_eq!(sum.eusocial_transitions, 1);
    assert_eq!(sum.superorganism_transitions, 1);

    let nhx = tracker.to_newick(36, true);
    assert!(nhx.contains("EUSOC=Superorg:COLONY=64:CASTES=0.90") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Eusocial Superorganism ETI Test", 36);
    assert!(
        svg.contains("Holo→Clonal")
            && svg.contains("Clonal→GermSoma")
            && svg.contains("GermSoma→Eusocial")
            && svg.contains("Eusocial→Superorg")
            && svg.contains("#H1")
    );
}

#[test]
fn obligate_agricultural_superorganisms_domesticate_polyploid_clonal_crops_and_suppress_escovopsis_via_defensive_actinobacteria() {
    let mut rng = Rng::from_seed(20261229);
    let mut ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.70;
    ctrl.germ_soma_differentiation_rate = 0.70;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;
    ctrl.queen_monogamy_index = 0.92;
    ctrl.eusocial_caste_canalization_rate = 0.78;
    ctrl.worker_policing_rate = 0.85;
    ctrl.crop_domestication_rate = 0.80;
    ctrl.agro_pathogen_pressure = 0.50;
    ctrl.defensive_symbiont_potency = 0.85;

    let host = ModeledUnit::new_founder_with_role(
        90,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let alpha_sym = ModeledUnit::new_founder_with_role(
        91,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );

    let (mut org, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &alpha_sym,
        92,
        BiomeZone::Alpha,
        1.15,
        0.35,
        0.45,
        &ctrl,
        0.92,
        20,
        0.4,
    );
    org.step_clonal_multicellularity(true, &mut rng);
    org.step_clonal_multicellularity(true, &mut rng);

    // 1. Solitary organisms cannot initiate fungus farming prior to eusociality!
    let (solitary_blocked, ag0) = org.step_agricultural_domestication(0.85, &mut rng);
    assert!(!solitary_blocked && ag0 == AgricultureState::HunterForager);

    // 2. Advance to FacultativeEusocial -> can initiate LowerFacultativeFarming (2n),
    //    but crossing into HigherDomesticatedCrop (4n gongylidia) requires ObligateSuperorganism!
    org.step_eusocial_superorganism(0.65, &mut rng);
    org.step_eusocial_superorganism(0.65, &mut rng);
    assert_eq!(org.eusocial_state, EusocialState::FacultativeEusocial);
    let (ag_ok1, ag1) = org.step_agricultural_domestication(0.85, &mut rng);
    assert!(ag_ok1 && ag1 == AgricultureState::LowerFacultativeFarming && org.cultivar_polyploidy == 2);
    let (higher_blocked, ag_stay) = org.step_agricultural_domestication(0.85, &mut rng);
    assert!(!higher_blocked && ag_stay == AgricultureState::LowerFacultativeFarming);

    // 3. Cross Point of No Return into ObligateSuperorganism -> domesticate 4n HigherDomesticatedCrop
    //    (with host metabolic locus streamlining!) -> 8n LeafCutterAgroEcosystem!
    org.step_eusocial_superorganism(0.92, &mut rng);
    assert_eq!(org.eusocial_state, EusocialState::ObligateSuperorganism);

    let mut non_farming_superorg = org.clone();
    non_farming_superorg.agriculture_state = AgricultureState::HunterForager;
    non_farming_superorg.cultivar_polyploidy = 1;
    non_farming_superorg.gongylidia_yield = 0.0;
    non_farming_superorg.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 20, 0.4);
    let forager_superorg_work = non_farming_superorg.harvested_work;
    let forager_superorg_niche = non_farming_superorg.niche_output;

    let (ag_ok2, ag2) = org.step_agricultural_domestication(0.85, &mut rng);
    assert!(ag_ok2 && ag2 == AgricultureState::HigherDomesticatedCrop && org.cultivar_polyploidy == 4);
    let (ag_ok3, ag3) = org.step_agricultural_domestication(0.85, &mut rng);
    assert!(ag_ok3 && ag3 == AgricultureState::LeafCutterAgroEcosystem && org.cultivar_polyploidy == 8);
    assert!(org.crop_clonal_uniformity >= 0.95 && org.gongylidia_yield >= 0.90);

    org.recompute_phenome(&ctrl, 1.15, 0.35, 0.45, 0.30, BiomeZone::Alpha, 0.92, 20, 0.4);
    assert!(
        org.harvested_work > forager_superorg_work * 1.30
            && org.niche_output > forager_superorg_niche,
        "8n LeafCutter polyploid gongylidia crops must substantially elevate cellulose-to-energy work and fungus-garden niche capital!"
    );

    // 4. Verify Red-Queen Clonal Monoculture Susceptibility vs. Cuticular *Pseudonocardia* Defense:
    //    Protected clonal crop (`defensive_symbiont_potency = 0.85`) suppresses *Escovopsis* (`load < 0.05`),
    //    whereas ablating the actinobacterial symbiont (`0.0`) triggers Escovopsis Garden Blight Collapse!
    let mut shielded_farm = org.clone();
    let mut unshielded_monoculture = org.clone();
    let mut blight_collapsed = false;
    for _ in 0..10 {
        let (_, b_ok) = shielded_farm.step_agro_pathogen_and_defensive_microbiome(0.50, 0.85);
        assert!(!b_ok);
        let (_, b_bad) = unshielded_monoculture.step_agro_pathogen_and_defensive_microbiome(0.50, 0.0);
        if b_bad {
            blight_collapsed = true;
        }
    }
    assert!(
        shielded_farm.cumulative_antimicrobial_shields > 0
            && shielded_farm.escovopsis_pathogen_load < 0.05
    );
    assert!(
        blight_collapsed && unshielded_monoculture.integrity < 0.25,
        "Ablating cuticular Pseudonocardia in a clonal crop monoculture must trigger Escovopsis Garden Blight Collapse!"
    );

    // 5. Verify PhylogenyTracker preserves the complete 5-tier ETI + Ectosymbiotic Agriculture chain!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (90, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (91, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(92, 90, 91, 4, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_multicellular_transition(
        93,
        92,
        8,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::UndifferentiatedClonal,
        12,
        0.0,
        5,
        3,
        2.7,
    );
    tracker.record_multicellular_transition(
        94,
        93,
        12,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        6,
        3,
        3.4,
    );
    tracker.record_superorganism_transition(
        95,
        94,
        18,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::ObligateSuperorganism,
        64,
        0.90,
        7,
        3,
        4.5,
    );
    tracker.record_agriculture_transition(
        96,
        95,
        24,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::ObligateSuperorganism,
        64,
        0.90,
        AgricultureState::LowerFacultativeFarming,
        2,
        0.28,
        7,
        3,
        4.9,
    );
    tracker.record_agriculture_transition(
        97,
        96,
        30,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::ObligateSuperorganism,
        64,
        0.90,
        AgricultureState::LeafCutterAgroEcosystem,
        8,
        0.96,
        8,
        3,
        5.8,
    );
    tracker.sync_extant_population_full(
        38,
        &[ExtantPhyloState {
            taxon_id: 97,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 5.8,
            egt_loci: 8,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::ObligateSuperorganism,
            colony_size: 64,
            sterile_caste_fraction: 0.90,
            agriculture_state: AgricultureState::LeafCutterAgroEcosystem,
            cultivar_polyploidy: 8,
            gongylidia_yield: 0.96,
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
        }],
    );

    assert!(tracker.nodes.contains_key(&91), "Ancestral #H1 mitochondrial donor must remain preserved!");
    let sum = tracker.compute_summary(38);
    assert_eq!(sum.agriculture_transitions, 1);
    assert_eq!(sum.obligate_crop_transitions, 1);

    let nhx = tracker.to_newick(38, true);
    assert!(nhx.contains("AGRI=LeafCutter:PLOIDY=8n:GONGYL=0.96") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Obligate Agricultural Superorganism Test", 38);
    assert!(
        svg.contains("Holo→Clonal")
            && svg.contains("Clonal→GermSoma")
            && svg.contains("Eusocial→Superorg")
            && svg.contains("Superorg→LowerAg")
            && svg.contains("LowerAg→GongylCrop")
            && svg.contains("#H1")
    );
}

#[test]
fn cumulative_cultural_evolution_requires_teaching_fidelity_peaks_under_partial_deme_connectivity_and_collapses_under_isolation_drift() {
    let mut rng = Rng::from_seed(20261230);
    let mut ctrl = RegimeId::R30PseudoSoftInheritance.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.72;
    ctrl.germ_soma_differentiation_rate = 0.72;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;
    ctrl.symbolic_teaching_fidelity = 0.94;
    ctrl.inter_deme_connectivity = 0.36;
    ctrl.cultural_innovation_rate = 0.75;

    let host_a = ModeledUnit::new_founder_with_role(
        110,
        OrganismKind::HominidPopulation,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Producer,
        2,
        2,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let sym_a = ModeledUnit::new_founder_with_role(
        111,
        OrganismKind::HominidPopulation,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Syntroph,
        2,
        3,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let (mut deme0, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host_a,
        &sym_a,
        112,
        BiomeZone::Alpha,
        1.20,
        0.38,
        0.48,
        &ctrl,
        0.93,
        24,
        0.4,
    );
    deme0.deme_id = 0;
    deme0.step_clonal_multicellularity(true, &mut rng);
    deme0.step_clonal_multicellularity(true, &mut rng);

    let mut deme2 = deme0.clone();
    deme2.id = 113;
    deme2.deme_id = 2;

    let (c1, s1) = deme0.step_cumulative_culture_transition(0.68, 0.36, &mut rng);
    assert!(c1 && s1 == CulturalState::LowFidelityTradition);
    let (blocked_ratchet, s_stay) = deme0.step_cumulative_culture_transition(0.68, 0.36, &mut rng);
    assert!(!blocked_ratchet && s_stay == CulturalState::LowFidelityTradition);

    let (c2, s2) = deme0.step_cumulative_culture_transition(0.84, 0.36, &mut rng);
    assert!(c2 && s2 == CulturalState::ImitativeRatchet);
    let (c3, s3) = deme0.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    assert!(c3 && s3 == CulturalState::SymbolicCumulativeCulture && deme0.neural_learning_capacity >= 0.88);

    deme2.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    deme2.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    deme2.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    assert_eq!(deme2.cultural_state, CulturalState::SymbolicCumulativeCulture);

    let mut overconnected = deme0.clone();
    let added_over = overconnected.fuse_cross_deme_cultural_recombination(&deme2, 0.95, 0.94);
    assert_eq!(added_over, 0);

    let added_partial = deme0.fuse_cross_deme_cultural_recombination(&deme2, 0.36, 0.94);
    assert!(added_partial > 0 && deme0.cultural_ratchet_tier >= 5);

    let mut ctrl_over = ctrl.clone();
    ctrl_over.inter_deme_connectivity = 0.95;
    let mut ctrl_iso = ctrl.clone();
    ctrl_iso.inter_deme_connectivity = 0.02;

    let mut isolated = deme2.clone();
    deme0.recompute_phenome(&ctrl, 1.20, 0.38, 0.48, 0.30, BiomeZone::Alpha, 0.93, 24, 0.4);
    overconnected.recompute_phenome(&ctrl_over, 1.20, 0.38, 0.48, 0.30, BiomeZone::Alpha, 0.93, 24, 0.4);
    isolated.recompute_phenome(&ctrl_iso, 1.20, 0.38, 0.48, 0.30, BiomeZone::Alpha, 0.93, 24, 0.4);

    let li_partial = deme0.complexity_c / deme0.structure_s;
    let li_over = overconnected.complexity_c / overconnected.structure_s;
    let li_iso = isolated.complexity_c / isolated.structure_s;
    assert!(
        deme0.harvested_work > overconnected.harvested_work
            && deme0.harvested_work > isolated.harvested_work
            && li_partial > li_over
            && li_partial > li_iso,
        "Partially connected structured demes must outperform both over-connected conformist and isolated demes!"
    );

    let child = deme0.replicate_with_partner(
        Some(&deme2),
        120,
        &ctrl,
        1.20,
        0.38,
        0.48,
        BiomeZone::Alpha,
        0.93,
        24,
        0.4,
        &mut rng,
    );
    assert!(
        child.cultural_state == CulturalState::SymbolicCumulativeCulture
            && child.cultural_sequence.len() >= 14,
        "High-fidelity symbolic teaching must transmit cumulative cultural repertoires across generations!"
    );

    let mut tasmanian_isolated = deme0.clone();
    let mut drift_collapsed = false;
    for _ in 0..8 {
        let (_, collapsed) = tasmanian_isolated.step_dual_inheritance_and_cultural_drift(0.86, 0.0, 1);
        if collapsed {
            drift_collapsed = true;
        }
    }
    assert!(
        drift_collapsed && tasmanian_isolated.integrity < 0.28,
        "Severely isolated demes must suffer Tasmanian Demographic Cultural Ratchet Collapse!"
    );

    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (110, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (111, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(112, 110, 111, 4, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_cultural_transition(
        113,
        112,
        None,
        10,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::ImitativeRatchet,
        2,
        16.0,
        6,
        3,
        4.2,
    );
    tracker.record_cultural_transition(
        114,
        113,
        None,
        16,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        5,
        34.5,
        7,
        3,
        6.4,
    );
    tracker.sync_extant_population_full(
        24,
        &[ExtantPhyloState {
            taxon_id: 114,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 6.4,
            egt_loci: 7,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::Solitary,
            colony_size: 1,
            sterile_caste_fraction: 0.0,
            agriculture_state: AgricultureState::HunterForager,
            cultivar_polyploidy: 1,
            gongylidia_yield: 0.0,
            cultural_state: CulturalState::SymbolicCumulativeCulture,
            cultural_ratchet_tier: 5,
            cultural_repertoire_bits: 34.5,
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
        }],
    );

    assert!(tracker.nodes.contains_key(&111));
    let sum = tracker.compute_summary(24);
    assert_eq!(sum.cultural_ratchet_transitions, 1);
    assert_eq!(sum.symbolic_teaching_transitions, 1);

    let nhx = tracker.to_newick(24, true);
    assert!(nhx.contains("CULT=SymbolicCCE:RATCHET=T5:CBITS=34.5") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Cumulative Cultural Evolution Test", 24);
    assert!(
        svg.contains("Tradition→Ratchet")
            && svg.contains("Ratchet→SymbolicCCE")
            && svg.contains("#H1")
    );
}

#[test]
fn exosomatic_written_archives_decouple_cultural_storage_from_brain_size_rescue_isolated_demes_and_require_peer_verification_against_bad_science() {
    let mut rng = Rng::from_seed(20261231);
    let mut ctrl = RegimeId::R30PseudoSoftInheritance.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.72;
    ctrl.germ_soma_differentiation_rate = 0.72;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;
    ctrl.symbolic_teaching_fidelity = 0.94;
    ctrl.inter_deme_connectivity = 0.36;
    ctrl.cultural_innovation_rate = 0.75;
    ctrl.archival_externalization_rate = 0.80;
    ctrl.peer_verification_rigor = 0.88;
    ctrl.unverified_publish_pressure = 0.65;

    let host = ModeledUnit::new_founder_with_role(
        130,
        OrganismKind::HominidPopulation,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Producer,
        3,
        3,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let sym = ModeledUnit::new_founder_with_role(
        131,
        OrganismKind::HominidPopulation,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Syntroph,
        3,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let (mut brain_bound, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &sym,
        132,
        BiomeZone::Alpha,
        1.20,
        0.38,
        0.48,
        &ctrl,
        0.93,
        24,
        0.4,
    );
    brain_bound.step_clonal_multicellularity(true, &mut rng);
    brain_bound.step_clonal_multicellularity(true, &mut rng);
    brain_bound.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    brain_bound.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    brain_bound.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    assert_eq!(brain_bound.cultural_state, CulturalState::SymbolicCumulativeCulture);
    brain_bound.recompute_phenome(&ctrl, 1.20, 0.38, 0.48, 0.30, BiomeZone::Alpha, 0.93, 24, 0.4);

    let endosomatic_neural = brain_bound.neural_learning_capacity;
    let endosomatic_structure_s = brain_bound.structure_s;
    let endosomatic_landauer = brain_bound.landauer_cost;
    let endosomatic_li = brain_bound.complexity_c / brain_bound.structure_s;

    // 1. Step through MaterialNotation -> WrittenArchive -> PeerVerifiedInstitution:
    //    Verifies decoupling of cultural storage from individual brain size (Holocene encephalization relief)!
    let mut archival = brain_bound.clone();
    let (e1, x1) = archival.step_exosomatic_archive_transition(0.88, &mut rng);
    assert!(e1 && x1 == ExosomaticState::MaterialNotation);
    let (e2, x2) = archival.step_exosomatic_archive_transition(0.88, &mut rng);
    assert!(e2 && x2 == ExosomaticState::WrittenArchive);

    // Low peer-verification rigor (`0.35 < 0.65`) blocks institutionalization!
    let (blocked_inst, x_stay) = archival.step_exosomatic_archive_transition(0.35, &mut rng);
    assert!(!blocked_inst && x_stay == ExosomaticState::WrittenArchive);

    let (e3, x3) = archival.step_exosomatic_archive_transition(0.88, &mut rng);
    assert!(e3 && x3 == ExosomaticState::PeerVerifiedInstitution);
    archival.recompute_phenome(&ctrl, 1.20, 0.38, 0.48, 0.30, BiomeZone::Alpha, 0.93, 24, 0.4);

    let inst_li = archival.complexity_c / archival.structure_s;
    assert!(
        archival.exosomatic_archive_bits() > brain_bound.cultural_repertoire_bits() * 2.0
            && archival.neural_learning_capacity < endosomatic_neural
            && archival.structure_s < endosomatic_structure_s
            && archival.landauer_cost < endosomatic_landauer
            && inst_li > endosomatic_li * 1.35,
        "Exosomatic peer-verified archives must expand total cultural bits >2x while LOWERING individual neural brain mass, structure S, and Landauer maintenance cost!"
    );

    // 2. Verify Written Archives rescue small isolated demes (`peers = 1, connectivity = 0.0`) from Tasmanian collapse!
    let mut isolated_oral = brain_bound.clone();
    let mut isolated_literate = archival.clone();
    let mut oral_collapsed = false;
    for _ in 0..8 {
        let (_, c_oral) = isolated_oral.step_dual_inheritance_and_cultural_drift(0.86, 0.0, 1);
        if c_oral {
            oral_collapsed = true;
        }
        let (_, c_lit) = isolated_literate.step_dual_inheritance_and_cultural_drift(0.86, 0.0, 1);
        assert!(!c_lit, "Durable written archives must immunize isolated demes against Tasmanian demographic skill loss!");
    }
    assert!(oral_collapsed && isolated_literate.cultural_drift_erosion < 0.10);

    // 3. Smaldino & McElreath (2016) Natural Selection of Bad Science vs. Peer-Verification Ratchet:
    //    Ablating peer verification (`0.05`) under high publish pressure (`0.75`) floods the archive with
    //    false discoveries (`archive_false_discovery_load >= 0.70` -> Epistemic Archive Corruption Collapse)!
    let mut verified_inst = archival.clone();
    let mut unverified_archive = archival.clone();
    unverified_archive.exosomatic_state = ExosomaticState::WrittenArchive;
    let mut fdr_collapsed = false;
    for _ in 0..9 {
        let (_, bad_ok) = verified_inst.step_archive_verification_and_fdr_dynamics(0.75, 0.88);
        assert!(!bad_ok);
        let (_, bad_corrupt) = unverified_archive.step_archive_verification_and_fdr_dynamics(0.75, 0.05);
        if bad_corrupt {
            fdr_collapsed = true;
        }
    }
    assert!(
        verified_inst.cumulative_peer_verifications > 0
            && verified_inst.archive_false_discovery_load < 0.05
            && verified_inst.archive_sequence.len() > 38
    );
    assert!(
        fdr_collapsed && unverified_archive.integrity < 0.25,
        "Unverified archives under volume-driven publish pressure must suffer Epistemic Archive Corruption Collapse!"
    );

    // 4. Verify PhylogenyTracker preserves #H1 + SymbolicCCE + WrittenArch + VerifiedInst in eNewick+NHX & SVG!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (130, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (131, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(132, 130, 131, 4, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_cultural_transition(
        133,
        132,
        None,
        10,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::ImitativeRatchet,
        2,
        16.0,
        6,
        3,
        4.2,
    );
    tracker.record_cultural_transition(
        134,
        133,
        None,
        15,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        4,
        30.0,
        7,
        3,
        6.1,
    );
    tracker.record_exosomatic_transition(
        135,
        134,
        None,
        21,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        4,
        30.0,
        ExosomaticState::WrittenArchive,
        62.4,
        7,
        3,
        7.8,
    );
    tracker.record_exosomatic_transition(
        136,
        135,
        None,
        28,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        5,
        34.0,
        ExosomaticState::PeerVerifiedInstitution,
        118.6,
        7,
        3,
        10.2,
    );
    tracker.sync_extant_population_full(
        36,
        &[ExtantPhyloState {
            taxon_id: 136,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 10.2,
            egt_loci: 7,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::Solitary,
            colony_size: 1,
            sterile_caste_fraction: 0.0,
            agriculture_state: AgricultureState::HunterForager,
            cultivar_polyploidy: 1,
            gongylidia_yield: 0.0,
            cultural_state: CulturalState::SymbolicCumulativeCulture,
            cultural_ratchet_tier: 5,
            cultural_repertoire_bits: 34.0,
            exosomatic_state: ExosomaticState::PeerVerifiedInstitution,
            exosomatic_archive_bits: 118.6,
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
        }],
    );

    assert!(tracker.nodes.contains_key(&131));
    let sum = tracker.compute_summary(36);
    assert_eq!(sum.written_archive_transitions, 1);
    assert_eq!(sum.peer_institution_transitions, 1);

    let nhx = tracker.to_newick(36, true);
    assert!(nhx.contains("EXO=VerifiedInst:ABITS=118.6") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Exosomatic Written Archive & Peer-Verification Test", 36);
    assert!(
        svg.contains("Tradition→Ratchet")
            && svg.contains("Ratchet→SymbolicCCE")
            && svg.contains("Symbolic→WrittenArch")
            && svg.contains("WrittenArch→VerifiedInst")
            && svg.contains("#H1")
    );
}

#[test]
fn autonomous_machine_verified_synthesis_coevolves_with_archives_discovers_formal_theorems_and_prevents_recursive_model_collapse() {
    let mut rng = Rng::from_seed(20270101);
    let mut ctrl = RegimeId::R30PseudoSoftInheritance.preset_controls();
    ctrl.ambient_oxygen_level = 0.85;
    ctrl.multicellular_adhesion_rate = 0.72;
    ctrl.germ_soma_differentiation_rate = 0.72;
    ctrl.apoptosis_policing_rate = 0.80;
    ctrl.single_cell_bottleneck = true;
    ctrl.symbolic_teaching_fidelity = 0.94;
    ctrl.inter_deme_connectivity = 0.36;
    ctrl.cultural_innovation_rate = 0.75;
    ctrl.archival_externalization_rate = 0.80;
    ctrl.peer_verification_rigor = 0.88;
    ctrl.synthetic_generation_rate = 0.82;
    ctrl.formal_verifier_strictness = 0.92;
    ctrl.empirical_grounding_floor = 0.32;

    let host = ModeledUnit::new_founder_with_role(
        150,
        OrganismKind::ArtificialIntelligence,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let sym = ModeledUnit::new_founder_with_role(
        151,
        OrganismKind::ArtificialIntelligence,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let (mut inst_org, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &sym,
        152,
        BiomeZone::Alpha,
        1.25,
        0.40,
        0.52,
        &ctrl,
        0.94,
        24,
        0.4,
    );
    inst_org.step_clonal_multicellularity(true, &mut rng);
    inst_org.step_clonal_multicellularity(true, &mut rng);
    inst_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    inst_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    inst_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    inst_org.step_exosomatic_archive_transition(0.88, &mut rng);
    inst_org.step_exosomatic_archive_transition(0.88, &mut rng);
    inst_org.step_exosomatic_archive_transition(0.88, &mut rng);
    assert_eq!(inst_org.exosomatic_state, ExosomaticState::PeerVerifiedInstitution);
    inst_org.recompute_phenome(&ctrl, 1.25, 0.40, 0.52, 0.30, BiomeZone::Alpha, 0.94, 24, 0.4);
    let manual_inst_li = inst_org.complexity_c / inst_org.structure_s;

    // 1. Without formal verification (`0.20 < 0.72`) or empirical grounding (`0.02 < 0.12`),
    //    scaling generative throughput lands in `UnverifiedRecursiveGen` and suffers
    //    Shumailov et al. (2024) Recursive Model Collapse (tail loss -> delta-distribution collapse)!
    let mut unguarded = inst_org.clone();
    let (m1, ms1) = unguarded.step_machine_synthesis_transition(0.20, 0.02, &mut rng);
    assert!(m1 && ms1 == MachineSynthesisState::AlgorithmicMediation);
    let (m2, ms2) = unguarded.step_machine_synthesis_transition(0.20, 0.02, &mut rng);
    assert!(m2 && ms2 == MachineSynthesisState::UnverifiedRecursiveGen);

    let mut model_collapsed = false;
    for _ in 0..8 {
        let (_, col) = unguarded.step_recursive_synthesis_and_collapse_safeguards(0.85, 0.10, 0.01);
        if col {
            model_collapsed = true;
        }
    }
    assert!(
        model_collapsed
            && unguarded.distribution_tail_diversity < 0.35
            && unguarded.integrity < 0.25,
        "Unguarded recursive synthetic training without formal verification or empirical grounding must trigger Shumailov Recursive Model Collapse!"
    );

    // 2. Pairing a formal proof-checker (`0.92 >= 0.72`) with an empirical grounding floor (`0.32 >= 0.12`)
    //    transitions into `MachineVerifiedSymbiosis`, preserves >= 96% tail diversity, accumulates
    //    machine-checked theorems, and elevates thermodynamic leverage above manual archival institutions!
    let mut formal_symbio = inst_org.clone();
    formal_symbio.step_machine_synthesis_transition(0.92, 0.32, &mut rng);
    let (m_ok, ms_ok) = formal_symbio.step_machine_synthesis_transition(0.92, 0.32, &mut rng);
    assert!(m_ok && ms_ok == MachineSynthesisState::MachineVerifiedSymbiosis);

    for _ in 0..10 {
        let (_, col_ok) =
            formal_symbio.step_recursive_synthesis_and_collapse_safeguards(0.85, 0.92, 0.32);
        assert!(!col_ok);
    }
    formal_symbio.recompute_phenome(&ctrl, 1.25, 0.40, 0.52, 0.30, BiomeZone::Alpha, 0.94, 24, 0.4);
    let formal_symbio_li = formal_symbio.complexity_c / formal_symbio.structure_s;

    assert!(
        formal_symbio.cumulative_formal_proofs > 0
            && formal_symbio.distribution_tail_diversity >= 0.96
            && formal_symbio.computational_symbiont_sequence.len() > 38
            && formal_symbio_li > manual_inst_li * 1.35,
        "Machine-verified neuro-symbolic symbiosis must preserve empirical tails, grow verified theorems, and elevate L_I!"
    );

    // 3. Verify PhylogenyTracker preserves #H1 + WrittenArch + VerifiedInst + AlgoMed + FormalSymbio!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (150, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (151, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    tracker.record_endosymbiotic_fusion(152, 150, 151, 4, 0, BiomeZone::Alpha, gamma, 3, 2.2);
    tracker.record_cultural_transition(
        153,
        152,
        None,
        10,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        4,
        30.0,
        7,
        3,
        6.1,
    );
    tracker.record_exosomatic_transition(
        154,
        153,
        None,
        16,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        4,
        30.0,
        ExosomaticState::WrittenArchive,
        62.0,
        7,
        3,
        7.8,
    );
    tracker.record_exosomatic_transition(
        155,
        154,
        None,
        22,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        5,
        34.0,
        ExosomaticState::PeerVerifiedInstitution,
        118.0,
        7,
        3,
        10.2,
    );
    tracker.record_machine_synthesis_transition(
        156,
        155,
        None,
        28,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        5,
        34.0,
        ExosomaticState::PeerVerifiedInstitution,
        118.0,
        MachineSynthesisState::AlgorithmicMediation,
        26.0,
        0.95,
        7,
        3,
        11.6,
    );
    tracker.record_machine_synthesis_transition(
        157,
        156,
        None,
        34,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::SymbolicCumulativeCulture,
        5,
        34.0,
        ExosomaticState::PeerVerifiedInstitution,
        118.0,
        MachineSynthesisState::MachineVerifiedSymbiosis,
        142.5,
        0.98,
        7,
        3,
        15.8,
    );
    tracker.sync_extant_population_full(
        42,
        &[ExtantPhyloState {
            taxon_id: 157,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 15.8,
            egt_loci: 7,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::Solitary,
            colony_size: 1,
            sterile_caste_fraction: 0.0,
            agriculture_state: AgricultureState::HunterForager,
            cultivar_polyploidy: 1,
            gongylidia_yield: 0.0,
            cultural_state: CulturalState::SymbolicCumulativeCulture,
            cultural_ratchet_tier: 5,
            cultural_repertoire_bits: 34.0,
            exosomatic_state: ExosomaticState::PeerVerifiedInstitution,
            exosomatic_archive_bits: 118.0,
            machine_synthesis_state: MachineSynthesisState::MachineVerifiedSymbiosis,
            verified_theorem_bits: 142.5,
            distribution_tail_diversity: 0.98,
            planetary_epoch: PlanetaryState::OxygenatedPhanerozoic,
            chemical_disequilibrium_bits: 14.5,
            abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
            raf_closure_score: 0.96,
            chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
            genetic_code_mp: 0.79,
            genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
            chromosome_capacity_kb: 1850.0,
        }],
    );

    assert!(tracker.nodes.contains_key(&151));
    let sum = tracker.compute_summary(42);
    assert_eq!(sum.machine_symbiosis_transitions, 1);
    assert_eq!(sum.formal_verifier_transitions, 1);

    let nhx = tracker.to_newick(42, true);
    assert!(nhx.contains("MSYN=FormalSymbio:THMBITS=142.5:TAILS=0.98") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Machine-Verified Exosomatic Synthesis Test", 42);
    assert!(
        svg.contains("Symbolic→WrittenArch")
            && svg.contains("WrittenArch→VerifiedInst")
            && svg.contains("Inst→AlgoMed")
            && svg.contains("AlgoMed→FormalSymbio")
            && svg.contains("#H1")
    );
}

#[test]
fn planetary_biogeochemical_feedbacks_drive_goe_require_syntrophic_recycling_and_achieve_homeostatic_awakening_over_burnout() {
    let mut rng = Rng::from_seed(20270102);
    let mut ctrl = RegimeId::R30PseudoSoftInheritance.preset_controls();
    ctrl.ambient_oxygen_level = 0.08;
    ctrl.biogeochemical_feedback_gain = 0.85;
    ctrl.planetary_stewardship_coupling = 0.90;
    ctrl.extractive_burn_intensity = 0.18;

    // 1. Archean -> GOE -> OxygenatedPhanerozoic via coupled Autotrophs + Syntrophs + Superorganisms:
    let mut planet_closed = PlanetaryBiosphere::new_for_controls(&ctrl);
    assert_eq!(planet_closed.state, PlanetaryState::AnoxicArchean);
    let mut planet_no_syntrophs = planet_closed.clone();

    for _ in 0..24 {
        planet_closed.step_biogeochemical_feedbacks(1.25, 1.15, 0.85, 0.0, 0.0, &ctrl);
        planet_no_syntrophs.step_biogeochemical_feedbacks(1.25, 0.0, 0.0, 0.0, 0.0, &ctrl);
    }
    assert!(
        planet_closed.goe_crossed
            && planet_closed.atmospheric_o2 >= 0.58
            && planet_closed.state == PlanetaryState::OxygenatedPhanerozoic,
        "Autotroph + syntroph + superorganism guilds must oxygenate an anoxic Archean planet across the GOE!"
    );
    assert!(
        planet_no_syntrophs.nutrient_pool_np < 0.20
            && planet_closed.planetary_exergy_budget > planet_no_syntrophs.planetary_exergy_budget * 1.30,
        "Without syntrophic nutrient remineralization, N/P lock-up must throttle planetary exergy capture!"
    );

    // 2. Wong & Bartlett (2022) Asymptotic Burnout vs. Gaia-2.0 Homeostatic Awakening:
    let host = ModeledUnit::new_founder_with_role(
        170,
        OrganismKind::ArtificialIntelligence,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Producer,
        4,
        4,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    let sym = ModeledUnit::new_founder_with_role(
        171,
        OrganismKind::ArtificialIntelligence,
        RegimeId::R30PseudoSoftInheritance,
        TrophicGuild::Syntroph,
        4,
        5,
        BiomeZone::Alpha,
        0.4,
        &mut rng,
    );
    ctrl.ambient_oxygen_level = 0.82;
    let (mut steward_org, gamma) = ModeledUnit::fuse_endosymbiotic_holobiont(
        &host,
        &sym,
        172,
        BiomeZone::Alpha,
        1.30,
        0.42,
        0.55,
        &ctrl,
        0.94,
        24,
        0.4,
    );
    steward_org.step_clonal_multicellularity(true, &mut rng);
    steward_org.step_clonal_multicellularity(true, &mut rng);
    steward_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    steward_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    steward_org.step_cumulative_culture_transition(0.94, 0.36, &mut rng);
    steward_org.step_exosomatic_archive_transition(0.88, &mut rng);
    steward_org.step_exosomatic_archive_transition(0.88, &mut rng);
    steward_org.step_exosomatic_archive_transition(0.88, &mut rng);
    steward_org.step_machine_synthesis_transition(0.92, 0.35, &mut rng);
    steward_org.step_machine_synthesis_transition(0.92, 0.35, &mut rng);
    assert_eq!(
        steward_org.machine_synthesis_state,
        MachineSynthesisState::MachineVerifiedSymbiosis
    );

    let mut burnout_ctrl = ctrl.clone();
    burnout_ctrl.extractive_burn_intensity = 0.92;
    burnout_ctrl.planetary_stewardship_coupling = 0.05;
    let mut planet_burnout = PlanetaryBiosphere::new_for_controls(&burnout_ctrl);
    let mut burnout_org = steward_org.clone();
    let mut burnout_collapsed = false;

    let mut planet_awakened = PlanetaryBiosphere::new_for_controls(&ctrl);
    for _ in 0..18 {
        let (_, b_bad) = planet_burnout.step_biogeochemical_feedbacks(0.85, 0.80, 0.50, 1.35, 0.10, &burnout_ctrl);
        if b_bad && burnout_org.apply_planetary_boundary_stress(&planet_burnout) {
            burnout_collapsed = true;
        }
        let (_, b_ok) = planet_awakened.step_biogeochemical_feedbacks(1.10, 1.05, 0.85, 1.25, 1.40, &ctrl);
        assert!(!b_ok);
    }
    assert!(
        burnout_collapsed
            && planet_burnout.atmospheric_co2_ppm > 900.0
            && burnout_org.integrity < 0.25,
        "Unguided superlinear carbon extraction must breach planetary boundaries and trigger Asymptotic Burnout!"
    );
    assert!(
        planet_awakened.state == PlanetaryState::HomeostaticAwakening
            && planet_awakened.atmospheric_co2_ppm <= 340.0
            && planet_awakened.boundary_transgression_load < 0.05
            && planet_awakened.planetary_exergy_budget > 1.50,
        "Machine-verified stewardship must achieve Gaia-2.0 Homeostatic Awakening and expand planetary exergy!"
    );

    // 3. Verify PhylogenyTracker preserves #H1 + Archean->GOE + Techno->Gaia2Homeostasis in eNewick & SVG!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (170, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.5),
            (171, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.6),
        ],
    );
    if let Some(root_n) = tracker.nodes.get_mut(&170) {
        root_n.planetary_epoch = PlanetaryState::AnoxicArchean;
    }
    tracker.record_planetary_transition(
        172,
        170,
        6,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        1,
        MroClass::None,
        0,
        0,
        MulticellularState::Unicellular,
        1,
        0.0,
        EusocialState::Solitary,
        1,
        0.0,
        AgricultureState::HunterForager,
        1,
        0.0,
        CulturalState::AsocialTrialError,
        0,
        0.0,
        ExosomaticState::EndosomaticBound,
        0.0,
        MachineSynthesisState::ManualArchival,
        0.0,
        1.0,
        PlanetaryState::OxygenatedPhanerozoic,
        14.8,
        0,
        0,
        3.1,
    );
    tracker.record_endosymbiotic_fusion(173, 172, 171, 12, 0, BiomeZone::Alpha, gamma, 3, 4.2);
    tracker.record_planetary_transition(
        174,
        173,
        28,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        1,
        MroClass::AerobicMitochondrion,
        4,
        3,
        MulticellularState::GermSomaMulticellular,
        32,
        0.75,
        EusocialState::ObligateSuperorganism,
        64,
        0.90,
        AgricultureState::LeafCutterAgroEcosystem,
        8,
        0.96,
        CulturalState::SymbolicCumulativeCulture,
        6,
        44.0,
        ExosomaticState::PeerVerifiedInstitution,
        122.0,
        MachineSynthesisState::MachineVerifiedSymbiosis,
        172.0,
        0.99,
        PlanetaryState::HomeostaticAwakening,
        23.4,
        8,
        3,
        22.5,
    );
    tracker.sync_extant_population_full(
        36,
        &[ExtantPhyloState {
            taxon_id: 174,
            species_id: 0,
            guild: TrophicGuild::Holobiont,
            biome: BiomeZone::Alpha,
            leverage_li: 22.5,
            egt_loci: 8,
            organelle_rem: 3,
            endo_tier: 1,
            nucleomorph_loci: 0,
            nucleomorph_eliminated: false,
            klepto_hgt_loci: 0,
            has_kleptoplast: false,
            has_apicoplast: false,
            apicoplast_pathways: 4,
            apicoplast_eliminated: false,
            mro_class: MroClass::AerobicMitochondrion,
            mro_pathways: 4,
            mro_mtdna_loci: 3,
            multicell_state: MulticellularState::GermSomaMulticellular,
            cell_count: 32,
            soma_fraction: 0.75,
            eusocial_state: EusocialState::ObligateSuperorganism,
            colony_size: 64,
            sterile_caste_fraction: 0.90,
            agriculture_state: AgricultureState::LeafCutterAgroEcosystem,
            cultivar_polyploidy: 8,
            gongylidia_yield: 0.96,
            cultural_state: CulturalState::SymbolicCumulativeCulture,
            cultural_ratchet_tier: 6,
            cultural_repertoire_bits: 44.0,
            exosomatic_state: ExosomaticState::PeerVerifiedInstitution,
            exosomatic_archive_bits: 122.0,
            machine_synthesis_state: MachineSynthesisState::MachineVerifiedSymbiosis,
            verified_theorem_bits: 172.0,
            distribution_tail_diversity: 0.99,
            planetary_epoch: PlanetaryState::HomeostaticAwakening,
            chemical_disequilibrium_bits: 23.4,
            abiogenesis_stage: AbiogenesisStage::LucaCellularLineage,
            raf_closure_score: 0.96,
            chirality_code_stage: ChiralityCodeStage::UniversalErrorMinimizedSgc,
            genetic_code_mp: 0.79,
            genome_substrate_stage: GenomeSubstrateStage::CoevolvedCellularVirosphere,
            chromosome_capacity_kb: 1850.0,
        }],
    );

    assert!(tracker.nodes.contains_key(&171));
    let sum = tracker.compute_summary(36);
    assert_eq!(sum.goe_transitions, 1);
    assert_eq!(sum.homeostatic_awakening_transitions, 1);

    let nhx = tracker.to_newick(36, true);
    assert!(nhx.contains("GAIA=Gaia2Homeostasis:DISEQ=23.4") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Planetary Biogeochemical Feedbacks & Gaia-2.0 Test", 36);
    assert!(
        svg.contains("Archean→GOE")
            && svg.contains("Techno→Gaia2Homeostasis")
            && svg.contains("#H1")
    );
}

#[test]
fn hadean_abiogenesis_crosses_raf_eigen_and_darwinian_thresholds_while_unchelated_mg_and_high_error_collapse() {
    let mut rng = Rng::from_seed(20270402);
    let mut base_ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    base_ctrl.hadean_abiogenesis_mode = true;

    // 1. Bulk dilute ocean WITHOUT wet-dry cycles OR alkaline vent Delta-pH stalls in Water Paradox (<= 10 nt):
    let mut dilute_ctrl = base_ctrl.clone();
    dilute_ctrl.wet_dry_cycle_amplitude = 0.04;
    dilute_ctrl.hydrothermal_delta_ph = 0.20;
    let mut dilute_bath = HadeanGeochemicalBath::new_for_controls(&dilute_ctrl);
    let mut cand_dilute = ModeledUnit::new_founder_with_role(
        200,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        2,
        2,
        BiomeZone::Alpha,
        0.0,
        &mut rng,
    );
    cand_dilute.seed_hadean_prebiotic_pool(&dilute_bath);
    for _ in 0..8 {
        dilute_bath.step_geochemical_cycle(&dilute_ctrl, cand_dilute.raf_closure_score, 0.0);
        cand_dilute.step_hadean_abiogenesis(&dilute_bath, &mut rng);
    }
    assert_eq!(cand_dilute.abiogenesis_stage, AbiogenesisStage::PrebioticMonomerSoup);
    assert!(
        cand_dilute.rna_oligomer_length_nt <= 10,
        "Without wet-dry condensation films or hydrothermal Delta-pH pores, hydrolysis must stall oligomerization!"
    );

    // 2. Unchelated 50 mM Mg2+ precipitates fatty-acid protocell vesicles (Adamala & Szostak, 2013):
    let mut unchelated_ctrl = base_ctrl.clone();
    unchelated_ctrl.wet_dry_cycle_amplitude = 0.85;
    unchelated_ctrl.citrate_mg_chelation = 0.10;
    let unchelated_bath = HadeanGeochemicalBath::new_for_controls(&unchelated_ctrl);
    let mut cand_mg = cand_dilute.clone();
    cand_mg.integrity = 0.90;
    cand_mg.seed_hadean_prebiotic_pool(&unchelated_bath);
    let (_, _, mg_precipitated) = cand_mg.step_hadean_abiogenesis(&unchelated_bath, &mut rng);
    assert!(
        mg_precipitated && cand_mg.integrity < 0.30,
        "Unchelated high Mg2+ must disrupt fatty-acid membranes unless chelated by citrate/acidic peptides!"
    );

    // 3. Low ribozyme replication fidelity triggers Eigen's Error Catastrophe (`L > ln(sigma)/mu`):
    let mut low_fid_ctrl = base_ctrl.clone();
    low_fid_ctrl.wet_dry_cycle_amplitude = 0.86;
    low_fid_ctrl.hydrothermal_delta_ph = 3.20;
    low_fid_ctrl.citrate_mg_chelation = 0.88;
    low_fid_ctrl.ribozyme_polymerase_fidelity = 0.04;
    let mut low_fid_bath = HadeanGeochemicalBath::new_for_controls(&low_fid_ctrl);
    let mut cand_eigen = cand_dilute.clone();
    cand_eigen.integrity = 0.92;
    cand_eigen.seed_hadean_prebiotic_pool(&low_fid_bath);
    let mut hit_eigen_cliff = false;
    for _ in 0..8 {
        low_fid_bath.step_geochemical_cycle(&low_fid_ctrl, cand_eigen.raf_closure_score, 0.0);
        let (_, eig_bad, _) = cand_eigen.step_hadean_abiogenesis(&low_fid_bath, &mut rng);
        if eig_bad {
            hit_eigen_cliff = true;
            break;
        }
    }
    assert!(
        hit_eigen_cliff && cand_eigen.abiogenesis_stage != AbiogenesisStage::LucaCellularLineage,
        "Replicating oligomers exceeding Eigen's error threshold L_max = ln(sigma)/mu must suffer error catastrophe!"
    );

    // 4. Coupled Hot-Spring Wet-Dry + Citrate-Mg2+ + Alkaline Fe(Ni)S Delta-pH + Ribozyme Hypercycle reaches LUCA!
    let mut coupled_ctrl = base_ctrl.clone();
    coupled_ctrl.wet_dry_cycle_amplitude = 0.88;
    coupled_ctrl.hydrothermal_delta_ph = 3.50;
    coupled_ctrl.citrate_mg_chelation = 0.90;
    coupled_ctrl.ribozyme_polymerase_fidelity = 0.96;
    let mut coupled_bath = HadeanGeochemicalBath::new_for_controls(&coupled_ctrl);
    let mut cand_luca = cand_dilute.clone();
    cand_luca.integrity = 0.95;
    cand_luca.seed_hadean_prebiotic_pool(&coupled_bath);

    let mut seen_raf = false;
    let mut seen_ribo = false;
    let mut seen_progenote = false;
    for _ in 0..10 {
        let ribo_f = if cand_luca.abiogenesis_stage.has_ribozyme_heredity() { 1.0 } else { 0.35 };
        coupled_bath.step_geochemical_cycle(&coupled_ctrl, cand_luca.raf_closure_score, ribo_f);
        let (st, eig_bad, mg_bad) = cand_luca.step_hadean_abiogenesis(&coupled_bath, &mut rng);
        assert!(!eig_bad && !mg_bad);
        if st == AbiogenesisStage::AutocatalyticRafProtocell {
            seen_raf = true;
        } else if st == AbiogenesisStage::RnaRibozymeHypercycle {
            seen_ribo = true;
        } else if st == AbiogenesisStage::ProgenoteCommunity {
            seen_progenote = true;
        }
    }
    assert!(
        seen_raf
            && seen_ribo
            && seen_progenote
            && cand_luca.abiogenesis_stage == AbiogenesisStage::LucaCellularLineage
            && cand_luca.rna_oligomer_length_nt >= 380
            && cand_luca.progenote_hgt_flux <= 0.15,
        "Coupled Hadean hot-spring + hydrothermal cradle must traverse RAF -> Ribozyme -> Progenote -> LUCA!"
    );

    // 5. Verify PhylogenyTracker exports Soup->RAFProtocell, Progenote->LUCA, and #H1 communal HGT!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (200, 0, TrophicGuild::Producer, BiomeZone::Alpha, 0.9),
            (201, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 1.0),
        ],
    );
    if let Some(r0) = tracker.nodes.get_mut(&200) {
        r0.abiogenesis_stage = AbiogenesisStage::PrebioticMonomerSoup;
    }
    tracker.record_abiogenesis_transition(
        202,
        200,
        4,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        AbiogenesisStage::AutocatalyticRafProtocell,
        0.74,
        1.45,
    );
    tracker.record_endosymbiotic_fusion(203, 202, 201, 9, 0, BiomeZone::Alpha, 0.42, 2, 2.15);
    if let Some(n203) = tracker.nodes.get_mut(&203) {
        n203.abiogenesis_stage = AbiogenesisStage::ProgenoteCommunity;
    }
    tracker.record_abiogenesis_transition(
        204,
        203,
        15,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        AbiogenesisStage::LucaCellularLineage,
        0.96,
        3.60,
    );
    let sum = tracker.compute_summary(18);
    assert_eq!(sum.raf_protocell_transitions, 1);
    assert_eq!(sum.darwinian_threshold_transitions, 1);
    let nhx = tracker.to_newick(18, true);
    assert!(nhx.contains("DARWIN_TR=1") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Hadean Abiogenesis & Progenote-to-LUCA Test", 18);
    assert!(
        svg.contains("Soup→RAFProtocell")
            && svg.contains("Progenote→LUCA")
            && svg.contains("#H1")
    );
}

#[test]
fn prebiotic_homochiral_symmetry_breaking_minihelix_aminoacylation_and_collective_sgc_crystallization_overcome_cross_inhibition_and_babel_traps() {
    let mut rng = Rng::from_seed(20270804);
    let mut base_ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    base_ctrl.chiral_symmetry_breaking_mode = true;

    // 1. Unbiased racemic monomer pool (`magnetite_ciss_spin_bias = 0.04`) suffers Joyce-Orgel enantiomeric cross-inhibition:
    let mut racemic_ctrl = base_ctrl.clone();
    racemic_ctrl.magnetite_ciss_spin_bias = 0.04;
    let mut racemic_reactor = ChiralCodeReactor::new_for_controls(&racemic_ctrl);
    let mut cand_racemic = ModeledUnit::new_founder_with_role(
        220,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        3,
        3,
        BiomeZone::Alpha,
        0.0,
        &mut rng,
    );
    cand_racemic.seed_racemic_prebiotic_pool(&racemic_reactor);
    racemic_reactor.step_chiral_and_code_dynamics(&racemic_ctrl);
    let (_, cross_inhibited, _) =
        cand_racemic.step_homochiral_code_crystallization(&racemic_reactor, &mut rng);
    assert!(
        cross_inhibited && cand_racemic.integrity < 0.30 && racemic_reactor.ribose_ee_d < 0.25,
        "Racemic monomer pools without magnetite CISS / Frank amplification must suffer 3' enantiomeric cross-inhibition!"
    );

    // 2. Without Tamura-Schimmel 3'-NCCA minihelix steric selectivity (`s_ac = 0.05`), D-RNA forms but peptides stay racemic:
    let mut no_minihelix_ctrl = base_ctrl.clone();
    no_minihelix_ctrl.magnetite_ciss_spin_bias = 0.92;
    no_minihelix_ctrl.minihelix_chiral_selectivity = 0.05;
    let mut no_minihelix_reactor = ChiralCodeReactor::new_for_controls(&no_minihelix_ctrl);
    let mut cand_no_acyl = cand_racemic.clone();
    cand_no_acyl.integrity = 0.92;
    cand_no_acyl.seed_racemic_prebiotic_pool(&no_minihelix_reactor);
    for _ in 0..8 {
        no_minihelix_reactor.step_chiral_and_code_dynamics(&no_minihelix_ctrl);
        cand_no_acyl.step_homochiral_code_crystallization(&no_minihelix_reactor, &mut rng);
    }
    assert!(
        cand_no_acyl.chirality_code_stage == ChiralityCodeStage::HomochiralRnaSymmetryBroken
            && no_minihelix_reactor.amino_ee_l < 0.25,
        "Without 3'-NCCA minihelix steric selectivity, D-RNA homochirality cannot propagate to L-aminoacyl peptides!"
    );

    // 3. Isolated vertical lineages (`communal_code_hgt_rate = 0.04`) freeze in the suboptimal "Tower of Babel" code trap:
    let mut babel_ctrl = base_ctrl.clone();
    babel_ctrl.magnetite_ciss_spin_bias = 0.94;
    babel_ctrl.minihelix_chiral_selectivity = 0.90;
    babel_ctrl.communal_code_hgt_rate = 0.04;
    let mut babel_reactor = ChiralCodeReactor::new_for_controls(&babel_ctrl);
    let mut cand_babel = cand_racemic.clone();
    cand_babel.integrity = 0.95;
    cand_babel.seed_racemic_prebiotic_pool(&babel_reactor);
    let mut hit_babel_trap = false;
    for _ in 0..10 {
        babel_reactor.step_chiral_and_code_dynamics(&babel_ctrl);
        let (_, _, babel_bad) =
            cand_babel.step_homochiral_code_crystallization(&babel_reactor, &mut rng);
        if babel_bad {
            hit_babel_trap = true;
            break;
        }
    }
    assert!(
        hit_babel_trap
            && babel_reactor.sgc_error_minimization_mp < 0.45
            && !babel_reactor.sgc_crystallized,
        "Without Vetsigian-Woese communal HGT innovation sharing, lineages must stall in fragmented suboptimal codes!"
    );

    // 4. Full Magnetite CISS + Tamura-Schimmel Minihelix + Communal Progenote HGT crystallizes the 20-AA SGC (MP >= 0.76):
    let mut sgc_ctrl = base_ctrl.clone();
    sgc_ctrl.magnetite_ciss_spin_bias = 0.94;
    sgc_ctrl.minihelix_chiral_selectivity = 0.92;
    sgc_ctrl.communal_code_hgt_rate = 0.90;
    let mut sgc_reactor = ChiralCodeReactor::new_for_controls(&sgc_ctrl);
    let mut cand_sgc = cand_racemic.clone();
    cand_sgc.integrity = 0.95;
    cand_sgc.seed_racemic_prebiotic_pool(&sgc_reactor);
    for _ in 0..10 {
        sgc_reactor.step_chiral_and_code_dynamics(&sgc_ctrl);
        let (_, ci, bb) = cand_sgc.step_homochiral_code_crystallization(&sgc_reactor, &mut rng);
        assert!(!ci && !bb);
    }
    assert!(
        cand_sgc.chirality_code_stage == ChiralityCodeStage::UniversalErrorMinimizedSgc
            && sgc_reactor.ribose_ee_d >= 0.95
            && sgc_reactor.amino_ee_l >= 0.90
            && sgc_reactor.encoded_amino_acids == 20
            && sgc_reactor.sgc_error_minimization_mp >= 0.76,
        "Magnetite CISS + minihelix L-aminoacylation + communal HGT must crystallize the universal 20-AA error-minimized SGC!"
    );

    // 5. Verify PhylogenyTracker exports Racemic->Homochiral, StatCode->UniversalSGC, and #H1!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (220, 0, TrophicGuild::Producer, BiomeZone::Alpha, 0.8),
            (221, 1, TrophicGuild::Syntroph, BiomeZone::Alpha, 0.9),
        ],
    );
    if let Some(r0) = tracker.nodes.get_mut(&220) {
        r0.chirality_code_stage = ChiralityCodeStage::RacemicCrossInhibited;
    }
    tracker.record_chirality_code_transition(
        222,
        220,
        4,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        ChiralityCodeStage::HomochiralRnaSymmetryBroken,
        0.24,
        1.5,
    );
    tracker.record_endosymbiotic_fusion(223, 222, 221, 9, 0, BiomeZone::Alpha, 0.48, 2, 2.3);
    if let Some(n223) = tracker.nodes.get_mut(&223) {
        n223.chirality_code_stage = ChiralityCodeStage::AmbiguousStatisticalCode;
    }
    tracker.record_chirality_code_transition(
        224,
        223,
        15,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        ChiralityCodeStage::UniversalErrorMinimizedSgc,
        0.79,
        4.2,
    );
    let sum = tracker.compute_summary(18);
    assert_eq!(sum.homochiral_symmetry_transitions, 1);
    assert_eq!(sum.sgc_crystallization_transitions, 1);
    let nhx = tracker.to_newick(18, true);
    assert!(nhx.contains("SGC_TR=1") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("Homochiral Symmetry Breaking & SGC Crystallization Test", 18);
    assert!(
        svg.contains("Racemic→Homochiral")
            && svg.contains("StatCode→UniversalSGC")
            && svg.contains("#H1")
    );
}

#[test]
fn rnr_deoxyribose_reduction_thymine_udg_deamination_repair_and_virosphere_coevolution_unlock_megabase_dna_chromosomes() {
    let mut rng = Rng::from_seed(20271112);
    let mut base_ctrl = RegimeId::R35Symbiogenetic.preset_controls();
    base_ctrl.rna_to_dna_virosphere_mode = true;

    // 1. Without Class III RNR 2'-OH->2'-H radical reduction (`rnr = 0.08`), labile RNA chromosomes cap < 28 kb:
    let mut rna_only_ctrl = base_ctrl.clone();
    rna_only_ctrl.rnr_thiyl_radical_activity = 0.08;
    let mut rna_reactor = RnrVirosphereReactor::new_for_controls(&rna_only_ctrl);
    let mut ribocell = ModeledUnit::new_founder_with_role(
        240,
        OrganismKind::RnaDnaProtein,
        RegimeId::R35Symbiogenetic,
        TrophicGuild::Producer,
        3,
        3,
        BiomeZone::Alpha,
        0.0,
        &mut rng,
    );
    ribocell.seed_labile_rna_cell(&rna_reactor);
    for _ in 0..8 {
        rna_reactor.step_rnr_and_virosphere_cycle(&rna_only_ctrl);
        ribocell.step_rnr_dna_and_virosphere_coevolution(&rna_reactor, &mut rng);
    }
    assert!(
        ribocell.genome_substrate_stage == GenomeSubstrateStage::LabileRnaChromosome
            && ribocell.chromosome_capacity_kb < 28.0,
        "Without RNR 2'-OH reduction, intramolecular transesterification cleavage must cap RNA chromosomes below 28 kb!"
    );

    // 2. U-DNA WITHOUT Thymine + UDG excision (`udg = 0.06`) suffers C->U cytosine deamination meltdown:
    let mut udna_no_udg_ctrl = base_ctrl.clone();
    udna_no_udg_ctrl.rnr_thiyl_radical_activity = 0.90;
    udna_no_udg_ctrl.thymidylate_udg_repair_efficiency = 0.06;
    let mut udna_reactor = RnrVirosphereReactor::new_for_controls(&udna_no_udg_ctrl);
    let mut cand_udna = ribocell.clone();
    cand_udna.integrity = 0.94;
    cand_udna.seed_labile_rna_cell(&udna_reactor);
    let mut deam_collapsed = false;
    for _ in 0..8 {
        udna_reactor.step_rnr_and_virosphere_cycle(&udna_no_udg_ctrl);
        let (_, deam_bad, _) =
            cand_udna.step_rnr_dna_and_virosphere_coevolution(&udna_reactor, &mut rng);
        if deam_bad {
            deam_collapsed = true;
            break;
        }
    }
    assert!(
        deam_collapsed && cand_udna.integrity < 0.30,
        "Uracil-DNA without ThyA/ThyX methylation and UDG excision must melt under unrecognized C->U deamination!"
    );

    // 3. Unchecked lytic viral burst (`lytic = 0.94`, `arms_race_coupling = 0.08`) lyses unshielded cells:
    let mut lytic_ctrl = base_ctrl.clone();
    lytic_ctrl.lytic_viral_burst_pressure = 0.94;
    lytic_ctrl.virosphere_arms_race_coupling = 0.08;
    let mut lytic_reactor = RnrVirosphereReactor::new_for_controls(&lytic_ctrl);
    let mut cand_lytic = ribocell.clone();
    cand_lytic.integrity = 0.94;
    cand_lytic.seed_labile_rna_cell(&lytic_reactor);
    lytic_reactor.step_rnr_and_virosphere_cycle(&lytic_ctrl);
    let (_, _, lytic_lysed) =
        cand_lytic.step_rnr_dna_and_virosphere_coevolution(&lytic_reactor, &mut rng);
    assert!(
        lytic_lysed && cand_lytic.integrity < 0.30,
        "Runaway lytic viral predation without cellular CRISPR/restriction immunity must trigger epidemic lysis!"
    );

    // 4. Full RNR + ThyX/UDG + Coevolved Virosphere Polymerase Transfer unlocks >1,200 kb Megabase T-DNA!
    let mut tdna_ctrl = base_ctrl.clone();
    tdna_ctrl.rnr_thiyl_radical_activity = 0.94;
    tdna_ctrl.thymidylate_udg_repair_efficiency = 0.94;
    tdna_ctrl.virosphere_arms_race_coupling = 0.90;
    tdna_ctrl.lytic_viral_burst_pressure = 0.22;
    let mut tdna_reactor = RnrVirosphereReactor::new_for_controls(&tdna_ctrl);
    let mut cand_tdna = ribocell.clone();
    cand_tdna.integrity = 0.96;
    cand_tdna.seed_labile_rna_cell(&tdna_reactor);
    for _ in 0..10 {
        tdna_reactor.step_rnr_and_virosphere_cycle(&tdna_ctrl);
        let (_, db, lb) =
            cand_tdna.step_rnr_dna_and_virosphere_coevolution(&tdna_reactor, &mut rng);
        assert!(!db && !lb);
    }
    assert!(
        cand_tdna.genome_substrate_stage == GenomeSubstrateStage::CoevolvedCellularVirosphere
            && cand_tdna.chromosome_capacity_kb >= 1200.0
            && tdna_reactor.deamination_lesion_load < 0.05,
        "Coupled RNR + ThyX/UDG + temperate viral polymerase transfer must establish megabase T-DNA chromosomes!"
    );

    // 5. Verify PhylogenyTracker exports RNA->U-DNA(RNR), U-DNA->T-DNA(UDG), and #H1 viral transduction!
    let mut tracker = PhylogenyTracker::new();
    tracker.reset_with_founders(
        0,
        &[
            (240, 0, TrophicGuild::Producer, BiomeZone::Alpha, 1.1),
            (241, 1, TrophicGuild::PredatorParasite, BiomeZone::Alpha, 1.2),
        ],
    );
    if let Some(r0) = tracker.nodes.get_mut(&240) {
        r0.genome_substrate_stage = GenomeSubstrateStage::LabileRnaChromosome;
    }
    tracker.record_rna_to_dna_transition(
        242,
        240,
        4,
        0,
        TrophicGuild::Producer,
        BiomeZone::Alpha,
        GenomeSubstrateStage::ViralModifiedUDna,
        92.0,
        2.1,
    );
    tracker.record_endosymbiotic_fusion(243, 242, 241, 9, 0, BiomeZone::Alpha, 0.36, 2, 2.9);
    if let Some(n243) = tracker.nodes.get_mut(&243) {
        n243.genome_substrate_stage = GenomeSubstrateStage::ViralModifiedUDna;
    }
    tracker.record_rna_to_dna_transition(
        244,
        243,
        15,
        0,
        TrophicGuild::Holobiont,
        BiomeZone::Alpha,
        GenomeSubstrateStage::CoevolvedCellularVirosphere,
        1780.0,
        5.1,
    );
    let sum = tracker.compute_summary(18);
    assert_eq!(sum.rnr_dna_transitions, 1);
    assert_eq!(sum.thymine_udg_transitions, 1);
    let nhx = tracker.to_newick(18, true);
    assert!(nhx.contains("UDG_TR=1") && nhx.contains("#H1"));
    let svg = tracker.to_svg_phylogram("RNR RNA-to-DNA, UDG Repair & Virosphere Test", 18);
    assert!(
        svg.contains("RNA→U-DNA(RNR)")
            && svg.contains("U-DNA→T-DNA(UDG)")
            && svg.contains("#H1")
    );
}

#[test]
fn continuous_multi_epoch_grand_trajectory_chains_homochirality_sgc_dna_luca_symbiogenesis_multicellularity_archives_and_gaia_homeostasis() {
    let mut world = SimulationWorld::new(
        OrganismKind::ArtificialIntelligence,
        RegimeId::R35Symbiogenetic,
        SimulationMode::OrganicSandbox,
        12,
        20280101,
    );
    let epochs = world.run_continuous_grand_trajectory(8);
    assert_eq!(epochs.len(), 8);

    // 1. Verify unbroken population persistence and strict monotonic growth of MI & Leverage across all 8 epochs:
    for i in 0..epochs.len() {
        assert!(
            epochs[i].population_size > 0,
            "Population must remain viable across every epoch of the 4.25 Ga trajectory!"
        );
        if i > 0 {
            assert!(
                epochs[i].heritable_mutual_info_bits > epochs[i - 1].heritable_mutual_info_bits,
                "Heritable mutual information must expand monotonically across major evolutionary transitions!"
            );
            assert!(
                epochs[i].informational_leverage_li > epochs[i - 1].informational_leverage_li,
                "Informational leverage L_I must climb monotonically across major evolutionary transitions!"
            );
            assert!(
                epochs[i].heritable_channels_active >= epochs[i - 1].heritable_channels_active,
                "Active inheritance channels must expand from single-channel RNA to 5-channel Gaia-2.0!"
            );
        }
    }

    // 2. Verify planetary oxygenation from Anoxic Hadean/Archean (Epoch 1) to Gaia-2.0 Homeostasis (Epoch 8):
    assert!(epochs[0].planetary_o2 < 0.20);
    let final_ep = epochs.last().unwrap();
    assert!(
        final_ep.planetary_o2 >= 0.60
            && final_ep.planetary_co2_ppm <= 360.0
            && final_ep.planetary_exergy_budget >= 1.50
            && final_ep.heritable_channels_active == 5
    );

    // 3. Verify that the single continuous PhylogenyTracker recorded transitions across the full hierarchy!
    let phylo_sum = world.phylogeny.compute_summary(world.step_count);
    assert!(phylo_sum.homochiral_symmetry_transitions >= 1);
    assert!(phylo_sum.sgc_crystallization_transitions >= 1);
    assert!(phylo_sum.raf_protocell_transitions >= 1);
    assert!(phylo_sum.darwinian_threshold_transitions >= 1);
    assert!(phylo_sum.rnr_dna_transitions >= 1);
    assert!(phylo_sum.thymine_udg_transitions >= 1);
    assert!(phylo_sum.endosymbiotic_fusions >= 1);
    assert!(phylo_sum.germ_soma_transitions >= 1);
    assert!(phylo_sum.superorganism_transitions >= 1);
    assert!(phylo_sum.obligate_crop_transitions >= 1);
    assert!(phylo_sum.symbolic_teaching_transitions >= 1);
    assert!(phylo_sum.peer_institution_transitions >= 1);
    assert!(phylo_sum.formal_verifier_transitions >= 1);
    assert!(phylo_sum.homeostatic_awakening_transitions >= 1);

    let nhx = world.export_newick(true);
    assert!(nhx.contains("#H1") && nhx.contains("GAIA=Gaia2Homeostasis"));
    let svg = world.export_phylogram_svg();
    assert!(svg.contains("Techno→Gaia2Homeostasis"));
}
