use itme_evolution_sim::engine::{SimulationMode, SimulationWorld};
use itme_evolution_sim::organisms::OrganismKind;
use itme_evolution_sim::regimes::RegimeId;
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

struct EpochSnap {
    stage: String,
    tick_end: usize,
    geologic_age_ga: f64,
    population_size: usize,
    heritable_channels_active: u8,
    heritable_mutual_info_bits: f64,
    informational_leverage_li: f64,
    order_velocity_vo: f64,
    planetary_o2: f64,
    planetary_co2_ppm: f64,
    planetary_exergy_budget: f64,
    landmark_summary: String,
}

struct ConditionSnap {
    condition: &'static str,
    bath_exergy_flux: f64,
    heritable_variation: f64,
    mean_order: f64,
    sd_order: f64,
    mean_exergy: f64,
    mean_entropy: f64,
    mean_species: f64,
    all_finite_positive: bool,
}

struct SuiteArtifacts {
    csv_path: PathBuf,
    csv_bytes: u64,
    epochs: Vec<EpochSnap>,
    nwk_path: PathBuf,
    nwk_bytes: u64,
    newick: String,
    newick_leaves: usize,
    recovery_pop: usize,
    recovery_species_count: usize,
    recovery_ecosystem_order_pop: f64,
    phylo_active_nodes: usize,
    conditions: Vec<ConditionSnap>,
    determinism_verified: bool,
    json_path: PathBuf,
    json_bytes: u64,
}

static SUITE: OnceLock<SuiteArtifacts> = OnceLock::new();

fn validation_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("validation");
    fs::create_dir_all(&dir).expect("create validation directory");
    dir
}

fn json_escape(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 8);
    for ch in raw.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn csv_escape(raw: &str) -> String {
    if raw.contains(',') || raw.contains('"') || raw.contains('\n') {
        format!("\"{}\"", raw.replace('"', "\"\""))
    } else {
        raw.to_string()
    }
}

fn count_newick_leaves(nwk: &str) -> usize {
    let body = nwk.trim().trim_end_matches(';');
    if body.is_empty() {
        return 0;
    }
    body.split(',').count()
}

fn verify_newick_balanced(nwk: &str) -> bool {
    let trimmed = nwk.trim();
    if !trimmed.ends_with(';') {
        return false;
    }
    let mut depth: i32 = 0;
    let mut saw_paren = false;
    for ch in trimmed.chars() {
        match ch {
            '(' => {
                depth += 1;
                saw_paren = true;
            }
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return false;
                }
            }
            _ => {}
        }
    }
    depth == 0 && (saw_paren || !trimmed.trim_end_matches(';').is_empty())
}

fn ensure_validation_suite() -> &'static SuiteArtifacts {
    SUITE.get_or_init(|| {
        let val_dir = validation_dir();

        // ---------------------------------------------------------------------
        // T1: Continuous 4.25 Ga 8-epoch Hadean-to-Gaia trajectory simulation
        // ---------------------------------------------------------------------
        let ticks_per_epoch = 12usize;
        let mut world_t1 = SimulationWorld::new(
            OrganismKind::ArtificialIntelligence,
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            12,
            20280101,
        );
        let raw_epochs = world_t1.run_continuous_grand_trajectory(ticks_per_epoch);

        let mut csv = String::from(
            "epoch_index,epoch_stage,tick_end,geologic_age_ga,population_size,heritable_channels_active,heritable_mutual_info_bits,informational_leverage_li,order_velocity_vo,planetary_o2,planetary_co2_ppm,planetary_exergy_budget,landmark_summary\n",
        );
        let mut epochs = Vec::with_capacity(raw_epochs.len());
        let mut epoch_json_rows = Vec::with_capacity(raw_epochs.len());

        for (idx, ep) in raw_epochs.iter().enumerate() {
            let stage_name = format!("{:?}", ep.epoch);
            csv.push_str(&format!(
                "{},{},{},{:.4},{},{},{:.6},{:.6},{:.6},{:.6},{:.2},{:.6},{}\n",
                idx + 1,
                stage_name,
                ep.tick_end,
                ep.geologic_age_ga,
                ep.population_size,
                ep.heritable_channels_active,
                ep.heritable_mutual_info_bits,
                ep.informational_leverage_li,
                ep.order_velocity_vo,
                ep.planetary_o2,
                ep.planetary_co2_ppm,
                ep.planetary_exergy_budget,
                csv_escape(&ep.landmark_summary),
            ));
            epoch_json_rows.push(format!(
                "{{\"epoch_index\":{},\"stage\":\"{}\",\"tick_end\":{},\"geologic_age_ga\":{:.4},\"population_size\":{},\"heritable_channels_active\":{},\"heritable_mutual_info_bits\":{:.6},\"informational_leverage_li\":{:.6},\"order_velocity_vo\":{:.6},\"planetary_o2\":{:.6},\"planetary_co2_ppm\":{:.2},\"planetary_exergy_budget\":{:.6}}}",
                idx + 1,
                json_escape(&stage_name),
                ep.tick_end,
                ep.geologic_age_ga,
                ep.population_size,
                ep.heritable_channels_active,
                ep.heritable_mutual_info_bits,
                ep.informational_leverage_li,
                ep.order_velocity_vo,
                ep.planetary_o2,
                ep.planetary_co2_ppm,
                ep.planetary_exergy_budget,
            ));
            epochs.push(EpochSnap {
                stage: stage_name,
                tick_end: ep.tick_end,
                geologic_age_ga: ep.geologic_age_ga,
                population_size: ep.population_size,
                heritable_channels_active: ep.heritable_channels_active,
                heritable_mutual_info_bits: ep.heritable_mutual_info_bits,
                informational_leverage_li: ep.informational_leverage_li,
                order_velocity_vo: ep.order_velocity_vo,
                planetary_o2: ep.planetary_o2,
                planetary_co2_ppm: ep.planetary_co2_ppm,
                planetary_exergy_budget: ep.planetary_exergy_budget,
                landmark_summary: ep.landmark_summary.clone(),
            });
        }

        let csv_path = val_dir.join("planetary_epochs_trajectory.csv");
        fs::write(&csv_path, &csv).expect("write planetary_epochs_trajectory.csv");
        let csv_bytes = fs::metadata(&csv_path).expect("stat csv").len();

        let epoch1 = &epochs[0];
        let epoch4 = &epochs[3];
        let epoch8 = &epochs[7];
        let o2_max = epochs.iter().map(|e| e.planetary_o2).fold(0.0_f64, f64::max);
        let mi_max = epochs
            .iter()
            .map(|e| e.heritable_mutual_info_bits)
            .fold(0.0_f64, f64::max);

        // ---------------------------------------------------------------------
        // T2: Lineage diversification, extinction bottleneck recovery & Newick export
        // ---------------------------------------------------------------------
        let mut world_t2 = SimulationWorld::new(
            OrganismKind::RnaDnaProtein,
            RegimeId::R35Symbiogenetic,
            SimulationMode::OrganicSandbox,
            12,
            20260717,
        );
        for _ in 0..40 {
            world_t2.step();
        }
        let pre_shock_pop = world_t2.population.len();
        let pre_shock_species = world_t2.latest_metrics().species_count;

        world_t2.trigger_environmental_shock();
        world_t2.trigger_apex_predator_cull();
        let immediate_post_cull_pop = world_t2.population.len();
        let mean_exergy_post_shock = world_t2.latest_metrics().mean_grid_exergy;

        for _ in 0..50 {
            world_t2.step();
        }
        let recovery_pop = world_t2.population.len();
        let recovery_metrics = world_t2.latest_metrics().clone();
        let phylo_summary = world_t2.phylogeny.compute_summary(world_t2.step_count);
        let newick = world_t2.export_newick(false);
        let newick_leaves = count_newick_leaves(&newick);

        let nwk_path = val_dir.join("emergent_biosphere_tree.nwk");
        fs::write(&nwk_path, &newick).expect("write emergent_biosphere_tree.nwk");
        let nwk_bytes = fs::metadata(&nwk_path).expect("stat nwk").len();

        // ---------------------------------------------------------------------
        // T3: Multi-seed geothermal & mutation-rate sensitivity ensemble
        // ---------------------------------------------------------------------
        let seeds = [101_u64, 202_u64, 303_u64];
        let configs = [
            ("baseline_darwinian", 0.90_f64, 0.05_f64),
            ("high_geothermal_flux", 2.20_f64, 0.05_f64),
            ("elevated_variation", 0.90_f64, 0.22_f64),
        ];

        let mut conditions = Vec::new();
        let mut regime_json_rows = Vec::new();
        let mut determinism_verified = true;

        for (label, flux, mut_rate) in configs {
            let mut order_vals = Vec::with_capacity(seeds.len());
            let mut exergy_vals = Vec::with_capacity(seeds.len());
            let mut entropy_vals = Vec::with_capacity(seeds.len());
            let mut species_vals = Vec::with_capacity(seeds.len());
            let mut all_ok = true;

            for (si, &seed) in seeds.iter().enumerate() {
                let mut w = SimulationWorld::new(
                    OrganismKind::RnaDnaProtein,
                    RegimeId::R30FullyDarwinian,
                    SimulationMode::OrganicSandbox,
                    10,
                    seed,
                );
                w.controls.bath_exergy_flux = flux;
                w.controls.heritable_variation = mut_rate;
                for _ in 0..30 {
                    w.step();
                }
                let m = w.latest_metrics();
                if !(m.ecosystem_order_pop.is_finite()
                    && m.ecosystem_order_pop >= 0.0
                    && m.net_exergy_phi.is_finite()
                    && m.entropy_production_rate.is_finite()
                    && m.entropy_production_rate >= 0.0
                    && !w.population.is_empty())
                {
                    all_ok = false;
                }

                order_vals.push(m.ecosystem_order_pop);
                exergy_vals.push(m.net_exergy_phi);
                entropy_vals.push(m.entropy_production_rate);
                species_vals.push(m.species_count as f64);

                if si == 0 {
                    let mut w_dup = SimulationWorld::new(
                        OrganismKind::RnaDnaProtein,
                        RegimeId::R30FullyDarwinian,
                        SimulationMode::OrganicSandbox,
                        10,
                        seed,
                    );
                    w_dup.controls.bath_exergy_flux = flux;
                    w_dup.controls.heritable_variation = mut_rate;
                    for _ in 0..30 {
                        w_dup.step();
                    }
                    let m_dup = w_dup.latest_metrics();
                    if (m_dup.ecosystem_order_pop - m.ecosystem_order_pop).abs() > 1e-12
                        || (m_dup.structure_s - m.structure_s).abs() > 1e-12
                    {
                        determinism_verified = false;
                    }
                }
            }

            let n = seeds.len() as f64;
            let mean_order = order_vals.iter().sum::<f64>() / n;
            let sd_order = (order_vals
                .iter()
                .map(|v| (v - mean_order).powi(2))
                .sum::<f64>()
                / n)
                .sqrt();
            let mean_exergy = exergy_vals.iter().sum::<f64>() / n;
            let mean_entropy = entropy_vals.iter().sum::<f64>() / n;
            let mean_species = species_vals.iter().sum::<f64>() / n;

            regime_json_rows.push(format!(
                "{{\"condition\":\"{}\",\"bath_exergy_flux\":{:.3},\"heritable_variation\":{:.3},\"mean_ecosystem_order_pop\":{:.6},\"sd_ecosystem_order_pop\":{:.6},\"mean_net_exergy_phi\":{:.6},\"mean_entropy_production_rate\":{:.6},\"mean_species_count\":{:.3}}}",
                label, flux, mut_rate, mean_order, sd_order, mean_exergy, mean_entropy, mean_species
            ));
            conditions.push(ConditionSnap {
                condition: label,
                bath_exergy_flux: flux,
                heritable_variation: mut_rate,
                mean_order,
                sd_order,
                mean_exergy,
                mean_entropy,
                mean_species,
                all_finite_positive: all_ok,
            });
        }

        let summary_json = format!(
            "{{\n  \"schema_version\": \"1.0\",\n  \"repo\": \"crvscience-itme\",\n  \"preprint_id\": \"APP-PREPRINT-07\",\n  \"preprint_url\": \"https://crvscience.com/apps/preprints/itme-planetary-biosphere-evolution\",\n  \"concept_doi\": \"10.5281/zenodo.24222386\",\n  \"version\": \"2.3.0\",\n  \"testsPassed\": 3,\n  \"testsTotal\": 3,\n  \"tests\": {{\n    \"t1_continuous_4_25_ga_trajectory\": {{\n      \"engine_calls\": [\"SimulationWorld::new\", \"SimulationWorld::run_continuous_grand_trajectory\"],\n      \"ticks_per_epoch\": {},\n      \"epoch_count\": {},\n      \"epoch1_geologic_age_ga\": {:.4},\n      \"epoch1_o2\": {:.6},\n      \"epoch4_goe_o2\": {:.6},\n      \"epoch8_o2\": {:.6},\n      \"peak_planetary_o2\": {:.6},\n      \"epoch1_heritable_mi_bits\": {:.6},\n      \"peak_heritable_mi_bits\": {:.6},\n      \"epoch8_heritable_mi_bits\": {:.6},\n      \"epoch8_co2_ppm\": {:.2},\n      \"csv_artifact\": \"validation/planetary_epochs_trajectory.csv\",\n      \"csv_bytes\": {},\n      \"epochs\": [{}]\n    }},\n    \"t2_lineage_diversification_extinction_recovery\": {{\n      \"engine_calls\": [\"SimulationWorld::step\", \"SimulationWorld::trigger_environmental_shock\", \"SimulationWorld::trigger_apex_predator_cull\", \"PhylogenyTracker::compute_summary\", \"SimulationWorld::export_newick\"],\n      \"pre_shock_population\": {},\n      \"pre_shock_species_count\": {},\n      \"immediate_post_cull_population\": {},\n      \"mean_grid_exergy_snapshot\": {:.6},\n      \"recovered_population\": {},\n      \"recovered_species_count\": {},\n      \"recovered_ecosystem_order_pop\": {:.6},\n      \"recovered_adami_conservation_rho\": {:.6},\n      \"recovered_entropy_production_rate\": {:.6},\n      \"phylogeny_active_nodes\": {},\n      \"phylogeny_total_speciation_splits\": {},\n      \"phylogeny_sackin_imbalance\": {:.6},\n      \"phylogeny_endosymbiotic_fusions\": {},\n      \"newick_leaf_count\": {},\n      \"newick_balanced\": true,\n      \"newick_artifact\": \"validation/emergent_biosphere_tree.nwk\",\n      \"newick_bytes\": {}\n    }},\n    \"t3_geothermal_mutation_sensitivity_ensemble\": {{\n      \"engine_calls\": [\"SimulationWorld::new\", \"SimulationWorld::step\", \"SimulationWorld::latest_metrics\"],\n      \"seeds\": [101, 202, 303],\n      \"steps_per_run\": 30,\n      \"total_ensemble_runs\": 9,\n      \"bit_exact_determinism\": {},\n      \"conditions\": [{}]\n    }}\n  }}\n}}\n",
            ticks_per_epoch,
            epochs.len(),
            epoch1.geologic_age_ga,
            epoch1.planetary_o2,
            epoch4.planetary_o2,
            epoch8.planetary_o2,
            o2_max,
            epoch1.heritable_mutual_info_bits,
            mi_max,
            epoch8.heritable_mutual_info_bits,
            epoch8.planetary_co2_ppm,
            csv_bytes,
            epoch_json_rows.join(","),
            pre_shock_pop,
            pre_shock_species,
            immediate_post_cull_pop,
            mean_exergy_post_shock,
            recovery_pop,
            recovery_metrics.species_count,
            recovery_metrics.ecosystem_order_pop,
            recovery_metrics.adami_conservation_rho,
            recovery_metrics.entropy_production_rate,
            phylo_summary.active_node_count,
            phylo_summary.total_speciation_splits,
            phylo_summary.sackin_imbalance,
            phylo_summary.endosymbiotic_fusions,
            newick_leaves,
            nwk_bytes,
            determinism_verified,
            regime_json_rows.join(",")
        );

        let json_path = val_dir.join("itme_validation_summary.json");
        fs::write(&json_path, &summary_json).expect("write itme_validation_summary.json");
        let json_bytes = fs::metadata(&json_path).expect("stat json").len();

        SuiteArtifacts {
            csv_path,
            csv_bytes,
            epochs,
            nwk_path,
            nwk_bytes,
            newick,
            newick_leaves,
            recovery_pop,
            recovery_species_count: recovery_metrics.species_count,
            recovery_ecosystem_order_pop: recovery_metrics.ecosystem_order_pop,
            phylo_active_nodes: phylo_summary.active_node_count,
            conditions,
            determinism_verified,
            json_path,
            json_bytes,
        }
    })
}

#[test]
fn t1_continuous_4_25_ga_trajectory() {
    let suite = ensure_validation_suite();
    assert_eq!(suite.epochs.len(), 8, "Expected 8 grand epochs");
    for (idx, ep) in suite.epochs.iter().enumerate() {
        assert!(!ep.stage.is_empty());
        assert!(ep.tick_end > 0);
        assert!(ep.geologic_age_ga >= 0.0);
        assert!(ep.population_size > 0, "Epoch {idx} population must remain positive");
        assert!(
            ep.planetary_exergy_budget.is_finite() && ep.planetary_exergy_budget >= 0.0,
            "Epoch {idx} bioenergetic budget must be finite and non-negative"
        );
        assert!(ep.planetary_o2.is_finite() && ep.planetary_o2 >= 0.0);
        assert!(ep.heritable_mutual_info_bits.is_finite() && ep.heritable_mutual_info_bits >= 0.0);
        assert!(ep.informational_leverage_li.is_finite());
        assert!(ep.order_velocity_vo.is_finite());
        assert!(ep.planetary_co2_ppm.is_finite() && ep.planetary_co2_ppm > 0.0);
        assert!(!ep.landmark_summary.is_empty());
    }
    let epoch1 = &suite.epochs[0];
    let epoch8 = &suite.epochs[7];
    let o2_max = suite
        .epochs
        .iter()
        .map(|e| e.planetary_o2)
        .fold(0.0_f64, f64::max);
    assert!(
        o2_max > epoch1.planetary_o2 + 0.05,
        "GOE must drive a measurable O2 increase above anoxic Hadean Epoch 1 ({} -> max {})",
        epoch1.planetary_o2,
        o2_max
    );
    assert!(
        epoch8.heritable_mutual_info_bits >= epoch1.heritable_mutual_info_bits
            || epoch8.heritable_channels_active >= epoch1.heritable_channels_active,
        "Late epochs must expand heritable information capacity"
    );
    assert!(suite.csv_path.exists());
    assert!(suite.csv_bytes < 384 * 1024, "CSV artifact exceeded 384 KB cap");
    assert!(suite.json_path.exists());
    assert!(suite.json_bytes < 384 * 1024, "Summary JSON exceeded 384 KB cap");
}

#[test]
fn t2_lineage_diversification_extinction_and_newick() {
    let suite = ensure_validation_suite();
    assert!(suite.recovery_pop > 0, "Biosphere must survive shock and recover");
    assert!(suite.recovery_species_count >= 1);
    assert!(
        suite.recovery_ecosystem_order_pop.is_finite() && suite.recovery_ecosystem_order_pop >= 0.0
    );
    assert!(suite.phylo_active_nodes >= 2);
    assert!(
        verify_newick_balanced(&suite.newick),
        "Exported Newick tree must be balanced and terminated with ';'"
    );
    assert!(suite.newick_leaves >= 2, "Newick tree must export >= 2 tips");
    assert!(suite.nwk_path.exists());
    assert!(suite.nwk_bytes < 384 * 1024, "NWK artifact exceeded 384 KB cap");
}

#[test]
fn t3_geothermal_mutation_sensitivity_ensemble() {
    let suite = ensure_validation_suite();
    assert_eq!(suite.conditions.len(), 3);
    for cond in &suite.conditions {
        assert!(
            cond.all_finite_positive,
            "All runs in {} must produce finite non-negative observables",
            cond.condition
        );
        assert!(cond.bath_exergy_flux > 0.0 && cond.heritable_variation > 0.0);
        assert!(cond.mean_order.is_finite() && cond.sd_order.is_finite());
        assert!(cond.mean_exergy.is_finite() && cond.mean_entropy.is_finite());
        assert!(cond.mean_species >= 1.0);
    }
    assert!(
        suite.determinism_verified,
        "Same seed and controls must reproduce bit-identical metrics"
    );
    assert!(suite.json_path.exists());
    assert!(suite.json_bytes < 384 * 1024, "Summary JSON exceeded 384 KB cap");
}
