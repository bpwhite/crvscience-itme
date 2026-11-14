//! Headless console mode (`--headless`) for the ITME v12.0 (`v1.2.0`) Anaerobic Mitochondrial
//! Reduction (`Hydrogenosomes` & `Mitosomes`), Apicoplast, Kleptoplasty & Serial Endosymbiosis Engine.

use crate::engine::{SimulationMode, SimulationWorld};
use crate::geometry::{render_ascii_3d, render_stage_atlas_svg, Camera3D};
use crate::model::TrophicGuild;
use crate::organisms::OrganismKind;
use crate::regimes::RegimeId;
use crate::spatial::render_ascii_spatial_grid;
use std::fs;
use std::io;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct ConsoleOptions {
    pub organisms: Vec<OrganismKind>,
    pub regimes: Vec<RegimeId>,
    pub steps_per_regime: usize,
    pub seed: u64,
    pub show_ascii_3d: bool,
    pub csv_path: Option<String>,
    pub svg_path: Option<String>,
    pub newick_path: Option<String>,
    pub phylo_svg_path: Option<String>,
    pub grand_trajectory: bool,
}

impl Default for ConsoleOptions {
    fn default() -> Self {
        Self {
            organisms: OrganismKind::ALL.to_vec(),
            regimes: RegimeId::ALL.to_vec(),
            steps_per_regime: 24,
            seed: 20250813,
            show_ascii_3d: true,
            csv_path: None,
            svg_path: None,
            newick_path: None,
            phylo_svg_path: None,
            grand_trajectory: false,
        }
    }
}

pub fn run_grand_trajectory_headless(opts: &ConsoleOptions) -> io::Result<()> {
    println!("===================================================================================================================================");
    println!(" ITME v23.0 (v2.3.0) — CONTINUOUS 4.25 Ga HADEAN-TO-GAIA MULTI-EPOCH GRAND TRAJECTORY BENCHMARK");
    println!("===================================================================================================================================");
    let mut world = SimulationWorld::new(
        OrganismKind::ArtificialIntelligence,
        RegimeId::R35Symbiogenetic,
        SimulationMode::OrganicSandbox,
        opts.steps_per_regime.max(12),
        opts.seed,
    );
    let summaries = world.run_continuous_grand_trajectory(opts.steps_per_regime.max(12));
    println!(
        " {:<22} {:>7} {:>6} {:>9} {:>10} {:>8} {:>9} {:>8}  {}",
        "Epoch", "Age(Ga)", "Chan", "MI(bits)", "Lev(L_I)", "O2", "CO2(ppm)", "Exergy", "Transition Landmark"
    );
    println!(" {}", "-".repeat(126));
    let mut csv_out = String::from(
        "epoch,age_ga,tick_end,population,channels,mutual_info_bits,leverage_li,o2,co2_ppm,exergy,landmark\n",
    );
    for s in &summaries {
        println!(
            " {:<22} {:>7.3} {:>6} {:>9.2} {:>10.3} {:>8.2} {:>9.0} {:>8.2}  {}",
            s.epoch.code(),
            s.geologic_age_ga,
            s.heritable_channels_active,
            s.heritable_mutual_info_bits,
            s.informational_leverage_li,
            s.planetary_o2,
            s.planetary_co2_ppm,
            s.planetary_exergy_budget,
            s.landmark_summary
        );
        csv_out.push_str(&format!(
            "{},{:.3},{},{},{},{:.3},{:.3},{:.3},{:.1},{:.3},\"{}\"\n",
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
            s.landmark_summary
        ));
    }
    if let Some(ref cp) = opts.csv_path {
        fs::write(Path::new(cp), csv_out)?;
    }
    if let Some(ref np) = opts.newick_path {
        fs::write(Path::new(np), world.export_newick(true))?;
    }
    if let Some(ref sp) = opts.phylo_svg_path {
        fs::write(Path::new(sp), world.export_phylogram_svg())?;
    }
    Ok(())
}

pub fn run_headless(opts: &ConsoleOptions) -> io::Result<()> {
    if opts.grand_trajectory {
        return run_grand_trajectory_headless(opts);
    }
    println!("=========================================================================================================================================================================");
    println!(" INFORMATION-THEORETIC MODEL OF EVOLUTION (ITME v23.0 / v2.3.0) — ANAEROBIC MROs (HYDROGENOSOMES & MITOSOMES), APICOPLASTS, KLEPTOPLASTY & 1°/2°/3° ENDOSYMBIOSIS");
    println!(
        " Organisms: {} | Regimes: {} | Steps/Regime: {} | PRNG Seed: {}",
        opts.organisms.len(),
        opts.regimes.len(),
        opts.steps_per_regime,
        opts.seed
    );
    println!("=========================================================================================================================================================================");

    let mut csv_buf = String::new();
    if opts.csv_path.is_some() {
        csv_buf.push_str("organism,regime_code,regime_name,step,S_raw_bits,S_eff_bits,Adami_rho_within,C_emergent,L_I,S_T_minus,Order_O,O_pop,v_O,Phi_net,Q_landauer,Grid_Exergy,Byproduct_Pool,Niche_Capital,Producers,Syntrophs,Predators,Aero_Mitos,Hydrogenosomes,Mitosomes,Amiteukaryotes,Mean_MRO_Paths,Hydro_Transitions,MitoS_Transitions,MRO_Losses,Apico_Parasites,Mean_Apico_Paths,Apico_Transitions,Apico_Losses,Klepto_Hosts,Mean_KP_Integ,Mean_Klepto_HGT,Klepto_HGT_Events,Klepto_Fixations,Mosaic_Donors,Holobionts,Enslaved_Holobionts,Secondary_Holobionts,Tertiary_Holobionts,Active_Nucleomorphs,Mean_NM_Len,NM_Loss_Events,Serial_Fusions,Mean_EGT_Frac,Mean_Organelle_Len,EGT_Events,Mean_TL,FoodWeb_H_bits,Trophic_Flux,CrossFeed_Flux,Species_Count,F_ST,TMRCA_steps,Sackin_Imbalance,Speciation_Splits,Guild_Transitions,Endosymbiotic_Fusions,Active_Reticulations,MI_bits,Price_Cov,Price_Trans,Pop,Emergent_Dc\n");
    }

    let mut newick_buf = String::new();
    let mut last_phylo_svg = String::new();
    let mut svg_panels = Vec::new();
    let cam = Camera3D::default();

    for &org in &opts.organisms {
        println!("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        println!(" MODEL ORGANISM: {}", org.title());
        println!(" Contextual Alphabet: {}", org.alphabet_description());
        println!("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        println!(
            " {:<10} {:<23} {:>6} {:>6} {:>11} {:>9} {:>7} {:>7} {:>5} {:>5}  {:<18}",
            "Regime",
            "Emergent Stage",
            "Order",
            "O_pop",
            "P/Syn/Hol/X",
            "Hy/Mi/Am0",
            "Apic/P",
            "KP/Fix",
            "Fuse",
            "MRO>0",
            "Emergent Dc"
        );
        println!(" {}", "-".repeat(162));

        let mut world = SimulationWorld::new(
            org,
            opts.regimes[0],
            SimulationMode::OrganicSandbox,
            opts.steps_per_regime,
            opts.seed,
        );

        for (r_idx, &regime) in opts.regimes.iter().enumerate() {
            let preserve = r_idx > 2;
            world.apply_regime_preset(regime, preserve);

            for _ in 0..opts.steps_per_regime {
                world.step();
                if opts.csv_path.is_some() {
                    let m = world.latest_metrics();
                    csv_buf.push_str(&format!(
                        "{},{},\"{}\",{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.5},{:.4},{:.4},{:.4},{:.4},{:.4},{},{},{},{},{},{},{},{:.2},{},{},{},{},{:.2},{},{},{},{:.3},{:.2},{},{},{:.2},{},{},{},{},{},{:.2},{},{},{:.4},{:.2},{},{:.3},{:.3},{:.4},{:.4},{},{:.4},{},{:.4},{},{},{},{},{:.3},{:.5},{:.5},{},\"{}\"\n",
                        org.slug(),
                        regime.code(),
                        regime.name(),
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
                        m.diagnostic
                    ));
                }
            }

            let m = world.latest_metrics().clone();
            let exact_match = world.classified_regime == regime;
            let dc_mark = if exact_match { "OK" } else { "~" };
            let guild_str = format!(
                "{}/{}/{}/{}",
                m.producer_count, m.syntroph_count, m.holobiont_count, m.predator_count
            );
            let mro_str = format!(
                "{}/{}/{}",
                m.hydrogenosome_count, m.mitosome_count, m.amiteukaryote_count
            );
            let apic_str = format!(
                "{}/{:.1}",
                m.apicoplast_parasite_count, m.mean_apicoplast_pathways
            );
            let kp_str = format!("{}/{}", m.kleptoplastic_host_count, m.klepto_fixation_events);

            println!(
                " {:<10} {:<23} {:>6.2} {:>6.2} {:>11} {:>9} {:>7} {:>7} {:>5} {:>5}  {} [{}]",
                world.classified_regime.code(),
                truncate_str(world.classified_regime.name(), 23),
                m.order_o,
                m.ecosystem_order_pop,
                guild_str,
                mro_str,
                apic_str,
                kp_str,
                m.endosymbiotic_fusions,
                m.mro_complete_loss_events,
                m.diagnostic,
                dc_mark
            );

            if opts.newick_path.is_some() && regime >= RegimeId::R30FullyDarwinian {
                newick_buf.push_str(&format!(
                    "# {} | {} | step={} | Fusions(#H)={} | Mito->Hydro={} | Hydro->MitoS={} | MRO->0={} | Holo->Apico={}\n{}\n\n",
                    org.title(),
                    regime.full_title(),
                    m.step,
                    m.endosymbiotic_fusions,
                    m.hydrogenosome_transition_events,
                    m.mitosome_transition_events,
                    m.mro_complete_loss_events,
                    m.apicoplast_transition_events,
                    world.export_newick(true)
                ));
            }
            if opts.phylo_svg_path.is_some()
                && (regime == RegimeId::R35Symbiogenetic
                    || (last_phylo_svg.is_empty() && regime >= RegimeId::R30FullyDarwinian))
            {
                last_phylo_svg = world.export_phylogram_svg();
            }

            let scene = world.current_scene();
            if opts.svg_path.is_some() {
                let panel_title = format!("{} | {}", org.short_title(), regime.code());
                let metric_line = format!(
                    "O={:.2} Hydro={} MitoS={} Apico={}",
                    m.order_o,
                    m.hydrogenosome_count,
                    m.mitosome_count,
                    m.apicoplast_parasite_count
                );
                svg_panels.push((panel_title, metric_line, scene.clone()));
            }

            let highlight_stage = opts.regimes.len() <= 3
                || matches!(
                    regime,
                    RegimeId::R00Static
                        | RegimeId::R20ProtoDarwinian
                        | RegimeId::R30FullyDarwinian
                        | RegimeId::R35Symbiogenetic
                        | RegimeId::R30PseudoSoftInheritance
                );

            if opts.show_ascii_3d && highlight_stage {
                println!(
                    "\n   [3D Structure, Spatial Web & MRO/Apicoplast Network] {} @ {} — {}",
                    org.short_title(),
                    regime.full_title(),
                    org.regime_narrative(regime)
                );
                if let Some(best) = world.best_unit() {
                    println!(
                        "   Champion @ ({}, {}) [{} {}° | MRO:{} | Sp {}]: {}",
                        best.grid_x,
                        best.grid_y,
                        best.guild.short_code(),
                        best.endo_tier,
                        best.mro_class.short_code(),
                        best.species_id,
                        best.formatted_blueprint(10)
                    );
                }
                let ascii = render_ascii_3d(&scene, &cam, 56, 13);
                for line in ascii.lines() {
                    println!("   {}", line);
                }
                if regime >= RegimeId::R30FullyDarwinian {
                    let coords: Vec<(usize, usize, usize, TrophicGuild)> = world
                        .population
                        .iter()
                        .map(|u| (u.grid_x, u.grid_y, u.species_id, u.guild))
                        .collect();
                    let map_ascii = render_ascii_spatial_grid(&world.grid, &coords);
                    print!("{}", map_ascii);
                    let nhx_preview = world.export_newick(true);
                    println!(
                        "   eNewick+MRO Network: {}",
                        truncate_str(&nhx_preview, 122)
                    );
                }
            }
        }

        if opts.newick_path.is_some() && newick_buf.is_empty() {
            newick_buf.push_str(&world.export_newick(true));
            newick_buf.push('\n');
        }
        if opts.phylo_svg_path.is_some() && last_phylo_svg.is_empty() {
            last_phylo_svg = world.export_phylogram_svg();
        }
    }

    if let Some(ref path) = opts.csv_path {
        write_file_creating_dirs(path, &csv_buf)?;
        println!(
            "\n[Export] Wrote ITME v12.0 anaerobic MRO, apicoplast & symbiogenesis CSV to: {}",
            path
        );
    }

    if let Some(ref path) = opts.svg_path {
        let svg_doc = render_stage_atlas_svg(
            "ITME v12.0 Anaerobic MROs (Hydrogenosomes & Mitosomes), Apicoplasts & Endosymbiosis Atlas",
            &svg_panels,
        );
        write_file_creating_dirs(path, &svg_doc)?;
        println!(
            "[Export] Wrote ITME v12.0 3D/2D structural stage atlas SVG to: {}",
            path
        );
    }

    if let Some(ref path) = opts.newick_path {
        write_file_creating_dirs(path, &newick_buf)?;
        println!(
            "[Export] Wrote Extended Newick (eNewick #H + MRO/APICO/KLEPTO NHX) networks to: {}",
            path
        );
    }

    if let Some(ref path) = opts.phylo_svg_path {
        write_file_creating_dirs(path, &last_phylo_svg)?;
        println!(
            "[Export] Wrote anastomosing SVG phylogram (Mito->Hydro, Hydro->MitoS, MRO->0 & Holo->Apico) to: {}",
            path
        );
    }

    println!("\nHeadless ITME v12.0 (v1.2.0) anaerobic MRO, apicoplast & symbiogenesis simulation completed.");
    Ok(())
}

fn write_file_creating_dirs(path: &str, contents: &str) -> io::Result<()> {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    fs::write(path, contents)
}

fn truncate_str(s: &str, max_chars: usize) -> String {
    if s.chars().count() <= max_chars {
        s.to_string()
    } else {
        let mut t: String = s.chars().take(max_chars.saturating_sub(2)).collect();
        t.push_str("..");
        t
    }
}
