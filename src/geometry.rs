//! 2D and 3D procedural geometry generation, 3D perspective camera projection,
//! terminal ASCII/ANSI rendering, and standalone SVG structural stage atlas export (ITME v12.0 / `v1.2.0`).

use crate::model::{AbiogenesisStage, AgricultureState, ChiralityCodeStage, CulturalState, EusocialState, ExosomaticState, GenomeSubstrateStage, MachineSynthesisState, MroClass, MulticellularState, PlanetaryState, TrophicGuild};
use crate::organisms::{ModeledUnit, OrganismKind};
use crate::regimes::RegimeId;
use std::f64::consts::{PI, TAU};
use std::fmt::Write as _;

pub type Rgb = [u8; 3];

#[derive(Copy, Clone, Debug, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }

    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }

    pub fn lerp(self, o: Self, t: f64) -> Self {
        Self::new(
            self.x + (o.x - self.x) * t,
            self.y + (o.y - self.y) * t,
            self.z + (o.z - self.z) * t,
        )
    }
}

#[derive(Clone, Debug)]
pub struct Node3D {
    pub pos: Vec3,
    pub radius: f64,
    pub color: Rgb,
    pub label: String,
    pub glow: f64,
}

#[derive(Clone, Debug)]
pub struct Edge3D {
    pub a: usize,
    pub b: usize,
    pub width: f64,
    pub color: Rgb,
    pub dashed: bool,
}

#[derive(Clone, Debug)]
pub struct Node2D {
    pub x: f64,
    pub y: f64,
    pub radius: f64,
    pub color: Rgb,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct Edge2D {
    pub a: usize,
    pub b: usize,
    pub width: f64,
    pub color: Rgb,
    pub directed: bool,
}

#[derive(Clone, Debug)]
pub struct StructureScene {
    pub kind: OrganismKind,
    pub regime: RegimeId,
    pub subtitle: String,
    pub nodes_3d: Vec<Node3D>,
    pub edges_3d: Vec<Edge3D>,
    pub nodes_2d: Vec<Node2D>,
    pub edges_2d: Vec<Edge2D>,
}

#[derive(Copy, Clone, Debug)]
pub struct Camera3D {
    pub yaw: f64,
    pub pitch: f64,
    pub zoom: f64,
    pub distance: f64,
}

impl Default for Camera3D {
    fn default() -> Self {
        Self {
            yaw: 0.58,
            pitch: 0.34,
            zoom: 1.0,
            distance: 4.8,
        }
    }
}

impl Camera3D {
    pub fn project(&self, p: Vec3, width: f64, height: f64) -> (f64, f64, f64, f64) {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();

        let x1 = p.x * cy + p.z * sy;
        let z1 = -p.x * sy + p.z * cy;
        let y2 = p.y * cp - z1 * sp;
        let z2 = p.y * sp + z1 * cp;

        let depth = (self.distance - z2).max(0.65);
        let focal = width.min(height) * 0.62 * self.zoom;
        let persp = focal / depth;

        let sx = width * 0.5 + x1 * persp;
        let sy_px = height * 0.5 - y2 * persp;
        (sx, sy_px, depth, persp)
    }
}

pub fn build_structure_scene(
    unit: &ModeledUnit,
    regime: RegimeId,
    population_size: usize,
    phase_time: f64,
) -> StructureScene {
    let mut scene = match unit.kind {
        OrganismKind::RnaDnaProtein => {
            build_rna_dna_protein_scene(unit, regime, population_size, phase_time)
        }
        OrganismKind::HominidPopulation => {
            build_hominid_scene(unit, regime, population_size, phase_time)
        }
        OrganismKind::ArtificialIntelligence => {
            build_ai_scene(unit, regime, population_size, phase_time)
        }
    };
    append_trophic_guild_ornaments(&mut scene, unit, phase_time);
    append_mro_ornaments(&mut scene, unit, phase_time);
    append_multicellular_ornaments(&mut scene, unit, phase_time);
    append_superorganism_ornaments(&mut scene, unit, phase_time);
    append_agriculture_ornaments(&mut scene, unit, phase_time);
    append_cultural_ornaments(&mut scene, unit, phase_time);
    append_exosomatic_ornaments(&mut scene, unit, phase_time);
    append_machine_synthesis_ornaments(&mut scene, unit, phase_time);
    append_planetary_biosphere_ornaments(&mut scene, unit, phase_time);
    append_abiogenesis_protocell_ornaments(&mut scene, unit, phase_time);
    append_homochiral_minihelix_ornaments(&mut scene, unit, phase_time);
    append_dna_double_helix_and_capsid_ornaments(&mut scene, unit, phase_time);
    scene
}

/// Renders 3D structural ornaments for **Mitochondrion-Related Organelles (`MRO`s)**:
/// - `Hydrogenosome`: double-membrane teal-emerald fermentative organelle (`mtDNA = 0`) emitting `H2` efflux spheres!
/// - `Mitosome`: compact double-membrane minimalist vesicle housing a glowing `[2Fe-2S]` `ISC` cluster!
/// - `Amiteukaryote`: cytosolic `SUF` `[Fe-S]` scaffold nodes (`MRO -> 0`, *Monocercomonoides*)!
/// Renders 3D **Clonal Multicellular Rosette / Volvox Spheroid** ornaments:
/// outer ring of motile somatic cells + inner sequestered quiescent germline (`gonidia`)
/// + apoptotic cheater-policing / propagule-scission sparks!
/// Renders 3D **Eusocial Superorganism Colony Architecture**:
/// inner royal **Queen Germline** chamber + orbital ring of sterile **Foraging Workers**
/// + outer armored **Soldier Caste** perimeter + **Worker Kin-Policing** arcs!
/// Renders 3D **Obligate Agricultural Fungus-Garden & Cuticular Actinobacterial Shield**:
/// subterranean polyploid **Gongylidia** crop clusters + cyan *Pseudonocardia* antimicrobial crypt halo!
/// Renders 3D **Cumulative Cultural Knowledge Halo & Symbolic Teaching Lattice**:
/// overhead cerulean **Neural Plasticity Crown** + orbital gold **Symbolic Repertoire Nodes** (`T1..T6`)!
/// Renders 3D **Exosomatic Written Archive Monoliths & Peer-Verification Ring**:
/// external crystalline archive tablets (`[20, 184, 166]` teal) + platinum replication audit beams!
/// Renders 3D **Neuro-Symbolic Computational Symbiont & Formal Proof-Checker Polyhedron**:
/// violet-cobalt **Sampler-Evaluator Octahedral Nodes** (`[139, 92, 246]`) + cyan empirical-anchor beams!
/// Renders 3D **Planetary Atmospheric O2/CO2 Biosignature & Gaia-2.0 Stewardship Shell**:
/// equatorial orbital planetary nodes (`[16, 185, 129]` emerald-cyan in `HomeostaticAwakening`,
/// `[249, 115, 22]` amber-crimson when planetary boundaries are transgressed)!
/// Renders 3D **Hadean Protocell Amphiphilic Micelle & RNA Ribozyme Hairpin** geometry
/// during pre-LUCA abiogenesis stages (`MonomerSoup`, `RAFProtocell`, `RibozymeHypercycle`, `ProgenotePool`).
/// Renders 3D **Right-Handed `D-RNA` + `L-Aminoacyl` Acceptor Minihelix & `SGC` Codon Prism**
/// during pre-SGC homochiral symmetry-breaking and code-crystallization transitions.
/// Renders 3D **B-Form Deoxyribose Double-Helix & Primordial Icosahedral Virion** geometry
/// during active `RNR` `RNA -> U-DNA -> T-DNA` and virosphere coevolution transitions!
fn append_dna_double_helix_and_capsid_ornaments(
    scene: &mut StructureScene,
    unit: &ModeledUnit,
    phase_time: f64,
) {
    if unit.genome_substrate_stage == GenomeSubstrateStage::CoevolvedCellularVirosphere {
        return;
    }
    let col = match unit.genome_substrate_stage {
        GenomeSubstrateStage::LabileRnaChromosome => [244, 114, 182],
        GenomeSubstrateStage::ViralModifiedUDna => [34, 211, 238],
        GenomeSubstrateStage::DeaminationRepairableTDna => [52, 211, 153],
        GenomeSubstrateStage::CoevolvedCellularVirosphere => [52, 211, 153],
    };
    let rungs = 6usize;
    let base_idx = scene.nodes_3d.len();
    for k in 0..rungs {
        let t = (k as f64) / ((rungs - 1) as f64);
        let ang = t * TAU * 0.85 + phase_time * 0.40;
        let y = -0.80 + 1.60 * t;
        let a_idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 0.72, y, ang.sin() * 0.72),
            radius: 0.065,
            color: col,
            label: if k == 0 {
                format!(
                    "{} ({:.0}kb)",
                    unit.genome_substrate_stage.short_code(),
                    unit.chromosome_capacity_kb
                )
            } else {
                String::new()
            },
            glow: 0.92,
        });
        let b_idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(-ang.cos() * 0.72, y, -ang.sin() * 0.72),
            radius: 0.065,
            color: col,
            label: String::new(),
            glow: 0.92,
        });
        scene.edges_3d.push(Edge3D {
            a: a_idx,
            b: b_idx,
            width: 1.1,
            color: col,
            dashed: unit.genome_substrate_stage == GenomeSubstrateStage::LabileRnaChromosome,
        });
        if k > 0 {
            scene.edges_3d.push(Edge3D {
                a: base_idx + 2 * (k - 1),
                b: a_idx,
                width: 1.25,
                color: col,
                dashed: false,
            });
        }
    }
}

fn append_homochiral_minihelix_ornaments(
    scene: &mut StructureScene,
    unit: &ModeledUnit,
    phase_time: f64,
) {
    if unit.chirality_code_stage == ChiralityCodeStage::UniversalErrorMinimizedSgc {
        return;
    }
    let col = match unit.chirality_code_stage {
        ChiralityCodeStage::RacemicCrossInhibited => [251, 113, 133],
        ChiralityCodeStage::HomochiralRnaSymmetryBroken => [192, 132, 252],
        ChiralityCodeStage::StereoselectiveMinihelixAcylated => [56, 189, 248],
        ChiralityCodeStage::AmbiguousStatisticalCode => [250, 204, 21],
        ChiralityCodeStage::UniversalErrorMinimizedSgc => [52, 211, 153],
    };
    let n_pts = 7usize;
    let base_idx = scene.nodes_3d.len();
    for k in 0..n_pts {
        let t = (k as f64) / ((n_pts - 1) as f64);
        let ang = t * TAU + phase_time * 0.45;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 0.95, -0.75 + 1.50 * t, ang.sin() * 0.95),
            radius: 0.075,
            color: col,
            label: if k == 0 {
                format!(
                    "{} (MP={:.0}%)",
                    unit.chirality_code_stage.short_code(),
                    unit.genetic_code_mp * 100.0
                )
            } else {
                String::new()
            },
            glow: 0.94,
        });
        if k > 0 {
            scene.edges_3d.push(Edge3D {
                a: base_idx + k - 1,
                b: idx,
                width: 1.35,
                color: col,
                dashed: unit.chirality_code_stage == ChiralityCodeStage::RacemicCrossInhibited,
            });
        }
    }
}

fn append_abiogenesis_protocell_ornaments(
    scene: &mut StructureScene,
    unit: &ModeledUnit,
    phase_time: f64,
) {
    if unit.abiogenesis_stage == AbiogenesisStage::LucaCellularLineage {
        return;
    }
    let col = match unit.abiogenesis_stage {
        AbiogenesisStage::PrebioticMonomerSoup => [244, 114, 182],
        AbiogenesisStage::AutocatalyticRafProtocell => [251, 191, 36],
        AbiogenesisStage::RnaRibozymeHypercycle => [56, 189, 248],
        AbiogenesisStage::ProgenoteCommunity => [52, 211, 153],
        AbiogenesisStage::LucaCellularLineage => [52, 211, 153],
    };
    let n_pts = 6usize;
    let base_idx = scene.nodes_3d.len();
    for k in 0..n_pts {
        let ang = (k as f64) * TAU / (n_pts as f64) - phase_time * 0.35;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 1.42, 0.32 * (ang * 2.0).sin(), ang.sin() * 1.42),
            radius: 0.08,
            color: col,
            label: if k == 0 {
                format!("{} ({}nt)", unit.abiogenesis_stage.short_code(), unit.rna_oligomer_length_nt)
            } else {
                String::new()
            },
            glow: 0.92,
        });
        scene.edges_3d.push(Edge3D {
            a: base_idx + ((k + n_pts - 1) % n_pts),
            b: idx,
            width: 1.2,
            color: col,
            dashed: unit.abiogenesis_stage == AbiogenesisStage::PrebioticMonomerSoup,
        });
    }
}

fn append_planetary_biosphere_ornaments(
    scene: &mut StructureScene,
    unit: &ModeledUnit,
    phase_time: f64,
) {
    if unit.planetary_epoch == PlanetaryState::OxygenatedPhanerozoic
        && unit.cumulative_stewardship_actions == 0
    {
        return;
    }
    let col = match unit.planetary_epoch {
        PlanetaryState::HomeostaticAwakening => [52, 211, 153],
        PlanetaryState::OxygenatedPhanerozoic => [56, 189, 248],
        PlanetaryState::AnoxicArchean => [168, 85, 247],
        PlanetaryState::AnthropoceneDisequilibrium => [249, 115, 22],
    };
    let ring_n = 8usize;
    let base_idx = scene.nodes_3d.len();
    for r in 0..ring_n {
        let ang = (r as f64) * TAU / (ring_n as f64) + phase_time * 0.22;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 3.18, 0.0, ang.sin() * 3.18),
            radius: 0.09,
            color: col,
            label: if r == 0 {
                format!("Biosphere ({})", unit.planetary_epoch.short_code())
            } else {
                String::new()
            },
            glow: 0.96,
        });
        scene.edges_3d.push(Edge3D {
            a: base_idx + ((r + ring_n - 1) % ring_n),
            b: idx,
            width: 1.25,
            color: col,
            dashed: unit.planetary_epoch == PlanetaryState::AnthropoceneDisequilibrium,
        });
    }
}

fn append_machine_synthesis_ornaments(
    scene: &mut StructureScene,
    unit: &ModeledUnit,
    phase_time: f64,
) {
    if !unit.machine_synthesis_state.has_computational_symbiont() {
        return;
    }
    let sym_n = match unit.machine_synthesis_state {
        MachineSynthesisState::MachineVerifiedSymbiosis => 6usize,
        MachineSynthesisState::UnverifiedRecursiveGen => 4usize,
        MachineSynthesisState::AlgorithmicMediation => 3usize,
        MachineSynthesisState::ManualArchival => 0usize,
    };
    let base_idx = scene.nodes_3d.len();
    for m in 0..sym_n {
        let ang = (m as f64) * TAU / (sym_n as f64) - phase_time * 0.92;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 2.25, 1.88 + 0.12 * (ang * 2.0).cos(), ang.sin() * 2.25),
            radius: 0.11,
            color: if unit.machine_synthesis_state.is_formally_verified() {
                [139, 92, 246]
            } else {
                [244, 63, 94]
            },
            label: if m == 0 {
                format!(
                    "{} ({:.1}b Proofs, {:.0}% Tails)",
                    unit.machine_synthesis_state.short_code(),
                    unit.verified_theorem_bits(),
                    unit.distribution_tail_diversity * 100.0
                )
            } else {
                String::new()
            },
            glow: 1.0,
        });
        if m > 0 {
            scene.edges_3d.push(Edge3D {
                a: base_idx + m - 1,
                b: idx,
                width: 1.45,
                color: [56, 189, 248],
                dashed: !unit.machine_synthesis_state.is_formally_verified(),
            });
        }
    }
}

fn append_exosomatic_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    if !unit.exosomatic_state.is_archival() {
        return;
    }
    let arch_n = match unit.exosomatic_state {
        ExosomaticState::PeerVerifiedInstitution => 8usize,
        ExosomaticState::WrittenArchive => 6usize,
        ExosomaticState::MaterialNotation => 3usize,
        ExosomaticState::EndosomaticBound => 0usize,
    };
    let base_idx = scene.nodes_3d.len();
    for a in 0..arch_n {
        let ang = (a as f64) * TAU / (arch_n as f64) + phase_time * 0.28;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 2.86, 1.05, ang.sin() * 2.86),
            radius: 0.105,
            color: [20, 184, 166],
            label: if a == 0 {
                format!(
                    "{} ({:.1}b Archive, -{:.0}% Brain Cost)",
                    unit.exosomatic_state.short_code(),
                    unit.exosomatic_archive_bits(),
                    unit.cognitive_offloading_fraction * 100.0
                )
            } else {
                String::new()
            },
            glow: 1.0,
        });
        if a > 0 {
            scene.edges_3d.push(Edge3D {
                a: base_idx + a - 1,
                b: idx,
                width: 1.4,
                color: if unit.exosomatic_state == ExosomaticState::PeerVerifiedInstitution {
                    [248, 250, 252]
                } else {
                    [20, 184, 166]
                },
                dashed: false,
            });
        }
    }
}

fn append_cultural_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    if !unit.cultural_state.is_cultural() {
        return;
    }
    let crown_idx = scene.nodes_3d.len();
    scene.nodes_3d.push(Node3D {
        pos: Vec3::new(0.0, 1.42, 0.0),
        radius: 0.145,
        color: [56, 189, 248],
        label: format!(
            "{} (T{}, {:.0}% Neural)",
            unit.cultural_state.short_code(),
            unit.cultural_ratchet_tier,
            unit.neural_learning_capacity * 100.0
        ),
        glow: 1.0,
    });

    let sym_n = (unit.cultural_ratchet_tier as usize * 2).clamp(3, 10);
    for s in 0..sym_n {
        let ang = (s as f64) * TAU / (sym_n as f64) - phase_time * 0.72;
        let s_idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 1.72, 1.32 + 0.10 * (ang * 2.0).sin(), ang.sin() * 1.72),
            radius: 0.088,
            color: [250, 204, 21],
            label: if s == 0 {
                format!(
                    "Cultural Corpus ({:.1} bits, Deme {})",
                    unit.cultural_repertoire_bits(),
                    unit.deme_id
                )
            } else {
                String::new()
            },
            glow: 0.98,
        });
        scene.edges_3d.push(Edge3D {
            a: crown_idx,
            b: s_idx,
            width: 1.3,
            color: [250, 204, 21],
            dashed: true,
        });
    }
}

fn append_agriculture_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    if !unit.agriculture_state.is_farming() {
        return;
    }
    let garden_n = match unit.agriculture_state {
        AgricultureState::LeafCutterAgroEcosystem => 8usize,
        AgricultureState::HigherDomesticatedCrop => 6usize,
        AgricultureState::LowerFacultativeFarming => 4usize,
        AgricultureState::HunterForager => 0usize,
    };
    let garden_base = scene.nodes_3d.len();
    for g in 0..garden_n {
        let ang = (g as f64) * TAU / (garden_n as f64) + phase_time * 0.42;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 1.48, -0.34 + 0.08 * (ang * 2.0).cos(), ang.sin() * 1.48),
            radius: 0.11,
            color: [163, 230, 53],
            label: if g == 0 {
                format!(
                    "{} ({}n, {:.0}% Gongylidia)",
                    unit.agriculture_state.short_code(),
                    unit.cultivar_polyploidy,
                    unit.gongylidia_yield * 100.0
                )
            } else {
                String::new()
            },
            glow: 0.96,
        });
        if g > 0 {
            scene.edges_3d.push(Edge3D {
                a: garden_base + g - 1,
                b: idx,
                width: 1.35,
                color: [163, 230, 53],
                dashed: false,
            });
        }
    }
    if unit.defensive_actinobacteria_cover >= 0.20 {
        for c in 0..4 {
            let ca = (c as f64) * TAU / 4.0 + phase_time * 0.68;
            let c_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(ca.cos() * 2.62, -0.16, ca.sin() * 2.62),
                radius: 0.085,
                color: [34, 211, 238],
                label: if c == 0 {
                    format!(
                        "Pseudonocardia Shield ({:.0}%)",
                        unit.defensive_actinobacteria_cover * 100.0
                    )
                } else {
                    String::new()
                },
                glow: 0.98,
            });
            scene.edges_3d.push(Edge3D {
                a: garden_base + (c % garden_n.max(1)),
                b: c_idx,
                width: 1.25,
                color: [34, 211, 238],
                dashed: true,
            });
        }
    }
}

fn append_superorganism_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    if !unit.eusocial_state.is_social() {
        return;
    }
    let queen_idx = scene.nodes_3d.len();
    scene.nodes_3d.push(Node3D {
        pos: Vec3::new(0.0, 0.36, 0.0),
        radius: 0.16,
        color: [245, 158, 11],
        label: format!("Queen Germline ({})", unit.eusocial_state.short_code()),
        glow: 1.0,
    });

    let worker_n = (unit.colony_size as usize / 8).clamp(4, 10);
    let worker_base = scene.nodes_3d.len();
    for k in 0..worker_n {
        let ang = (k as f64) * TAU / (worker_n as f64) + phase_time * 0.55;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 2.05, 0.25 + 0.12 * (ang * 3.0).sin(), ang.sin() * 2.05),
            radius: 0.095,
            color: [16, 185, 129],
            label: if k == 0 {
                format!(
                    "Sterile Workers ({:.0}%)",
                    unit.sterile_worker_fraction * 100.0
                )
            } else {
                String::new()
            },
            glow: 0.92,
        });
        scene.edges_3d.push(Edge3D {
            a: queen_idx,
            b: idx,
            width: 1.2,
            color: [245, 158, 11],
            dashed: true,
        });
    }

    if unit.soldier_fraction >= 0.05 {
        for s in 0..4 {
            let sa = (s as f64) * TAU / 4.0 - phase_time * 0.35;
            let s_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(sa.cos() * 2.38, 0.52, sa.sin() * 2.38),
                radius: 0.12,
                color: [239, 68, 68],
                label: if s == 0 {
                    format!("Soldier Caste ({:.0}%)", unit.soldier_fraction * 100.0)
                } else {
                    String::new()
                },
                glow: 0.95,
            });
            scene.edges_3d.push(Edge3D {
                a: worker_base + (s % worker_n),
                b: s_idx,
                width: 1.5,
                color: [239, 68, 68],
                dashed: false,
            });
        }
    }
}

fn append_multicellular_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    if !unit.multicell_state.is_multicellular() {
        return;
    }
    let outer_n = (unit.cell_count as usize).clamp(6, 12);
    let soma_col: Rgb = if unit.multicell_state == MulticellularState::GermSomaMulticellular {
        [56, 189, 248]
    } else if unit.multicell_state == MulticellularState::UndifferentiatedClonal {
        [129, 140, 248]
    } else {
        [249, 115, 22]
    };

    let soma_base = scene.nodes_3d.len();
    for k in 0..outer_n {
        let ang = (k as f64) * TAU / (outer_n as f64) + phase_time * 0.42;
        let idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * 1.72, -0.18 + 0.16 * (ang * 2.0).cos(), ang.sin() * 1.72),
            radius: 0.105,
            color: soma_col,
            label: if k == 0 {
                format!(
                    "{} ({} cells, {:.0}% Soma)",
                    unit.multicell_state.short_code(),
                    unit.cell_count,
                    unit.soma_fraction * 100.0
                )
            } else {
                String::new()
            },
            glow: 0.92,
        });
        if k > 0 {
            scene.edges_3d.push(Edge3D {
                a: idx - 1,
                b: idx,
                width: 1.5,
                color: soma_col,
                dashed: false,
            });
        }
    }
    scene.edges_3d.push(Edge3D {
        a: soma_base + outer_n - 1,
        b: soma_base,
        width: 1.5,
        color: soma_col,
        dashed: false,
    });

    if unit.multicell_state == MulticellularState::GermSomaMulticellular {
        for g in 0..3 {
            let ga = (g as f64) * TAU / 3.0 - phase_time * 0.30;
            let g_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(ga.cos() * 0.28, 0.05, ga.sin() * 0.28),
                radius: 0.125,
                color: [255, 214, 102],
                label: if g == 0 {
                    "Protected Germline (Low-ROS)".to_string()
                } else {
                    String::new()
                },
                glow: 1.0,
            });
            scene.edges_3d.push(Edge3D {
                a: g_idx,
                b: soma_base + (g * 2) % outer_n,
                width: 1.3,
                color: [255, 214, 102],
                dashed: true,
            });
        }
    }

    if unit.cumulative_apoptotic_culls > 0 {
        let ap_idx = scene.nodes_3d.len();
        scene.nodes_3d.push(Node3D {
            pos: Vec3::new(1.85, 0.55, 0.0),
            radius: 0.088,
            color: [255, 92, 138],
            label: format!("Mito-Apoptosis Culls ({})", unit.cumulative_apoptotic_culls),
            glow: 1.0,
        });
        scene.edges_3d.push(Edge3D {
            a: soma_base,
            b: ap_idx,
            width: 1.4,
            color: [255, 92, 138],
            dashed: true,
        });
    }
}

fn append_mro_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    match unit.mro_class {
        MroClass::None | MroClass::AerobicMitochondrion | MroClass::AnaerobicMitochondrion => {}
        MroClass::Hydrogenosome => {
            let core_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(-0.62, -0.42, 0.45),
                radius: 0.125,
                color: [6, 214, 160],
                label: "Hydrogenosome (PFO/HydA H2 + ISC, mtDNA=0)".to_string(),
                glow: 1.0,
            });
            for k in 0..3 {
                let ang = (k as f64) * TAU / 3.0 + phase_time * 1.1;
                let h2_idx = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(
                        -0.62 + ang.cos() * 0.28,
                        -0.25 + 0.10 * (phase_time + k as f64).sin(),
                        0.45 + ang.sin() * 0.28,
                    ),
                    radius: 0.068,
                    color: [167, 243, 208],
                    label: if k == 0 { "H2".to_string() } else { String::new() },
                    glow: 0.90,
                });
                scene.edges_3d.push(Edge3D {
                    a: core_idx,
                    b: h2_idx,
                    width: 1.3,
                    color: [6, 214, 160],
                    dashed: true,
                });
            }
        }
        MroClass::Mitosome => {
            let isc_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(-0.52, -0.38, 0.38),
                radius: 0.11,
                color: [255, 209, 102],
                label: "Mitosome ([2Fe-2S] ISC Core Only)".to_string(),
                glow: 1.0,
            });
            if !scene.nodes_3d.is_empty() {
                scene.edges_3d.push(Edge3D {
                    a: isc_idx,
                    b: 0,
                    width: 1.4,
                    color: [255, 209, 102],
                    dashed: true,
                });
            }
        }
        MroClass::Amiteukaryote => {
            let suf_idx = scene.nodes_3d.len();
            scene.nodes_3d.push(Node3D {
                pos: Vec3::new(-0.48, -0.35, 0.35),
                radius: 0.095,
                color: [56, 189, 248],
                label: "Cytosolic SUF (MRO->0)".to_string(),
                glow: 0.95,
            });
            if !scene.nodes_3d.is_empty() {
                scene.edges_3d.push(Edge3D {
                    a: suf_idx,
                    b: 0,
                    width: 1.3,
                    color: [56, 189, 248],
                    dashed: true,
                });
            }
        }
    }
}

fn append_trophic_guild_ornaments(scene: &mut StructureScene, unit: &ModeledUnit, phase_time: f64) {
    match unit.guild {
        TrophicGuild::Producer => {}
        TrophicGuild::Syntroph => {
            let base = scene.nodes_3d.len();
            for k in 0..5 {
                let ang = (k as f64) * TAU / 5.0 + phase_time * 0.6;
                let idx = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(ang.cos() * 1.45, 1.15, ang.sin() * 1.45),
                    radius: 0.085,
                    color: [80, 245, 165],
                    label: if k == 0 {
                        "Syntroph B-Uptake".to_string()
                    } else {
                        String::new()
                    },
                    glow: 0.85,
                });
                if k > 0 {
                    scene.edges_3d.push(Edge3D {
                        a: idx - 1,
                        b: idx,
                        width: 1.4,
                        color: [80, 245, 165],
                        dashed: true,
                    });
                }
            }
            scene.edges_3d.push(Edge3D {
                a: base + 4,
                b: base,
                width: 1.4,
                color: [80, 245, 165],
                dashed: true,
            });
        }
        TrophicGuild::Holobiont => {
            let endo_count = unit.endosymbiont_sequence.len().clamp(3, 10);
            let inner_r = 0.36 + 0.03 * (endo_count as f64);
            let inner_base = scene.nodes_3d.len();
            for k in 0..endo_count {
                let ang = (k as f64) * TAU / (endo_count as f64) - phase_time * 0.75;
                let idx = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(
                        ang.cos() * inner_r,
                        0.20 * (ang * 2.0).sin(),
                        ang.sin() * inner_r,
                    ),
                    radius: 0.092,
                    color: [80, 245, 165],
                    label: if k == 0 {
                        format!("CoRR-Org ({} loci)", unit.endosymbiont_sequence.len())
                    } else {
                        String::new()
                    },
                    glow: 0.95,
                });
                if k > 0 {
                    scene.edges_3d.push(Edge3D {
                        a: idx - 1,
                        b: idx,
                        width: 1.8,
                        color: [80, 245, 165],
                        dashed: false,
                    });
                }
            }
            scene.edges_3d.push(Edge3D {
                a: inner_base + endo_count - 1,
                b: inner_base,
                width: 1.8,
                color: [80, 245, 165],
                dashed: false,
            });

            if !unit.nucleomorph_sequence.is_empty() {
                let nm_vis = unit.nucleomorph_sequence.len().min(6);
                let nm_center_ang = phase_time * 0.48;
                let nm_cx = nm_center_ang.cos() * 0.92;
                let nm_cz = nm_center_ang.sin() * 0.92;
                let nm_base = scene.nodes_3d.len();
                for k in 0..nm_vis {
                    let sub_a = (k as f64) * TAU / (nm_vis.max(1) as f64) + phase_time;
                    let idx = scene.nodes_3d.len();
                    scene.nodes_3d.push(Node3D {
                        pos: Vec3::new(
                            nm_cx + sub_a.cos() * 0.16,
                            0.52 + 0.08 * sub_a.sin(),
                            nm_cz + sub_a.sin() * 0.16,
                        ),
                        radius: 0.078,
                        color: [255, 214, 102],
                        label: if k == 0 {
                            format!("Nucleomorph ({} loci)", unit.nucleomorph_sequence.len())
                        } else {
                            String::new()
                        },
                        glow: 1.0,
                    });
                    if k > 0 {
                        scene.edges_3d.push(Edge3D {
                            a: idx - 1,
                            b: idx,
                            width: 1.5,
                            color: [255, 183, 3],
                            dashed: false,
                        });
                    }
                }
                scene.edges_3d.push(Edge3D {
                    a: nm_base,
                    b: inner_base,
                    width: 1.4,
                    color: [255, 214, 102],
                    dashed: true,
                });
            }

            let total_relocated =
                (unit.egt_transferred_loci + unit.nucleomorph_transferred_loci).min(8);
            for k in 0..total_relocated {
                let ang = (k as f64) * TAU / (total_relocated.max(1) as f64) + phase_time * 0.22;
                let idx = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(ang.cos() * 1.12, -0.55 + 0.12 * (k as f64), ang.sin() * 1.12),
                    radius: 0.088,
                    color: [0, 245, 212],
                    label: if k == 0 {
                        format!(
                            "Nuc-EGT (+{})",
                            unit.egt_transferred_loci + unit.nucleomorph_transferred_loci
                        )
                    } else {
                        String::new()
                    },
                    glow: 0.95,
                });
                scene.edges_3d.push(Edge3D {
                    a: idx,
                    b: inner_base + (k % endo_count),
                    width: 1.4,
                    color: [255, 209, 102],
                    dashed: true,
                });
            }

            let ring_layers = unit.endo_tier.clamp(1, 3) as usize;
            let layer_colors: [Rgb; 3] = [[199, 125, 255], [255, 159, 28], [247, 37, 133]];
            for layer in 0..ring_layers {
                let ring_r = 1.46 + 0.16 * (layer as f64);
                let ring_y = 0.95 + 0.14 * (layer as f64);
                let col = layer_colors[layer];
                let ring_base = scene.nodes_3d.len();
                for k in 0..6 {
                    let ang = (k as f64) * TAU / 6.0 + phase_time * (0.32 - 0.06 * (layer as f64));
                    let idx = scene.nodes_3d.len();
                    scene.nodes_3d.push(Node3D {
                        pos: Vec3::new(ang.cos() * ring_r, ring_y, ang.sin() * ring_r),
                        radius: 0.10,
                        color: col,
                        label: if k == 0 && layer + 1 == ring_layers {
                            format!(
                                "{}° Envelope ({}m)",
                                unit.endo_tier.max(1),
                                unit.organelle_membranes.max(2)
                            )
                        } else {
                            String::new()
                        },
                        glow: 0.90,
                    });
                    if k > 0 {
                        scene.edges_3d.push(Edge3D {
                            a: idx - 1,
                            b: idx,
                            width: 1.6,
                            color: col,
                            dashed: true,
                        });
                    }
                }
                scene.edges_3d.push(Edge3D {
                    a: ring_base + 5,
                    b: ring_base,
                    width: 1.6,
                    color: col,
                    dashed: true,
                });
            }
        }
        TrophicGuild::PredatorParasite => {
            for k in 0..4 {
                let ang = (k as f64) * TAU / 4.0 - phase_time * 0.8;
                let tip = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(ang.cos() * 1.55, 0.95, ang.sin() * 1.55),
                    radius: 0.10,
                    color: [255, 82, 118],
                    label: if k == 0 {
                        if unit.has_apicoplast() {
                            "Apical Invasion Conoid".to_string()
                        } else if unit.has_active_kleptoplast() {
                            "Mixotroph Phago-Hook".to_string()
                        } else {
                            "Parasite/Predator Hook".to_string()
                        }
                    } else {
                        String::new()
                    },
                    glow: 0.95,
                });
                if k < scene.nodes_3d.len().saturating_sub(5) {
                    scene.edges_3d.push(Edge3D {
                        a: k,
                        b: tip,
                        width: 1.6,
                        color: [255, 82, 118],
                        dashed: false,
                    });
                }
            }

            if unit.has_apicoplast() {
                let core_idx = scene.nodes_3d.len();
                scene.nodes_3d.push(Node3D {
                    pos: Vec3::new(0.0, 0.42, 0.0),
                    radius: 0.13,
                    color: [255, 209, 102],
                    label: format!(
                        "Apicoplast 4m (MEP/SUF P:{}/4)",
                        unit.apicoplast_pathway_count()
                    ),
                    glow: 1.0,
                });
                let env_pts = 5usize;
                let env_base = scene.nodes_3d.len();
                for k in 0..env_pts {
                    let ang = (k as f64) * TAU / (env_pts as f64) + phase_time * 0.55;
                    let idx = scene.nodes_3d.len();
                    scene.nodes_3d.push(Node3D {
                        pos: Vec3::new(ang.cos() * 0.48, 0.42, ang.sin() * 0.34),
                        radius: 0.078,
                        color: [199, 125, 255],
                        label: String::new(),
                        glow: 0.88,
                    });
                    scene.edges_3d.push(Edge3D {
                        a: core_idx,
                        b: idx,
                        width: 1.4,
                        color: [255, 209, 102],
                        dashed: true,
                    });
                    if k > 0 {
                        scene.edges_3d.push(Edge3D {
                            a: idx - 1,
                            b: idx,
                            width: 1.5,
                            color: [199, 125, 255],
                            dashed: false,
                        });
                    }
                }
                scene.edges_3d.push(Edge3D {
                    a: env_base + env_pts - 1,
                    b: env_base,
                    width: 1.5,
                    color: [199, 125, 255],
                    dashed: false,
                });
            } else if unit.has_active_kleptoplast() {
                let kp_vis = unit.kleptoplast_sequence.len().clamp(3, 6);
                let integ = unit.kleptoplast_integrity.clamp(0.10, 1.0);
                let kp_base = scene.nodes_3d.len();
                for k in 0..kp_vis {
                    let ang = (k as f64) * TAU / (kp_vis as f64) + phase_time * 0.95;
                    let idx = scene.nodes_3d.len();
                    scene.nodes_3d.push(Node3D {
                        pos: Vec3::new(
                            ang.cos() * 0.56,
                            0.35 + 0.10 * (ang * 2.0).sin(),
                            ang.sin() * 0.56,
                        ),
                        radius: 0.065 + 0.040 * integ,
                        color: [163, 255, 87],
                        label: if k == 0 {
                            format!("Kleptoplast ({:.0}%)", integ * 100.0)
                        } else {
                            String::new()
                        },
                        glow: 0.45 + 0.55 * integ,
                    });
                    if k > 0 {
                        scene.edges_3d.push(Edge3D {
                            a: idx - 1,
                            b: idx,
                            width: 1.5,
                            color: [163, 255, 87],
                            dashed: true,
                        });
                    }
                }
                scene.edges_3d.push(Edge3D {
                    a: kp_base + kp_vis - 1,
                    b: kp_base,
                    width: 1.5,
                    color: [163, 255, 87],
                    dashed: true,
                });

                let hgt_vis = unit.klepto_hgt_loci.min(5);
                for h in 0..hgt_vis {
                    let ang = (h as f64) * TAU / (hgt_vis.max(1) as f64) - phase_time * 0.35;
                    let h_idx = scene.nodes_3d.len();
                    scene.nodes_3d.push(Node3D {
                        pos: Vec3::new(
                            ang.cos() * 1.05,
                            -0.35 + 0.14 * (h as f64),
                            ang.sin() * 1.05,
                        ),
                        radius: 0.082,
                        color: [0, 245, 212],
                        label: if h == 0 {
                            format!("Shopping-Bag HGT ({}/4)", unit.klepto_hgt_loci)
                        } else {
                            String::new()
                        },
                        glow: 0.92,
                    });
                    scene.edges_3d.push(Edge3D {
                        a: h_idx,
                        b: kp_base + (h % kp_vis),
                        width: 1.4,
                        color: [0, 245, 212],
                        dashed: true,
                    });
                }
            }
        }
    }
}

fn symbol_color_molecular(sym: u16) -> Rgb {
    match sym {
        0 => [76, 201, 240],
        1 => [247, 127, 0],
        2 => [67, 170, 139],
        3 => [249, 65, 68],
        5 | 8 => [255, 209, 102],
        6 | 7 => [239, 71, 111],
        _ => [168, 130, 255],
    }
}

fn build_rna_dna_protein_scene(
    unit: &ModeledUnit,
    regime: RegimeId,
    population_size: usize,
    phase_time: f64,
) -> StructureScene {
    let symbols = unit.kind.alphabet_symbols();
    let n = unit.sequence.len().max(4);
    let mut nodes_3d = Vec::with_capacity(n * 2 + 40);
    let mut edges_3d = Vec::with_capacity(n * 3 + 32);
    let mut nodes_2d = Vec::with_capacity(n);
    let mut edges_2d = Vec::with_capacity(n * 2);

    for (i, &sym) in unit.sequence.iter().enumerate() {
        let frac = (i as f64) / ((n - 1).max(1) as f64);
        let mem = unit.lifetime_memory.get(i).copied().unwrap_or(0.0);
        let base_col = symbol_color_molecular(sym);
        let label = symbols.get(sym as usize).copied().unwrap_or("?").to_string();

        let pos = match regime {
            RegimeId::R00Static => {
                let gx = ((i % 4) as f64 - 1.5) * 0.65;
                let gy = (((i / 4) % 3) as f64 - 1.0) * 0.65;
                let gz = ((i / 12) as f64 - 0.5) * 0.65;
                Vec3::new(gx, gy, gz)
            }
            RegimeId::R05Decaying => {
                let angle = frac * TAU * 1.7 + (i as f64) * 0.4;
                let drift =
                    1.15 + (1.0 - unit.integrity) * 1.15 + 0.12 * (phase_time + i as f64).sin();
                Vec3::new(
                    angle.cos() * drift,
                    (frac - 0.5) * 2.3 + 0.18 * (i as f64 * 1.9).cos(),
                    angle.sin() * drift,
                )
            }
            RegimeId::R10AdaptiveDissipation => {
                let theta = frac * TAU + phase_time * 0.8;
                let minor = 0.35 * (frac * TAU * 3.0 + phase_time * 1.6).cos();
                let r = 1.25 + minor;
                Vec3::new(
                    theta.cos() * r,
                    0.42 * (frac * TAU * 3.0 + phase_time * 1.6).sin(),
                    theta.sin() * r,
                )
            }
            RegimeId::R15Autoadaptive | RegimeId::R20AutoadaptiveReplicative => {
                let fold_tight = 1.15 - 0.42 * mem.clamp(0.0, 1.0);
                let theta = frac * PI * 2.4 - PI * 0.2;
                let y = (1.0 - (2.0 * frac - 1.0).abs()) * 1.35 - 0.65;
                Vec3::new(
                    theta.cos() * fold_tight,
                    y,
                    theta.sin() * fold_tight * 0.75,
                )
            }
            RegimeId::R15Replicative | RegimeId::R20InteractiveReplicative => {
                let turns = frac * TAU * 2.1;
                let radius = 0.88;
                Vec3::new(turns.cos() * radius, (frac - 0.5) * 2.5, turns.sin() * radius)
            }
            RegimeId::R20ProtoDarwinian => {
                let turns = frac * TAU * 2.2;
                let jitter = 0.24 * ((phase_time * 2.7 + i as f64 * 1.7).sin());
                let radius = 0.95 + jitter;
                Vec3::new(
                    turns.cos() * radius,
                    (frac - 0.5) * 2.4 + jitter * 0.5,
                    turns.sin() * radius,
                )
            }
            RegimeId::R30FullyDarwinian
            | RegimeId::R30PseudoDarwinian
            | RegimeId::R30SoftInheritance
            | RegimeId::R35Symbiogenetic
            | RegimeId::R30PseudoSoftInheritance => {
                if sym <= 3 {
                    let turns = frac * TAU * 2.45;
                    let r = 0.70;
                    Vec3::new(turns.cos() * r, (frac - 0.5) * 2.35, turns.sin() * r)
                } else {
                    let pocket_angle = frac * TAU * 3.2 + 0.7;
                    let r = if matches!(
                        regime,
                        RegimeId::R30PseudoDarwinian | RegimeId::R30PseudoSoftInheritance
                    ) {
                        1.15
                    } else {
                        1.28
                    };
                    Vec3::new(
                        pocket_angle.cos() * r,
                        0.85 * (frac * PI * 2.0 - 1.0).sin(),
                        pocket_angle.sin() * r,
                    )
                }
            }
        };

        nodes_3d.push(Node3D {
            pos,
            radius: if sym <= 3 { 0.11 } else { 0.135 },
            color: base_col,
            label: label.clone(),
            glow: mem.clamp(0.0, 1.0) * 0.75,
        });

        if i > 0 && (regime != RegimeId::R05Decaying || i % 3 != 0) && regime != RegimeId::R00Static
        {
            edges_3d.push(Edge3D {
                a: i - 1,
                b: i,
                width: 2.2,
                color: [170, 195, 220],
                dashed: false,
            });
        }

        let ang_2d = frac * TAU - PI * 0.5;
        nodes_2d.push(Node2D {
            x: ang_2d.cos() * 0.76,
            y: ang_2d.sin() * 0.76,
            radius: 0.055,
            color: base_col,
            label,
        });
        if i > 0 && regime != RegimeId::R00Static && (regime != RegimeId::R05Decaying || i % 3 != 0)
        {
            edges_2d.push(Edge2D {
                a: i - 1,
                b: i,
                width: 1.8,
                color: [120, 145, 175],
                directed: false,
            });
        }
    }

    if matches!(
        regime,
        RegimeId::R15Replicative
            | RegimeId::R20InteractiveReplicative
            | RegimeId::R30FullyDarwinian
            | RegimeId::R30PseudoDarwinian
            | RegimeId::R30SoftInheritance
            | RegimeId::R35Symbiogenetic
            | RegimeId::R30PseudoSoftInheritance
    ) {
        let base_count = unit.sequence.len().min(14);
        let offset = nodes_3d.len();
        for k in 0..base_count {
            let p = nodes_3d[k].pos;
            let comp_pos = Vec3::new(-p.x * 0.88, p.y, -p.z * 0.88);
            let comp_idx = nodes_3d.len();
            nodes_3d.push(Node3D {
                pos: comp_pos,
                radius: 0.09,
                color: [95, 225, 195],
                label: String::new(),
                glow: 0.25,
            });
            edges_3d.push(Edge3D {
                a: k,
                b: comp_idx,
                width: 1.2,
                color: [110, 240, 210],
                dashed: true,
            });
            if k > 0 {
                edges_3d.push(Edge3D {
                    a: offset + k - 1,
                    b: comp_idx,
                    width: 1.6,
                    color: [85, 180, 165],
                    dashed: false,
                });
            }
        }
    }

    if matches!(
        regime,
        RegimeId::R30SoftInheritance
            | RegimeId::R35Symbiogenetic
            | RegimeId::R30PseudoSoftInheritance
    ) {
        for m_idx in (0..unit.sequence.len().min(16)).step_by(4) {
            let anchor = nodes_3d[m_idx].pos;
            let tag_idx = nodes_3d.len();
            nodes_3d.push(Node3D {
                pos: anchor.scale(1.26),
                radius: 0.075,
                color: [80, 250, 190],
                label: "Me".to_string(),
                glow: 0.9,
            });
            edges_3d.push(Edge3D {
                a: m_idx,
                b: tag_idx,
                width: 1.3,
                color: [80, 250, 190],
                dashed: false,
            });
        }
    }

    if matches!(
        regime,
        RegimeId::R15Autoadaptive
            | RegimeId::R20AutoadaptiveReplicative
            | RegimeId::R30FullyDarwinian
            | RegimeId::R30PseudoDarwinian
            | RegimeId::R30SoftInheritance
            | RegimeId::R35Symbiogenetic
            | RegimeId::R30PseudoSoftInheritance
    ) {
        let mg_idx = nodes_3d.len();
        let engineered = matches!(
            regime,
            RegimeId::R30PseudoDarwinian | RegimeId::R30PseudoSoftInheritance
        );
        let cofactor_color = if engineered {
            [255, 95, 210]
        } else {
            [140, 255, 120]
        };
        let cofactor_label = if engineered { "Ligand*" } else { "Mg2+" };
        nodes_3d.push(Node3D {
            pos: Vec3::new(0.0, 0.0, 0.0),
            radius: 0.18,
            color: cofactor_color,
            label: cofactor_label.to_string(),
            glow: 0.95,
        });
        for i in (0..unit.sequence.len()).step_by(5) {
            edges_3d.push(Edge3D {
                a: i,
                b: mg_idx,
                width: 1.1,
                color: cofactor_color,
                dashed: true,
            });
        }
    }

    for &(u, v) in &unit.wiring {
        if u < unit.sequence.len() && v < unit.sequence.len() {
            edges_3d.push(Edge3D {
                a: u,
                b: v,
                width: 1.4,
                color: [255, 214, 102],
                dashed: true,
            });
            edges_2d.push(Edge2D {
                a: u,
                b: v,
                width: 1.5,
                color: [255, 205, 86],
                directed: false,
            });
        }
    }

    append_population_satellites(&mut nodes_3d, population_size, [76, 201, 240], phase_time);

    StructureScene {
        kind: unit.kind,
        regime,
        subtitle: format!(
            "3D RNP Fold | {} ({}°) | MRO:{} | Apico:{} | S:{:.1}b",
            unit.guild.short_code(),
            unit.endo_tier,
            unit.mro_class.short_code(),
            unit.apicoplast_pathway_count(),
            unit.structure_s
        ),
        nodes_3d,
        edges_3d,
        nodes_2d,
        edges_2d,
    }
}

fn build_hominid_scene(
    unit: &ModeledUnit,
    regime: RegimeId,
    population_size: usize,
    phase_time: f64,
) -> StructureScene {
    let symbols = unit.kind.alphabet_symbols();
    let mut nodes_3d = Vec::new();
    let mut edges_3d = Vec::new();
    let mut nodes_2d = Vec::new();
    let mut edges_2d = Vec::new();

    let rings = 5usize;
    let radial_pts = 6usize;
    let symmetry = match regime {
        RegimeId::R00Static => 0.15,
        RegimeId::R05Decaying => 0.05,
        RegimeId::R10AdaptiveDissipation => 0.25,
        RegimeId::R15Autoadaptive => 0.52,
        RegimeId::R15Replicative | RegimeId::R20InteractiveReplicative => 0.48,
        RegimeId::R20AutoadaptiveReplicative => 0.62,
        RegimeId::R20ProtoDarwinian => 0.35,
        RegimeId::R30FullyDarwinian => 0.88,
        RegimeId::R30PseudoDarwinian => 0.94,
        RegimeId::R30SoftInheritance
        | RegimeId::R35Symbiogenetic
        | RegimeId::R30PseudoSoftInheritance => 1.0,
    };

    for r_idx in 0..rings {
        let v = (r_idx as f64) / ((rings - 1) as f64);
        let y = (v - 0.45) * 1.95;
        let cobble_profile = (1.0 - (2.0 * v - 1.0).powi(2)).max(0.08).sqrt() * 0.78;
        let biface_width = (1.0 - v.powf(1.25)) * 0.85 + 0.07;
        let rx = cobble_profile * (1.0 - symmetry) + biface_width * symmetry;
        let rz = cobble_profile * (1.0 - 0.70 * symmetry);

        for k in 0..radial_pts {
            let theta = (k as f64) * TAU / (radial_pts as f64);
            let scar_jitter = if regime == RegimeId::R20ProtoDarwinian {
                0.16 * ((theta * 3.0 + r_idx as f64 * 2.1 + phase_time * 2.0).sin())
            } else if regime == RegimeId::R05Decaying {
                0.28 * ((k + r_idx) as f64 * 1.7).cos()
            } else {
                0.0
            };

            let pos = Vec3::new(
                theta.cos() * (rx + scar_jitter),
                y,
                theta.sin() * (rz + scar_jitter * 0.6),
            );
            let stone_col = if regime >= RegimeId::R30SoftInheritance {
                [115, 235, 195]
            } else if regime >= RegimeId::R30FullyDarwinian {
                [224, 159, 62]
            } else if regime == RegimeId::R20ProtoDarwinian {
                [198, 112, 92]
            } else {
                [148, 156, 168]
            };
            nodes_3d.push(Node3D {
                pos,
                radius: 0.065,
                color: stone_col,
                label: if r_idx == rings - 1 && k == 0 && symmetry > 0.7 {
                    "Biface Tip".to_string()
                } else {
                    String::new()
                },
                glow: symmetry * 0.40,
            });

            let curr = r_idx * radial_pts + k;
            let next_in_ring = r_idx * radial_pts + ((k + 1) % radial_pts);
            if regime != RegimeId::R05Decaying || (k + r_idx) % 2 == 0 {
                edges_3d.push(Edge3D {
                    a: curr,
                    b: next_in_ring,
                    width: 1.5,
                    color: stone_col,
                    dashed: false,
                });
                if r_idx > 0 {
                    edges_3d.push(Edge3D {
                        a: (r_idx - 1) * radial_pts + k,
                        b: curr,
                        width: 1.4,
                        color: [185, 145, 95],
                        dashed: false,
                    });
                }
            }
        }
    }

    if regime >= RegimeId::R30PseudoDarwinian
        || (regime == RegimeId::R30FullyDarwinian && unit.sequence.contains(&5))
    {
        let top_haft = nodes_3d.len();
        nodes_3d.push(Node3D {
            pos: Vec3::new(0.0, -0.85, 0.0),
            radius: 0.10,
            color: [168, 112, 64],
            label: "Haft Binding".to_string(),
            glow: 0.5,
        });
        let base_haft = nodes_3d.len();
        nodes_3d.push(Node3D {
            pos: Vec3::new(0.0, -1.75, 0.0),
            radius: 0.09,
            color: [142, 92, 50],
            label: "Shaft".to_string(),
            glow: 0.2,
        });
        edges_3d.push(Edge3D {
            a: top_haft,
            b: base_haft,
            width: 3.2,
            color: [168, 112, 64],
            dashed: false,
        });
    }

    let tribe_count = match regime {
        RegimeId::R00Static | RegimeId::R05Decaying => 0,
        RegimeId::R10AdaptiveDissipation | RegimeId::R15Autoadaptive => 1,
        _ => population_size.clamp(3, 9),
    };
    let hearth_start = nodes_3d.len();
    for t in 0..tribe_count {
        let ang = (t as f64) * TAU / (tribe_count.max(1) as f64) + phase_time * 0.18;
        let r = 1.72;
        nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * r, -0.45 + 0.12 * (t as f64).sin(), ang.sin() * r),
            radius: 0.13,
            color: [255, 123, 84],
            label: format!("Band {}", t + 1),
            glow: 0.75,
        });
        if t > 0 && regime >= RegimeId::R20InteractiveReplicative {
            edges_3d.push(Edge3D {
                a: hearth_start + t - 1,
                b: hearth_start + t,
                width: 1.4,
                color: [255, 180, 105],
                dashed: true,
            });
        }
        if regime >= RegimeId::R30SoftInheritance && t % 2 == 0 {
            edges_3d.push(Edge3D {
                a: hearth_start + t,
                b: (rings - 1) * radial_pts,
                width: 1.3,
                color: [115, 245, 200],
                dashed: true,
            });
        }
    }
    if tribe_count >= 3 && regime >= RegimeId::R30FullyDarwinian {
        edges_3d.push(Edge3D {
            a: hearth_start + tribe_count - 1,
            b: hearth_start,
            width: 1.4,
            color: [255, 180, 105],
            dashed: true,
        });
    }

    for (i, &sym) in unit.sequence.iter().enumerate() {
        let col_idx = i % 6;
        let row_idx = i / 6;
        let x = -0.75 + (col_idx as f64) * 0.30;
        let y = 0.55 - (row_idx as f64) * 0.36;
        let col = if sym <= 5 {
            [224, 159, 62]
        } else {
            [92, 200, 255]
        };
        nodes_2d.push(Node2D {
            x: x.clamp(-0.9, 0.9),
            y: y.clamp(-0.9, 0.9),
            radius: 0.058,
            color: col,
            label: symbols.get(sym as usize).copied().unwrap_or("?").to_string(),
        });
        if i > 0 && regime != RegimeId::R00Static {
            edges_2d.push(Edge2D {
                a: i - 1,
                b: i,
                width: 1.8,
                color: [160, 175, 195],
                directed: true,
            });
        }
    }
    for &(u, v) in &unit.wiring {
        if u < nodes_2d.len() && v < nodes_2d.len() {
            edges_2d.push(Edge2D {
                a: u,
                b: v,
                width: 1.3,
                color: [255, 140, 95],
                directed: true,
            });
        }
    }

    StructureScene {
        kind: unit.kind,
        regime,
        subtitle: format!(
            "3D Lithic Biface | {} ({}°) | MRO:{} | Apico:{} | S:{:.1}b",
            unit.guild.short_code(),
            unit.endo_tier,
            unit.mro_class.short_code(),
            unit.apicoplast_pathway_count(),
            unit.structure_s
        ),
        nodes_3d,
        edges_3d,
        nodes_2d,
        edges_2d,
    }
}

fn build_ai_scene(
    unit: &ModeledUnit,
    regime: RegimeId,
    population_size: usize,
    phase_time: f64,
) -> StructureScene {
    let symbols = unit.kind.alphabet_symbols();
    let n = unit.sequence.len().max(2);
    let mut nodes_3d = Vec::with_capacity(n + 40);
    let mut edges_3d = Vec::with_capacity(n * 3 + 28);
    let mut nodes_2d = Vec::with_capacity(n);
    let mut edges_2d = Vec::with_capacity(n * 2);

    let layers = 5usize;
    for (i, &sym) in unit.sequence.iter().enumerate() {
        let layer = if i == 0 {
            0
        } else if i + 1 == n {
            layers - 1
        } else {
            1 + ((i - 1) * (layers - 2)) / (n.saturating_sub(2).max(1))
        };

        let x_layer = -1.55 + (layer as f64) * (3.10 / ((layers - 1) as f64));
        let ring_angle = (i as f64) * 1.85
            + if regime == RegimeId::R10AdaptiveDissipation {
                phase_time * 1.1
            } else {
                0.0
            };
        let radial_spread = if i == 0 || i + 1 == n {
            0.0
        } else if regime == RegimeId::R00Static {
            0.55
        } else {
            0.85
        };

        let jitter = if regime == RegimeId::R05Decaying {
            0.45 * (1.0 - unit.integrity) * ((i as f64 * 2.3).sin())
        } else if regime == RegimeId::R20ProtoDarwinian {
            0.22 * ((phase_time * 3.1 + i as f64).cos())
        } else {
            0.0
        };

        let pos = Vec3::new(
            x_layer + jitter,
            ring_angle.cos() * radial_spread,
            ring_angle.sin() * radial_spread,
        );

        let mem = unit.lifetime_memory.get(i).copied().unwrap_or(0.0);
        let op_color = match sym {
            0 | 9 => [255, 214, 102],
            4 => [0, 245, 212],
            5 | 6 => [241, 91, 181],
            7 | 8 => [155, 93, 229],
            _ => [72, 149, 239],
        };
        let dimmed = if regime == RegimeId::R00Static {
            [op_color[0] / 2, op_color[1] / 2, op_color[2] / 2]
        } else {
            op_color
        };

        let label = symbols.get(sym as usize).copied().unwrap_or("op").to_string();
        nodes_3d.push(Node3D {
            pos,
            radius: if sym == 0 || sym == 9 { 0.16 } else { 0.115 },
            color: dimmed,
            label: label.clone(),
            glow: if regime == RegimeId::R00Static {
                0.0
            } else {
                (0.25 + mem * 0.65).clamp(0.1, 1.0)
            },
        });

        let y_2d = if i == 0 || i + 1 == n {
            0.0
        } else {
            (((i % 3) as f64) - 1.0) * 0.52
        };
        let x_2d = -0.84 + 1.68 * ((i as f64) / ((n - 1).max(1) as f64));
        nodes_2d.push(Node2D {
            x: x_2d,
            y: y_2d,
            radius: 0.055,
            color: dimmed,
            label,
        });

        if i > 0 && (regime != RegimeId::R05Decaying || i % 2 == 1) {
            let flow_col = if regime == RegimeId::R00Static {
                [80, 90, 105]
            } else {
                [110, 200, 255]
            };
            edges_3d.push(Edge3D {
                a: i - 1,
                b: i,
                width: 2.0,
                color: flow_col,
                dashed: false,
            });
            edges_2d.push(Edge2D {
                a: i - 1,
                b: i,
                width: 1.8,
                color: flow_col,
                directed: true,
            });
        }
    }

    if regime >= RegimeId::R30SoftInheritance {
        let teacher_idx = nodes_3d.len();
        nodes_3d.push(Node3D {
            pos: Vec3::new(0.0, 1.45, 0.0),
            radius: 0.17,
            color: [115, 250, 180],
            label: "Teacher Distill".to_string(),
            glow: 0.95,
        });
        for k in (1..n.saturating_sub(1)).step_by(4) {
            edges_3d.push(Edge3D {
                a: teacher_idx,
                b: k,
                width: 1.3,
                color: [115, 250, 180],
                dashed: true,
            });
        }
    }

    for &(u, v) in &unit.wiring {
        if u < nodes_3d.len() && u < n && v < nodes_3d.len() && v < n {
            edges_3d.push(Edge3D {
                a: u,
                b: v,
                width: 1.5,
                color: [0, 245, 212],
                dashed: true,
            });
            if u < nodes_2d.len() && v < nodes_2d.len() {
                edges_2d.push(Edge2D {
                    a: u,
                    b: v,
                    width: 1.4,
                    color: [0, 245, 212],
                    directed: true,
                });
            }
        }
    }

    append_population_satellites(&mut nodes_3d, population_size, [0, 245, 212], phase_time);

    StructureScene {
        kind: unit.kind,
        regime,
        subtitle: format!(
            "3D Neural DAG | {} ({}°) | MRO:{} | Apico:{} | S:{:.1}b",
            unit.guild.short_code(),
            unit.endo_tier,
            unit.mro_class.short_code(),
            unit.apicoplast_pathway_count(),
            unit.structure_s
        ),
        nodes_3d,
        edges_3d,
        nodes_2d,
        edges_2d,
    }
}

fn append_population_satellites(
    nodes_3d: &mut Vec<Node3D>,
    population_size: usize,
    color: Rgb,
    phase_time: f64,
) {
    if population_size <= 1 {
        return;
    }
    let count = population_size.min(12);
    for s in 0..count {
        let ang = (s as f64) * TAU / (count as f64) - phase_time * 0.25;
        let r = 1.92;
        nodes_3d.push(Node3D {
            pos: Vec3::new(ang.cos() * r, -1.25, ang.sin() * r),
            radius: 0.048,
            color,
            label: String::new(),
            glow: 0.3,
        });
    }
}

pub fn render_ascii_3d(scene: &StructureScene, camera: &Camera3D, cols: usize, rows: usize) -> String {
    let w = cols.max(32);
    let h = rows.max(12);
    let mut chars = vec![vec![' '; w]; h];
    let mut zbuf = vec![vec![f64::INFINITY; w]; h];

    let proj: Vec<(f64, f64, f64, f64)> = scene
        .nodes_3d
        .iter()
        .map(|n| camera.project(n.pos, w as f64, (h * 2) as f64))
        .collect();

    for edge in &scene.edges_3d {
        if edge.a >= proj.len() || edge.b >= proj.len() {
            continue;
        }
        let (x0, y0_raw, z0, _) = proj[edge.a];
        let (x1, y1_raw, z1, _) = proj[edge.b];
        let y0 = y0_raw * 0.5;
        let y1 = y1_raw * 0.5;
        let steps = ((x1 - x0).abs().max((y1 - y0).abs()) as usize).clamp(2, 120);
        for s in 0..=steps {
            if edge.dashed && s % 2 == 1 {
                continue;
            }
            let t = (s as f64) / (steps as f64);
            let xi = (x0 + (x1 - x0) * t).round() as isize;
            let yi = (y0 + (y1 - y0) * t).round() as isize;
            let z = z0 + (z1 - z0) * t;
            if xi >= 0 && (xi as usize) < w && yi >= 0 && (yi as usize) < h {
                let (ux, uy) = (xi as usize, yi as usize);
                if z < zbuf[uy][ux] {
                    zbuf[uy][ux] = z;
                    chars[uy][ux] = if edge.dashed { ':' } else { '.' };
                }
            }
        }
    }

    let mut order: Vec<usize> = (0..scene.nodes_3d.len()).collect();
    order.sort_by(|&a, &b| proj[b].2.total_cmp(&proj[a].2));
    for idx in order {
        let node = &scene.nodes_3d[idx];
        let (x, y_raw, z, _) = proj[idx];
        let xi = x.round() as isize;
        let yi = (y_raw * 0.5).round() as isize;
        if xi >= 0 && (xi as usize) < w && yi >= 0 && (yi as usize) < h {
            let (ux, uy) = (xi as usize, yi as usize);
            if z <= zbuf[uy][ux] + 0.25 {
                zbuf[uy][ux] = z;
                chars[uy][ux] = if node.radius >= 0.15 {
                    '@'
                } else if node.glow > 0.45 {
                    '*'
                } else if node.radius >= 0.10 {
                    'O'
                } else {
                    'o'
                };
            }
        }
    }

    let mut out = String::with_capacity((w + 4) * (h + 3));
    let border: String = "-".repeat(w);
    let _ = writeln!(out, "+{}+", border);
    for row in chars {
        let line: String = row.into_iter().collect();
        let _ = writeln!(out, "|{}|", line);
    }
    let _ = writeln!(out, "+{}+", border);
    out
}

pub fn render_stage_atlas_svg(
    title: &str,
    panels: &[(String, String, StructureScene)],
) -> String {
    let cols = 6usize;
    let rows = (panels.len() + cols - 1) / cols;
    let cell_w = 275.0;
    let cell_h = 250.0;
    let pad = 24.0;
    let header_h = 64.0;
    let total_w = pad * 2.0 + (cols as f64) * cell_w;
    let total_h = header_h + pad * 2.0 + (rows.max(1) as f64) * cell_h;

    let mut svg = String::with_capacity(72_000);
    let _ = writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w:.0} {total_h:.0}" width="{total_w:.0}" height="{total_h:.0}">"#
    );
    let _ = writeln!(
        svg,
        r##"<rect width="100%" height="100%" fill="#0b1020"/>"##
    );
    let _ = writeln!(
        svg,
        r##"<text x="{pad:.0}" y="38" fill="#f0f4ff" font-family="monospace, sans-serif" font-size="19" font-weight="bold">{}</text>"##,
        escape_xml(title)
    );

    let cam = Camera3D {
        yaw: 0.62,
        pitch: 0.34,
        zoom: 1.05,
        distance: 4.9,
    };

    for (idx, (panel_title, metric_line, scene)) in panels.iter().enumerate() {
        let col = idx % cols;
        let row = idx / cols;
        let ox = pad + (col as f64) * cell_w;
        let oy = header_h + (row as f64) * cell_h;
        let box_w = cell_w - 12.0;
        let box_h = cell_h - 12.0;

        let _ = writeln!(
            svg,
            r##"<rect x="{ox:.1}" y="{oy:.1}" width="{box_w:.1}" height="{box_h:.1}" rx="8" fill="#131b31" stroke="#283659" stroke-width="1.2"/>"##
        );
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" fill="#e2ebff" font-family="monospace, sans-serif" font-size="11" font-weight="bold">{}</text>"##,
            ox + 10.0,
            oy + 19.0,
            escape_xml(panel_title)
        );
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" fill="#8fa8d6" font-family="monospace, sans-serif" font-size="9">{}</text>"##,
            ox + 10.0,
            oy + 33.0,
            escape_xml(metric_line)
        );

        let vp_w = box_w - 74.0;
        let vp_h = box_h - 46.0;
        let vp_x = ox + 6.0;
        let vp_y = oy + 40.0;

        let proj: Vec<(f64, f64, f64, f64)> = scene
            .nodes_3d
            .iter()
            .map(|n| {
                let (sx, sy, d, s) = cam.project(n.pos, vp_w, vp_h);
                (vp_x + sx, vp_y + sy, d, s)
            })
            .collect();

        for e in &scene.edges_3d {
            if e.a < proj.len() && e.b < proj.len() {
                let (x1, y1, _, _) = proj[e.a];
                let (x2, y2, _, _) = proj[e.b];
                let [r, g, b] = e.color;
                let dash = if e.dashed {
                    r#" stroke-dasharray="3,3""#
                } else {
                    ""
                };
                let _ = writeln!(
                    svg,
                    r#"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="rgb({r},{g},{b})" stroke-opacity="0.72" stroke-width="{:.1}"{dash}/>"#,
                    e.width * 0.85
                );
            }
        }

        let mut order: Vec<usize> = (0..scene.nodes_3d.len()).collect();
        order.sort_by(|&a, &b| proj[b].2.total_cmp(&proj[a].2));
        for ni in order {
            let node = &scene.nodes_3d[ni];
            let (cx, cy, _, persp) = proj[ni];
            if cx < vp_x || cx > vp_x + vp_w || cy < vp_y || cy > vp_y + vp_h {
                continue;
            }
            let rad = (node.radius * persp).clamp(1.8, 10.5);
            let [r, g, b] = node.color;
            if node.glow > 0.35 {
                let _ = writeln!(
                    svg,
                    r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{:.1}" fill="rgb({r},{g},{b})" fill-opacity="0.22"/>"#,
                    rad * 1.85
                );
            }
            let _ = writeln!(
                svg,
                r##"<circle cx="{cx:.1}" cy="{cy:.1}" r="{rad:.1}" fill="rgb({r},{g},{b})" stroke="#0a0e1a" stroke-width="0.7"/>"##
            );
        }

        let in_size = 60.0;
        let in_x = ox + box_w - in_size - 8.0;
        let in_y = oy + box_h - in_size - 8.0;
        let _ = writeln!(
            svg,
            r##"<rect x="{in_x:.1}" y="{in_y:.1}" width="{in_size:.1}" height="{in_size:.1}" rx="5" fill="#0d1426" stroke="#2c3e66" stroke-width="0.9"/>"##
        );
        for e2 in &scene.edges_2d {
            if e2.a < scene.nodes_2d.len() && e2.b < scene.nodes_2d.len() {
                let na = &scene.nodes_2d[e2.a];
                let nb = &scene.nodes_2d[e2.b];
                let x1 = in_x + in_size * (0.5 + na.x * 0.42);
                let y1 = in_y + in_size * (0.5 - na.y * 0.42);
                let x2 = in_x + in_size * (0.5 + nb.x * 0.42);
                let y2 = in_y + in_size * (0.5 - nb.y * 0.42);
                let [r, g, b] = e2.color;
                let _ = writeln!(
                    svg,
                    r#"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="rgb({r},{g},{b})" stroke-opacity="0.65" stroke-width="0.9"/>"#
                );
            }
        }
        for n2 in &scene.nodes_2d {
            let cx = in_x + in_size * (0.5 + n2.x * 0.42);
            let cy = in_y + in_size * (0.5 - n2.y * 0.42);
            let [r, g, b] = n2.color;
            let _ = writeln!(
                svg,
                r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="1.9" fill="rgb({r},{g},{b})"/>"#
            );
        }
    }

    let _ = writeln!(svg, "</svg>");
    svg
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
