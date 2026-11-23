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

//! ITME v7.0 Multi-Trophic & Holobiont Spatial Ecology Engine:
//! - 2D Reaction-Diffusion Primary Exergy Grid `E(x, y, t)`
//! - Diffusible Secondary Metabolite / Cross-Feeding Public-Good Field `B(x, y, t)`
//! - Persistent Ecological-Inheritance / Niche-Construction Field `N(x, y, t)`
//! - Trophic-Aware Chemotactic Foraging (`Producer`, `Syntroph`, `Holobiont`, `PredatorParasite`)
//! - Assortative Genotypic Compatibility, Wright's `F_ST` & Within-Species Adami Conservation

use crate::model::TrophicGuild;
use std::fmt::Write as _;

/// Default lattice width & height (`16 x 16 = 256` spatial patches).
pub const DEFAULT_GRID_DIM: usize = 16;

/// Ecological biome zone on the spatial grid, defining divergent selective targets.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BiomeZone {
    /// North-West hotspot (Acidic Hydrothermal Vent / Highland Chert Quarry / Dense Tensor Cluster)
    Alpha = 0,
    /// South-East hotspot (Alkaline Cold Seep / Riverine Woodland / Edge Sensor Mesh)
    Beta = 1,
}

/// 2D reaction-diffusion lattice coupling primary exergy `E(x, y)`, secondary cross-feeding
/// byproducts `B(x, y)`, and constructed niche capital `N(x, y)`.
#[derive(Clone, Debug)]
pub struct SpatialGrid {
    pub width: usize,
    pub height: usize,
    pub exergy: Vec<f64>,
    pub byproduct: Vec<f64>,
    pub niche: Vec<f64>,
    pub source_profile: Vec<f64>,
    pub biomes: Vec<BiomeZone>,
    scratch: Vec<f64>,
}

impl SpatialGrid {
    pub fn new(width: usize, height: usize, base_flux: f64) -> Self {
        let w = width.max(6);
        let h = height.max(6);
        let len = w * h;
        let mut source_profile = vec![0.0; len];
        let mut biomes = vec![BiomeZone::Alpha; len];

        let ax = (w as f64) * 0.26;
        let ay = (h as f64) * 0.28;
        let bx = (w as f64) * 0.74;
        let by = (h as f64) * 0.72;
        let sigma2 = ((w.min(h) as f64) * 0.24).powi(2).max(1.5);

        for y in 0..h {
            for x in 0..w {
                let idx = y * w + x;
                let d_alpha = ((x as f64 - ax).powi(2) + (y as f64 - ay).powi(2)) / sigma2;
                let d_beta = ((x as f64 - bx).powi(2) + (y as f64 - by).powi(2)) / sigma2;

                let g_alpha = (-0.5 * d_alpha).exp();
                let g_beta = (-0.5 * d_beta).exp();
                source_profile[idx] = (0.32 + 1.18 * (g_alpha + g_beta)).clamp(0.25, 1.55);
                biomes[idx] = if x + y < (w + h) / 2 {
                    BiomeZone::Alpha
                } else {
                    BiomeZone::Beta
                };
            }
        }

        let safe_flux = if base_flux.is_finite() { base_flux } else { 0.0 };
        let init_e = if safe_flux.abs() <= 0.03 {
            0.0
        } else {
            safe_flux
        };
        let exergy = source_profile
            .iter()
            .map(|&s| if init_e <= 0.0 { init_e } else { init_e * s })
            .collect();

        Self {
            width: w,
            height: h,
            exergy,
            byproduct: vec![0.0; len],
            niche: vec![0.0; len],
            source_profile,
            biomes,
            scratch: vec![0.0; len],
        }
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        (y % self.height) * self.width + (x % self.width)
    }

    #[inline]
    fn toroidal_stencil(&self, x: usize, y: usize) -> (usize, usize, usize, usize, usize) {
        let w = self.width;
        let h = self.height;
        let c = self.idx(x, y);
        let left = self.idx((x + w - 1) % w, y);
        let right = self.idx((x + 1) % w, y);
        let up = self.idx(x, (y + h - 1) % h);
        let down = self.idx(x, (y + 1) % h);
        (c, left, right, up, down)
    }

    #[inline]
    pub fn local_exergy(&self, x: usize, y: usize) -> f64 {
        self.exergy[self.idx(x, y)]
    }

    #[inline]
    pub fn local_byproduct(&self, x: usize, y: usize) -> f64 {
        self.byproduct[self.idx(x, y)]
    }

    #[inline]
    pub fn local_niche(&self, x: usize, y: usize) -> f64 {
        self.niche[self.idx(x, y)]
    }

    #[inline]
    pub fn biome_at(&self, x: usize, y: usize) -> BiomeZone {
        self.biomes[self.idx(x, y)]
    }

    pub fn apply_trophic_flux(
        &mut self,
        x: usize,
        y: usize,
        consumed_e: f64,
        produced_b: f64,
        consumed_b: f64,
        deposited_n: f64,
    ) {
        let i = self.idx(x, y);
        if self.exergy[i] > 0.0 && consumed_e > 0.0 {
            self.exergy[i] = (self.exergy[i] - consumed_e * 0.22).max(0.04);
        }
        self.byproduct[i] =
            (self.byproduct[i] + produced_b - consumed_b * 0.65).clamp(0.0, 2.2);
        self.niche[i] = (self.niche[i] + deposited_n).clamp(0.0, 2.5);
    }

    pub fn step_reaction_diffusion(
        &mut self,
        bath_exergy_flux: f64,
        exergy_diffusion: f64,
        niche_decay: f64,
        niche_diffusion: f64,
    ) {
        if !bath_exergy_flux.is_finite()
            || !exergy_diffusion.is_finite()
            || !niche_decay.is_finite()
            || !niche_diffusion.is_finite()
        {
            return;
        }
        if bath_exergy_flux.abs() <= 0.03 {
            self.exergy.fill(0.0);
            for b in &mut self.byproduct {
                *b *= 0.90;
            }
            for n in &mut self.niche {
                *n *= 0.95;
            }
            return;
        }

        if bath_exergy_flux < -0.03 {
            self.exergy.fill(bath_exergy_flux);
            for b in &mut self.byproduct {
                *b *= 0.78;
            }
            for n in &mut self.niche {
                *n *= (1.0 - niche_decay * 2.0).clamp(0.5, 0.96);
            }
            return;
        }

        let d_e = exergy_diffusion.clamp(0.0, 0.24);
        let replenish = 0.26;
        let w = self.width;
        let h = self.height;

        for y in 0..h {
            for x in 0..w {
                let (c, left, right, up, down) = self.toroidal_stencil(x, y);

                let lap_e = self.exergy[left] + self.exergy[right] + self.exergy[up]
                    + self.exergy[down]
                    - 4.0 * self.exergy[c];
                let target_e = bath_exergy_flux * self.source_profile[c];
                self.scratch[c] = (self.exergy[c]
                    + replenish * (target_e - self.exergy[c])
                    + d_e * lap_e)
                    .clamp(0.03, 3.0);
            }
        }
        self.exergy.copy_from_slice(&self.scratch);

        let d_b = (d_e * 0.85).clamp(0.015, 0.20);
        let decay_b = 0.085;
        for y in 0..h {
            for x in 0..w {
                let (c, left, right, up, down) = self.toroidal_stencil(x, y);

                let lap_b = self.byproduct[left] + self.byproduct[right] + self.byproduct[up]
                    + self.byproduct[down]
                    - 4.0 * self.byproduct[c];
                self.scratch[c] =
                    (self.byproduct[c] * (1.0 - decay_b) + d_b * lap_b).clamp(0.0, 2.2);
            }
        }
        self.byproduct.copy_from_slice(&self.scratch);

        let d_n = niche_diffusion.clamp(0.0, 0.15);
        let decay = niche_decay.clamp(0.005, 0.35);
        for y in 0..h {
            for x in 0..w {
                let (c, left, right, up, down) = self.toroidal_stencil(x, y);

                let lap_n = self.niche[left] + self.niche[right] + self.niche[up]
                    + self.niche[down]
                    - 4.0 * self.niche[c];
                self.scratch[c] =
                    (self.niche[c] * (1.0 - decay) + d_n * lap_n).clamp(0.0, 2.5);
            }
        }
        self.niche.copy_from_slice(&self.scratch);
    }

    pub fn mean_exergy(&self) -> f64 {
        if self.exergy.is_empty() {
            0.0
        } else {
            self.exergy.iter().sum::<f64>() / (self.exergy.len() as f64)
        }
    }

    pub fn mean_byproduct(&self) -> f64 {
        if self.byproduct.is_empty() {
            0.0
        } else {
            self.byproduct.iter().sum::<f64>() / (self.byproduct.len() as f64)
        }
    }

    pub fn mean_niche(&self) -> f64 {
        if self.niche.is_empty() {
            0.0
        } else {
            self.niche.iter().sum::<f64>() / (self.niche.len() as f64)
        }
    }

    /// Trophic-guild-specific chemotactic/cognitive foraging across the `3 x 3` Moore neighborhood:
    /// - `Producer`: climbs `E(x, y) + 0.65 N(x, y) - 0.28 predator_density`
    /// - `Syntroph`: climbs `1.45 B(x, y) + 0.55 N(x, y) + 0.25 host_density`
    /// - `Holobiont`: climbs `1.15 E(x, y) + 0.85 B(x, y) + 0.80 N(x, y) - 0.22 predator_density`
    /// - `PredatorParasite`: climbs `0.85 host_density + 0.28 E(x, y) - 0.22 predator_density`
    pub fn best_guild_neighbor(
        &self,
        x: usize,
        y: usize,
        guild: TrophicGuild,
        host_counts: &[u16],
        predator_counts: &[u16],
    ) -> (usize, usize) {
        let mut best = (x, y);
        let eval_cell = |cx: usize, cy: usize| -> f64 {
            let i = self.idx(cx, cy);
            let h_c = host_counts.get(i).copied().unwrap_or(0) as f64;
            let p_c = predator_counts.get(i).copied().unwrap_or(0) as f64;
            match guild {
                TrophicGuild::Producer => self.exergy[i] + 0.65 * self.niche[i] - 0.28 * p_c,
                TrophicGuild::Syntroph => {
                    1.45 * self.byproduct[i] + 0.55 * self.niche[i] + 0.25 * h_c - 0.20 * p_c
                }
                TrophicGuild::Holobiont => {
                    1.15 * self.exergy[i]
                        + 0.85 * self.byproduct[i]
                        + 0.80 * self.niche[i]
                        - 0.22 * p_c
                }
                TrophicGuild::PredatorParasite => {
                    0.85 * h_c + 0.28 * self.exergy[i] + 0.25 * self.byproduct[i] - 0.22 * p_c
                }
            }
        };

        let mut best_val = eval_cell(x, y);
        let w = self.width;
        let h = self.height;
        for dy in [h - 1, 0, 1] {
            for dx in [w - 1, 0, 1] {
                let nx = (x + dx) % w;
                let ny = (y + dy) % h;
                let score = eval_cell(nx, ny);
                if score > best_val + 0.015 {
                    best_val = score;
                    best = (nx, ny);
                }
            }
        }
        best
    }
}

pub fn normalized_sequence_divergence(a: &[u16], b: &[u16]) -> f64 {
    let max_l = a.len().max(b.len());
    if max_l == 0 {
        return 0.0;
    }
    let min_l = a.len().min(b.len());
    let mut diffs = max_l - min_l;
    for i in 0..min_l {
        if a[i] != b[i] {
            diffs += 1;
        }
    }
    (diffs as f64) / (max_l as f64)
}

#[derive(Clone, Debug)]
pub struct SpeciationSummary {
    pub species_count: usize,
    pub unit_species: Vec<usize>,
    pub fst_index: f64,
    pub within_species_adami_rho: f64,
}

pub fn analyze_spatial_speciation(
    sequences: &[&[u16]],
    biomes: &[BiomeZone],
    alphabet_size: usize,
    theta_spec: f64,
) -> SpeciationSummary {
    let n = sequences.len();
    if n == 0 {
        return SpeciationSummary {
            species_count: 0,
            unit_species: Vec::new(),
            fst_index: 0.0,
            within_species_adami_rho: 1.0,
        };
    }
    if n == 1 {
        return SpeciationSummary {
            species_count: 1,
            unit_species: vec![0],
            fst_index: 0.0,
            within_species_adami_rho: 1.0,
        };
    }

    let mut labels = vec![usize::MAX; n];
    let mut species_count = 0usize;
    let mut queue = Vec::with_capacity(n);

    for start in 0..n {
        if labels[start] != usize::MAX {
            continue;
        }
        let sp_id = species_count;
        species_count += 1;
        labels[start] = sp_id;
        queue.clear();
        queue.push(start);
        let mut head = 0usize;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            for v in 0..n {
                if labels[v] == usize::MAX
                    && normalized_sequence_divergence(sequences[u], sequences[v]) <= theta_spec
                {
                    labels[v] = sp_id;
                    queue.push(v);
                }
            }
        }
    }

    let mut alpha_seqs = Vec::new();
    let mut beta_seqs = Vec::new();
    for (idx, &seq) in sequences.iter().enumerate() {
        match biomes.get(idx).copied().unwrap_or(BiomeZone::Alpha) {
            BiomeZone::Alpha => alpha_seqs.push(seq),
            BiomeZone::Beta => beta_seqs.push(seq),
        }
    }

    let h_total = mean_site_entropy(sequences, alphabet_size);
    let fst_index = if !alpha_seqs.is_empty() && !beta_seqs.is_empty() && h_total > 1e-5 {
        let w_a = (alpha_seqs.len() as f64) / (n as f64);
        let w_b = (beta_seqs.len() as f64) / (n as f64);
        let h_sub = w_a * mean_site_entropy(&alpha_seqs, alphabet_size)
            + w_b * mean_site_entropy(&beta_seqs, alphabet_size);
        ((h_total - h_sub) / h_total).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let mut weighted_rho_sum = 0.0;
    let mut weight_total = 0.0;
    for sp in 0..species_count {
        let cluster: Vec<&[u16]> = sequences
            .iter()
            .enumerate()
            .filter_map(|(i, &s)| (labels[i] == sp).then_some(s))
            .collect();
        if cluster.len() >= 2 {
            let h_norm = normalized_ensemble_entropy(&cluster, alphabet_size);
            let rho_sp = (1.0 - h_norm).clamp(0.08, 1.0);
            let w = cluster.len() as f64;
            weighted_rho_sum += w * rho_sp;
            weight_total += w;
        }
    }

    let within_species_adami_rho = if weight_total > 0.0 {
        (weighted_rho_sum / weight_total).clamp(0.08, 1.0)
    } else {
        (1.0 - normalized_ensemble_entropy(sequences, alphabet_size)).clamp(0.08, 1.0)
    };

    SpeciationSummary {
        species_count,
        unit_species: labels,
        fst_index,
        within_species_adami_rho,
    }
}

fn mean_site_entropy(sequences: &[&[u16]], alphabet_size: usize) -> f64 {
    if sequences.len() < 2 || alphabet_size < 2 {
        return 0.0;
    }
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    if max_len == 0 {
        return 0.0;
    }
    let mut counts = vec![0usize; alphabet_size];
    let mut h_sum = 0.0;
    let mut loci = 0usize;
    for pos in 0..max_len {
        counts.fill(0);
        let mut pres = 0usize;
        for s in sequences {
            if let Some(&sym) = s.get(pos) {
                counts[(sym as usize) % alphabet_size] += 1;
                pres += 1;
            }
        }
        if pres >= 2 {
            loci += 1;
            let nf = pres as f64;
            for &c in &counts {
                if c > 0 {
                    let p = (c as f64) / nf;
                    h_sum -= p * p.log2();
                }
            }
        }
    }
    if loci == 0 {
        0.0
    } else {
        h_sum / (loci as f64)
    }
}

fn normalized_ensemble_entropy(sequences: &[&[u16]], alphabet_size: usize) -> f64 {
    if sequences.len() < 2 || alphabet_size < 2 {
        return 0.0;
    }
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    if max_len == 0 {
        return 0.0;
    }
    let h_max = (alphabet_size as f64).log2();
    let mut counts = vec![0usize; alphabet_size];
    let mut norm_sum = 0.0;
    let mut loci = 0usize;

    for pos in 0..max_len {
        counts.fill(0);
        let mut pres = 0usize;
        for s in sequences {
            if let Some(&sym) = s.get(pos) {
                counts[(sym as usize) % alphabet_size] += 1;
                pres += 1;
            }
        }
        if pres >= 2 {
            loci += 1;
            let nf = pres as f64;
            let mut h_pos = 0.0;
            for &c in &counts {
                if c > 0 {
                    let p = (c as f64) / nf;
                    h_pos -= p * p.log2();
                }
            }
            let cap = h_max.min(nf.log2());
            if cap > 1e-9 {
                norm_sum += (h_pos / cap).clamp(0.0, 1.0);
            }
        }
    }
    if loci == 0 {
        0.0
    } else {
        norm_sum / (loci as f64)
    }
}

/// Renders an ASCII multi-trophic spatial map showing Composite Holobionts (`@`), Producers (`P`),
/// Syntrophs (`s`), Predators/Parasites (`X`), Byproduct pools (`~`), and Niche halos (`#`).
pub fn render_ascii_spatial_grid(
    grid: &SpatialGrid,
    unit_coords_species_guild: &[(usize, usize, usize, TrophicGuild)],
) -> String {
    let mut occ = vec![None::<(usize, TrophicGuild)>; grid.width * grid.height];
    for &(x, y, sp, guild) in unit_coords_species_guild {
        if x < grid.width && y < grid.height {
            let idx = y * grid.width + x;
            let replace = match (occ[idx], guild) {
                (None, _) => true,
                (Some((_, prev_g)), g) => (g as u8) > (prev_g as u8),
            };
            if replace {
                occ[idx] = Some((sp, guild));
            }
        }
    }

    let mut out = String::new();
    let border = "-".repeat(grid.width * 2);
    let _ = writeln!(
        out,
        "   +{}+  (Legend: @=Holobiont | P=Producer | s=Syntroph | X=Predator | #=Niche | ~=Byproduct)",
        border
    );
    for y in 0..grid.height {
        let _ = write!(out, "   |");
        for x in 0..grid.width {
            let idx = y * grid.width + x;
            if let Some((_sp, guild)) = occ[idx] {
                let ch = match guild {
                    TrophicGuild::Producer => 'P',
                    TrophicGuild::Syntroph => 's',
                    TrophicGuild::PredatorParasite => 'X',
                    TrophicGuild::Holobiont => '@',
                };
                let _ = write!(out, "{ch} ");
            } else if grid.niche[idx] >= 0.35 {
                let _ = write!(out, "# ");
            } else if grid.byproduct[idx] >= 0.25 {
                let _ = write!(out, "~ ");
            } else if grid.exergy[idx] >= 0.85 {
                let _ = write!(out, ": ");
            } else if grid.exergy[idx] >= 0.35 {
                let _ = write!(out, ". ");
            } else {
                let _ = write!(out, "  ");
            }
        }
        let _ = writeln!(out, "|");
    }
    let _ = writeln!(out, "   +{}+", border);
    out
}
