# Simulate GA-Selected vs. Truncation-Selected Founders Over Generations

Runs recurrent-selection forward simulation for two founder sets –
typically the output of
[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md)
and
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md)
– and tracks realised genetic gain (mean/max GEBV of the breeding
population) over `n_generations`. This is how HapSelect's own
`localGEBV_vs_TS_simulation()`/ `Haplotype_vs_TS_simulation()`
demonstrate that a GA-optimised founder set outperforms plain truncation
selection, and this function is a wrapper around the same underlying
engine HapSelect uses for that comparison: the genomicSimulation package
(<https://github.com/vllrs/genomicSimulation>). See the file header of
`R/forward_simulation.R` for how the wrapper is structured (each scheme
runs as its own isolated genomicSimulation session, since
genomicSimulation itself has only one active simulation at a time).

## Usage

``` r
ga_vs_ts_simulation(
  hap1,
  hap2,
  snp_info,
  snp_effects,
  ga_selected = NULL,
  ts_selected = NULL,
  schemes = NULL,
  blocks = NULL,
  n_generations = 10L,
  pop_size = 100L,
  recomb_rate = 1e-08,
  selection_intensity = 0.2,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- hap1, hap2:

  Numeric matrices (SNPs x individuals), values 0/1 – e.g.
  `read_phased_vcf()$hap1` / `$hap2`. Column names are individual IDs;
  row names/order must match `snp_info`.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`, same row order as
  `hap1`/`hap2`.

- snp_effects:

  Named numeric vector of per-SNP additive effects (from
  [`backsolve_snp_effects`](https://FAkohoue.github.io/HapBlockR/reference/backsolve_snp_effects.md)
  or
  [`estimate_marker_effects`](https://FAkohoue.github.io/HapBlockR/reference/estimate_marker_effects.md)),
  names matching `snp_info$SNP`. GEBVs are calculated by
  genomicSimulation on the same centred-dosage scale used elsewhere in
  HapBlockR (\\(x - 2p) \cdot \alpha\\, via
  [`genomicSimulation::change.eff.set.centres()`](https://rdrr.io/pkg/genomicSimulation/man/change.eff.set.centres.html)),
  so trajectories are directly comparable to
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/HapBlockR/reference/run_haplotype_prediction.md)
  output.

- ga_selected:

  Character vector of founder IDs (from `select_parents_ga()$selected`),
  must be a subset of `colnames(hap1)`. Ignored (with a message) if
  `schemes` is supplied. Together with `ts_selected`, this is the
  original 2-scheme calling convention – kept for backward
  compatibility; it builds
  `schemes = list(GA = list(founders = ga_selected, mating_scheme = "truncation"), TS = list(founders = ts_selected, mating_scheme = "truncation"))`
  internally.

- ts_selected:

  Character vector of founder IDs (from
  `truncation_selection()$selected`), same requirement as `ga_selected`.

- schemes:

  Optional named list, one entry per scheme to simulate (overrides
  `ga_selected`/`ts_selected` when supplied; each list name is used as
  the scheme's label in the returned `summary`). Each entry is itself a
  list with:

  `founders`

  :   Required. Character vector of founder IDs (subset of
      `colnames(hap1)`, \>= 2 individuals).

  `mating_scheme`

  :   One of `"truncation"` (default if omitted – random-mate the
      current generation, keep the top `selection_intensity` fraction by
      GEBV; the original, only behaviour of this function), `"ocs"`
      (each generation,
      [`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
      picks an optimum-contribution mating plan from the current
      population, executed via
      [`genomicSimulation::make.targeted.crosses()`](https://rdrr.io/pkg/genomicSimulation/man/make.targeted.crosses.html);
      requires `blocks`), or `"uc"` (each generation,
      [`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
      ranks every candidate pair by predicted progeny usefulness and the
      top-ranked pairs are crossed; requires `blocks`).

  Any other named element

  :   Passed through to the underlying per-generation call:
      [`select_parents_ocs()`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md)
      for `"ocs"` (e.g. `n_crosses`, `max_contrib_per_parent`,
      `target_degree`), or
      [`usefulness_criterion()`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md)
      for `"uc"` (e.g. `selected_proportion`), plus the `"uc"`-only
      `n_crosses` (how many top-ranked pairs are actually crossed each
      generation; default `pop_size`, one offspring per pair).

  See the "OCS/UC-informed rapid-cycling extension" section of the
  `R/forward_simulation.R` file header for the full design and its
  "Verification status" note.

- blocks:

  LD-block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/HapBlockR/reference/run_Big_LD_all_chr.md)
  (or
  [`tune_LD_params`](https://FAkohoue.github.io/HapBlockR/reference/tune_LD_params.md)).
  Required if any scheme uses `mating_scheme = "uc"` (each generation's
  per-block local GEBVs are computed via
  [`compute_local_gebv`](https://FAkohoue.github.io/HapBlockR/reference/compute_local_gebv.md),
  which needs the block table). Not needed by `"ocs"` (which builds a
  whole-genome GRM directly from dosage, no block decomposition
  involved) or `"truncation"`. Validated up front, before any simulation
  work is done, whenever at least one scheme requests `"uc"`.

- n_generations:

  Integer. Number of generations to simulate. Default `10L`.

- pop_size:

  Integer. Offspring produced per generation, per scheme. Default
  `100L`.

- recomb_rate:

  Numeric. Used to convert `snp_info$POS` (bp) into the genetic map (cM)
  that genomicSimulation's own crossing functions consume:
  `pos_cM = POS_bp * recomb_rate * 100`. Default `1e-8` (~1 cM/Mb, a
  common genome-wide average approximation). Rescale for a species with
  a known genetic map:
  `recomb_rate = map_length_cM / 100 / genome_length_bp`.

- selection_intensity:

  Numeric in (0,1\] or `NULL`. Fraction of each generation's offspring
  retained as parents for the next generation (recurrent truncation
  selection within each scheme, via
  [`genomicSimulation::break.group.by.GEBV()`](https://rdrr.io/pkg/genomicSimulation/man/break.group.by.GEBV.html)).
  Default `0.2` (top 20%). Set `NULL` to mate every offspring generation
  in full (no within-scheme truncation).

- seed:

  Integer or `NULL`. Random seed for reproducibility
  ([`set.seed()`](https://rdrr.io/r/base/Random.html);
  genomicSimulation's own crossing functions use R's random number
  generator).

- verbose:

  Logical. Print per-generation progress. Default `FALSE`.

## Value

Named list:

- `summary`:

  Data frame: `scheme` (scheme label – "GA"/ "TS" under the legacy
  2-scheme convention, or the `schemes` list's own names), `generation`
  (0 = founders themselves), `mean_gebv`, `max_gebv`, `sd_gebv`,
  `n_pop`.

- `final`:

  Named list (one element per scheme, same names as `summary$scheme`) of
  lists with `hap1`/`hap2` (reconstructed from genomicSimulation's
  internal storage; markers x individuals, row names = marker names in
  genomicSimulation's own marker order, which may differ from
  `snp_info`'s input row order) and `gebv` (named numeric vector) for
  the final generation's breeding population of each scheme, for further
  analysis (e.g. diversity metrics).

- `ga_final`, `ts_final`:

  Only present when called via the legacy `ga_selected`/`ts_selected`
  convention (i.e. `schemes` not supplied) – identical content to
  `final$GA`/`final$TS`, kept as top-level names unchanged from previous
  releases so existing calling code does not break.

## Installation of genomicSimulation

genomicSimulation is not on CRAN. Install it from a source release:

1.  Download the `.tar.gz` from
    <https://github.com/vllrs/genomicSimulation/releases>.

2.  `install.packages("path/to/genomicSimulation_x.y.z.tar.gz", repos = NULL)`

On a shared HPC system, compile on the same CPU architecture/node type
you will run jobs on (see your cluster's documentation for details –
e.g. UQ Bunya requires compiling on an "epyc3" node to avoid "illegal
instruction" errors at runtime on other nodes).

**Requires phased haplotypes.** Unphased 0/1/2 dosage cannot support a
block-preserving meiosis simulation – without cis/trans information at
heterozygous sites, simulated recombination would not respect the LD
blocks the founder sets were selected to stack. Use
[`read_phased_vcf`](https://FAkohoue.github.io/HapBlockR/reference/read_phased_vcf.md)
or another phased source for `hap1`/ `hap2`.

## See also

[`select_parents_ga`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ga.md),
[`truncation_selection`](https://FAkohoue.github.io/HapBlockR/reference/truncation_selection.md),
[`select_parents_ocs`](https://FAkohoue.github.io/HapBlockR/reference/select_parents_ocs.md),
[`usefulness_criterion`](https://FAkohoue.github.io/HapBlockR/reference/usefulness_criterion.md),
[`plot_ga_vs_ts_simulation`](https://FAkohoue.github.io/HapBlockR/reference/plot_ga_vs_ts_simulation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
phased <- read_phased_vcf("mydata_phased.vcf.gz")
res    <- run_haplotype_prediction(phased$dosage, phased$snp_info, blocks,
                                   blues = blues_vec)
top    <- select_top_blocks(res$block_importance, n = 15)
vmat   <- res$local_gebv[, top$block_id, drop = FALSE]
ga     <- select_parents_ga(vmat, n_founders = 20, seed = 1)
ts     <- truncation_selection(res$gebv, n_founders = 20)

# Legacy 2-scheme convention (unchanged from previous releases):
sim <- ga_vs_ts_simulation(
  phased$hap1, phased$hap2, phased$snp_info, res$snp_effects,
  ga_selected = ga$selected, ts_selected = ts$selected,
  n_generations = 10, seed = 1
)
plot_ga_vs_ts_simulation(sim)

# New: 4-way comparison including OCS- and UC-informed rapid cycling.
G <- compute_haplotype_grm(phased$dosage)
sim4 <- ga_vs_ts_simulation(
  phased$hap1, phased$hap2, phased$snp_info, res$snp_effects,
  schemes = list(
    GA  = list(founders = ga$selected, mating_scheme = "truncation"),
    TS  = list(founders = ts$selected, mating_scheme = "truncation"),
    OCS = list(founders = ts$selected, mating_scheme = "ocs",
               n_crosses = 20),
    UC  = list(founders = ts$selected, mating_scheme = "uc",
               selected_proportion = 0.1)
  ),
  blocks = blocks, n_generations = 10, seed = 1
)
plot_ga_vs_ts_simulation(sim4)
} # }
```
