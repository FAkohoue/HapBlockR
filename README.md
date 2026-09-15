# HapBlockR

[![R-CMD-check](https://github.com/FAkohoue/HapBlockR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/HapBlockR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

HapBlockR detects linkage-disequilibrium blocks and supports haplotype
analysis, genomic prediction, parent selection, optimum-contribution
selection, genomic mating, and core-collection design.

The package is under active development. HapBlockR converts the data,
breeding objectives and constraints supplied by the user into reproducible
parent, cross and mating recommendations. Breeder-facing results include
quality control, provenance, validation and uncertainty information so that
the basis of each recommendation is explicit. Recommendation relevance is
therefore driven principally by the quality and representativeness of the
input data and by how accurately the specified parameters express the
breeding programme's objectives.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("FAkohoue/HapBlockR")
```

Required package dependencies install automatically. Some methods use
optional R packages or external tools and fail explicitly when the requested
engine is unavailable.

## Ten-minute workflow

The included deterministic values represent externally adjusted genotype
means with equal precision. They illustrate target preparation, internal
prediction and parent selection:

```r
library(HapBlockR)

data(ldx_geno)
data(ldx_snp_info)
data(ldx_blocks)
data(ldx_blues)

targets <- prepare_breeding_targets(
  data.frame(
    id = ldx_blues$id,
    trait = "YLD",
    value = ldx_blues$YLD,
    precision = 1
  ),
  input_type = "adjusted_mean",
  precision_col = "precision"
)

prediction <- run_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info = ldx_snp_info,
  blocks = ldx_blocks,
  blues = targets,
  seed = 42,
  min_reliability = 0.30,
  verbose = FALSE
)

top_blocks <- select_top_blocks(
  prediction$block_importance,
  n = min(10L, nrow(prediction$block_importance))
)

parent_values <- prediction$local_gebv[
  , top_blocks$block_id, drop = FALSE
]

parents <- select_parents_ga(
  value_matrix = parent_values,
  n_founders = 8,
  seed = 42,
  n_reps = 5,
  verbose = FALSE
)

parents$selected
validate(parents)
summary(parents)
```

This call uses coverage-only GA. When whole-genome merit must remain active
while complementary blocks are assembled, use the explicit GA+TS tool:

```r
hybrid_parents <- select_parents_ga_ts(
  value_matrix = parent_values,
  n_founders = 8,
  merit_score = prediction$gebv,
  merit_priority = 50,
  seed = 42,
  n_reps = 5,
  verbose = FALSE
)
```

Both GA tools run five independent searches by default, check hard constraints
for every result, and automatically return the feasible replicate with the
greatest complete objective. The breeder reviews the returned stability
summary; individual runs do not require manual selection.

For leakage-aware assessment before selection:

```r
cv <- cv_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info = ldx_snp_info,
  blocks = ldx_blocks,
  blues = targets,
  k = 5,
  seed = 42,
  verbose = FALSE
)

cv$pa_pooled
cv$gebv_all
validate(cv)
```

Use grouped or forward validation when families, populations, sites, years,
or prediction timing could leak information between training and validation
sets.

## Capability map

| Task | Main functions |
| --- | --- |
| Import and identity-safe caching | `read_geno()`, `read_geno_bigmemory()` |
| LD and block detection | `compute_r2()`, `compute_rV2()`, `Big_LD()`, `run_Big_LD_all_chr()` |
| Statistical phasing | `phase_with_beagle()`, `read_phased_vcf()` |
| Haplotype extraction and diversity | `infer_block_haplotypes()`, `extract_haplotypes()`, `compute_haplotype_diversity()` |
| Association and epistasis | `test_block_haplotypes()`, `estimate_diplotype_effects()`, `scan_block_epistasis()` |
| Genomic prediction | `run_haplotype_prediction()`, `cv_haplotype_prediction()` |
| Multi-trait and G×E models | `fit_multitrait_gblup()`, `build_selection_index()`, `fit_gxe_gblup()` |
| Externally analysed breeding targets | `breeding_target_types()`, `prepare_breeding_targets()` |
| Local breeding values | `estimate_marker_effects()`, `compute_local_gebv()` |
| Parent shortlisting | `truncation_selection()`, coverage-only `select_parents_ga()`, joint GA+TS `select_parents_ga_ts()`, and `select_parents_by_family()`; both GA tools support Optimal Haplotype Selection (OHS) and Optimal Population Value (OPV) strategies |
| Cross and contribution decisions | `usefulness_criterion()`, `select_parents_ocs()`, `select_parents_pareto()` |
| Plan validation and diversity | `validate_crosses_exact()`, `select_core_collection()` |
| Operational feasibility | `screen_candidate_crosses()`, `certify_mating_plan()` |
| Standards-oriented exchange | `validate_breeding_metadata()`, `build_breeding_exchange()` |
| End-to-end LD workflow | `run_ldx_pipeline()` |

The package website contains the full function reference and task-oriented
vignettes.

## Breeding-target input contract

HapBlockR starts from genotype-level estimates produced by a field-trial or
genetic-evaluation analysis outside the package. It does not accept raw plot
records, fit replicate, block, row or column effects, or calculate an
external selection index. Units are optional. Every supplied estimate must
have a standard error, posterior standard deviation, precision weight or full
sampling covariance; random-effect predictions instead require reliability
or prediction error variance and the relevant genetic variance.

| `input_type` | Meaning and treatment |
| --- | --- |
| `"adjusted_mean"` | Model-adjusted entry mean, including a Bayesian posterior adjusted mean when the analysis does not distinguish BLUE from BLUP; not a raw arithmetic mean |
| `"BLUE"` | Best Linear Unbiased Estimate of an entry fitted as fixed; retained without deregression |
| `"BLUP_identity"` | Random genotype or entry effect with covariance `I × genetic variance`; deregressed from reliability or prediction error variance |
| `"PBLUP"` | Pedigree Best Linear Unbiased Prediction with covariance `A × additive genetic variance`; requires the named numerator relationship matrix `A` |
| `"BV"` | Additive, transmissible breeding value; the fixed, identity or pedigree estimation basis must be declared |
| `"GCA"` | General combining ability in a declared tester or mate population; fixed GCA is retained and random GCA is deregressed |
| `"TGV"` | Total genetic value; additive and dominance components must be supplied separately |

External genomic Best Linear Unbiased Predictions (GBLUPs), genomic estimated
breeding values and selection-index values are rejected as target types.
HapBlockR uses the accepted external summaries as responses in its internal
marker, haplotype, block, multi-trait or genotype-by-environment models. The
resulting internal predictions feed `build_selection_index()` and become the
directionally aligned `merit_score` used in parent selection.

`prepare_breeding_targets()` reverses lower-is-better traits once, so every
downstream `model_value`, selection-index score and `merit_score` follows the
same rule: larger is better. For estimates with standard error `SE`, relative
precision is proportional to `1 / SE^2` and is normalised within each
trait-analysis group. For a pedigree BLUP, reliability is
`1 - PEV / (A_ii × additive genetic variance)`.

`build_selection_index()` provides four internal methods:

- `"smith_hazel"` for economic weights;
- `"pesek_baker"` for a deterministic desired-gain index;
- `"dgsi"` for replicated Desired-Gain Selection Index optimisation through
  DesiredGainR; and
- `"qgsi"` for a Quadratic Genomic Selection Index through DesiredGainR using
  an explicit, symmetric matrix of squared and cross-product weights.

Trait direction and objective magnitude have one unambiguous contract.
`directions` declares which traits increase or decrease, while
`economic_weights` and `desired_gains` contain non-negative magnitudes after
orientation to the favourable direction. HapBlockR rejects negative
objectives instead of silently changing their signs.

For DGSI, desired gains are expressed in candidate standard deviations,
whether or not trait scaling is requested. For Pesek-Baker, they are expressed
in original trait units. Divide an original-unit DGSI target by the candidate
standard deviation of that trait before passing it to `desired_gains`.

The DGSI engine selects its best replicate automatically using its declared
holdout or validation rule and reports coefficient, rank and selected-set
stability. HapBlockR uses DesiredGainR's model-expected transmitted response in
original trait units and reports the realised selected-set differential
separately in candidate standard-deviation units. DGSI and QGSI use the
selection intensity for the number of candidates actually selected. QGSI
model-expected gains use the total linear-plus-quadratic index variance and are
converted from DesiredGainR's analysis space back to the original trait units.
QGSI has no single global coefficient vector: HapBlockR reports its linear
weights, quadratic-weight matrix and candidate-specific contributions
separately. All four methods return the same coefficient-table columns, with
`NA` only where a field does not apply to that method.

`dgsi_control` and `qgsi_control` accept exact DesiredGainR argument names;
partial names and overrides of HapBlockR's structural arguments are rejected.
QGSI's `Gamma` is the covariance of genomic predictions, not an automatic
substitute for the supplied genetic covariance. Supply it in original trait
units through `qgsi_control`, or let DesiredGainR estimate it from the candidate
or reference genomic predictions, optionally using a relationship matrix.
QGSI weights must refer to the oriented and, if requested, scaled trait space.

For DGSI, `coefficients_original_units` includes the direction and scale
conversion needed to combine original-unit marker, haplotype or block effects.
`score_intercept` accounts for reference centring when reconstructing candidate
scores. The unmodified `engine_result` remains available for DesiredGainR's
comparison and diagnostic tools.

When `prepare_breeding_targets()` supplies a full sampling covariance,
`fit_multitrait_gblup()` uses it as the complete record-error covariance by
default. Set `sampling_covariance_mode = "sampling_plus_residual"` only when
the model requires a separate residual nugget in addition to known sampling
error. Record keys and diagonal precision must agree exactly enough to prevent
silent loss or duplication of uncertainty.

## Common result contract

Decision-critical methods return objects inheriting from
`hapblockr_result`. The versioned contract records:

- method, call, normalised parameters, and seed;
- immutable sample and variant IDs and SHA-256 input hashes;
- transformation history and software versions;
- quality-control gates, warnings, fallbacks, and exclusions;
- decision and uncertainty tables; and
- machine-readable validation status.

Use:

```r
validate(result)
print(result)
summary(result)
as.data.frame(result)
plot(result)
```

A failed validation gate must be resolved before a result is promoted to a
breeding recommendation.

## Beagle phasing

HapBlockR integrates a user-supplied Beagle 5.x JAR, allowing the programme to
retain explicit control of the external-tool version and licence. Download it
from the official Beagle site and configure the path:

```r
options(HapBlockR.beagle_jar = "/absolute/path/to/beagle.jar")

phased <- phase_with_beagle(
  input_vcf = "genotypes.vcf.gz",
  out_prefix = "results/genotypes_phased",
  nthreads = 2,
  seed = 42,
  min_genotype_concordance = 0.99,
  return_details = TRUE
)

phased$quality_control
phased$provenance
```

The path may instead be supplied through `beagle_jar`,
`HAPBLOCKR_BEAGLE_JAR`, or a file named `beagle.jar` beside
`out_prefix`. HapBlockR restricts the integration to one or two threads,
records Java and JAR provenance, verifies Beagle 5.x, and checks sample,
variant, allele, and observed-genotype identity after phasing.

Official Beagle download and licence:
<https://faculty.washington.edu/browning/beagle/beagle.html>

## Data scale and memory

Memory behaviour depends on the source and backend:

- GDS, PLINK BED, and `bigmemory` pathways can provide file-backed or
  subset-oriented access;
- text dosage, HapMap, and ordinary VCF pathways may materialise substantial
  objects in memory; and
- downstream algorithms can still require dense matrices even when import is
  streamed.

Measure wall time and peak resident memory on data representative of the
intended programme. HapBlockR provides file-backed, subset-oriented and dense
analysis pathways; the selected pathway determines whether a full matrix is
materialised for a particular operation.

## Interpreting HapBlockR results

All breeder-facing merit inputs follow one direction: larger values mean
greater breeding merit. Declare lower-is-better traits during target
preparation or index construction; do not reverse their signs a second time.
Standardised superior-candidate filtering is available through
`min_sel_mode = "sd_above_mean"`; `min_sel_value = 0` retains candidates at
or above the mean and `1` requires at least one SD superiority.
Use `min_sel_mode = "relaxed_pool"` only when the programme deliberately
admits some candidates below the mean for complementarity or diversity. The
former name `"sd_below_mean"` remains as a deprecated alias.

- HapBlockR accepts phased data and also provides optional statistical phasing
  through `phase_with_beagle()`. Use phasing when the breeding question
  depends on the parental chromosome carrying an allele; dosage-based
  workflows remain available when that distinction is unnecessary.
- Haplotype inference, `hap1`/`hap2` representations, Beagle integration, and
  compiled r-squared or rV-squared kernels operate on diploid, biallelic data.
  Dosage-centred relationship and marker-effect calculations can accept
  alternative ploidy values.
- `run_haplotype_prediction()` supports single- and multiple-trait inputs.
  Use `fit_multitrait_gblup()` when genetic and residual covariance should be
  estimated jointly, and `fit_gxe_gblup()` for reaction-norm predictions.
  Both models report REML likelihoods with residual degrees of freedom and
  prediction error variances adjusted for fitted fixed effects. Compare REML
  likelihoods only between models with the same fixed-effect design.
- HapBlockR controls population structure and genomic relatedness in its
  adjusted association models, applies multiple-testing procedures, compares
  effects across populations, and provides within- and between-block
  fine-mapping tools. A reproducible haplotype association can identify a
  segment that harbours a functional gene or causal variant. The strength of
  a causal conclusion then depends on positional resolution, replication and
  functional evidence.
- Prediction relevance is determined by the phenotype quality, genotype and
  haplotype representation, training population, target environments,
  validation design and parameters supplied by the user. Grouped and forward
  validation quantify transfer to new families or cycles.
- HapBlockR supports representative future populations through
  family- or genetic-cluster selection, GA selection of complementary
  favourable haplotypes, core-collection design, relationship control and
  environment-aware modelling.
- Reliability, uncertainty and validation gates show how strongly the data
  support each reported recommendation.

See the vignettes for method-specific assumptions and diagnostics.

## Breeder's guide

The versioned, source-backed guide covers nine decision tools, their
variants, a worked numerical decision, data-quality and feasibility checks,
interpretation guidance, references, and a sign-off template. It is
distributed in accessible PDF and editable Word formats:

```r
open_breeder_guide(format = "pdf")
open_breeder_guide(format = "docx")
```

Its source is
`inst/guide/HapBlockR_Breeder_Guide.md`. The locked, reproducible builder is
`tools/build_breeder_guide_accessible.cjs`:

```sh
corepack enable
corepack prepare pnpm@11.9.0 --activate
pnpm --dir tools install --frozen-lockfile
pnpm --dir tools exec playwright install chromium
node tools/build_breeder_guide_accessible.cjs
python tools/build_breeder_guide_docx.py
```

The build produces a tagged A4 PDF with a document outline, semantic table
headers, internal links, British-English language metadata and page numbers.
`tools/build_breeder_guide_docx.py` produces the corresponding editable Word
edition with real headings, lists, table headers and page numbering.
Continuous integration rebuilds and verifies both editions.

## Reproducibility

For an auditable analysis:

```r
set.seed(42)
packageVersion("HapBlockR")
sessionInfo()
```

Retain the result object, source data release, input hashes, external-tool
checksums, exclusions, warnings, validation output, and signed decision
record. Do not place confidential germplasm or phenotype records in public
issues or repositories.

## Citation

Use the installed citation so that the package version is current:

```r
citation("HapBlockR")
```

Machine-readable citation metadata are provided in `CITATION.cff`. A DOI has
not yet been assigned; do not cite an invented or provisional DOI.

## Contributing and support

Read [CONTRIBUTING.md](CONTRIBUTING.md), the
[code of conduct](CODE_OF_CONDUCT.md), [security policy](SECURITY.md), and
[support guidance](SUPPORT.md). Report reproducible bugs through the
[GitHub issue tracker](https://github.com/FAkohoue/HapBlockR/issues).

## Licence

HapBlockR is distributed under the MIT licence. External tools and optional
dependencies retain their own licences.
