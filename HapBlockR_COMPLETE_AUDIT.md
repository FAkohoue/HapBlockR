# Complete audit of HapBlockR

**Package version:** 0.3.12.9000  
**Audit date:** 25 July 2026  
**Audit scope:** source code, public functions, compiled code, data input, statistical methods, parent and mating decisions, documentation, six vignettes, breeder's guide, examples, tests, continuous integration, package metadata and release readiness  
**Audit standard:** decision-grade use in an operational plant-breeding programme

## Executive assessment

HapBlockR is an unusually broad development-stage package. It connects linkage disequilibrium (LD), haplotype definition, association testing, genomic prediction, parent selection, optimal contribution selection (OCS), mate allocation, diversity management and forward simulation. The package contains 72 exported functions, 1,069 `testthat` cases and several compiled kernels. The available test paths completed without a failed expectation.

The package is nevertheless **not yet ready to be presented as a decision-grade breeding platform**. Its current breadth exceeds the strength of its validation, provenance controls and release engineering. Several public claims are stronger than the corresponding implementation:

- `select_parents_ga()` does not select repeated genetic-algorithm runs using the documented complete objective. It also treats the requested relatedness ceiling and parent count as soft penalties.
- the `optiSel` path can return an empty mating plan after solver failures while the tests still pass;
- cross-validation is random individual-level cross-validation, and the documented out-of-fold genomic estimated breeding values (GEBVs) are not returned;
- genotype cache invalidation and sample relabelling are unsafe for some backends;
- the advertised out-of-core behaviour is format-dependent rather than universal;
- OpenMP kernels hard-code four threads, contrary to the Comprehensive R Archive Network (CRAN) limit of two, and cannot be interrupted;
- three of six vignettes disable evaluation globally, and 28 test cases were skipped in the audited environment;
- the source package is approximately 27.6 MB, contains Java archives without the required source layout, and is not currently CRAN-ready;
- the breeder's guide is attractive and accessible but is not reproducibly built, is versionless, contains internal inconsistencies and provides no worked numerical case.

The recommended course is to stabilise the decision contract before adding more methods. First correct the selection objectives, enforce hard constraints, expose validation status, repair data provenance and make every claimed integration fail loudly in dedicated integration tests. Then add breeder-facing uncertainty, operational crossing constraints and environment-aware validation.

## Audit boundary and evidence

The audit combined static source inspection, package inventory, full `testthat` execution, package checking, documentation comparison and visual inspection of every page of the breeder's guide.

| Item | Observed |
|---|---:|
| Exported functions | 72 |
| Registered S3 methods | 9 |
| Top-level R functions | 201 |
| R source files | 27 |
| R source lines | 24,659 |
| C/C++ source files | 7 |
| Test files | 29 |
| Test lines | 14,073 |
| `test_that()` cases | 1,069 |
| Expectations | 2,203 |
| Vignettes | 6 |
| Vignette source lines | 3,168 |
| Vignette code chunks | 121 |
| Rd help pages | 85 |
| Rd pages with examples | 54 |
| Rd pages without examples | 31 |
| Breeder's guide pages | 20 |

The dynamic results apply to the installed dependencies and Windows environment used for this audit. Beagle end-to-end operation, all optional prediction engines and whole-genome-scale performance were not demonstrated because the applicable tests were skipped or no acceptance-scale fixture was supplied. A zero-failure unit-test result must therefore not be interpreted as validation of those paths.

The final package check reported **0 errors, 0 warnings and 4 notes**. Two notes are attributable to package content: CRAN incoming feasibility identified the development version, `Remotes` and redistributed Java archives without a top-level `java` directory; `estimate_diplotype_effects()` took 9.25 seconds in its example. Two notes were audit-environment artefacts: the audit's temporary top-level directory and a residual check-library lock. The package tests took 402 seconds, and vignette rebuilding took 137 seconds.

## Priority system

- **P0 — release blocker:** can produce an incorrect breeding decision, violates a documented hard guarantee, corrupts provenance or prevents a lawful/reliable release.
- **P1 — high:** materially limits scientific validity, operational use, reproducibility or user trust.
- **P2 — medium:** impairs maintainability, usability, performance evidence or documentation quality.
- **P3 — low:** consistency, presentation or governance improvement.

## Findings register

### P0 — Correct the genetic-algorithm objective and hard constraints

**Evidence.** In [`R/parent_selection.R`](R/parent_selection.R), the internal run returns a coverage-only value as `fitness`, although the optimisation can also include whole-genome merit, coancestry, a relatedness-ceiling penalty and a cardinality penalty. The outer function chooses the best replicate using that incomplete value. Consequently, with the default repeated runs, the selected run need not be the run that best satisfies the configured objective.

The same function documents `target_degree` as a relatedness ceiling. The implementation converts violations into a finite quadratic penalty. An arbitrarily small violation can therefore remain preferable to a feasible solution. The required number of parents is also penalised rather than repaired or enforced; the function can warn and return the wrong count.

**Required correction.**

1. Define one canonical objective function and return its complete, decomposed value from every run.
2. Rank repeated runs by that complete objective after a deterministic feasibility check.
3. Treat `n_founders` and a documented relationship ceiling as hard constraints. Repair infeasible chromosomes during the search or reject infeasible final solutions.
4. Return `feasible`, `constraint_violations`, `objective_components`, `termination_reason`, `seed`, `run_id` and a candidate-level audit table.
5. Add property tests asserting the exact parent count, the ceiling within a stated tolerance and consistency between the selected run and the complete objective.
6. Amend the README and help pages immediately if the ceiling remains soft.

**Release gate.** No P0 selection function may return a result labelled successful when a documented hard constraint is violated.

### P0 — Make the `optiSel` integration fail safely

**Evidence.** In [`R/ocs.R`](R/ocs.R), the phenotype table supplies missing sex values. During the audit, `optiSel::noffspring()` failed with `missing value where TRUE/FALSE needed`, and `optiSel::matings()` rejected missing sex values. The wrapper substituted rounded contributions and then returned an empty mating plan. The corresponding assertions in [`tests/testthat/test-ocs.R`](tests/testthat/test-ocs.R) require only that `ok` is logical, not that it is `TRUE`, and check a non-empty mating plan only conditionally.

`n_parents_max` is applied by zeroing small contributions after optimisation, without re-solving or consistently renormalising the constrained problem. Selfing and repeated-cross restrictions are also filtered after allocation without refilling the plan. These operations can invalidate contribution totals and the optimality or feasibility of the returned solution.

**Required correction.**

1. Supply valid sex categories or use the documented `optiSel` mode that does not require sex.
2. Treat solver, offspring-allocation and mating-allocation failures as errors by default.
3. Permit an explicit `allow_heuristic_fallback = TRUE` mode only if the returned object is marked `approximate`, contains the original error and passes all postconditions.
4. Integrate parent-count, selfing, repeated-cross and capacity constraints into the optimisation or re-solve after any restriction.
5. Assert contribution sum, non-negativity, parent count, cross count, sex compatibility, repeat limits, selfing policy and relatedness target in tests.
6. Add an integration job that installs a fixed supported `optiSel` version and requires a non-empty valid plan.

### P0 — Repair genotype identity and cache provenance

**Evidence.** In [`R/read_geno.R`](R/read_geno.R), user-supplied `sample_ids` replace the backend sample identifiers. The Genomic Data Structure (GDS) reader later uses those labels as physical `sample.id` keys. Relabelling a GDS input can therefore break or misdirect retrieval. The VCF and HapMap cache check verifies that an existing GDS file opens but does not fingerprint the source path, size, modification time or content. A changed source can silently reuse stale genotype data.

**Required correction.**

1. Store immutable physical identifiers separately from display labels.
2. Persist a cache manifest containing canonical source path, size, modification time, strong content hash, import options, package version and backend version.
3. Refuse a mismatched cache or rebuild it explicitly.
4. Record sample and variant order before and after every transformation.
5. Validate uniqueness of sample and variant identifiers, finite and ordered positions, chromosome labels, dosage range and ploidy.
6. Return a machine-readable import report listing discarded, relabelled and reordered records.

### P0 — Resolve CRAN and redistribution blockers

**Evidence.** The source-package check reports:

- a source tarball of approximately 27,564,376 bytes;
- large version components in `0.3.12.9000`;
- an unknown `DESCRIPTION` field, `Remotes`;
- Java `.jar`/`.class` content without a top-level `java` source directory.

The repository also tracks compiled objects and a package dynamic library in [`src`](src), two Beagle Java archives under [`inst/extdata`](inst/extdata), and very large graphics, including a schematic Scalable Vector Graphics (SVG) file of approximately 13.5 MB. No adjacent Beagle licence, source bundle or checksum manifest was found. The two archives differ, and the documentation refers to an older 2021 release.

CRAN policy requires source packages to avoid executable binary content, requires source for redistributed free and open-source Java components in a top-level `java` directory, asks packages to remain below approximately 10 MB and limits parallel code to no more than two cores. See the [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html).

**Required correction.**

1. Decide whether Beagle is an external prerequisite or redistributed component.
2. Prefer a secure, opt-in downloader that verifies an approved URL and checksum, records the version and honours licensing; otherwise provide the complete corresponding source and notices in the CRAN-prescribed layout.
3. Add `inst/COPYRIGHTS` or equivalent third-party notices and document licence compatibility.
4. Remove tracked `.o` and `.dll` artefacts from source control and ignore generated build products.
5. optimise or regenerate large SVG assets without embedded raster content.
6. remove `Remotes` from the CRAN build or move remote-only dependencies behind a documented development profile.
7. use a CRAN-compatible development version and keep the built source package below 10 MB.
8. replace stale claims in [`cran-comments.md`](cran-comments.md), including the statement that no binary executable content is shipped.

### P0 — Fix compiled parallel execution

**Evidence.** [`src/ld_core.cpp`](src/ld_core.cpp) contains three `num_threads(4)` directives. They ignore user configuration and exceed CRAN's normal two-core limit. Long loops do not call `Rcpp::checkUserInterrupt()`. The configure logic assumes OpenMP availability by operating system instead of compiling and linking a feature probe.

**Required correction.**

1. default to one thread and cap automatic use at two;
2. accept a validated `n_threads` parameter and honour standard environment controls;
3. add periodic user-interrupt checks outside parallel regions;
4. use a compile-and-link OpenMP probe and a tested serial fallback;
5. add deterministic equality tests across one and two threads;
6. benchmark nested parallelism and prevent oversubscription.

### P1 — Return the promised cross-validation predictions

**Evidence.** [`R/haplotype_analysis.R`](R/haplotype_analysis.R) documents `gebv_all` as out-of-fold GEBVs. The function constructs an out-of-fold vector but omits it from the returned object. The accuracy summary takes an unweighted mean of fold-specific Pearson correlations. Fit failures can be converted to `NULL` and silently omitted.

**Required correction.**

1. return a row-level out-of-fold table with sample, repetition, fold, observation, prediction and exclusion reason;
2. report pooled out-of-fold correlation and error metrics, with confidence intervals;
3. if fold-wise correlations are retained, use an explicitly justified weighting or Fisher transformation;
4. fail when a requested model does not fit unless the caller selects a documented partial-results policy;
5. add family-, cluster-, population-, year-, site- and generation-aware fold assignment;
6. add forward prediction, CV1/CV2 and nested tuning modes.

### P1 — Distinguish truly out-of-core and in-memory formats

**Evidence.** The numeric reader preallocates a full genotype matrix; HapMap import uses a full `fread`; fallback VCF and phased-VCF paths read the whole file. GDS, PLINK/BED and `bigmemory` paths provide the stronger streaming behaviour. The documentation currently generalises the stronger behaviour to all formats.

**Required correction.**

1. publish a backend capability table covering memory model, random access, multiallelic handling, compression, thread support and conversion cost;
2. label in-memory paths clearly and estimate required memory before import;
3. reject jobs above a configurable memory budget;
4. implement iterator-based fallbacks or require conversion to a supported out-of-core format;
5. include peak resident set size in performance tests.

### P1 — Correct multiallelic and genotype parsing semantics

**Evidence.** The VCF documentation says that the first alternative allele is used. The parser counts every non-zero allele as alternative, thereby collapsing all alternative alleles. The automatic GDS path instead retains biallelic records only. The same source can therefore produce backend-dependent results.

**Required correction.**

1. define one explicit policy: reject, split or select multiallelic variants;
2. apply it identically across backends;
3. retain reference/alternative allele identity, ploidy, phasing and missingness;
4. test phased, unphased, haploid, diploid, polyploid and multiallelic records;
5. make any lossy transformation opt-in and report every affected variant.

### P1 — Remove session-wide side effects

**Evidence.** The `bigmemory` backend sets `options(bigmemory.typecast.warning = FALSE)` without restoring it in one execution path.

**Required correction.** Use `withr::local_options()` or an `on.exit()` restoration in every path. Add a test that snapshots and compares all modified options before and after both successful and failing calls.

### P1 — Repair epistasis fine mapping and eliminate skip-on-error tests

**Evidence.** The least absolute shrinkage and selection operator (LASSO) path can pass degenerate columns to `glmnet::cv.glmnet()` and fail with `non-conformable arrays`. The tests catch the error and call `skip()`, converting a regression into an apparent pass. Six block-by-block epistasis tests also skip because their own fixture contains no valid haplotype columns.

**Required correction.**

1. remove zero-variance and aliased features before scaling;
2. calculate an effective fold count from the usable sample size and response distribution;
3. return an explicit `not_estimable` status only for scientifically valid insufficiency;
4. make the fixtures satisfy their stated preconditions;
5. never skip after an unexpected model error;
6. test coefficient recovery and false-positive behaviour on simulated data.

### P1 — Do not disguise degenerate inference as adjusted inference

**Evidence.** `test_block_haplotypes()` can replace a failed genomic relationship matrix (GRM) with an identity matrix and continue after a warning. Too-small or monomorphic subsets can also be skipped.

**Required correction.**

1. fail by default when the requested population-structure adjustment cannot be estimated;
2. offer an explicit `fallback = "unadjusted"` mode;
3. mark every row with `qc_status`, `model_fitted`, `adjustment_used`, effective sample size and exclusion reason;
4. distinguish a scientific null result from a model that was not estimable.

### P1 — Correct the core-collection distance guarantee

**Evidence.** [`R/core_collection.R`](R/core_collection.R) derives `D = diag(G) + diag(G) - 2G`, which is squared Euclidean distance, and does not take the square root. Squared Euclidean distance is not generally a metric. The documented Gonzalez two-approximation guarantee therefore does not follow. User-supplied distance matrices are not comprehensively checked for symmetry, non-negativity, finiteness or the triangle inequality.

**Required correction.**

1. validate and, where justified, bend the GRM before deriving distance;
2. use Euclidean distance, `sqrt(pmax(D, 0))`, when the theoretical guarantee is invoked;
3. validate distance-matrix properties and expose the tolerance;
4. describe the method as a heuristic unless its assumptions are verified;
5. test the guarantee on metric fixtures and reject invalid inputs.

### P1 — Separate gametic haplotypes from unphased diplotypes

**Evidence.** Several workflows encode unphased dosage strings and subsequently call them haplotypes. These are diplotype summaries, not observed gametic haplotypes. Population differentiation calculated by treating each unphased observation as one copy is not conventional haplotype F-statistic estimation.

**Required correction.**

1. introduce separate classes for phased haplotypes, inferred haplotypes, diplotypes and dosage patterns;
2. require a compatible class for recombination- or phase-dependent methods;
3. rename statistics that use dosage-pattern frequencies or require phased input;
4. propagate phasing probability or uncertainty instead of treating inferred states as certain;
5. provide a conversion audit and prevent silent class coercion.

### P1 — Validate Beagle as an external scientific dependency

**Evidence.** Eleven end-to-end phasing tests were skipped because Beagle was not functional. The package contains two different archives, but the default function does not consistently discover the installed package copy. Output does not provide switch-error, imputation accuracy, missingness-change or allele-concordance metrics, and it does not capture the exact Java and Beagle version/checksum.

The official Beagle site currently lists Beagle 5.5, dated 27 February 2025, with Java 8 requirements, GNU General Public License version 3 terms and source availability: [Beagle software](https://faculty.washington.edu/browning/beagle/beagle.html).

**Required correction.**

1. support one declared Beagle version range and verify it in a dedicated integration job;
2. search only explicit, documented locations and report the selected binary;
3. capture executable path, checksum, command, Java version, Beagle version, reference panel and random seed;
4. validate sample/variant/allele identity before and after phasing;
5. add truth-set tests reporting switch error and dosage/imputation accuracy;
6. return a structured quality-control report and fail below user-defined thresholds.

### P1 — Add uncertainty to breeder-facing decisions

**Evidence.** Most selection outputs are point estimates. They do not consistently report prediction error variance (PEV), reliability, bootstrap or posterior intervals, probability of favourable state, sensitivity to thresholds, or probability that one cross dominates another.

**Required correction.**

- return uncertainty for local GEBVs, haplotype effects, usefulness criterion and cross rankings;
- propagate effect and phasing uncertainty through cross simulation;
- display expected gain together with downside risk, not only rank;
- provide threshold-sensitivity and leave-one-environment/population analyses;
- require minimum reliability and effective sample size before an output can be promoted to a recommendation.

### P1 — Strengthen multi-trait and genotype-by-environment modelling

**Evidence.** The multi-trait workflow fits traits independently. It does not estimate a multivariate genetic covariance structure, economic selection index or desired-gain index. Environment stability is largely represented by Finlay–Wilkinson summaries; spatial field structure and factor-analytic genotype-by-environment (G×E) models are absent.

**Required correction.**

1. add multivariate genomic best linear unbiased prediction (GBLUP) with covariance diagnostics;
2. add economic and desired-gain indices with unit and direction checks;
3. support environment-aware cross-validation and reaction-norm or factor-analytic models;
4. integrate phenotype quality control, trial design, spatial correction and environmental covariates;
5. report environment-specific and overall recommendation stability.

### P1 — Add operational mating constraints

**Evidence.** Current selection criteria cover merit, relatedness and selected mating restrictions, but not the principal feasibility constraints of an operational programme.

**Required correction.** Support sex and role, flowering synchrony, fertility, reciprocal effects, seed availability, crossing capacity by parent and period, forbidden/required pairs, heterotic groups, subpopulation quotas, quarantine constraints and minimum family size. Every plan should include a feasibility certificate and a list of binding constraints.

### P1 — Build a common decision and provenance object

**Evidence.** Outputs are heterogeneous lists and data frames. They do not share a schema, schema version, input hashes, call, package/external-tool versions, seed, session information, validation status or complete warning/error ledger.

**Required correction.** Introduce a versioned `hapblockr_result` contract with:

- method and schema version;
- immutable sample and variant identifiers;
- normalised parameters and random seed;
- input hashes and transformation history;
- software and external-tool versions;
- quality-control gates and their results;
- warnings, fallbacks and excluded records;
- decision table, uncertainty table and machine-readable provenance;
- `validate()`, `print()`, `summary()`, `plot()` and `as.data.frame()` methods.

### P1 — State the limits of polyploid support

**Evidence.** Some dosage centring accepts alternative ploidy values, but haplotype inference, LD interpretation, phasing and several kernels remain diploid in their assumptions. The package description can be read as broader polyploid support.

**Required correction.** Publish a function-by-function ploidy matrix. Reject unsupported combinations. Add polyploid-specific simulations and reference datasets before claiming operational support.

### P1 — Rebuild performance claims as reproducible evidence

**Evidence.** The documentation makes whole-genome-scale, ten-million-marker and large speed-up claims without a benchmark repository, fixed hardware specification, peak-memory measurement, package-version lock or correctness comparison. Continuous integration contains no large-data performance gate.

**Required correction.**

1. add a versioned benchmark suite with simulated and public reference data;
2. record wall time, central processing unit time, peak resident memory, threads, compiler, hardware and input shape;
3. compare numerical output with a trusted scalar implementation;
4. publish medians and variability from repeated runs;
5. separate measured results from projections;
6. enforce non-regression thresholds on representative medium-scale fixtures.

### P1 — Make vignettes executable specifications

**Evidence.** Three of six vignettes set `eval = FALSE` globally: the large-scale, phasing and end-to-end workflow vignettes. Across all vignettes, 121 chunks are present, but important external and large-scale paths are narratives rather than continuously verified workflows.

**Required correction.**

1. make a small deterministic version of every primary workflow executable;
2. move expensive demonstrations to scheduled integration jobs;
3. show expected result schemas and quality-control interpretation;
4. make every printed table derive from current code;
5. include session information and package version;
6. fail documentation builds on stale code or missing output.

### P1 — Turn optional integrations into explicit test lanes

**Evidence.** The full suite completed with 28 skips and five warnings. Skipped areas included Beagle, epistasis and optional backends. The current workflow installs `covr` but does not run it. No minimum-R, sanitizer, `rchk`, Valgrind, spelling/linting or coverage-threshold job was found.

**Required correction.**

- create separate required jobs for Beagle, `optiSel`, `glmnet`, GDS, PLINK, `bigmemory` and each prediction engine;
- test R 4.2 because it is the declared minimum;
- add AddressSanitizer/UndefinedBehaviorSanitizer and Valgrind or `rchk` coverage for compiled code;
- calculate test coverage and set a decision-critical branch threshold;
- run spelling under British English and a maintained word list;
- pin external actions and remote dependencies to immutable revisions;
- run a scheduled full reverse-dependency and benchmark job.

### P2 — Reduce function size and separate responsibilities

**Evidence.** Several exported functions are difficult to review and test in isolation:

| Function | Approximate lines | Branch tokens |
|---|---:|---:|
| `test_block_haplotypes()` | 887 | 77 |
| `Big_LD()` | 825 | 93 |
| `compute_ld_decay()` | 677 | 68 |
| `select_parents_by_family()` | 642 | 113 |
| `run_ldx_pipeline()` | 593 | 36 |
| `compare_gwas_effects()` | 406 | 52 |
| `select_parents_ga()` | 378 | 46 |
| `usefulness_criterion()` | 329 | 46 |

**Required correction.** Divide input validation, transformation, model fitting, optimisation, diagnostics and presentation into internal units. Use pure functions for objective calculation and postconditions. Keep public wrappers stable while testing each unit directly.

### P2 — Rationalise the public application programming interface

**Evidence.** Naming mixes `Big_LD`, `run_Big_LD_all_chr`, `CLQD`, `ldx_*` and snake case. Both British and American spellings occur: `summarise_blocks()` and `summarize_parent_haplotypes()`, for example. `ldx` names remain after the package identity changed.

**Required correction.**

1. choose snake case for new public functions;
2. choose one spelling policy for prose while preserving stable function names;
3. publish lifecycle badges and a deprecation schedule;
4. retain old names as documented aliases for at least two minor releases;
5. define semantic-versioning and breaking-change policy;
6. audit argument names, result columns and units for consistent terminology.

### P2 — Improve examples and discoverability

**Evidence.** Thirty-one of 85 help pages have no examples. Twenty-four pages use `\dontrun{}` and 28 use `\donttest{}`. Sixteen exports are absent from every vignette, including marker-effect, local-GEBV, stability, visualisation and writer functions.

**Required correction.**

- add a fast executable example for every public function;
- use `\donttest{}` only for genuinely expensive or external operations;
- create task-oriented index pages for data import, haplotypes, prediction, parent selection, mating and reporting;
- cross-link every exported function to one executable workflow;
- provide explicit input/output schemas and units.

### P2 — Repair the README, NEWS and citation path

**Evidence.** [`README.md`](README.md) is approximately 2,034 lines and 112 KB. [`NEWS.md`](NEWS.md) exceeds 3,100 lines and 189 KB, repeats development-version headings and contains a very long current section. The README citation still identifies version 0.3.1 and 2025, and contribution guidance points to a non-existent `tests/testthat/test-core.R`. No `inst/CITATION`, `CITATION.cff`, archived digital object identifier (DOI), `CONTRIBUTING.md`, code of conduct or security policy was found.

**Required correction.**

1. reduce the README to installation, a ten-minute workflow, capability map, limitations and links;
2. move detailed method catalogues to the website;
3. consolidate NEWS by released version and put unreleased items under one heading;
4. generate citations dynamically and add `inst/CITATION` and `CITATION.cff`;
5. archive releases and provide a DOI;
6. add contribution, conduct, security and support policies.

### P2 — Establish British-English spelling quality control

**Evidence.** The package has no `Language` field, so `spelling` defaults to `en-US`. Prose and comments mix programme/program and optimisation/optimization families. [`R/breeder_guide.R`](R/breeder_guide.R) contains “program managers”, “program-shape” and the typo “rather than than”. The package-specific vocabulary also produces a large uncurated spelling report.

**Required correction.**

1. add `Language: en-GB`;
2. create and review `inst/WORDLIST`;
3. correct genuine typographical errors in source, regenerated Rd files and the guide;
4. use British English in prose while retaining established code identifiers;
5. run spelling in continuous integration and fail on new unrecognised prose terms.

### P2 — Correct package-level side claims and metadata

**Evidence.** `DESCRIPTION` is approximately 12.5 KB and its description is a marketing catalogue. It says that the breeder's guide covers seven tools with worked examples; the PDF covers eight tools and contains no worked numerical example. Related R documentation and the breeding-decisions vignette still say seven. The README says eight.

**Required correction.** Make `DESCRIPTION` concise and verifiable. State the eight tools consistently. Replace “worked examples” with an accurate description until the guide contains one. Add automated consistency tests for version, tool count, guide title and citation.

### P2 — Improve interoperability and data governance

**Evidence.** No native Breeding API (BrAPI), Minimum Information About a Plant Phenotyping Experiment (MIAPPE) or Crop Ontology import/export contract was found. These are important for traceable integration with breeding-management systems. BrAPI provides a standard interface for breeding phenotype and genotype databases and aligns with related standards: [BrAPI](https://brapi.org/).

**Required correction.**

- provide canonical mapping tables rather than embedding remote-service logic in modelling functions;
- validate ontology identifiers, units, study/trial/environment identifiers and germplasm identity;
- support export of recommendations and provenance to stable interoperable tables;
- add round-trip tests with representative BrAPI/MIAPPE records.

### P3 — Improve repository hygiene

**Evidence.** Generated graphics appeared in `tests/testthat/Rplots.pdf` during the test run, indicating that at least one plotting test opens the default device. Build artefacts are tracked, and the ignore rules do not cover standard compiled files.

**Required correction.** Use a temporary graphics device in every plotting test; ignore and remove `.o`, `.dll`, temporary plots, check directories and local executables; add a clean-tree assertion after tests.

## Breeder's guide audit

### What works

The 20-page guide is visually clean, readable and consistently styled. It is tagged, contains page numbers and has functional internal table-of-contents links. Its repeated structure—what the tool is, the question it answers, analogy, when to use it, when to skip it and how it works in practice—is appropriate for a non-programming audience. The decision table and glossary are particularly useful. The guide also distinguishes OCS from SimpleMating and acknowledges that the genetic algorithm is stochastic.

### What must improve

1. **Reproducibility.** The source Word document is intentionally excluded, and the PDF is converted manually. The guide cannot be diffed, reviewed or regenerated from the repository. Store the source and build script, or rewrite it in Quarto/R Markdown and render it in continuous integration.
2. **Version control.** The PDF has no package version, release date, compatibility statement, document identifier or change history. Its metadata author is “Un-named”.
3. **Internal consistency.** The guide describes eight tools while package-level documentation says seven. Table-of-contents section 1.7 wraps its page number onto a separate line, and section numbering under Part 2 continues as 1.5 rather than 2.1.
4. **Scientific traceability.** There is no reference list, method citation, source link or explanation of which package versions implement the described behaviour.
5. **Worked decision.** Despite the `DESCRIPTION` claim, the guide contains no worked numerical case. Add one complete example from candidates to recommendation, including an input table, quality-control results, objective components, uncertainty, binding constraints and final decision.
6. **Limits.** Add prominent boxes for phased versus unphased data, diploid versus polyploid support, training-population relevance, G×E limits, minimum sample sizes and the distinction between prediction and causal inference.
7. **Audit use.** Add a one-page sign-off template recording data version, model, seed, constraints, excluded candidates, reviewer and approval.
8. **Language.** Apply British English consistently and correct the source typo before regenerating the PDF.
9. **Accessibility.** Add document bookmarks and meaningful external links, and verify reading order, table headers and alternative text with an accessibility checker.

## Tests and continuous-integration assessment

### Positive evidence

- The test suite is extensive in count and covers many input validations.
- The available paths completed without a failed expectation.
- The continuous-integration matrix includes Linux, macOS and Windows across release, development and old-release R.
- Optional dependencies are commonly guarded, which helps minimal installation.
- Compiled-code checks completed successfully in the audited build.

### Why the current green result is insufficient

The suite completed with 28 skips and five warnings. Some skips represent unavailable external software, but others are triggered by unexpected errors or unsuitable fixtures. More importantly, several assertions check only output type or conditional structure. A solver can fail, `ok` can be `FALSE`, the mating plan can be empty, and the test can still pass. This is the central testing weakness: **the suite often verifies that the wrapper returns something, not that the breeding guarantee was met**.

For decision-critical functions, every test should assert:

1. preconditions and fixture validity;
2. estimator or optimiser convergence;
3. hard-constraint satisfaction;
4. numerical postconditions and invariants;
5. expected behaviour under deliberately infeasible input;
6. a scientifically meaningful result on a truth-known simulation;
7. deterministic reproducibility under a fixed seed;
8. provenance completeness.

## Statistical validation programme

Before a 1.0 release, establish a documented validation battery.

### Simulation truth sets

- founder haplotypes with known block boundaries and recombination;
- causal haplotype, additive, dominance and epistatic architectures;
- population structure, relatives and admixed populations;
- missingness, genotyping error, imputation error and phasing error;
- multi-environment and multi-trait genetic covariance;
- diploid and explicitly supported polyploid scenarios;
- selection across generations, tracking gain, inbreeding and effective population size.

### Required metrics

- block boundary precision/recall and switch error;
- effect bias, root-mean-square error, confidence-interval coverage and false discovery rate;
- prediction correlation, calibration slope, bias, mean squared error and reliability;
- realised versus predicted cross mean and variance;
- constraint-violation rate and optimisation gap on small exact problems;
- genetic gain, diversity loss, inbreeding increment and effective population size across generations;
- wall time, peak memory and numerical agreement across backends and thread counts.

### External validation

Use at least two public crop datasets with contrasting mating systems and population structures. Pre-register the analysis, freeze the data split, compare against transparent baselines and publish all scripts. Do not use the same individuals to discover favourable haplotypes and evaluate the resulting selection advantage.

## Recommended target architecture

```text
Validated input
  ├─ immutable germplasm, sample and variant identity
  ├─ genotype/phase/phenotype/environment quality control
  └─ provenance manifest
          ↓
Typed genetic representation
  ├─ dosage
  ├─ phased haplotype
  ├─ inferred haplotype with uncertainty
  └─ diplotype pattern
          ↓
Model and validation
  ├─ leakage-safe training and tuning
  ├─ family/time/environment-aware validation
  ├─ uncertainty and calibration
  └─ model card
          ↓
Decision engine
  ├─ canonical objective components
  ├─ hard operational constraints
  ├─ feasibility certificate
  └─ sensitivity and trade-off frontier
          ↓
Auditable recommendation
  ├─ parents and mating plan
  ├─ predicted gain and risk
  ├─ diversity/inbreeding trajectory
  └─ signed, versioned decision record
```

## Phased improvement roadmap

### Phase 1 — 0 to 30 days: stop unsafe claims

- fix the genetic-algorithm replicate objective and hard constraints;
- make `optiSel` errors fail safely and strengthen its tests;
- return cross-validation predictions;
- add cache fingerprints and immutable physical identifiers;
- cap threads at two and add interrupts;
- correct multiallelic documentation and backend inconsistency;
- remove tracked compiled artefacts and decide the Beagle redistribution model;
- correct the seven/eight-tool discrepancy, citation, typo and stale CRAN comments;
- label unvalidated integrations and performance claims as experimental.

### Phase 2 — 31 to 90 days: establish release engineering

- introduce the common result/provenance schema;
- create required integration lanes for all optional engines;
- make minimal versions of every vignette executable;
- add coverage, sanitizer, spelling and minimum-R jobs;
- publish the backend capability and ploidy matrices;
- build reproducible benchmarks;
- create `CONTRIBUTING.md`, security policy, citation files and release checklist;
- regenerate the breeder's guide from tracked source.

### Phase 3 — 3 to 6 months: scientific validation

- implement family-, population-, time- and environment-aware validation;
- add uncertainty propagation and reliability gates;
- validate haplotype, association, prediction and optimisation modules on truth-known simulations;
- compare heuristic selections with exact solutions on small problems;
- validate Beagle with truth sets;
- publish two external crop case studies.

### Phase 4 — 6 to 12 months: operational breeding platform

- add multivariate and G×E prediction;
- add economic/desired-gain indices and operational crossing constraints;
- add pedigree-genomic integration and multi-generation diversity monitoring;
- implement BrAPI/MIAPPE/Crop Ontology mappings;
- add model cards, decision sign-off and versioned export bundles;
- release 1.0 only after every decision-critical release gate passes.

## Proposed 1.0 release gates

HapBlockR should not be labelled production- or decision-grade until all of the following are true:

- `R CMD check --as-cran` has no package-caused ERROR, WARNING or NOTE on supported platforms;
- the source tarball is below 10 MB and all redistributed components have verified licensing/source;
- every documented hard constraint is tested as an invariant;
- no unexpected error is converted to a skip;
- all supported external integrations pass required jobs;
- all primary vignettes execute on deterministic fixtures;
- input cache, identity and provenance tests pass across backends;
- prediction validation includes leakage-resistant family/time/environment splits;
- uncertainty and model status accompany every breeder recommendation;
- numerical equality is established across supported backends and one/two-thread execution;
- public performance claims are linked to reproducible benchmarks;
- the breeder's guide is generated from tracked source and matches the current package version;
- the application programming interface, deprecation and semantic-version policies are published.

## Strengths to retain

The package should not be reduced to a narrow collection of unrelated functions. Its principal asset is the integrated path from genomic data to breeding action. The following qualities merit preservation:

- a coherent ambition to connect haplotype analysis with parent and mating decisions;
- broad format support and compiled performance work;
- extensive input checking and many informative warnings;
- a large existing test base and multi-platform continuous integration;
- example data and task-oriented vignettes;
- explicit attention to diversity, relatedness and multiple selection strategies;
- a breeder's guide designed for readers who do not use R;
- candid documentation in several places about approximation, stochastic search and method limits.

The next development cycle should convert those strengths into enforceable contracts. A sophisticated breeding tool is not defined by the number of methods it exposes; it is defined by whether every recommendation is feasible, calibrated, reproducible, scientifically qualified and traceable to the exact data and software that produced it.

## Authoritative external references

- [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
- [Official Beagle software page](https://faculty.washington.edu/browning/beagle/beagle.html)
- [Breeding API](https://brapi.org/)
