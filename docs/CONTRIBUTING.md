# Contributing to HapBlockR

HapBlockR welcomes focused bug reports, documentation corrections,
tests, and scientifically justified method improvements.

Before opening an issue, please search the existing issues and reduce
the problem to a reproducible example. Include the HapBlockR version, R
version, operating system, relevant optional-package versions, random
seed, and non-sensitive input dimensions. Do not attach confidential
germplasm, genotype, phenotype, pedigree, or trial data.

For a code contribution:

1.  open an issue for substantial behavioural or interface changes;
2.  create a branch from the current development branch;
3.  preserve existing public interfaces unless a deprecation has been
    agreed;
4.  add tests covering the corrected behaviour and its failure modes;
5.  run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
    [`devtools::test()`](https://devtools.r-lib.org/reference/test.html),
    and `devtools::check(args = "--as-cran")`;
6.  use British English in prose while preserving established function
    and argument names; and
7.  describe scientific assumptions, data exclusions, fallbacks, and any
    change to a breeder-facing recommendation.

Tests belong in the most closely related file under `tests/testthat/`;
create a clearly named new file when no suitable file exists. Never
commit external executables, JAR archives, credentials, confidential
data, generated check directories, or compiled objects.

By contributing, you agree that your contribution is distributed under
the package licence and that you will follow
[CODE_OF_CONDUCT.md](https://FAkohoue.github.io/HapBlockR/CODE_OF_CONDUCT.md).
