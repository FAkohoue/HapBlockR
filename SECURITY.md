# Security and responsible disclosure

## Supported versions

Security and data-integrity fixes are applied to the current development
version and the latest released minor version where practical.

## Reporting

Do not open a public issue for a vulnerability that could expose data,
execute untrusted content, corrupt scientific results, or compromise an
external tool invocation. Report it privately to `akohoue.f@gmail.com` with:

- the affected HapBlockR and R versions;
- a minimal reproduction using synthetic, non-sensitive data;
- the expected and observed behaviour;
- the likely impact; and
- any proposed mitigation.

Do not include credentials, private repository links, confidential
germplasm identifiers, or real genotype and phenotype records.

HapBlockR invokes external tools only through paths and arguments supplied or
explicitly configured by the user. Users remain responsible for obtaining
those tools from authoritative sources, verifying licences and checksums, and
protecting the input and output directories.
