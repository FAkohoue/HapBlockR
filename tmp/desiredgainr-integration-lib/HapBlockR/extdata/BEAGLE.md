# Beagle is an external dependency

HapBlockR does not distribute a Beagle JAR archive.

The supported integration targets Beagle 5.x. Download Beagle and its
corresponding source from the official Beagle website, review its GNU General
Public License terms, and provide the local JAR path through one of these
explicit mechanisms:

1. `beagle_jar = "/absolute/path/to/beagle.jar"`;
2. `options(HapBlockR.beagle_jar = "/absolute/path/to/beagle.jar")`;
3. the `HAPBLOCKR_BEAGLE_JAR` environment variable; or
4. `beagle.jar` in the directory containing `out_prefix`.

Official download and licence information:
<https://faculty.washington.edu/browning/beagle/beagle.html>

The required integration job currently tests
`beagle.27Feb25.75f.jar` (Beagle 5.5) and verifies SHA-256
`7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70`
before execution. This checksum identifies the tested file; it is not an
endorsement or a substitute for reviewing the upstream licence and source.

For a truth-known validation set, pass `truth_vcf` to
`phase_with_beagle()` and set the switch-error, dosage-accuracy,
allele-concordance, and call-rate thresholds appropriate to the study.
`assess_phasing_accuracy()` can also compare two phased VCFs directly.
