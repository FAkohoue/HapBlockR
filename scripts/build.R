
# Build favicons from the maintained PNG asset.
pkgdown::build_favicons(overwrite = TRUE)

options(timeout = 3000)  # 5 minutes
pkgdown::build_favicons(overwrite = TRUE)

# Regenerate favicons from the correct logo
pkgdown::build_favicons(overwrite = TRUE)
pkgdown::build_site()
pkgdown::build_site(override = list(template = list(favicon = FALSE)))

library(magick)

# Load your logo
img <- image_read("man/figures/logo.png")

# Create favicon sizes
sizes <- c(16, 32, 48, 64, 180, 192, 512)

dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)

for (s in sizes) {
  resized <- image_resize(img, paste0(s, "x", s))
  image_write(resized, paste0("pkgdown/assets/favicon-", s, ".png"))
}

# Create favicon.ico (multi-size)
ico <- image_resize(img, "64x64")
image_write(ico, "pkgdown/assets/favicon.ico")

writeLines(
  c(
    '<link rel="icon" type="image/svg+xml" href="logo.svg">',
    '<link rel="icon" type="image/png" sizes="32x32" href="favicon-32.png">',
    '<link rel="apple-touch-icon" href="favicon-180.png">'
  ),
  "pkgdown/assets/favicon.html"
)


# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")


#############################################################################

# Historical raster-to-SVG wrapper generation is disabled because it produced
# a 13.5 MB SVG with embedded raster content. The website uses the maintained
# PNG directly.
if (FALSE) {

library(rsvg)
library(magick)
library(base64enc)

png_file <- "man/figures/HapBlockR_schematic.png"
svg_file <- "man/figures/HapBlockR_schematic.svg"

# Maximum Base64 length for each embedded PNG strip.
# This remains safely below the approximately 10 MB XML buffer limit.
max_base64_chars <- 7.5e6

# -------------------------------------------------------------------
# 1. Validate and read the original PNG
# -------------------------------------------------------------------

if (!file.exists(png_file)) {
  stop("Input PNG does not exist: ", png_file)
}

png_file <- normalizePath(
  png_file,
  winslash = "/",
  mustWork = TRUE
)

img  <- magick::image_read(png_file)
info <- magick::image_info(img)[1, ]

width_px  <- as.integer(info$width)
height_px <- as.integer(info$height)

message(
  "Original dimensions: ",
  width_px, " × ", height_px, " pixels"
)

# -------------------------------------------------------------------
# 2. Create temporary working directory
# -------------------------------------------------------------------

tmp_dir <- tempfile("rsvg_png_to_svg_")
dir.create(tmp_dir, recursive = TRUE)

on.exit(
  unlink(tmp_dir, recursive = TRUE, force = TRUE),
  add = TRUE
)

# Estimate the required number of strips from the original PNG size.
estimated_base64_size <- file.info(png_file)$size * 4 / 3

n_strips <- max(
  1L,
  ceiling(estimated_base64_size / max_base64_chars)
)

strip_height <- ceiling(height_px / n_strips)

# -------------------------------------------------------------------
# 3. Split the PNG until every encoded strip is below the XML limit
# -------------------------------------------------------------------

repeat {

  strips <- list()
  strip_too_large <- FALSE

  y_positions <- seq.int(
    from = 0L,
    to   = height_px - 1L,
    by   = strip_height
  )

  for (i in seq_along(y_positions)) {

    y <- y_positions[i]
    h <- min(strip_height, height_px - y)

    tile <- magick::image_crop(
      img,
      geometry = sprintf(
        "%dx%d+0+%d",
        width_px,
        h,
        y
      ),
      repage = TRUE
    )

    tile_file <- file.path(
      tmp_dir,
      sprintf("strip_%04d.png", i)
    )

    magick::image_write(
      tile,
      path   = tile_file,
      format = "png"
    )

    encoded <- base64enc::base64encode(tile_file)

    encoded_size <- nchar(
      encoded,
      type = "bytes"
    )

    if (encoded_size > max_base64_chars) {
      strip_too_large <- TRUE
      break
    }

    strips[[i]] <- list(
      y       = y,
      height  = h,
      base64  = encoded
    )
  }

  if (!strip_too_large) {
    break
  }

  if (strip_height <= 1L) {
    stop(
      "A one-pixel-high strip still exceeds the XML buffer limit."
    )
  }

  strip_height <- max(
    1L,
    floor(strip_height / 2)
  )
}

message(
  "Number of embedded strips: ",
  length(strips)
)

# -------------------------------------------------------------------
# 4. Construct the self-contained intermediate SVG
# -------------------------------------------------------------------

wrapper_file <- file.path(
  tmp_dir,
  "HapBlockR_schematic_wrapper.svg"
)

con <- file(
  wrapper_file,
  open = "wb"
)

on.exit(
  try(close(con), silent = TRUE),
  add = TRUE
)

writeLines(
  sprintf(
    paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<svg xmlns="http://www.w3.org/2000/svg"\n',
      '     width="%d"\n',
      '     height="%d"\n',
      '     viewBox="0 0 %d %d">\n'
    ),
    width_px,
    height_px,
    width_px,
    height_px
  ),
  con = con,
  useBytes = TRUE
)

for (strip in strips) {

  node <- sprintf(
    paste0(
      '  <image x="0" y="%d"\n',
      '         width="%d" height="%d"\n',
      '         preserveAspectRatio="none"\n',
      '         href="data:image/png;base64,%s" />\n'
    ),
    strip$y,
    width_px,
    strip$height,
    strip$base64
  )

  writeLines(
    node,
    con = con,
    useBytes = TRUE
  )
}

writeLines(
  "</svg>",
  con = con,
  useBytes = TRUE
)

close(con)

# -------------------------------------------------------------------
# 5. Use rsvg to generate the final SVG
#
# Do not specify width or height: the native dimensions declared in
# the wrapper are retained.
# -------------------------------------------------------------------

rsvg::rsvg_svg(
  svg  = wrapper_file,
  file = svg_file
)

if (!file.exists(svg_file) || file.info(svg_file)$size == 0) {
  stop("The final SVG was not created correctly.")
}
}



######################################################################################################
#######################################################################################################

# =============================================================================
# HapBlockR - clean build script
# Run SESSION A top-to-bottom, let it restart R, then run SESSION B.
# Both sessions must start with the working directory set to the package root
# (the folder containing DESCRIPTION).
# =============================================================================

# ------------------------------- SESSION A -----------------------------------
# Goal: wipe every stale artifact so the next session starts from zero.

# 1. Remove the installed package so no old DLL can be loaded by mistake.
remove.packages("HapBlockR")

unlink(file.path(.libPaths()[1], "00LOCK-HapBlockR"), recursive = TRUE)

# 2. Delete the compiled DLL from the source tree (src/*.so / src/*.dll).
#    Do this BEFORE the restart so there is nothing to unload conflicts with.
#devtools::clean_dll()

# 3. Restart R.  All loaded DLLs are released, file locks are cleared.
.rs.restartR()
# -- after restart, continue in SESSION B --------------------------------------



# ------------------------------- SESSION B -----------------------------------
# Goal: rebuild everything from source in the correct dependency order.

# 1. Regenerate example data (.rda files in data/ and flat files in
#    inst/extdata/).  Must come first because devtools::document() will
#    try to lazy-load data/ when it parses roxygen @examples.
source("data-raw/generate_example_data.R")

# 2. First document pass: generates NAMESPACE and man/ from roxygen tags in
#    R/*.R.  RcppExports.R does not exist yet so its tags are skipped.
devtools::document()

# 3. Regenerate RcppExports.R (in R/) and src/RcppExports.cpp from the
#    [[Rcpp::export]] annotations in src/ld_core.cpp.
#    This is the canonical source of truth for the C++ symbol table.
Rcpp::compileAttributes()

# 4. Second document pass: picks up the roxygen tags now present in
#    the freshly written R/RcppExports.R.
devtools::document()

# 5. Compile C++ and install into the library.
#    upgrade = FALSE  - do not touch other packages.
devtools::install()

roxygen2::roxygenise()

devtools::load_all()

# 6. Run the test suite.  All C++ symbols are now registered in the
#    installed DLL, so load_all() will find them.
devtools::test()

devtools::test(filter = "association")
devtools::test(filter = "haplotypes")
devtools::test(filter = "full-pipeline")
devtools::test(filter = "breeder-guide")

# 7. Full CRAN check (run after tests pass).
devtools::check()


# 9. Build vignettes
#options(pkgdown.internet = FALSE)
library(HapBlockR)

# Build everything except home, then build home separately
#pkgdown::build_favicons(overwrite = TRUE)

#pkgdown::check_pkgdown()


pkgdown::build_reference()
pkgdown::build_articles()
pkgdown::build_news()

# Build home with network disabled at the curl level
httr2_mock <- function(...) stop("no network", call. = FALSE)
pkgdown::build_home()

#unloadNamespace("OptSLDP")


#pkgdown::clean_site(force = TRUE)
pkgdown::build_site()

# 10. Build package
devtools::build()





while (!is.null(dev.list())) dev.off()

dev.new(width = 7, height = 7)

testthat::test_file(
  "tests/testthat/test-extensions.R"
)

# Delete the cached rendered output for this vignette
unlink("docs/articles/HapBlockR-large-scale.html")
unlink("docs/articles/HapBlockR-large-scale_files", recursive = TRUE)

# Also clear any pkgdown article cache
unlink("vignettes/cache", recursive = TRUE)
unlink("vignettes/HapBlockR-large-scale_cache", recursive = TRUE)

# Now rebuild just that article from scratch
pkgdown::build_article("HapBlockR-large-scale")


file <- "R/haplotype_association.R"

tools::showNonASCIIfile(file)

txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("->", "->", txt, fixed = TRUE)
txt <- gsub("x", "x", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("...", "...", txt, fixed = TRUE)

writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)


txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)

writeLines(txt, file, useBytes = TRUE)

txt <- readLines(file, encoding = "UTF-8")

txt <- gsub("->", "->", txt, fixed = TRUE)
txt <- gsub("x", "x", txt, fixed = TRUE)
txt <- gsub("-", "-", txt, fixed = TRUE)
txt <- gsub("...", "...", txt, fixed = TRUE)
txt <- gsub("union", "union", txt, fixed = TRUE)
txt <- gsub("<=", "<=", txt, fixed = TRUE)
txt <- gsub(">=", ">=", txt, fixed = TRUE)
writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)




txt <- readLines(file, encoding = "UTF-8")

replacements <- c(
  "->"="->",
  "x"="x",
  "-"="-",
  "-"="-",
  "-"="-",
  "..."="...",
  "union"="union",
  "<="="<=",
  ">="=">=",
  "^3"="^3",
  "Lambda"="Lambda"
)

for (sym in names(replacements)) {
  txt <- gsub(sym, replacements[[sym]], txt, fixed = TRUE)
}

writeLines(txt, file, useBytes = TRUE)

tools::showNonASCIIfile(file)
