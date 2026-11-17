# make_country_adjacency.R — country adjacency registry for the brms_kl
#   spatial term, from the taxdat::afr_sf country polygons
#
# Usage:
#   Rscript analysis/utils/make_country_adjacency.R
#   Rscript analysis/utils/make_country_adjacency.R --tol_km 1 --out <name>
#
# Details:
#   Writes analysis/groupings/<out>.csv (iso3_a,iso3_b; one undirected edge per
#   row, iso3_a < iso3_b). Reads package data only, nothing under data/.
#   Read back with utils/config.R::load_adjacency(); consumed by
#   utils/dlnm_kernel.R::kl_adjacency_data() when
#   occurrence.kernel.spatial.enable. Isolated countries (COM, STP absent from
#   afr_sf; MDG, YEM with no neighbour): the ICAR field pins them near 0, DAGAR
#   gives them N(0, tau).

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
})
args <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--tol_km",
      default = 1,
      type = "double",
      help = "distance tolerance (km) for neighbouring polygons [default 1]"
    ),
    make_option(
      "--out",
      default = "country_adjacency",
      type = "character",
      help = "registry name under analysis/groupings/ [default country_adjacency]"
    )
  )
))

e <- new.env()
utils::data("afr_sf", package = "taxdat", envir = e)
x <- sf::st_make_valid(janitor::clean_names(e$afr_sf))
# Neighbours = polygons within --tol_km. The default 1 km recovers the
# Botswana-Zambia border at Kazungula (~150 m, missed by a strict touch on the
# simplified polygons) and adds no sea crossings (10 km adds none either;
# 50 km starts adding Gulf of Aden / Lake Tanganyika links).
near <- sf::st_is_within_distance(
  x,
  x,
  dist = units::set_units(args$tol_km, "km"),
  sparse = FALSE
)
diag(near) <- FALSE
p <- which(near, arr.ind = TRUE)
p <- p[p[, 1] < p[, 2], , drop = FALSE]
a <- x$country[p[, 1]]
b <- x$country[p[, 2]]
edges <- data.frame(
  iso3_a = pmin(a, b),
  iso3_b = pmax(a, b),
  stringsAsFactors = FALSE
)
edges <- unique(edges[order(edges$iso3_a, edges$iso3_b), ])

path <- file.path("analysis/groupings", paste0(args$out, ".csv"))
writeLines(
  c(
    sprintf(
      "# country adjacency: taxdat::afr_sf polygons within %g km (utils/make_country_adjacency.R)",
      args$tol_km
    ),
    "iso3_a,iso3_b",
    paste(edges$iso3_a, edges$iso3_b, sep = ",")
  ),
  path
)
iso <- sort(unique(x$country))
deg <- table(factor(c(edges$iso3_a, edges$iso3_b), levels = iso))
cat(sprintf(
  "Wrote %s: %d edges over %d countries; isolated: %s\n",
  path,
  nrow(edges),
  length(iso),
  paste(names(deg)[deg == 0], collapse = ", ")
))
