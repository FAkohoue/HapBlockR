## demo_pipeline.R
## ---------------------------------------------------------------------------
## An end-to-end demonstration of DesiredGainR on the example maize programme.
##
## The sequence follows the order in which the decisions actually have to be
## made: state the objective, test whether it is attainable, build the index,
## evaluate it, and only then rank candidates.
##
## Run with: source(system.file("examples", "demo_pipeline.R",
##                              package = "DesiredGainR"))
## ---------------------------------------------------------------------------

library(DesiredGainR)

data(dgr_traits)
data(dgr_G)
data(dgr_P)
data(dgr_candidates)
data(dgr_gebv)
data(dgr_history)

traits <- dgr_traits$trait
lower_is_better <- dgr_traits$trait[dgr_traits$direction == "decrease"]

message("Traits: ", paste(traits, collapse = ", "))
message("To be decreased: ", paste(lower_is_better, collapse = ", "))

## -- 1. Inspect the covariance matrices before trusting anything -------------
## Index coefficients come from inverting these matrices, so their conditioning
## determines whether the coefficients mean anything.

print(matrix_diagnostics(dgr_G, "G")[
  c("condition_number", "positive_definite")
])
print(matrix_diagnostics(dgr_P, "P")[
  c("condition_number", "positive_definite")
])

## Most of that conditioning comes from the trait scales rather than from the
## correlations. Standardising removes it.
message(
  "Condition number of the genetic correlation matrix: ",
  format(matrix_diagnostics(stats::cov2cor(dgr_G))$condition_number,
         digits = 4)
)

## -- 2. Translate between the two ways of stating an objective ---------------
## A breeder who cannot state economic weights can state desired gains, and
## each implies the other exactly.

## Desired gains are stated as improvements, and lower_is_better tells the
## function which traits improve by falling. Without it, G and P would be read
## in the raw trait direction and a positive gain for plant height would be a
## request to grow taller.

desired_gains <- c(GY = 1.0, PHT = 0.4, AD = 0.6, ASI = 0.5,
                   EPP = 0.4, GLS = 0.6)

implied <- implied_economic_weights(
  desired_gains, dgr_G, dgr_P,
  lower_is_better = lower_is_better, gain_units = "genetic_sd"
)
print(round(implied, 4))

## Two of these weights are negative, for traits the objective asks to improve.
## That is correct rather than an error, and it does not mean those traits
## should decrease: lower_is_better has already oriented everything, so larger
## is better throughout. Once oriented, grain yield correlates +0.55 with the
## anthesis-silking interval and +0.45 with ears per plant, so selecting for a
## full standard deviation of yield would carry both past the smaller gains
## requested for them. To deliver the stated ratio the index must hold them
## back, and their weights turn negative.

## The magnitudes mislead more than the signs do. A weight carries inverse
## trait units, so a trait on a small scale receives a large number for the
## same emphasis. Rescaling by the genetic standard deviation makes them
## comparable, and the picture changes: the largest number in the vector turns
## out to be one of the smaller effects.

genetic_sd <- stats::setNames(dgr_traits$genetic_sd, dgr_traits$trait)
print(round(
  sort(implied * genetic_sd[names(implied)], decreasing = TRUE), 3
))

## The translation is exactly invertible when the units match, which is the
## check worth making before trusting either direction.
print(round(
  implied_desired_gains(
    implied, dgr_G, dgr_P,
    lower_is_better = lower_is_better, gain_units = "genetic_sd"
  ),
  4
))

## -- 3. Ask whether the objective is attainable at all -----------------------
## Grain yield and anthesis date are positively correlated, yet the objective
## requires yield up and anthesis date down, so not every target is reachable.

feasibility <- gain_feasibility(
  desired_gains = desired_gains,
  G = dgr_G, P = dgr_P,
  n_candidates = nrow(dgr_candidates),
  n_select = 20L,
  lower_is_better = lower_is_better,
  gain_units = "genetic_sd"
)
print(feasibility)

## When the target is out of reach, the attainable response is the requested
## direction rescaled to the intensity actually applied. Read it in the units
## the target was stated in, alongside the original trait units.
print(round(feasibility$attainable_response_input_units, 3))
print(round(feasibility$attainable_response, 3))

## -- 4. Recover what the programme has been selecting for --------------------
## Where an index has never been used, the historical selection differentials
## already encode the objective.

recovered <- retrospective_weights(
  selected_values = dgr_candidates[dgr_history$selected, traits],
  population_values = dgr_candidates[, traits],
  trait_cols = traits
)
print(recovered)
print(round(recovered$selection_differential_sd, 3))

## -- 5. Build and compare index families -------------------------------------
## Traits are standardised, following Covarrubias-Pazaran (2021), so that a
## desired gain of one means one standard deviation of progress.

economic_weights <- c(GY = 1.0, PHT = 0.2, AD = 0.5, ASI = 0.4,
                      EPP = 0.3, GLS = 0.5)

smith_hazel <- selection_index(
  dgr_candidates, traits, method = "smith_hazel",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
base_index <- selection_index(
  dgr_candidates, traits, method = "base",
  G = dgr_G, P = dgr_P, economic_weights = economic_weights,
  lower_is_better = lower_is_better, n_select = 20L, main_trait = "GY"
)
rank_sum <- selection_index(
  dgr_candidates, traits, method = "mulamba_mock",
  lower_is_better = lower_is_better, n_select = 20L
)

print(smith_hazel)
print(base_index)

## Guimaraes et al. (2021) found the rank-sum index competitive despite
## requiring neither weights nor covariance matrices, so it is worth comparing.
## Report the rank correlation alongside the overlap in the selected set: the
## two must be consistent with each other, and a high correlation paired with a
## near-empty intersection would indicate a fault rather than a finding.

comparison <- merge(
  smith_hazel$ranking[, c("id", "score")],
  rank_sum$ranking[, c("id", "score")],
  by = "id", suffixes = c("_sh", "_mm")
)
message(
  "Spearman correlation between the two index scores: ",
  round(stats::cor(comparison$score_sh, comparison$score_mm,
                   method = "spearman"), 3)
)
message(
  "Candidates selected by both Smith-Hazel and the rank-sum index: ",
  length(intersect(smith_hazel$selected$id, rank_sum$selected$id)), " of 20"
)

## -- 6. Check which traits the index is really acting on ---------------------
## An index coefficient is not comparable across traits on different scales.

print(smith_hazel$effective_weights)

## -- 7. Ask how much the decision depends on the weights ---------------------
## Weights are estimates. If the selected set survives perturbing them, further
## refinement will not change the decision.

sensitivity <- weight_sensitivity(
  economic_weights, dgr_candidates, dgr_G, dgr_P,
  n_select = 20L, trait_cols = traits, n_draws = 100L
)
print(sensitivity)

## -- 8. Fit the optimised desired-gain index ---------------------------------
## run_dgsi() searches for a desired-gain vector whose realised response in the
## selected set approaches the stated target.

dgsi <- run_dgsi(
  init_data = dgr_candidates[, c("GenoID", "Family")],
  cand_data = dgr_candidates,
  trait_cols = traits,
  dg = desired_gains,
  G = dgr_G,
  P = dgr_P,
  lower_is_better = lower_is_better,
  scale_traits = TRUE,
  n_select = 20L,
  n_iter = 200L,
  n_rep = 5L,
  seed = 42L
)

print(round(dgsi$realised_response, 3))
print(dgsi$replicate_diagnostics)

## Read this against step 3. The feasibility check reported that only part of
## the requested gain is attainable *in the requested proportion*. The
## optimisation is not bound to that proportion, so it distributes the
## shortfall unevenly: it overshoots grain yield and undershoots anthesis date,
## rather than delivering a uniformly scaled version of the target. Whether
## that trade is acceptable is a breeding decision, not a numerical one.

print(round(dgsi$realised_response - desired_gains, 3))

## The classical solution, without iteration, is returned as a comparator.
## Note the anthesis-silking interval: the non-iterated index moves it in the
## wrong direction, which is precisely the failure Joukhadar et al. (2024)
## reported for the unoptimised method.
print(round(dgsi$non_iterated$realised_response, 3))

## -- 9. Score genomic estimated breeding values with the quadratic index -----

quadratic_weights <- diag(c(
  GY = 0.05, PHT = -0.03, AD = -0.03, ASI = -0.04, EPP = 0.02, GLS = -0.03
))
dimnames(quadratic_weights) <- list(traits, traits)

## scale_traits = TRUE matters more here than anywhere else in this script.
## Economic weights multiply the genomic estimated breeding values directly, so
## on the original scale plant height, whose standard deviation is a hundred
## times that of ears per plant, would dominate the index whatever weight it
## was given, and grain yield would go backwards. Standardising puts the stated
## weights in charge of the objective.

qgsi <- run_qgsi(
  init_data = dgr_gebv["GenoID"],
  gebv_data = dgr_gebv,
  trait_cols = traits,
  linear_weights = economic_weights,
  W = quadratic_weights,
  lower_is_better = lower_is_better,
  scale_traits = TRUE,
  n_select = 20L
)

print(qgsi)
print(qgsi$expected_gain_per_trait)

## The same call on the original scale warns, and the warning is worth reading
## rather than suppressing.
unscaled <- tryCatch(
  run_qgsi(
    init_data = dgr_gebv["GenoID"],
    gebv_data = dgr_gebv,
    trait_cols = traits,
    linear_weights = economic_weights,
    W = quadratic_weights,
    lower_is_better = lower_is_better,
    scale_traits = FALSE,
    n_select = 20L
  ),
  warning = function(w) {
    message("Expected warning: ", conditionMessage(w))
    suppressWarnings(run_qgsi(
      init_data = dgr_gebv["GenoID"],
      gebv_data = dgr_gebv,
      trait_cols = traits,
      linear_weights = economic_weights,
      W = quadratic_weights,
      lower_is_better = lower_is_better,
      scale_traits = FALSE,
      n_select = 20L
    ))
  }
)
message(
  "Expected grain yield gain, standardised: ",
  round(qgsi$expected_gain_per_trait$Expected_Genetic_Gain[1L], 4),
  "; unstandardised: ",
  round(unscaled$expected_gain_per_trait$Expected_Genetic_Gain[1L], 4)
)

message("Demonstration complete.")
