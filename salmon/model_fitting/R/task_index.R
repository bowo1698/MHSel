# SLURM array index -> (iteration, scenario); 6 scenarios per iteration.
SCENARIOS <- c("S1_equal", "S2_unequal", "S3_equal", "S4_unequal", "S5_equal", "S6_unequal")

task_to_iter_scenario <- function(task, scenarios = SCENARIOS) {
  list(iter = (task - 1) %/% length(scenarios) + 1,
       scenario = scenarios[(task - 1) %% length(scenarios) + 1])
}

# A scenario is done when both models wrote results_all.rds.
should_skip <- function(out_dir, force = FALSE) {
  !force && all(file.exists(file.path(out_dir, c("gblup", "bayesR"), "results_all.rds")))
}
