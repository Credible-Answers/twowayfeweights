#' Internal helpers (renaming, transformation, filtering, summary).
#'
#' @importFrom data.table as.data.table is.data.table setnames set setDT
#' @importFrom stats as.formula coef vcov
#' @noRd

printf <- function(...) print(sprintf(...))

fn_ctrl_rename            <- function(x) paste("ctrl",   x, sep = "_")
fn_treatment_rename       <- function(x) paste("OT",     x, sep = "_")
fn_treatment_weight_rename <- function(x) paste0("weight_", x)
fn_random_weight_rename   <- function(x) paste("RW",     x, sep = "_")

get_controls_rename       <- function(controls)   unlist(lapply(controls,   fn_ctrl_rename))
get_treatments_rename     <- function(treatments) unlist(lapply(treatments, fn_treatment_rename))
get_random_weight_rename  <- function(ws)         unlist(lapply(ws,         fn_random_weight_rename))


# -----------------------------------------------------------------------------
# Rename input columns to internal canonical names.
# -----------------------------------------------------------------------------
twowayfeweights_rename_var <- function (df, Y, G, T, D, D0, controls, treatments, random_weights){
  
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  
  # Initialisation
  Y_df              <- NULL
  G_df              <- NULL
  T_df              <- NULL
  D_df              <- NULL

  D0_df             <- NULL
  random_weight_df  <- NULL
  controls_df       <- NULL
  treatments_df     <- NULL

  controls_rename       <- NULL
  treatments_rename     <- NULL
  random_weight_rename  <- NULL

  # To solve the devtools::check() NOTE
  ..cols <- NULL
  ..id.col <- NULL

  # Generate a random character value that
  # is not contained in the column names
  # to act as the identification column
  id.col = paste(sample(LETTERS, 10, TRUE), collapse = '')
  while(id.col %in% colnames(df)){
    id.col = paste(sample(LETTERS, 10, TRUE), collapse = '')
  }
  df[, (id.col) := 1:nrow(df)]

  # Create one dataframe for the 1-length vector variables
  core.variables <- c(Y, G, T, D, D0)
  names.variables <- c("Y", "G", "T", "D", "D0")
  for(x in seq_along(core.variables)){
    if(!is.null(core.variables[x])){
      name.df <- paste0(names.variables[x], "_df")
      cols <- c(core.variables[x], id.col)
      assign(name.df, df[, ..cols])

      new.names <- c(names.variables[x], id.col)
      data.table::setnames(eval(as.symbol(name.df)), new.names)
    }
  }

  if(!is.null(controls)){
    controls_rename     <- get_controls_rename(controls)
    cols <- c(controls, (id.col))
    controls_df = df[, ..cols]

    new.names.controls  = c(controls_rename, id.col)
    data.table::setnames(controls_df, new.names.controls)
  }

  if(!is.null(treatments)){
    treatments_rename   <- get_treatments_rename(treatments)
    cols <- c(treatments, (id.col))
    treatments_df = df[, ..cols]

    new.names.treatments  = c(treatments_rename, id.col)
    data.table::setnames(treatments_df, (new.names.treatments))
  }

  if(length(random_weights) > 0) {
    random_weight_rename <- get_random_weight_rename(random_weights)
    cols <- c(random_weights, (id.col))
    random_weight_df = df[, ..cols]

    new.names.random_weights   = c(random_weight_rename, id.col)
    data.table::setnames(random_weight_df, (new.names.random_weights))
  }

  out <- df[, ..id.col]
  final.variables <- c("Y", "G", "T", "D", "D0")
  for(x in final.variables){
    df.to.merge <- eval(as.symbol(paste0(x, "_df")))
    if(!is.null(df.to.merge)){
      out <- out[df.to.merge, on = (id.col)]
    }
  }

  # Merge the other variables
  # complementary.variables <- c(controls_rename, treatments_rename, random_weight_rename)
  complementary.variables <- c("controls_df", "treatments_df", "random_weight_df")
  for(x in seq_along(complementary.variables)){
    if(!is.null(complementary.variables[x])){
      df.to.merge <- eval(as.symbol(paste0(complementary.variables[x])))
      if(!is.null(df.to.merge)){
        out <- out[df.to.merge, on = (id.col)]
      }
    }
  }

  # Remove the id column
  out[, (id.col) := NULL][]

  out
}


# -----------------------------------------------------------------------------
# Normalize within-cell variation, attach weights, build T factors.
# -----------------------------------------------------------------------------
twowayfeweights_transform <- function(df, controls, weights, treatments) {

  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  ret <- twowayfeweights_normalize_var(df, "D")
  if (ret$retcode) {
    df <- ret$df
    printf("The treatment variable in the regression varies within some group * period cells.")
    printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
    printf("with a group * period level treatment.")
    printf("The command will replace the treatment by its average value in each group * period.")
    printf("The results below apply to the two-way fixed effects regression with that treatment variable.")
  }

  for (control in controls) {
    ret <- twowayfeweights_normalize_var(df, control)
    if (ret$retcode) {
      df <- ret$df
      printf("The control variable %s in the regression varies within some group * period cells.", control)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with controls apply to group * period level controls.")
      printf("The command will replace replace control variable %s by its average value in each group * period.", control)
      printf("The results below apply to the regression with control variable %s averaged at the group * period level.", control)
    }
  }

  for (treatment in treatments) {
    ret <- twowayfeweights_normalize_var(df, treatment)
    if (ret$retcode) {
      df <- ret$df
      printf("The other treatment variable %s in the regression varies within some group * period cells.", treatment)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with several treatments apply to group * period level controls.")
      printf("The command will replace replace other treatment variable %s by its average value in each group * period.", treatment)
      printf("The results below apply to the regression with other treatment variable %s averaged at the group * period level.", treatment)
    }
  }

  if (is.null(weights)) {
    data.table::set(df, j = "weights", value = 1)
  } else {
    data.table::set(df, j = "weights", value = weights)
  }

  df[, Tfactor := factor(T)]
  df[, TFactorNum := as.numeric(Tfactor)]

  df
}


# -----------------------------------------------------------------------------
# Drop rows with missing values according to the chosen estimation type.
# Note: at this point in the pipeline `df` already has the canonical column
# names "Y", "G", "T", "D" (and "D0" if set), so the original Y/G/T/D/D0
# parameter values are unused -- we keep them in the signature for backwards
# compatibility with the call site.
# -----------------------------------------------------------------------------
twowayfeweights_filter <- function(df, Y, G, T, D, D0, cmd_type, controls,
                                   treatments) {

  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  na_count <- function(dt, cols) {
    if (length(cols) == 0) return(rep(0L, nrow(dt)))
    Reduce(`+`, lapply(cols, function(cc) as.integer(is.na(dt[[cc]]))))
  }

  if (cmd_type != "fdTR") {
    cols <- c("Y", "G", "T", "D", controls, treatments)
    df   <- df[na_count(df, cols) == 0]
  } else {
    tag1 <- na_count(df, c("D", "T", "Y"))
    tag2 <- na_count(df, "D0")
    keep <- (tag1 == 0) | (tag2 == 0)
    df   <- df[keep]
    tag1 <- tag1[keep]

    if (length(controls) > 0) {
      tag3 <- na_count(df, controls)
      df   <- df[(tag1 == 1) | (tag3 == 0)]
    }
  }

  df
}


# -----------------------------------------------------------------------------
# Tally positive / negative weights.
# -----------------------------------------------------------------------------
twowayfeweights_summarize_weights <- function(df, var_weight) {

  w <- df[[var_weight]]
  ok <- !is.na(w)

  weight_plus  <- w[ok & w > 0]
  weight_minus <- w[ok & w < 0]

  list(
    nr_plus    = length(weight_plus),
    nr_minus   = length(weight_minus),
    nr_weights = length(weight_plus) + length(weight_minus),
    sum_plus   = sum(weight_plus),
    sum_minus  = sum(weight_minus)
  )
}


# -----------------------------------------------------------------------------
# Test correlation between candidate variables and the weights.
# -----------------------------------------------------------------------------
twowayfeweights_test_random_weights <- function(df, random_weights) {

  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  mat <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mat) <- c("Coef", "SE", "t-stat", "Correlation")

  df_filtered <- df[is.finite(W) & nat_weight != 0]

  for (v in random_weights) {
    fml <- stats::as.formula(sprintf("%s ~ W", v))
    rw_lm <- fixest::feols(fml = fml, data = df_filtered,
                           weights = ~nat_weight, vcov = ~G)
    beta <- stats::coef(rw_lm)[["W"]]
    se   <- sqrt(diag(stats::vcov(rw_lm)))[["W"]]
    r2   <- fixest::r2(rw_lm)["r2"]
    mat[v, ] <- c(beta, se, beta / se,
                  if (beta > 0) sqrt(r2) else -sqrt(r2))
  }

  mat
}
