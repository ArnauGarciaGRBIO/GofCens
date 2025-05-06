print.cumhazPlot <- function(x, ...) {
  if (x$prnt && x$outp == "list") {
    cat("Parameter estimates\n")
    cat(" ", "\n")
    for (dist in x$distributions) {
      if (!is.null(x$params[[dist]])) {
        cat(dist, "\n", sep = "")
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(x$params[[dist]])) {
            cat("Location (se): ", round(x$params[[dist]][1], x$degs), " " ,
                "(", round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
          }
          if ("shape" %in% names(x$params[[dist]])) {
            cat("   Shape (se): ", round(x$params[[dist]][1], x$degs), " " ,
                "(",round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
          }
          cat("   Scale (se): ", round(x$params[[dist]][2], x$degs), " " ,
              "(", round(x$se[[dist]][2], x$degs), ")", "\n", sep = "")
        } else {
          if (dist == "exponential") {
            cat("   Scale (se): ", round(x$params[[dist]], x$degs), " " ,
                "(", round(x$se[[dist]], x$degs), ")", "\n", sep = "")
          } else {
            cat("  Shape1 (se): ", round(x$params[[dist]]$parameters[1], x$degs), " " ,
                "(", round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
            cat("  Shape2 (se): ", round(x$params[[dist]]$parameters[2], x$degs), " " ,
                "(", round(x$se[[dist]][2], x$degs), ")", "\n", sep = "")
            cat("  Domain:", round(x$params[[dist]]$domain[1], x$degs), "-",
                round(x$params[[dist]]$domain[2],x$degs), "\n")
          }
        }
        if(x$print.AIC){
          cat( "   AIC:", round(x$aic[[dist]], x$degs), "\n")
        }
        if(x$print.BIC){
          cat( "   BIC:", round(x$bic[[dist]], x$degs), "\n")
        }
        cat("\n")
      }
    }
  } else if (x$prnt && x$outp == "table"){

    column_list <- list()
    count <- 1
    for (dist in x$distributions) {
      if (!is.null(x$params[[dist]])) {
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic", "beta")) {

          param1 <- paste0(round(x$params[[dist]][1], x$degs), " " ,
                        "(",round(x$se[[dist]][1], x$degs), ")", sep = "")
          param2 <- paste0(round(x$params[[dist]][2], x$degs), " " ,
                        "(",round(x$se[[dist]][2], x$degs), ")", sep = "")
          if(x$print.AIC){
            aic_prt <- round(x$aic[[dist]], x$degs)
          } else {
            aic_prt <- NA
          }
          if(x$print.BIC){
            bic_prt <- round(x$bic[[dist]], x$degs)
          } else {
            bic_prt <- NA
          }
          column_list[[count]] <- c(dist, param1, param2, aic_prt, bic_prt)

        } else {

          param1 <- paste0(round(x$params[[dist]], x$degs), " " ,
                        "(",round(x$se[[dist]], x$degs), ")", sep = "")
          if(x$print.AIC){
            aic_prt <- round(x$aic[[dist]], x$degs)
          } else {
            aic_prt <- NA
          }
          if(x$print.BIC){
            bic_prt <- round(x$bic[[dist]], x$degs)
          } else {
            bic_prt <- NA
          }
          column_list[[count]] <- c(dist, param1, NA, aic_prt, bic_prt)
        }
        count <- count + 1
      }
    }

    df_prt <- as.data.frame(column_list)
    names(df_prt) <- as.character(df_prt[1, ])
    df_prt <- df_prt[-1, ]
    rownames(df_prt) <- c("Param 1 (s.e.)", "Param 2 (s.e.)", "AIC", "BIC")
    if(any(is.na(df_prt["AIC",]))){
      df_prt <- df_prt[rownames(df_prt) != "AIC", ]
    }
    if(any(is.na(df_prt["BIC",]))){
      df_prt <- df_prt[rownames(df_prt) != "BIC", ]
    }

    print(df_prt)

  }
  if ("beta" %in% x$distributions) {
    x$params$beta <- x$params$beta$parameters
  }
  invisible(lapply(x$params, round, x$degs))
}
