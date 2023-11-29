## Table with prediction scores
scores_table <- function(calib, country, type, dir_out = "Output/"){
  if (country == "all"){
    if (type == "case") scores_obj <- rbindlist(lapply(calib,"[[","scores")) else scores_obj <- rbindlist(lapply(calib,"[[","scores_death"))
  } else {
    if (type == "case") scores_obj <- calib$scores else scores_obj <- calib$scores_death
  }
  
  quantile_na <- function(x, probs, ...) quantile(x, probs = probs, na.rm = T,  ...)
  cols <- c("crps", "dss", "se_mean")
  quant_scores <- scores_obj[, unlist(
    lapply(.SD, function(x) list(q95l = quantile_na(x, probs = 0.025),
                                 q50l = quantile_na(x, probs = 0.25),
                                 med = quantile_na(x, probs = 0.5),
                                 q50u = quantile_na(x, probs = 0.75),
                                 q95u = quantile_na(x, probs = 0.975))),
    recursive = F), .SDcols = cols, by = .(country, horizon)
  ]
  tbl_scores <- quant_scores[,.(
    Country = country,
    `Forecast horizon (weeks ahead)` = horizon,
    RPS = med_and_CI(crps.med, crps.q95l, crps.q95u, d = 3, method = "signif"),
    DSS = med_and_CI(dss.med, dss.q95l, dss.q95u, d = 3, method = "signif"),
    SES = med_and_CI(se_mean.med, se_mean.q95l, se_mean.q95u, d = 3, method = "signif")
  )]
  write.csv(tbl_scores, paste0(dir_out, "table_scores_", country, "_", type,".csv"), row.names = F)
  quant_scores_reg <- scores_obj[, unlist(
    lapply(.SD, function(x) list(q95l = quantile_na(x, probs = 0.025),
                                 q50l = quantile_na(x, probs = 0.25),
                                 med = quantile_na(x, probs = 0.5),
                                 q50u = quantile_na(x, probs = 0.75),
                                 q95u = quantile_na(x, probs = 0.975))),
    recursive = F), .SDcols = cols, by = .(country, horizon, reg_nb)
  ]
  tbl_scores_reg <- quant_scores_reg[,.(
    Country = country,
    Region = reg_nb,
    `Forecast horizon (weeks ahead)` = horizon,
    RPS = med_and_CI(crps.med, crps.q95l, crps.q95u, d = 3, method = "signif"),
    DSS = med_and_CI(dss.med, dss.q95l, dss.q95u, d = 3, method = "signif"),
    SES = med_and_CI(se_mean.med, se_mean.q95l, se_mean.q95u, d = 3, method = "signif")
  )]
  write.csv(tbl_scores_reg, paste0(dir_out, "table_scores_reg_", country, "_", type, ".csv"), row.names = F)
  quant_scores_age <- scores_obj[, unlist(
    lapply(.SD, function(x) list(q95l = quantile_na(x, probs = 0.025),
                                 q50l = quantile_na(x, probs = 0.25),
                                 med = quantile_na(x, probs = 0.5),
                                 q50u = quantile_na(x, probs = 0.75),
                                 q95u = quantile_na(x, probs = 0.975))),
    recursive = F), .SDcols = cols, by = .(country, horizon, age)
  ]
  tbl_scores_age <- quant_scores_age[,.(
    Country = country,
    Age = age,
    `Forecast horizon (weeks ahead)` = horizon,
    RPS = med_and_CI(crps.med, crps.q95l, crps.q95u, d = 3, method = "signif"),
    DSS = med_and_CI(dss.med, dss.q95l, dss.q95u, d = 3, method = "signif"),
    SES = med_and_CI(se_mean.med, se_mean.q95l, se_mean.q95u, d = 3, method = "signif")
  )]
  write.csv(tbl_scores_age, paste0(dir_out, "table_scores_age_", country, "_", type, ".csv"), row.names = F)
}

## Plot prediction scores
scores_figure <- function(calib, country, type, which = c("crps", "dss", "se_mean"), total = F, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  if(country == "all") loop <- 3 else loop <- 1
  
  ylbls <- c(crps = "RPS", dss = "DSS", se_mean = "SES")
  y <- sym(which)
  
  for(j in seq_len(loop)){
    if (country == "all"){
      country_j <- countries[j]
      calib_j <- calib[[j]]
    } else {
      country_j <- country
      calib_j <- calib
    }
    
    if (type == "case") scores_obj <- calib_j$scores else scores_obj <- calib_j$scores_death
    
    if (total){
      scores_obj <- scores_obj[, lapply(.SD, mean), .SDcols = c("crps", "dss", "se_mean"), by = .(country, horizon,reg_nb)]  
    }
    
    p <- ggplot(scores_obj, aes(x = as.factor(horizon), y = !!y)) + 
      geom_boxplot() +
      labs(x = "Forecast horizon (weeks ahead)", y = ylbls[which]) + 
      theme_cowplot()
    if(!total){
      p <- p + facet_wrap(~age, scales = "free")
    }
    png(paste0(dir_out, "scores_", country_j, "_", type, "_", 
               if(total) "total" else "age",".png"), width = 1200, height = 900)
    print(p)
    dev.off()
  }
}

## Generate PIT histograms
pit_figure <- function(calib, country, type, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  # Check if the plots are generated for one country or all of them
  if(country == "all") loop <- 3 else loop <- 1
  
  for(i in seq_len(loop)){
    ## Extract the pit objects from calib
    if(country == "all"){
      if(type == "case") list_p <- calib[[i]]$pit else list_p <- calib[[i]]$pit_death
      country_i <- countries[i]
    } else{
      if(type == "case") list_p <- calib$pit else list_p <- calib$pit_death
      country_i <- country
    }
    
    ## Extract the proportion of simulations below x and xm1
    px <- list_p[["px"]]
    pxm1 <- list_p[["pxm1"]]
    
    ## Function to generate Probability Integral Transform
    Fbar1 <- function (u, Px, Pxm1){
      F_u <- punif(u, Pxm1, Px)
      mean(F_u)
    }
    ## Define the number of bins in PIT histogram
    J <- 10
    breaks <- (0:J)/J
    max_iter <- nrow(px)
    ## If the type is "case", then daily forecasts, if "death": weekly forecasts
    if(type == "case") nb_days <- 7 else nb_days <- 1
    png(paste0(dir_out, "pit_", country_i, "_", type, ".png"), width = 600, height = 450)
    par(mfrow = c(2,2), mar = c(3, 4, 2, 0.5), mar = c(3, 4, 1, 0),
        oma = c(2, 2, 0, 2), las = 1, bty = "l")
    for(j in seq_len(max_iter / nb_days)){
      ## Extract values of px and pxm1 in the jth week
      px_j <- c(px[(j-1) * nb_days + seq_len(nb_days),,])
      pxm1_j <- c(pxm1[(j-1) * nb_days + seq_len(nb_days),,])
      ## Remove NA values
      pxm1_j <- pxm1_j[!is.na(px_j)]
      px_j <- px_j[!is.na(px_j)]
      Fbar_seq_prop <- vapply(X = breaks, FUN = Fbar1, FUN.VALUE = 0, Px = px_j,Pxm1 = pxm1_j, USE.NAMES = FALSE)
      ## Plot PIT histogram
      plot((diff(Fbar_seq_prop) * J) ~ breaks[-1], type = "h", lwd = 10, lend = "butt", ylim = c(0,4), xlab = "", 
           ylab = "")
      label <- paste(type, "PIT histogram\n", country_i, j, ifelse(j > 1, "weeks ahead", "week ahead"), sep = " ")
      title(main = label, line = -2)
      abline(h = 1, lty = 2)
      abline(h = c(1.5, .5), lty = 2, col = "red")
    }
    title(xlab = "Probability Integral Transform", line = 0, outer = TRUE, cex.lab = 1.5)
    title(ylab = "Relative Frequency", line = 0, outer = TRUE, cex.lab = 1.5)
    dev.off()
    ## Generate PIT by broad age group + week
    if(country_i == "CZ" || country_i == "FR"){
      png(paste0(dir_out, "pit_per_age_", country_i, "_", type, ".png"), width = 900, height = 600)
      par(mfrow = c(4,4), mar = c(3, 4, 2, 0.5), mar = c(3, 4, 1, 0),
          oma = c(2, 2, 0, 2), las = 1, bty = "l")
      ## Generate PIT by broad age groups
      broad_groups <- list(c("0-9", "10-19"), c("20-29", "30-39", "40-49", "50-59"),
                           c("60-69", "70-79"), "80+")
      lab_broad <- c("0-19", "20-59", "60-79", "80+")
      for(i in seq_len(nrow(px) / nb_days)){
        for(j in seq_along(broad_groups)){
          age_groups <- broad_groups[[j]] 
          cols <- which(is.element(sub(".*[.]", "", colnames(px)), age_groups))
          px_i <- c(px[(i-1) * nb_days + seq_len(nb_days), cols,])
          pxm1_i <- c(pxm1[(i-1) * nb_days + seq_len(nb_days), cols,])
          pxm1_i <- pxm1_i[!is.na(px_i)]
          px_i <- px_i[!is.na(px_i)]
          breaks <- (0:J)/J
          Fbar_seq_prop <- vapply(X = breaks, FUN = Fbar1, FUN.VALUE = 0, Px = px_i,Pxm1 = pxm1_i,USE.NAMES = FALSE)
          breaks <- breaks - .05
          barplot((diff(Fbar_seq_prop) * J) ~ breaks[-1], ylim = c(0,4), xlab = "", ylab = "",
                  border = NA, col = "black", space = 0, xlim = c(0, 9.5))
          abline(h = 1, lty = 2)
          abline(h = c(1.5, .5), lty = 2, col = "red")
          label <- paste(lab_broad[j], "years old\n", country_i, i, 
                         ifelse(i > 1, "weeks ahead", "week ahead"), 
                         sep = " ")
          title(main = label, line = -3)
        }
      }
      title(xlab = "Probability Integral Transform", line = 0, outer = TRUE, cex.lab = 1.5)
      title(ylab = "Relative Frequency", line = 0, outer = TRUE, cex.lab = 1.5)
      dev.off()
    }
  }
}

## Plot comparison with ensemble forecasts
compar_figure <- function(calib, country, type, dt_ensemble, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  if(country == "all") loop <- 3 else loop <- 1
  
  ens <- dt_ensemble[grepl(type, target),]
  
  for(j in seq_len(loop)){
    ## Extract calib object
    if(country == "all"){
      country_j <- countries[j]
      calib_j <- calib[[j]] 
    } else{
      country_j <- country
      calib_j <- calib
    }
    
    png(paste0(dir_out, "compar_", country_j, "_", type, ".png"), width = 600, height = 450)
    par(mfrow = c(4,1), mar = c(3, 4, 2, 0.5), mar = c(3, 4, 1, 0),
        oma = c(2, 2, 0, 2), las = 1, bty = "l")
    ## Extract prediction from calib
    if(type == "case") calib_week <- calib_j$calib else if(type == "death") calib_week <- calib_j$calib_death
    for(i in seq_len(4)){
      ## Extract observed data and predictions from our model
      obs_i <- calib_week[[1]][i:length(calib_week[[1]])]
      pred_i <- calib_week[[2]][,i][seq_along(obs_i)]
      up_i <- calib_week[["up"]][,i][seq_along(obs_i)]
      low_i <- calib_week[["low"]][,i][seq_along(obs_i)]
      names(obs_i) <- names(pred_i)
      ## Extract predictions from the ensemble model
      ens_i <- dcast(ens[location == country_j & quantile %in% c(0.025, 0.5, 0.975) & 
                           (as.character(as.Date(forecast_date) - 2) %in% names(pred_i)) &
                           grepl(i, target),], forecast_date ~ quantile)
      ens_i <- as.matrix(ens_i[, !"forecast_date"], rownames.value = ens_i[, as.character(as.Date(forecast_date) - 2)])
      
      if(!all(diff(as.Date(rownames(ens_i))) == 7)){
        missing_dates <- names(pred_i)[!is.element(names(pred_i), rownames(ens_i))]
        missing_matrix <- matrix(NA, nrow = length(missing_dates), ncol = ncol(ens_i))
        rownames(missing_matrix) <- missing_dates
        ens_i <- rbind(ens_i, missing_matrix)
        ens_i <- ens_i[order(as.Date(rownames(ens_i))),]
      }
      
      max_y <- max(c(obs_i, pred_i, ens_i), na.rm = T)
      max_y <- min(max_y, max(obs_i, na.rm = T) * 4)
      
      min_y <- 0
      log_axis <- ""
      if(type == "case"){
        min_y <- min(c(obs_i, pred_i, ens_i), na.rm = T)
        ens_i[ens_i == 0] <- 1
        obs_i[obs_i == 0] <- 1
        pred_i[pred_i == 0] <- 1
        if(min_y == 0) min_y <- 1
        if(country_j == "IT") min_y <- 1e3
        if(country_j == "CZ") min_y <- 1e2
        if(country_j == "FR") min_y <- 1e3
        log_axis <- "y"
      }
      ## Plot observed data
      date_plot <- as.Date(names(pred_i))
      plot(obs_i ~ date_plot, type = "l", ylim = c(min_y, max_y), 
           xlab = "", ylab = "", lwd = 2, log = log_axis)
      ## Add 95% prediction interval for the ensemble model, removing points where
      ## no forecast is available
      if(any(is.na(ens_i[,1]))){
        index_i <- seq_len(which(is.na(ens_i[,1]))[1] -1)
        polygon(x = c(date_plot[index_i], rev(date_plot[index_i])), 
                y = c(ens_i[index_i, 1], rev(ens_i[index_i, 3])), 
                col = transp("#66c2a5", .4), border = NA)
        for(k in which(is.na(ens_i[, 1]))){
          if(any(is.na(ens_i[-seq_len(k), 1]))){
            index_i <- seq(k + 1, which(is.na(ens_i[,1]) & seq_along(ens_i[,1])>k)[1] - 1)
          } else index_i <- seq(k + 1, nrow(ens_i))
          if(all(diff(index_i) > 0) | length(index_i) == 1)
            polygon(x = c(date_plot[index_i], rev(date_plot[index_i])), 
                    y = c(ens_i[index_i, 1], rev(ens_i[index_i, 3])), 
                    col = transp("#66c2a5", .4), border = NA)
        }
      } else 
        polygon(x = c(date_plot, rev(date_plot)), 
                y = c(ens_i[, 1], rev(ens_i[, 3])), 
                col = transp("#66c2a5", .4), border = NA)
      ## Plot median ensemble forecasts
      lines(ens_i[, 2] ~ date_plot, col = "#66c2a5", lwd = 1, lty = 2)
      ## Add 95% prediction interval of our model
      polygon(x = c(date_plot, rev(date_plot)), y = c(up_i, rev(low_i)), 
              col = transp("#fc8d62", .4), border = NA)
      ## Plot our predictions
      lines(pred_i ~ date_plot, col = "#fc8d62", lwd = 1, lty = 2)
      
      label <- paste(country_j, i, ifelse(i > 1, "weeks ahead", "week ahead"), sep = " ")
      title(main = label, line = 0)
      
      if(i == 1) legend("topright", legend = c("ensemble", "EE model", "data"), 
                        col = c("#66c2a5", "#fc8d62", "black"), lwd = 1, bty = "n")
      # max_y <- max(abs((abs(pred_i - obs_i) - abs(ens_i[,2] - obs_i))))
      # 
      # diff_pred_i <- abs(pred_i - obs_i)
      # diff_ens_i <- abs(ens_i[,2] - obs_i)
      # 
      # if(type == "case"){
      #   if(country_j == "FR" || country_j == "IT") thresh <- 5000
      #   if(country_j == "CZ") thresh <- 1000
      # } else if(type == "death"){ 
      #   if(country_j == "FR" || country_j == "IT") thresh <- 50
      #   if(country_j == "CZ") thresh <- 0
      # }
      # 
      # print(c(
      #   sum((diff_pred_i - diff_ens_i) < -thresh, na.rm = T), # nb when pred_i is much better
      #   sum((diff_pred_i - diff_ens_i) > thresh, na.rm = T)) # nb when ens_i is much better
      # )
    }
    title(xlab = "Time", line = 0, outer = TRUE, cex.lab = 1.5)
    title(ylab = paste0("Number of ", type, "s"), line = 1, outer = TRUE, cex.lab = 1.5)
    dev.off()
  }
}

compar_scores <- function(calib, country, type, quants, dt_ensemble, countries = c("IT", "CZ", "FR")){
  if(country == "all") loop <- 3 else loop <- 1
  list_EE <- vector("list", loop)
  list_obs <- vector("list", loop)
  # Loop over countries to extract EE model median and quantile forecasts
  for(j in seq_len(loop)){
    ## Extract calib object
    if(country == "all"){
      country_j <- countries[j]
      calib_j <- calib[[j]]
    } else{
      country_j <- country
      calib_j <- calib
    }
    
    if(type == "case") calib_week <- calib_j$calib else calib_week <- calib_j$calib_death
    
    obs_j <- as.data.table(calib_week$obs, keep.rownames = T)
    names(obs_j) <- c("date", "true_value")
    obs_j[, date := as.Date(date)]
    
    dt_EE_j <- calib_week$q_pred[round(quantile,3) %in% round(quants,3)]
    setnames(dt_EE_j, "date", "forecast_date")
    dt_EE_j[, horizon := as.integer(horizon)]
    dt_EE_j[, target_end_date := forecast_date + 7 * horizon]
    dt_EE_j[, location := country_j]
    dt_EE_j[, date := target_end_date - 7]
    dt_EE_j <- merge(dt_EE_j, obs_j, by = "date")
    dt_EE_j[, date := NULL]
    dt_EE_j[, model := "EE"]
    
    list_EE[[j]] <- dt_EE_j
    list_obs[[j]] <- obs_j
  }
  
  dt_EE <- rbindlist(list_EE)
  obs <- rbindlist(list_obs, idcol = "id")
  obs[,location := countries[id]]
  obs[,id := NULL]
  
  # Process European COVID-19 Forecast Hub ensemble forecasts
  dt_ensemble1 <- copy(dt_ensemble)
  # Select countries
  dt_ensemble1 <- dt_ensemble1[location %in% countries]
  # Select case or death forecasts based on type
  dt_ensemble1 <- dt_ensemble1[grep(type, target)]
  # Convert quantile column to numeric
  dt_ensemble1[, quantile := as.numeric(quantile)]
  # Select quantiles to use when calculating score
  dt_ensemble1 <- dt_ensemble1[quantile %in% quants]
  setnames(dt_ensemble1, "value", "prediction")
  dt_ensemble1[, horizon := as.integer(sub("^([0-9]+).*", "\\1", target))]
  dt_ensemble1[, forecast_date := as.Date(forecast_date)]
  dt_ensemble1[, target_end_date := as.Date(target_end_date)] 
  # dt_ensemble1[, date := target_end_date + 2 - 7]
  dt_ensemble1[, date := target_end_date - 7]
  dt_ensemble1 <- merge(obs, dt_ensemble1, by = c("location","date"))
  dt_ensemble1[, `:=`(date = NULL, target = NULL)]
  dt_ensemble1[, model := "ensemble"]
  
  # Bind EE and ensemble forecasts
  dt_EE_ensemble <- rbind(dt_EE, dt_ensemble1)
  # Remove dates without forecasts for both models
  dt_EE_ensemble <- dt_EE_ensemble[
    target_end_date %in% 
      intersect(dt_EE[, target_end_date], dt_ensemble1[, target_end_date])]
  # Calculate weighted interval score for quantile forecasts
  scores_quantile <- score(dt_EE_ensemble, metrics = "interval_score")
  setorder(scores_quantile, model, location, horizon)
  
  dt_EE_ensemble_point <- dt_EE_ensemble[quantile == 0.5]
  dt_EE_ensemble_point[, quantile := NA]
  # Calculate scores for point (median) forecasts
  scores_point <- score(dt_EE_ensemble_point)
  scores_point[, which(sapply(scores_point, function(x) all(is.na(x)))) := NULL]
  setorder(scores_point, model, location, horizon)
  
  return(list(scores_quantile = scores_quantile, scores_point = scores_point))
}

compar_scores_table <- function(scores_obj, country, type, forecast_type, dir_out = "Output/"){
  nms <- names(scores_obj)
  cols <- nms[!(nms %in% c("forecast_date", "horizon", "target_end_date", "location", "model"))]
  
  tbl_scores <- scores_obj[, 
    lapply(.SD, function(x) 
      med_and_CI(
        quantile(x, probs = 0.5),
        quantile(x, probs = 0.025),
        quantile(x, probs = 0.975),
        d = 3,
        method = "signif"
      )), .SDcols = cols, by = .(model, location, horizon)
  ]
  
  setnames(tbl_scores, c("location", "horizon"), c("Country", "Forecast horizon (weeks)"))
  
  write.csv(tbl_scores, paste0(dir_out, "compar_scores_", country, "_", type, "_", forecast_type, ".csv"), row.names = F)
}

compar_scores_figure <- function(scores_obj, country, type, which, ylbl, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  if(country == "all") loop <- 3 else loop <- 1
  y <- sym(which)
  
  for(i in seq_len(loop)){
    country_i <- countries[i]
    scores_obj_i <- scores_obj[location == country_i,]
    
    png(filename = 
          paste0(dir_out, "compar_scores_", country_i, "_", type, "_", 
                 if(which == "interval_score") "quantile" else "point", ".png"), 
        width = 1200, height = 300)
    print(ggplot(scores_obj_i, aes(x = target_end_date,y = !!y,group = model, color = model)) +
      geom_line() + 
      labs(x = "Date", y = ylbl, color = "Model") +
      theme_cowplot() + 
      facet_wrap(~horizon, nrow = 1))
    dev.off()
  }
}

## Comparison data and median predictions
compar_median <- function(calib, country, type, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  # Check if the plots are generated for one country or all of them
  if(country == "all") loop <- 3 else loop <- 1
  
  ## Extract the pit objects from calib
  for(i in seq_len(loop)){
    if(country == "all"){
      if(type == "case") compar <- calib[[i]]$compar_cases else compar <- calib[[i]]$compar_deaths
      country_i <- countries[i]
    } else{
      if(type == "case") compar <- calib$compar_cases else compar <- calib$compar_deaths
      country_i <- country
    }
    
    mat_prop <- data.frame()
    for(j in 1:4){
      if(type == "death"){
        obs_vec <- compar[[j]][1,]
        med_pred_vec <- compar[[j]][2,]
      } else {
        obs_vec <- compar[[j]][1,]
        med_pred_vec <- compar[[j]][2,]
      }
      if(type == "case"){
        if(countries[i] == "IT"){
          groups <- c(-1, 50, 100, 200, 400, 1000, 2000, 5e6, 1e7)
          groups_lab <- c("0-50", "50-100", "100-200", "200-400", "400-1000", "1000-2000", "2000+")
        } else if(countries[i] == "CZ") { 
          groups <- c(-1, 10, 25, 50, 75, 100, 150, 250, 500, 1000, 1e7)
          groups_lab <- c("0-10", "10-25", "25-50", "50-75", "75-100", "100-150", "150-250", "500-1000", "1000+")
        } else if(countries[i] == "FR"){
          groups <- c(-1, 25, 50, 100, 250, 500, 1000, 5e6, 1e7)
          groups_lab <- c("0-25", "25-50", "50-100", "100-250", "250-500", "500-1000", "1000+")
        }
      } else {
        groups <- c(-1, 1, 5, 10, 20, 40, 100, 200, 5e6, 1e7)
        groups_lab <- c("0-1", "1-5", "5-10", "10-20", "20-40", "40-100", "100-200", "200+")
      }
      groups_lab <- groups_lab[-which(groups > max(obs_vec))]
      groups <- groups[-which(groups > max(obs_vec))[-1]]
      med_pred_group <- cut(med_pred_vec, groups, groups_lab)
      obs_group <- cut(obs_vec, groups, groups_lab)
      
      nb_group <- table(med_pred_group, obs_group)
      prop_group <- t(nb_group)/as.numeric(table(obs_group))
      
      melt_prop <- melt(prop_group)
      melt_prop$nb <- melt(t(nb_group))$value
      melt_prop$week <- paste(type, country_i, j, ifelse(j > 1, "weeks ahead", "week ahead"), sep = " ")
      mat_prop <- rbind(mat_prop, melt_prop)
    }
    
    png(filename = paste0(dir_out, country_i, "_", type, "heatmap_prop.png"), width = 600, height = 600)
    gg_prop <- ggplot(mat_prop, aes(obs_group, med_pred_group)) + geom_tile(aes(fill = value)) + 
      facet_wrap(~week) + theme_bw() + labs(x = "Data", y = "Predicted number of cases", 
                                            fill = "Proportion") + 
      scale_fill_gradient(low = "white", high = "#132B43")
    plot(gg_prop)
    dev.off()
    png(filename = paste0(dir_out, country_i, "_", type, "heatmap_nb.png"), width = 600, height = 600)
    gg_nb <- ggplot(mat_prop, aes(obs_group, med_pred_group)) + geom_tile(aes(fill = nb)) + 
      facet_wrap(~week) + theme_bw() + labs(x = "Data", y = "Predicted number of cases", 
                                            fill = "Number") + 
      scale_fill_gradient(low = "white", high = "#132B43")
    plot(gg_nb)
    dev.off()
    
  }
}
