#### Plots of prediction scores
figure_scores <- function(calib, type, which = c("crps", "dss", "se_mean"), countries = c("IT", "CZ", "FR")){
  ylbls <- c(crps = "RPS", dss = "DSS", se_mean = "SES")
  if (length(countries) > 1){ # there's more than one country in calibration output
    if (type == "case") scores_obj <- rbindlist(lapply(calib,"[[","scores")) else scores_obj <- rbindlist(lapply(calib,"[[","scores_death"))  
  } else {
    if (type == "case") scores_obj <- calib$scores else scores_obj <- calib$scores_death 
  }
  
  y <- sym(which)
  ggplot(scores_obj, aes(x = as.factor(horizon), y = !!y)) + 
    geom_boxplot() +
    labs(x = "Forecast horizon (weeks ahead)", y = ylbls[which]) + 
    facet_wrap(~country, nrow = 1, scales = "free") + 
    theme_cowplot()
}

#### Table of prediction scores for full model and empty model (no covariates or seasonality)
table_scores <- function(all_total, type){
  path <- paste0(if(all_total) "total"  else "age_structured")
  tbl_scores_full <- fread(paste0("Output/", path, "/table_scores_all_", type, ".csv"))
  tbl_scores_full[, model := "Full"]
  tbl_scores_empty <- fread(paste0("Output/", path, "_empty/table_scores_all_", type, ".csv"))
  tbl_scores_empty[, model := "Empty"]
  
  tbl_scores <- rbind(tbl_scores_full, tbl_scores_empty)
  tbl_scores <- dcast(tbl_scores, Country + `Forecast horizon (weeks ahead)` ~ model, value.var = c("RPS", "DSS", "SES"))
  
  write.csv(tbl_scores, paste0("Output/table_scores_", type, "_", if(all_total) "total"  else "age_structured", ".csv"), row.names = F)
}

#### Comparison with ensemble forecasts
figure_country <- function(calib, type, dt_ensemble, min_date, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  loop <- length(countries)
  
  ens <- dt_ensemble[grepl(type, target),]
  
  png(paste0(dir_out, "compar_all_", type, ".png"), width = 900, height = 600)
  par(mfrow = c(4,loop), mar = c(3, 4, 2, 0.5), mar = c(3, 4, 1, 0),
      oma = c(2, 2, 0, 2), las = 1, bty = "l")
  # Plot the weekly country-wide forecast for each forecast horizon and country
  for(i in seq_len(4)){
    for(j in seq_len(loop)){
      ## Extract calib object
      country_j <- countries[j]
      calib_j <- calib[[j]] 

      ## Extract prediction from calib
      if(type == "case") calib_week <- calib_j$calib else if(type == "death") calib_week <- calib_j$calib_death
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
      
      obs_i <- obs_i[names(obs_i) > min_date]
      up_i <- up_i[names(pred_i) > min_date]
      low_i <- low_i[names(pred_i) > min_date]
      pred_i <- pred_i[names(pred_i) > min_date]
      ens_i <- ens_i[rownames(ens_i) > min_date,]

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
      max_y <- max(abs((abs(pred_i - obs_i) - abs(ens_i[,2] - obs_i))))
      
      diff_pred_i <- abs(pred_i - obs_i)
      diff_ens_i <- abs(ens_i[,2] - obs_i)
      
      if(type == "case"){
        if(country_j == "FR" || country_j == "IT") thresh <- 5000
        if(country_j == "CZ") thresh <- 1000
      } else if(type == "death"){ 
        if(country_j == "FR" || country_j == "IT") thresh <- 50
        if(country_j == "CZ") thresh <- 0
      }
    }
  }
  title(xlab = "Prediction date (weeks)", line = 0, outer = TRUE, cex.lab = 1.5)
  title(ylab = paste0("Number of ", type, "s"), line = 1, outer = TRUE, cex.lab = 1.5)
  dev.off()
}

#### Comparison of scores with those of ensemble forecasts
figure_compar_scores <- function(scores_obj, which, ylbl){
  y <- sym(which)
  ggplot(scores_obj, aes(x = target_end_date, y = !!y, group = model, color = model)) +
    geom_line() +
    facet_grid(location ~ horizon, scales = "free") +
    labs(x = "Date", y = ylbl, color = "Model") +
    theme_cowplot() #+
    # theme(strip.background = element_blank())
}

#### Comparison of scores with those of ensemble forecasts
table_compar_scores <- function(all_total, empty, type){
  path <- paste0(if(all_total) "total"  else "age_structured", if(empty) "_empty", "", "/")
  wis <- fread(paste0("Output/", path, "compar_scores_all_", type, "_quantile.csv"))
  wis <- dcast(wis, Country + `Forecast horizon (weeks)` ~ model, value.var = "interval_score")
  
  se_median <- fread(paste0("Output/", path, "compar_scores_all_", type, "_point.csv"))
  se_median <- dcast(se_median, Country + `Forecast horizon (weeks)` ~ model, value.var = "se_point")
  
  tbl_scores <- merge(wis, se_median, by = c("Country","Forecast horizon (weeks)"), suffixes = c(" WIS"," SE median"))
  write.csv(tbl_scores, paste0("Output/", path, "compar_scores_", type, ".csv"), row.names = F)
}

#### Comparison data and median predictions 
figure_median <- function(calib, type, min_date, countries = c("IT", "CZ", "FR"), prop = T){
  ggplot_list <- list()
  ## Generate a heatmap for each country and horizon forecast
  for(i in seq_len(3)){
    mat_prop <- data.frame()
    for(j in 1:4){
      # Extract death or case data and forecasts
      if(type == "death"){
        obs_vec <- calib[[i]]$compar_deaths[[j]][1,]
        med_pred_vec <- calib[[i]]$compar_deaths[[j]][2,]
      } else {
        obs_vec <- calib[[i]]$compar_cases[[j]][1,]
        med_pred_vec <- calib[[i]]$compar_cases[[j]][2,]
      }
      # Create groups of observations (depending on the country and type)
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
        if(all(obs_vec >= 1)){
          groups <- c(-1, 5, 10, 20, 40, 100, 200, 5e6, 1e7)
          groups_lab <- c("0-5", "5-10", "10-20", "20-40", "40-100", "100-200", "200+")
        }
      }
      # Exclude groups with no observation
      groups_lab <- groups_lab[-which(groups > max(obs_vec))]
      groups <- groups[-which(groups > max(obs_vec))[-1]]
      med_pred_group <- cut(med_pred_vec, groups, groups_lab)
      obs_group <- cut(obs_vec, groups, groups_lab)
      
      # Compute the number of observation in each data / forecast group
      nb_group <- table(med_pred_group, obs_group)
      prop_group <- t(nb_group)/as.numeric(table(obs_group))
      melt_prop <- melt(prop_group)
      melt_prop$nb <- melt(t(nb_group))$value
      if(countries[i] == "IT"){
        melt_prop$week <- paste(countries[i], j, ifelse(j > 1, "weeks ahead,", "week ahead,"), 
                                "per region", sep = " ")
      } else melt_prop$week <- paste(countries[i], j, ifelse(j > 1, "weeks ahead,", "week ahead,"), 
                                     "per region and age group", sep = " ")
      mat_prop <- rbind(mat_prop, melt_prop)
    }
    if(!prop){
      mat_prop$value <- mat_prop$nb
      lim <- c(0, max(mat_prop$value))
      lab <- "Number of \nobservations"
    } else{
      lim <- c(0, 1)
      lab <- "Proportion"
    }
    # Create the heat map per country
    gg_prop <- ggplot(mat_prop, aes(obs_group, med_pred_group)) + geom_tile(aes(fill = value)) + 
      facet_wrap(week ~ ., nrow = 4, scales = "free") +
      theme_bw() + 
      labs(fill = lab) + 
      scale_fill_gradient(low = "white", high = "#132B43", limits = lim)
    ggplot_list[[i]] <- gg_prop
  }
  # Merge the heatmap for each country in the same grid
  plot <- plot_grid(
    ggplot_list[[1]] + theme(legend.position = "none", axis.title = element_blank()), 
    ggplot_list[[2]] + theme(legend.position = "none", axis.title = element_blank()),
    ggplot_list[[3]] + theme(legend.position = "none", axis.title = element_blank()),
    get_legend(ggplot_list[[1]]), nrow = 1, rel_widths = c(3,3,3,1))  
  # Add axis titles
  y.grob <- textGrob(paste0("Predicted number of ", type, "s"), 
                     gp=gpar(fontface="bold", fontsize=15), rot=90)
  
  x.grob <- textGrob("Data", 
                     gp=gpar(fontface="bold", fontsize=15))
  
  #add to plot
  
  grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
  
  
}

#### PIT histogram 
figure_pit <- function(calib, type, min_date, countries = c("IT", "CZ", "FR"), dir_out = "Output/"){
  length_country <- 3
  png(paste0(dir_out, "pit_all_", type, ".png"), width = 600, height = 400)
  par(mfrow = c(4,3), mar = c(3, 4, 1, 0), oma = c(2, 2, 0, 2), las = 1, bty = "l")
  ## If the type is "case", then daily forecasts, if "death": weekly forecasts
  if(type == "case") nb_days <- 7 else nb_days <- 1
  for(j in seq_len(4)){
    ## Extract the pit objects from calib
    for(i in seq_len(length_country)){
      if(country == "all"){
        if(type == "case") list_p <- calib[[i]]$pit else list_p <- calib[[i]]$pit_death
        country_i <- countries[i]
      }
      which_dates <- which(rownames(calib[[1]]$calib$pred_week) >= min_date)
      
      ## Extract the proportion of simulations below x and xm1
      px <- list_p[["px"]][,,which_dates]
      pxm1 <- list_p[["pxm1"]][,,which_dates]
      
      ## Function to generate Probability Integral Transform
      Fbar1 <- function (u, Px, Pxm1){
        F_u <- punif(u, Pxm1, Px)
        mean(F_u)
      }
      ## Define the number of bins in PIT histogram
      J <- 10
      breaks <- (0:J)/J
      ## Extract values of px and pxm1 in the jth week
      px_j <- c(px[(j-1) * nb_days + seq_len(nb_days),,])
      pxm1_j <- c(pxm1[(j-1) * nb_days + seq_len(nb_days),,])
      ## Remove NA values
      pxm1_j <- pxm1_j[!is.na(px_j)]
      px_j <- px_j[!is.na(px_j)]
      if(type == "death") ymax <- 2 else ymax <- 3
      Fbar_seq_prop <- vapply(X = breaks, FUN = Fbar1, FUN.VALUE = 0, Px = px_j,Pxm1 = pxm1_j, USE.NAMES = FALSE)
      breaks <- breaks - .05
      ## Plot PIT histogram
      barplot((diff(Fbar_seq_prop) * J) ~ breaks[-1], ylim = c(0, ymax), xlab = "", ylab = "",
              border = NA, col = "black", space = 0, xlim = c(0, 9.5))
      label <- paste(country_i, j, ifelse(j > 1, "weeks ahead", "week ahead"), sep = " ")
      title(main = label, line = -1)
      abline(h = 1, lty = 2)
      abline(h = c(1.5, .5), lty = 2, col = "red")
    }
    title(xlab = "Probability Integral Transform", line = 0, outer = TRUE, cex.lab = 1.5)
    title(ylab = "Relative Frequency", line = 0, outer = TRUE, cex.lab = 1.5)
  }
  dev.off()
}
