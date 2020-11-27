

require(assertthat)
require(tidyverse)
require(dplyr)
require(reshape2)

VERSION <- '0.2.2'






mix_wells <- function(assay.df, sample.df){
  ## Mix wells from sample and assay together ##
  
  # check the column headers of the assay
  assay_headers <- c("assay_well", "assay_enzyme", 
                     "assay_enzyme_conc", "assay_enzyme_unit", "repeat_num")
  
  if (!are_equal(colnames(assay.df), assay_headers)){
    stop(sprintf('The column headers in the assay input file must be: %s', paste(assay_headers, collapse=', ')))
  }
  
  # check the column headers of the sample
  sample_headers <- c("sample_well", "sample_substrate", "sample_substrate_conc", 
                      "sample_substrate_unit")

  if (!are_equal(colnames(sample.df), sample_headers)){
    stop(sprintf('The column headers in the sample input file must be: %s', paste(assay_headers, collapse=', ')))
  }
    
  # make all combinations of wells
  mixed.df <- expand.grid(sample.df$sample_well, assay.df$assay_well)
  
  # fix the colnames
  colnames(mixed.df) <- c('sample_well', 'assay_well')
  
  # add in the original data
  mixed.df <- merge(mixed.df, sample.df, by='sample_well')
  mixed.df <- merge(mixed.df, assay.df, by='assay_well')
  
  
  # calculate final concentrations of things
  mixed.df$final_well <- paste(mixed.df$sample_well, mixed.df$assay_well, sep='-')
  
  mixed.df$final_substrate_name <- mixed.df$sample_substrate_name
  mixed.df$final_substrate_conc <- mixed.df$sample_substrate_conc
  mixed.df$final_substrate_unit <- mixed.df$sample_substrate_unit
  
  mixed.df$final_enzyme_name <- mixed.df$assay_enzyme_name
  mixed.df$final_enzyme_conc <- mixed.df$assay_enzyme_conc
  mixed.df$final_enzyme_unit <- mixed.df$assay_enzyme_unit
  
  return(mixed.df)
}



combine_data <- function(mixed.df, rox.df, fam.df){
  ## Combine data from the different sources ##
  
  # convert to long format
  rox.df.long <- melt(rox.df, id=c('Chamber.ID')) 
  colnames(rox.df.long) <- c('final_well', 'cycle', 'rox_value')
  rox.df.long$cycle <- as.integer(gsub('X', '', rox.df.long$cycle))
  
  fam.df.long <- melt(fam.df, id=c('Chamber.ID')) 
  colnames(fam.df.long) <- c('final_well', 'cycle', 'fam_value')
  fam.df.long$cycle <- as.integer(gsub('X', '', fam.df.long$cycle))
  
  # Add the data into my frame with concentration info
  all.df <- merge(fam.df.long, rox.df.long, by=c('final_well', 'cycle'))
  all.df <- merge(all.df, mixed.df, by=c('final_well'))
  
  # Change units of the protein input to mg/ml, if nessecary
  if (unique(all.df$assay_enzyme_unit) == 'mg/ml') {
    # everything is fine, do nothing
    placeholder <- 0
    
  } else if (unique(all.df$assay_enzyme_unit) == 'ug/ml') {
    all.df$final_enzyme_conc <- all.df$final_enzyme_conc / 1000
    all.df$final_enzyme_unit <- 'mg/ml'
 
  } else if (unique(all.df$assay_enzyme_unit) == 'ng/ml') {
    all.df$final_enzyme_conc <- all.df$final_enzyme_conc / 1000000
    all.df$final_enzyme_unit <- 'mg/ml'
       
  } else {
    stop(print('The protein concentration in the assay input file must be either in ng/ml, ug/ml or mg/ml'))
    
  }
  
  
  return(all.df)
}





# ----- Define a function for plotting a matrix ----- #
# from http://www.phaget4.org/R/image_matrix.html
myImagePlot <- function(x, min, max, ...){
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  num_cols = 256
  ColorRamp <- colorRampPalette(c('#222222', '#EEEEEE'))(n = num_cols)
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Add a separate color for removed samples, by adding another "step" to the color palette
  zstep <- (max - min) / num_cols # step in the color palette
  newz.na <- max + zstep # the value for na will be one step above the current max
  max <- max + zstep # extend top limit to include the two new values above and na
  na.color = 'red'
  
  # find NA and prelace these with the new value
  x[which(is.na(x))] <- newz.na # same for newz.na
  
  # add the NA color to the color ramp
  col <- c(ColorRamp, na.color) # we construct the new color range by including: na.color and na.outside
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=col, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  ColorLevels <- seq(min, max, length=length(col))
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
        col=col,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #




plot_raw_data_heatmaps <- function(df, min=0, max=65000){
  ## Make heatmaps of the raw data ##
  
  for (cycle in c(1, max(df$cycle))){
    # make a subset of the data frames for plotting
    rox_data_first_cycle.df <- df[df$cycle == cycle, which(names(df) %in% c('assay_well', 'sample_well', 'rox_value'))]
    fam_data_first_cycle.df <- df[df$cycle == cycle, which(names(df) %in% c('assay_well', 'sample_well', 'fam_value'))]
    
    # make wide format
    rox_data_wide.df <- dcast(rox_data_first_cycle.df, assay_well ~ sample_well, value.var='rox_value')
    fam_data_wide.df <- dcast(fam_data_first_cycle.df, assay_well ~ sample_well, value.var='fam_value')
    
    # use first column as row names
    rox.mat_data <- data.matrix(rox_data_wide.df[, 2:length(names(rox_data_wide.df))]) 
    rnames <- rox_data_wide.df[,'assay_well']   # assign labels in the "File" column to "rnames"
    rownames(rox.mat_data) <- rnames 
    
    fam.mat_data <- data.matrix(fam_data_wide.df[, 2:length(names(fam_data_wide.df))]) 
    rnames <- fam_data_wide.df[,'assay_well']   # assign labels in the "File" column to "rnames"
    rownames(fam.mat_data) <- rnames 
    
    
    # Now make the heatmap
    myImagePlot(rox.mat_data, min=0, max=68000, title=(sprintf('ROX channel, cycle %s', cycle)))
    myImagePlot(fam.mat_data, min=0, max=68000, title=(sprintf('FAM channel, cycle %s', cycle)))
  }
}



subtract_no_enz_fluorescence <- function(df){
  ## subtract fluorescence of no enzyme control  ##
  
  # first get the values from control
  control.df <- df[df$final_enzyme_conc==0, c('cycle', 'assay_enzyme', 'sample_substrate', 'sample_substrate_conc', 'repeat_num', 'rox_value')]
  
  # rename the rox value
  colnames(control.df)[colnames(control.df) == 'rox_value'] <- 'rox_value_control'
  
  # do a left join to add on the controls to the samples
  df.combined <- merge(df, control.df, by=(c('cycle', 'assay_enzyme', 'sample_substrate', 'sample_substrate_conc', 'repeat_num')), all.x=TRUE)
  
  # now calculate the normalized values
  df.combined$rox_value_norm <- df.combined$rox_value - df.combined$rox_value_control
  
  # remove the control samples, they are not needed anymore as the data is now normalized
  df.combined <- df.combined[!df.combined$final_enzyme_conc==0, ]
  df.combined$final_substrate_norm_unit <- unique(df.combined$final_substrate_unit)
  
  # drop unneeded columns
  small.df <- df.combined[, -which(names(df.combined) %in% c('rox_value_control'))]
  
  # re-order the data
  small.df <- small.df %>% 
    select(assay_well, sample_well, assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, minutes, repeat_num, rox_value_norm, final_substrate_norm_unit)
  
  return(small.df)
}



trim_data_low_range <- function(df, tolerance_at_zero){
  # remove samples where no activity is seen, within a certain tolerance #

  # what's the maximum value observed in the measurements
  max_val = max(concentration.df$resorufin_value_norm)

  # find which enzyme-substrate concentrations give no perceptibel activity (within a tolerance level)
  no_activity.df <- concentration.df %>%
    group_by(assay_well, sample_well) %>%
    summarize(assay_enzyme = list(assay_enzyme),
              final_enzyme_conc = list(final_enzyme_conc),
              final_enzyme_unit = list(final_enzyme_unit),
              sample_substrate = list(sample_substrate),
              final_substrate_conc = list(final_substrate_conc),
              final_substrate_unit = list(final_substrate_unit),
              minutes = list(minutes),
              repeat_num = list(repeat_num),
              final_substrate_norm_unit = list(final_substrate_norm_unit),
              resorufin_value_norm = list(resorufin_value_norm),
              resorufin_value_norm_unit = list(resorufin_value_norm_unit)
    ) 
  
  # now actually filter out the bad samples
  without_no_activity.df <- no_activity.df[sapply(no_activity.df$resorufin_value_norm, max) > 0 + tolerance_at_zero * max_val, ]
  
  # unbundle
  unbundled.df <- unnest(without_no_activity.df)
  
  return(unbundled.df)
}


fluorescence_to_concentration <- function(df, chip_type){
  # convert measured fluorescence to product concentration using a standard curve
  
  # y = kx + m
  # x = (y-m)/k
  
  if (chip_type == 'flex_six'){
    # y = 1384x + 2614
    k = 1247
    m = 199
    
  } else if (chip_type == '96_96') {
    stop(print('No standard curve has been defined for this chip.'))
    
  } else if (chip_type == '48_48') {  # only the Gene Expression chip! Not Genotyping!
    stop(print('No standard curve has been defined for this chip.'))
    
  } else if (chip_type == '192_24') {
    stop(print('No standard curve has been defined for this chip.'))
    
  } else if (chip_type == 'plate') {
    # y = 4460x + 152
    k = 4460
    m = 152
    
  } else {
    stop(sprintf('Unknown chip: %s', chip_type))
  }
  
  concentration.df <- cbind(df, resorufin_value_norm = mapply(function(y, m, k) {(y-m)/k}, df$rox_value_norm, m, k))
  concentration.df$resorufin_value_norm_unit <- 'uM'
  
  # drop unneeded columns
  small.df <- concentration.df[, -which(names(concentration.df) %in% c('cycle', 'final_well', 'fam_value', 'rox_value_norm', 'assay_enzyme_conc', 'assay_enzyme_unit', 'sample_substrate_conc', 'sample_substrate_unit'))]
  
  return(small.df)
}



compute_slopes <- function(df){
  ## fit a linear regression to data to get the slope (how the signal changes with time), extract slope and the r2 ##
  ## I also compute the spearman correlation ##

  # first I need to bundle all the measurements from one well
  df.grouped <- df %>%
    group_by(assay_well, sample_well, assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num) %>%
    summarize(minutes=list(minutes), 
              resorufin_value_norm=list(resorufin_value_norm),
              resorufin_value_norm_unit=list(resorufin_value_norm_unit)) %>%
    filter(sapply(resorufin_value_norm, length) != 0) %>%
    mutate(slope = coef( lm(unlist(resorufin_value_norm)~unlist(minutes)) )[2],
           slope_unit = sprintf('%s min^-1', unique(unlist(resorufin_value_norm_unit))),
           r2 = summary( lm(unlist(resorufin_value_norm)~unlist(minutes)) )$r.squared,
           corr = cor.test(x=unlist(minutes), y=unlist(resorufin_value_norm), method='spearman')$estimate)
  
  # drop unneeded columns
  small.df <- df.grouped[, -which(names(df.grouped) %in% c('resorufin_value_norm', 'resorufin_value_norm_unit', 'minutes'))]
  
  return(small.df)
}



plot_slope_data_heatmaps <- function(df){
  ## Plot the slope of the linear models as well as the R2 values, which show the quality of the fit ##
  
  # make a subset of the data frames for plotting
  slope_data.df <- df[, which(names(df) %in% c('assay_well', 'sample_well', 'slope'))]
  
  # make wide format
  slope_data_wide.df <- dcast(slope_data.df, assay_well ~ sample_well, value.var='slope')
  
  # use first column as row names
  mat_data <- data.matrix(slope_data_wide.df[, 2:length(names(slope_data_wide.df))]) 
  rnames <- slope_data_wide.df[,'assay_well']   # assign labels in the "File" column to "rnames"
  rownames(mat_data) <- rnames 
  
  
  # Now make the heatmap to show the slope
  myImagePlot(mat_data, max=max(mat_data, na.rm=TRUE), min=0, title=('Slope (linear model)'))
  
  
  ### Plotting R2 of linear model (to assess the quality of the fit) ###
  
  # make a subset of the data frames for plotting 
  subset.df <- df[, which(names(df) %in% c('assay_well', 'sample_well', 'r2'))]
  
  
  # make wide format
  df.wide <- dcast(subset.df, assay_well ~ sample_well, value.var='r2')
  
  # use first column as row names
  mat_data <- data.matrix(df.wide[, 2:length(names(df.wide))]) 
  rnames <- df.wide[,'assay_well']   # assign labels in the "File" column to "rnames"
  rownames(mat_data) <- rnames 
  
  # Now make the heatmap
  myImagePlot(mat_data, max=max(mat_data, na.rm=TRUE), min=0, title=('R2 of slope (linear model fit)'))
  
  
  ### plotting spearman correlation ###
  
  # make a subset of the data frames for plotting 
  subset.df <- df[, which(names(df) %in% c('assay_well', 'sample_well', 'corr'))]
  
  
  # make wide format
  df.wide <- dcast(subset.df, assay_well ~ sample_well, value.var='corr')
  
  # use first column as row names
  mat_data <- data.matrix(df.wide[, 2:length(names(df.wide))]) 
  rnames <- df.wide[,'assay_well']   # assign labels in the "File" column to "rnames"
  rownames(mat_data) <- rnames 
  
  # Now make the heatmap
  myImagePlot(mat_data, max=max(mat_data, na.rm=TRUE), min=0, title=('Spearman correlation between time and signal'))
  
}



filter_r2 <- function(df, r2_cutoff, corr_cutoff){

  # First I need to bundle all the measurements from each enzyme concentration.
  # Also remove any samples with less than 5 observations (5 enzyme-substrate combinations),
  # as this seems like a minium number of data points for a kinetic curve.
  bundled.slopes.df <- df %>%
    group_by(assay_well, final_enzyme_conc, repeat_num, assay_enzyme, sample_substrate) %>%
    filter(sapply(list(slope), length) >= 5) %>%
    summarize(final_enzyme_unit=list(final_enzyme_unit),
              final_substrate_conc=list(final_substrate_conc),
              final_substrate_unit=list(final_substrate_unit),
              slope=list(slope),
              slope_unit=list(slope_unit),
              r2=list(r2),
              corr=list(corr))

  
  # now check, for each enzyme concentration, whether any of the R2 scores are below the cutoff
  # a bad R2 score indicates a poor fit of the data to the linear model
  bundled.slopes.nonegs.df <- bundled.slopes.df %>%
    group_by(assay_well, final_enzyme_conc, repeat_num, assay_enzyme, sample_substrate) %>%
    mutate(bad_r2 = any(unlist(r2) < r2_cutoff),
           bad_corr = any(unlist(corr) < corr_cutoff))
  
  # now actually filter the data
  bundled.slopes.nonegs.df <- bundled.slopes.nonegs.df[bundled.slopes.nonegs.df$bad_r2 == FALSE | bundled.slopes.nonegs.df$bad_corr == FALSE, ]
  
  # re-order the data
  small.df <- bundled.slopes.nonegs.df %>% 
    select(assay_well, assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num, slope, slope_unit)
  
  
  return(small.df)
}



calculate_slope_and_substrate_correlation <- function(df){
  # get the correlation between slope (reaction rate) and substrate concentration
  # this will give me a measure with which to evaluate which enzyme concentration gave the best data
  bundled.slopes.nonegs.spearman.df <- df %>%
    select(assay_well, final_enzyme_conc, final_enzyme_unit, repeat_num, assay_enzyme, sample_substrate, final_substrate_conc, final_substrate_unit, slope, slope_unit) %>%
    group_by(assay_well, final_enzyme_conc, repeat_num, assay_enzyme, sample_substrate) %>%
    mutate(slope_spearman_pval = cor.test(unlist(final_substrate_conc), unlist(slope),  method='spearman')$p.value,
           slope_spearman_rho = cor.test(unlist(final_substrate_conc), unlist(slope), method='spearman')[[4]])
  
  # re-order the data
  small.df <- bundled.slopes.nonegs.spearman.df %>% 
    select(assay_well, assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num, slope, slope_unit, slope_spearman_pval, slope_spearman_rho)
  
  
  return(small.df)
}



select_best_samples <- function(df){
  # choose the enzyme concentration that gave the best data
  
  
  ## The best one as a maximum sum of all repeats (use best concentration for all repeats) ##
  
  # get sum of spearman rhos
  sum.df <- df %>%
    group_by(assay_enzyme, final_enzyme_conc, sample_substrate) %>%
    summarise(slope_sum = sum(slope_spearman_rho))

  # figure out which enzyme concentration has the maximum sum
  best_concentration <- sum.df[which.max(sum.df$slope_sum), ]$final_enzyme_conc

  # slice the original data frame
  max.df <- df[df$final_enzyme_conc == best_concentration, ]
  
  return(max.df)
}


unbundle_data <- function(df){
  # separate out the values again
  filtered.df <- df[, -which(names(df) %in% c('slope_spearman_rho', 'bad_r2', 'slope_spearman_pval'))]
  
  unbundled.df <- unnest(filtered.df)
  
  # drop unneeded columns
  small.df <- unbundled.df[, -which(names(unbundled.df) %in% c('assay_well'))]
  
  # re-order the data
  small.df <- small.df %>% 
    select(assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num, slope, slope_unit)
  
  return(small.df)
}


subtract_no_substrate_slope <- function(df){
  # subtract the no substrate control from the others
  
  # first get the no substrate data
  slope.of.control.df <- df[df$final_substrate_conc==0, c('assay_enzyme', 'sample_substrate', 'repeat_num', 'slope')]
  
  # rename the slope value column
  colnames(slope.of.control.df)[length(colnames(slope.of.control.df))] <- 'slope_control'
  
  # do a left join to add on the controls to the samples
  final.data.df <- merge(df, slope.of.control.df, by=(c('assay_enzyme', 'sample_substrate', 'repeat_num')), all.x=TRUE)
  
  # now calculate the normalized
  final.data.df$activity <- final.data.df$slope - final.data.df$slope_control
  final.data.df$activity_unit <- final.data.df$slope_unit
  
  # re-order the data
  small.df <- final.data.df %>% 
    select(assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num, activity, activity_unit)
  
  return(small.df)
}


concentration_to_moles <- function(df, chip_type){
  # convert molar substrate to moles
  # n = c*v

  if (chip_type == 'flex_six'){
    # Sample volume = 7.9 nL
    #	Assay volume = 0.97 nL
    reaction_vol <- 8.87*10^-9   # 7.9 nL + 0.97 nL
    
  } else if (chip_type == '96_96') {
    # Sample volume = 6.7 nL
    #	Assay volume = 0.49 nL
    reaction_vol <- 7.19*10^-9   # 6.7 nL + 0.49 nL
    
  } else if (chip_type == '48_48') {  # only the Gene Expression chip! Not Genotyping!
    # Sample volume = 9.1 nL
    #	Assay volume = 1 nL
    reaction_vol <- 10.1*10^-9   # 9.1 nL + 1 nL
    
  } else if (chip_type == '192_24') {
    # Sample volume = 8 nL
    #	Assay volume = 1 nL
    reaction_vol <- 9*10^-9   # 8 nL + 1 nL
    
  } else if (chip_type == 'plate') {
    reaction_vol <- 20*10^-6 # 20 ul
    
  } else {
    stop(sprintf('Unknown chip: %s', chip_type))
  }
  
  moles.df <- cbind(df, activity_moles = mapply(function(c) c*reaction_vol, df$slope))
  moles.df$activity_moles_unit <- 'umol min^-1'

  return(moles.df)
}



normalize_to_protein <- function(df, chip_type){
  # take the rate values and normalize to per mg protein
  # the protein unit is in ng/ml so I need to convert to mg/l and then multiply by the volume used (in L)
  
  if (chip_type == 'flex_six'){
    # Sample volume = 7.9 nL
    #	Assay volume = 0.97 nL
    reaction_vol <- 8.87*10^-9   # 7.9 nL + 0.97 nL
    
  } else if (chip_type == '96_96') {
    # Sample volume = 6.7 nL
    #	Assay volume = 0.49 nL
    reaction_vol <- 7.19*10^-9   # 6.7 nL + 0.49 nL
    
  } else if (chip_type == '48_48') {  # only the Gene Expression chip! Not Genotyping!
    # Sample volume = 9.1 nL
    #	Assay volume = 1 nL
    reaction_vol <- 10.1*10^-9   # 9.1 nL + 1 nL
    
  } else if (chip_type == '192_24') {
    # Sample volume = 8 nL
    #	Assay volume = 1 nL
    reaction_vol <- 9*10^-9   # 8 nL + 1 nL
    
  } else if (chip_type == 'plate') {
    reaction_vol <- 20*10^-6 # 20 ul
      
  } else {
    stop(sprintf('Unknown chip: %s', chip_type))
  }
  
  df$activity_moles_protein <- df$activity_moles / ( (unique(df$final_enzyme_conc)*1000) * reaction_vol ) # * 1000 because the values are in mL and I want it in L
  df$activity_moles_protein_unit <- 'umol min^-1 mg^-1'
    
  return(df)
}


trim_data_top_range <- function(df, tolerance_at_max){

  # drop unnessecary columns
  #smaller.df <- df[, -which(names(df) %in% c('slope', 'slope_unit', 'activity', 'activity_unit', 'activity_moles', 'activity_moles_unit'))]
  
  # determine the substrate concentration at which the maximal slope is obtained
  max.df <- df %>%
    group_by(assay_enzyme, sample_substrate, repeat_num) %>%
    slice(which.max(activity_moles_protein))
  
  # drop the rox data and rename the substrate concentration column
  max.df <- max.df[, -which(names(max.df) %in% c('activity_moles', 'activity_moles_unit', 'activity_moles_protein', 'activity_moles_protein_unit'))]
  names(max.df)[grep("final_substrate_conc", colnames(max.df))] <- 'final_substrate_conc_at_max'
  
  # join the data frames
  final.data.filtered.df <- merge(df, max.df, by=(c('assay_enzyme', 'sample_substrate', 'final_substrate_unit', 'repeat_num', 'final_enzyme_conc', 'final_enzyme_unit')), all.x=TRUE)
  
  # now remove values on the plateau if they are lower than the maximum rate, within some tolerance
  final.data.filtered.both.df <- final.data.filtered.df %>%
    group_by(final_enzyme_conc, final_enzyme_unit, assay_enzyme, sample_substrate, repeat_num) %>%
    filter(final_substrate_conc <= final_substrate_conc_at_max |
             final_substrate_conc > final_substrate_conc_at_max & activity_moles_protein >= max(activity_moles_protein) - tolerance_at_max * max(activity_moles_protein))
  
  # drop the max column
  small.df <- final.data.filtered.both.df[, -which(names(final.data.filtered.both.df) %in% c('final_substrate_conc_at_max'))]
  
  # re-order the data
  small.df <- small.df %>% 
    select(assay_enzyme, final_enzyme_conc, final_enzyme_unit, sample_substrate, final_substrate_conc, final_substrate_unit, repeat_num, activity_moles, activity_moles_unit, activity_moles_protein, activity_moles_protein_unit)
  
  return(small.df)
}



calculate_kinetics <- function(df){

  ## Need to generalize this for many enzymes.... ##
  
  
  # "velocity = Vmax times S divided by (Km plus S)", stored in MMcurve
  MMcurve <- formula(v ~ Vmax*S/(Km+S))
  
  
  # fit the equation to data using non-linear least squares
  kinetic_data <- df[, which(names(df) %in% c('assay_enzyme', 'activity_moles_protein', 'final_substrate_conc'))]
  
  # rename to make it work with the formula
  slope_index <- grep("activity_moles_protein", colnames(kinetic_data))
  names(kinetic_data)[slope_index] <- 'v'
  
  substrate_index <- grep("final_substrate_conc", colnames(kinetic_data))
  names(kinetic_data)[substrate_index] <- 'S'
  
  
  # estimate Vmax
  vmax_estimate <- max(kinetic_data$v)
  
  
  # estimate km at concentration closest to half of Vmax
  index_of_half_vmax <- which.min(abs(kinetic_data$v - max(kinetic_data$v)/2))
  km_estimate = kinetic_data$S[index_of_half_vmax]
  
  
  # perform non-linear least squares regression
  bestfit <- nls(MMcurve, kinetic_data, start=list(Vmax=vmax_estimate, Km=km_estimate))
  
  
  # generate new x-values that can be used to predict y using the model (for making the smooth curve)
  xFit <- seq(0, max(df$final_substrate_conc), length.out = 100)
  yFit <- predict(bestfit, newdata=list(S = xFit))
  
  return(c(data.frame(xFit=xFit, yFit=yFit), coef(bestfit)['Km'], coef(bestfit)['Vmax']))
}