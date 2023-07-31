# ppp_to_spatial = function(ppp_list, phenotype_names){
#   simulated_cell_positive = lapply(names(ppp_list), function(core){
#     #convert ppp into long formatted dataframe
#     whole_core = data.frame(x = ppp_list[[core]]$x,
#                             y = ppp_list[[core]]$y,
#                             marks = ppp_list[[core]]$marks) %>%
#       dplyr::mutate(label = "Tumor",
#                     pos = 1) %>%
#       tidyr::spread(key = marks, value = pos) %>%
#       dplyr::mutate(background_pos = 1) %>%
#       dplyr::mutate_at(.vars = grep(paste0(phenotype_names,collapse = "|"), colnames(.), value = T), .funs = ~ ifelse(is.na(.x), 0, 1)) %>%
#       dplyr::mutate(dplyr::across(grep(paste0(phenotype_names, collapse = "|"), colnames(.), value = T),
#                                   list(intensity = function(x){
#                                     ifelse(x == 1,
#                                            rnorm(dplyr::n(), 10, 2),
#                                            rnorm(dplyr::n(), 3, 1))
#                                   }),
#                                   .names = "{gsub('pos', 'intensity', {col}, fixed = TRUE)}"))
#
#     #spread_wide
#     return(whole_core)
#   })
#
#   names(simulated_cell_positive) = names(ppp_list)
#   return(simulated_cell_positive)
# }

#kernal for generating distributions around random points for k centers
gaussian_kernel = function(k = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, sdmin = 1/2, sdmax = 2){
  #simulate centers centers of groups
  center_peak = data.frame(x = stats::runif(k,xmin,xmax), y = stats::runif(k,ymin,ymax),
                           sd.x = stats::runif(k,sdmin,sdmax), sd.y = stats::runif(k, sdmin, sdmax))
  #spearman correlations for how the distribution falls in 2D space
  center_peak$rho = sapply(1:k, function(a){
    stats::runif(1,-center_peak$sd.x[a]*center_peak$sd.y[a], center_peak$sd.x[a]*center_peak$sd.y[a])
  })
  return(center_peak)
}

#randomly generates area to remove based on hole_prob for min and max
generate_holes = function(xmin = 0, xmax = 10, ymin = 0,
                          ymax = 10, sdmin = 1/2, sdmax = 2, hole_prob = c(0.2, 0.35)){
  #Random select the percent to total area missing
  percent_area_missing = stats::runif(1,hole_prob[1],hole_prob[2])
  #Random select the number of holes, can range from 1 to floor(percent_area_missing*10)
  #for example if percent_area_missing = 0.213974, then there can be either 1 or 2 holes
  num_holes = sample(2:floor(percent_area_missing*10),1)
  if(percent_area_missing < 0.2 | num_holes == 1){
    num_holes = 1
    area_allocation = percent_area_missing
  }else{
    #RandVec is from the Surrogate package this function can create a vector such that
    #it's sum is s, but each value in the vector is in the interval [a,b]
    #This is a key function the package will have to be a dependency if we include
    #this simulation code. The Surrogate package was not easy to install
    area_allocation = generate_sum_vector(num_holes, 0.1, percent_area_missing, percent_area_missing)
  }
  #generate table for hole centers
  center_hole = data.frame(x = stats::runif(num_holes,xmin,xmax), y = stats::runif(num_holes,ymin,ymax), rho = 0,
                           sd.x = stats::runif(num_holes,sdmin,sdmax), sd.y = stats::runif(num_holes,sdmin,sdmax))
  #add random spearman correlations again guessing for direction of skewness but sd should also provide that
  center_hole$rho = sapply(1:num_holes, function(a){
    stats::runif(1,-center_hole$sd.x[a]*center_hole$sd.y[a], center_hole$sd.x[a]*center_hole$sd.y[a])
  })
  #The max distance was derived by setting the area_allocation = area_circle/total_area (side^2) and solving
  #for r.
  return(cbind.data.frame(center_hole, area = area_allocation,
                          max_dist = sqrt(area_allocation * (xmax - xmin) *(ymax - ymin) / pi)))
}

CalculateGrid = function(grid, gauss_tab, cores){
  parallel::mclapply(1:nrow(grid), function(val){
    val = grid[val,] %>% as.numeric()
    sapply(1:nrow(gauss_tab), function(a){
      diff = c((val[1] - gauss_tab$x[a]), (val[2] - gauss_tab$y[a]))
      sigma = matrix(c(gauss_tab$sd.x[a]^2, rep(gauss_tab$rho[a], 2),
                       gauss_tab$sd.y[a]^2), nrow = 2, ncol = 2)
      z = t(diff) %*% solve(sigma) %*% diff
      closest = exp(-z)
    }) %>% max()
  }, mc.cores = cores) %>% unlist()
}

CalculateGridHoles = function(grid, gauss_tab, cores){
  parallel::mclapply(1:nrow(grid), function(val){
    val = grid[val,] %>% as.numeric()
    mins = sapply(1:nrow(gauss_tab), function(a){
      diff = c((val[1] - gauss_tab$x[a]), (val[2] - gauss_tab$y[a]))
      sigma = matrix(c(gauss_tab$sd.x[a]^2, rep(gauss_tab$rho[a], 2),
                       gauss_tab$sd.y[a]^2), nrow = 2, ncol = 2)
      z = t(diff) %*% solve(sigma) %*% diff #solve finds in the inverse of a single matrix without the righthand vector
    })
    data.frame(`Closest Hole` = which.min(mins),
      `Hole Z` = min(mins), check.names = FALSE)
  }, mc.cores = cores) %>% do.call(dplyr::bind_rows, .)
}

is.empty <- function(obj, slot.name) {
  if (slot.name %in% methods::slotNames(class(obj))) {
    return(is.null(methods::slot(obj, slot.name)) || length(methods::slot(obj, slot.name)) == 0)
  } else {
    stop(paste0("Object does not have a '", slot.name, "' slot."))
  }
}

scale_probs = function(x, probs){
  s = (x - min(x))/(max(x) - min(x))
  (probs[2] - probs[1]) * s + probs[1]
}

replace_na <- function(x, y) {
  ifelse(is.na(y), x, y)
}

generate_sum_vector <- function(num_vals, min_val, max_val, sum_val) {
  vals <- numeric(num_vals)
  while (TRUE) {
    vals <- stats::runif(num_vals, min_val, max_val)
    if (sum(vals) < sum_val) {
      break
    }
  }
  dif = sum_val - sum(vals)
  vals[1] = vals[1] + dif
  return(vals)
}




make_dis = function(spat, positive_mean, negative_mean, positive_sd, negative_sd){
  cells = grep("Cell", colnames(spat), value = TRUE)
  dat = spat %>% dplyr::arrange(get(cell))
  #ugh for loop
  for(cell_n in seq(cells)){
    cell = cells[cell_n]
    counts = data.frame(table(spat[[cell]]))
    negative_cells = stats::rnorm(counts$Freq[1], negative_mean, negative_sd)
    positive_cells = stats::rnorm(counts$Freq[2], positive_mean, positive_sd)
    dat[paste0("Cell ", cell_n, " Var")] = c(negative_cells, positive_cells)
  }

  return(dat)
}

