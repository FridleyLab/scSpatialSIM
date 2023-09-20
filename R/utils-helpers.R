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
    val = grid[val,1:2] %>% as.numeric()
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
  #cat("running\n")
  cn = 1
  while (TRUE) {
    vals <- stats::runif(num_vals, min_val, max_val)
    if (sum(vals) < sum_val) {
      break
    }
    #have to slowly slowly decrease the threshold otherwise samples get into a situation where,
    #for example, trying to randomly generate 3 numbers between 0.1 and 0.35, that add up to 0.3
    #and unles the numbers are identical it'll never achieve it
    if(cn == 10000){
      #cat("tick\n")
      min_val = min_val - (min_val * 0.1)
      cn = 1
    }
    cn = cn + 1
  }
  dif = sum_val - sum(vals)
  vals[1] = vals[1] + dif
  return(vals)
}

make_dis = function(spat, positive_mean, negative_mean, positive_sd, negative_sd){
  cells = grep("Cell", colnames(spat), value = TRUE)
  dat = spat %>% dplyr::arrange(get(cells))
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

gaussian_kernel_shift = function(kern, shift, win_limits){
  if(nrow(kern) == 1){
    #only single center
    #move towards opposite side?
    n_x = win_limits[2] - (kern$x[1] - win_limits[1])
    n_y = win_limits[4] - (kern$y[1] - win_limits[2])
    #shift towards the new point by shift amount
    return(kern %>%
             dplyr::mutate(x = (n_x) + (x - n_x) * (1-shift),
                           y = (n_y) + (y - n_y) * (1-shift)))
  } else { #if(nrow(kern) <= 3)
    #going to go with a different method for introducing a bit more randomness
    #number of kernels is less than or equal to 5, find all intersections and sample nros
    #perform the dirichlet boundary calculation
    pp_obj = spatstat.geom::ppp(x = kern$x, y = kern$y, window =
                                  spatstat.geom::owin(win_limits[1:2], win_limits[3:4]))
    dirichlet_boundary = spatstat.geom::dirichlet(pp_obj)
    #pull out only the boundary and vertex information
    dirich_boundary_extract = lapply(dirichlet_boundary$tiles, function(x){
      x$bdry[[1]]
    })
    #find places that are for more than 1 point (probably all)
    dirichlet_intersect = lapply(seq(dirich_boundary_extract), function(point_num){
      others = lapply(dirich_boundary_extract[-point_num], data.frame) %>%
        do.call(dplyr::bind_rows, .) %>%
        dplyr::distinct() #%>%
      #dplyr::filter(!(x %in% c(0, 10)), !( y %in% c(0,10)))
      point_bounds = dirich_boundary_extract[[point_num]] %>% data.frame()
      dplyr::inner_join(point_bounds, others, by = dplyr::join_by(x, y))
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::distinct()
    #sample them to kern rows
    dirichlet_intersect = dirichlet_intersect[sample(seq(nrow(dirichlet_intersect)), nrow(kern)),]
    #find nearest real points
    distance_matrix = proxy::dist(x = kern[,1:2], y = dirichlet_intersect)
    #get minimum distance
    x_diff = apply(distance_matrix, 1, which.min)
    y_diff = apply(distance_matrix, 1, which.min)
    #calculate new spot?
    vec_dat = kern %>%
      dplyr::select(x, y) %>%
      dplyr::mutate(xend = dirichlet_intersect[x_diff, "x"],
                    yend = dirichlet_intersect[y_diff, "y"])
    vec_dat2 = vec_dat %>%
      dplyr::mutate(xnew = (xend) + (x - xend) * (1-shift),
                    ynew = (yend) + (y - yend) * (1-shift))
    kern2 = kern
    kern2$x = vec_dat2$xnew
    kern2$y = vec_dat2$ynew
    return(kern2)
  }
  # else {
  #   #over 5, just use the resulting non-border locations since they *should* be spread out enough
  #   #perform the dirichlet boundary calculation
  #   pp_obj = spatstat.geom::ppp(x = kern$x, y = kern$y, window =
  #                                 spatstat.geom::owin(win_limits[1:2], win_limits[3:4]))
  #   dirichlet_boundary = spatstat.geom::dirichlet(pp_obj)
  #   #pull out only the boundary and vertex information
  #   dirich_boundary_extract = lapply(dirichlet_boundary$tiles, function(x){
  #     x$bdry[[1]]
  #   })
  #   #find places that are for more than 1 point (probably all)
  #   dirichlet_intersect = lapply(seq(dirich_boundary_extract), function(point_num){
  #     others = lapply(dirich_boundary_extract[-point_num], data.frame) %>%
  #       do.call(dplyr::bind_rows, .) %>%
  #       dplyr::distinct() %>%
  #       dplyr::filter(!(x %in% win_limits[1:2]), !( y %in% win_limits[3:4]))
  #     point_bounds = dirich_boundary_extract[[point_num]] %>% data.frame()
  #     dplyr::inner_join(point_bounds, others, by = dplyr::join_by(x, y))
  #   }) %>%
  #     do.call(dplyr::bind_rows, .) %>%
  #     dplyr::distinct()
  #   #find nearest real points
  #   distance_matrix = proxy::dist(x = kern[,1:2], y = dirichlet_intersect)
  #   #get minimum distance
  #   x_diff = apply(distance_matrix, 1, which.min)
  #   y_diff = apply(distance_matrix, 1, which.min)
  #   #calculate new spot?
  #   vec_dat = kern %>%
  #     dplyr::select(x, y) %>%
  #     dplyr::mutate(xend = dirichlet_intersect[x_diff, "x"],
  #                   yend = dirichlet_intersect[y_diff, "y"])
  #   vec_dat2 = vec_dat %>%
  #     dplyr::mutate(xnew = (xend) + (x - xend) * (1-shift),
  #                   ynew = (yend) + (y - yend) * (1-shift))
  #   kern2 = kern
  #   kern2$x = vec_dat2$xnew
  #   kern2$y = vec_dat2$ynew
  #   return(kern2)
  # }
}
