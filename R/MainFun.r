#' iNterpolation and EXTrapolation for Beta diversity 
#' 
#' \code{iNEXTBetaDiv}: Interpolation and extrapolation of Beta diversity with order q
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-assemblages abundance matrix.\cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-sampling units.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being \code{0} (non-detection) or \code{1} (detection).
#' @param base Sample-sized-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}) or coverage-based rarefaction and extrapolation for gamma, alpha and beta diversity (\code{base = "coverage"}). Default is \code{base = "coverage"}.
#' @param level A numerical vector specifying the particular value of sample coverage (between 0 and 1 when \code{base = "coverage"}) or sample size (\code{base = "size"}). \code{level = 1} (\code{base = "coverage"}) means complete coverage (the corresponding diversity represents asymptotic diversity).\cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the gamma and alpha diversity estimates up to double the reference sample size. If \code{base = "coverage"} and \code{level = NULL}, then this function computes the gamma and alpha diversity estimates up to one (for \code{q = 1, 2}) or up to the coverage of double the reference sample size (for \code{q = 0});\cr 
#' the corresponding beta diversity is computed up to the same maximum coverage as the alpha diversity.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{20}. If more accurate results are required, set \code{nboot = 100} (or \code{200}).
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' 
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import colorRamps
#' @import iNEXT.3D
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tidytree
#' 
#' @return A list of seven data frames, three for diversity and four for dissimilarity. Also, if more than one dataset is included, the list will contain one component for each dataset, and each component will be the list of seven data frames.
#' 
#' @examples
#' ## abundance data & coverage-based
#' data(beetle_abu)
#' output1 = iNEXTBetaDiv(data = beetle_abu, datatype = 'abundance',
#'                        nboot = 20, conf = 0.95)
#' output1
#' 
#' 
#' ## incidence data & coverage-based
#' data(beetle_inc)
#' output2 = iNEXTBetaDiv(data = beetle_inc, datatype = 'incidence_raw',
#'                        nboot = 20, conf = 0.95)
#' output2
#' 
#' 
#' ## abundance data & size-based
#' data(beetle_abu)
#' output3 = iNEXTBetaDiv(data = beetle_abu, datatype = 'abundance', base = 'size',
#'                        nboot = 20, conf = 0.95)
#' output3
#' 
#' 
#' ## incidence data & size-based
#' data(beetle_inc)
#' output4 = iNEXTBetaDiv(data = beetle_inc, datatype = 'incidence_raw', base = 'size',
#'                        nboot = 20, conf = 0.95)
#' output4
#' 
#' 
#' @export
iNEXTBetaDiv = function(data, q = c(0, 1, 2), datatype = 'abundance', base = "coverage", level = NULL, nboot = 20, conf = 0.95) {
  max_alpha_coverage = F
  if (datatype == 'abundance') {
    
    if( class(data) == "data.frame" | class(data) == "matrix" ) data = list(Region_1 = data)
    
    if(class(data) == "list"){
      if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
    }
    
  }
  
  if (datatype == 'incidence_raw') {
    
    if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
    Ns = sapply(data, length)
    data_list = data
    
  }
  
  
  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  
  trunc = ifelse(is.null(level), T, F)
  if ( is.null(level) & base == 'coverage' ) level = seq(0.5, 1, 0.025) else if ( base == 'size' ) {
    if ( is.null(level) ) {
      
      if (datatype == "abundance") {
        endpoint <- sapply(data_list, function(x) 2*sum(x))
      } else if (datatype == "incidence_raw") {
        endpoint <- sapply(data_list, function(x) 2*ncol(x[[1]]))
      }
      
      level <- lapply(1:length(data_list), function(i) {
        
        if(datatype == "abundance") {
          ni <- sum(data_list[[i]])
        }else if(datatype == "incidence_raw"){
          ni <- ncol(data_list[[i]][[1]])
        }
        
        mi <- floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint[i], length.out = 20)))
      })
      
    } else {
      
      if (class(level) == "numeric" | class(level) == "integer" | class(level) == "double") {
        level <- list(level = level)
      }
      
      if (length(level) != length(data_list)) level <- lapply(1:length(data_list), function(x) level[[1]])
      
      level <- lapply(1:length(data_list), function(i) {
        
        if (datatype == "abundance") {
          ni <- sum(data_list[[i]])
        } else if (datatype == "incidence_raw"){
          ni <- ncol(data_list[[i]][[1]])
        }
        
        if( sum(level[[i]] == ni) == 0 ) mi <- sort(c(ni, level[[i]])) else mi <- level[[i]]
        unique(mi)
      })
    }
  }
  
  for_each_region = function(data, region_name, N) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
      
      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level<1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
      
      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level < 1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq, i, datatype='incidence_freq'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq, i, datatype='incidence_freq'))
      
    }
    
    
    
    if (datatype == 'abundance') {
      
      gamma = lapply(1:length(level), function(i){
        estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
      }) %>% do.call(rbind,.)
      
      alpha = lapply(1:length(level), function(i){
        estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
      }) %>% do.call(rbind,.)
      
    }
    
    if (datatype == 'incidence_raw') {
      
      gamma = lapply(1:length(level), function(i){
        estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
      }) %>% do.call(rbind,.)
      
      alpha = lapply(1:length(level), function(i){
        estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
      }) %>% do.call(rbind,.)
      
      
      
    }
    
    gamma = (cbind(level = rep(level, each=length(q)), gamma[,-c(1,2,8,9)]) %>% 
               mutate(Method = ifelse(level>=ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
    )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
    
    # for (i in 0:2) gamma$Order[gamma$Order==paste0('q = ', i)] = i
    # gamma$Order = as.numeric(gamma$Order)
    
    if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level > 0
    gamma = gamma[under_max_alpha,]
    
    
    
    alpha = (cbind(level = rep(level, each = length(q)), alpha[,-c(1,2,8,9)]) %>% 
               mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
    )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))
    
    alpha$Estimate = alpha$Estimate / N
    
    # for (i in 0:2) alpha$Order[alpha$Order == paste0('q = ', i)] = i
    # alpha$Order = as.numeric(alpha$Order)
    
    alpha = alpha[under_max_alpha,]
    
    
    
    beta = alpha
    beta$Estimate = gamma$Estimate/alpha$Estimate
    beta[beta == "Observed"] = "Observed_alpha"
    beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate), 
                                   Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))
    
    C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
    U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
    V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
    S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))
    
    if(nboot>1){
      
      # cl = makeCluster(cluster_numbers)
      # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
      #                     'datatype', 'data_2D'))
      # clusterEvalQ(cl, library(tidyverse, magrittr))
      
      # plan(sequential)
      # plan(multiprocess)
      
      # se = parSapply(cl, 1:nboot, function(i){
      
      # start = Sys.time()
      se = future_lapply(1:nboot, function(i){
        
        if (datatype == 'abundance') {
          
          bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
          bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
          
          bootstrap_data_gamma = rowSums(bootstrap_sample)
          bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
          bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
          bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
          
          gamma = lapply(1:length(level), function(i){
            estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
          }) %>% do.call(rbind,.)
          
          alpha = lapply(1:length(level), function(i){
            estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
          }) %>% do.call(rbind,.)
          
          beta_obs = (obs3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) / 
                        (obs3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) / N)) %>% unlist()
          
        }
        
        if (datatype == 'incidence_raw') {
          
          bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
          
          raw = lapply(1:ncol(bootstrap_population), function(j){
            
            lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
            
          })
          
          gamma = Reduce('+', raw)
          gamma[gamma > 1] = 1
          bootstrap_data_gamma_freq = c(n, rowSums(gamma))
          
          bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
          
          bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
          bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
          
          gamma = lapply(1:length(level), function(i){
            estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
          }) %>% do.call(rbind,.)
          
          alpha = lapply(1:length(level), function(i){
            estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
          }) %>% do.call(rbind,.)
          
          beta_obs = (obs3D(as.numeric(bootstrap_data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) / 
                        (obs3D(as.numeric(bootstrap_data_alpha_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) / N)) %>% unlist()
          
        }
        
        gamma = gamma[,c(6,3,7)]$qD[under_max_alpha]
        
        alpha = alpha[,c(6,3,7)]$qD[under_max_alpha]
        alpha = alpha / N
        
        beta = gamma/alpha
        
        gamma = c(gamma, rep(0, length(q)))
        alpha = c(alpha, rep(0, length(q)))
        beta = c(beta, beta_obs)
        
        order = rep(q, length(level) + 1)[under_max_alpha]
        
        beta = data.frame(Estimate=beta, order)
        
        C = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
        U = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
        V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
        S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
        
        beta = beta$Estimate
        
        cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
        
        # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
      }) %>% abind(along = 3) %>% apply(1:2, sd)
      # end = Sys.time()
      # end - start
      
      # stopCluster(cl)
      # plan(sequential)
      
    } else {
      
      se = matrix(0, ncol = 7, nrow = nrow(beta))
      colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      se = as.data.frame(se)
      
    }
    
    
    if (nboot>0 & datatype == 'abundance') {
      gamma.se = estimate3D(data_gamma, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
      
      alpha.se = estimate3D(data_alpha, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
      alpha.se = alpha.se / N
      
    } else if (nboot>0 & datatype == 'incidence_raw') {
      
      gamma.se = estimate3D(data_gamma_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
      
      alpha.se = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
      alpha.se = alpha.se / N
      
    }
    
    se[1:( length(level) * length(q) ), 'gamma'] = gamma.se
    
    se[1:( length(level) * length(q) ), 'alpha'] = alpha.se
    
    
    se = as.data.frame(se)
    
    gamma = gamma %>% mutate(s.e. = se$gamma[1:(nrow(se) - length(q))],
                             LCL = Estimate - tmp * se$gamma[1:(nrow(se) - length(q))],
                             UCL = Estimate + tmp * se$gamma[1:(nrow(se) - length(q))],
                             Region = region_name)
    
    alpha = alpha %>% mutate(s.e. = se$alpha[1:(nrow(se) - length(q))],
                             LCL = Estimate - tmp * se$alpha[1:(nrow(se) - length(q))],
                             UCL = Estimate + tmp * se$alpha[1:(nrow(se) - length(q))],
                             Region = region_name)
    
    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Region = region_name)
    
    C = C %>% mutate(        s.e. = se$C,
                             LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Region = region_name)
    
    
    U = U %>% mutate(        s.e. = se$U,
                             LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Region = region_name)
    
    V = V %>% mutate(        s.e. = se$V,
                             LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Region = region_name)
    
    S = S %>% mutate(        s.e. = se$S,
                             LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Region = region_name)
    
    if (trunc) {
      
      gamma = gamma %>% filter(!(Order==0 & round(Size)>2*n))
      
      alpha = alpha %>% filter(!(Order==0 & round(Size)>2*n))
      
      beta  = beta  %>% filter(!(Order==0 & round(Size)>2*n))
      
       C    =  C    %>% filter(!(Order==0 & round(Size)>2*n))
      
       U    =  U    %>% filter(!(Order==0 & round(Size)>2*n))
      
       V    =  V    %>% filter(!(Order==0 & round(Size)>2*n))
      
       S    =  S    %>% filter(!(Order==0 & round(Size)>2*n))
      
    }
    
    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)
    
  }
  
  for_each_region.size = function(data, region_name, N, level) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    if (datatype == 'abundance') {
      
      gamma = estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
      
      alpha = estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
      
    }
    
    if (datatype == 'incidence_raw') {
      
      gamma = estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
      
      alpha = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
      
    }
    
    se = cbind(gamma$s.e., alpha$s.e. / N)
    colnames(se) = c("gamma", "alpha")
    se = as.data.frame(se)
    se[is.na(se)] = 0
    
    gamma = (cbind(Size = rep(level, each=length(q)), gamma[,-c(1,2,8,9)]) %>% 
               mutate(Method = ifelse(Size>=ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
    )[,c(5,3,2,4,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))
    
    
    alpha = (cbind(Size = rep(level, each = length(q)), alpha[,-c(1,2,8,9)]) %>% 
               mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
    )[,c(5,3,2,4,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))
    
    alpha$Estimate = alpha$Estimate / N
    
    # if(nboot>1){
    #   
    #   se = future_lapply(1:nboot, function(i){
    #     
    #     if (datatype == 'abundance') {
    #       
    #       bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
    #       bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
    #       
    #       bootstrap_data_gamma = rowSums(bootstrap_sample)
    #       bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
    #       bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
    #       bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
    #       
    #       gamma = lapply(1:length(level), function(i){
    #         estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level[i], nboot = 0)
    #       }) %>% do.call(rbind,.)
    #       
    #       alpha = lapply(1:length(level), function(i){
    #         estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level[i], nboot = 0)
    #       }) %>% do.call(rbind,.)
    #       
    #     }
    #     
    #     if (datatype == 'incidence_raw') {
    #       
    #       bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
    #       
    #       raw = lapply(1:ncol(bootstrap_population), function(j){
    #         
    #         lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
    #         
    #       })
    #       
    #       gamma = Reduce('+', raw)
    #       gamma[gamma > 1] = 1
    #       bootstrap_data_gamma_freq = c(n, rowSums(gamma))
    #       
    #       bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
    #       
    #       bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
    #       bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
    #       
    #       gamma = lapply(1:length(level), function(i){
    #         estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level[i], nboot = 0)
    #       }) %>% do.call(rbind,.)
    #       
    #       alpha = lapply(1:length(level), function(i){
    #         estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level[i], nboot = 0)
    #       }) %>% do.call(rbind,.)
    #       
    #     }
    #     
    #     gamma = gamma$qD
    #     
    #     alpha = alpha$qD
    #     alpha = alpha / N
    #     
    #     cbind(gamma, alpha) %>% as.matrix
    #     
    #   }) %>% abind(along = 3) %>% apply(1:2, sd)
    #   
    # } else {
    #   
    #   se = matrix(0, ncol = 2, nrow = nrow(gamma))
    #   colnames(se) = c("gamma", "alpha")
    #   se = as.data.frame(se)
    #   
    # }
    
    se = as.data.frame(se)
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Region = region_name)
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Region = region_name)
    
    if (datatype == 'incidence_raw') {
      colnames(gamma)[colnames(gamma) == 'Size'] = 'nT'
      colnames(alpha)[colnames(alpha) == 'Size'] = 'nT'
    }
    
    list(gamma = gamma, alpha = alpha)
    
  }
  
  if (base == 'coverage') output = lapply(1:length(data_list), function(i) for_each_region(data = data_list[[i]], region_name = region_names[i], N = Ns[i]))
  if (base == 'size') output = lapply(1:length(data_list), function(i) for_each_region.size(data = data_list[[i]], region_name = region_names[i], N = Ns[i], level = level[[i]]))
  names(output) = region_names
  
  return(output)
  
}



#' ggplot for Beta diversity
#' 
#' \code{ggiNEXTBetaDiv}: ggplot for Interpolation and extrapolation of Beta diversity with order q
#' 
#' @param output the output from iNEXTBetaDiv
#' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ; \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' @param scale Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.
#' @param main The title of the plot.
#' @param transp a value between 0 and 1 controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.
#' 
#' @return a figure for Beta diversity or dissimilarity diversity.
#' 
#' @examples
#' ## abundance data & coverage-based
#' data(beetle_abu)
#' output1 = iNEXTBetaDiv(data = beetle_abu, datatype = 'abundance', 
#'                        nboot = 20, conf = 0.95)
#' 
#' ggiNEXTBetaDiv(output1, type = 'B', scale = 'free', main = NULL, transp = 0.4)
#' ggiNEXTBetaDiv(output1, type = 'D', scale = 'free', main = NULL, transp = 0.4)
#' 
#' 
#' ## incidence data & coverage-based
#' data(beetle_inc)
#' output2 = iNEXTBetaDiv(data = beetle_inc, datatype = 'incidence_raw', 
#'                        nboot = 20, conf = 0.95)
#' 
#' ggiNEXTBetaDiv(output2, type = 'B', scale = 'free', main = NULL, transp = 0.4)
#' ggiNEXTBetaDiv(output2, type = 'D', scale = 'free', main = NULL, transp = 0.4)
#' 
#' 
#'  ## abundance data & size-based
#' data(beetle_abu)
#' output3 = iNEXTBetaDiv(data = beetle_abu, datatype = 'abundance', base = 'size', 
#'                        nboot = 20, conf = 0.95)
#' ggiNEXTBetaDiv(output3, scale = 'free', main = NULL, transp = 0.4)
#' 
#' 
#' ## incidence data & size-based
#' data(beetle_inc)
#' output4 = iNEXTBetaDiv(data = beetle_inc, datatype = 'incidence_raw', base = 'size', 
#'                        nboot = 20, conf = 0.95)
#' ggiNEXTBetaDiv(output4, scale = 'free', main = NULL, transp = 0.4)
#' 
#' 
#' @export
ggiNEXTBetaDiv = function(output, type = 'B', scale = 'free', main = NULL, transp = 0.4){
  
  if (length(output[[1]]) == 7) {
    if (type == 'B'){
      
      gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
      alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
      beta =  lapply(output, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
      beta = beta %>% filter(Method != 'Observed')
      beta[beta == 'Observed_alpha'] = 'Observed'
      
      # # Dropping out the points extrapolated over double reference size
      # gamma1 = data.frame() ; alpha1 = data.frame() ; beta1 = data.frame()
      # 
      # for(i in 1:length(unique(gamma$Region))){
      #   
      #   Gamma <- gamma %>% filter(Region==unique(gamma$Region)[i]) ; ref_size = unique(Gamma[Gamma$Method=="Observed",]$Size)
      #   Gamma = Gamma %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   
      #   Alpha <- alpha %>% filter(Region==unique(gamma$Region)[i]) ; Alpha = Alpha %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   Beta <- beta %>% filter(Region==unique(gamma$Region)[i]) ; Beta = Beta %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   
      #   gamma1 = rbind(gamma1,Gamma) ; alpha1 = rbind(alpha1,Alpha) ; beta1 = rbind(beta1,Beta)
      #   
      # }
      # 
      # gamma = gamma1 ; alpha = alpha1 ; beta= beta1
      
      df = rbind(gamma, alpha, beta)
      for (i in unique(gamma$Order)) df$Order[df$Order == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
      
      id_obs = which(df$Method == 'Observed')
      
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$level = new$level - 0.0001
        new$Method = 'Interpolated'
        
        newe = df[id_obs[i],]
        newe$level = newe$level + 0.0001
        newe$Method = 'Extrapolated'
        
        df = rbind(df, new, newe)
        
      }
      
      ylab = "Taxonomic diversity"
      
    }
    
    if (type == 'D'){
      
      C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
      U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
      V = lapply(output, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
      S = lapply(output, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
      C = C %>% filter(Method != 'Observed')
      U = U %>% filter(Method != 'Observed')
      V = V %>% filter(Method != 'Observed')
      S = S %>% filter(Method != 'Observed')
      C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = V[V == 'Observed_alpha'] = S[S == 'Observed_alpha'] = 'Observed'
      
      # # Dropping out the points extrapolated over double reference size
      # c1 = data.frame() ; u1 = data.frame() ; v1 = data.frame() ; s1 = data.frame()
      # 
      # for(i in 1:length(unique(C$Region))){
      #   
      #   CC <- C %>% filter(Region==unique(C$Region)[i]) ; ref_size = unique(CC[CC$Method=="Observed",]$Size)
      #   CC = CC %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   
      #   UU <- U %>% filter(Region==unique(C$Region)[i]) ; UU = UU %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   VV <- V %>% filter(Region==unique(C$Region)[i]) ; VV = VV %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   SS <- S %>% filter(Region==unique(C$Region)[i]) ; SS = SS %>% filter(!(Order==0 & round(Size)>2*ref_size))
      #   
      #   c1 = rbind(c1,CC) ; u1 = rbind(u1,UU) ; v1 = rbind(v1,VV) ; s1 = rbind(s1,SS)
      #   
      # }
      # 
      # C = c1 ; U = u1 ; V = v1 ; S = s1
      
      
      df = rbind(C, U, V, S)
      for (i in unique(C$Order)) df$Order[df$Order == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
      
      id_obs = which(df$Method == 'Observed')
      
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$level = new$level - 0.0001
        new$Method = 'Interpolated'
        
        newe = df[id_obs[i],]
        newe$level = newe$level + 0.0001
        newe$Method = 'Extrapolated'
        
        df = rbind(df, new, newe)
        
      }
      
      ylab = "Taxonomic dissimilarity"
      
    }
    
    lty = c(Interpolated = "solid", Extrapolated = "dashed")
    df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))
    
    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolated" & round(Size) %in% double_size)
    
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    
    ggplot(data = df, aes(x = level, y = Estimate, col = Region)) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) + 
      geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) + scale_linetype_manual(values = lty) +
      # geom_line(lty=2) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5)+
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) + 
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) + 
      scale_colour_manual(values = cbPalette) + 
      scale_fill_manual(values = cbPalette) + 
      facet_grid(div_type ~ Order, scales = scale) +
      theme_bw() + 
      theme(legend.position = "bottom", legend.title = element_blank()) +
      labs(x = 'Sample coverage', y = ylab, title = main)
    
  } else if (length(output[[1]]) == 2){
    
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    
    gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
    
    if ('nT' %in% colnames(gamma)) {
      xlab = 'Number of sampling units'
      colnames(gamma)[colnames(gamma) == 'nT'] = 'Size'
      colnames(alpha)[colnames(alpha) == 'nT'] = 'Size'
    } else xlab = 'Number of individuals'
    
    df = rbind(gamma, alpha)
    for (i in unique(gamma$Order)) df$Order[df$Order == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha"))
    
    id_obs = which(df$Method == 'Observed')
    
    for (i in 1:length(id_obs)) {
      
      new = df[id_obs[i],]
      new$Size = new$Size - 0.0001
      new$Method = 'Interpolated'
      
      newe = df[id_obs[i],]
      newe$Size = newe$Size + 0.0001
      newe$Method = 'Extrapolated'
      
      df = rbind(df, new, newe)
      
    }
    
    ylab = "Taxonomic diversity"
    
    lty = c(Interpolated = "solid", Extrapolated = "dashed")
    df$Method = factor(df$Method, levels = c('Interpolated', 'Extrapolated', 'Observed'))
    
    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolated" & round(Size) %in% double_size)
    
    ggplot(data = df, aes(x = Size, y = Estimate, col = Region)) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) + 
      geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) + scale_linetype_manual(values = lty) +
      # geom_line(lty=2) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5)+
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) + 
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) + 
      scale_colour_manual(values = cbPalette) + 
      scale_fill_manual(values = cbPalette) + 
      facet_grid(div_type ~ Order, scales = scale) +
      theme_bw() + 
      theme(legend.position = "bottom", legend.title = element_blank()) +
      labs(x = xlab, y = ylab, title = main)
  }
  
}


coverage_to_size = function(x, C, datatype = 'abundance'){
  
  if (datatype == 'abundance'){
    
    n <- sum(x)
    refC <- iNEXT.3D:::Coverage(x, 'abundance', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'abundance', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  } else {
    
    m <- NULL
    n <- max(x)
    refC <- iNEXT.3D:::Coverage(x, 'incidence_freq', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'incidence_freq', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  }
  
  return(mm)
}

bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
  
  if (datatype == 'abundance'){
    
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
    
    output = apply(data, 2, function(x){
      
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      
      if(length(p_i_hat) != length(x)){
        
        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)
        
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
        
        p_i_hat
        
      } else {
        
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
        
      }
    })
    
  }
  
  if (datatype == 'incidence'){
    
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
    
    output = apply(data, 2, function(x){
      
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      
      if(length(pi_i_hat) != (length(x) - 1)){
        
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs
        
        pi_i_hat
        
      } else {
        
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
        
      }
    })
    
  }
  
  return(output)
  
}

