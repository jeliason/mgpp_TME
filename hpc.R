library(ppjsdm)
`%>%` <- magrittr::`%>%`
library(dplyr)
library(spatstat)

ROOT_PATH = here::here()

DATA_PATH <- paste0(ROOT_PATH,"/data/")

df_raw <- readr::read_csv(paste0(DATA_PATH,"CRC_cleaned.csv")) %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)

#### FIT ALL MODELS
region_encoding <- function(df) {
  window_ss <- owin(c(min(df$X),max(df$X)),c(min(df$Y),max(df$Y)))
  pat <- ppp(df$X,df$Y,marks = df$type,window = window_ss)

  cvx_window <- convexhull(pat)
  dens <- density(unmark(pat),sigma = bw.diggle)
  dens_df <- dens %>%
    as.data.frame() %>%
    mutate(value = inside.owin(x,y,cvx_window)) %>%
    mutate(value = ifelse(value,"inside","outside"))

  cvx_pat <- pat
  Window(cvx_pat) <- cvx_window
  cvx_dens <- density(unmark(cvx_pat),sigma = bw.diggle)

  cvx_dens_df <- (cvx_dens > intensity(unmark(cvx_pat))) %>%
    as.data.frame() %>%
    mutate(value = ifelse(value,"high_dens","low_dens"))

  dens_df <- dens_df %>%
    left_join(cvx_dens_df,by=c("x","y")) %>%
    mutate(value = ifelse(is.na(value.y),value.x,value.y)) %>%
    select(-value.x,-value.y)


  hot_enc <- fastDummies::dummy_cols(dens_df,remove_first_dummy = TRUE,select_columns = "value",remove_selected_columns = TRUE) %>%
    as.data.frame()

  low_dens_var <- hot_enc %>%
    select(x,y,value_low_dens) %>%
    rename(value = value_low_dens) %>%
    as.im()

  outside_var <- hot_enc %>%
    select(x,y,value_outside) %>%
    rename(value = value_outside) %>%
    as.im()

  list(outside_var=outside_var,low_dens_var=low_dens_var)
}

df_raw %>%
  dplyr::group_by(Spot) %>%
  dplyr::group_map(~{
    df = .x %>%
      dplyr::mutate(type = forcats::fct_lump_min(type,20)) %>%
      # mutate(type = fct_lump_n(type,)) %>%
      dplyr::filter(type != "Other") %>%
      droplevels()

    config <- Configuration(df$X,df$Y,df$type)
    window <- Rectangle_window(x_range = c(min(df$X),max(df$X)), y_range = c(min(df$Y),max(df$Y)))

    region_var <- region_encoding(.x)

    list(config=config,window=window,outside_var=region_var$outside_var,spot=.y$Spot)
  }) -> cw


fit_configs <- function(obj) {

  get_papangelous <- function(fit,nthreads=1) {
    config <- fit$configuration_list[[1]]
    window <- fit$window
    types <- levels(config$types)

    papangelous <- setNames(lapply(types, function(sp) {
      print(sp)
      plot_papangelou(fit,
                      window = window,
                      configuration = config,
                      type = sp,
                      use_log = FALSE,
                      drop_type_from_configuration = TRUE,
                      return_papangelou = TRUE,
                      grid_steps = c(200, 200),
                      nthreads=nthreads)
    }), nm = types)
  }

  auc_papangelous <- function(papangelous,fit) {
    config <- fit$configuration_list[[1]]
    window <- fit$window
    types <- names(papangelous)
    aucs <- sapply(types, function(sp) {
      intensity <- papangelous[[sp]]
      conf <- ppp(x = config$x[config$types == sp],
                  y = config$y[config$types == sp],
                  window = owin(window$x_range, window$y_range))
      spatstat.explore::auc(X = conf, covariate = intensity)
    })
  }

  config <- obj$config
  window <- obj$window
  outside_var <- obj$outside_var
  spot <- obj$spot


  number_of_species <- length(levels(config$types))
  short_range <- matrix(30, number_of_species, number_of_species)
  medium_range <- matrix(70, number_of_species, number_of_species)
  long_range <- matrix(150, number_of_species, number_of_species)
  saturation <- 5

  covariates <- list()
  covariates$outside = outside_var
  # need to add in nndist-vasculature here as well

  args = list(configuration = config,
              window = window,
              covariates = covariates,
              model = "exponential",
              medium_range_model = "exponential",
              short_range = short_range,
              medium_range = medium_range,
              long_range = long_range,
              fitting_package = "glmnet",
              dummy_distribution = "stratified",
              saturation = saturation,
              min_dummy = 5e3,
              max_dummy = 5e3,
              debug = TRUE,
              nthreads = 2L)
  tm <- Sys.time()
  fit <- do.call(gibbsm,args)
  print(Sys.time() - tm)

  s <- tryCatch({
    summary(fit,nthreads = 2L)
  },
  error=\(e) NULL)

  aucs <- tryCatch({
    papangelous <- get_papangelous(fit,nthreads = 2)
    auc_papangelous(papangelous,fit)
  },
  error=\(e) NULL)

  list(summary=s,aucs=aucs,spot=spot)
}

rslurm::slurm_map(cw,fit_configs,nodes=length(cw),cpus_per_node = 1,slurm_options=list(time='30:00',mem='10g'),submit=FALSE)
