load_pt_data <- function(DATA_PATH = paste0(here::here(),"/data/")) {
  pt_data <<- readr::read_csv(paste0(DATA_PATH,"CRC_pt_metadata.csv"))
}

data_load_CRC <- function() {
  ROOT_PATH <<- here::here()
  DATA_PATH <<- paste0(ROOT_PATH,"/data/")

  df_raw1 <- readr::read_csv(paste0(DATA_PATH,"first_half_CRC_cleaned.csv"))

  df_raw2 <- readr::read_csv(paste0(DATA_PATH,"second_half_CRC_cleaned.csv"))

  df_raw <<- bind_rows(df_raw1, df_raw2) %>%
    dplyr::mutate(type = as.factor(type))

  load_pt_data(DATA_PATH = DATA_PATH)

}

nndist_cell_type <- function(df,type) {
  window_ss <- owin(c(min(df$X),max(df$X)),c(min(df$Y),max(df$Y)))
  pat <- ppp(df$X,df$Y,marks = df$type,window = window_ss)

  dens <- density(unmark(pat),sigma = bw.diggle) %>%
    as.data.frame()

  pat_dens <- ppp(dens$x,dens$y,window = window_ss)
  nnd <- nncross(pat_dens,subset(pat,marks==type),what="dist")

  dens$value <- nnd

  as.im(dens)
}

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

plot_coef <- function(sum_fit,coef="alpha",return_df=FALSE) {
  nms <- rownames(sum_fit$se$alpha[[1]])

  nms_tib <- tibble(type=nms,num_type=1:length(nms))
  sum_fit$coefficients %>%
    as_tibble(rownames = "coef_name") %>%
    filter(stringr::str_detect(coef_name,coef)) %>%
    separate(coef_name,into = c(coef,"type1","type2"),sep="_") %>%
    select(-!!sym(coef)) %>%
    mutate(type1=as.numeric(type1),
           type2=as.numeric(type2)) %>%
    left_join(nms_tib,by=c("type1"="num_type")) %>%
    select(-type1) %>%
    rename(type1=type) %>%
    left_join(nms_tib,by=c("type2"="num_type")) %>%
    select(-type2) %>%
    rename(type2=type) -> df

  df %>%
    ggplot(aes(type1,type2,fill=coefficients)) +
    geom_tile() +
    anglex() +
    sig_stars(p_values="Pval") +
    scale_fill_gradient2() -> p
  if(return_df) {
    return(list(p=p,df=df))
  }
  else {
    return(p)
  }
}

anglex <- function() {
  list(ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1,vjust=1)))
}


sig_stars <- function(alpha=0.05,p_values="p.adj") {
  list(
    ggplot2::geom_point(ggplot2::aes(shape=ifelse(!!sym(p_values) < alpha, "dot", "no_dot"),
                                     alpha=ifelse(!!sym(p_values) < alpha, "dot", "no_dot")),
                        color='red',size=2),
    ggplot2::scale_shape_manual(values=c(dot=8, no_dot=0), guide="none"),
    # Rather than setting the shape to be NA, which throws a warning, I
    # set the alpha based on significance to make non-significant stars
    # invisible
    ggplot2::scale_alpha_manual(values=c(dot=1, no_dot=0), guide="none")
  )
}

expand.grid.unique <- function(x, y, include.equals=TRUE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

plot_all_potentials <- function(fit,types=fit$type_names) {
  combs <- expand.grid.unique(types,types) %>%
    as.data.frame()

  names(combs) <- c("type1","type2")

  pots <- combs %>%
    purrr::pmap(\(type1,type2) {
      df <- plot_potentials(fit,type1,type2,return_df = T)
      df$type1 = type1
      df$type2 = type2
      df
    }) %>%
    bind_rows()

  pots %>%
    ggplot2::ggplot(ggplot2::aes(x,value,color=potential,linetype=dashed)) +
    ggplot2::geom_line(linewidth=1) +
    ggplot2::guides(linetype="none") +
    ggplot2::facet_grid(type1~type2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,hjust=1,vjust=1),
                   strip.text = ggplot2::element_text(size = 7)) -> p

  # print(p)
  p
  # return(pots)
}

plot_potentials <- function(fit,type1,type2,return_df=FALSE) {
  pot <- ppjsdm::potentials(fit,type1,type2)
  lr <- max(pot$long_range[[1]],pot$short_range[[1]][[1]])
  xrange <- seq(0.01,lr,1)
  df <- tibble::tibble(x=xrange,
                       short=pot$short[[1]][[1]](xrange),
                       medium=pot$medium[[1]](xrange),
                       overall=pot$overall[[1]](xrange)) %>%
    tidyr::pivot_longer(-x) %>%
    dplyr::mutate(dashed=ifelse(name=="overall",FALSE,TRUE)) %>%
    dplyr::rename(potential=name)
  if(return_df)
    return(df)

  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x,value,color=potential,linetype=dashed)) +
    ggplot2::geom_line(linewidth=1) +
    ggplot2::guides(linetype="none") +
    ggplot2::ggtitle(paste0("Interaction potential between ", type1, " and ", type2))

  p
  # print(p)
}

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
                    drop_type_from_configuration = FALSE,
                    return_papangelou = TRUE,
                    grid_steps = c(200, 200),
                    nthreads=nthreads)
  }), nm = types)
}

auc_papangelous <- function(papangelous,fit=NULL,config=NULL,window=NULL) {
  if(!is.null(fit) & is.null(config) & is.null(window)) {
    config <- fit$configuration_list[[1]]
    window <- fit$window
  }
  types <- names(papangelous)
  aucs <- sapply(types, function(sp) {
    intensity <- papangelous[[sp]]
    conf <- ppp(x = config$x[config$types == sp],
                y = config$y[config$types == sp],
                window = owin(window$x_range, window$y_range))
    spatstat.explore::auc(X = conf, covariate = intensity)
  }) %>%
    enframe() %>%
    rename(type=name,AUC=value)
}

plot_papangelous <- function(papangelous,fit,types) {
  config <- fit$configuration_list[[1]]
  aucs <- auc_papangelous(papangelous,fit)
  lapply(types,\(type) {
    df_papa <- papangelous[[type]] %>%
      as.data.frame()

    df_ppp <- tibble(x = config$x[config$types == type],
                     y = config$y[config$types == type])

    ggplot() +
      geom_tile(data=df_papa,aes(x,y,fill=log(value))) +
      geom_point(data=df_ppp,aes(x,y)) +
      ggtitle(paste0("Log-conditional intensity for ",type, ", AUC: ", round(aucs$AUC[aucs$type == type],3))) +
      scico::scale_fill_scico(palette = "tokyo",direction = 1) -> p
  }) -> p
}

coefs_from_summary <- function(summary_fit,types) {
  coef_num <- 1:length(types)
  tb_match <- tibble(coef_num=as.character(coef_num),type=types)
  summary_fit$coefficients %>%
    as_tibble(rownames = "coef_name") %>%
    select(coefficients,coef_name) %>%
    filter(str_detect(coef_name,"alpha") | str_detect(coef_name,"gamma")) %>%
    separate(coef_name,c("name","type1","type2"),sep="_") %>%
    drop_na() %>%
    left_join(tb_match,by=c("type1"="coef_num")) %>%
    select(-type1) %>%
    rename(type1=type) %>%
    left_join(tb_match,by=c("type2"="coef_num")) %>%
    select(-type2) %>%
    rename(type2=type)
}
