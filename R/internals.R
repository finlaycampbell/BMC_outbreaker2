##===== Setting up =====#

## Load all libraries passed as arguments
libs <- function(...) {

  invisible(lapply(unlist(list(...)), require, character.only = TRUE))

}

## Connect to cluster
set.up <- function(clust = 'mrc') {

  setwd("/home/fc1915/mnt/fc1915/ob2")

  options(didehpc.cluster = "fi--dideclusthn",
          didehpc.credentials = "~/.smbcredentials")

  if(clust == "mrc") options(didehpc.cluster = "fi--didemrchnb")

  didehpc::didehpc_config_global(temp = didehpc::path_mapping('tmp',
                                                              '/home/fc1915/mnt/tmp',
                                                              '//fi--didef3.dide.ic.ac.uk/tmp',
                                                              'T:'))

  didehpc::web_login()

  ## Remove outbreaker2 if only using phybreak
  our.pkgs <- c('ape', 'outbreaker2', 'TransPhylo', 'o2mod.transphylo',
                'distcrete', 'phybreak', 'tidyr', 'dplyr')

  pkg <- provisionr::package_sources(local = c("~/mnt/fc1915/ob2/outbreaker2_1.0-0.tar.gz",
                                               "~/mnt/fc1915/ob2/o2mod.transphylo_0.0.0.9000.tar.gz"),
                                     github = c("xavierdidelot/TransPhylo", "reconhub/distcrete"))

  our.sources <- "~/mnt/fc1915/ob2/functions.R"

  ctx <- context::context_save("contexts_2", packages = our.pkgs,
                               sources = our.sources,
                               package_sources = pkg)

  obj <- didehpc::queue_didehpc(ctx)

  return(obj)
}

## Cluster on windows machine
set.win.up <- function(clust = 'mrc') {

  setwd("Q:\\ob2")

  options(didehpc.cluster = "fi--dideclusthn",
          didehpc.credentials = "~\\.smbcredentials",
          didehpc.home = "q:\\ob2")

  if(clust == "mrc") options(didehpc.cluster = "fi--didemrchnb")

  #didehpc::didehpc_config_global(temp = didehpc::path_mapping('tmp',
  #                                                            '/home/fc1915/mnt/tmp',
  #                                                            '//fi--didef3.dide.ic.ac.uk/tmp',
                                        #                                                            'T:'))

  didehpc::web_login()

  ## Remove outbreaker2 if only using phybreak
  our.pkgs <- c('ape', 'outbreaker2', 'TransPhylo', 'o2mod.transphylo',
                'distcrete', 'phybreak', 'tidyr', 'dplyr')

  pkg <- provisionr::package_sources(local = c("Q:\\ob2/outbreaker2_1.0-0.tar.gz",
                                               "Q:\\ob2/o2mod.transphylo_0.0.0.9000.tar.gz"),
                                     github = c("xavierdidelot/TransPhylo", "reconhub/distcrete"))

  our.sources <- "Q:\\ob2\\functions.R"

  ctx <- context::context_save("contexts_3", packages = our.pkgs,
                               sources = our.sources,
                               package_sources = pkg)

  obj <- didehpc::queue_didehpc(ctx)

  return(obj)

}

##===== Analysis =====##

## Simulate a phybreak outbreak
run.sim <- function(size) {

  R0 <- 1.8
  w.mean <- 14.4
  w.sd <- 8.9
  mu <- 3.10e-6
  seql <- 18958

  shape.gen <- w.mean^2/w.sd^2
  scale.gen <- w.mean/shape.gen

  si <- distcrete("gamma", shape = shape.gen, scale = scale.gen,
                  interval = 1L, w = 0)
  w <- si$d(1:(w.mean*10))


  sim <- sim.phybreak(obsize = size,
                      R0 = R0,
                      mean.gen = w.mean,
                      shape.gen = shape.gen,
                      mean.sample = w.mean,
                      shape.sample = shape.gen,
                      sequence.length = seql,
                      mu = mu,
                      wh.model = 3,
                      wh.slope = 1)
  sim$w <- w
  sim$w.mean <- w.mean
  sim$shape <- shape.gen
  sim$scale <- scale.gen

  return(sim)

}

## Runs o2mod.transphylo on a phybreak sim
run.o2mod <- function(sim) {

  dates <- sim$sample.times
  phy <- sim$sim.tree
  nsam <- length(dates)
  phy$tip.label <- 1:nsam
  data <- outbreaker_data(dates = dates, w_dens = sim$w,
                          dna = as.DNAbin(matrix('A', nsam, 1)))
  data$ptree <- ptreeFromPhylo(phy, max(dates))
  data$neg <- sim$w.mean
  config <- list(n_iter = 1e3,
                 sample_every = 5)
  res <- o2mod.transphylo(data, config)

}

## Runs transphylo on a phybreak sim
run.trans <- function(sim) {

  times <- sim$sample.times*50
  phy <- sim$sim.tree
  phy$edge.length <- phy$edge.length*50
  phy$root.edge <- phy$root.edge*50
  phy$tip.label <- 1:length(times)

  ptree <- ptreeFromPhylo(phy, max(sim$sample.times))

  # Everything needs to be in unit days, not years
  # Neg = 0.25 is in years (from paper ~100 days)
  # Therefore multiply by 365
  # Our shape and scale are already in days
  # Scale up by 50 to improve mixing

  res <- inferTTree(ptree, w.shape = sim$shape, w.scale = sim$scale*50,
                    startPi = 1, updatePi = F, startNeg = sim$w.mean*50,
                    updateNeg = F, startOff.r = 1.8, updateOff.r = F,
                    startOff.p = 0.5, updateOff.p = F, thinning = 50,
                    mcmcIterations = 1e4, dateT = Inf)

  return(res)

}

## Runs o2 on a phybreak sim
run.o2 <- function(sim) {

  dates <- sim$sample.times
  dna <- as.DNAbin(sim$sequences)
  rownames(dna) <- seq_len(nrow(dna))
  data <- outbreaker_data(dates = dates, w_dens = sim$w, dna = dna)
  config <- list(n_iter = 1e3,
                 sample_every = 5)
  res <- outbreaker2::outbreaker(data, config)

}

## Returns the accuracy of a outbreak run on a phybreak sim
get.o2mod.acc <- function(sim, res, burnin = 0.1) {

  cons <- get.o2mod.cons(res, burnin = burnin)

  ind <- cons > 0
  cons[ind] <- paste0('host.', cons[ind])
  cons[!ind] <- 'index'

  acc <- mean(cons == sim$sim.infectors)

  return(acc)

}

## Returns the accuracy of a transphylo run on a phybreak sim
get.trans.acc <- function(sim, res, burnin = 0.1) {

  cons <- get.trans.cons(res, burnin = burnin)

  ind <- cons > 0
  cons[ind] <- paste0('host.', cons[ind])
  cons[!ind] <- 'index'

  acc <- mean(cons == sim$sim.infectors)

  return(acc)

}

## Get the consensus ancestors
get.o2mod.cons <- function(res, burnin = 0.1) {

  cons <- summary(res, burnin = burnin*length(res$mu))$tree$from
  cons[is.na(cons)] <- 0
  return(cons)

}

## Get the consensus ancestors
get.trans.cons <- function(res, burnin = 0.1) {

  mat <- computeMatWIW(res, burnin)
  # this counts the 'index' assigments
  ind <- 1 - apply(mat, 2, sum)
  mat <- rbind(ind, mat)

  cons <- apply(mat, 2, which.max) - 1

}

## Returns all ancestries with their posterior probability
get.o2mod.anc <- function(res, burnin) {

    res <- res[res$step > max(res$step)*burnin,,drop = FALSE]

    alpha <- as.matrix(res[,grep("alpha", names(res))])
    colnames(alpha) <- seq_len(ncol(alpha))
    from <- as.vector(alpha)
    to <- as.vector(col(alpha))
    from[is.na(from)] <- 0
    out <- data.frame(xyTable(from,to))
    names(out) <- c("from", "to", "frequency")

    ## Calculate proportion among ancestries
    get.prop <- function(i) {
      ind <- which(out$to == out$to[i])
      out[[3]][i]/sum(out[[3]][ind])
    }

    out[3] <- vapply(seq_along(out[[3]]), get.prop, 1)
    out$mod <- 'o2mod'

    return(out)

}

## Returns all ancestries with their posterior probability
get.trans.anc <- function(res, burnin) {

  mat <- computeMatWIW(res, burnin)
  ind <- 1 - apply(mat, 2, sum)
  mat <- rbind(ind, mat)

  w <- which(mat > 0, arr.ind = T, useNames = F)
  from <- w[,1] - 1
  to <- w[,2]

  out <- data.frame(xyTable(from, to))
  out$number <- NULL
  names(out) <- c('from', 'to')

  mat <- t(mat)
  out$frequency <- mat[which(mat>0)]
  out$mod <- 'trans'

  return(out)

}

## How similar are o2mod and trans in their consensus trees
get.simil <- function(o2mod.res, trans.res) {

  mean(get.o2mod.cons(o2mod.res) == get.trans.cons(trans.res))

}

## How similar are o2 and trans in their consensus trees
get.o2.simil <- function(res1, res2) {

  mean(get.o2mod.cons(res1) == get.o2mod.cons(res2))

}

## Runs o2.mod / transphylo on the cluster
run.analysis <- function(runs, size) {

  store <- list()
  store$sim <- store$o2mod.res <- store$trans.res <- store$o2.res <- list()
  summ <- data.frame(o2mod = numeric(runs), trans = numeric(runs), o2 = numeric(runs))
  times <- data.frame(o2mod = numeric(runs), trans = numeric(runs), o2 = numeric(runs))

  # Wrapper function to allow random errors to be discarded and try again
  attempt <- function() {

    sim <- run.sim(size = size)
    o2mod.time <- system.time(o2mod.res <- run.o2mod(sim))[3]
    trans.time <- system.time(trans.res <- run.trans(sim))[3]
    o2.time <- system.time(o2.res <- run.o2(sim))[3]

    out <- list(sim = sim, o2mod.res = o2mod.res, trans.res = trans.res, o2.res = o2.res,
                o2mod.time = o2mod.time, trans.time = trans.time, o2.time = o2.time)

    return(out)

  }

  for(i in 1:runs) {

    print(i)

    while(TRUE){
       tmp <- try(attempt(), silent=FALSE)
       if(!is(tmp, 'try-error')) break
    }

    summ$o2mod[i] <- get.o2mod.acc(tmp$sim, tmp$o2mod.res)
    summ$trans[i] <- get.trans.acc(tmp$sim, tmp$trans.res)
    summ$o2[i] <- get.o2mod.acc(tmp$sim, tmp$o2.res)

    times$o2mod[i] <- tmp$o2mod.time
    times$trans[i] <- tmp$trans.time
    times$o2[i] <- tmp$o2.time

    store$sim[[i]] <- tmp$sim
    store$o2mod.res[[i]] <- tmp$o2mod.res
    store$trans.res[[i]] <- tmp$trans.res
    store$o2.res[[i]] <- tmp$o2.res

  }

  summ$diff <- summ$o2mod - summ$trans
  store$summ <- summ
  store$times <- times

  return(store)

}


##===== Result collection =====#

## Add loaded objects to store
create.store <- function(obj, bundle.name, dir, load, dl, store = NULL) {

  files <- list.files(dir)
  rem <- grep("store", files)
  if(length(rem) > 0) files <- files[-grep("store", files)]
  n.files <- length(files)

  if(load) {
    pb <- txtProgressBar(min = 1, max = length(files), style = 3)
    for(i in seq_along(files)) {
      setTxtProgressBar(pb, i)
      load(paste0(dir, files[i]))
      store <- mk.summary(store, r)
    }
  }

  if(dl) {
    task_bundle <- obj$task_bundle_get(bundle.name)
    ids <- task_bundle$ids
    ids <- ids[task_bundle$status() == "COMPLETE"]
    pb <- txtProgressBar(min = 1, max = length(ids), style = 3)
    for(i in seq_along(ids)) {
      setTxtProgressBar(pb, i)
      task <- obj$task_get(ids[i])
      r <- task$result()
      save(r, file = paste0(dir, "r.", n.files + i, ".RData"))
      store <- mk.summary(store, r)
    }
  }

  return(store)

}

## Add one run object to store
mk.summary <- function(store = NULL, r) {

  if(is.null(store)) {

    store <- list(acc = data.frame(o2mod = numeric(),
                                   trans = numeric(),
                                   o2 = numeric()),
                  simil = data.frame(o2mod = numeric(),
                                     o2 = numeric()),
                  times = data.frame(o2mod = numeric(),
                                     trans = numeric(),
                                     o2 = numeric()))

  }
  
  simil <- data.frame(o2mod = sapply(seq_along(r$o2mod.res),
                                     function(i) get.simil(r$o2mod.res[[i]],
                                                           r$trans.res[[i]])),
                      o2 = sapply(seq_along(r$o2mod.res),
                                  function(i) get.simil(r$o2.res[[i]],
                                                        r$trans.res[[i]])))

  store$acc <- rbind(store$acc, r$summ)
  store$times <- rbind(store$times, r$times)
  store$simil <- rbind(store$simil, simil)

  return(store)

}

## Load the i'th file in a directory
load.ifile <- function(dir, i = 1) {

  load(list.files(dir, full.names = T)[i], .GlobalEnv)

}


##===== Plots =====#

## Provides a label for facetting
create.flab <- function() {
  as_labeller(c(o2mod = 'o2mod.TransPhylo',
                trans = 'TransPhylo'))
}

## Axis labels
create.axlab <- function(a) {

  def <- list(acc = 'Accuracy of outbreak reconstruction',
              simil = 'Proportion of identical ancestry assignments',
              diff = 'Difference in accuracy of outbreak reconstruction')
  def[[a]]

}

## Uniform theme
create.theme <- function() {
  theme_light(base_size = 18)
}

## Visualise density plot of accuracies
vis.dens <- function(store) {

  ## Convert to long format
  store$acc %>%
    gather(mod, acc, -diff, -simil) %>%
    ggplot(aes(mod, acc, fill = mod)) +
    geom_violin() +
    viridis::scale_fill_viridis(discrete = TRUE) +
    scale_x_discrete(labels = create.flab()) +
    xlab(NULL) +
    ylab('Accuracy of outbreak reconstruction') +
    theme(legend.position = 'none')

}

## Compare o2mod and trans accuracies
vis.rel <- function(store) {

  df <- store$acc %>%
    arrange(trans, diff) %>%
    mutate(i = seq_along(diff))
           #o2mod = o2mod - trans,
           #trans = 0)

  df2 <- gather(df, mod, acc, -diff, -i, -simil)

  ggplot(df) +
    geom_segment(aes(x = i, y = trans, xend = i, yend = o2mod, colour = diff), size = 1) +
    geom_point(data = df2, aes(i, acc, fill = mod), shape = 21, size = 1.5) +
    create.theme() +
    scale_color_viridis(guide = 'none') +
    scale_fill_manual(values = c('white', 'black'),
                    labels = create.flab(),
                    name = NULL) +
    labs(y = 'Accuracy of outbreak reconstruction',
         x = 'Run') +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal') +
    ylim(0, 1)

}

## Compare o2mod and trans ancestries
vis.anc <- function(o2mod.res, trans.res, burnin = 0.1, nbreaks = 10) {

  df <- rbind(get.o2mod.anc(o2mod.res, burnin),
              get.trans.anc(trans.res, burnin))

  ggplot(df) +
    geom_point(aes(x = to, y = from, size = frequency, color = factor(from))) +
    facet_grid(mod ~ ., labeller = create.flab()) +
    scale_size_area(name = 'Posterior\nfrequency') +
    viridis::scale_color_viridis(discrete = TRUE) +
    guides(colour = FALSE) +
    labs(x = 'Infectee', y = 'Infector') +
    scale_x_continuous(breaks = seq(1, max(df$to),
                                    round(max(df$to)/nbreaks, 0))) +
    scale_y_continuous(breaks = seq(0, max(df$from),
                                    round(max(df$from)/nbreaks, 0))) +
    create.theme() +
    theme(strip.text.y = element_text(colour = "black"),
          strip.background = element_rect(colour = "darkgrey", fill = "grey90"),
          panel.grid.minor = element_blank())

}

## Histogram of similarity in ancestry assignments
vis.hist <- function(store, a, bw = 0.07) {

  ggplot(store$summ, aes_string(x = a)) +
    geom_histogram(binwidth = bw) +
    create.theme() +
    labs(x = create.axlab(a),
         y = 'Count')

}

## Compare o2mod and trans chains
vis.chains <- function(o2mod.res, trans.res, burnin = 1, burnend = length(o2mod.res$like)) {

  to.keep <- burnin:burnend

  o2mod.chain <- o2mod.res$like[to.keep]
  trans.chain <- sapply(trans.res[to.keep],
                        function(x) x$pTTree+x$pPTree)

  o2mod <- data.frame(iter = to.keep,
                      like = o2mod.chain,
                      mod = 'o2mod')

  trans = data.frame(iter = to.keep,
                     like = trans.chain,
                     mod = 'trans')

  df <- rbind(o2mod, trans)

  ggplot(df, aes(iter, like, group = mod)) +
    geom_line() +
    facet_grid(mod ~  ., scale = 'free', labeller = create.flab()) +
    create.theme() +
    labs(y = 'Posterior log likelihood', x = 'MCMC iteration') +
    theme(legend.position = 'none',
          strip.text.y = element_text(colour = "black"),
          strip.background = element_rect(colour = "darkgrey", fill = "grey90"))

}

## Save plots to ob2/figs directory
vis.save <- function(p, name, ext = 'svg', dpi = 500,
                     width = 10, height = 10, ...) {

  ggsave(paste0("~/Dropbox/phd/ob2/figs/", name, ".", ext),
         p, dpi = dpi, width = width, height = height, ...)

}

## Update the plots in the final submission dired
vis.update <- function() {

  p1 <- vis.chains(chains$o2mod.res[[1]], chains$trans.res[[1]])
  vis.save(p1, 'chains', 'png', width = 10, height = 10)

  p2 <- vis.anc(ances$o2mod.res[[1]], ances$trans.res[[1]])
  vis.save(p2, 'ances', 'png', width = 10, height = 10)

  p3 <- vis.rel(store)
  vis.save(p3, 'rel', 'png', width = 15, height = 10, dpi = 200)

  #3 <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  #is.save(p3, 'join', 'png', width = 20, height = 10, dpi = 500)

}
