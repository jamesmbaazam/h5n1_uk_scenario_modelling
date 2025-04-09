# script to plot the outbreak size (and length) distribution for H5N1 scenarios

library(epichains)
library(ggplot2)
library(cowplot)

# Simulating an outbreak size distribution

statistic <- "size"
offspring_dist <- "rpois"
R <- seq(0.1, 1.1, 0.1)

# parameter space
scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  stringsAsFactors = FALSE
)

n_chains <- 1e5

breaks <- c(0, 1, 2, 5, 10, 20, 50, Inf)

outbreak_list <- vector(mode = "list", length = nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains, 
    statistic = scenarios[i, "statistic"], 
    offspring_dist = offspring_dist_fun,
    lambda = scenarios[i, "R"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# remove index case
outbreak_list <- lapply(outbreak_list, function(x) x - 1)

intervals <- lapply(
  outbreak_list, 
  cut, 
  breaks = breaks, 
  right = FALSE, 
  include.lowest = TRUE
)
prop <- lapply(intervals, function(interval) table(interval) / sum(table(interval)))
outbreak_size_list <- lapply(prop, as.data.frame)
for (i in seq_len(nrow(scenarios))) {
  outbreak_size_list[[i]]$R <- scenarios[i, "R"]
  outbreak_size_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_size_list[[i]]$statistic <- scenarios[i, "statistic"]
}
outbreak_size <- do.call(rbind, outbreak_size_list)
head(outbreak_size)

pois_size <- ggplot2::ggplot(data = outbreak_size) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(
    name = "Reproduction number (R)",
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak size\n(secondary cases)", 
    palette = "Spectral"
  ) + 
  ggplot2::labs(title = "Poisson distribution") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = element_text(hjust = 0.5))

# summary statistics from outbreak size

# proportion of outbreaks that do not have any secondary transmission
no_secondary <- outbreak_size[outbreak_size$interval == "[0,1)", ]

# R values where most outbreaks (>50%) do not cause any secondary transmission
no_secondary[no_secondary$Freq > 0.5, "R"]
# so the maximum value is the max R where the majority of outbreaks are only index cases

# proportion of outbreaks that contain >=20 cases
gt_20_secondary <- outbreak_size[outbreak_size$interval %in% c("[20,50)", "[50,Inf]"), ]

gt_20_secondary |> dplyr::summarise(Freq = sum(Freq), .by = "R")

# Simulating an outbreak size distribution with a Negative binomial distribution
statistic <- "size"
offspring_dist <- "rnbinom"
R <- seq(0.1, 1.1, 0.1)
k <- c(0.1, 0.5, 5, 1000)

scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  k = k,
  stringsAsFactors = FALSE
)

outbreak_list <- vector(mode = "list", length = nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains, 
    statistic = scenarios[i, "statistic"], 
    offspring_dist = offspring_dist_fun,
    mu = scenarios[i, "R"],
    size = scenarios[i, "k"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# remove index case
outbreak_list <- lapply(outbreak_list, function(x) x - 1)

intervals <- lapply(
  outbreak_list, 
  cut, 
  breaks = breaks, 
  right = FALSE, 
  include.lowest = TRUE
)
prop <- lapply(intervals, function(interval) table(interval) / sum(table(interval)))
outbreak_size_list <- lapply(prop, as.data.frame)
for (i in seq_len(nrow(scenarios))) {
  outbreak_size_list[[i]]$R <- scenarios[i, "R"]
  outbreak_size_list[[i]]$k <- scenarios[i, "k"]
  outbreak_size_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_size_list[[i]]$statistic <- scenarios[i, "statistic"]
}
outbreak_size <- do.call(rbind, outbreak_size_list)
head(outbreak_size)

nbinom_size <- ggplot2::ggplot(data = outbreak_size) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(
    name = "Reproduction number (R)",
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak size\n(secondary cases)", 
    palette = "Spectral"
  ) + 
  ggplot2::facet_wrap(
    facets = c("k"), 
    labeller = ggplot2::label_both
  ) +
  ggplot2::labs(title = "Negative Binomial distribution") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  )


# summary statistics from outbreak size

# proportion of outbreaks that cause >= 50 secondary cases when k = 0.5
outbreak_size[
  outbreak_size$k == 0.5 &
    outbreak_size$interval %in% "[50,Inf]", 
]
# less than 10% of outbreaks for R = 1.1 and k = 0.5 exceed 50 secondary cases

size <- cowplot::plot_grid(
  pois_size + ggplot2::theme(legend.position = "none"), 
  nbinom_size, 
  nrow = 1,
  rel_widths = c(0.5, 1),
  labels = c("A", "B")
)

ggplot2::ggsave(
  filename = file.path("plots", "outbreak_size.png"), 
  plot = size,
  device = "png", 
  width = 250, 
  height = 150,
  units = "mm",
  dpi = 300
)

# Simulating an outbreak length distribution

statistic <- "length"
offspring_dist <- "rpois"
scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  stringsAsFactors = FALSE
)

outbreak_list <- vector(mode = "list", length = nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains, 
    statistic = scenarios[i, "statistic"], 
    offspring_dist = offspring_dist_fun,
    lambda = scenarios[i, "R"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# remove index case
outbreak_list <- lapply(outbreak_list, function(x) x - 1)

intervals <- lapply(
  outbreak_list, 
  cut, 
  breaks = breaks, 
  right = FALSE, 
  include.lowest = TRUE
)
prop <- lapply(intervals, function(interval) table(interval) / sum(table(interval)))
outbreak_length_list <- lapply(prop, as.data.frame)
for (i in seq_len(nrow(scenarios))) {
  outbreak_length_list[[i]]$R <- scenarios[i, "R"]
  outbreak_length_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_length_list[[i]]$statistic <- scenarios[i, "statistic"]
}
outbreak_length <- do.call(rbind, outbreak_length_list)
head(outbreak_length)

pois_length <- ggplot2::ggplot(data = outbreak_length) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(
    name = "Reproduction number (R)",
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak length\n(Secondary generations)", 
    palette = "Spectral"
  ) + 
  ggplot2::labs(title = "Poisson distribution") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = element_text(hjust = 0.5))


# Simulating an outbreak length distribution with a Negative binomial distribution
statistic <- "length"
offspring_dist <- "rnbinom"
R <- seq(0.1, 1.1, 0.1)
k <- c(0.1, 0.5, 5, 1000)

scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  k = k,
  stringsAsFactors = FALSE
)

outbreak_list <- vector(mode = "list", length = nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains, 
    statistic = scenarios[i, "statistic"], 
    offspring_dist = offspring_dist_fun,
    mu = scenarios[i, "R"],
    size = scenarios[i, "k"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# remove index case
outbreak_list <- lapply(outbreak_list, function(x) x - 1)

intervals <- lapply(
  outbreak_list, 
  cut, 
  breaks = breaks, 
  right = FALSE, 
  include.lowest = TRUE
)
prop <- lapply(intervals, function(interval) table(interval) / sum(table(interval)))
outbreak_length_list <- lapply(prop, as.data.frame)
for (i in seq_len(nrow(scenarios))) {
  outbreak_length_list[[i]]$R <- scenarios[i, "R"]
  outbreak_length_list[[i]]$k <- scenarios[i, "k"]
  outbreak_length_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_length_list[[i]]$statistic <- scenarios[i, "statistic"]
}
outbreak_length <- do.call(rbind, outbreak_length_list)
head(outbreak_length)

nbinom_length <- ggplot2::ggplot(data = outbreak_size) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(
    name = "Reproduction number (R)",
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak length\n(secondary generations)", 
    palette = "Spectral"
  ) + 
  ggplot2::facet_wrap(
    facets = c("k"), 
    labeller = ggplot2::label_both
  ) +
  ggplot2::labs(title = "Negative Binomial distribution") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  )

length <- cowplot::plot_grid(
  pois_length + ggplot2::theme(legend.position = "none"), 
  nbinom_length, 
  nrow = 1,
  rel_widths = c(0.5, 1),
  labels = c("A", "B")
)

ggplot2::ggsave(
  filename = file.path("plots", "outbreak_length.png"), 
  plot = length,
  device = "png", 
  width = 250, 
  height = 150,
  units = "mm",
  dpi = 300
)

