#' Infection trajectory visualisation
#'
#' @param infection_traj simulated trajectories
#' @param n_shown how many to visualise
#' @param right_truncate how many days to truncate at the end
#'
#' @returns plot
#'
#' @importFrom dplyr bind_rows mutate filter
#' @importFrom tidyr pivot_longer any_of starts_with
#' @importFrom ggplot2 ggplot facet_wrap geom_line scale_x_date theme
#' @importFrom cowplot theme_cowplot panel_border
#'
#' @export
plot_infection_traj <- function (infection_traj,
                                 n_shown = 10,
                                 right_truncate = 14) {

  # combine into one data frame
  infection_traj_df <- dplyr::bind_rows(infection_traj) %>%
    # pivot for plotting
    tidyr::pivot_longer(cols = tidyr::any_of(tidyr::starts_with("draw")),
                        names_to = "draw",
                        values_to = "value") %>%
    # turn draw label back to numeric
    dplyr::mutate(draw = gsub("draw", "", draw)) %>%
    dplyr::mutate(draw = as.numeric(draw),
           date = as.Date(date))

  # sample a few trajectory, if wanted
  if (!is.null(n_shown)) {
    sampled_traj <- sample(max(infection_traj_df$draw),
                           n_shown)
    infection_traj_df <- infection_traj_df %>%
      dplyr::filter(draw %in% sampled_traj)
  }

  # remove most recent days, if wanted
  if (!is.null(right_truncate)) {
    date_cutoff <- (max(infection_traj_df$date) - right_truncate)
    infection_traj_df <- infection_traj_df %>%
      filter(date <= date_cutoff)
  }

  # plot
  infection_traj_df %>%
    dplyr::mutate(draw = as.factor(draw)) %>%
    ggplot2::ggplot(ggplot2::aes(x = date, y = value, colour = draw)) +
    ggplot2::facet_wrap(~jurisdiction) +
    ggplot2::geom_line(alpha = 0.5) +
    cowplot::theme_cowplot() +
    cowplot::panel_border(remove = TRUE) +
    ggplot2::scale_x_date(date_breaks = "1 month",
                          date_labels = "%b") +
    ggplot2::theme(legend.position = "none")

}
