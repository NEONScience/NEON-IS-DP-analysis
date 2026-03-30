# install.packages(c("ggplot2", "dplyr", "forcats", "scales"))  # if needed
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)

wind_rose <- function(
  data,
  dir_col = "dir",          # name of direction column (degrees 0..360)
  spd_col = "spd",          # name of speed column
  dir_bin = 22.5,           # sector width (e.g., 22.5 -> 16 sectors)
  spd_breaks = c(0, 1, 3, 5, 8, 11, 15, Inf), # Beaufort-ish default (customize freely)
  spd_labels = NULL,        # optional custom labels (length = length(spd_breaks) - 1)
  calm_threshold = 0.2,     # <= this is treated as "Calm"
  missing_dir_policy = c("drop", "impute"), # typically drop NA directions
  palette = scales::viridis_pal(option = "C")(length(spd_breaks) - 1),
  title = "Wind Rose",
  subtitle = NULL
) {
  missing_dir_policy <- match.arg(missing_dir_policy)

  stopifnot(dir_col %in% names(data), spd_col %in% names(data))

  df <- data %>%
    select(dir = all_of(dir_col), spd = all_of(spd_col)) %>%
    # sanitize input
    mutate(
      dir = as.numeric(dir),
      spd = as.numeric(spd)
    )

  # Handle NA directions
  if (missing_dir_policy == "drop") {
    df <- df %>% filter(!is.na(dir))
  } else {
    # If impute: clamp to [0,360) and replace missing with 0 (rarely desired)
    df <- df %>%
      mutate(dir = ifelse(is.na(dir), 0, dir))
  }

  # Normalize directions to [0, 360)
  df <- df %>%
    mutate(dir = (dir %% 360))

  # Mark calms
  df <- df %>%
    mutate(is_calm = !is.na(spd) & spd <= calm_threshold)

  # Define speed bins (excluding calms)
  if (is.null(spd_labels)) {
    # Create nice labels like "0–1", "1–3", ..., "15+"
    spd_labels <- paste0(
      head(spd_breaks, -1),
      "–",
      replace(tail(spd_breaks, -1), is.infinite(tail(spd_breaks, -1)), "+")
    )
  }

  df <- df %>%
    mutate(
      spd_bin = case_when(
        is_calm ~ "Calm",
        TRUE    ~ as.character(cut(spd, breaks = spd_breaks, labels = spd_labels, right = FALSE))
      ),
      # Build direction sectors centered on standard bearings.
      # Example with 22.5° sectors: -11.25..11.25 => N, 11.25..33.75 => NNE, etc.
      # We'll compute sector index by floor((dir + dir_bin/2)/dir_bin).
      sector_idx = floor((dir + dir_bin / 2) / dir_bin),
      sector_deg = (sector_idx %% round(360 / dir_bin)) * dir_bin
    )

  # Create human-readable sector labels (N, NNE, NE, ...)
  # For standard sector names up to 16, we can map by index.
  n_sectors <- round(360 / dir_bin)
  stopifnot(360 %% dir_bin == 0, n_sectors >= 4)
  # Standard 16-wind names; if you use other dir_bin values, we'll fallback to degrees.
  standard_16 <- c("N","NNE","NE","ENE","E","ESE","SE","SSE",
                   "S","SSW","SW","WSW","W","WNW","NW","NNW")
  sector_names <- if (n_sectors == 16) {
    standard_16
  } else if (n_sectors == 8) {
    c("N","NE","E","SE","S","SW","W","NW")
  } else if (n_sectors == 4) {
    c("N","E","S","W")
  } else {
    # Fallback: label by degrees
    paste0((0:(n_sectors - 1)) * dir_bin, "°")
  }

  # Map sector_deg to factor labels
  df <- df %>%
    mutate(
      sector_order = (sector_deg / dir_bin) + 1L,
      sector_lab = factor(
        sector_names[pmax(1, pmin(n_sectors, sector_order))],
        levels = sector_names
      ),
      spd_bin = factor(spd_bin, levels = c("Calm", spd_labels))
    )

  # Aggregate to frequencies
  counts <- df %>%
    count(sector_lab, spd_bin, name = "n")

  total_n <- sum(counts$n, na.rm = TRUE)
  if (total_n == 0) {
    stop("No observations after preprocessing. Check your inputs/filters.")
  }

  freq <- counts %>%
    group_by(sector_lab) %>%
    mutate(
      sector_total = sum(n, na.rm = TRUE),
      pct_of_total = (n / total_n) * 100
    ) %>%
    ungroup()

  # Keep "Calm" separate (as a center annotation) if desired.
  calm_n <- freq %>% filter(spd_bin == "Calm") %>% summarise(n = sum(n)) %>% pull(n)
  calm_pct <- if (!is.null(calm_n) && length(calm_n) == 1 && total_n > 0) 100 * calm_n / total_n else 0

  # Plot stacked bars by speed bin in polar coords
  p <- ggplot(
    data = freq %>% filter(spd_bin != "Calm"),
    aes(x = sector_lab, y = pct_of_total, fill = spd_bin)
  ) +
    geom_col(width = 0.95, color = "grey20", size = 0.2) +
    coord_polar(start = -dir_bin/2*pi/180, direction = 1) +  # start at North, clockwise
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("Calm" = "grey80", palette)) +
    guides(fill = guide_legend(title = "Speed")) +
    labs(
      title = title,
      subtitle = if (is.null(subtitle)) sprintf("Calm: %s%%", round(calm_pct, 1)) else subtitle,
      x = NULL,
      y = "Percent of all observations"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_text(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  return(p)
}
