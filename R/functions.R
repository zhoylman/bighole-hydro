######################################
#      DEFINE HELPER FUNCTIONS       #
######################################

get_gridmet_daily_tmmx_tmmn_sf = function(sites_sf,
                                          out_units = c("C","K","F"),
                                          keep_cols = c("network","site_no","station_nm")) {
  # ------------------------------------------------------------------------------
  # Title:    Extract Daily gridMET Max/Min Temperature for Multiple Sites
  # Author:   Zachary H. Hoylman
  # Date:     8-17-2025
  #
  # Description:
  #   This function extracts daily maximum (tmmx) and minimum (tmmn) temperature
  #   values from the gridMET dataset using the THREDDS OPeNDAP service.
  #   It accepts an sf object containing POINT geometries and associated site
  #   metadata, identifies the nearest gridMET cell for each point, and returns
  #   a tidy tibble of daily temperatures in user-specified units (C, K, or F).
  #
  # Input:
  #   sites_sf   : An sf object with POINT geometries (WGS84 CRS).
  #                Must contain at least one unique identifier column if you
  #                want metadata preserved. Common fields: "network", "site_no",
  #                "station_nm".
  #   out_units  : Units for output temperature. One of:
  #                  - "C" = degrees Celsius
  #                  - "K" = Kelvin
  #                  - "F" = degrees Fahrenheit
  #   keep_cols  : Character vector of column names in sites_sf to preserve
  #                in the output (e.g., c("network","site_no","station_nm")).
  #
  # Output:
  #   A tibble with columns:
  #     - id columns (from keep_cols, if present)
  #     - lon, lat   (decimal degrees, WGS84)
  #     - date       (daily Date stamp)
  #     - tmmx       (daily maximum temperature, in out_units)
  #     - tmmn       (daily minimum temperature, in out_units)
  #
  # Notes:
  #   * gridMET tmmx and tmmn values are in Kelvin; these are converted on the fly.
  #   * Daily slices correspond approximately to calendar days ending at midnight MST.
  #   * The function opens the netCDF files once and reuses the connection for
  #     efficiency. Each site is matched to its nearest gridMET cell.
  #   * If duplicate coordinates exist in the sf object, they are handled gracefully
  #     and replicated in the output.
  #   * Requires packages: sf, ncdf4, dplyr, purrr, tidyr, lubridate, tibble
  #
  # Example:
  #   daily_temps = get_gridmet_daily_tmmx_tmmn_sf(all_site_meta, out_units = "C")
  #   dplyr::glimpse(daily_temps)
  #
  # ------------------------------------------------------------------------------
  out_units = match.arg(out_units)
  
  # ---- sanity & coords ----------------------------------------------------
  stopifnot(inherits(sites_sf, "sf"))
  # ensure WGS84 lon/lat
  if (sf::st_is_longlat(sites_sf) == FALSE ||
      is.na(sf::st_crs(sites_sf))) {
    sites_sf = sf::st_set_crs(sites_sf, 4326)
  } else if (sf::st_crs(sites_sf)$epsg != 4326) {
    sites_sf = sf::st_transform(sites_sf, 4326)
  }
  
  coords = sf::st_coordinates(sites_sf)  # columns: X (lon), Y (lat)
  sites_tbl =
    sites_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(lon = coords[, "X"], lat = coords[, "Y"]) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(keep_cols), as.character, .names = "{col}")) |>
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # choose id columns to carry through (only those present)
  id_cols = intersect(keep_cols, names(sites_tbl))
  if (length(id_cols) == 0L) id_cols = ".row_id"
  
  # de-duplicate exact same lon/lat pairs to avoid redundant reads
  points_unique =
    sites_tbl |>
    dplyr::distinct(lon, lat, .keep_all = TRUE) |>
    dplyr::select(dplyr::all_of(id_cols), lon, lat)
  
  # ---- open gridMET once --------------------------------------------------
  urls = list(
    tmmx = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc",
    tmmn = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmn_1979_CurrentYear_CONUS.nc"
  )
  nc_x = ncdf4::nc_open(urls$tmmx); on.exit(try(ncdf4::nc_close(nc_x), silent = TRUE), add = TRUE)
  nc_n = ncdf4::nc_open(urls$tmmn); on.exit(try(ncdf4::nc_close(nc_n), silent = TRUE), add = TRUE)
  
  lon_vals = nc_x$dim$lon$vals
  lat_vals = nc_x$dim$lat$vals
  time_len = nc_x$dim$day$len
  time_num = ncdf4::ncvar_get(nc_x, "day", start = 1, count = time_len)
  dates = as.Date(time_num, origin = "1900-01-01")
  
  # pick the 3D data variables (robust to dim order/names)
  .pick_var = function(nc) {
    vname =
      names(nc$var)[sapply(nc$var, function(v) {
        length(v$dim) == 3 &&
          all(c("lon","lat","day") %in% vapply(v$dim, function(d) d$name, character(1)))
      })][1]
    if (is.na(vname) || !nzchar(vname)) stop("No 3D data variable found.")
    vname
  }
  v_x = .pick_var(nc_x)
  v_n = .pick_var(nc_n)
  
  # read one point from one nc file
  .read_point = function(nc, vname, lat_in, lon_in) {
    ix_lon = which.min(abs(lon_vals - lon_in))
    ix_lat = which.min(abs(lat_vals - lat_in))
    
    dims = vapply(nc$var[[vname]]$dim, function(d) d$name, character(1))
    start = count = rep(1L, 3)
    start[match("lon", dims)] = ix_lon; count[match("lon", dims)] = 1L
    start[match("lat", dims)] = ix_lat; count[match("lat", dims)] = 1L
    start[match("day", dims)] = 1L;     count[match("day", dims)] = time_len
    
    vals = ncdf4::ncvar_get(nc, vname, start = start, count = count)
    miss = ncdf4::ncatt_get(nc, vname, "_FillValue")$value
    if (is.null(miss)) miss = ncdf4::ncatt_get(nc, vname, "missing_value")$value
    if (!is.null(miss) && is.finite(miss)) vals[vals == miss] = NA_real_
    
    as.numeric(vals)
  }
  
  # unit conversion from Kelvin
  .K_to = function(x, units = out_units) {
    dplyr::case_when(
      units == "K" ~ x,
      units == "C" ~ x - 273.15,
      units == "F" ~ (x - 273.15) * 9/5 + 32
    )
  }
  
  # ---- pull both series for each unique point ----------------------------
  read_one = function(lon, lat) {
    tmmx_K = .read_point(nc_x, v_x, lat, lon)
    tmmn_K = .read_point(nc_n, v_n, lat, lon)
    tibble::tibble(
      date = dates,
      tmmx = .K_to(tmmx_K),
      tmmn = .K_to(tmmn_K)
    )
  }
  
  unique_results =
    points_unique |>
    dplyr::mutate(data = purrr::pmap(list(lon = lon, lat = lat), read_one)) |>
    tidyr::unnest(data)
  
  # join back to all rows (so duplicate coordinates replicate correctly)
  out =
    sites_tbl |>
    dplyr::select(dplyr::all_of(id_cols), lon, lat) |>
    dplyr::left_join(unique_results, by = c(id_cols, "lon","lat")) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(id_cols)), date)
  
  out
}

fit_per_station_gam = function(df, by = 0.1) {
  # ------------------------------------------------------------------------------
  # Title:    Per-Station Temporal Smoothers for Thermal Exceedance Hours
  # Author:   Zachary H. Hoylman
  # Date:     8-19-2025
  #
  # Description:
  #   This function fits site-level generalized additive models (GAMs) to estimate
  #   smooth nonlinear trends in the cumulative annual hours exceeding thermal
  #   thresholds (e.g., 21 °C, 25 °C). The response is modeled with a
  #   quasipoisson(log) family, ensuring non-negative predictions and robustness
  #   to overdispersion. Basis dimension (k) is kept modest to avoid overfitting
  #   when few years of data are available. For sites with all-zero observations,
  #   the function returns flat zero predictions. Predictions are returned on a
  #   fine year grid with fitted values and 95% confidence intervals.
  #
  # Justification:
  #   A GAM framework is used because exceedance hours often vary nonlinearly
  #   through time—declining in some periods, rising in others—making linear GLMs
  #   inadequate. The smooth spline penalty balances flexibility with parsimony,
  #   allowing site-specific curves to adapt to available data while avoiding
  #   spurious wiggles. Using a log link maintains interpretability on the odds/
  #   rate scale and ensures predictions remain ≥ 0. Compared to GLMs, GAMs better
  #   capture nonlinear trends; compared to more complex GLMMs/GAMMs, this model
  #   is appropriate for annual aggregates where serial autocorrelation is minimal.
  #   Importantly, this is used simply as a visual aid. 
  #
  # Input:
  #   dat_i : A tibble for a single site with columns:
  #           - year (numeric)
  #           - cumulative_hours (non-negative integer)
  #
  # Output:
  #   Tibble with station name, prediction year grid, fitted mean,
  #   and lower/upper 95% CI bounds (on the response scale).
  #
  # ------------------------------------------------------------------------------
  df %>%
    dplyr::group_by(station_nm) %>%
    dplyr::group_modify(~{
      dat_i = .x %>% dplyr::arrange(year)
      ny    = dplyr::n_distinct(dat_i$year)
      
      # prediction grid across observed span
      new_year = tibble::tibble(year = seq(min(dat_i$year), max(dat_i$year), by = by))
      
      # ---- all zeros → return flat zero ----
      if (all(dat_i$cumulative_hours == 0, na.rm = TRUE)) {
        return(tibble::tibble(
          station_nm = dat_i$station_nm[1],
          year = new_year$year,
          fit  = 0,
          lwr  = 0,
          upr  = 0
        ))
      }
      
      # modest k to avoid wiggles
      k_use = max(3, min(5, floor(ny/3) + 1))
      
      mod = mgcv::gam(
        cumulative_hours ~ s(year, k = k_use),
        data   = dat_i,
        family = quasipoisson(link = "log"),
        method = "REML"
      )
      
      pr = predict(mod, newdata = new_year, se.fit = TRUE, type = "link")
      
      tibble::tibble(
        station_nm = dat_i$station_nm[1],
        year = new_year$year,
        fit  = exp(pr$fit),
        lwr  = exp(pr$fit - 1.96 * pr$se.fit),
        upr  = exp(pr$fit + 1.96 * pr$se.fit)
      ) %>%
        dplyr::mutate(dplyr::across(c(fit, lwr, upr), ~ pmax(., 0)))
    }) %>%
    dplyr::ungroup()
}

plot_exceedance = function(df,
                           y_col = c("cumulative_hours_21","cumulative_hours_25"),
                           smooth_df = NULL,
                           response_col = NULL,   # if smooth_df is NULL, compute via fit_per_station_gam()
                           facet_var = "station_nm",
                           title = NULL,
                           y_lab = NULL,
                           caption = NULL,
                           nrow_facets = 3,
                           width = 10, height = 8,
                           save_path = NULL) {
  # ------------------------------------------------------------------------------
  # Title:    Plot Thermal Exceedance Hours with Site-Level GAM Smooths
  # Author:   Zachary H. Hoylman
  # Date:     8-19-2025
  #
  # Description:
  #   This function generates site-level plots of cumulative thermal exceedance
  #   hours (e.g., >21 °C or >25 °C) through time. Raw annual exceedance counts
  #   are shown as points, with smoothed trajectories estimated using per-station
  #   GAM fits (quasipoisson log-link, ensuring non-negative predictions).
  #   95% confidence ribbons are drawn from model-based uncertainty, and each
  #   site is displayed in a separate facet.
  #
  # Justification:
  #   GAM smoothing provides a flexible yet controlled description of nonlinear
  #   temporal patterns in exceedance hours, capturing declines and increases
  #   without imposing linearity. The quasipoisson family accounts for overdispersed
  #   counts and ensures fitted values remain ≥ 0. Faceting by station highlights
  #   site-level differences while maintaining common structure. This visualization
  #   is intended as a descriptive tool rather than an inferential test.
  #
  # Input:
  #   df           : Data frame with at least `year`, the response column
  #                  (e.g., cumulative_hours_21), and a station identifier.
  #   y_col        : Name of the exceedance column to plot ("cumulative_hours_21"
  #                  or "cumulative_hours_25").
  #   smooth_df    : Optional pre-computed smoother tibble from fit_per_station_gam().
  #   response_col : If smoother not provided, this specifies which response column
  #                  to pass to fit_per_station_gam().
  #   facet_var    : Station ID variable for faceting (default = "station_nm").
  #   title, y_lab, caption : Optional plot labels (auto-filled if NULL).
  #   nrow_facets  : Number of rows in facet layout.
  #   width,height : Plot dimensions for saved output.
  #   save_path    : Optional file path to save figure (PNG, PDF, etc).
  #
  # Output:
  #   A ggplot object with raw exceedance points, fitted GAM smooth, and 95% CI band.
  #   Optionally, the plot is saved to disk.
  # ------------------------------------------------------------------------------
  
  y_col = rlang::as_name(rlang::ensym(y_col))
  
  # sanity
  if (!("year" %in% names(df))) stop("`df` must have a `year` column.")
  if (!(y_col %in% names(df))) stop("`df` must contain the column '", y_col, "'.")
  if (!(facet_var %in% names(df))) stop("`facet_var` '", facet_var, "' not found in `df`.")
  
  # build smoother if not supplied
  if (is.null(smooth_df)) {
    if (is.null(response_col)) response_col = y_col
    smooth_df = fit_per_station_gam(df, response_col = response_col, group_var = facet_var)
  }
  
  # choose default labels if not provided
  if (is.null(title)) {
    thr = if (grepl("25", y_col)) "25" else "21"
    title = paste0("Cumulative Hours > ", thr, "°C")
  }
  if (is.null(y_lab)) {
    thr = if (grepl("25", y_col)) "25" else "21"
    y_lab = paste0("Cumulative Hours Above ", thr, "°C")
  }
  if (is.null(caption)) {
    caption = stringr::str_wrap("Data: USGS & MT DNRC", width = 140)
  }
  
  # x breaks (one per observed year)
  x_breaks = sort(unique(df$year))
  
  # tidy-eval for aesthetics
  facet_sym = rlang::sym(facet_var)
  y_sym     = rlang::sym(y_col)
  
  p = ggplot() +
    geom_point(
      data = df,
      aes(x = year, y = !!y_sym),
      color = "steelblue", size = 2, alpha = 0.7
    ) +
    geom_ribbon(
      data = smooth_df,
      aes(x = year, ymin = pmax(lwr, 0), ymax = pmax(upr, 0)),
      alpha = 0.2
    ) +
    geom_line(
      data = smooth_df,
      aes(x = year, y = pmax(fit, 0)),
      color = "black", linewidth = 1
    ) +
    facet_wrap(vars(!!facet_sym), scales = "free_y", nrow = nrow_facets) +
    labs(
      title = title,
      x = "Year",
      y = y_lab,
      caption = caption
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_continuous(breaks = x_breaks)
  
  if (!is.null(save_path)) {
    ggsave(p, file = save_path, width = width, height = height)
  }
  
  p
}


make_exceedance_fixed_effects_table = function(
    fit,
    threshold   = 21,
    labels      = list(
      intercept = "Intercept",
      flow      = "Streamflow (per +1 CFS)",
      tmmx      = "Max Air Temperature (per +1 °C)",
      tmmn      = "Min Air Temperature (per +1 °C)",
      period    = "Period: pre-2017 vs post-2017"
    ),
    period_term = "periodpre2017",
    interp_text = NULL,
    html_path   = NULL,
    png_path    = NULL
) {
  # ------------------------------------------------------------------------------
  # Title:    Fixed-Effects Table (AR(1) Binomial GAMM) with ORs & 95% CI
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-20
  #
  # Description:
  #   Given a fitted mgcv::bam object from your exceedance model (with AR(1),
  #   site random intercepts, and site-specific tmmx slopes), this function:
  #     * extracts parametric (fixed) effects,
  #     * computes odds ratios and 95% CI,
  #     * builds a gt table with clean labels and a short interpretation note,
  #     * optionally saves HTML/PNG outputs.
  #
  # Args:
  #   fit         : mgcv::bam model object (binomial link=logit)
  #   threshold   : numeric; exceedance threshold for subtitle (default 21)
  #   labels      : list of pretty labels used in the table & interpretation
  #                 (see defaults below)
  #   period_term : character; the coefficient name for the pre/post effect in `summary(fit)$p.table`
  #                 (default "periodpre2017")
  #   html_path   : optional file path to save the gt table as HTML
  #   png_path    : optional file path to save the gt table as PNG (needs webshot2)
  #
  # Returns:
  #   list(
  #     param_tab = <tidy fixed-effects data.frame with ORs>,
  #     gt_table  = <gt object>,
  #     interp    = <character interpretation string>
  #   )
  # ------------------------------------------------------------------------------
  
  # ---- extract parametric table safely ----
  pt = as.data.frame(summary(fit)$p.table)
  pt$term = rownames(pt)
  
  # normalize column names (some locales print "Std. Error" with a space)
  names(pt) = gsub("\\s+", "_", names(pt))
  stopifnot(all(c("Estimate","Std._Error","z_value","Pr(>|z|)","term") %in% names(pt)))
  
  param_tab = pt |>
    tibble::as_tibble() |>
    dplyr::rename(
      estimate = Estimate,
      se       = Std._Error,
      z        = z_value,
      p        = `Pr(>|z|)`
    ) |>
    dplyr::mutate(
      OR    = exp(estimate),
      OR_lo = exp(estimate - 1.96 * se),
      OR_hi = exp(estimate + 1.96 * se),
      Term  = dplyr::case_when(
        term == "(Intercept)"             ~ labels$intercept,
        term == "daily_median_discharge"  ~ labels$flow,
        term == "tmmx"                    ~ labels$tmmx,
        term == "tmmn"                    ~ labels$tmmn,
        term == period_term               ~ labels$period,
        TRUE                              ~ term
      )
    ) |>
    dplyr::select(Term, estimate, se, z, p, OR, OR_lo, OR_hi)
  
  # ---- pull rows for interpretation (if present) ----
  get_row = function(lbl) dplyr::filter(param_tab, Term == lbl)
  r_flow  = get_row(labels$flow)
  r_tmx   = get_row(labels$tmmx)
  r_tmn   = get_row(labels$tmmn)
  r_per   = get_row(labels$period)
  
  fmt_p = function(x) ifelse(is.na(x), NA_character_,
                             ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x)))
  
  # You can expand this to a full narrative if desired; leaving placeholder per your draft
  interp = paste0("**Interpretation.** ", interp_text)
  
  # ---- build gt table ----
  fixed_gt = param_tab |>
    dplyr::mutate(
      dplyr::across(c(estimate, se, z), ~ round(., 3)),
      p        = format.pval(p, digits = 3, eps = .001),
      OR       = round(OR, 2),
      `95% CI (OR)` = sprintf("[%0.2f, %0.2f]", round(OR_lo, 2), round(OR_hi, 2))
    ) |>
    dplyr::select(Term, estimate, se, z, p, OR, `95% CI (OR)`) |>
    gt::gt() |>
    gt::tab_header(
      title    = gt::md("**Exceedance Model (AR(1) GAMM): Fixed Effects**"),
      subtitle = paste0("Outcome: Pr(stream daily max > ", threshold, " °C)")
    ) |>
    gt::cols_label(
      Term     = "Predictor",
      estimate = "Log-odds",
      se       = "SE",
      z        = "z",
      p        = "p",
      OR       = "Odds ratio",
      `95% CI (OR)` = "95% CI"
    ) |>
    gt::tab_source_note(gt::md(interp))
  
  # ---- optional saves ----
  if (!is.null(html_path)) {
    gt::gtsave(fixed_gt, filename = html_path)
  }
  if (!is.null(png_path)) {
    # requires webshot2; catch and message if missing
    try({
      gt::gtsave(fixed_gt, filename = png_path, expand = 5, zoom = 1.25)
    }, silent = TRUE)
  }
  
  list(
    param_tab = param_tab,
    gt_table  = fixed_gt,
    interp    = interp
  )
}

# ------------------------------------------------------------------------------
# Simple HUC8 Sites Map (NHD remote fetch with robust TNM fallback, no flowlines)
# ------------------------------------------------------------------------------

plot_big_hole_sites_nhd = function(
    sites,                         # sf POINTS with cols: network, site_no, station_nm
    huc8      = "10020004",        # Big Hole River HUC8
    out_crs   = 5070,              # NAD83 / Conus Albers (planar)
    label     = TRUE,
    label_col = "site_no",
    palette   = c(USGS = "#1f77b4", StAGE = "#2ca02c"),
    title     = "Monitoring Sites in the Big Hole Watershed",
    caption   = "Data: USGS & MT DNRC StAGE",
    save_path = NULL,
    width = 9, height = 7, dpi = 300,
    # NEW: optional ArcGIS REST layer for NHD flowlines
    flowlines_url = NULL           # e.g. "https://<...>/FeatureServer/3" or ".../MapServer/3"
) {
  # ------------------------------------------------------------------------------
  # Title:    Big Hole Watershed Sites Map (HUC8 from NHD/TNM)
  # Author:   Zachary H. Hoylman
  # Date:     8-20-2025
  #
  # Description:
  #   Produces a publication-ready map of monitoring sites within the Big Hole 
  #   River watershed. The HUC8 boundary is fetched remotely from the NHDPlus/TNM
  #   services (with multiple fallbacks to ensure robustness). Sites are colored 
  #   by network (e.g., USGS, StAGE) and labeled if requested. 
  #
  # Justification:
  #   Visualizing site locations relative to the watershed provides important 
  #   spatial context for interpreting exceedance analyses and temperature 
  #   modeling. Using an authoritative NHD/TNM boundary ensures consistency with 
  #   hydrologic datasets. Stream flowlines are intentionally omitted here for 
  #   simplicity and to avoid brittle remote calls.
  #
  # Inputs:
  #   sites   : sf POINTS with required columns:
  #               - network   (factor or character; e.g., USGS, StAGE)
  #               - site_no   (unique site identifier)
  #               - station_nm (human-readable site name)
  #   huc8    : character; HUC8 code to fetch (default = "10020004").
  #   out_crs : integer EPSG; projection for output (default = 5070).
  #   label   : logical; whether to annotate sites.
  #   label_col : character; column in `sites` used for labels.
  #   palette : named vector of colors keyed by network values.
  #   title, caption : plot text.
  #   save_path : optional file path to save the figure.
  #   width,height,dpi : graphics device settings if saving.
  #
  # Output:
  #   A ggplot2 object showing watershed boundary and monitoring sites,
  #   optionally saved as a static figure (PNG/PDF/etc).
  # ------------------------------------------------------------------------------
  stopifnot(inherits(sites, "sf"))
  need = c("network","site_no","station_nm")
  miss = setdiff(need, names(sites))
  if (length(miss)) stop("`sites` is missing: ", paste(miss, collapse = ", "))
  if (!label_col %in% names(sites)) stop("`label_col` not found in `sites`.")
  
  # -- helper: robust HUC8 fetch via nhdplusTools then TNM GeoJSON fallback
  fetch_huc8 = function(huc8, crs_out) {
    bh = NULL
    if (requireNamespace("nhdplusTools", quietly = TRUE)) {
      bh = try(nhdplusTools::get_huc(id = huc8, type = "huc08"), silent = TRUE)
      if (inherits(bh, "try-error") || !inherits(bh, "sf") || nrow(bh) == 0) {
        bh = try(nhdplusTools::subset_wbd(huc = huc8), silent = TRUE)
        if (is.list(bh) && !inherits(bh, "sf")) {
          idx = which(vapply(bh, inherits, logical(1), what = "sf"))
          if (length(idx)) bh = bh[[idx[1]]]
        }
      }
    }
    if (!inherits(bh, "sf") || nrow(bh) == 0) {
      base  = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"
      layers = c(8, 6, 4) # try a few common HU8 layers
      for (lyr in layers) {
        url = sprintf(
          "%s/%d/query?where=HUC8%%3D%%27%s%%27&outFields=*&returnGeometry=true&f=geojson",
          base, lyr, utils::URLencode(huc8, reserved = TRUE)
        )
        bh_try = try(sf::read_sf(url), silent = TRUE)
        if (!inherits(bh_try, "try-error") && inherits(bh_try, "sf") && nrow(bh_try) > 0) {
          bh = bh_try
          break
        }
      }
    }
    if (!inherits(bh, "sf") || nrow(bh) == 0)
      stop("Failed to fetch HUC8 boundary for '", huc8, "'.")
    if (is.na(sf::st_crs(bh))) sf::st_crs(bh) = sf::st_crs(4326)
    sf::st_transform(bh, crs_out)
  }
  
  # -- NEW helper: fetch flowlines from ArcGIS REST for the HUC8 bbox, then clip
  # Fetch lighter-weight flowlines from ArcGIS REST, then clip to a polygon
  fetch_flowlines_arcgis = function(
    flowlines_layer_url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/3",
    aoi_sf,
    out_crs = 5070,            # project to meters so maxAllowableOffset is in meters
    min_order = 3,             # keep only StreamOrde >= 3 (tweak as needed)
    exclude_artificial = TRUE, # drop artificial paths (FCode 55800)
    max_offset_m = 100,        # server-side geometry generalization (meters)
    fields = c("GNIS_ID","GNIS_NAME","StreamOrde","FCode") # keep it lean
  ) {
    stopifnot(inherits(aoi_sf, "sf"))
    aoi_wgs84 = sf::st_transform(aoi_sf, 4326)
    bb = sf::st_bbox(aoi_wgs84)
    env = paste(bb$xmin, bb$ymin, bb$xmax, bb$ymax, sep = ",")
    
    # Attribute filter
    where_parts = c(sprintf("StreamOrde >= %d", min_order))
    if (exclude_artificial) where_parts = c(where_parts, "FCode <> 55800")
    where = paste(where_parts, collapse = " AND ")
    
    # Build query (ask for generalization + projected outSR for meter units)
    base = sub("\\/+$", "", flowlines_layer_url)
    q = sprintf(
      "%s/query?where=%s&geometry=%s&geometryType=esriGeometryEnvelope&inSR=4326&spatialRel=esriSpatialRelIntersects&outFields=%s&outSR=%s&returnGeometry=true&f=geojson&maxAllowableOffset=%s",
      base,
      utils::URLencode(where, reserved = TRUE),
      utils::URLencode(env,   reserved = TRUE),
      utils::URLencode(paste(fields, collapse=","), reserved = TRUE),
      as.integer(out_crs),
      as.numeric(max_offset_m)
    )
    
    fl = try(sf::read_sf(q), silent = TRUE)
    if (inherits(fl, "try-error") || !inherits(fl, "sf") || nrow(fl) == 0) return(NULL)
    
    # Transform & hard-clip to your watershed polygon
    fl = sf::st_transform(fl, sf::st_crs(aoi_sf))
    bh = aoi_sf
    bh_valid = try(sf::st_make_valid(bh), silent = TRUE); if (inherits(bh_valid, "try-error")) bh_valid = bh
    fl_crop = suppressWarnings(sf::st_crop(fl, sf::st_bbox(bh_valid)))
    fl_clip = try(suppressWarnings(sf::st_intersection(fl_crop, bh_valid)), silent = TRUE)
    if (inherits(fl_clip, "try-error")) fl_clip = fl_crop
    
    # Optional extra simplification (client-side) if still too detailed
    # fl_clip = sf::st_simplify(fl_clip, dTolerance = 50)  # ~50 m
    
    fl_clip
  }
  
  # HUC8 boundary & sites to map CRS
  big_hole = fetch_huc8(huc8, out_crs)
  sites    = sf::st_transform(sites, sf::st_crs(big_hole))
  
  # After you’ve fetched big_hole and transformed sites:
  flowlines = fetch_flowlines_arcgis(
    aoi_sf = big_hole,
    out_crs = 5070,        # meters
    min_order = 1,         # try 3–4 for cleaner network
    exclude_artificial = TRUE,
    max_offset_m = 150     # bigger number => simpler lines
  )
  
  # Optional: flowlines
  flows = NULL
  if (!is.null(flowlines_url)) {
    flows = fetch_flowlines_arcgis(flowlines_url, big_hole, sf::st_crs(big_hole))
  }
  
  # Colors (add any unexpected networks)
  nets  = sort(unique(as.character(sites$network)))
  extra = setdiff(nets, names(palette))
  if (length(extra)) {
    auto_cols = grDevices::hcl.colors(length(extra), palette = "Dark3")
    names(auto_cols) = extra
    palette = c(palette, auto_cols)
  }
  
  # Label coords (planar)
  label_df = NULL
  if (label) {
    coords = sf::st_coordinates(sites)
    label_df = cbind(sf::st_drop_geometry(sites), as.data.frame(coords))
  }
  
  # bbox padding
  bb  = sf::st_bbox(big_hole)
  pad = 0.01 * max(bb$xmax - bb$xmin, bb$ymax - bb$ymin)
  xlim = c(bb$xmin - pad, bb$xmax + pad)
  ylim = c(bb$ymin - pad, bb$ymax + pad)
  
  # Plot
  p = ggplot2::ggplot() +
    ggplot2::geom_sf(data = big_hole, fill = scales::alpha("grey80", 0.25),
                     color = "grey30", linewidth = 0.4) +
    { if (!is.null(flows))
      ggplot2::geom_sf(data = flows, color = "steelblue", linewidth = 0.3, alpha = 0.75)
      else NULL } +
    ggplot2::geom_sf(data = sites, ggplot2::aes(color = network), size = 2.8, alpha = 0.95) +
    {
      if (label && !is.null(label_df)) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(x = X, y = Y, label = .data[[label_col]]),
            size = 3, seed = 42, box.padding = 0.4,
            segment.color = "grey40", max.overlaps = Inf
          )
        } else {
          ggplot2::geom_text(
            data = label_df,
            ggplot2::aes(x = X, y = Y, label = .data[[label_col]]),
            nudge_y = 2000, size = 3, color = "grey20"
          )
        }
      } else NULL
    } +
    ggplot2::scale_color_manual(values = palette, name = "Network") +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(title = title, caption = caption, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption     = ggplot2::element_text(hjust = 0) # left-justify caption
    )
  
  if (!is.null(flowlines)) {
    p = p + geom_sf(data = flowlines, color = "steelblue", linewidth = 0.3, alpha = 0.4)
  }
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  p
}

fit_exceedance_bam = function(data,
                              threshold  = 21,
                              break_date = as.Date("2017-01-01")) {
  # ------------------------------------------------------------------------------
  # Title:    Fit AR(1) Binomial GAMM for >21°C Exceedance (site RE + tmmx slopes)
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-20
  #
  # Description:
  #   Given a daily data set with stream temperature, air temperatures (tmmx, tmmn),
  #   discharge, site, and date, this function:
  #     (1) constructs a binary exceedance indicator (> threshold °C),
  #     (2) builds a pre/post period factor (default breakpoint = 2017-01-01),
  #     (3) de-duplicates (site_no, date), orders within site, and flags AR(1) starts,
  #     (4) estimates an AR(1) rho from a baseline GLM’s Pearson residual ACF,
  #     (5) fits mgcv::bam with binomial(logit), random site intercepts, and
  #         site-specific slopes for tmmx:  over21 ~ flow + tmmx + tmmn + period
  #                                         + s(site_no, bs="re")
  #                                         + s(site_no, by = tmmx, bs="re")
  #     (6) returns model object, rho, prepared data, and a small param/OR table.
  #
  # Inputs (columns must exist in `data`):
  #   data      : tibble/data.frame with
  #               - date (Date)
  #               - site_no (id)
  #               - daily_max_stream_temp (numeric, °C)
  #               - daily_median_discharge (numeric, CFS)
  #               - tmmx, tmmn (numeric, °C)
  #   threshold : exceedance threshold in °C (default 21)
  #   break_date: Date for period breakpoint (default "2017-01-01")
  #
  # Output:
  #   A list with:
  #     - dat_bam   : modeling data used (with AR_start flag)
  #     - rho_guess : scalar AR(1) correlation used
  #     - fit       : mgcv::bam model object
  #     - params    : data frame of fixed effects with odds ratios and 95% CI
  #
  # Notes:
  #   * Predictors are used on their original scales (no centering) to match your
  #     script’s specification, but you can center upstream if desired.
  #   * AR_start resets correlation when (date_t != date_{t-1} + 1) within site.
  #   * Requires: dplyr, tidyr, mgcv, stats.
  # ------------------------------------------------------------------------------
  # --- sanity checks ----------------------------------------------------------
  need = c("date", "site_no", "daily_max_stream_temp",
           "daily_median_discharge", "tmmx", "tmmn")
  miss = setdiff(need, names(data))
  if (length(miss)) stop("`data` is missing required columns: ",
                         paste(miss, collapse = ", "))
  
  if (!inherits(data$date, "Date")) {
    stop("`date` must be class Date.")
  }
  
  # --- prep: binary outcome, period, dedupe, AR.start -------------------------
  dat_bam = data |>
    dplyr::mutate(
      over21 = as.integer(.data$daily_max_stream_temp > threshold),
      period = factor(dplyr::if_else(.data$date < break_date, glue::glue("pre{year(break_date)}"), glue::glue("post{year(break_date)}"))),
      site_no = as.factor(.data$site_no)
    ) |>
    tidyr::drop_na(over21, daily_median_discharge, tmmx, tmmn, period, site_no) |>
    dplyr::arrange(site_no, date) |>
    dplyr::group_by(site_no, date) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::arrange(site_no, date) |>
    dplyr::group_by(site_no) |>
    dplyr::mutate(AR_start = dplyr::row_number() == 1 | date != dplyr::lag(date) + 1) |>
    dplyr::ungroup() |>
    dplyr::select(date, site_no, over21, daily_median_discharge, tmmx, tmmn, period, AR_start)
  
  # --- rough AR(1) rho guess from baseline GLM --------------------------------
  base_glm = stats::glm(
    over21 ~ daily_median_discharge + tmmx + tmmn + period,
    data = dat_bam,
    family = stats::binomial()
  )
  rp = stats::residuals(base_glm, type = "pearson")
  ac = stats::acf(rp, plot = FALSE, lag.max = 1)$acf
  rho_guess = pmin(pmax(ifelse(length(ac) >= 2, ac[2], 0.4), 0.05), 0.95)
  
  # --- fit BAM with site RE + site-specific slope for tmmx --------------------
  fit = mgcv::bam(
    over21 ~ daily_median_discharge + tmmx + tmmn + period +
      s(site_no, bs = "re") +
      s(site_no, by = tmmx, bs = "re"),
    data    = dat_bam,
    family  = binomial(),
    method  = "fREML",
    discrete = TRUE,
    rho     = rho_guess,
    AR.start = dat_bam$AR_start
  )
  
  # --- tidy fixed effects with ORs --------------------------------------------
  pt = as.data.frame(summary(fit)$p.table)
  pt$term = rownames(pt)
  names(pt) = sub(" ", "_", names(pt), fixed = TRUE)
  
  # build dynamic period term and label
  period_year  = lubridate::year(break_date)
  period_term  = paste0("periodpre", period_year)
  period_label = sprintf("Period: pre-%d vs post-%d", period_year, period_year)
  
  # base mapping + dynamic entry
  map_vec = c(
    "(Intercept)"            = "Intercept",
    "daily_median_discharge" = "Streamflow (per +1 CFS)",
    "tmmx"                   = "Max Air Temperature (per +1 °C)",
    "tmmn"                   = "Min Air Temperature (per +1 °C)"
  )
  map_vec = c(map_vec, stats::setNames(period_label, period_term))  # or rlang::set_names(period_label, period_term)
  
  params = pt |>
    dplyr::rename(
      estimate = Estimate,
      se       = Std._Error,
      z        = z_value,
      p        = `Pr(>|z|)`
    ) |>
    dplyr::mutate(
      OR     = exp(estimate),
      OR_lo  = exp(estimate - 1.96 * se),
      OR_hi  = exp(estimate + 1.96 * se),
      Term   = dplyr::recode(term, !!!map_vec, .default = term)
    ) |>
    dplyr::select(Term, estimate, se, z, p, OR, OR_lo, OR_hi)
  
  # return
  list(
    dataset   = dat_bam,
    rho_guess = rho_guess,
    fit       = fit,
    params    = params
  )
}

interpret_exceedance_gamm = function(fit, scale_flow = 100, break_date = NULL, threshold = 21) {
  # ------------------------------------------------------------------------------
  # Title:    Automated Interpretation of Exceedance GAMM Results (dynamic period)
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-21
  #
  # Description:
  #   Extracts fixed effects from a binomial GAMM/BAM fit and produces a standardized
  #   interpretation paragraph. Discharge effects are reported per +`scale_flow` CFS.
  #   The pre/post "period" term/year is detected from the model (e.g., periodpre2017),
  #   or can be supplied via `break_date`.
  #
  # Args:
  #   fit        : mgcv::bam/gam object (binomial logit) with terms daily_median_discharge,
  #                tmmx, tmmn, and a two-level factor "period" coded as pre/post-YEAR.
  #   scale_flow : numeric, CFS increment to report for discharge (default 100).
  #   break_date : Date, optional; if provided, forces the label to that year.
  #   threshold  : exceedance threshold (°C) used in the model, for text only.
  #
  # Returns:
  #   A single character string with the interpretation paragraph.
  # ------------------------------------------------------------------------------
  
  stopifnot(inherits(fit, "gam") || inherits(fit, "bam"))
  param_tab  = as.data.frame(summary(fit)$p.table)
  smooth_tab = as.data.frame(summary(fit)$s.table)
  
  # Helper to compute OR, CI, p for a named coefficient with optional scaling
  get_OR = function(term, scale = 1) {
    if (!term %in% rownames(param_tab))
      return(list(OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
    est = param_tab[term, "Estimate"]
    se  = param_tab[term, "Std. Error"]
    p   = param_tab[term, "Pr(>|z|)"]
    list(
      OR = exp(est * scale),
      lo = exp((est - 1.96 * se) * scale),
      hi = exp((est + 1.96 * se) * scale),
      p  = p
    )
  }
  
  # --- Identify the period coefficient name and year ---------------------------
  period_term = NULL
  if (!is.null(break_date)) {
    # expect names like "periodpreYYYY"
    y = lubridate::year(break_date)
    candidate = paste0("periodpre", y)
    if (candidate %in% rownames(param_tab)) period_term = candidate
  }
  if (is.null(period_term)) {
    # auto-detect first coefficient that looks like periodpre#### (robust to other factors)
    hits = grep("^periodpre\\d{4}$", rownames(param_tab), value = TRUE)
    if (length(hits) >= 1) period_term = hits[1]
  }
  period_year = if (!is.null(period_term)) as.integer(sub("^periodpre(\\d{4})$", "\\1", period_term)) else NA_integer_
  
  # --- Pull fixed effects ------------------------------------------------------
  flow = get_OR("daily_median_discharge", scale = scale_flow)
  tmx  = get_OR("tmmx")
  tmn  = get_OR("tmmn")
  per  = if (!is.null(period_term)) get_OR(period_term) else list(OR=NA, lo=NA, hi=NA, p=NA)
  
  # --- Random effect significance (if present) ---------------------------------
  re_intercept_p = tryCatch(smooth_tab["s(site_no)", "p-value"],     error = function(e) NA_real_)
  re_slope_p     = tryCatch(smooth_tab["s(site_no):tmmx", "p-value"], error = function(e) NA_real_)
  
  re_text = ""
  if (!is.na(re_intercept_p) || !is.na(re_slope_p)) {
    if (!is.na(re_intercept_p) && re_intercept_p < 0.05 && !is.na(re_slope_p) && re_slope_p < 0.05) {
      re_text = "Site-specific random intercepts and slopes were significant, indicating heterogeneity in both baseline risk and thermal sensitivity across sites. "
    } else if (!is.na(re_intercept_p) && re_intercept_p < 0.05 && (is.na(re_slope_p) || re_slope_p >= 0.05)) {
      re_text = "Random intercepts were significant, capturing baseline differences among sites; random slopes for air temperature were not significant. "
    } else if ((is.na(re_intercept_p) || re_intercept_p >= 0.05) && !is.na(re_slope_p) && re_slope_p < 0.05) {
      re_text = "Random slopes for air temperature were significant, indicating heterogeneous thermal sensitivity; baseline differences among sites were not significant. "
    } else {
      re_text = "Random effects were not significant, suggesting broadly consistent baseline risk and thermal sensitivity across sites. "
    }
  }
  
  # --- Fit summary stats -------------------------------------------------------
  summ = summary(fit)
  r2   = if (!is.null(summ$r.sq)) round(summ$r.sq, 2) else NA_real_
  dev  = if (!is.null(summ$dev.expl)) round(summ$dev.expl * 100, 1) else NA_real_
  
  # --- Formatting helpers ------------------------------------------------------
  fmt_p = function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "< 0.001", sprintf("%.2f", p)))
  pct   = function(x) sprintf("%d", round(100 * x))
  
  # Percent interpretations (guard NA)
  flow_pct  = if (is.na(flow$OR)) "n/a" else pct(1 - flow$OR)
  tmx_pct   = if (is.na(tmx$OR))  "n/a" else pct(tmx$OR - 1)
  tmn_pct   = if (is.na(tmn$OR))  "n/a" else pct(tmn$OR - 1)
  per_OR    = if (is.na(per$OR))  "n/a" else sprintf("%.2f", per$OR)
  
  # --- Build the paragraph -----------------------------------------------------
  period_clause = if (!is.na(period_year)) {
    paste0("The pre-", period_year, " vs post-", period_year, " period effect was ",
           ifelse(!is.na(per$p) && per$p < 0.05, "significant", "not significant"),
           " (OR ≈ ", per_OR, "; p ", fmt_p(per$p), "). ")
  } else {
    ""
  }
  
  ar_text = "An AR(1) residual structure accounted for temporal clustering. "
  fit_text = if (!is.na(dev) && !is.na(r2)) {
    paste0("Overall, the model explained ~", dev, "% of deviance (adj. R² ≈ ", r2, ").")
  } else {
    ""
  }
  
  glue::glue(
    "Exceedance risk (> {threshold} °C) was strongly governed by air temperature and streamflow. ",
    "Each +{scale_flow} CFS increase in discharge lowered the odds by ~{flow_pct}% (p {fmt_p(flow$p)}), ",
    "while a +1 °C increase in maximum and minimum air temperature raised the odds by ~{tmx_pct}% ",
    "and ~{tmn_pct}%, respectively (max p {fmt_p(tmx$p)}; min p {fmt_p(tmn$p)}). ",
    period_clause,
    re_text,
    ar_text,
    fit_text
  ) |> as.character()
}

# helpers for StAGE data specifically
# -- helper to extract "41D 07900" from a header like
# "Water Temp.water temp@41D 07900 (degC)"
.extract_site_from_header = function(colname) {
  colname |>
    stringr::str_extract("(?<=@)[^\\(]+") |>
    stringr::str_squish()
}

# -- cleaner for one tibble
clean_stage_tbl = function(df) {
  # find the temperature column by pattern
  temp_col = names(df)[stringr::str_detect(names(df), "(?i)^\\s*water\\s*temp")]
  flow_col = names(df)[stringr::str_detect(names(df), "(?i)^\\s*discharge")]
  if (length(temp_col) != 1) {
    stop("Could not uniquely identify the temperature column. Found: ",
         paste(temp_col, collapse = ", "))
  }
  site_name = .extract_site_from_header(temp_col)
  
  df |>
    dplyr::mutate(site_no = site_name) |>
    dplyr::rename(temp = tidyselect::all_of(temp_col),
                  discharge = tidyselect::all_of(flow_col)) |>
    dplyr::relocate(site_no, .after = "Timestamp")
}

plot_exceedance_vs_tmmx = function(
    fit, dat_bam,
    mode = c("by_site", "prepost_by_site", "fixed", "fixed_prepost"),
    flow_ref   = NULL,
    tmmn_ref   = NULL,
    period_ref = NULL,
    tmmx_probs = c(0, 1),
    n_pts = 120,
    facet_nrow = 3,
    palette = c(pre2017 = "black", post2017 = "steelblue"),
    title = NULL,
    caption = NULL,
    x_lab = "Daily max air temperature (°C)",
    y_lab = "Probability of >21°C",
    save_path = NULL, width = 10, height = 8, dpi = 300
) {
  # ------------------------------------------------------------------------------
  # Title:    Exceedance Probability vs. Air Max (tmmx) Plotter
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-20 (updated)
  #
  # Key updates:
  #   * "by_site": site-specific medians for flow/tmmn; EXCLUDE random effects;
  #                curves limited to each site's observed tmmx range.
  #   * "prepost_by_site": INCLUDE random effects; curves limited to each site's
  #                        observed tmmx range.
  #
  # Other modes:
  #   * "fixed": single fixed-effects curve (exclude RE)
  #   * "fixed_prepost": fixed-effects curves for pre vs post (exclude RE)
  # ------------------------------------------------------------------------------
  
  mode = match.arg(mode)
  
  # sanity checks
  stopifnot(all(c("site_no","tmmx","tmmn","period","daily_median_discharge") %in% names(dat_bam)))
  if (!is.factor(dat_bam$site_no)) dat_bam$site_no = factor(dat_bam$site_no)
  if (!is.factor(dat_bam$period))  dat_bam$period  = factor(dat_bam$period)
  
  # reference values (global) used in some modes
  if (is.null(flow_ref))  flow_ref  = stats::median(dat_bam$daily_median_discharge, na.rm = TRUE)
  if (is.null(tmmn_ref))  tmmn_ref  = stats::median(dat_bam$tmmn, na.rm = TRUE)
  if (is.null(period_ref)) period_ref = names(sort(table(dat_bam$period), decreasing = TRUE))[1]
  period_ref = factor(period_ref, levels = levels(dat_bam$period))
  
  site_levels = levels(dat_bam$site_no)
  
  # helper: predict response with optional exclusion of random effects
  .pred_resp = function(newdata, exclude = NULL) {
    as.numeric(predict(fit, newdata = newdata, type = "response", exclude = exclude))
  }
  
  if (mode %in% c("by_site", "prepost_by_site")) {
    
    # -------- site-specific ranges & medians ----------
    ranges <- dat_bam |>
      dplyr::group_by(site_no) |>
      dplyr::summarise(
        t_lo = as.numeric(stats::quantile(tmmx, tmmx_probs[1], na.rm = TRUE)),
        t_hi = as.numeric(stats::quantile(tmmx, tmmx_probs[2], na.rm = TRUE)),
        flow_ref_site = stats::median(daily_median_discharge, na.rm = TRUE),
        tmmn_ref_site = stats::median(tmmn, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        t_seq = purrr::map2(t_lo, t_hi, ~ seq(.x, .y, length.out = n_pts))
      )
    
    if (mode == "by_site") {
      # ---- site medians; exclude REs; site-specific x range ----
      new_df <- ranges |>
        tidyr::unnest(t_seq) |>
        dplyr::transmute(
          site_no,
          tmmx = t_seq,
          daily_median_discharge = flow_ref_site,
          tmmn   = tmmn_ref_site,
          period = period_ref,
          AR_start = FALSE
        )
      
      new_df$prob <- .pred_resp(new_df, exclude = c("s(site_no)", "s(site_no):tmmx"))
      
      if (is.null(title))   title   = "Exceedance Probability (>21°C) vs Daily Max Air Temperature"
      if (is.null(caption)) caption = stringr::str_wrap(
        "Curves show fixed-effects predictions by site with streamflow and minimum air temperature held at each site's medians. Random effects are excluded, so differences among panels reflect covariate medians. Each facet's curve is limited to that site's observed air-temperature range.",
        width = 140
      )
      
      p <- ggplot2::ggplot(new_df, ggplot2::aes(x = tmmx, y = prob)) +
        ggplot2::geom_line(linewidth = 1, color = "steelblue") +
        ggplot2::facet_wrap(~ site_no, scales = "fixed", nrow = facet_nrow) +
        ggplot2::labs(title = title, x = x_lab, y = y_lab, caption = caption) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = 11),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(hjust = 0)
        ) +
        ggplot2::ylim(0,1)
      
    } else { # prepost_by_site
      # ---- global medians for flow/tmmn; include REs; site-specific x range ----
      new_df <- ranges |>
        tidyr::unnest(t_seq) |>
        tidyr::expand_grid(period = factor(c("pre2017","post2017"),
                                           levels = levels(dat_bam$period))) |>
        dplyr::transmute(
          site_no,
          tmmx = t_seq,
          period,
          daily_median_discharge = flow_ref,
          tmmn = tmmn_ref,
          AR_start = FALSE
        )
      
      new_df$prob <- .pred_resp(new_df)  # includes random effects
      
      if (is.null(title))   title   = "Exceedance Probability (>21°C) vs Air Max: Pre vs Post 2017"
      if (is.null(caption)) caption = stringr::str_wrap(
        "Predicted exceedance probability comparing pre- vs post-2017 by site. Lines include site random effects and hold streamflow and minimum air temperature at global medians. Each facet's curves are limited to that site's observed air-temperature range.",
        width = 140
      )
      
      p <- ggplot2::ggplot(new_df, ggplot2::aes(x = tmmx, y = prob, color = period)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_wrap(~ site_no, scales = "fixed", nrow = facet_nrow) +
        ggplot2::scale_color_manual(values = palette, name = "Period") +
        ggplot2::labs(title = title, x = x_lab, y = y_lab, caption = caption) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold", size = 11),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          legend.position = "top",
          plot.caption = ggplot2::element_text(hjust = 0)
        ) +
        ggplot2::ylim(0,1)
    }
    
  } else {
    # ---------- Single-panel fixed-effects curves (exclude random terms) ----------
    exclude_terms = c("s(site_no)", "s(site_no):tmmx")
    dummy_site    = site_levels[1]  # for factor integrity only
    
    # global tmmx sequence
    tmmx_lo  = as.numeric(stats::quantile(dat_bam$tmmx, tmmx_probs[1], na.rm = TRUE))
    tmmx_hi  = as.numeric(stats::quantile(dat_bam$tmmx, tmmx_probs[2], na.rm = TRUE))
    tmmx_seq = seq(tmmx_lo, tmmx_hi, length.out = n_pts)
    
    if (mode == "fixed") {
      new_df = tibble::tibble(
        site_no = factor(dummy_site, levels = site_levels),
        tmmx = tmmx_seq,
        daily_median_discharge = flow_ref,
        tmmn = tmmn_ref,
        period = period_ref,
        AR_start = FALSE
      )
      new_df$prob = .pred_resp(new_df, exclude = exclude_terms)
      
      if (is.null(title))   title   = "Fixed-effects Prediction: Pr(>21°C) vs Air Max"
      if (is.null(caption)) caption = stringr::str_wrap(
        "Fixed-effects prediction of exceedance probability >21 °C versus maximum air temperature. Random intercepts and slopes are excluded, so the curve reflects the average relationship across all sites. Streamflow and minimum air temperature are fixed at global medians.",
        width = 140
      )
      
      p <- ggplot2::ggplot(new_df, ggplot2::aes(x = tmmx, y = prob)) +
        ggplot2::geom_line(linewidth = 1, color = "black") +
        ggplot2::labs(title = title, x = x_lab, y = y_lab, caption = caption) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(hjust = 0)
        ) +
        ggplot2::ylim(0,1)
      
    } else { # fixed_prepost
      new_df <- tidyr::expand_grid(
        tmmx = tmmx_seq,
        period = factor(c("pre2017","post2017"), levels = levels(dat_bam$period))
      ) |>
        dplyr::mutate(
          site_no = factor(dummy_site, levels = site_levels),
          daily_median_discharge = flow_ref,
          tmmn = tmmn_ref,
          AR_start = FALSE
        )
      new_df$prob <- .pred_resp(new_df, exclude = exclude_terms)
      
      if (is.null(title))   title   = "Fixed-effects Prediction: Pre vs Post 2017"
      if (is.null(caption)) caption = stringr::str_wrap(
        "Fixed-effects predictions of exceedance probability >21 °C versus maximum air temperature, contrasting pre- vs post-2017. Random effects are excluded to highlight the overall period effect. Streamflow and minimum air temperature fixed at global medians.",
        width = 140
      )
      
      p <- ggplot2::ggplot(new_df, ggplot2::aes(x = tmmx, y = prob, color = period)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_color_manual(values = palette, name = "Period") +
        ggplot2::labs(title = title, x = x_lab, y = y_lab, caption = caption) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          legend.position = "top",
          plot.caption = ggplot2::element_text(hjust = 0)
        ) +
        ggplot2::ylim(0,1)
    }
  }
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  p
}

summarize_exceedance_gamm = function(fit, dataset,
                                      tmmx_vals = c(15, 20, 25),
                                      delta_temps = c(5, 10),
                                      flow_scale = 100) {
  # ------------------------------------------------------------------------------
  # Title:    Summarize Exceedance GAMM Effects (ORs and Probabilities)
  # Author:   Zachary H. Hoylman
  # Date:     8-21-2025
  #
  # Description:
  #   Given a fitted binomial GAMM (mgcv::bam) of exceedance probability,
  #   this function:
  #     * extracts odds ratios per +1 °C (tmmx, tmmn) and per +100 CFS (flow),
  #     * computes odds ratios for multi-degree increments (e.g., +5°C, +10°C),
  #     * predicts probabilities at representative values of predictors,
  #       holding other covariates at their medians / most common level.
  #
  # Inputs:
  #   fit          : mgcv::bam object (binomial link)
  #   dataset      : the dataset used to fit (must include site_no, tmmx, tmmn,
  #                  daily_median_discharge, period)
  #   tmmx_vals    : numeric vector of tmmx values to compute probabilities at
  #   delta_temps  : numeric vector of increments (°C) to compute ORs for
  #   flow_scale   : numeric; scaling for flow OR (default 100 CFS)
  #
  # Output:
  #   list with:
  #     - OR_table : odds ratios (per-unit and multi-unit increments)
  #     - probs    : predicted exceedance probabilities at ref conditions
  # ------------------------------------------------------------------------------
  stopifnot(inherits(fit, "bam"), "tmmx" %in% names(dataset))
  
  # Extract parametric coefficients
  coefs = coef(fit)
  
  # --- ORs for main effects ---
  OR_flow = exp(coefs["daily_median_discharge"] * flow_scale)
  OR_tmmx = exp(coefs["tmmx"])
  OR_tmmn = exp(coefs["tmmn"])
  
  OR_multi = tibble::tibble(
    term = c(paste0("+", flow_scale, " CFS flow"),
             paste0("+1 °C max air temp"),
             paste0("+1 °C min air temp")),
    OR   = c(OR_flow, OR_tmmx, OR_tmmn)
  )
  
  # Add multi-degree increments
  for (d in delta_temps) {
    OR_multi = dplyr::bind_rows(
      OR_multi,
      tibble::tibble(
        term = c(paste0("+", d, " °C max air temp"),
                 paste0("+", d, " °C min air temp")),
        OR   = c(OR_tmmx^d, OR_tmmn^d)
      )
    )
  }
  
  # --- Representative values ---
  flow_ref   = median(dataset$daily_median_discharge, na.rm = TRUE)
  tmmn_ref   = median(dataset$tmmn, na.rm = TRUE)
  period_ref = names(sort(table(dataset$period), decreasing = TRUE))[1]
  site_ref   = levels(dataset$site_no)[1]  # dummy level for fixed effects
  
  # Build newdata for fixed-effects predictions (exclude REs)
  new_df = expand.grid(
    site_no = factor(site_ref, levels = levels(dataset$site_no)),
    tmmx = tmmx_vals,
    daily_median_discharge = flow_ref,
    tmmn   = tmmn_ref,
    period = factor(period_ref, levels = levels(dataset$period)),
    AR_start = FALSE
  )
  
  new_df$prob = as.numeric(predict(
    fit, newdata = new_df, type = "response",
    exclude = c("s(site_no)", "s(site_no):tmmx")
  ))
  
  probs = tibble::as_tibble(new_df) |>
    dplyr::select(tmmx, prob)
  
  list(
    OR_table = OR_multi,
    probs    = probs
  )
}

plot_airtemp_trends = function(
    daily_temps,
    months      = 5:10,          # restrict to May–Oct by default
    end_year    = 2024,          # last full year to include
    palette_fun = viridis::turbo,
    save_path   = NULL,
    width = 10, height = 7, dpi = 300
) {
  # ------------------------------------------------------------------------------
  # Title:    Annual Median Air Temperature Trends (gridMET) with Site-specific Significance
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-21
  #
  # Description:
  #   Computes annual (May–Oct) median daily max air temperature (tmmx) per site,
  #   orders sites by their long-term average, and plots site-colored time series
  #   with linear fits only if slope ≠ 0 (p < 0.05). A pooled regression with site
  #   fixed effects estimates the average warming rate across sites in °C per decade.
  #   Caption text reflects pooled slope and whether all site-specific slopes are significant.
  # ------------------------------------------------------------------------------
  
  stopifnot(all(c("site_no","date","tmmx") %in% names(daily_temps)))
  
  # --- Ensure date column is proper Date ---
  if (!inherits(daily_temps$date, "Date")) {
    daily_temps = daily_temps |>
      dplyr::mutate(date = as.Date(date))
  }
  
  # --- Annual medians May–Oct ---
  climate_trends =
    daily_temps |>
    dplyr::filter(lubridate::month(date) %in% months,
                  lubridate::year(date) <= end_year) |>
    dplyr::group_by(site_no, year = lubridate::year(date)) |>
    dplyr::summarise(median_max = median(tmmx, na.rm = TRUE), .groups = "drop")
  
  # --- Site order (warmest to coolest long-term) ---
  site_order =
    climate_trends |>
    dplyr::group_by(site_no) |>
    dplyr::summarise(avg_med = mean(median_max, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(avg_med)) |>
    dplyr::pull(site_no)
  climate_trends$site_no = factor(climate_trends$site_no, levels = site_order)
  
  # --- Pooled warming rate (°C per decade) ---
  pooled_lm = stats::lm(median_max ~ year + site_no, data = climate_trends)
  slope_decade = unname(stats::coef(pooled_lm)["year"] * 10)
  ci_decade    = stats::confint(pooled_lm)["year", ] * 10
  slope_txt = sprintf("%.2f", slope_decade)
  ci_lo_txt  = sprintf("%.2f", ci_decade[1])
  ci_hi_txt  = sprintf("%.2f", ci_decade[2])
  
  # --- Per-site regressions for significance check ---
  site_slopes =
    climate_trends |>
    dplyr::group_by(site_no) |>
    dplyr::group_map(~{
      mod = lm(median_max ~ year, data = .x)
      broom::tidy(mod) |>
        dplyr::filter(term == "year") |>
        dplyr::mutate(
          slope_decade = estimate * 10,
          site_no = unique(.y$site_no)
        )
    }) |>
    dplyr::bind_rows()
  
  # Which sites significant?
  site_slopes = site_slopes |>
    dplyr::mutate(sig = p.value < 0.05)
  
  all_sig = all(site_slopes$sig)
  
  # --- Caption text ---
  caption_txt = stringr::str_wrap(
    paste0(
      "Lines show annual May–Oct medians. ",
      "Pooled warming across sites: ",
      slope_txt, " °C/decade (95% CI ", ci_lo_txt, " to ", ci_hi_txt, "). ",
      if (all_sig) {
        "All site-specific warming trends are significant."
      } else {
        "Regression lines shown only for sites with significant (p < 0.05) warming trends."
      }
    ),
    width = 120   # tweak width until it looks good
  )
  
  # --- Palette ---
  n_sites = length(levels(climate_trends$site_no))
  cols    = palette_fun(max(n_sites, 3)) |> rev()
  
  # --- Plot ---
  p =
    ggplot2::ggplot(climate_trends, ggplot2::aes(x = year, y = median_max, color = site_no, group = site_no)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.3) +
    {
      if (all_sig) {
        ggplot2::geom_smooth(method = "lm", se = FALSE)
      } else {
        sig_sites = site_slopes$site_no[site_slopes$sig]
        ggplot2::geom_smooth(
          data = climate_trends |> dplyr::filter(site_no %in% sig_sites),
          method = "lm", se = FALSE
        )
      }
    } +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title    = "Annual Median Daily Max Air Temperature (gridMET)",
      subtitle = "May–October medians by site in the Big Hole watershed",
      x = "Year",
      y = "Temperature (°C)",
      caption = caption_txt
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5),
      plot.caption     = ggplot2::element_text(hjust = 0)
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  # Attach slope info
  attr(p, "pooled_slope_C_per_decade") = slope_decade
  attr(p, "pooled_slope_CI_C_per_decade") = ci_decade
  attr(p, "site_slopes") = site_slopes
  return(p)
}

plot_streamflow_trends = function(
    merged,
    months      = 5:10,
    end_year    = 2024,
    palette_fun = viridis::turbo,
    save_path   = NULL,
    width = 10, height = 7, dpi = 300
) {
  # ------------------------------------------------------------------------------
  # Title:    Annual Median Streamflow Trends (May–Oct) with Site-specific Significance
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-21
  #
  # Description:
  #   Computes annual May–Oct median daily streamflow (CFS) per site from merged
  #   USGS and StAGE records. Sites are ordered by their long-term average flows
  #   (highest to lowest). A pooled linear regression with site fixed effects is
  #   used to estimate the average decadal trend in streamflow (CFS/decade).
  #
  #   Site-specific regressions are also fit. Trend lines are drawn only for sites
  #   with statistically significant slopes (p < 0.05). If all site trends are
  #   significant, this is noted in the caption. Caption also reports the pooled
  #   slope and its 95% confidence interval.
  #
  # Inputs:
  #   merged      : tibble/data.frame with columns
  #                   - site_no (character/factor)
  #                   - date (Date)
  #                   - daily_median_discharge (numeric; CFS)
  #   months      : integer vector of months to include (default May–Oct)
  #   end_year    : last year to include in the analysis
  #   palette_fun : function(n) -> vector of n colors (default viridis::turbo)
  #   save_path   : optional file path to save figure
  #   width,height,dpi : graphics device parameters for saving
  #
  # Output:
  #   ggplot object of site-colored annual median flows with trend lines.
  #   Attributes attached:
  #     - pooled_slope_CFS_per_decade
  #     - pooled_slope_CI_CFS_per_decade
  #     - site_slopes (per-site slope and p-value)
  #
  # Notes:
  #   * Annual median = median of daily median discharge between May–Oct.
  #   * Trend lines drawn only if p < 0.05 for site regressions.
  #   * Caption dynamically reflects pooled slope, confidence interval,
  #     and whether all site-specific slopes are significant.
  # ------------------------------------------------------------------------------
  stopifnot(all(c("site_no","date","daily_median_discharge") %in% names(merged)))
  if (!inherits(merged$date, "Date")) {
    merged = merged |> dplyr::mutate(date = as.Date(date))
  }
  
  flow_trends =
    merged |>
    dplyr::filter(lubridate::month(date) %in% months,
                  lubridate::year(date) <= end_year) |>
    dplyr::group_by(site_no, year = lubridate::year(date)) |>
    dplyr::summarise(
      median_flow = median(daily_median_discharge, na.rm = TRUE),
      .groups = "drop"
    )
  
  site_order =
    flow_trends |>
    dplyr::group_by(site_no) |>
    dplyr::summarise(avg_med = mean(median_flow, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(avg_med)) |>
    dplyr::pull(site_no)
  flow_trends$site_no = factor(flow_trends$site_no, levels = site_order)
  
  pooled_lm   = stats::lm(median_flow ~ year + site_no, data = flow_trends)
  slope_dec   = unname(stats::coef(pooled_lm)["year"] * 10)
  ci_dec      = stats::confint(pooled_lm)["year", ] * 10
  slope_txt   = sprintf("%.1f", slope_dec)
  ci_lo_txt   = sprintf("%.1f", ci_dec[1])
  ci_hi_txt   = sprintf("%.1f", ci_dec[2])
  
  # --- FIXED: use group_map() and add site_no afterward
  site_slopes =
    flow_trends |>
    dplyr::group_by(site_no) |>
    dplyr::group_map(~{
      mod  = stats::lm(median_flow ~ year, data = .x)
      tt   = broom::tidy(mod)
      out  = tt[tt$term == "year", , drop = FALSE]
      dplyr::tibble(
        slope_CFS_y      = out$estimate,
        slope_CFS_decade = out$estimate * 10,
        p_value          = out$p.value
      )
    }) |>
    dplyr::bind_rows(.id = "site_no") |>
    dplyr::mutate(
      site_no = levels(flow_trends$site_no)[as.integer(site_no)],
      sig = p_value < 0.05
    )
  
  all_sig = all(site_slopes$sig)
  
  caption_txt = stringr::str_wrap(
    paste0(
      "Lines show annual May–Oct medians. ",
      "Pooled trend across sites: ", slope_txt, " CFS/decade (95% CI ",
      ci_lo_txt, " to ", ci_hi_txt, "). ",
      if (all_sig) "All site-specific trends are significant (p < 0.05)." else
        "Regression lines are drawn only for sites with significant (p < 0.05) trends."
    ),
    width = 110
  )
  
  n_sites = length(levels(flow_trends$site_no))
  cols    = palette_fun(max(n_sites, 3)) |> rev()
  
  p =
    ggplot2::ggplot(flow_trends,
                    ggplot2::aes(x = year, y = median_flow, color = site_no, group = site_no)
    ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.35) +
    {
      if (all_sig) {
        ggplot2::geom_smooth(method = "lm", se = FALSE)
      } else {
        sig_sites = site_slopes$site_no[site_slopes$sig]
        ggplot2::geom_smooth(
          data = flow_trends |> dplyr::filter(site_no %in% sig_sites),
          method = "lm", se = FALSE
        )
      }
    } +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title    = "Annual Median Streamflow (May–Oct)",
      subtitle = "Daily median discharge aggregated to seasonal medians by site",
      x = "Year", y = "Discharge (CFS)", caption = caption_txt
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.caption  = ggplot2::element_text(hjust = 0),
      legend.title  = ggplot2::element_text(face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  attr(p, "pooled_slope_CFS_per_decade")     = slope_dec
  attr(p, "pooled_slope_CI_CFS_per_decade")  = ci_dec
  attr(p, "site_slopes")                     = site_slopes
  return(p)
}

plot_exceedance_vs_airtemp_from_bases = function(
    daily_temps,
    exceedance_hours,
    response    = c("cumulative_hours_21","cumulative_hours_25"),
    months      = 5:10,          # May–Oct by default
    end_year    = 2024,          # last full year to include
    palette_fun = viridis::turbo,
    save_path   = NULL,
    width = 10, height = 7, dpi = 300
) {
  # ------------------------------------------------------------------------------
  # Title:    Exceedance vs Air Temperature (build air_stream internally)
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-21
  #
  # Description:
  #   Computes annual (May–Oct) median daily max air temperature (tmmx) per site
  #   from `daily_temps`, joins to `exceedance_hours` to create `air_stream`,
  #   then plots seasonal exceedance hours vs seasonal median air temp by site.
  #   Per-site linear fits are drawn only if significant (p < 0.05).
  #   Also reports a pooled slope (site fixed effects) in the caption.
  #
  # Inputs:
  #   daily_temps      : tibble with columns site_no, date (Date), tmmx (°C)
  #   exceedance_hours : tibble with site_no, year, cumulative_hours_21/25
  #   response         : "cumulative_hours_21" or "cumulative_hours_25"
  #   months           : months to include (default 5:10)
  #   end_year         : last year to include (default 2024)
  # ------------------------------------------------------------------------------
  
  response = match.arg(response)
  
  # --- sanity checks -----------------------------------------------------------
  need_daily = c("site_no","date","tmmx")
  miss_daily = setdiff(need_daily, names(daily_temps))
  if (length(miss_daily)) stop("`daily_temps` missing: ", paste(miss_daily, collapse=", "))
  
  need_ex = c("site_no","year", "cumulative_hours_21","cumulative_hours_25")
  miss_ex = setdiff(need_ex, names(exceedance_hours))
  if (length(miss_ex)) stop("`exceedance_hours` missing: ", paste(miss_ex, collapse=", "))
  
  # Ensure date class
  if (!inherits(daily_temps$date, "Date")) {
    daily_temps = dplyr::mutate(daily_temps, date = as.Date(date))
  }
  
  # --- compute annual_temps as specified --------------------------------------
  annual_temps =
    daily_temps |>
    dplyr::filter(lubridate::month(date) %in% months,
                  lubridate::year(date) <= end_year) |>
    dplyr::group_by(site_no, year = lubridate::year(date)) |>
    dplyr::summarise(median_max = median(tmmx, na.rm = TRUE), .groups = "drop")
  
  # --- build air_stream via left_join -----------------------------------------
  air_stream =
    exceedance_hours |>
    dplyr::filter(year <= end_year) |>
    dplyr::left_join(annual_temps, by = c("site_no","year"))
  
  # drop rows without the needed pieces
  df = air_stream |>
    dplyr::select(site_no, year, y = !!rlang::sym(response), x = median_max) |>
    tidyr::drop_na(y, x)
  
  # --- order sites by long-term exceedance level ------------------------------
  site_order =
    df |>
    dplyr::group_by(site_no) |>
    dplyr::summarise(avg_y = mean(y, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(avg_y)) |>
    dplyr::pull(site_no)
  df$site_no = factor(df$site_no, levels = site_order)
  
  # --- pooled regression with site fixed effects ------------------------------
  pooled_lm = stats::lm(y ~ x + site_no, data = df)
  slope     = unname(stats::coef(pooled_lm)["x"])
  ci        = stats::confint(pooled_lm)["x", ]
  
  slope_txt = sprintf("%.2f", slope)
  ci_lo_txt = sprintf("%.2f", ci[1])
  ci_hi_txt = sprintf("%.2f", ci[2])
  
  # --- per-site regressions; only draw if significant -------------------------
  site_slopes =
    df |>
    dplyr::group_by(site_no) |>
    dplyr::group_map(~{
      mod = stats::lm(y ~ x, data = .x)
      broom::tidy(mod) |>
        dplyr::filter(term == "x") |>
        dplyr::mutate(site_no = unique(.x$site_no))
    }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(sig = p.value < 0.05)
  
  all_sig = all(site_slopes$sig)
  
  # --- caption ----------------------------------------------------------------
  cap = stringr::str_wrap(
    paste0(
      "Seasonal exceedance hours (≥21 °C) vs median May–Oct air temperature. ",
      "Pooled slope across sites: ", slope_txt,
      " hours per +1 °C (95% CI ", ci_lo_txt, " to ", ci_hi_txt, "). ",
      if (all_sig) "All site-specific regressions are significant."
      else "Regression lines shown only for sites with significant (p < 0.05) relationships."
    ),
    width = 110
  )
  
  # --- palette ----------------------------------------------------------------
  n_sites = length(levels(df$site_no))
  cols    = palette_fun(max(n_sites, 3)) |> rev()
  
  # --- plot -------------------------------------------------------------------
  p = ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = site_no)) +
    ggplot2::geom_point(alpha = 0.55) +
    {
      if (all_sig) {
        ggplot2::geom_smooth(method = "lm", se = FALSE)
      } else {
        ggplot2::geom_smooth(
          data = df |> dplyr::filter(site_no %in% site_slopes$site_no[site_slopes$sig]),
          method = "lm", se = FALSE
        )
      }
    } +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title    = "Exceedance vs Seasonal Air Temperature",
      subtitle = paste0("Response: ", response, " | Driver: Median May–Oct tmmx"),
      x = "Median May–Oct Daily Max Air Temperature (°C)",
      y = paste0("Seasonal ", response),
      caption = cap
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.caption  = ggplot2::element_text(hjust = 0)
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  attr(p, "pooled_slope_hours_per_C")    = slope
  attr(p, "pooled_slope_CI_hours_per_C") = ci
  attr(p, "site_slopes")                 = site_slopes
  attr(p, "air_stream")                  = air_stream
  return(p)
}

analyze_seasonal_exceedance = function(
    daily_temps,
    daily_max_complete,
    exceedance_hours,
    response      = c("cumulative_hours_21","cumulative_hours_25"),
    months        = 5:10,     # May–Oct by default
    end_year      = 2024,
    hold_at       = c("overall_median","site_median"), # partials reference
    palette_fun   = viridis::turbo,
    save_prefix   = NULL,     # e.g., "~/bighole-hydro/figs/seasonal_exceedance"
    caption_width = 90,      # <<< wrap width for captions
    width = 8, height = 8, dpi = 300
){
  # ------------------------------------------------------------------------------
  # Title:  Seasonal Exceedance vs. Climate & Flow (Annual LM + Partial Effects)
  # Author: Zachary H. Hoylman
  # Date:   2025-08-21
  #
  # Description:
  #   Uses seasonal (May–Oct) medians of daily max air temperature (gridMET tmmx)
  #   and daily median streamflow to explain interannual variability in seasonal
  #   exceedance hours. Fits lm:
  #       cumulative_hours_XX ~ median_max_temp + median_flow
  #   and returns partial-effects plots (with interpretation as plot captions)
  #   for temperature (holding flow fixed) and for flow (holding temp fixed).
  # ------------------------------------------------------------------------------
  
  response = match.arg(response)
  hold_at  = match.arg(hold_at)
  
  stopifnot(all(c("site_no","date","tmmx") %in% names(daily_temps)))
  stopifnot(all(c("site_no","date","daily_median_discharge") %in% names(daily_max_complete)))
  stopifnot(all(c("site_no","year", response) %in% names(exceedance_hours)))
  
  if (!inherits(daily_temps$date, "Date"))         daily_temps$date = as.Date(daily_temps$date)
  if (!inherits(daily_max_complete$date, "Date"))  daily_max_complete$date = as.Date(daily_max_complete$date)
  
  # Seasonal medians per site-year
  annual_temps =
    daily_temps |>
    dplyr::filter(lubridate::month(date) %in% months,
                  lubridate::year(date) <= end_year) |>
    dplyr::group_by(site_no, year = lubridate::year(date)) |>
    dplyr::summarise(median_max_temp = median(tmmx, na.rm = TRUE), .groups = "drop")
  
  annual_flow  =
    daily_max_complete |>
    dplyr::filter(lubridate::month(date) %in% months,
                  lubridate::year(date) <= end_year) |>
    dplyr::group_by(site_no, year = lubridate::year(date)) |>
    dplyr::summarise(median_flow = median(daily_median_discharge, na.rm = TRUE), .groups = "drop")
  
  air_stream_flow =
    exceedance_hours |>
    dplyr::left_join(annual_temps, by = c("site_no","year")) |>
    dplyr::left_join(annual_flow,  by = c("site_no","year")) |>
    dplyr::filter(!is.na(.data[[response]]),
                  !is.na(median_max_temp),
                  !is.na(median_flow))
  
  # Model
  fml = stats::as.formula(paste(response, "~ median_max_temp + median_flow"))
  fit = stats::lm(fml, data = air_stream_flow)
  
  tb  = broom::tidy(fit)
  gl  = broom::glance(fit)
  
  # Coefs for interpretation
  b_t   = tb$estimate[tb$term == "median_max_temp"]
  p_t   = tb$p.value[tb$term == "median_max_temp"]
  b_q   = tb$estimate[tb$term == "median_flow"]
  p_q   = tb$p.value[tb$term == "median_flow"]
  r2    = gl$r.squared
  r2adj = gl$adj.r.squared
  per100 = 100 * b_q  # hours per +100 CFS
  
  # Prediction grids (overall/site)
  if (hold_at == "overall_median") {
    flow_ref = stats::median(air_stream_flow$median_flow, na.rm = TRUE)
    tmmx_ref = stats::median(air_stream_flow$median_max_temp, na.rm = TRUE)
    
    temp_grid = data.frame(
      median_max_temp = seq(min(air_stream_flow$median_max_temp, na.rm = TRUE),
                            max(air_stream_flow$median_max_temp, na.rm = TRUE),
                            length.out = 200),
      median_flow = flow_ref
    )
    flow_grid = data.frame(
      median_max_temp = tmmx_ref,
      median_flow = seq(min(air_stream_flow$median_flow, na.rm = TRUE),
                        max(air_stream_flow$median_flow, na.rm = TRUE),
                        length.out = 200)
    )
    
    pr_t = stats::predict(fit, newdata = temp_grid, se.fit = TRUE)
    pred_temp = cbind(temp_grid, .fitted = as.numeric(pr_t$fit), .se.fit = as.numeric(pr_t$se.fit)) |>
      dplyr::mutate(lwr = .fitted - 1.96*.se.fit, upr = .fitted + 1.96*.se.fit)
    
    pr_q = stats::predict(fit, newdata = flow_grid, se.fit = TRUE)
    pred_flow = cbind(flow_grid, .fitted = as.numeric(pr_q$fit), .se.fit = as.numeric(pr_q$se.fit)) |>
      dplyr::mutate(lwr = .fitted - 1.96*.se.fit, upr = .fitted + 1.96*.se.fit)
    
    pred_temp$site_no = NA
    pred_flow$site_no = NA
    
  } else {
    # per-site clamped ranges
    ranges = air_stream_flow |>
      dplyr::group_by(site_no) |>
      dplyr::summarise(
        t_lo = min(median_max_temp, na.rm = TRUE),
        t_hi = max(median_max_temp, na.rm = TRUE),
        q_lo = min(median_flow,      na.rm = TRUE),
        q_hi = max(median_flow,      na.rm = TRUE),
        t_ref = stats::median(median_max_temp, na.rm = TRUE),
        q_ref = stats::median(median_flow,      na.rm = TRUE),
        .groups = "drop"
      )
    
    temp_grid = ranges |>
      tidyr::unnest(tibble::tibble(
        median_max_temp = seq(t_lo, t_hi, length.out = pmax(2L, ceiling((t_hi - t_lo) * 20)))
      )) |>
      dplyr::transmute(site_no, median_max_temp, median_flow = q_ref)
    
    flow_grid = ranges |>
      tidyr::unnest(tibble::tibble(
        median_flow = seq(q_lo, q_hi, length.out = pmax(2L, ceiling((q_hi - q_lo) * 0.4)))
      )) |>
      dplyr::transmute(site_no, median_max_temp = t_ref, median_flow)
    
    pr_t = stats::predict(fit, newdata = as.data.frame(temp_grid), se.fit = TRUE)
    pred_temp = cbind(temp_grid, .fitted = as.numeric(pr_t$fit), .se.fit = as.numeric(pr_t$se.fit)) |>
      dplyr::mutate(lwr = .fitted - 1.96*.se.fit, upr = .fitted + 1.96*.se.fit)
    
    pr_q = stats::predict(fit, newdata = as.data.frame(flow_grid), se.fit = TRUE)
    pred_flow = cbind(flow_grid, .fitted = as.numeric(pr_q$fit), .se.fit = as.numeric(pr_q$se.fit)) |>
      dplyr::mutate(lwr = .fitted - 1.96*.se.fit, upr = .fitted + 1.96*.se.fit)
  }
  
  # Interpretation (to embed as plot caption)
  interp = glue::glue(
    "Across sites and years, seasonal exceedance hours are strongly coupled to maximum temperature and streamflow. ",
    "Holding streamflow constant, each +1 °C increase in median May–Oct daily max air temperature is associated with ~{round(b_t,1)} additional hours ",
    if (response == 'cumulative_hours_21') "≥21 °C" else "≥25 °C",
    " (p = {ifelse(p_t < 0.001, '<0.001', sprintf('%.3f', p_t))}). ",
    "Holding air temperature constant, each +100 CFS increase in seasonal median flow is associated with ~{round(per100,1)} fewer hours ",
    if (response == 'cumulative_hours_21') "≥21 °C" else "≥25 °C",
    " (p = {ifelse(p_q < 0.001, '<0.001', sprintf('%.3f', p_q))}). ",
    "The model explains R² = {sprintf('%.2f', r2)}, adj. R² = {sprintf('%.2f', r2adj)}."
  )
  caption_text = stringr::str_wrap(as.character(interp), width = caption_width)
  
  # Colors
  n_sites = dplyr::n_distinct(air_stream_flow$site_no)
  cols    = palette_fun(max(n_sites, 3))
  
  # Plot: Temperature partial
  p_temp = ggplot2::ggplot(
    air_stream_flow,
    ggplot2::aes(x = median_max_temp, y = .data[[response]], color = site_no)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title = "Marginal Effect of Air Temperature",
      subtitle = "Median May–Oct Daily Max Air Temperature (holding streamflow constant)",
      x = "Median May–Oct Daily Max Air Temperature (°C)",
      y = if (response == "cumulative_hours_21") "Seasonal Hours ≥ 21 °C" else "Seasonal Hours ≥ 25 °C",
      caption = caption_text
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  if ("site_no" %in% names(pred_temp) && any(!is.na(pred_temp$site_no))) {
    p_temp = p_temp +
      ggplot2::geom_ribbon(
        data = pred_temp,
        ggplot2::aes(x = median_max_temp, ymin = lwr, ymax = upr),
        inherit.aes = FALSE, alpha = 0.12, fill = "grey70"
      ) +
      ggplot2::geom_line(
        data = pred_temp,
        ggplot2::aes(x = median_max_temp, y = .fitted, group = site_no),
        inherit.aes = FALSE, color = "black", linewidth = 0.9
      )
  } else {
    p_temp = p_temp +
      ggplot2::geom_line(
        data = pred_temp,
        ggplot2::aes(x = median_max_temp, y = .fitted),
        inherit.aes = FALSE, color = "black", linewidth = 1
      ) +
      ylim(0,800)
  }
  
  # Plot: Flow partial
  p_flow = ggplot2::ggplot(
    air_stream_flow,
    ggplot2::aes(x = median_flow, y = .data[[response]], color = site_no)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title = "Marginal Effect of Streamflow",
      subtitle = "Seasonal exceedance hours vs. median May–Oct streamflow (holding air temperature constant)",
      x = "Median May–Oct streamflow (CFS)",
      y = if (response == "cumulative_hours_21") "Seasonal hours ≥ 21 °C" else "Seasonal hours ≥ 25 °C",
      caption = caption_text
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  if ("site_no" %in% names(pred_flow) && any(!is.na(pred_flow$site_no))) {
    p_flow = p_flow +
      ggplot2::geom_line(
        data = pred_flow,
        ggplot2::aes(x = median_flow, y = .fitted, group = site_no),
        inherit.aes = FALSE, color = "black", linewidth = 0.9
      )
  } else {
    p_flow = p_flow +
      ggplot2::geom_ribbon(
        data = pred_flow,
        ggplot2::aes(x = median_flow, ymin = lwr, ymax = upr),
        inherit.aes = FALSE, alpha = 0.12, fill = "grey70"
      ) +
      ggplot2::geom_line(
        data = pred_flow,
        ggplot2::aes(x = median_flow, y = .fitted),
        inherit.aes = FALSE, color = "black", linewidth = 1
      )
  }
  
  if (!is.null(save_prefix)) {
    ggplot2::ggsave(paste0(save_prefix, "_partial_temp.png"), p_temp,  width = width, height = height, dpi = dpi)
    ggplot2::ggsave(paste0(save_prefix, "_partial_flow.png"), p_flow,  width = width, height = height, dpi = dpi)
  }
  
  list(
    data           = air_stream_flow,
    model          = fit,
    coefs          = tb,
    pooled_glance  = gl,
    p_temp         = p_temp,
    p_flow         = p_flow,
    interpretation = as.character(interp)
  )
}

plot_z_anomalies_by_site <- function(
    df,
    title    = "Standardized Anomalies (z-scores) by Site",
    subtitle = "Within-site z-scores; streamflow sign reversed (positive = higher thermal risk)",
    save_path = NULL,
    width = 9, height = 6.5, dpi = 300
) {
  # ------------------------------------------------------------------------------
  # Title:    Standardized Exceedance, Air Temp, and Flow Anomalies by Site
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-24
  #
  # Description:
  #   Converts exceedance hours (≥21 °C), median May–Oct air temperature, and
  #   median May–Oct streamflow to within-site z-scores to allow direct comparison
  #   on a common scale. Streamflow z-scores are sign-reversed so that lower flows
  #   (higher thermal risk) appear in the same direction as higher exceedance
  #   hours and higher air temperatures.
  #
  # Inputs:
  #   df         : tibble with columns site_no, year, cumulative_hours_21,
  #                median_max_temp, median_flow
  #   title      : main plot title
  #   subtitle   : secondary text for plot context
  #   save_path  : optional file path for saving output (PNG)
  #   width/height/dpi : figure save dimensions
  #
  # Output:
  #   Facetted ggplot with standardized anomalies (z-scores) for each variable,
  #   shown across years and facetted by site.
  # ------------------------------------------------------------------------------
  stopifnot(all(c("site_no","year","cumulative_hours_21","median_max_temp","median_flow") %in% names(df)))
  
  df_z <- df %>%
    select(site_no, year, cumulative_hours_21, median_max_temp, median_flow) %>%
    filter(!is.na(cumulative_hours_21), !is.na(median_max_temp), !is.na(median_flow)) %>%
    group_by(site_no) %>%
    mutate(
      z_hours = as.numeric(scale(cumulative_hours_21)),
      z_temp  = as.numeric(scale(median_max_temp)),
      z_flow  = -as.numeric(scale(median_flow))   # reverse sign so lower flow => higher value
    ) %>%
    ungroup() %>%
    pivot_longer(
      cols = c(z_hours, z_temp, z_flow),
      names_to = "series", values_to = "z"
    ) %>%
    mutate(
      series = recode(series,
                      z_temp  = "Air temp (°C)",
                      z_hours = "Hours ≥21°C",
                      z_flow  = "Streamflow (reversed)"
      )
    )
  
  p <- ggplot(df_z, aes(x = year, y = z, group = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, color = "grey55") +
    
    # Hours ≥21 °C (full opacity)
    geom_line(
      data = dplyr::filter(df_z, series == "Hours ≥21°C"),
      aes(color = series),
      linewidth = 0.9, alpha = 1
    ) +
    geom_point(
      data = dplyr::filter(df_z, series == "Hours ≥21°C"),
      aes(color = series),
      size = 1.6, alpha = 1
    ) +
    
    # Air temp (lower alpha)
    geom_line(
      data = dplyr::filter(df_z, series == "Air temp (°C)"),
      aes(color = series),
      linewidth = 0.9, alpha = 0.5
    ) +
    geom_point(
      data = dplyr::filter(df_z, series == "Air temp (°C)"),
      aes(color = series),
      size = 1.6, alpha = 0.5
    ) +
    
    # Streamflow (lower alpha)
    geom_line(
      data = dplyr::filter(df_z, series == "Streamflow (reversed)"),
      aes(color = series),
      linewidth = 0.9, alpha = 0.6
    ) +
    geom_point(
      data = dplyr::filter(df_z, series == "Streamflow (reversed)"),
      aes(color = series),
      size = 1.6, alpha = 0.6
    ) +
    
    facet_wrap(~ site_no, scales = "free_x") +
    scale_color_manual(
      values = c(
        "Hours ≥21°C"           = "darkred",
        "Air temp (°C)"         = "forestgreen",
        "Streamflow (reversed)" = "steelblue"
      ),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Year",
      y = "Z-score (centered & scaled within site)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  p
}

plot_streamflow_trends_usgs = function(
    usgs_raw,
    discharge_col = "X_00060_00000",   # USGS discharge column (CFS)
    site_col      = "site_no",
    time_col      = "dateTime",
    months        = 5:10,               # May–Oct
    end_year      = 2024,
    palette_fun   = viridis::turbo,
    save_path     = NULL,
    width = 10, height = 7, dpi = 300
){
  # ------------------------------------------------------------------------------
  # Title:    Annual Median Streamflow Trends (May–Oct) from raw USGS subhourly data
  # Author:   Zachary H. Hoylman
  # Date:     2025-08-21
  #
  # Description:
  #   Uses raw USGS subhourly points to compute daily median discharge, then
  #   aggregates to May–Oct median discharge by site–year. Fits a pooled linear
  #   model with site fixed effects to estimate the basin-wide trend (CFS/decade),
  #   and per-site simple trends; draws site trend lines only where p < 0.05.
  #
  # Inputs:
  #   usgs_raw     : data.frame/tibble with at least:
  #                  - site_col (e.g., "site_no")
  #                  - time_col (e.g., "dateTime", POSIXt or parseable)
  #                  - discharge_col (e.g., "X_00060_00000", numeric CFS)
  #   discharge_col: name of discharge column in CFS
  #   site_col     : name of site id column
  #   time_col     : name of datetime column (POSIXt or string)
  #   months       : months to include (default May–Oct)
  #   end_year     : last year to include
  #   palette_fun  : function(n) -> n colors (default viridis::turbo)
  #   save_path    : optional file path to save figure
  #   width,height,dpi : graphics params for saving
  #
  # Output:
  #   ggplot object with attributes:
  #     - pooled_slope_CFS_per_decade
  #     - pooled_slope_CI_CFS_per_decade
  #     - site_slopes (per-site slope & p-value)
  #
  # Notes:
  #   * Daily median discharge = median of all subhourly points for a calendar day.
  #   * Seasonal median = median of daily medians within May–Oct of a year.
  #   * Pooled model: median_flow ~ year + site_no (site fixed effects).
  #   * Site lines drawn only if site p < 0.05 in simple lm(median_flow ~ year).
  # ------------------------------------------------------------------------------
  
  # --- basic checks & standardize columns
  stopifnot(all(c(site_col, time_col, discharge_col) %in% names(usgs_raw)))
  df <- usgs_raw
  names(df)[match(c(site_col, time_col, discharge_col), names(df))] <- c("site_no","dateTime","discharge")
  
  # coerce datetime
  if (!inherits(df$dateTime, "POSIXt")) {
    df$dateTime <- as.POSIXct(df$dateTime, tz = "UTC")
  }
  
  # keep finite discharge
  df <- df |>
    dplyr::filter(is.finite(discharge))
  
  # --- daily median discharge per site/day
  daily <- df |>
    dplyr::mutate(date = as.Date(lubridate::floor_date(dateTime, unit = "day"))) |>
    dplyr::group_by(site_no, date) |>
    dplyr::summarise(daily_median_discharge = stats::median(discharge, na.rm = TRUE),
                     .groups = "drop")
  
  # --- May–Oct seasonal median by site–year
  flow_trends <- daily |>
    dplyr::mutate(year = lubridate::year(date),
                  mon  = lubridate::month(date)) |>
    dplyr::filter(mon %in% months, year <= end_year) |>
    dplyr::group_by(site_no, year) |>
    dplyr::summarise(median_flow = stats::median(daily_median_discharge, na.rm = TRUE),
                     .groups = "drop")
  
  # order facets/colors by long-term average flow
  site_order <- flow_trends |>
    dplyr::group_by(site_no) |>
    dplyr::summarise(avg_med = mean(median_flow, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(avg_med)) |>
    dplyr::pull(site_no)
  flow_trends$site_no <- factor(flow_trends$site_no, levels = site_order)
  
  # --- pooled site-fixed-effects regression (CFS/decade)
  pooled_lm <- stats::lm(median_flow ~ year + site_no, data = flow_trends)
  slope_dec <- unname(stats::coef(pooled_lm)["year"] * 10)
  ci_dec    <- stats::confint(pooled_lm)["year", ] * 10
  
  slope_txt <- sprintf("%.1f", slope_dec)
  ci_lo_txt <- sprintf("%.1f", ci_dec[1])
  ci_hi_txt <- sprintf("%.1f", ci_dec[2])
  
  # --- per-site simple trends
  site_slopes <- flow_trends |>
    dplyr::group_by(site_no) |>
    dplyr::group_map(~{
      mod <- stats::lm(median_flow ~ year, data = .x)
      tt  <- broom::tidy(mod)
      out <- tt[tt$term == "year", , drop = FALSE]
      dplyr::tibble(
        slope_CFS_y      = out$estimate,
        slope_CFS_decade = out$estimate * 10,
        p_value          = out$p.value
      )
    }) |>
    dplyr::bind_rows(.id = "site_no") |>
    dplyr::mutate(
      site_no = levels(flow_trends$site_no)[as.integer(site_no)],
      sig     = p_value < 0.05
    )
  
  all_sig <- all(site_slopes$sig)
  
  caption_txt <- stringr::str_wrap(
    paste0(
      "Lines show annual May–Oct medians from raw USGS subhourly discharge. ",
      "Pooled trend (site fixed effects): ", slope_txt, " CFS/decade (95% CI ",
      ci_lo_txt, " to ", ci_hi_txt, "). ",
      if (all_sig) "All site-specific trends are significant (p < 0.05)."
      else "Site trend lines are drawn only where p < 0.05."
    ),
    width = 110
  )
  
  n_sites <- length(levels(flow_trends$site_no))
  cols    <- palette_fun(max(n_sites, 3)) |> rev()
  
  p <- ggplot2::ggplot(
    flow_trends,
    ggplot2::aes(x = year, y = median_flow, color = site_no, group = site_no)
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.35) +
    {
      if (all_sig) {
        ggplot2::geom_smooth(method = "lm", se = FALSE)
      } else {
        sig_sites <- site_slopes$site_no[site_slopes$sig]
        ggplot2::geom_smooth(
          data = dplyr::filter(flow_trends, site_no %in% sig_sites),
          method = "lm", se = FALSE
        )
      }
    } +
    ggplot2::scale_color_manual(values = cols, name = "Site") +
    ggplot2::labs(
      title    = "Annual Median Streamflow (May–Oct) for USGS Sites",
      subtitle = "Daily medians from subhourly USGS discharge; seasonal medians by site–year",
      x = "Year", y = "Discharge (CFS)", caption = caption_txt
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.caption  = ggplot2::element_text(hjust = 0),
      legend.title  = ggplot2::element_text(face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  attr(p, "pooled_slope_CFS_per_decade")    <- slope_dec
  attr(p, "pooled_slope_CI_CFS_per_decade") <- ci_dec
  attr(p, "site_slopes")                    <- site_slopes
  return(p)
}
