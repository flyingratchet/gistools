# CORE GIS FUNCTIONS ------------------------------------------------------

#' Function which filters coordinates in a data frame based on whether they
#' fall within a polygon
#' @param df a data frame of occurrence records that contains lat and decimal_longitude columns
#' @param poly a shape files with polygons
#' @param lat decimal latitude field from the data from
#' @param decimal_longitude decimal longitude field from the data from
#' @param results_col a character string representing the column where you like the label to be stored that indicates selected points
#' @param tag a character string with which to label selected points
#' @param crs a projection system
#' @return a data frame with records filtered based on falling inside the focal polygon
#' @export
tag_by_poly <- function(df, poly, decimal_latitude = "decimal_latitude", decimal_longitude = "decimal_longitude", crs, results_col = "results", tag){


    df_no_coords <- df %>% filter(is.na(decimal_latitude) | is.na(decimal_longitude))
    df %<>% filter(!is.na(decimal_latitude) | !is.na(decimal_longitude))

    # convert data frame to simple feature (sf) format
    df_points <- df %>% sf::st_as_sf(coords = c(x = {{decimal_longitude}}, y = {{decimal_latitude}}), crs = crs)

    # query which points are within polygon and record results in focal column
    # final lengths > 0 creates a TRUE/FALSE vector from the st_within output
    indicator <- sf::st_within(df_points, poly) %>% lengths > 0

    # note that the last argument in if_else needs to be the content of the column and not the column itself (confusingly)
    df[{{results_col}}] <- if_else(indicator, tag, df[[{{results_col}}]])

    df <- bind_rows(df, df_no_coords)

    return(df)
}





#' Build a unified, time-sorted tibble of GPX trackpoints
#'
#' Recursively finds all .gpx files under `path`, parses `<trkpt>` nodes,
#' and returns a tibble of latitude, longitude, optional elevation, and timestamps.
#' Namespace-safe (uses local-name()) and robust to empty/bad files.
#'
#' @param path Character. Parent folder containing GPX files.
#' @param glob File pattern. Default `"*.gpx"`.
#' @param lat_attr Attribute name for latitude on `<trkpt>`. Default `"lat"`.
#' @param lon_attr Attribute name for longitude on `<trkpt>`. Default `"lon"`.
#' @param ele_node Child node name (local-name) for elevation. Default `"ele"`.
#' @param time_node Child node name (local-name) for time. Default `"time"`.
#' @param out_lat,out_lon,out_ele,out_time Output column names.
#' @param tz Timezone for parsed times. Default `"UTC"`.
#' @param quiet Suppress non-fatal parse warnings. Default `TRUE`.
#' @param latlon_digits Integer digits to round lat/lon to. Set `NULL` to skip rounding. Default `6`.
#' @param as_character Logical: if `TRUE`, output lat/lon as fixed-width strings with trailing zeros
#'   (using `latlon_digits`). Default `FALSE`.
#' @return A tibble with columns named per `out_*`. Elevation may be `NA`.
#' @export
build_coordinate_table <- function(
    path,
    glob            = "*.gpx",
    lat_attr        = "lat",
    lon_attr        = "lon",
    ele_node        = "ele",
    time_node       = "time",
    out_lat         = "decimal_latitude",
    out_lon         = "decimal_longitude",
    out_ele         = "elevation_in_meters",
    out_time        = "date_time",
    tz              = "UTC",
    quiet           = TRUE,
    latlon_digits   = 6,
    as_character    = FALSE
) {
  requireNamespace("fs", quietly = TRUE)
  requireNamespace("xml2", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)
  requireNamespace("lubridate", quietly = TRUE)

  files <- fs::dir_ls(path, recurse = TRUE, glob = glob)

  empty_out <- tibble::tibble(
    !!out_lat  := vector("double"),
    !!out_lon  := vector("double"),
    !!out_ele  := vector("double"),
    !!out_time := as.POSIXct(character())
  )
  if (length(files) == 0L) return(empty_out)

  parse_one_file <- function(f) {
    doc <- tryCatch(xml2::read_xml(f), error = function(e) NULL)
    if (is.null(doc)) {
      return(tibble::tibble(.lat = numeric(), .lon = numeric(), .ele = numeric(), .time = character()))
    }

    pts <- xml2::xml_find_all(doc, ".//*[local-name()='trkpt']")
    if (length(pts) == 0L) {
      return(tibble::tibble(.lat = numeric(), .lon = numeric(), .ele = numeric(), .time = character()))
    }

    get_child_text <- function(nodeset, child_local) {
      xml2::xml_text(xml2::xml_find_first(nodeset, paste0("./*[local-name()='", child_local, "']")))
    }

    tibble::tibble(
      .lat  = suppressWarnings(as.numeric(xml2::xml_attr(pts, lat_attr))),
      .lon  = suppressWarnings(as.numeric(xml2::xml_attr(pts, lon_attr))),
      .ele  = suppressWarnings(as.numeric(get_child_text(pts, ele_node))),
      .time = get_child_text(pts, time_node)
    )
  }

  if (quiet) {
    old <- getOption("warn"); on.exit(options(warn = old), add = TRUE); options(warn = -1)
  }

  dat <- files %>%
    purrr::map(parse_one_file) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  needed <- c(".lat", ".lon", ".ele", ".time")
  if (nrow(dat) == 0L || !all(needed %in% names(dat))) return(empty_out)

  dat %>%
    dplyr::mutate(
      .time = lubridate::ymd_hms(.time, quiet = TRUE, tz = tz),
      .ele  = round(.ele),
      .lat  = if (!is.null(latlon_digits)) round(.lat, latlon_digits) else .lat,
      .lon  = if (!is.null(latlon_digits)) round(.lon, latlon_digits) else .lon
    ) %>%
    dplyr::mutate(
      .lat = if (as_character && !is.null(latlon_digits))
        sprintf(paste0("%.", latlon_digits, "f"), .lat) else .lat,
      .lon = if (as_character && !is.null(latlon_digits))
        sprintf(paste0("%.", latlon_digits, "f"), .lon) else .lon
    ) %>%
    dplyr::rename(
      !!out_lat  := .lat,
      !!out_lon  := .lon,
      !!out_ele  := .ele,
      !!out_time := .time
    ) %>%
    dplyr::arrange(!!rlang::sym(out_time))
}








#' Populate elevation using the geonames server based on
#' lat and decimal_longitude columns in a dataframe
#' I obtained my username (which is 'roverso') from here to validate the elevation lookup request
#' http://www.geonames.org/export/web-services.html
#' @param df a database in SYMBIOTA format
#' @return a dataframe with rounded decimal_latitude decimal_longitude values
#' @export
find_elevation <- function(df){
    df_coords <- df %>% select(decimal_latitude, decimal_longitude)
    df_coords$decimal_latitude %<>% as.numeric()
    df_coords$decimal_longitude %<>% as.numeric()
    df_coords %<>% rename(decimalLatitude = decimal_latitude, decimalLongitude = decimal_longitude)
    elevation_query <- rgbif::elevation(df_coords, username = 'roverso')
    df$elevation_in_meters <- elevation_query$elevation_geonames %>% as.character()
    return(df)
}




#' Attach a polygon attribute to a data frame by spatial join
#'
#' @param df A Symbiota-like data.frame with `decimal_latitude` and `decimal_longitude`
#' @param shapefile An sf polygon layer
#' @param field_name The single attribute in `shapefile` to bring over
#' @param new_name Name for the new/updated column on `df`
#' @return `df` with `new_name` column filled where points intersect polygons
#' @export
reverse_geocode <- function(df, shapefile, field_name, new_name) {

    # little helper for message sprintf above
    `%||%` <- function(a, b) if (is.null(a)) b else a

    # --- input checks ----------------------------------------------------------
    if (!inherits(shapefile, "sf")) {
        stop("`shapefile` must be an sf object (e.g., read with sf::st_read()).")
    }
    if (!field_name %in% names(shapefile)) {
        stop(sprintf("`field_name` '%s' not found in `shapefile`.", field_name))
    }
    if (is.na(sf::st_crs(shapefile))) {
        stop("`shapefile` has no CRS set. Assign the correct CRS first (e.g., sf::st_set_crs(shapefile, 4326)) before calling this function.")
    }

    # --- clean coordinate columns ---------------------------------------------
    # coerce to numeric safely; keep originals intact
    decimal_latitude <- suppressWarnings(as.numeric(df$decimal_latitude))
    decimal_longitude <- suppressWarnings(as.numeric(df$decimal_longitude))

    # rows with usable coords
    has_coords <- !is.na(decimal_latitude) & !is.na(decimal_longitude)
    df_no_coords <- df[!has_coords, , drop = FALSE]
    df_coords    <- df[ has_coords, , drop = FALSE]

    if (nrow(df_coords) == 0) {
        warning("No rows with valid numeric coordinates; returning `df` unchanged.")
        return(df)
    }

    # --- build sf points (WGS84) ----------------------------------------------
    df_points <- sf::st_as_sf(
        df_coords,
        coords = c(x = "decimal_longitude", y = "decimal_latitude"),
        crs = 4326
    )

    # --- make polygons valid & align CRS --------------------------------------
    shapefile <- sf::st_make_valid(shapefile)

    if (sf::st_crs(shapefile) != sf::st_crs(df_points)) {
        # Informative message rather than failing hard
        message(
            sprintf(
                "Transforming shapefile from EPSG:%s to EPSG:%s for the join.",
                sf::st_crs(shapefile)$epsg %||% "unknown",
                sf::st_crs(df_points)$epsg %||% "unknown"
            )
        )
        shapefile <- sf::st_transform(shapefile, sf::st_crs(df_points))
    }

    # --- spatial join (points <- polygons) ------------------------------------
    # keep just the requested attribute (plus geometry)
    shp_min <- shapefile[, field_name, drop = FALSE]

    joined <- sf::st_join(
        df_points,
        shp_min,
        left = TRUE,
        join = sf::st_intersects
    )

    # pull the attribute back out
    results <- joined[[field_name]]

    # --- write into `new_name` safely -----------------------------------------
    if (new_name %in% names(df_coords)) {
        # only fill where target is NA/blank
        target <- df_coords[[new_name]]
        is_blank <- function(x) is.na(x) | (is.character(x) & trimws(x) == "")
        df_coords[[new_name]] <- dplyr::coalesce(
            dplyr::if_else(is_blank(target), NA, target),
            results
        )
    } else {
        df_coords[[new_name]] <- results
    }

    # --- stitch back rows without coords (unchanged) ---------------------------
    out <- dplyr::bind_rows(df_coords, df_no_coords)

    # drop sf class if it slipped in (it shouldn't here)
    if (inherits(out, "sf")) out <- sf::st_drop_geometry(out)

    out
}



#' Rounds decimal_latitude decimal_longitude values based on sensible precision and manage formats
#' based on accuracy information in an coordinate uncertainty field
#' in meters
#' @param df a database in SYMBIOTA format
#' @return a dataframe with rounded decimal_latitude decimal_longitude values
#' @export
geography_coordinate_rounder <- function(df){
    # Ensure decimal_latitude/decimal_longitude and coordinate_uncertainty are numeric
    df <- df %>%
        mutate(
            decimal_latitude = as.numeric(decimal_latitude),
            decimal_longitude = as.numeric(decimal_longitude),
            coordinate_uncertainty = as.numeric(coordinate_uncertainty)
        )

    # Separate rows based on the presence of numeric uncertainty
    df_valid <- df %>%
        filter(!is.na(decimal_latitude) & grepl("\\d", coordinate_uncertainty))

    df_no_coords <- df %>%
        filter(is.na(decimal_latitude) | !grepl("\\d", coordinate_uncertainty))

    # Determine the number of decimal places based on coordinate_uncertainty
    df_valid <- df_valid %>%
        mutate(
            roundNum = case_when(
                coordinate_uncertainty >= 1000 ~ 2,
                coordinate_uncertainty >= 100 & coordinate_uncertainty < 1000 ~ 3,
                coordinate_uncertainty >= 10 & coordinate_uncertainty < 100 ~ 4,
                coordinate_uncertainty < 10 ~ 5,
                TRUE ~ NA_real_
            )
        )

    # Round decimal_latitude/decimal_longitude values and format them to preserve trailing zeros
    df_valid <- df_valid %>%
        rowwise() %>%
        mutate(
            decimal_latitude = formatC(round(decimal_latitude, roundNum), format = "f", digits = roundNum),
            decimal_longitude = formatC(round(decimal_longitude, roundNum), format = "f", digits = roundNum)
        ) %>%
        ungroup()

    # Convert decimal_latitude/decimal_longitude to character to ensure consistency
    df_valid <- df_valid %>%
        mutate(
            decimal_latitude = as.character(decimal_latitude),
            decimal_longitude = as.character(decimal_longitude)
        )

    # Ensure df_no_coords decimal_latitude/decimal_longitude are also characters
    df_no_coords <- df_no_coords %>%
        mutate(
            decimal_latitude = as.character(decimal_latitude),
            decimal_longitude = as.character(decimal_longitude)
        )

    # Combine all parts of the data frame back together
    df <- bind_rows(df_valid, df_no_coords)

    # Sort by a relevant column if it exists
    if("SpecimenCode" %in% colnames(df)){
        df <- df %>% arrange(SpecimenCode)
    } else if("CollectionCode" %in% colnames(df)){
        df <- df %>% arrange(CollectionCode)
    }

    return(df)
}


#' DEPRACATED
#' Rounds decimal_latitude decimal_longitude values based on sensible precision and manage formats
#' based on accuracy information in an coordinate uncertainty field
#' in meters
#' @param df a database in SYMBIOTA format
#' @return a dataframe with rounded decimal_latitude decimal_longitude values
#' @export
round_coords_cc <- function(df){
    # Make sure decimal_latitude/decimal_longitude associated columns are numeric
    df$decimal_latitude <- as.numeric(df$decimal_latitude)
    df$decimal_longitude <- as.numeric(df$decimal_longitude)
    df$coordinate_uncertainty <- as.numeric(df$coordinate_uncertainty)

    # Temporarily remove records with "NA" or no numbers in the coordinate_uncertainty field so they don't get processed
    dfnoNumUncertainty <- df[grepl("^([^0-9]*)$", df$coordinate_uncertainty, perl = TRUE),]
    dfNAUncertainty <- df[is.na(df$coordinate_uncertainty),]
    df <- df[grepl("\\d", df$coordinate_uncertainty, perl = TRUE),]
    # remove records without coords temporarily
    dfNoCoords <- df %>% filter(is.na(decimal_latitude))
    dfCoords <- df %>% filter(!is.na(decimal_latitude))

    # ROUND decimal_latitude/decimal_longitude TO APPROPRIATE VALUES BASED ON ACCURACY
    dfCoords['roundNum'] <- NA
    dfCoords$roundNum[dfCoords$coordinate_uncertainty >= 1000] <- 2
    dfCoords$roundNum[dfCoords$coordinate_uncertainty >= 100 & dfCoords$coordinate_uncertainty < 1000] <- 3
    dfCoords$roundNum[dfCoords$coordinate_uncertainty >= 10 & dfCoords$coordinate_uncertainty < 100] <- 4
    dfCoords$roundNum[dfCoords$coordinate_uncertainty < 10] <- 5
    dfCoords$roundNum <- as.numeric(dfCoords$roundNum)

    # Overwrite rounded decimal_latitude/decimal_longitude results with original values
    results <- NA
    for(i in 1:nrow(dfCoords)){
        results[[i]] <- formatC(dfCoords$decimal_latitude[i], dfCoords$roundNum[i], format = "f")
    }
    dfCoords$decimal_latitude <- results

    results <- NA
    for(i in 1:nrow(dfCoords)){
        results[[i]] <- formatC(dfCoords$decimal_longitude[i], dfCoords$roundNum[i], format = "f")
    }
    dfCoords$decimal_longitude <- results

    # Delete columns that aren't needed
    dfCoords$roundNum <- NULL
    dfCoords$decimal_latitude <- as.numeric(dfCoords$decimal_latitude)
    dfCoords$decimal_longitude <- as.numeric(dfCoords$decimal_longitude)

    # Add back records that were removed pre-processing above and then re-sort database
    df <- bind_rows(dfCoords, dfNoCoords, dfNAUncertainty, dfnoNumUncertainty)

    # Format as character values so no funny business happens to decimal_latitude/decimal_longitude
    df$decimal_latitude <- as.character(df$decimal_latitude)
    df$decimal_longitude <- as.character(df$decimal_longitude)

    # Sort whatever column is necessary before returning output since I cut and pasted the df so much
    if("SpecimenCode" %in% colnames(df)){
        df <- df[with(df, order(SpecimenCode)),] # sort cc.new;
    } else{ if("CollectionCode" %in% colnames(df)){
        df <- df[with(df, order(CollectionCode)),] # sort cc.new
    }
    }
    return(df)
}






#' Geoencode records with GPX using per-record UTC offsets
#'
#' @param records Data frame with event_date, event_time, decimal_latitude, decimal_longitude, and timezone columns
#' @param time_coords GPS coordinates with UTC timestamps (data.frame or CSV path).
#'        Must have columns: date_time, decimal_latitude, decimal_longitude. Elevation column optional.
#' @param time_zone Column name in `records` containing UTC offsets per record (e.g., "UTC-7", "UTC+5.5")
#' @param max_gap_mins Maximum gap in minutes between record and nearest GPX point (default: 15)
#' @param elevation_col Name of the elevation column to write in `records` (default: "elevation_in_meters")
#' @param gps_elevation_col Name of the elevation column in the GPX data (default: "elevation_in_meters")
#' @param overwrite_latlon_only_if_na If TRUE, only fill decimal_latitude/lon when NA; else overwrite on match (default: FALSE)
#' @param overwrite_elev_only_if_na   If TRUE, only fill elevation when NA; else overwrite on match (default: FALSE)
#' @param verbose Print progress messages (default: TRUE)
#' @return `records` with updated decimal_latitude, decimal_longitude, and (optionally) `elevation_col`.
#'         Adds a helper column `time_diff_mins` indicating match gap (mins).
#' @export
#'
#' @importFrom data.table as.data.table setkey
#' @importFrom data.table .I .N
#' @importFrom data.table `:=`
extract_coordinates <- function(
    records,
    time_coords,
    time_zone,
    max_gap_mins = 15,
    elevation_col = "elevation_in_meters",
    gps_elevation_col = "elevation_in_meters",
    overwrite_latlon_only_if_na = FALSE,
    overwrite_elev_only_if_na   = FALSE,
    verbose = TRUE
) {
  require(dplyr, quietly = TRUE)
  require(data.table, quietly = TRUE)
  require(lubridate, quietly = TRUE)
  require(purrr, quietly = TRUE)

  # --- validate & snapshot types so we can restore later ---------------------
  stopifnot(is.data.frame(records))
  req_cols <- c("event_date", "event_time", "decimal_latitude", "decimal_longitude")
  miss <- setdiff(req_cols, names(records))
  if (length(miss)) stop("Missing required columns in `records`: ", paste(miss, collapse = ", "))

  if (!time_zone %in% names(records)) {
    stop(sprintf("Column '%s' not found in records.", time_zone))
  }

  col_classes <- purrr::map_chr(records, ~ class(.x)[1])
  orig_cols   <- names(records)

  # Coerce decimal_latitude/decimal_longitude to numeric safely (if present as character)
  records$decimal_latitude  <- as.numeric_safe(records$decimal_latitude)
  records$decimal_longitude <- as.numeric_safe(records$decimal_longitude)

  # Ensure elevation target column exists (numeric by default)
  if (!elevation_col %in% names(records)) records[[elevation_col]] <- NA_real_

  # --- parse local datetime + timezone to UTC --------------------------------
  records$local_datetime <- create_local_datetime(records$event_date, records$event_time)

  # Separate those with usable times
  with_time <- dplyr::filter(records, !is.na(local_datetime) & !is.na(.data[[time_zone]]))
  no_time   <- dplyr::filter(records,  is.na(local_datetime) |  is.na(.data[[time_zone]]))

  # Parse per-record time zone offsets (e.g., "UTC-7", "UTC+5.5")
  tz_info <- check_time_zone_strings(with_time[[time_zone]])
  if (any(!tz_info$ok)) {
    bad <- unique(as.character(with_time[[time_zone]][!tz_info$ok]))
    stop("Invalid timezone format found: ", paste(bad, collapse = ", "),
         "\nUse 'UTC±N' forms (e.g., 'UTC-7', 'UTC+5', 'UTC+5.5').")
  }
  with_time$offset_hours <- tz_info$offset_hours
  with_time$utc_time     <- with_time$local_datetime - lubridate::dhours(with_time$offset_hours)
  attr(with_time$utc_time, "tzone") <- "UTC"

  # --- load & normalize GPX ---------------------------------------------------
  gps_dt <- load_gps_data(time_coords, gps_elevation_col = gps_elevation_col) # adds gpx_elev

  # --- nearest-time match (rolling join) -------------------------------------
  rec_dt <- as.data.table(with_time)[, rec_idx := .I]
  matched <- gps_dt[rec_dt, on = .(date_time = utc_time), roll = "nearest",
                    .(
                      rec_idx        = i.rec_idx,
                      time_diff_mins = abs(as.numeric(difftime(x.date_time, i.utc_time, units = "mins"))),
                      gpx_lat        = x.decimal_latitude,
                      gpx_lon        = x.decimal_longitude,
                      gpx_elev       = x.gpx_elev
                    )
  ]

  within_gap <- matched$time_diff_mins <= max_gap_mins & !is.na(matched$time_diff_mins)

  # --- apply updates ----------------------------------------------------------
  updated_with_time <- update_from_matches(
    with_time, matched, within_gap, elevation_col,
    overwrite_latlon_only_if_na = overwrite_latlon_only_if_na,
    overwrite_elev_only_if_na   = overwrite_elev_only_if_na
  )

  # Recombine with no-time rows
  out <- dplyr::bind_rows(updated_with_time, no_time) %>%
    dplyr::arrange(event_date, event_time)

  # Restore original column classes for original columns only
  out <- recoerce_columns(out, col_classes)

  if (verbose) {
    cat("\n", sum(within_gap, na.rm = TRUE), " records updated (≤", max_gap_mins, " mins).\n", sep = "")
    n_unmatched <- nrow(with_time) - sum(within_gap, na.rm = TRUE)
    if (n_unmatched > 0) cat(n_unmatched, " records were unmatched and not updated\n", sep = "")
    if (nrow(no_time) > 0) cat("SKIPPED (no usable time): ", nrow(no_time), "\n", sep = "")
  }

  # Keep original cols plus updated fields + time_diff_mins if present
  keep <- unique(c(orig_cols, "decimal_latitude", "decimal_longitude", elevation_col, "time_diff_mins"))
  keep <- intersect(keep, names(out))
  dplyr::select(out, dplyr::all_of(keep))
}


#' @noRd
#' @keywords internal
as.numeric_safe <- function(x) {
  suppressWarnings(as.numeric(ifelse(is.na(x) | x == "", NA, x)))
}

#' @noRd
#' @keywords internal
create_local_datetime <- function(date_col, time_col) {
  datetime_str <- paste(date_col, time_col)
  # Try with seconds (YYYY-mm-dd HH:MM:SS)
  local_dt <- as.POSIXct(datetime_str, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  # Fallback to minutes (YYYY-mm-dd HH:MM)
  need_min <- is.na(local_dt) & !is.na(time_col)
  if (any(need_min)) {
    local_dt[need_min] <- as.POSIXct(datetime_str[need_min], format = "%Y-%m-%d %H:%M", tz = "UTC")
  }
  local_dt
}

#' @noRd
#' @keywords internal
check_time_zone_strings <- function(x) {
  x_char <- as.character(x)

  # Match exactly "UTC" followed by +/- and digits (with optional decimals)
  pattern <- "^UTC([+-])(\\d+(?:\\.\\d+)?)$"
  matches <- regexec(pattern, x_char, ignore.case = TRUE)

  offset_hours <- rep(NA_real_, length(x_char))

  for (i in seq_along(matches)) {
    match_data <- regmatches(x_char[i], matches[i])[[1]]
    if (length(match_data) == 3) {  # Full match + 2 capture groups
      sign <- ifelse(match_data[2] == "-", -1, 1)
      hours <- as.numeric(match_data[3])
      offset_hours[i] <- sign * hours
    }
  }

  data.frame(
    input = x_char,
    offset_hours = offset_hours,
    ok = !is.na(offset_hours),
    stringsAsFactors = FALSE
  )
}

#' @noRd
#' @keywords internal
load_gps_data <- function(time_coords, gps_elevation_col = "elevation_in_meters") {
  if (is.character(time_coords) && file.exists(time_coords)) {
    time_coords <- read.csv(time_coords, stringsAsFactors = FALSE)
  }

  required_cols <- c("date_time", "decimal_latitude", "decimal_longitude")
  missing_cols <- setdiff(required_cols, names(time_coords))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in time_coords: ", paste(missing_cols, collapse = ", "))
  }

  gps_dt <- as.data.table(time_coords)
  gps_dt[, `:=`(
    decimal_latitude  = as.numeric(decimal_latitude),
    decimal_longitude = as.numeric(decimal_longitude),
    date_time         = parse_utc_datetime(date_time)
  )]

  # Normalize elevation column to gpx_elev for reliable use in data.table joins
  if (gps_elevation_col %in% names(gps_dt)) {
    gps_dt[, gpx_elev := suppressWarnings(as.numeric(get(gps_elevation_col)))]
  } else {
    gps_dt[, gpx_elev := NA_real_]
  }

  setkey(gps_dt, date_time)
  gps_dt
}

#' @noRd
#' @keywords internal
parse_utc_datetime <- function(datetime_col) {
  if (inherits(datetime_col, "POSIXct")) return(lubridate::with_tz(datetime_col, "UTC"))
  parsed <- lubridate::ymd_hms(as.character(datetime_col), tz = "UTC", quiet = TRUE)
  if (any(is.na(parsed))) {
    parsed2 <- as.POSIXct(as.character(datetime_col), format = "%Y-%m-%d %H:%M", tz = "UTC")
    parsed[is.na(parsed)] <- parsed2[is.na(parsed)]
  }
  parsed
}

#' @noRd
#' @keywords internal
update_from_matches <- function(records_with_time, matched, within_gap, elevation_col,
                                overwrite_latlon_only_if_na = FALSE,
                                overwrite_elev_only_if_na   = FALSE) {
  if (!"time_diff_mins" %in% names(records_with_time)) records_with_time$time_diff_mins <- NA_real_

  idx <- matched$rec_idx[within_gap]
  if (length(idx) == 0) return(records_with_time)

  # Map the matched vectors to the TRUEs in within_gap
  lat_new  <- matched$gpx_lat[within_gap]
  lon_new  <- matched$gpx_lon[within_gap]
  elev_new <- matched$gpx_elev[within_gap]
  gap_new  <- matched$time_diff_mins[within_gap]

  # LAT
  if (overwrite_latlon_only_if_na) {
    write_lat <- is.na(records_with_time$decimal_latitude[idx])
  } else {
    write_lat <- rep(TRUE, length(idx))
  }
  records_with_time$decimal_latitude[idx[write_lat]] <- lat_new[write_lat]

  # LON
  if (overwrite_latlon_only_if_na) {
    write_lon <- is.na(records_with_time$decimal_longitude[idx])
  } else {
    write_lon <- rep(TRUE, length(idx))
  }
  records_with_time$decimal_longitude[idx[write_lon]] <- lon_new[write_lon]

  # ELEVATION
  if (overwrite_elev_only_if_na) {
    write_elev <- is.na(records_with_time[[elevation_col]][idx])
  } else {
    write_elev <- rep(TRUE, length(idx))
  }
  # If elevation column is character, keep type; else numeric
  if (is.character(records_with_time[[elevation_col]])) {
    records_with_time[[elevation_col]][idx[write_elev]] <- as.character(elev_new[write_elev])
  } else {
    records_with_time[[elevation_col]][idx[write_elev]] <- suppressWarnings(as.numeric(elev_new[write_elev]))
  }

  # Bookkeeping
  records_with_time$time_diff_mins[idx] <- gap_new
  records_with_time
}
#' @noRd
#' @keywords internal
recoerce_columns <- function(df, ref_classes) {
  for (col in intersect(names(df), names(ref_classes))) {
    df[[col]] <- switch(
      ref_classes[[col]],
      character = as.character(df[[col]]),
      numeric   = as.numeric(df[[col]]),
      double    = as.double(df[[col]]),
      integer   = as.integer(df[[col]]),
      factor    = as.factor(df[[col]]),
      df[[col]]  # fallback/no-op
    )
  }
  df
}





































# NEW REFACTORED GPX AND KML PIPELINE -------------------------------------

#' Convert GPX Files to KML with Labeled Points and Trails (Auto-Timezone Detection)
#'
#' This function processes all `.gpx` files in a given trip folder and converts
#' them into `.kml` files. Each output KML file includes:
#' - A styled `LineString` showing the GPS trail
#' - Individual timestamped `Point` placemarks with optional local time labels
#'
#' The function automatically detects the local timezone based on GPS coordinates
#' using the `lutz` package, eliminating the need for manual timezone configuration.
#' For multi-timezone trips, it uses the timezone of the first trackpoint.
#'
#' The function also generates a combined `all_tracks.kml` file for the trip folder
#' when multiple GPX files are present.
#'
#' @param trip_folder A character path to a folder containing one or more `.gpx` files.
#' @param hover_labels If TRUE, labels only appear on hover/click (default: FALSE - always visible)
#' @param point_frequency Show every Nth trackpoint (default: 1 - all points)
#' @param line_color KML color code for track lines (default: "ff0000ff" - red)
#' @param line_width Width of track lines in pixels (default: 3)
#' @param point_scale Scale factor for point icons (default: 0.6)
#' @param label_scale Scale factor for labels (default: 0.6)
#'
#' @return Writes `.kml` files and optionally a combined `all_tracks.kml` into a `kml/` subfolder.
#' @export

convert_dir_to_kml <- function(trip_folder,
                               hover_labels = FALSE,
                               point_frequency = 1,
                               line_color = "ff0000ff",
                               line_width = 3,
                               point_scale = 0.6,
                               label_scale = 0.6) {

    message("🚀 Converting: ", trip_folder)
    output_dir <- file.path(trip_folder, "kml")

    # Main processing
    gpx_files <- find_gpx_files(trip_folder)
    if (length(gpx_files) == 0) {
        message("⚠️ No GPX files found in: ", trip_folder)
        return(invisible(NULL))
    }

    # Process each GPX file
    for (gpx_file in gpx_files) {
        process_single_gpx(
            gpx_file = gpx_file,
            output_dir = output_dir,
            hover_labels = hover_labels,
            point_frequency = point_frequency,
            line_color = line_color,
            line_width = line_width,
            point_scale = point_scale,
            label_scale = label_scale
        )
    }

    # Create combined file if multiple GPX files
    if (length(gpx_files) > 1) {
        create_combined_kml(output_dir)
    } else {
        message("📝 Only one GPX file found - skipping combined file creation")
    }
}

#' Null-coalescing operator helper
#' @param x First value
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Find GPX files in a directory
#' @param folder Path to search for GPX files
#' @return Vector of GPX file paths
find_gpx_files <- function(folder) {
    gpx_files <- dir_ls(folder, regexp = "\\.gpx$", recurse = TRUE, type = "file")
    gpx_files[!grepl("/kml/", gpx_files)]  # Exclude files in kml subdirectories
}

#' Extract and process trackpoint data from GPX file
#' @param gpx_file Path to GPX file
#' @return List with decimal_latitude, lon, time_utc, time_local, local_tz
extract_trackpoint_data <- function(gpx_file) {
    xml <- read_xml(gpx_file)
    xml_clean <- xml_ns_strip(xml)
    trkpts <- xml_find_all(xml_clean, ".//trkpt | .//wpt | .//rtept")

    if (length(trkpts) == 0) {
        message("⚠️ No trackpoints found in: ", gpx_file)
        return(NULL)
    }

    # Extract coordinates and time
    decimal_latitude <- xml_attr(trkpts, "decimal_latitude")
    decimal_longitude <- xml_attr(trkpts, "decimal_longitude")
    time_raw <- map_chr(trkpts, ~ xml_text(xml_find_first(.x, "time")))
    time_parsed <- suppressWarnings(as.POSIXct(time_raw, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))

    # Filter valid data
    valid_idx <- which(!is.na(decimal_latitude) & !is.na(decimal_longitude) & !is.na(time_parsed))
    if (length(valid_idx) == 0) {
        message("⚠️ No valid data in: ", gpx_file)
        return(NULL)
    }

    decimal_latitude <- decimal_latitude[valid_idx]
    decimal_longitude <- decimal_longitude[valid_idx]
    time_utc <- time_parsed[valid_idx]

    # Detect and convert to local timezone
    local_tz <- detect_timezone(decimal_latitude, decimal_longitude)
    time_local <- as.POSIXct(time_utc, tz = "UTC")
    attr(time_local, "tzone") <- local_tz

    list(
        decimal_latitude = decimal_latitude,
        decimal_longitude = decimal_longitude,
        time_utc = time_utc,
        time_local = time_local,
        local_tz = local_tz
    )
}

#' Detect timezone from GPS coordinates
#' @param decimal_latitude Vector of latitudes
#' @param decimal_longitude Vector of longitudes
#' @return Timezone string
detect_timezone <- function(lat, lon) {
    tryCatch({
        if (length(lat) > 0 && length(lon) > 0 && !is.na(lat[1]) && !is.na(lon[1])) {
            tz_name <- tz_lookup_coords(lat = as.numeric(lat[1]),
                                        lon = as.numeric(lon[1]),
                                        method = "accurate")
            message("📍 Detected timezone: ", tz_name)
            return(tz_name)
        }
    }, error = function(e) {
        message("⚠️ Timezone detection failed: ", e$message)
    })

    message("⚠️ Using UTC as fallback")
    return("UTC")
}

#' Calculate bearing between two points
#' @param lat1,lon1 First point coordinates
#' @param lat2,lon2 Second point coordinates
#' @return Bearing in degrees
calculate_heading <- function(lat1, lon1, lat2, lon2) {
    rad <- pi / 180
    dLon <- (lon2 - lon1) * rad
    lat1 <- lat1 * rad
    lat2 <- lat2 * rad

    y <- sin(dLon) * cos(lat2)
    x <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon)
    heading <- atan2(y, x) * (180 / pi)
    (heading + 360) %% 360
}

#' Generate KML style definitions
#' @param hover_labels Whether labels should only appear on hover
#' @param line_color KML color code for lines
#' @param line_width Width of lines
#' @param point_scale Scale of point icons
#' @param label_scale Scale of labels
#' @return KML style string
build_kml_styles <- function(hover_labels = FALSE,
                             line_color = "ff0000ff",
                             line_width = 3,
                             point_scale = 0.6,
                             label_scale = 0.6) {

    line_style <- glue('
    <Style id="lineStyle">
      <LineStyle>
        <color>{line_color}</color>
        <width>{line_width}</width>
      </LineStyle>
    </Style>')

    if (hover_labels) {
        # StyleMap for hover-only labels
        point_styles <- glue('
    <StyleMap id="arrowStyle">
      <Pair>
        <key>normal</key>
        <styleUrl>#arrowNormal</styleUrl>
      </Pair>
      <Pair>
        <key>highlight</key>
        <styleUrl>#arrowHighlight</styleUrl>
      </Pair>
    </StyleMap>

    <Style id="arrowNormal">
      <IconStyle>
        <scale>{point_scale}</scale>
        <Icon>
          <href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href>
        </Icon>
      </IconStyle>
      <LabelStyle>
        <scale>0</scale>
      </LabelStyle>
    </Style>

    <Style id="arrowHighlight">
      <IconStyle>
        <scale>{point_scale * 1.2}</scale>
        <Icon>
          <href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href>
        </Icon>
      </IconStyle>
      <LabelStyle>
        <scale>{label_scale}</scale>
      </LabelStyle>
    </Style>')
    } else {
        # Always-visible labels
        point_styles <- glue('
    <Style id="arrowStyle">
      <IconStyle>
        <scale>{point_scale}</scale>
        <Icon>
          <href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href>
        </Icon>
      </IconStyle>
      <LabelStyle>
        <scale>{label_scale}</scale>
      </LabelStyle>
    </Style>')
    }

    paste(line_style, point_styles, sep = "\n")
}

#' Create formatted time label with timezone
#' @param time_local Local time (POSIXct)
#' @param time_utc UTC time (POSIXct)
#' @param local_tz Timezone string
#' @return Formatted label string
format_time_label <- function(time_local, time_utc, local_tz) {
    local_time_formatted <- format(time_local, "%Y-%m-%d %H:%M")
    tz_abbrev <- format(time_local, "%Z")

    # If timezone abbreviation is same as timezone name, use offset instead
    if (tz_abbrev == local_tz || tz_abbrev == "" || is.na(tz_abbrev)) {
        offset_hours <- as.numeric(difftime(time_local, time_utc, units = "hours"))
        tz_display <- if (offset_hours >= 0) {
            paste0("UTC+", offset_hours)
        } else {
            paste0("UTC", offset_hours)
        }
    } else {
        tz_display <- tz_abbrev
    }

    paste0(local_time_formatted, " (", tz_display, ")")
}

#' Generate line placemark KML
#' @param data Trackpoint data list
#' @param filename GPX filename for labeling
#' @return KML string for line placemark
create_line_placemark <- function(data, filename) {
    coordinates_str <- paste(glue("{data$lon},{data$lat},0"), collapse = "\n")

    glue('
      <Placemark>
        <name>{basename(filename)}</name>
        <styleUrl>#lineStyle</styleUrl>
        <LineString>
          <tessellate>1</tessellate>
          <coordinates>
            {coordinates_str}
          </coordinates>
        </LineString>
      </Placemark>')
}

#' Generate point placemarks KML
#' @param data Trackpoint data list
#' @param point_frequency Show every Nth point
#' @return KML string for all point placemarks
create_point_placemarks <- function(data, point_frequency = 1) {
    # Apply frequency filter
    indices <- seq(1, length(data$lat), by = point_frequency)

    map_chr(indices, function(i) {
        label <- format_time_label(data$time_local[i], data$time_utc[i], data$local_tz)
        time_utc_formatted <- format(data$time_utc[i], "%Y-%m-%dT%H:%M:%SZ")

        # Calculate heading for arrow direction
        heading_tag <- ""
        if (i < length(data$lat)) {
            next_idx <- min(i + 1, length(data$lat))  # Handle edge case
            h <- calculate_heading(as.numeric(data$lat[i]), as.numeric(data$lon[i]),
                                   as.numeric(data$lat[next_idx]), as.numeric(data$lon[next_idx]))
            heading_tag <- glue("<heading>{round(h, 2)}</heading>")
        }

        glue('
      <Placemark>
        <name>{label}</name>
        <styleUrl>#arrowStyle</styleUrl>
        <TimeStamp><when>{time_utc_formatted}</when></TimeStamp>
        <description>{label}</description>
        <Point><coordinates>{data$lon[i]},{data$lat[i]},0</coordinates></Point>
        <Style>
          <IconStyle>
            {heading_tag}
          </IconStyle>
        </Style>
      </Placemark>')
    }) %>% paste(collapse = "\n")
}

#' Process a single GPX file into KML
#' @param gpx_file Path to GPX file
#' @param output_dir Output directory
#' @param ... Style parameters
process_single_gpx <- function(gpx_file, output_dir, ...) {
    # Collect all forwarded args once
    dots <- list(...)

    # Extract only the style-related args for build_kml_styles()
    styles <- build_kml_styles(
        hover_labels = dots$hover_labels %||% FALSE,
        line_color   = dots$line_color   %||% "ff0000ff",
        line_width   = dots$line_width   %||% 3,
        point_scale  = dots$point_scale  %||% 0.6,
        label_scale  = dots$label_scale  %||% 0.6
    )

    # Handle point_frequency **only** here (don’t pass it to styles)
    pfreq <- dots$point_frequency %||% 1

    # Extract trackpoint data
    data <- extract_trackpoint_data(gpx_file)
    if (is.null(data)) return(invisible(NULL))

    # Ensure output directory exists
    if (!dir_exists(output_dir)) dir_create(output_dir, recurse = TRUE)

    # Build placemarks
    line_placemark   <- create_line_placemark(data, gpx_file)
    point_placemarks <- create_point_placemarks(data, pfreq)

    # Assemble KML
    kml_doc <- glue('<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>{basename(gpx_file)}</name>
    {styles}
    {line_placemark}
    {point_placemarks}
  </Document>
</kml>')

    # Write file
    out_file <- file.path(output_dir, paste0(path_ext_remove(path_file(gpx_file)), ".kml"))
    writeLines(kml_doc, con = out_file)
    message("✅ Wrote: ", out_file)
}


#' Create combined KML file from individual files
#' @param output_dir Directory containing individual KML files
create_combined_kml <- function(output_dir) {
    kml_files <- dir_ls(output_dir, glob = "*.kml")

    placemark_bodies <- map_chr(kml_files, function(file) {
        content <- readLines(file, warn = FALSE)
        start <- which(grepl("<Document>", content))
        end   <- which(grepl("</Document>", content))
        if (length(start) == 1 && length(end) == 1 && end > start) {
            paste(content[(start + 1):(end - 1)], collapse = "\n")
        } else {
            warning("⚠️ Could not extract KML from: ", basename(file))
            return("")
        }
    })

    master_kml <- glue('<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>All Tracks</name>
    {paste(placemark_bodies, collapse = "\n")}
  </Document>
</kml>')

    output_file <- file.path(output_dir, "all_tracks.kml")
    writeLines(master_kml, output_file)
    message("📦 Master KML created: ", output_file)
}















