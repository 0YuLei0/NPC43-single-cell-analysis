# =============================================================================
# 4. add image  — AICL FFPE LN H&E (Seurat image "roi")
# =============================================================================
# Source this after `seu` is in the workspace (named list of Seurat objects).
# Attaches histology as image "roi" and leaves the original SlideSeq "image"
# unchanged. Samples with no matching file are skipped.
#
# Default image folder (Bouchet):
#   /home/ly385/project_pi_rf273/shared/ly385/AICL/03.HE
# Top-level HE.tif seen there: LN14734, LN18427, LN21720, LN22630, LN32298,
# LN3470, LN4737. LN3524 and LN5034 are not at the top level — the finder
# also searches v2/ and v3/ (and any other subfolder) for those.
#
# These TIFFs are ~1–3 GB. Install magick and run on a compute node, not
# login1. Images are downsampled to HE_MAX_PX on the long edge before they
# are stored on the Seurat object.
#
# Required: Seurat, ggplot2, patchwork
# Image I/O: magick (preferred) or tiff / jpeg / png
# =============================================================================

AICL_SAMPLE_IDS <- c(
  "LN14734", "LN18427", "LN21720", "LN22630", "LN32298",
  "LN3470", "LN3524", "LN4737", "LN5034"
)

if (!exists("HE_MAX_PX")) HE_MAX_PX <- 2000L

default_he_path <- function() {
  cand <- c(
    if (exists("PATH_HE_IMG")) PATH_HE_IMG else NULL,
    "/home/ly385/project_pi_rf273/shared/ly385/AICL/03.HE",
    file.path(getwd(), "03.HE"),
    if (identical(basename(getwd()), "03.HE")) getwd() else NULL
  )
  cand <- unique(as.character(cand[!vapply(cand, is.null, logical(1))]))
  hit <- cand[dir.exists(cand)]
  if (length(hit)) hit[[1]] else cand[[1]]
}

he_sid_token <- function(sid) {
  paste0("(^|[^0-9A-Za-z])", sid, "([^0-9]|$)")
}

list_he_image_files <- function(root) {
  if (!dir.exists(root)) return(character(0))
  list.files(
    root,
    pattern = "\\.(tif|tiff|jpg|jpeg|png)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

list_he_candidates <- function(sid, root = default_he_path()) {
  files <- list_he_image_files(root)
  if (!length(files)) return(character(0))
  keep <- grepl(he_sid_token(sid), basename(files), ignore.case = TRUE)
  files[keep]
}

he_candidate_rank <- function(path, sid) {
  bn <- basename(path)
  if (grepl(paste0("^FFPE_", sid, "\\.HE\\."), bn, ignore.case = TRUE)) return(1L)
  if (grepl(paste0("^FFPE_ALCL_", sid, "\\.HE\\."), bn, ignore.case = TRUE)) return(2L)
  if (grepl(paste0("^FFPE_", sid, "\\."), bn, ignore.case = TRUE)) return(3L)
  if (grepl("HE", bn, ignore.case = TRUE)) return(4L)
  if (grepl("ROI", bn, ignore.case = TRUE)) return(5L)
  6L
}

# Pick the best H&E / ROI file for one Seurat sample id (e.g. "LN14734").
# Preference: FFPE_{sid}.HE.* at a shallower path, then ALCL-prefixed HE,
# then any *HE* / *ROI* match, then any image whose filename contains sid.
# LN3470 will not match LN34700.
find_he_file <- function(sid, root = default_he_path()) {
  hits <- list_he_candidates(sid, root)
  if (!length(hits)) return(NA_character_)
  depth <- vapply(
    hits,
    function(p) length(strsplit(p, "/", fixed = TRUE)[[1]]),
    integer(1)
  )
  rank <- vapply(hits, he_candidate_rank, integer(1), sid = sid)
  hits[[order(rank, depth, basename(hits))[1]]]
}

find_he_files <- function(sids, root = default_he_path()) {
  setNames(vapply(sids, find_he_file, character(1), root = root), sids)
}

# Edit per sample: clockwise rotate in {0, 90, 180, 270}; flips apply to the
# image only. Start at identity and change after comparing images "image" vs "roi".
default_roi_orient <- function(sids) {
  setNames(lapply(sids, function(sid) {
    list(rotate = 0L, flip_x = FALSE, flip_y = FALSE)
  }), sids)
}

he_array_stats <- function(img) {
  sprintf(
    "dim=%s range=[%.4f, %.4f] mean=%.4f chroma=%.4f",
    paste(dim(img), collapse = "x"),
    min(img, na.rm = TRUE), max(img, na.rm = TRUE), mean(img, na.rm = TRUE),
    he_colorfulness(img)
  )
}

is_nearly_black <- function(img, mean_max = 0.02, q_max = 0.05) {
  mu <- mean(img, na.rm = TRUE)
  q <- suppressWarnings(as.numeric(stats::quantile(img, 0.99, na.rm = TRUE)))
  is.finite(mu) && mu < mean_max && is.finite(q) && q < q_max
}

# Real H&E is pink/purple. A failed JPEG/YCbCr decode is almost grayscale
# (R≈G≈B) with striping — reject those and try another page/backend.
he_colorfulness <- function(img) {
  if (length(dim(img)) != 3L || dim(img)[3] < 3L) return(0)
  r <- img[, , 1]
  g <- img[, , 2]
  b <- img[, , 3]
  mean(abs(r - g) + abs(g - b) + abs(r - b), na.rm = TRUE) / 2
}

is_gray_decode <- function(img, max_color = 0.025) {
  he_colorfulness(img) < max_color
}

he_looks_usable <- function(img) {
  !is.null(img) && !is_nearly_black(img) && !is_gray_decode(img)
}

# 8-bit / 16-bit / unit-interval. 16-bit scanner TIFFs are often 0–65535,
# or 8-bit data stored in the high byte (multiples of 256) — both look
# black if they are treated as 0–1 or divided by 255 and then clipped.
scale_to_unit <- function(img) {
  if (is.raw(img)) img <- as.integer(img)
  storage.mode(img) <- "double"
  mx <- max(img, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(pmin(pmax(img, 0), 1))
  if (mx > 255.5) {
    img <- img / 65535
  } else if (mx > 1.5) {
    img <- img / 255
  }
  q <- suppressWarnings(as.numeric(stats::quantile(img, 0.995, na.rm = TRUE)))
  if (is.finite(q) && q > 1e-6 && q < 0.15) img <- img / q
  pmin(pmax(img, 0), 1)
}

normalize_rgb_array <- function(img) {
  if (is.list(img) && !is.array(img)) {
    stop("Got a list of TIFF pages; pass a single page array")
  }
  if (length(dim(img)) == 2L) img <- replicate(3L, img)
  if (length(dim(img)) != 3L) stop("Expected a 2D or 3D image array")
  if (dim(img)[3] > 3L) img <- img[, , 1:3, drop = FALSE]
  scale_to_unit(img)
}

downsample_array <- function(img, max_px = 2000L) {
  max_px <- as.integer(max_px)
  h <- dim(img)[1]
  w <- dim(img)[2]
  long <- max(h, w)
  if (is.na(max_px) || max_px < 1L || long <= max_px) return(img)
  scale <- max_px / long
  nr <- max(1L, as.integer(round(h * scale)))
  nc <- max(1L, as.integer(round(w * scale)))
  ri <- as.integer(round(seq(1, h, length.out = nr)))
  ci <- as.integer(round(seq(1, w, length.out = nc)))
  img[ri, ci, , drop = FALSE]
}

magick_to_array <- function(im) {
  im <- magick::image_convert(im, colorspace = "sRGB", depth = 8)
  probe <- as.numeric(magick::image_data(magick::image_scale(im, "64x64"), channels = "rgb"))
  if (length(probe) && mean(probe, na.rm = TRUE) < 8) {
    im <- magick::image_normalize(im)
    im <- magick::image_convert(im, colorspace = "sRGB", depth = 8)
  }
  data <- magick::image_data(im, channels = "rgb")
  arr <- as.numeric(data) / 255
  dim(arr) <- dim(data)
  aperm(arr, c(3L, 2L, 1L))
}

he_identify_pages <- function(path) {
  fmt <- "%w %h %[colorspace] %C\\n"
  out <- character()
  if (nzchar(Sys.which("identify"))) {
    out <- suppressWarnings(system2(
      "identify", c("-quiet", "-format", fmt, path),
      stdout = TRUE, stderr = FALSE
    ))
  } else if (nzchar(Sys.which("magick"))) {
    out <- suppressWarnings(system2(
      "magick", c("identify", "-quiet", "-format", fmt, path),
      stdout = TRUE, stderr = FALSE
    ))
  }
  out <- out[nzchar(trimws(out))]
  if (!length(out)) return(NULL)
  rows <- strsplit(trimws(out), "[[:space:]]+")
  if (any(vapply(rows, length, 1L) < 2L)) return(NULL)
  data.frame(
    width = as.integer(vapply(rows, `[[`, "", 1L)),
    height = as.integer(vapply(rows, `[[`, "", 2L)),
    colorspace = vapply(rows, function(z) if (length(z) >= 3L) z[[3]] else "sRGB", ""),
    compression = vapply(rows, function(z) if (length(z) >= 4L) z[[4]] else "", ""),
    stringsAsFactors = FALSE
  )
}

pick_magick_page <- function(info, max_px = 2000L) {
  long <- pmax(info$width, info$height)
  short <- pmin(info$width, info$height)
  aspect <- long / pmax(short, 1)
  ok <- short >= 64 & aspect <= 4
  if ("colorspace" %in% names(info)) {
    rgbish <- grepl("RGB|sRGB|CMYK", info$colorspace, ignore.case = TRUE)
    if (any(ok & rgbish)) ok <- ok & rgbish
  }
  if (!any(ok)) ok <- rep(TRUE, nrow(info))
  target <- max(as.integer(max_px) * 2L, 1600L)
  which.min(abs(long - target) + ifelse(ok, 0, 1e9))
}

read_one_magick_page <- function(path, page, max_px = 2000L) {
  im <- magick::image_read(sprintf("%s[%d]", path, as.integer(page)))
  info1 <- magick::image_info(im)
  long <- max(info1$width[[1]], info1$height[[1]])
  if (long > max_px) {
    geom <- if (info1$width[[1]] >= info1$height[[1]]) {
      sprintf("%dx", max_px)
    } else {
      sprintf("x%d", max_px)
    }
    im <- magick::image_scale(im, geom)
  }
  normalize_rgb_array(magick_to_array(im))
}

# ImageMagick CLI resizes on disk and applies TIFF predictor/tiles correctly.
# The R magick binding's image_data() on a 29k×32k tiled TIFF produced the
# gray wavy 800×746 preview (right aspect ratio, wrong pixels).
read_he_imagemagick_cli <- function(path, max_px = 2000L) {
  bin <- Sys.which("magick")
  if (!nzchar(bin)) bin <- Sys.which("convert")
  if (!nzchar(bin)) return(NULL)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  geom <- sprintf("%dx%d>", as.integer(max_px), as.integer(max_px))
  message("ImageMagick CLI resize ", geom, " via ", bin)
  st <- suppressWarnings(system2(bin, c(path, "-resize", geom, tmp),
                                 stdout = TRUE, stderr = TRUE))
  if (!file.exists(tmp) || isTRUE(file.info(tmp)$size < 200)) {
    message(paste(st, collapse = "\n"))
    return(NULL)
  }
  arr <- read_preview_image_file(tmp)
  if (is.null(arr)) return(NULL)
  normalize_rgb_array(arr)
}

read_he_magick <- function(path, max_px = 2000L) {
  info <- he_identify_pages(path)
  if (!is.null(info) && nrow(info) == 1L &&
      max(info$width[[1]], info$height[[1]]) > 8000L) {
    message("magick R: skip in-memory decode of ",
            info$width[[1]], "x", info$height[[1]], " TIFF")
    return(NULL)
  }
  pages <- 0L
  if (!is.null(info) && nrow(info)) {
    long <- pmax(info$width, info$height)
    # Do not decode a 20k full-res page just to thumbnail it.
    keep <- which(long >= 64 & long <= 8000)
    if (!length(keep)) keep <- pick_magick_page(info, max_px)
    pages <- keep - 1L
    message("magick: probing pages ", paste(pages, collapse = ","))
  }
  best <- NULL
  best_score <- -Inf
  for (page in unique(pages)) {
    arr <- tryCatch(read_one_magick_page(path, page, max_px), error = function(e) NULL)
    if (!he_looks_usable(arr)) next
    sc <- he_colorfulness(arr)
    message("  page ", page, " colorfulness=", round(sc, 4))
    if (sc > best_score) {
      best <- arr
      best_score <- sc
    }
  }
  if (is.null(best)) {
    tryCatch(read_one_magick_page(path, pages[[1]], max_px), error = function(e) NULL)
  } else {
    best
  }
}

read_preview_image_file <- function(tmp) {
  ext <- tolower(tools::file_ext(tmp))
  if (ext %in% c("jpg", "jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
    return(jpeg::readJPEG(tmp, native = FALSE))
  }
  if (ext == "png" && requireNamespace("png", quietly = TRUE)) {
    return(png::readPNG(tmp, native = FALSE))
  }
  if (requireNamespace("magick", quietly = TRUE)) {
    return(magick_to_array(magick::image_read(tmp)))
  }
  NULL
}

vips_thumbnail_one <- function(src, max_px) {
  vips <- Sys.which("vips")
  if (!nzchar(vips)) return(NULL)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  args <- c("thumbnail", src, tmp, as.character(as.integer(max_px)))
  suppressWarnings(system2(vips, args, stdout = TRUE, stderr = TRUE))
  if (!file.exists(tmp) || isTRUE(file.info(tmp)$size < 200)) return(NULL)
  arr <- read_preview_image_file(tmp)
  if (is.null(arr)) return(NULL)
  normalize_rgb_array(arr)
}

read_he_vips <- function(path, max_px = 2000L) {
  if (!nzchar(Sys.which("vips"))) return(NULL)
  srcs <- path
  info <- he_identify_pages(path)
  if (!is.null(info) && nrow(info) > 1L) {
    long <- pmax(info$width, info$height)
    keep <- which(long >= 64 & long <= 8000)
    if (!length(keep)) keep <- seq_len(min(nrow(info), 6L))
    srcs <- c(
      sprintf("%s[page=%d]", path, keep - 1L),
      path
    )
  }
  best <- NULL
  best_score <- -Inf
  for (src in unique(srcs)) {
    arr <- tryCatch(vips_thumbnail_one(src, max_px), error = function(e) NULL)
    if (!he_looks_usable(arr)) next
    sc <- he_colorfulness(arr)
    message("vips ", src, " colorfulness=", round(sc, 4))
    if (sc > best_score) {
      best <- arr
      best_score <- sc
    }
  }
  best
}

he_thumbnail_py <- function() {
  cand <- c(
    "Code/he_thumbnail.py",
    file.path(getwd(), "Code/he_thumbnail.py"),
    file.path(getwd(), "he_thumbnail.py")
  )
  ofiles <- unlist(lapply(sys.frames(), function(fr) fr$ofile), use.names = FALSE)
  if (length(ofiles)) {
    cand <- c(file.path(dirname(ofiles), "he_thumbnail.py"), cand)
  }
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[[1]], mustWork = FALSE))
  NULL
}

read_he_python <- function(path, max_px = 2000L) {
  py <- Sys.which(c("python3", "python"))
  py <- py[nzchar(py)]
  if (!length(py)) return(NULL)
  script <- he_thumbnail_py()
  if (is.null(script)) return(NULL)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  st <- suppressWarnings(system2(
    py[[1]], c(script, path, tmp, as.character(as.integer(max_px))),
    stdout = TRUE, stderr = TRUE
  ))
  if (!file.exists(tmp) || isTRUE(file.info(tmp)$size < 200)) return(NULL)
  read_preview_image_file(tmp)
}

read_he_tiff_pkg <- function(path, max_px = 2000L) {
  if (!requireNamespace("tiff", quietly = TRUE)) return(NULL)
  info <- he_identify_pages(path)
  if (!is.null(info) && max(pmax(info$width, info$height), na.rm = TRUE) > 8000) {
    message("tiff R: skip full in-memory read of huge TIFF")
    return(NULL)
  }
  pages <- tryCatch(
    tiff::readTIFF(path, native = FALSE, all = TRUE),
    error = function(e) NULL
  )
  if (is.null(pages)) return(NULL)
  if (is.array(pages) || !is.list(pages)) pages <- list(pages)
  scored <- vapply(pages, function(p) {
    if (is.null(dim(p))) return(-Inf)
    arr <- normalize_rgb_array(p)
    if (is_nearly_black(arr)) return(-Inf)
    mean(arr, na.rm = TRUE) * max(dim(arr)[1:2])
  }, numeric(1))
  if (!any(is.finite(scored))) {
    pages[[1]]
  } else {
    pages[[which.max(scored)]]
  }
}

# Whole-slide HE.tif files are tiled / pyramidal / often 16-bit. The R tiff
# package and a naive magick read frequently return the first IFD (a black
# label, mask, or single tile). Prefer vips/PIL thumbnail, then magick with
# a mid-pyramid page, then tiff.
read_he_array <- function(path, max_px = 2000L) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop("Image not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  img <- NULL
  used <- NULL
  if (ext %in% c("jpg", "jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
    img <- jpeg::readJPEG(path, native = FALSE)
    used <- "jpeg"
  } else if (ext == "png" && requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(path, native = FALSE)
    used <- "png"
  } else {
    try_backend <- function(fun, name) {
      got <- tryCatch(fun(), error = function(e) NULL)
      if (is.null(got)) return(NULL)
      got <- downsample_array(normalize_rgb_array(got), max_px = max_px)
      if (is_nearly_black(got)) {
        message(name, " returned a near-black array; trying the next reader")
        return(NULL)
      }
      if (is_gray_decode(got)) {
        message(
          name, " returned grayscale (colorfulness=",
          round(he_colorfulness(got), 4), "); trying the next reader"
        )
        return(NULL)
      }
      list(img = got, used = name)
    }
    info <- he_identify_pages(path)
    if (!is.null(info)) {
      message(
        "identify: ", nrow(info), " page(s); ",
        paste0(info$width, "x", info$height, " ", info$colorspace, collapse = "; ")
      )
    }
    hit <- try_backend(function() read_he_vips(path, max_px), "vips")
    if (is.null(hit)) {
      hit <- try_backend(function() read_he_imagemagick_cli(path, max_px), "magick-cli")
    }
    if (is.null(hit)) hit <- try_backend(function() read_he_python(path, max_px), "python")
    if (is.null(hit) && requireNamespace("magick", quietly = TRUE)) {
      hit <- try_backend(function() read_he_magick(path, max_px = max_px), "magick")
    }
    if (is.null(hit) && ext %in% c("tif", "tiff")) {
      hit <- try_backend(function() read_he_tiff_pkg(path, max_px), "tiff")
    }
    if (!is.null(hit)) {
      img <- hit$img
      used <- hit$used
    }
  }
  if (is.null(img)) {
    stop(
      "Cannot read ", path,
      ". Install libvips (`vips thumbnail`) or: pip install --user pillow tifffile"
    )
  }
  if (used %in% c("jpeg", "png")) {
    img <- downsample_array(normalize_rgb_array(img), max_px = max_px)
  }
  message("H&E via ", used, " — ", he_array_stats(img))
  if (is_nearly_black(img)) {
    warning(
      "H&E array is nearly black (", he_array_stats(img), "). ",
      "This TIFF is likely a pyramid/WSI; install vips or `pip install --user pillow tifffile`.",
      call. = FALSE
    )
  }
  img
}

rotate_roi_jpg <- function(img, rotate = 0L, flip_x = FALSE, flip_y = FALSE) {
  rotate <- as.integer(rotate) %% 360L
  if (!rotate %in% c(0L, 90L, 180L, 270L)) {
    stop("rotate must be 0, 90, 180, or 270 (clockwise)")
  }
  h <- dim(img)[1]
  w <- dim(img)[2]
  if (isTRUE(flip_x)) img <- img[, w:1, , drop = FALSE]
  if (isTRUE(flip_y)) img <- img[h:1, , , drop = FALSE]
  h <- dim(img)[1]
  w <- dim(img)[2]
  if (rotate == 90L) {
    img <- aperm(img, c(2L, 1L, 3L))[, h:1, , drop = FALSE]
  } else if (rotate == 180L) {
    img <- img[h:1, w:1, , drop = FALSE]
  } else if (rotate == 270L) {
    img <- aperm(img, c(2L, 1L, 3L))[w:1, , , drop = FALSE]
  }
  img
} # end rotate_roi_jpg

parse_dbit_barcodes <- function(cells) {
  n <- length(cells)
  row <- rep(NA_integer_, n)
  col <- rep(NA_integer_, n)
  m <- regexec("([0-9]+)x([0-9]+)$", cells)
  hit <- vapply(m, function(z) z[[1]] != -1L, logical(1))
  if (any(hit)) {
    mm <- regmatches(cells, m)
    row[hit] <- as.integer(vapply(mm[hit], `[[`, "", 2L))
    col[hit] <- as.integer(vapply(mm[hit], `[[`, "", 3L))
  }
  miss <- is.na(row)
  if (any(miss)) {
    m2 <- regexec("([0-9]+)_([0-9]+)$", cells[miss])
    hit2 <- vapply(m2, function(z) z[[1]] != -1L, logical(1))
    if (any(hit2)) {
      mm2 <- regmatches(cells[miss], m2)
      idx <- which(miss)[hit2]
      row[idx] <- as.integer(vapply(mm2[hit2], `[[`, "", 2L))
      col[idx] <- as.integer(vapply(mm2[hit2], `[[`, "", 3L))
    }
  }
  data.frame(row = row, col = col, row.names = cells, check.names = FALSE)
}

infer_ndim <- function(cells, default = 50L) {
  rc <- parse_dbit_barcodes(cells)
  ok <- !is.na(rc$row) & !is.na(rc$col)
  if (!any(ok)) return(as.integer(default))
  as.integer(max(c(rc$row[ok], rc$col[ok])))
}

resolve_ndim <- function(sid, cells) {
  if (exists("seu_params") && is.data.frame(seu_params) &&
      all(c("sid", "ndim") %in% names(seu_params))) {
    hit <- seu_params$ndim[seu_params$sid == sid]
    if (length(hit) && !is.na(hit[[1]])) return(as.integer(hit[[1]]))
  }
  if (exists("HE_NDIM") && sid %in% names(HE_NDIM) && !is.na(HE_NDIM[[sid]])) {
    return(as.integer(HE_NDIM[[sid]]))
  }
  infer_ndim(cells)
}

# Map DBiT / SlideSeq barcodes onto a VisiumV1 coordinate frame so
# SpatialDimPlot(images = "roi") matches images = "image".
dbit_roi_coords <- function(cells, img, ndim) {
  rc <- parse_dbit_barcodes(cells)
  ok <- !is.na(rc$row) & !is.na(rc$col)
  if (!any(ok)) {
    stop("No rowxcol (or row_col) barcodes — cannot place spots on the H&E")
  }
  row <- rc$row[ok]
  col <- rc$col[ok]
  n_px <- 2 * ndim - 1
  pixel_h <- dim(img)[1] / n_px
  pixel_w <- dim(img)[2] / n_px
  spot_px <- min(pixel_h, pixel_w)
  coords <- data.frame(
    tissue = 1L, row = row, col = col,
    imagecol = 2 * (row - 1) * pixel_w + pixel_w / 2,
    imagerow = 2 * (ndim - col) * pixel_h + pixel_h / 2,
    row.names = cells[ok]
  )
  list(coords = coords, spot_px = spot_px)
}

# Add histology as image "roi"; leave the original SlideSeq "image" unchanged.
add_dbit_roi_image <- function(seu, roi_jpg, assay = "long", ndim = 50L,
                               rotate = 0L, flip_x = FALSE, flip_y = FALSE,
                               max_px = 2000L, img = NULL) {
  if (is.null(img)) {
    if (is.na(roi_jpg) || !nzchar(roi_jpg) || !file.exists(roi_jpg)) {
      warning("No H&E image — skip", call. = FALSE)
      return(seu)
    }
    img <- read_he_array(roi_jpg, max_px = max_px)
  }
  img <- rotate_roi_jpg(img, rotate = rotate, flip_x = flip_x, flip_y = flip_y)
  placed <- dbit_roi_coords(Seurat::Cells(seu), img, ndim)
  seu[["roi"]] <- methods::new(
    Class = "VisiumV1",
    image = img,
    scale.factors = Seurat::scalefactors(
      spot = placed$spot_px, fiducial = placed$spot_px, hires = 1, lowres = 1
    ),
    coordinates = placed$coords,
    spot.radius = placed$spot_px / max(dim(img)[1:2]),
    assay = assay,
    key = "roi_"
  )
  seu
}

rel_to_root <- function(path, root) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(NA_character_)
  if (!is.na(root) && nzchar(root) && startsWith(path, root)) {
    return(sub("^/+", "", substr(path, nchar(root) + 1L, nchar(path))))
  }
  path
}

print_he_inventory <- function(sids, he_img, root) {
  missing <- is.na(he_img) | !nzchar(as.character(he_img))
  status <- ifelse(missing, "missing", "found")
  rel <- vapply(he_img, rel_to_root, character(1), root = root, USE.NAMES = FALSE)
  tab <- data.frame(sid = sids, status = status, file = rel, stringsAsFactors = FALSE)
  print(tab, right = FALSE, row.names = FALSE)
  files <- list_he_image_files(root)
  used <- he_img[status == "found"]
  extra <- setdiff(
    normalizePath(files, mustWork = FALSE),
    normalizePath(used, mustWork = FALSE)
  )
  if (length(extra)) {
    message("Unmatched image files (not assigned to a sample):")
    print(basename(extra))
  }
  invisible(tab)
}

add_all_he_images <- function(seu,
                              root = default_he_path(),
                              he_img = NULL,
                              roi_orient = NULL,
                              max_px = HE_MAX_PX,
                              plot = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to attach H&E images")
  }
  sids <- names(seu)
  if (is.null(he_img)) he_img <- find_he_files(sids, root)
  if (is.null(roi_orient)) roi_orient <- default_roi_orient(sids)
  message("PATH_HE_IMG = ", root)
  print_he_inventory(sids, he_img, root)

  for (sid in sids) {
    o <- roi_orient[[sid]]
    if (is.null(o)) o <- list(rotate = 0L, flip_x = FALSE, flip_y = FALSE)
    ndim <- resolve_ndim(sid, Seurat::Cells(seu[[sid]]))
    assay_use <- if ("long" %in% Seurat::Assays(seu[[sid]])) {
      "long"
    } else {
      Seurat::DefaultAssay(seu[[sid]])
    }
    message(sid, " ndim=", ndim, " assay=", assay_use, " file=", he_img[[sid]])
    seu[[sid]] <- add_dbit_roi_image(
      seu[[sid]], roi_jpg = he_img[[sid]], assay = assay_use, ndim = ndim,
      rotate = o$rotate, flip_x = o$flip_x, flip_y = o$flip_y, max_px = max_px
    )
  }

  if (isTRUE(plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE) ||
        !requireNamespace("patchwork", quietly = TRUE)) {
      warning("ggplot2 + patchwork needed for SpatialDimPlot wrap; skip plot")
      return(seu)
    }
    if (!exists("sci_cell_12")) {
      sci_cell_12 <- Seurat::DiscretePalette(12, palette = "alphabet")
    }
    has_roi <- vapply(seu, function(x) "roi" %in% Seurat::Images(x), logical(1))
    if (any(has_roi)) {
      print(patchwork::wrap_plots(lapply(sids[has_roi], function(sid) {
        Seurat::SpatialDimPlot(seu[[sid]], images = "roi", cols = sci_cell_12) +
          ggplot2::ggtitle(sid)
      }), ncol = 3))
    }
    if (any(!has_roi)) {
      message("No roi image: ", paste(sids[!has_roi], collapse = ", "))
    }
  }
  seu
}

plot_he_array <- function(img, main = "") {
  if (length(dim(img)) == 2L) img <- replicate(3L, img)
  h <- dim(img)[1]
  w <- dim(img)[2]
  stats <- he_array_stats(img)
  if (!nzchar(main)) main <- stats else main <- paste0(main, "\n", stats)
  ras <- grDevices::as.raster(pmin(pmax(img, 0), 1))
  op <- graphics::par(mar = c(1, 1, 3.2, 1))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(
    0, type = "n", xlim = c(0, w), ylim = c(0, h),
    xlab = "", ylab = "", axes = FALSE, asp = 1, main = main
  )
  graphics::rasterImage(ras, 0, 0, w, h, interpolate = TRUE)
  invisible(img)
}

pick_preview_sid <- function(he_img, preferred = c("LN18427", "LN4737", "LN14734")) {
  have <- names(he_img)[!is.na(he_img) & nzchar(as.character(he_img))]
  if (!length(have)) stop("No H&E file is available to preview")
  hit <- preferred[preferred %in% have]
  if (length(hit)) hit[[1]] else have[[1]]
}

# Read one H&E (800 px) and plot it. If `seu` is present, attach only that
# sample as image "roi" and draw SpatialDimPlot(roi) next to the original
# SlideSeq image so orientation can be checked.
quick_he_preview <- function(seu = get0("seu"),
                             sid = NULL,
                             root = default_he_path(),
                             he_img = NULL,
                             max_px = 800L,
                             attach = TRUE) {
  sids <- if (is.list(seu) && length(seu)) names(seu) else AICL_SAMPLE_IDS
  if (is.null(he_img)) he_img <- find_he_files(sids, root)
  if (is.null(sid)) sid <- pick_preview_sid(he_img)
  path <- he_img[[sid]]
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    stop("No H&E file for ", sid)
  }
  message("Preview ", sid, "\n  ", path)
  img <- read_he_array(path, max_px = max_px)
  plot_he_array(img, main = paste0(sid, " H&E"))
  if (is_gray_decode(img) || is_nearly_black(img)) {
    warning(
      "Decoded H&E does not look like histology. Run: identify ", path,
      call. = FALSE
    )
  }

  if (isTRUE(attach) && is.list(seu) && sid %in% names(seu)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      warning("Seurat not available — raw H&E only")
      return(invisible(list(sid = sid, path = path, image = img, seu = seu)))
    }
    o <- list(rotate = 0L, flip_x = FALSE, flip_y = FALSE)
    if (exists("roi_orient") && sid %in% names(roi_orient)) o <- roi_orient[[sid]]
    ndim <- resolve_ndim(sid, Seurat::Cells(seu[[sid]]))
    assay_use <- if ("long" %in% Seurat::Assays(seu[[sid]])) {
      "long"
    } else {
      Seurat::DefaultAssay(seu[[sid]])
    }
    message("Attaching roi for ", sid, " (overwrite, ndim=", ndim, ")")
    seu[[sid]] <- add_dbit_roi_image(
      seu[[sid]], roi_jpg = path, assay = assay_use, ndim = ndim,
      rotate = o$rotate, flip_x = o$flip_x, flip_y = o$flip_y, max_px = max_px,
      img = img
    )
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        requireNamespace("patchwork", quietly = TRUE) &&
        "roi" %in% Seurat::Images(seu[[sid]])) {
      cols <- if (exists("sci_cell_12")) sci_cell_12 else NULL
      p_he <- Seurat::SpatialDimPlot(seu[[sid]], images = "roi", cols = cols) +
        ggplot2::ggtitle(paste(sid, "H&E roi"))
      plots <- list(p_he)
      if ("image" %in% Seurat::Images(seu[[sid]])) {
        plots <- c(plots, list(
          Seurat::SpatialDimPlot(seu[[sid]], images = "image", cols = cols) +
            ggplot2::ggtitle(paste(sid, "SlideSeq image"))
        ))
      }
      print(patchwork::wrap_plots(plots, ncol = 2))
    }
  }
  invisible(list(sid = sid, path = path, image = img, seu = seu))
}

# -----------------------------------------------------------------------------
# Driver: map files, then preview one H&E. Set AICL_HE_ATTACH_ALL <- TRUE
# before source() to attach every sample instead.
# -----------------------------------------------------------------------------
if (exists("seu") && is.list(seu) && length(seu) &&
    !isTRUE(get0(".AICL_HE_SKIP_ATTACH", ifnotfound = FALSE))) {
  PATH_HE_IMG <- default_he_path()
  SAMPLE_IDS <- names(seu)
  HE_IMG <- find_he_files(SAMPLE_IDS, PATH_HE_IMG)
  print(HE_IMG)
  print_he_inventory(SAMPLE_IDS, HE_IMG, PATH_HE_IMG)

  roi_orient <- default_roi_orient(SAMPLE_IDS)
  # Examples after QC against SpatialDimPlot(..., images = "image"):
  # roi_orient$LN14734$flip_x <- TRUE
  # roi_orient$LN14734$flip_y <- TRUE
  # roi_orient$LN18427$rotate <- 90L

  if (isTRUE(get0("AICL_HE_ATTACH_ALL", ifnotfound = FALSE))) {
    seu <- add_all_he_images(
      seu, root = PATH_HE_IMG, he_img = HE_IMG,
      roi_orient = roi_orient, max_px = HE_MAX_PX, plot = TRUE
    )
  } else {
    preview <- quick_he_preview(
      seu, sid = pick_preview_sid(HE_IMG), root = PATH_HE_IMG,
      he_img = HE_IMG, max_px = 800L, attach = TRUE
    )
    seu <- preview$seu
    message(
      "Preview only. After orientation looks right, attach every slide with:\n",
      "  AICL_HE_ATTACH_ALL <- TRUE\n",
      "  seu <- add_all_he_images(seu)"
    )
  }
}
