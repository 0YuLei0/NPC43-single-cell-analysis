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

normalize_rgb_array <- function(img) {
  if (length(dim(img)) == 2L) img <- replicate(3L, img)
  if (length(dim(img)) != 3L) stop("Expected a 2D or 3D image array")
  if (dim(img)[3] > 3L) img <- img[, , 1:3, drop = FALSE]
  if (is.raw(img) || is.integer(img) || max(img, na.rm = TRUE) > 1.5) {
    img <- img / 255
  }
  storage.mode(img) <- "double"
  pmin(pmax(img, 0), 1)
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

read_he_magick <- function(path, max_px = 2000L) {
  im <- magick::image_read(path)
  info <- magick::image_info(im)
  if (nrow(info) > 1L) {
    long <- pmax(info$width, info$height)
    ge <- which(long >= max_px)
    idx <- if (length(ge)) ge[which.min(long[ge])] else which.max(long)
    im <- im[idx]
    info <- magick::image_info(im)
  }
  long <- max(info$width[[1]], info$height[[1]])
  if (long > max_px) {
    geom <- if (info$width[[1]] >= info$height[[1]]) {
      sprintf("%dx", max_px)
    } else {
      sprintf("x%d", max_px)
    }
    im <- magick::image_scale(im, geom)
  }
  im <- magick::image_convert(im, colorspace = "sRGB")
  data <- magick::image_data(im, channels = "rgb")
  arr <- as.numeric(data) / 255
  dim(arr) <- dim(data)
  # magick::image_data is [channel, width, height] -> [row, col, channel]
  aperm(arr, c(3L, 2L, 1L))
}

read_he_array <- function(path, max_px = 2000L) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop("Image not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (requireNamespace("magick", quietly = TRUE)) {
    return(normalize_rgb_array(read_he_magick(path, max_px = max_px)))
  }
  img <- NULL
  if (ext %in% c("jpg", "jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
    img <- jpeg::readJPEG(path, native = FALSE)
  } else if (ext == "png" && requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(path, native = FALSE)
  } else if (ext %in% c("tif", "tiff") && requireNamespace("tiff", quietly = TRUE)) {
    img <- tiff::readTIFF(path, native = FALSE)
  } else {
    stop(
      "Cannot read ", path,
      ". Install magick (recommended for these ~GB TIFFs) or tiff/jpeg/png."
    )
  }
  downsample_array(normalize_rgb_array(img), max_px = max_px)
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
                               max_px = 2000L) {
  if (is.na(roi_jpg) || !nzchar(roi_jpg) || !file.exists(roi_jpg)) {
    warning("No H&E image — skip", call. = FALSE)
    return(seu)
  }
  img <- rotate_roi_jpg(
    read_he_array(roi_jpg, max_px = max_px),
    rotate = rotate, flip_x = flip_x, flip_y = flip_y
  )
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

# -----------------------------------------------------------------------------
# Driver: runs when this file is sourced in a session that already has `seu`.
# -----------------------------------------------------------------------------
if (exists("seu") && is.list(seu) && length(seu) &&
    !isTRUE(get0(".AICL_HE_SKIP_ATTACH", ifnotfound = FALSE))) {
  PATH_HE_IMG <- default_he_path()
  SAMPLE_IDS <- names(seu)
  HE_IMG <- find_he_files(SAMPLE_IDS, PATH_HE_IMG)
  print(HE_IMG)

  roi_orient <- default_roi_orient(SAMPLE_IDS)
  # Examples after QC against SpatialDimPlot(..., images = "image"):
  # roi_orient$LN14734$flip_x <- TRUE
  # roi_orient$LN14734$flip_y <- TRUE
  # roi_orient$LN18427$rotate <- 90L

  seu <- add_all_he_images(
    seu, root = PATH_HE_IMG, he_img = HE_IMG,
    roi_orient = roi_orient, max_px = HE_MAX_PX, plot = TRUE
  )
}
