# Base-R tests for AICL H&E file matching and geometry helpers.
# Run: Rscript Code/test_AICL_add_HE_image.R

.AICL_HE_SKIP_ATTACH <- TRUE
src <- if (file.exists("Code/AICL_add_HE_image.R")) {
  "Code/AICL_add_HE_image.R"
} else if (file.exists("AICL_add_HE_image.R")) {
  "AICL_add_HE_image.R"
} else {
  stop("Cannot find AICL_add_HE_image.R")
}
sys.source(src, envir = environment(), keep.source = TRUE)

parsed <- parse(src, keep.source = TRUE)
stopifnot(length(parsed) >= 1L)

n_ok <- 0L
check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  n_ok <<- n_ok + 1L
  message("OK ", msg)
}

# --- mock the 03.HE listing from login1.bouchet ---
root <- tempfile("aicl-he-")
dir.create(file.path(root, "v2"), recursive = TRUE)
dir.create(file.path(root, "v3"), recursive = TRUE)
known_tif <- c(
  "FFPE_LN14734.HE.tif",
  "FFPE_LN18427.HE.tif",
  "FFPE_LN21720.HE.tif",
  "FFPE_LN22630.HE.tif",
  "FFPE_LN32298.HE.tif",
  "FFPE_LN3470.HE.tif",
  "FFPE_LN4737.HE.tif"
)
scripts <- c(
  "FFPE_ALCL_LN18427.work.sh",
  "FFPE_ALCL_LN22630.work.sh",
  "FFPE_LN14734.work.sh",
  "FFPE_LN21720.work.sh",
  "FFPE_LN32298.work.sh",
  "FFPE_LN3470.work.sh",
  "FFPE_LN4737.work.sh"
)
invisible(file.create(file.path(root, c(known_tif, scripts, "LN34700.HE.tif"))))

sids <- AICL_SAMPLE_IDS
found <- find_he_files(sids, root)

check(identical(unname(found[["LN14734"]]), file.path(root, "FFPE_LN14734.HE.tif")),
      "LN14734 maps to FFPE_LN14734.HE.tif")
check(identical(unname(found[["LN18427"]]), file.path(root, "FFPE_LN18427.HE.tif")),
      "LN18427 maps to FFPE_LN18427.HE.tif")
check(identical(unname(found[["LN3470"]]), file.path(root, "FFPE_LN3470.HE.tif")),
      "LN3470 does not pick LN34700.HE.tif")
check(is.na(found[["LN3524"]]), "LN3524 missing at top level")
check(is.na(found[["LN5034"]]), "LN5034 missing at top level")
check(sum(!is.na(found)) == 7L, "seven top-level HE TIFFs assigned")
inv <- print_he_inventory(sids, found, root)
check(identical(inv$status[inv$sid == "LN3524"], "missing") &&
        identical(inv$status[inv$sid == "LN5034"], "missing"),
      "inventory prints missing LN3524/LN5034 without error")
check(identical(inv$status[inv$sid == "LN14734"], "found"),
      "inventory marks LN14734 found")

# Recover missing slides from v3 / v2 (the live folders on Bouchet).
invisible(file.create(c(
  file.path(root, "v3", "LN3524.HE.tif"),
  file.path(root, "v2", "FFPE_LN5034.HE.jpg")
)))
found2 <- find_he_files(sids, root)
check(identical(unname(found2[["LN3524"]]), file.path(root, "v3", "LN3524.HE.tif")),
      "LN3524 recovered from v3")
check(identical(unname(found2[["LN5034"]]), file.path(root, "v2", "FFPE_LN5034.HE.jpg")),
      "LN5034 recovered from v2")

# Exact top-level HE wins over a nested extra.
invisible(file.create(file.path(root, "v3", "LN14734_preview.jpg")))
check(identical(unname(find_he_file("LN14734", root)),
                file.path(root, "FFPE_LN14734.HE.tif")),
      "exact FFPE_*.HE.tif beats nested preview")

# ALCL-prefixed HE is used only when the plain FFPE_*.HE file is absent.
unlink(file.path(root, "FFPE_LN18427.HE.tif"))
invisible(file.create(file.path(root, "FFPE_ALCL_LN18427.HE.tif")))
check(identical(unname(find_he_file("LN18427", root)),
                file.path(root, "FFPE_ALCL_LN18427.HE.tif")),
      "ALCL HE used when the plain HE is missing")

# work.sh must never be chosen
check(!any(grepl("\\.sh$", list_he_candidates("LN18427", root))),
      "scripts are not image candidates")

# --- geometry ---
img <- array(0, dim = c(4, 6, 3))
img[1, 1, 1] <- 1
r90 <- rotate_roi_jpg(img, rotate = 90L)
check(identical(dim(r90), c(6L, 4L, 3L)), "90° rotation swaps height/width")
check(isTRUE(all.equal(r90[1, 4, 1], 1)), "90° clockwise moves (1,1) to (1,4)")
r180 <- rotate_roi_jpg(img, rotate = 180L)
check(isTRUE(all.equal(r180[4, 6, 1], 1)), "180° moves (1,1) to (h,w)")
fx <- rotate_roi_jpg(img, flip_x = TRUE)
check(isTRUE(all.equal(fx[1, 6, 1], 1)), "flip_x mirrors columns")

tiny <- array(runif(20 * 10 * 3), dim = c(20, 10, 3))
ds <- downsample_array(tiny, max_px = 10L)
check(dim(ds)[1] == 10L && dim(ds)[2] == 5L, "downsample keeps aspect ratio")
check(identical(downsample_array(tiny, 50L), tiny), "no upsample")

rc <- parse_dbit_barcodes(c("1x2", "LN14734_50x3", "12_8", "bad"))
check(identical(rc$row[1:3], c(1L, 50L, 12L)) && identical(rc$col[1:3], c(2L, 3L, 8L)),
      "parse rowxcol, prefix_rowxcol, and row_col")
check(is.na(rc$row[4]), "unparseable barcode is NA")
check(identical(infer_ndim(c("1x2", "48x3", "12x50")), 50L), "infer ndim from max row/col")
check(identical(infer_ndim("nope"), 50L), "infer ndim fallback")

cells <- c("1x50", "50x1")
img50 <- array(0, dim = c(99, 99, 3))
placed <- dbit_roi_coords(cells, img50, ndim = 50L)
check(nrow(placed$coords) == 2L, "coords for both spots")
# col=50 -> imagerow = 2*(50-50)*h + h/2 = pixel_h/2 (top)
# row=50 -> imagecol = 2*(50-1)*w + w/2 (right)
check(placed$coords["1x50", "imagerow"] < placed$coords["50x1", "imagerow"],
      "high barcode col maps toward the top of the H&E")
check(placed$coords["50x1", "imagecol"] > placed$coords["1x50", "imagecol"],
      "high barcode row maps toward the right of the H&E")

unlink(root, recursive = TRUE)
message("All ", n_ok, " checks passed")
