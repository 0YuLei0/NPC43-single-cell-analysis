#!/usr/bin/env python3
"""Downsample a large / pyramidal HE TIFF to an RGB JPEG thumbnail."""
from __future__ import annotations

import sys


def _save_rgb(im, dst: str, max_px: int) -> None:
    im = im.convert("RGB")
    im.thumbnail((max_px, max_px))
    im.save(dst, format="JPEG", quality=92)


def _score_rgb(im) -> float:
    thumb = im.convert("RGB").resize((64, 64))
    px = list(thumb.getdata())
    mean = sum(sum(p) for p in px) / (len(px) * 3.0 * 255.0)
    w, h = im.size
    if min(w, h) < 32:
        return -1.0
    return mean * float(min(w, h))


def via_openslide(src: str, dst: str, max_px: int) -> bool:
    try:
        import openslide
    except Exception:
        return False
    slide = openslide.OpenSlide(src)
    try:
        dims = slide.level_dimensions
        # lowest-res level whose long edge is still >= max_px, else the last
        pick = len(dims) - 1
        for i, (w, h) in enumerate(dims):
            if max(w, h) >= max_px:
                pick = i
        w, h = dims[pick]
        im = slide.read_region((0, 0), pick, (w, h)).convert("RGB")
        _save_rgb(im, dst, max_px)
        return True
    finally:
        slide.close()


def via_tifffile(src: str, dst: str, max_px: int) -> bool:
    try:
        import numpy as np
        import tifffile
        from PIL import Image
    except Exception:
        return False
    with tifffile.TiffFile(src) as tif:
        series = max(tif.series, key=lambda s: int(np.prod(s.shape)))
        levels = list(getattr(series, "levels", [series])) or [series]
        level = levels[-1]
        for lv in reversed(levels):
            shape = tuple(int(x) for x in lv.shape)
            spatial = [d for d in shape if d > 4]
            if spatial and max(spatial) >= max_px:
                level = lv
                break
        arr = np.asarray(level.asarray())
        if arr.ndim == 2:
            arr = np.stack([arr, arr, arr], axis=-1)
        if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
            arr = np.moveaxis(arr, 0, -1)
        if arr.shape[-1] > 3:
            arr = arr[..., :3]
        im = Image.fromarray(arr)
    _save_rgb(im, dst, max_px)
    return True


def via_pillow(src: str, dst: str, max_px: int) -> bool:
    try:
        from PIL import Image
    except Exception:
        return False
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(src)
    n = int(getattr(im, "n_frames", 1) or 1)
    best = None
    best_score = -1.0
    for i in range(min(n, 16)):
        try:
            im.seek(i)
            frame = im.convert("RGB")
        except Exception:
            continue
        score = _score_rgb(frame)
        if score > best_score:
            best = frame.copy()
            best_score = score
    if best is None:
        best = im.convert("RGB")
    _save_rgb(best, dst, max_px)
    return True


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        sys.stderr.write("usage: he_thumbnail.py SRC DST [MAX_PX]\n")
        return 2
    src, dst = argv[1], argv[2]
    max_px = int(argv[3]) if len(argv) > 3 else 800
    for fn in (via_openslide, via_tifffile, via_pillow):
        try:
            if fn(src, dst, max_px):
                return 0
        except Exception as exc:
            sys.stderr.write(f"{fn.__name__}: {exc}\n")
    sys.stderr.write("no reader could decode the TIFF\n")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
