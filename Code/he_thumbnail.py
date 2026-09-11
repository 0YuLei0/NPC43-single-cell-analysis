#!/usr/bin/env python3
"""Downsample a large / pyramidal HE TIFF to an RGB thumbnail.

Prefers OpenSlide's get_thumbnail (correct Aperio/JPEG tiles + YCbCr).
Scores candidates by colorfulness so a gray JPEG-decode failure loses to
real pink/purple H&E.
"""
from __future__ import annotations

import sys


def _save_rgb(im, dst: str, max_px: int) -> None:
    im = im.convert("RGB")
    im.thumbnail((max_px, max_px))
    ext = dst.lower().rsplit(".", 1)[-1]
    if ext == "png":
        im.save(dst, format="PNG")
    else:
        im.save(dst, format="JPEG", quality=92)


def _colorfulness(im) -> float:
    thumb = im.convert("RGB").resize((64, 64))
    px = list(thumb.getdata())
    if not px:
        return -1.0
    acc = 0.0
    for r, g, b in px:
        acc += abs(r - g) + abs(g - b) + abs(r - b)
    return acc / (len(px) * 2.0 * 255.0)


def _score_rgb(im) -> float:
    w, h = im.size
    if min(w, h) < 32:
        return -1.0
    chroma = _colorfulness(im)
    # Reject near-grayscale JPEG/YCbCr failures (typical chroma ~ 0).
    if chroma < 0.025:
        return -1.0
    return chroma * float(min(w, h))


def via_openslide(src: str, dst: str, max_px: int) -> bool:
    try:
        import openslide
    except Exception:
        return False
    slide = openslide.OpenSlide(src)
    try:
        im = slide.get_thumbnail((max_px, max_px)).convert("RGB")
        if _score_rgb(im) < 0:
            return False
        _save_rgb(im, dst, max_px)
        return True
    finally:
        slide.close()


def _ycbcr_to_rgb(arr):
    import numpy as np

    arr = arr.astype(np.float32)
    y = arr[..., 0]
    cb = arr[..., 1] - 128.0
    cr = arr[..., 2] - 128.0
    r = y + 1.402 * cr
    g = y - 0.344136 * cb - 0.714136 * cr
    b = y + 1.772 * cb
    return np.clip(np.stack([r, g, b], axis=-1), 0, 255).astype("uint8")


def _as_hwc_rgb(arr, photometric=None):
    import numpy as np

    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = np.stack([arr, arr, arr], axis=-1)
    if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
        arr = np.moveaxis(arr, 0, -1)
    if arr.ndim != 3:
        raise ValueError(f"unexpected shape {arr.shape}")
    if arr.shape[-1] > 3:
        arr = arr[..., :3]
    name = ""
    if photometric is not None:
        name = str(photometric).lower()
    if "ycbcr" in name or photometric in (6, "YCBCR"):
        if arr.dtype != "uint8":
            arr = np.clip(arr, 0, 255).astype("uint8")
        arr = _ycbcr_to_rgb(arr)
    return arr


def via_tifffile(src: str, dst: str, max_px: int) -> bool:
    try:
        import numpy as np
        import tifffile
        from PIL import Image
    except Exception:
        return False
    best = None
    best_score = -1.0
    with tifffile.TiffFile(src) as tif:
        for series in tif.series:
            levels = list(getattr(series, "levels", [series])) or [series]
            photo = getattr(series, "photometric", None)
            for lv in reversed(levels):
                shape = tuple(int(x) for x in lv.shape)
                spatial = [d for d in shape if d > 4]
                if spatial and max(spatial) > 12000:
                    continue
                try:
                    arr = _as_hwc_rgb(lv.asarray(), photometric=photo)
                except Exception:
                    continue
                if arr.dtype != np.uint8:
                    mx = float(arr.max()) if arr.size else 1.0
                    if mx > 255.5:
                        arr = np.clip(arr / 65535.0 * 255.0, 0, 255).astype("uint8")
                    elif mx > 1.5:
                        arr = np.clip(arr, 0, 255).astype("uint8")
                    else:
                        arr = np.clip(arr * 255.0, 0, 255).astype("uint8")
                im = Image.fromarray(arr, mode="RGB")
                score = _score_rgb(im)
                if score > best_score:
                    best = im
                    best_score = score
    if best is None or best_score < 0:
        return False
    _save_rgb(best, dst, max_px)
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
    if best is None or best_score < 0:
        return False
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
    sys.stderr.write("no colorful H&E page found (decode may be YCbCr/gray)\n")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
