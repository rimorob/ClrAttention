# Cheap fidelity test: locate each matrix via its grey self-diagonal cells, sample the 14x14 cell colours,
# compare cell-wise (darkness = 1 - luminance, rescaled 0-1) with the real reference matrix.
import sys, numpy as np
from PIL import Image
from scipy import ndimage as ndi
def grey_boxes(im, region=None):
    a = np.asarray(im.convert('RGB')).astype(int)
    mx, mn = a.max(2), a.min(2); mean = a.mean(2)
    m = ((mx - mn) < 40) & (mean > 150) & (mean < 243)
    if region: 
        x0,y0,x1,y1 = region; mm = np.zeros_like(m); mm[y0:y1,x0:x1] = m[y0:y1,x0:x1]; m = mm
    lab, n = ndi.label(m); out = []
    for k, sl in enumerate(ndi.find_objects(lab), 1):
        h = sl[0].stop-sl[0].start; w = sl[1].stop-sl[1].start
        if 15 < h < 400 and 0.7 < w/h < 1.3 and (lab[sl]==k).sum() > 0.7*h*w: out.append((sl[1].start, sl[0].start, sl[1].stop, sl[0].stop))
    return out
def grid_from_boxes(b, n=14):
    b = sorted(b, key=lambda t: t[0])
    cw = np.median([t[2]-t[0] for t in b]); ch = np.median([t[3]-t[1] for t in b])
    x0 = min(t[0] for t in b); y0 = min(t[1] for t in b)
    return x0, y0, cw, ch, len(b)
def sample(im, x0, y0, cw, ch, n=14):
    a = np.asarray(im.convert('RGB')).astype(float)
    M = np.zeros((n, n))
    for r in range(n):
        for c in range(n):
            cx, cy = int(x0+(c+.5)*cw), int(y0+(r+.5)*ch)
            p = a[cy-2:cy+3, cx-2:cx+3].reshape(-1,3).mean(0)
            M[r, c] = 1 - (0.299*p[0]+0.587*p[1]+0.114*p[2])/255
    return M
if __name__ == '__main__':
    ref, pb = sys.argv[1], sys.argv[2]
    reg = [int(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else None
    for nm, p, rg in (('ref', ref, None), ('pb', pb, reg)):
        im = Image.open(p); b = grey_boxes(im, rg); print(nm, im.size, 'grey squares', len(b))
        if b: print(' grid', grid_from_boxes(b))
