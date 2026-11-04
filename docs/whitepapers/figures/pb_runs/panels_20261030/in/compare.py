import sys, numpy as np
sys.path.insert(0, 'in')
from cells import sample
from PIL import Image
def stats(ref, pb, rbox, pbox, n=14, I=3, J=4):
    ra = sample(Image.open(ref), rbox[0], rbox[1], rbox[2], rbox[3], n)
    pa = sample(Image.open(pb), pbox[0], pbox[1], pbox[2], pbox[3], n)
    keep = np.ones((n, n), bool); keep[I, :] = False; keep[:, J] = False; np.fill_diagonal(keep, False)
    a, b = ra[keep], pa[keep]
    # rescale to 0-1
    a = (a-a.min())/(a.max()-a.min()); b = (b-b.min())/(b.max()-b.min())
    r = np.corrcoef(a, b)[0, 1]
    # cells in the ref's top-quartile "dark" set that are also dark in PB
    qa = a >= np.quantile(a, .75); qb = b >= np.quantile(b, .75)
    return r, (qa & qb).sum()/qa.sum(), np.abs(a-b).mean()
if __name__ == '__main__':
    run = lambda r: 'runs/run_20261030_164706_%s/final_output.png' % r
    T = [('P1 D', 'ref/panel1.png', (180, 294, 49.5, 49.5), run(sys.argv[1]), tuple(map(float, sys.argv[2].split(',')))),
         ('P2 S', 'ref/panel2.png', (588, 522, 30, 30), run(sys.argv[3]), tuple(map(float, sys.argv[4].split(',')))),
         ('P3 S', 'ref/panel3.png', (156, 309, 24, 24), run(sys.argv[5]), tuple(map(float, sys.argv[6].split(',')))),
         ('P3 A', 'ref/panel3.png', (672, 309, 24, 24), run(sys.argv[5]), tuple(map(float, sys.argv[7].split(',')))),
         ('P4 t1', 'ref/panel4.png', (63, 309, 21.6, 21.6), run(sys.argv[8]), tuple(map(float, sys.argv[9].split(',')))),
         ('P4 t3', 'ref/panel4.png', (399, 309, 21.6, 21.6), run(sys.argv[8]), tuple(map(float, sys.argv[10].split(',')))),
         ('P4 t9', 'ref/panel4.png', (735, 309, 21.6, 21.6), run(sys.argv[8]), tuple(map(float, sys.argv[11].split(','))))]
    for nm, ref, rb, pb, pbx in T:
        x0, y0, x1, y1 = pbx
        r, ov, mad = stats(ref, pb, rb, (x0, y0, (x1-x0)/14, (y1-y0)/14))
        print(f'{nm}: cell-wise Pearson r = {r:.2f}; dark-quartile overlap = {ov:.2f}; mean |diff| = {mad:.2f}')
