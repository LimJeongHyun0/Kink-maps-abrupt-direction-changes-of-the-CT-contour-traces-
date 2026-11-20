"""Effect of chaining across sections: for every set, compute the kink statistics with the
elements chained within each section (kink_maps.py --chain section) and across all sections of
the set (--chain set, version 1.0), and list the joints that exist only in the latter, i.e.
joints created by joining end points of different sections.

Input:  <work>/traces_circ.json (one slide per section)
Output: <out>/chain_mode_comparison.csv   per set: statistics (i)-(iii) in both modes, number of
                                          extra joints and extra kinks (>= --kink-deg)
        <out>/cross_section_joints.csv    position (unit disc), radius and turning angle of each
                                          extra joint

Usage:  python compare_chain_modes.py [--work work] [--out out_kink] [--kink-deg 20] [--intensity-deg 20]
"""
import os, ast, csv, json, argparse
import numpy as np

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--work', default='work')
ap.add_argument('--out', default='out_kink')
ap.add_argument('--kink-deg', type=float, default=20.0)
ap.add_argument('--intensity-deg', type=float, default=20.0)
args = ap.parse_args()
KD, ID = args.kink_deg, args.intensity_deg
os.makedirs(args.out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))

mod = ast.parse(open(os.path.join(HERE, 'kink_maps.py')).read())
keep = [n for n in mod.body if (isinstance(n, ast.FunctionDef) and n.name in
        ('chain_polylines', 'polyline_kinks', 'sample_polyline', 'analyse')) or
        (isinstance(n, ast.Assign) and all(isinstance(t, ast.Name) for t in n.targets) and
         n.targets[0].id in ('KMAX', 'CHAIN_TOL', 'MIN_SEG', 'SIGMA_S', 'DS'))]
KM = {'np': np}
exec(compile(ast.Module(body=keep, type_ignores=[]), 'kink_maps.py', 'exec'), KM)

TOL, MIN_SEG = KM['CHAIN_TOL'], KM['MIN_SEG']

def chain_with_source(segs, srcs, tol=TOL):
    """kink_maps.chain_polylines (same greedy order and tests) that also returns, for every
    edge of every polyline, the section of the element that created it"""
    segs = [np.array([[s[0], s[1]], [s[2], s[3]]]) for s in segs]
    used = [False] * len(segs)
    out = []
    for i in range(len(segs)):
        if used[i]:
            continue
        used[i] = True
        pts = [segs[i][0], segs[i][1]]; esrc = [srcs[i]]
        grown = True
        while grown:
            grown = False
            for j in range(len(segs)):
                if used[j]:
                    continue
                a, b = segs[j]
                if np.hypot(*(pts[-1] - a)) < tol:
                    pts.append(b); esrc.append(srcs[j]); used[j] = True; grown = True
                elif np.hypot(*(pts[-1] - b)) < tol:
                    pts.append(a); esrc.append(srcs[j]); used[j] = True; grown = True
                elif np.hypot(*(pts[0] - b)) < tol:
                    pts.insert(0, a); esrc.insert(0, srcs[j]); used[j] = True; grown = True
                elif np.hypot(*(pts[0] - a)) < tol:
                    pts.insert(0, b); esrc.insert(0, srcs[j]); used[j] = True; grown = True
        P = np.array(pts)
        keep = [0]
        for k in range(1, len(P)):
            if np.hypot(*(P[k] - P[keep[-1]])) >= MIN_SEG or k == len(P) - 1:
                keep.append(k)
        # edge e of the merged polyline spans original edges keep[e] .. keep[e+1]-1
        src_in = [esrc[keep[e + 1] - 1] for e in range(len(keep) - 1)]     # last original edge
        src_out = [esrc[keep[e]] for e in range(len(keep) - 1)]            # first original edge
        P = P[keep]
        if len(P) >= 2 and np.hypot(*(P[-1] - P[-2])) < 1e-9:
            P = P[:-1]; src_in = src_in[:-1]; src_out = src_out[:-1]
        out.append((P, src_in, src_out))
    return out

tr = json.load(open(os.path.join(args.work, 'traces_circ.json')))
rows, joints = [], []
for name in sorted(tr):
    S = tr[name]
    _, _, _, _, st_sec = KM['analyse'](S, 'section', KD, ID)
    _, _, _, _, st_set = KM['analyse'](S, 'set', KD, ID)
    segs = [x for rec in S['slides'] for x in rec['segs_circ']]
    srcs = [k for k, rec in enumerate(S['slides']) for _ in rec['segs_circ']]
    extra = []
    for P, s_in, s_out in chain_with_source(segs, srcs):
        _, _, kappa = KM['polyline_kinks'](P)
        for v in range(1, len(P) - 1):           # joint v: incoming edge v-1, outgoing edge v
            if s_in[v - 1] != s_out[v]:
                x, y, kk = float(P[v, 0]), float(P[v, 1]), float(kappa[v - 1])
                extra.append((x, y, kk))
                joints.append(dict(set=name, section_in=s_in[v - 1] + 1, section_out=s_out[v] + 1,
                                   X=round(x, 4), Y=round(y, 4), r=round(float(np.hypot(x, y)), 4),
                                   turning_deg=round(kk, 1)))
    extra = np.array(extra).reshape(-1, 3)
    rows.append(dict(set=name, sections=len(S['slides']),
                     joints_section=st_sec['n_joints'], joints_set=st_set['n_joints'],
                     cross_section_joints=len(extra), cross_section_kinks=int(np.sum(extra[:, 2] >= KD)),
                     cross_section_joints_r_ge095=int(np.sum(np.hypot(extra[:, 0], extra[:, 1]) >= 0.95)),
                     kinks_per_R_section=round(st_sec['kinks_per_length'], 3),
                     kinks_per_R_set=round(st_set['kinks_per_length'], 3),
                     turning_deg_per_R_section=round(st_sec['turning_per_length'], 1),
                     turning_deg_per_R_set=round(st_set['turning_per_length'], 1),
                     frac_len_I_gt_section=round(st_sec['frac_len_I_gt'], 3),
                     frac_len_I_gt_set=round(st_set['frac_len_I_gt'], 3),
                     kink_deg=KD, intensity_deg=ID))
    r = rows[-1]
    print('%-8s cross-section joints %2d (kinks >=%g deg: %d, at r>=0.95: %d) | kinks/R %.2f -> %.2f | '
          'turning %.0f -> %.0f deg/R | frac(I>%g) %.3f -> %.3f   (section -> set)' % (
              name, r['cross_section_joints'], KD, r['cross_section_kinks'], r['cross_section_joints_r_ge095'],
              r['kinks_per_R_section'], r['kinks_per_R_set'], r['turning_deg_per_R_section'],
              r['turning_deg_per_R_set'], ID, r['frac_len_I_gt_section'], r['frac_len_I_gt_set']))
for fn, data in (('chain_mode_comparison.csv', rows), ('cross_section_joints.csv', joints)):
    with open(os.path.join(args.out, fn), 'w', newline='') as f:
        if data:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys())); w.writeheader(); w.writerows(data)
print('saved', os.path.join(args.out, 'chain_mode_comparison.csv'))
