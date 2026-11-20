"""Step 4 (replicate design): composite kink maps and between-group statistics for
groups of replicate sets named <Group>_<n> (e.g. High_1..High_3, Low_1..Low_3).
Use this instead of kink_composites.py, whose panel layout is fixed to the T/D series.

Input:  <out>/kink_metrics.csv and <out>/lines/, filled/, filled_faded/   (from kink_maps.py)
        <work>/traces_circ.json                                          (for the sensitivity check)
Output: <out>/replicates/
          composite_kink_lines.png, composite_kink_filled.png, composite_kink_filled_faded.png
          kink_metrics_replicates.png    statistics (i)-(iii): replicates, mean +/- SD
                                         (kink and I(s) thresholds read from kink_metrics.csv)
          replicate_stats.csv            per group: values, mean, SD, length-weighted pooled value;
                                         between groups: ratio, Welch t, exact Mann-Whitney U, Hedges g
          sensitivity_edge_band.csv      (i)-(iii) recomputed with joints in the outer FOV band
                                         (r >= 0.95 R) ignored

        --exclude SET ...: the sets are left out of the statistics and composites and drawn as open
        circles in the chart; output folder <out>/replicates_excl_<SET>_<SET>/

Usage:  python kink_replicates.py [--out out_kink] [--work work] [--groups High Low] [--chain section|set]
                                  [--exclude High_2 Low_1]
"""
import os, re, ast, csv, json, argparse, math
import numpy as np
from PIL import Image
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--out', default='out_kink', help='output folder of kink_maps.py')
ap.add_argument('--work', default='work', help='folder holding traces_circ.json')
ap.add_argument('--groups', nargs='*', default=None, help='group order (default: as found, sorted)')
ap.add_argument('--chain', choices=('section', 'set'), default=None,
                help='chaining for the sensitivity check (default: the mode recorded in kink_metrics.csv)')
ap.add_argument('--exclude', nargs='*', default=[], help='sets left out of the statistics (shown as open circles)')
args = ap.parse_args()
OUT = args.out
RDIR = os.path.join(OUT, 'replicates' + ('_excl_' + '_'.join(sorted(args.exclude)) if args.exclude else ''))
os.makedirs(RDIR, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))

N = 1000; CEN, RPX = N / 2.0, 0.45 * N
KMAX, SAT, VAL, HUE_MAX = 60.0, 0.80, 1.00, 225.0
GROUP_COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']        # categorical slots 1-4
INK, INK2, INK3 = '#0b0b0b', '#52514e', '#8a8984'

def val_to_rgb(v):
    v = np.asarray(v, dtype=float)
    hue = HUE_MAX * (1.0 - np.clip(v, 0, KMAX) / KMAX) / 360.0
    return mcolors.hsv_to_rgb(np.stack([hue, np.full_like(hue, SAT), np.full_like(hue, VAL)], -1))
CMAP = mcolors.ListedColormap(val_to_rgb(np.linspace(0, KMAX, 256)))
NORM = mcolors.Normalize(vmin=0, vmax=KMAX)

# ------------------------------------------------------------------ metrics and groups
metrics = {}
with open(os.path.join(OUT, 'kink_metrics.csv'), newline='') as f:
    for row in csv.DictReader(f):
        metrics[row['set']] = row
for r in metrics.values():          # kinks per R from the counts, not the 2-decimal rounded column
    if 'n_kinks' in r and float(r['total_length_R']) > 0:
        r['kinks_per_R'] = repr(int(r['n_kinks']) / float(r['total_length_R']))
unknown = [x for x in args.exclude if x not in metrics]
if unknown:
    raise SystemExit('unknown set(s) to exclude: %s' % ', '.join(unknown))
groups, excluded = {}, {}
for name in metrics:
    m = re.match(r'^([A-Za-z]+)_(\w+)$', name)
    if m:
        (excluded if name in args.exclude else groups).setdefault(m.group(1), []).append(name)
order = args.groups if args.groups else sorted(groups)
order = [g for g in order if g in groups]
for g in order:
    groups[g].sort(key=lambda s: (len(s), s))
if not order:
    raise SystemExit('no <Group>_<n> sets in kink_metrics.csv')
ncol = max(len(groups[g]) for g in order)

first = next(iter(metrics.values()))
if 'kinks_per_R' not in first:
    raise SystemExit('kink_metrics.csv is from kink_maps.py version 1.0; run the current kink_maps.py first')
KDEG, IDEG, CHAIN_USED = float(first['kink_deg']), float(first['intensity_deg']), first['chain']
if any((r['kink_deg'], r['intensity_deg'], r['chain']) != (first['kink_deg'], first['intensity_deg'], first['chain'])
       for r in metrics.values()):
    raise SystemExit('kink_metrics.csv mixes thresholds or chaining modes')
if args.chain is None:
    args.chain = CHAIN_USED
KEYS = [('kinks_per_R', '(i) Kinks (θ ≥ %g°) per unit trace length' % KDEG, '1/R'),
        ('turning_deg_per_R', '(ii) Total turning per unit trace length', '°/R'),
        ('frac_len_I_gt', '(iii) Length fraction with I(s) > %g°' % IDEG, '')]

def pooled(names, key):
    L = np.array([float(metrics[n]['total_length_R']) for n in names])
    if key == 'kinks_per_R':
        return float(sum(int(metrics[n]['n_kinks']) for n in names) / L.sum())
    v = np.array([float(metrics[n][key]) for n in names])
    return float(np.sum(v * L) / L.sum())                       # length-weighted

def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    d = (np.mean(b) - np.mean(a)) / sp if sp > 0 else np.nan
    return float(d * (1 - 3 / (4 * (na + nb) - 9)))

# ------------------------------------------------------------------ composites
def composite(subdir, prefix, fname):
    fig = plt.figure(figsize=(3.3 * ncol + 1.6, 3.55 * len(order) + 0.7), dpi=150)
    gs = fig.add_gridspec(len(order), ncol, wspace=0.06, hspace=0.30,
                          left=0.12, right=0.86, top=1 - 0.55 / (3.55 * len(order) + 0.7), bottom=0.04)
    for ri, g in enumerate(order):
        for ci in range(ncol):
            ax = fig.add_subplot(gs[ri, ci])
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if ci < len(groups[g]):
                name = groups[g][ci]
                p = os.path.join(OUT, subdir, '%s_%s.png' % (prefix, name))
                ax.imshow(Image.open(p))
                m = metrics[name]
                ax.set_xlabel('kinks ≥%g°: %s   ·   turning %.0f°/R\nlength with I(s) > %g°: %d%%' % (
                    KDEG, m['n_kinks'], float(m['turning_deg_per_R']), IDEG,
                    round(100 * float(m['frac_len_I_gt']))), fontsize=8, color='0.25', labelpad=2)
                ax.text(0.5, 1.01, name, transform=ax.transAxes, ha='center', va='bottom', fontsize=9, color='0.35')
            else:
                ax.imshow(np.ones((N, N, 3)))
                ax.add_patch(matplotlib.patches.Circle((CEN, CEN), RPX, fill=False, ls='--', lw=1.2, ec='0.55'))
            if ci == 0:
                ax.text(-0.08, 0.5, g, transform=ax.transAxes, ha='right', va='center',
                        fontsize=13, fontweight='bold')
    cax = fig.add_axes([0.885, 0.22, 0.02, 0.56])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=CMAP, norm=NORM, ticks=[0, 15, 30, 45, 60])
    cb.set_label('Local turning intensity I(s) of contour trace (°)', fontsize=10)
    cb.ax.set_yticklabels(['0', '15', '30', '45', '≥60'])
    cb.ax.tick_params(labelsize=9)
    fig.text(0.895, 0.80, 'strong turning', ha='center', fontsize=8, color='0.4')
    fig.text(0.895, 0.19, 'straight', ha='center', fontsize=8, color='0.4')
    fig.savefig(os.path.join(RDIR, fname))
    plt.close(fig)

composite('lines', 'kink_lines_plain', 'composite_kink_lines.png')
composite('filled', 'kink_fill', 'composite_kink_filled.png')
composite('filled_faded', 'kink_fill_faded', 'composite_kink_filled_faded.png')

# ------------------------------------------------------------------ statistics table
rows = []
for key, title, unit in KEYS:
    vals = {g: np.array([float(metrics[n][key]) for n in groups[g]]) for g in order}
    for g in order:
        v = vals[g]
        rows.append(dict(metric=key, group=g, n=len(v), values=' / '.join('%g' % x for x in v),
                         mean=round(float(v.mean()), 4), sd=round(float(v.std(ddof=1)), 4) if len(v) > 1 else '',
                         pooled_length_weighted=round(pooled(groups[g], key), 4)))
    if len(order) == 2:
        a, b = vals[order[0]], vals[order[1]]
        t = stats.ttest_ind(b, a, equal_var=False)
        u = stats.mannwhitneyu(b, a, alternative='two-sided', method='exact')
        rows.append(dict(metric=key, group='%s vs %s' % (order[1], order[0]),
                         ratio_of_means=round(float(b.mean() / a.mean()), 3),
                         welch_t=round(float(t.statistic), 3), welch_df=round(float(t.df), 2),
                         welch_p=round(float(t.pvalue), 4), mannwhitney_U=float(u.statistic),
                         mannwhitney_exact_p=round(float(u.pvalue), 4), hedges_g=round(hedges_g(a, b), 3)))
fields = ['metric', 'group', 'n', 'values', 'mean', 'sd', 'pooled_length_weighted', 'ratio_of_means',
          'welch_t', 'welch_df', 'welch_p', 'mannwhitney_U', 'mannwhitney_exact_p', 'hedges_g']
with open(os.path.join(RDIR, 'replicate_stats.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ chart
fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.1), dpi=150, gridspec_kw=dict(wspace=0.42))
for ax, (key, title, unit) in zip(axes, KEYS):
    ymax = 0
    for gi, g in enumerate(order):
        names = groups[g]
        v = np.array([float(metrics[n][key]) for n in names])
        col = GROUP_COL[gi % len(GROUP_COL)]
        mu, sd = v.mean(), (v.std(ddof=1) if len(v) > 1 else 0)
        ax.bar(gi, mu, width=0.56, color=col, alpha=0.22, edgecolor='none', zorder=1)
        ax.errorbar(gi - 0.13, mu, yerr=sd, fmt='none', ecolor=INK2, elinewidth=1.4, capsize=6, capthick=1.4, zorder=3)
        ax.plot([gi - 0.28, gi - 0.01], [mu, mu], color=col, lw=2.2, solid_capstyle='butt', zorder=2)
        offs = np.linspace(0.02, 0.18, len(v)) if len(v) > 1 else [0.08]
        for k, (x, y) in enumerate(zip(gi + np.array(offs), v)):
            ax.plot(x, y, 'o', ms=8, mfc=col, mec='white', mew=1.6, zorder=4)
            ax.text(x + 0.05, y, ' %s' % names[k].split('_', 1)[1], fontsize=7.5, color=INK3,
                    va='center', ha='left', zorder=5)
        ymax = max(ymax, v.max(), mu + sd)
        for name in excluded.get(g, []):
            y = float(metrics[name][key])
            ax.plot(gi + 0.24, y, 'o', ms=8, mfc='white', mec=INK3, mew=1.4, zorder=4)
            ax.text(gi + 0.29, y, ' %s (excl.)' % name.split('_', 1)[1], fontsize=7.5, color=INK3,
                    va='center', ha='left', zorder=5)
            ymax = max(ymax, y)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=10, color=INK)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, ymax * 1.30)
    ax.set_title(title, fontsize=9.5, color=INK, loc='left')
    ax.set_ylabel(unit if unit else 'fraction', fontsize=9, color=INK2)
    ax.tick_params(axis='y', labelsize=8.5, colors=INK2)
    ax.tick_params(axis='x', length=0)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#c9c8c3')
    ax.grid(axis='y', color='#ecebe8', lw=0.8, zorder=0); ax.set_axisbelow(True)
    if len(order) == 2:
        r = [r_ for r_ in rows if r_['metric'] == key and 'vs' in str(r_['group'])][0]
        ax.text(0.5, 0.97, 'Welch t: p = %.3f   ·   Mann–Whitney: p = %.2f\n%s / %s = %.1f×' % (
            r['welch_p'], r['mannwhitney_exact_p'], order[1], order[0], r['ratio_of_means']),
            transform=ax.transAxes, ha='center', va='top', fontsize=7.8, color=INK2)
ns = [len(groups[g]) for g in order]
note = 'Bars: mean ± SD; dots: individual sets (n = %s). Two-sided tests' % ' / '.join(str(n) for n in ns)
if len(order) == 2:
    note += '; with n = %d / %d the smallest attainable exact Mann–Whitney p is %.2f' % (
        ns[0], ns[1], min(1.0, 2.0 / math.comb(ns[0] + ns[1], ns[0])))
note += '.'
if args.exclude:
    note += '\nOpen circles: excluded from mean, SD and tests (%s).' % ', '.join(sorted(args.exclude))
fig.text(0.5, -0.04 if not args.exclude else -0.07, note, ha='center', fontsize=8, color=INK3)
fig.suptitle('Kink statistics of CT contour traces — replicate sets', fontsize=11, color=INK, y=1.02)
fig.savefig(os.path.join(RDIR, 'kink_metrics_replicates.png'), bbox_inches='tight')
plt.close(fig)

# ------------------------------------------------------------------ sensitivity: outer FOV band excluded
# Traces end on the FOV boundary, and end points of different traces lying within CHAIN_TOL of
# each other there are chained by kink_maps.chain_polylines into one polyline with a reversal.
# Check: recompute (i)-(iii) with the turning angle of every joint at r >= EDGE_R set to 0
# (same polylines, same trace length; only joints in the outer band are ignored).
EDGE_R = 0.95
tr_path = os.path.join(args.work, 'traces_circ.json')
if os.path.exists(tr_path):
    src = open(os.path.join(HERE, 'kink_maps.py')).read()
    mod = ast.parse(src)
    keep = []
    for node in mod.body:
        if isinstance(node, ast.FunctionDef) and node.name in ('chain_polylines', 'polyline_kinks', 'sample_polyline'):
            keep.append(node)
        elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) for t in node.targets) and \
                node.targets[0].id in ('KMAX', 'CHAIN_TOL', 'MIN_SEG', 'SIGMA_S', 'DS'):
            keep.append(node)
    KM = {'np': np}
    exec(compile(ast.Module(body=keep, type_ignores=[]), 'kink_maps.py', 'exec'), KM)

    def analyse_band(S, rmax):
        if args.chain == 'set':
            polys = KM['chain_polylines']([x for rec in S['slides'] for x in rec['segs_circ']])
        else:
            polys = [P for rec in S['slides'] for P in KM['chain_polylines'](rec['segs_circ'])]
        n45, L, turn, I_all, n_band = 0, 0.0, 0.0, [], 0
        for P in polys:
            s_, sv, kappa = KM['polyline_kinks'](P)
            if kappa.size:
                r = np.hypot(P[1:-1, 0], P[1:-1, 1])
                n_band += int(np.sum((r >= rmax) & (kappa >= KDEG)))
                kappa = np.where(r >= rmax, 0.0, kappa)
            _, I = KM['sample_polyline'](P, s_, sv, kappa)
            I_all.append(I); L += s_[-1]; turn += kappa.sum(); n45 += int(np.sum(kappa >= KDEG))
        I_all = np.concatenate(I_all)
        return dict(n45=n45, n45_in_band=n_band, cusps=n45 / L, turning=turn / L, frac=float(np.mean(I_all > IDEG)))

    tr = json.load(open(tr_path))
    srows, band = [], {}
    for g in order:
        for name in groups[g]:
            st = analyse_band(tr[name], EDGE_R)
            band[name] = st
            m = metrics[name]
            srows.append(dict(set=name, group=g, n_kinks=m['n_kinks'],
                              kinks_in_outer_band=st['n45_in_band'],
                              kinks_per_R=m['kinks_per_R'], kinks_per_R_excl=round(st['cusps'], 2),
                              turning_deg_per_R=m['turning_deg_per_R'], turning_deg_per_R_excl=round(st['turning'], 1),
                              frac_len_I_gt=m['frac_len_I_gt'], frac_len_I_gt_excl=round(st['frac'], 3)))
    if len(order) == 2:
        for k2, lab in (('cusps', 'kinks_per_R_excl'), ('turning', 'turning_deg_per_R_excl'), ('frac', 'frac_len_I_gt_excl')):
            a = np.array([band[n][k2] for n in groups[order[0]]])
            b = np.array([band[n][k2] for n in groups[order[1]]])
            t = stats.ttest_ind(b, a, equal_var=False)
            u = stats.mannwhitneyu(b, a, alternative='two-sided', method='exact')
            srows.append({'set': '%s vs %s: %s' % (order[1], order[0], lab),
                          'group': 'mean %s %.3g±%.3g, %s %.3g±%.3g; ratio %.2f; Welch p=%.4f; MW exact p=%.3f' % (
                              order[0], a.mean(), a.std(ddof=1), order[1], b.mean(), b.std(ddof=1),
                              b.mean() / a.mean(), t.pvalue, u.pvalue)})
    with open(os.path.join(RDIR, 'sensitivity_edge_band.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0].keys())); w.writeheader(); w.writerows(srows)

for r_ in rows:
    print({k: v for k, v in r_.items() if v not in ('', None)})
print('replicate composites, chart and tables in', RDIR)
