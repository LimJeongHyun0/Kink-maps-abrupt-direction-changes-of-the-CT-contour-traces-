"""Sensitivity of the kink statistics (i)-(iii) to their criteria, and the turning-angle
distribution used to choose the kink threshold.

The polylines, vertices and turning angles theta_k are those of kink_maps.py (same functions;
chaining within each section by default); only the criterion of each statistic is varied:

  (i)   kinks per unit trace length            kink threshold     theta_kink in KINK_T   (adopted 20 deg; v1.0 45 deg)
  (ii)  total turning per unit trace length    smallest angle     theta_min in TMIN      (adopted 0 deg)
        summed (vertices below it are left out as tracing noise)
  (iii) length fraction with I(s) > theta_I    threshold          theta_I in INT_T       (adopted 20 deg; v1.0 30 deg)
        and Gaussian spread                    sigma in SIGMAS                           (0.06 R)

For every criterion value: per-specimen values, group mean +/- SD, ratio of the means, Welch t,
exact two-sided Mann-Whitney U, Hedges g, and whether every specimen of the second group lies
above every specimen of the first (complete separation).

Output (<out>/sensitivity/):
  sensitivity_values.csv / sensitivity_summary.csv
  sensitivity_grid.png         bar chart of every criterion value; adopted value framed in black, v1.0 in grey
  sensitivity_trend.png        value, ratio and p against the criterion for (i)-(iii)
  sensitivity_iii_sigma.png    (iii) against theta_I for every sigma
  charts/<metric>_<value>.png  one bar chart per criterion value
  angle_distribution.png       exceedance per specimen and turning-angle distribution per unit trace length
  angle_bins.csv               vertices per unit trace length in 5-deg bins, per group (pooled) and specimen
  band_rates.csv               vertices per unit trace length in the bands 0-10, 10-20, >=20 deg per specimen,
                               with the group comparison of each band
  threshold_justification.png  (a) distribution per unit length with the three bands, (b) ratio per 5-deg bin,
                               (c) specimen values per band
  joint_density.csv            per specimen: vertices per unit trace length, median and mean angle

Usage:  python kink_sensitivity.py [--work work] [--out out_kink] [--groups High Low] [--chain section|set]
                                   [--labels "High=D3 (High)" "Low=D1 (Low)"]
"""
import os, re, ast, csv, json, math, argparse
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--work', default='work')
ap.add_argument('--out', default='out_kink')
ap.add_argument('--groups', nargs='*', default=None)
ap.add_argument('--chain', choices=('section', 'set'), default='section')
ap.add_argument('--labels', nargs='*', default=[], help='display names, e.g. "High=D3 (High)"')
args = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
SDIR = os.path.join(args.out, 'sensitivity')
os.makedirs(os.path.join(SDIR, 'charts'), exist_ok=True)

KINK_T = [10, 15, 20, 30, 45, 60, 90]
TMIN = [0, 5, 10, 15, 20, 30]
INT_T = [10, 15, 20, 30, 45, 60]
SIGMAS = [0.03, 0.06, 0.12]
ADOPTED = {'i': 20, 'ii': 0, 'iii': 20}
V10 = {'i': 45, 'ii': None, 'iii': 30}
BANDS = [(0, 10), (10, 20), (20, 181)]
SIGMA0 = 0.06
GROUP_COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
INK, INK2, INK3 = '#0b0b0b', '#52514e', '#8a8984'
LAB = dict(l.split('=', 1) for l in args.labels)

# ---- kink_maps.py functions and constants
mod = ast.parse(open(os.path.join(HERE, 'kink_maps.py')).read())
keep = [n for n in mod.body if (isinstance(n, ast.FunctionDef) and n.name in ('chain_polylines', 'polyline_kinks')) or
        (isinstance(n, ast.Assign) and all(isinstance(t, ast.Name) for t in n.targets) and
         n.targets[0].id in ('CHAIN_TOL', 'MIN_SEG', 'SIGMA_S', 'DS'))]
KM = {'np': np}
exec(compile(ast.Module(body=keep, type_ignores=[]), 'kink_maps.py', 'exec'), KM)
DS = KM['DS']

def intensity(s, sv, kappa, sigma):
    """kink_maps.sample_polyline with the spread sigma as a parameter"""
    total = s[-1]
    if total < 1e-9:
        return np.zeros(0)
    ss = np.r_[np.arange(0, total, DS), total]
    I = np.zeros_like(ss)
    for sk, kk in zip(sv, kappa):
        I += kk * np.exp(-0.5 * ((ss - sk) / sigma) ** 2)
    return I

tr = json.load(open(os.path.join(args.work, 'traces_circ.json')))
groups = {}
for name in tr:
    m = re.match(r'^([A-Za-z]+)_(\w+)$', name)
    if m:
        groups.setdefault(m.group(1), []).append(name)
order = [g for g in (args.groups or sorted(groups)) if g in groups]
for g in order:
    groups[g].sort(key=lambda s: (len(s), s))
lab = lambda g: LAB.get(g, g)

spec = {}
for g in order:
    for name in groups[g]:
        S = tr[name]
        if args.chain == 'set':
            polys = KM['chain_polylines']([x for rec in S['slides'] for x in rec['segs_circ']])
        else:
            polys = [P for rec in S['slides'] for P in KM['chain_polylines'](rec['segs_circ'])]
        kap, L, Is = [], 0.0, {sg: [] for sg in SIGMAS}
        for P in polys:
            s, sv, kappa = KM['polyline_kinks'](P)
            L += s[-1]; kap.append(kappa)
            for sg in SIGMAS:
                Is[sg].append(intensity(s, sv, kappa, sg))
        spec[name] = dict(kappa=np.concatenate(kap) if kap else np.zeros(0), L=L,
                          I={sg: np.concatenate(v) for sg, v in Is.items()})

def value(name, metric, p, sigma=SIGMA0):
    d = spec[name]
    if metric == 'i':
        return float(np.sum(d['kappa'] >= p) / d['L'])
    if metric == 'ii':
        k = d['kappa']
        return float(np.sum(k[k >= p]) / d['L'])
    return float(np.mean(d['I'][sigma] > p))

# ---- check against kink_metrics.csv when present
mpath = os.path.join(args.out, 'kink_metrics.csv')
if os.path.exists(mpath):
    with open(mpath) as f:
        km = {r['set']: r for r in csv.DictReader(f)}
    if km and 'kinks_per_R' in next(iter(km.values())):
        worst = 0.0
        for name in spec:
            if name in km:
                r = km[name]
                kd, idg = float(r['kink_deg']), float(r['intensity_deg'])
                worst = max(worst, abs(value(name, 'i', kd) - float(r['kinks_per_R'])) / 0.005,
                            abs(value(name, 'ii', 0) - float(r['turning_deg_per_R'])) / 0.05,
                            abs(value(name, 'iii', idg) - float(r['frac_len_I_gt'])) / 0.0005)
        print('criteria of kink_metrics.csv reproduced: %s' % ('yes, within rounding' if worst <= 1.0001 else
                                                              'NO (%.2f rounding units)' % worst))

def compare(a, b):
    out = dict(ratio=b.mean() / a.mean() if a.mean() > 0 else np.inf,
               welch_p=float(stats.ttest_ind(b, a, equal_var=False).pvalue) if (a.std() > 0 or b.std() > 0) else np.nan,
               mw_p=float(stats.mannwhitneyu(b, a, alternative='two-sided', method='exact').pvalue),
               separated=bool(b.min() > a.max()))
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    out['g'] = (b.mean() - a.mean()) / sp * (1 - 3 / (4 * (len(a) + len(b)) - 9)) if sp > 0 else np.nan
    return out

def group_stats(metric, p, sigma=SIGMA0):
    v = {g: np.array([value(n, metric, p, sigma) for n in groups[g]]) for g in order}
    out = dict(values=v)
    if len(order) == 2:
        out.update(compare(v[order[0]], v[order[1]]))
    return out

# metric: (title, mathtext symbol, ascii symbol, unit, y label)
LABEL = {'i': ('(i) Kinks per unit trace length', r'$\theta_{\mathrm{kink}}$', 'theta_kink', '1/R',
               r'kinks ($\theta_k \geq \theta_{\mathrm{kink}}$) per R'),
         'ii': ('(ii) Total turning per unit trace length', r'$\theta_{\min}$', 'theta_min', '°/R',
                r'$\Sigma\theta_k$ ($\theta_k \geq \theta_{\min}$) per R'),
         'iii': ('(iii) Length fraction with $I(s) > \\theta_I$', r'$\theta_I$', 'theta_I', '',
                 'fraction of trace length')}
GRID = {'i': KINK_T, 'ii': TMIN, 'iii': INT_T}

def tag(metric, p):
    return ' (adopted)' if p == ADOPTED[metric] else (' (v1.0)' if p == V10[metric] else '')

# ---- tables
vrows, srows = [], []
for metric, plist in GRID.items():
    for sg in (SIGMAS if metric == 'iii' else [None]):
        for p in plist:
            st = group_stats(metric, p, sg if sg else SIGMA0)
            crit = '%s=%g deg' % (LABEL[metric][2], p) + (', sigma=%gR' % sg if sg else '')
            adopted = int(p == ADOPTED[metric] and sg in (None, SIGMA0))
            for g in order:
                for n, x in zip(groups[g], st['values'][g]):
                    vrows.append(dict(metric=metric, criterion=crit, adopted=adopted, group=g, specimen=n,
                                      value=round(x, 5)))
            r = dict(metric=metric, criterion=crit, adopted=adopted)
            for g in order:
                vv = st['values'][g]
                r['%s_values' % g] = ' / '.join('%.4g' % x for x in vv)
                r['%s_mean' % g] = round(float(vv.mean()), 5)
                r['%s_sd' % g] = round(float(vv.std(ddof=1)), 5)
            if len(order) == 2:
                r.update(ratio=round(st['ratio'], 3) if np.isfinite(st['ratio']) else 'inf',
                         welch_p=round(st['welch_p'], 4), mannwhitney_exact_p=round(st['mw_p'], 4),
                         hedges_g=round(st['g'], 3), complete_separation=int(st['separated']))
            srows.append(r)
for fn, rows in (('sensitivity_values.csv', vrows), ('sensitivity_summary.csv', srows)):
    with open(os.path.join(SDIR, fn), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

def style(ax, small=False):
    ax.spines[['top', 'right']].set_visible(False); ax.spines[['left', 'bottom']].set_color('#c9c8c3')
    ax.tick_params(labelsize=7 if small else 8.5, colors=INK2)
    ax.grid(axis='y', color='#ecebe8', lw=0.8, zorder=0); ax.set_axisbelow(True)

# ---- bar chart of one criterion value
def bar_panel(ax, metric, p, sigma=SIGMA0, title=None, small=False):
    st = group_stats(metric, p, sigma)
    ymax = 0
    for gi, g in enumerate(order):
        v = st['values'][g]; col = GROUP_COL[gi]
        mu, sd = v.mean(), v.std(ddof=1)
        ax.bar(gi, mu, width=0.56, color=col, alpha=0.22, edgecolor='none', zorder=1)
        ax.errorbar(gi - 0.13, mu, yerr=sd, fmt='none', ecolor=INK2, elinewidth=1.2, capsize=4, capthick=1.2, zorder=3)
        ax.plot([gi - 0.28, gi - 0.01], [mu, mu], color=col, lw=2.0, solid_capstyle='butt', zorder=2)
        offs = np.linspace(0.02, 0.18, len(v))
        for k, (x, y) in enumerate(zip(gi + offs, v)):
            ax.plot(x, y, 'o', ms=5.5 if small else 8, mfc=col, mec='white', mew=1.2, zorder=4)
            if not small:
                ax.text(x + 0.05, y, ' %s' % groups[g][k].split('_', 1)[1], fontsize=7.5, color=INK3, va='center')
        ymax = max(ymax, v.max(), mu + sd)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([lab(g) for g in order], fontsize=7.5 if small else 10, color=INK)
    ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(0, (ymax if ymax > 0 else 1) * 1.38)
    ax.tick_params(axis='x', length=0)
    style(ax, small)
    if len(order) == 2:
        rt = '%.1f×' % st['ratio'] if np.isfinite(st['ratio']) else '∞'
        ax.text(0.5, 0.98, 'Welch p = %.3f · MW p = %.2f\n%s/%s = %s' % (st['welch_p'], st['mw_p'], order[1], order[0], rt),
                transform=ax.transAxes, ha='center', va='top', fontsize=6.6 if small else 7.8, color=INK2)
    if title:
        ax.set_title(title, fontsize=8.5 if small else 9.5, color=INK, loc='left',
                     fontweight='bold' if (p == ADOPTED[metric] and sigma == SIGMA0) else 'normal')

for metric, plist in GRID.items():
    for p in plist:
        fig, ax = plt.subplots(figsize=(3.6, 3.6), dpi=150)
        bar_panel(ax, metric, p, title='%s\n%s = %g°%s' % (LABEL[metric][0], LABEL[metric][1], p, tag(metric, p)))
        ax.set_ylabel(LABEL[metric][4] + (' (%s)' % LABEL[metric][3] if LABEL[metric][3] else ''), fontsize=8.5, color=INK2)
        fig.savefig(os.path.join(SDIR, 'charts', '%s_%s%02d.png' % (metric, LABEL[metric][2], p)), bbox_inches='tight')
        plt.close(fig)

ncol = max(len(v) for v in GRID.values())
fig, axes = plt.subplots(3, ncol, figsize=(2.25 * ncol + 0.6, 8.8), dpi=150, gridspec_kw=dict(wspace=0.42, hspace=0.62))
for ri, (metric, plist) in enumerate(GRID.items()):
    for ci in range(ncol):
        ax = axes[ri, ci]
        if ci >= len(plist):
            ax.axis('off'); continue
        p = plist[ci]
        bar_panel(ax, metric, p, title='%s = %g°%s' % (LABEL[metric][1], p, tag(metric, p)), small=True)
        if p == ADOPTED[metric] or p == V10[metric]:
            for s_ in ax.spines.values():
                s_.set_visible(True); s_.set_linewidth(1.3)
                s_.set_color(INK if p == ADOPTED[metric] else '#b5b4af')
        if ci == 0:
            ax.set_ylabel(LABEL[metric][4] + (' (%s)' % LABEL[metric][3] if LABEL[metric][3] else ''), fontsize=7.5, color=INK2)
    fig.text(0.01, axes[ri, 0].get_position().y1 + 0.035, LABEL[metric][0] + ('' if metric != 'iii' else '  (σ = 0.06 R)'),
             fontsize=10.5, color=INK, fontweight='bold')
fig.text(0.5, 0.005, 'Bars: mean ± SD; dots: specimens (n = %s). Black frame: adopted criterion; grey frame: version 1.0. '
         'Two-sided tests; smallest attainable exact Mann–Whitney p at this n: %.2f.' % (
             ' / '.join(str(len(groups[g])) for g in order),
             2.0 / math.comb(sum(len(groups[g]) for g in order), len(groups[order[0]]))),
         ha='center', fontsize=8, color=INK3)
fig.savefig(os.path.join(SDIR, 'sensitivity_grid.png'), bbox_inches='tight')
plt.close(fig)

# ---- trend figure
def trend_column(axv, axr, axp, metric, plist, sigma=SIGMA0):
    P = np.array(plist, float)
    res = [group_stats(metric, p, sigma) for p in plist]
    for gi, g in enumerate(order):
        col = GROUP_COL[gi]
        V = np.array([r['values'][g] for r in res])
        for k in range(V.shape[1]):
            axv.plot(P, V[:, k], color=col, lw=0.9, alpha=0.45, zorder=2)
            axv.text(P[-1] + (P[-1] - P[0]) * 0.02, V[-1, k], groups[g][k].split('_', 1)[1], fontsize=6.5, color=col, va='center')
        axv.errorbar(P + (gi - 0.5) * (P[-1] - P[0]) * 0.012, V.mean(1), yerr=V.std(1, ddof=1), color=col, lw=2.0,
                     marker='o', ms=5, capsize=3, zorder=3, label=lab(g))
    for ax in (axv, axr, axp):
        ax.axvline(ADOPTED[metric], color=INK, lw=1.0, ls='--', zorder=1)
        if V10[metric] is not None:
            ax.axvline(V10[metric], color=INK3, lw=0.9, ls=':', zorder=1)
    if len(order) == 2:
        ratio = np.array([r['ratio'] if np.isfinite(r['ratio']) else np.nan for r in res])
        axr.plot(P, ratio, color=INK, lw=1.8, marker='o', ms=4.5)
        for x, r in zip(P, res):
            if not np.isfinite(r['ratio']):
                axr.annotate('∞', (x, np.nanmax(ratio)), ha='center', fontsize=9, color=INK2)
        axr.axhline(1, color=INK3, lw=0.8)
        wp = np.array([r['welch_p'] for r in res]); mp = np.array([r['mw_p'] for r in res])
        axp.plot(P, wp, color=INK, lw=1.8, marker='o', ms=4.5, label='Welch t')
        axp.plot(P, mp, color=INK3, lw=1.2, marker='s', ms=4, ls=':', label='Mann–Whitney (exact)')
        axp.axhline(0.05, color='#c9c8c3', lw=0.9)
        axp.set_ylim(0, max(0.45, np.nanmax(np.r_[wp, mp]) * 1.1))
        for x, r in zip(P, res):
            if r['separated']:
                axp.plot(x, -0.02 * axp.get_ylim()[1], marker='^', color=INK, ms=4, clip_on=False)
    for ax in (axv, axr, axp):
        style(ax); ax.set_xticks(P)

fig, axes = plt.subplots(3, 3, figsize=(13.5, 9.6), dpi=150, gridspec_kw=dict(hspace=0.42, wspace=0.28, height_ratios=[1.5, 0.8, 0.8]))
for ci, (metric, plist) in enumerate(GRID.items()):
    trend_column(axes[0, ci], axes[1, ci], axes[2, ci], metric, plist)
    axes[0, ci].set_title(LABEL[metric][0] + ('  (σ = 0.06 R)' if metric == 'iii' else ''), fontsize=10, color=INK, loc='left')
    axes[0, ci].set_ylabel(LABEL[metric][4] + (' (%s)' % LABEL[metric][3] if LABEL[metric][3] else ''), fontsize=8.5, color=INK2)
    axes[2, ci].set_xlabel('%s (°)' % LABEL[metric][1], fontsize=9, color=INK2)
    axes[0, ci].set_ylim(bottom=0)
axes[1, 0].set_ylabel('%s / %s (ratio of means)' % (order[1], order[0]), fontsize=8.5, color=INK2)
axes[2, 0].set_ylabel('p (two-sided)', fontsize=8.5, color=INK2)
axes[0, 0].legend(fontsize=8, frameon=False, loc='upper right')
axes[2, 0].legend(fontsize=7.5, frameon=False, loc='upper left')
fig.text(0.5, 0.02, 'Top: group mean ± SD (thick) and specimens (thin, numbered). Dashed: adopted criterion; dotted: version 1.0. '
         '▲: every %s specimen above every %s specimen; grey line: p = 0.05.' % (order[1], order[0]),
         ha='center', fontsize=8, color=INK3)
fig.savefig(os.path.join(SDIR, 'sensitivity_trend.png'), bbox_inches='tight')
plt.close(fig)

fig, axes = plt.subplots(3, len(SIGMAS), figsize=(4.4 * len(SIGMAS), 8.6), dpi=150,
                         gridspec_kw=dict(hspace=0.42, wspace=0.28, height_ratios=[1.5, 0.8, 0.8]))
for ci, sg in enumerate(SIGMAS):
    trend_column(axes[0, ci], axes[1, ci], axes[2, ci], 'iii', INT_T, sigma=sg)
    axes[0, ci].set_title('(iii) σ = %g R%s' % (sg, '  (adopted)' if sg == SIGMA0 else ''), fontsize=10, color=INK, loc='left')
    axes[0, ci].set_ylim(0, 1)
    axes[2, ci].set_xlabel(r'$\theta_I$ (°)', fontsize=9, color=INK2)
axes[0, 0].set_ylabel(r'fraction of trace length with $I(s) > \theta_I$', fontsize=8.5, color=INK2)
axes[1, 0].set_ylabel('%s / %s' % (order[1], order[0]), fontsize=8.5, color=INK2)
axes[2, 0].set_ylabel('p (two-sided)', fontsize=8.5, color=INK2)
axes[0, 0].legend(fontsize=8, frameon=False, loc='upper right'); axes[2, 0].legend(fontsize=7.5, frameon=False, loc='upper left')
fig.savefig(os.path.join(SDIR, 'sensitivity_iii_sigma.png'), bbox_inches='tight')
plt.close(fig)

# ---- turning-angle distribution per unit trace length
bins = np.arange(0, 185, 5)
Lg = {g: sum(spec[n]['L'] for n in groups[g]) for g in order}
Kg = {g: np.concatenate([spec[n]['kappa'] for n in groups[g]]) for g in order}
Hg = {g: np.histogram(Kg[g], bins)[0] / Lg[g] for g in order}
with open(os.path.join(SDIR, 'angle_bins.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    names = [n for g in order for n in groups[g]]
    w.writerow(['bin_deg'] + ['%s_pooled_per_R' % g for g in order] +
               (['%s_over_%s' % (order[1], order[0])] if len(order) == 2 else []) + ['%s_per_R' % n for n in names])
    Hs = {n: np.histogram(spec[n]['kappa'], bins)[0] / spec[n]['L'] for n in names}
    for i in range(len(bins) - 1):
        if all(Hg[g][i] == 0 for g in order):
            continue
        row = ['%d-%d' % (bins[i], bins[i + 1])] + ['%.3f' % Hg[g][i] for g in order]
        if len(order) == 2:
            a, b = Hg[order[0]][i], Hg[order[1]][i]
            row.append('%.2f' % (b / a) if a > 0 else 'inf')
        row += ['%.3f' % Hs[n][i] for n in names]
        w.writerow(row)

band_rows, band_vals = [], {}
for lo, hi in BANDS:
    key = '%d-%s' % (lo, '' if hi > 180 else hi)
    v = {g: np.array([np.sum((spec[n]['kappa'] >= lo) & (spec[n]['kappa'] < hi)) / spec[n]['L'] for n in groups[g]]) for g in order}
    band_vals[(lo, hi)] = v
    r = dict(band_deg=('>=%d' % lo) if hi > 180 else '%d-%d' % (lo, hi))
    for g in order:
        r['%s_values' % g] = ' / '.join('%.2f' % x for x in v[g])
        r['%s_mean' % g] = round(float(v[g].mean()), 3); r['%s_sd' % g] = round(float(v[g].std(ddof=1)), 3)
    if len(order) == 2:
        c = compare(v[order[0]], v[order[1]])
        r.update(ratio=round(c['ratio'], 3), welch_p=round(c['welch_p'], 4), mannwhitney_exact_p=round(c['mw_p'], 4),
                 complete_separation=int(c['separated']))
    band_rows.append(r)
with open(os.path.join(SDIR, 'band_rates.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(band_rows[0].keys())); w.writeheader(); w.writerows(band_rows)

with open(os.path.join(SDIR, 'joint_density.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['specimen', 'group', 'trace_length_R', 'vertices', 'vertices_per_R', 'median_theta_deg', 'mean_theta_deg'])
    for g in order:
        for n in groups[g]:
            k = spec[n]['kappa']
            w.writerow([n, g, round(spec[n]['L'], 3), k.size, round(k.size / spec[n]['L'], 2),
                        round(float(np.median(k)), 1), round(float(k.mean()), 1)])
    allk = np.concatenate(list(Kg.values()))
    w.writerow(['all', 'pooled', round(sum(Lg.values()), 3), allk.size, round(allk.size / sum(Lg.values()), 2),
                round(float(np.median(allk)), 1), round(float(allk.mean()), 1)])

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), dpi=150, gridspec_kw=dict(wspace=0.25))
th = np.arange(0, 181, 1.0)
for gi, g in enumerate(order):
    col = GROUP_COL[gi]
    for k, n in enumerate(groups[g]):
        ex = np.array([np.sum(spec[n]['kappa'] >= t) for t in th]) / spec[n]['L']
        axes[0].plot(th, ex, color=col, lw=1.2, alpha=0.85, label=lab(g) if k == 0 else None)
        last = int(np.max(np.nonzero(ex > 0)[0])) if np.any(ex > 0) else 0
        axes[0].text(th[last] + 1.5, ex[last], n.split('_', 1)[1], fontsize=7, color=col, va='center')
    axes[1].hist(Kg[g], bins=bins, weights=np.full(Kg[g].size, 1.0 / Lg[g]), color=col, alpha=0.45,
                 label='%s (%d vertices, %.1f R of trace)' % (lab(g), Kg[g].size, Lg[g]), edgecolor='white', linewidth=0.5)
axes[0].set_yscale('log'); axes[0].set_ylim(0.02, None); axes[0].set_xlim(-3, 140); axes[1].set_xlim(-3, 140)
axes[0].set_xlabel(r'threshold $\theta$ (°)', fontsize=9, color=INK2)
axes[0].set_ylabel(r'vertices with $\theta_k \geq \theta$ per unit trace length (1/R)', fontsize=9, color=INK2)
axes[0].set_title(r'Exceedance per specimen ($\theta$ = 0: all vertices)', fontsize=10, color=INK, loc='left')
axes[1].set_xlabel(r'turning angle $\theta_k$ at a vertex (°), 5° bins', fontsize=9, color=INK2)
axes[1].set_ylabel('vertices per unit trace length (1/R)', fontsize=9, color=INK2)
axes[1].set_title('Turning-angle distribution (group, per unit trace length)', fontsize=10, color=INK, loc='left')
for ax in axes:
    ax.axvline(ADOPTED['i'], color=INK, lw=0.9, ls='--'); ax.axvline(V10['i'], color=INK3, lw=0.8, ls=':')
    style(ax); ax.grid(False); ax.legend(fontsize=8, frameon=False)
fig.savefig(os.path.join(SDIR, 'angle_distribution.png'), bbox_inches='tight')
plt.close(fig)

# ---- figure for the choice of the kink threshold
if len(order) == 2:
    a_, b_ = order
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), dpi=150, gridspec_kw=dict(wspace=0.3, width_ratios=[1.25, 1, 1]))
    ax = axes[0]
    centers = bins[:-1] + 2.5
    for gi, g in enumerate(order):
        ax.bar(centers + (gi - 0.5) * 2.1, Hg[g], width=2.0, color=GROUP_COL[gi], label=lab(g), zorder=2)
    ax.axvspan(0, 10, color='#f1f0ec', zorder=0); ax.axvspan(10, 20, color='#e6e5e0', zorder=0)
    for x, y, name in ((5, 1.01, 'common'), (15, 1.065, 'transition')):
        ax.text(x, y, name, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8, color=INK2)
    ax.text(22, 0.97, r'kink ($\theta_k \geq$ %g°)' % ADOPTED['i'], transform=ax.get_xaxis_transform(), ha='left', va='top',
            fontsize=8.5, color=INK, zorder=4)
    ax.axvline(ADOPTED['i'], color=INK, lw=1.0, ls='--')
    ax.set_xlim(0, 100); ax.set_xticks(range(0, 101, 10))
    ax.set_xlabel(r'turning angle $\theta_k$ at a vertex (°), 5° bins', fontsize=9, color=INK2)
    ax.set_ylabel('vertices per unit trace length (1/R)', fontsize=9, color=INK2)
    ax.set_title('(a) Turning-angle distribution', fontsize=10, color=INK, loc='left', pad=30)
    style(ax); ax.legend(fontsize=8, frameon=False, loc='upper right')

    ax = axes[1]
    ra, rb = Hg[a_], Hg[b_]
    na = np.histogram(Kg[a_], bins)[0]; nb = np.histogram(Kg[b_], bins)[0]
    ok = (ra > 0) & (centers < 100)
    few = ok & ((na < 5) | (nb < 5))
    ax.plot(centers[ok], rb[ok] / ra[ok], color=INK, lw=1.2, zorder=2)
    ax.plot(centers[ok & ~few], (rb / np.where(ra > 0, ra, 1))[ok & ~few], 'o', ms=5, color=INK, zorder=3,
            label='≥ 5 vertices in each group')
    ax.plot(centers[few], (rb / np.where(ra > 0, ra, 1))[few], 'o', ms=5, mfc='white', mec=INK, zorder=3,
            label='< 5 vertices in a group')
    lg = ax.legend(fontsize=7.5, frameon=True, loc='lower right', facecolor='white', edgecolor='none', framealpha=1)
    zero = (ra == 0) & (rb > 0) & (centers < 100)
    if zero.any():
        top = 2 ** np.ceil(np.log2(np.nanmax(rb[ok] / ra[ok]) * 1.5))
        ax.plot(centers[zero], np.full(zero.sum(), top), marker='^', ls='none', color=INK2, ms=5)
        ax.text(centers[zero][0], top * 1.25, '%s = 0' % a_, fontsize=7.5, color=INK2, ha='left')
    ax.set_yscale('log', base=2); ax.axhline(1, color=INK3, lw=0.9)
    ax.axvspan(0, 10, color='#f1f0ec', zorder=0); ax.axvspan(10, 20, color='#e6e5e0', zorder=0)
    ax.axvline(ADOPTED['i'], color=INK, lw=1.0, ls='--')
    ax.set_xlim(0, 100); ax.set_xticks(range(0, 101, 10))
    ax.set_xlabel(r'$\theta_k$ (°), 5° bins', fontsize=9, color=INK2)
    ax.set_ylabel('%s / %s (vertices per unit length)' % (b_, a_), fontsize=9, color=INK2)
    ax.set_title('(b) Ratio per bin (pooled)', fontsize=10, color=INK, loc='left', pad=30)
    style(ax)

    ax = axes[2]
    for bi, (lo, hi) in enumerate(BANDS):
        v = band_vals[(lo, hi)]
        for gi, g in enumerate(order):
            x = bi + (gi - 0.5) * 0.36
            ax.bar(x, v[g].mean(), width=0.32, color=GROUP_COL[gi], alpha=0.22, zorder=1)
            ax.plot(x + np.linspace(-0.08, 0.08, len(v[g])), v[g], 'o', ms=6, mfc=GROUP_COL[gi], mec='white', mew=1.1, zorder=3,
                    label=lab(g) if bi == 0 else None)
        c = compare(v[a_], v[b_])
        ax.text(bi, max(v[a_].max(), v[b_].max()) * 1.08 + 0.15, 'MW p = %.2f%s' % (c['mw_p'], '\nseparated' if c['separated'] else ''),
                ha='center', va='bottom', fontsize=7.5, color=INK2)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(['0–10°', '10–20°', '≥ 20°'], fontsize=9, color=INK)
    ax.set_ylabel('vertices per unit trace length (1/R)', fontsize=9, color=INK2)
    ax.set_title('(c) Specimens per angle band', fontsize=10, color=INK, loc='left', pad=30)
    ax.set_ylim(0, max(max(v.max() for v in band_vals[b].values()) for b in BANDS) * 1.35)
    style(ax); ax.tick_params(axis='x', length=0); ax.legend(fontsize=8, frameon=False, loc='upper left')
    fig.savefig(os.path.join(SDIR, 'threshold_justification.png'), bbox_inches='tight')
    plt.close(fig)

for r in srows:
    if 'sigma' not in r['criterion'] or 'sigma=0.06' in r['criterion']:
        print('%-4s %-26s %s' % (r['metric'], r['criterion'], '  '.join('%s=%s' % (k, r[k]) for k in
              ('%s_mean' % order[0], '%s_mean' % order[1], 'ratio', 'welch_p', 'mannwhitney_exact_p', 'complete_separation')
              if k in r)))
for r in band_rows:
    print('band', r)
print('saved in', SDIR)
