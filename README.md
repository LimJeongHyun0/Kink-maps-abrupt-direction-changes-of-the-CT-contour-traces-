[README.md](https://github.com/user-attachments/files/35117764/README.md)
[README.md](https://github.com/user-attachments/files/35117778/README.md)
# Kink maps of CT wrinkle-contour traces

Analysis code for the wrinkle-contour kink statistics (Figure 4; Supporting Notes S1 and S2 and the related Supporting Information) in:

> J. Lim *et al.*, "Porosity-Mediated Conformality of Nanomesh Electrodes to Skin Roughness and Curvature," 

The scripts take wrinkle contours traced manually as straight line elements on synchrotron X-ray micro-CT sections and chain them into polylines. They then measure the turning angle at every joint and report how often and how strongly the contours change direction. The output is three pooled statistics per specimen and a kink-intensity map.

## Version 1.1.0

Changes from version 1.0.0. The changes affect computed values; the version 1.0 behaviour can still be selected (see *Reproducing version 1.0*).

* **Chaining within each section.** Traced elements are now chained into polylines only within the section (slide) they were traced on (`--chain section`, default). Version 1.0 chained the elements of all ten sections of a specimen together. End points of different sections that happened to lie within `CHAIN_TOL` of each other were then joined, typically two contours ending at the same point of the field-of-view edge, and created spurious reversals of 150–175°. `compare_chain_modes.py` lists these joints.
* **Kinks.** Statistic (i), formerly *sharp cusps* (θ ≥ 45°), is now *kinks*: joints with θ ≥ θ_kink. Statistic (iii) uses I(s) > θ_I. The thresholds are options, `--kink-deg` and `--intensity-deg`, both 20° by default (version 1.0: 45° and 30°). The choice of 20° and the sensitivity of the statistics to all criteria are documented in Supporting Note S1 and produced by `kink_sensitivity.py`.
* **Metrics table.** `kink_metrics.csv` columns were renamed: `n_cusps_ge45` → `n_kinks`, `cusps45_per_R` → `kinks_per_R`, `frac_len_I_gt30` → `frac_len_I_gt`; `n_cusps_ge30` was removed. The criteria of the run are written to every row (`kink_deg`, `intensity_deg`, `chain`), and rows computed with other criteria are not mixed into the table. Values are written with more decimals.
* **Maps.** `kink_lines_marked_<set>.png` rings every kink (version 1.0: joints ≥ 60°).
* **New scripts.** `split_pptx_sets.py` (several specimens in one PowerPoint file), `kink_replicates.py` (replicate specimens: composites and group statistics), `kink_sensitivity.py` (sensitivity analysis and turning-angle distribution), `compare_chain_modes.py` (effect of the chaining mode).
* `parse_traces.py` and `circle_map.py` are unchanged.

## Specimens

Each set is one specimen with ten CT sections. The set names follow the input file names.

| Data | Sets (code) | Label in the manuscript |
|---|---|---|
| Density series, one specimen per mesh | `Concave_T1D3`, `Concave_T1D2`, `Concave_T1D1` | D3 (High), D2 (Mid), D1 (Low) |
| Replicate specimens, three per mesh | `High_1`–`High_3`, `Low_1`–`Low_3` | D3 (High), D1 (Low) |

The scripts also process the other sets of the density and diameter series (`Convex_T1D1`–`T1D3`, `Concave/Convex_T2D1`, `T3D1`); these are not part of the published results.

## Method

All lengths are normalized by the radius R of the circular field of view of the CT reconstruction, which is detected in each section (unit-disc coordinates). Parameters are defined at the top of `kink_maps.py` or set on its command line.

| Step | Parameter | Value |
|---|---|---|
| Chain traced elements of one section into polylines when end points coincide | `CHAIN_TOL` | 0.02 R |
| Merge joints closer than | `MIN_SEG` | 0.01 R |
| Turning angle θ at each interior joint | — | 0° = straight, 90° = right angle, 180° = reversal |
| Kink | `--kink-deg` | θ ≥ 20° |
| Resampling step along each polyline | `DS` | 0.004 R |
| Local turning intensity I(s) = Σ θ_k exp[−(s − s_k)² / (2σ²)] | `SIGMA_S` | σ = 0.06 R |
| Threshold of statistic (iii) | `--intensity-deg` | I(s) > 20° |
| Color scale of I(s) (HSB, S 0.8, B 1.0) | `KMAX` | 0° blue (hue 225°) → ≥60° red (hue 0°) |

Pooled statistics per set (all ten sections of one specimen), as defined in Supporting Note S1:

| Manuscript statistic | Column in `kink_metrics.csv` |
|---|---|
| (i) Kinks per unit trace length (1/R) | `kinks_per_R` |
| (ii) Total turning per unit trace length (°/R) | `turning_deg_per_R` |
| (iii) Length fraction with I(s) > θ_I | `frac_len_I_gt` |

The CSV also contains `n_polylines`, `n_joints`, `n_kinks`, `total_length_R`, `mean_I` and the criteria of the run (`kink_deg`, `intensity_deg`, `chain`).

The kink-intensity surface is an adaptive-kernel local mean of I(s) over the disc. It is an interpolated visualization and is not used for quantitative comparison. Statistics (i)–(iii) are. All values depend on the granularity of the manual tracing.

## Pipeline

Run the scripts in this order from one folder. Paths are relative to the current directory; the defaults are shown and can be changed with the options. Each script prints its usage with `-h`.

```
python split_pptx_sets.py traces.pptx --out raw/unpacked --map D1=Low D3=High   # only for several specimens in one file
python parse_traces.py    --input raw/unpacked --work work
python circle_map.py      --work work
python kink_maps.py       --work work --out out_kink
python kink_composites.py --out out_kink                                       # density / diameter series layout
python kink_replicates.py --work work --out out_kink --groups High Low          # replicate specimens
python kink_sensitivity.py --work work --out out_kink --groups High Low --labels "High=D3 (High)" "Low=D1 (Low)"
python compare_chain_modes.py --work work --out out_kink
```

`kink_maps.py` also accepts set names to process only those sets, e.g. `python kink_maps.py Concave_T1D1`. Rows of `kink_metrics.csv` for the processed sets are replaced and the other rows are kept, so reruns do not create duplicate rows.

| Script | Input | Output |
|---|---|---|
| `split_pptx_sets.py` | One PowerPoint file with several specimens (see Input) | One input folder per specimen for `parse_traces.py` |
| `parse_traces.py` | Annotated section slides (see Input) | `traces.json`: traced line elements and the section image of each slide |
| `circle_map.py` | `traces.json` | `traces_circ.json`: elements in unit-disc coordinates; overlay sheets in `work/check/` (traced lines in red, detected field of view in green) for visual checking |
| `kink_maps.py` | `traces_circ.json` | `kink_metrics.csv`, `lines/`, `filled/`, `filled_faded/`, `qa/` |
| `kink_composites.py` | Outputs of `kink_maps.py` | Composite figures and `kink_metrics_chart.png` for the density and diameter series |
| `kink_replicates.py` | Outputs of `kink_maps.py`, `traces_circ.json` | `replicates/`: composites, statistics chart and table for replicate specimens |
| `kink_sensitivity.py` | `traces_circ.json` | `sensitivity/`: statistics for a range of criteria and the turning-angle distribution |
| `compare_chain_modes.py` | `traces_circ.json` | Statistics with chaining within sections and across sections; joints created between sections |

### Output files of `kink_maps.py`

| Folder | Files | Content |
|---|---|---|
| `lines/` | `kink_lines_<set>.png`, `kink_lines_plain_<set>.png`, `kink_lines_marked_<set>.png` | Polylines colored by I(s): annotated, plain, and with rings at kinks |
| `filled/` | `kink_fill_<set>.png` | Fully interpolated disc |
| `filled_faded/` | `kink_fill_faded_<set>.png` | Interpolated disc faded to white away from traces (colors there are extrapolated) |
| `qa/` | `qa_kink_fill_<set>.png` | Interpolated disc with the polylines overlaid |

### Replicate specimens: `kink_replicates.py`

Sets named `<Group>_<n>` (e.g. `High_1`) are grouped by `<Group>`. Output in `out_kink/replicates/`:

* `composite_kink_lines.png`, `composite_kink_filled.png`, `composite_kink_filled_faded.png`: rows = groups, columns = specimens.
* `kink_metrics_replicates.png`: statistics (i)–(iii), specimen values with mean ± SD.
* `replicate_stats.csv`: per group the specimen values, mean, SD (n − 1) and length-weighted pooled value. Between two groups: ratio of the means, Welch t, exact two-sided Mann–Whitney U and Hedges g.
* `sensitivity_edge_band.csv`: (i)–(iii) recomputed with joints at r ≥ 0.95 R ignored.

`--exclude SET ...` leaves sets out of the statistics; they are drawn as open circles and the output goes to `replicates_excl_<sets>/`.

### Sensitivity to the criteria: `kink_sensitivity.py`

The polylines and turning angles are those of `kink_maps.py`; only the criterion of each statistic is varied. Output in `out_kink/sensitivity/`.

| Statistic | Criterion varied | Values |
|---|---|---|
| (i) | kink threshold θ_kink | 10, 15, 20, 30, 45, 60, 90° |
| (ii) | smallest angle summed θ_min (smaller joints left out) | 0, 5, 10, 15, 20, 30° |
| (iii) | threshold θ_I × Gaussian spread σ | 10, 15, 20, 30, 45, 60° × 0.03, 0.06, 0.12 R |

* `sensitivity_grid.png`, `charts/`: bar chart of every criterion value (adopted criterion framed in black, version 1.0 in grey).
* `sensitivity_trend.png`, `sensitivity_iii_sigma.png`: values, ratio of the group means and p against the criterion.
* `sensitivity_values.csv`, `sensitivity_summary.csv`: specimen values and group statistics for every criterion.
* `threshold_justification.png`, `angle_bins.csv`, `band_rates.csv`, `angle_distribution.png`, `joint_density.csv`: turning-angle distribution per unit trace length, per 5° bin and per angle band (0–10°, 10–20°, ≥ 20°).

At run time the script checks that it reproduces `kink_metrics.csv` at the criteria recorded there.

## Input

The traced sections are not included in this repository. They are available from the corresponding authors upon reasonable request, as stated in the Data Availability Statement of the manuscript.

**One file per specimen.** Each set is a PowerPoint file named `[Curvature][TxDy]_raw.pptx` (e.g. `[Concave][T1D1]_raw.pptx`). The first slide is a title slide, followed by ten slides. Each of those slides holds one CT section image and the wrinkle contour traced over it as straight-line connector shapes. Only red (`FF0000`) straight-line connector shapes are read as traces; other shapes are reported and ignored. Unzip each file into its own folder of the same name under the input directory (default `raw/unpacked/`) before running `parse_traces.py`:

```
raw/unpacked/[Concave][T1D1]_raw/ppt/...
raw/unpacked/[Concave][T1D2]_raw/ppt/...
```

**Several specimens in one file.** For each specimen, a title slide with the specimen code (e.g. `D1_2`) is followed by its section slides (one CT section picture with the red traces each). Slides with traces but no picture (e.g. an overview of all sections) are skipped. `split_pptx_sets.py` writes one input folder per specimen; `--map D1=Low D3=High` turns `D1_2` into `[Low][2]_raw` (set `Low_2`).

## Reproducing version 1.0

```
python kink_maps.py --work work --out out_kink_v10 --chain set --kink-deg 45 --intensity-deg 30
```

This gives the statistics and maps of version 1.0. The column names are those of version 1.1: `kinks_per_R` holds the former `cusps45_per_R`, and `frac_len_I_gt` the former `frac_len_I_gt30`. The plain and filled maps are identical to version 1.0. The annotated maps differ in their text, and the marked maps ring joints ≥ 45° instead of ≥ 60°. This was checked on the traced data.

## Requirements

Python 3 with NumPy, SciPy, scikit-image, Matplotlib, and Pillow.

```
pip install numpy scipy scikit-image matplotlib pillow
```

## Notes

The scripts were written with the assistance of a generative AI tool (Claude, Anthropic) and were reviewed and verified by the authors, as disclosed in the Acknowledgements of the manuscript.

## Citation

If you use this code, please cite the article above and the archived version of this repository:

> J. Lim, *Kink-map analysis of CT contour traces*, version 1.1.0, Zenodo (2026). https://doi.org/10.5281/zenodo.XXXXXXX

## License

MIT License. See `LICENSE`.
