# Kink maps of CT wrinkle-contour traces

Analysis code for the wrinkle-contour kink statistics (Figure 4b,c; Figures S5 and S7; Supporting Notes S1 and S2) in:

> J. Lim *et al.*, "Porosity-Mediated Conformality of Nanomesh Electrodes to Skin Roughness and Curvature," *Advanced Materials Interfaces* (submitted, 2026).

The scripts take wrinkle contours traced manually as straight line elements on synchrotron X-ray micro-CT sections, chain them into polylines, measure the turning angle at every joint, and report how often and how sharply the contours change direction. The output is three pooled statistics per specimen and a kink-intensity map.

## Specimens

Each set is one specimen on one jig curvature, with ten CT sections. The set names in the code follow the internal `TxDy` labels.

| Set (code) | Label in the manuscript | Reported in the manuscript |
|---|---|---|
| `[Concave][T1D3]` | D3 (High junction density) | Yes |
| `[Concave][T1D2]` | D2 (Mid) | Yes |
| `[Concave][T1D1]` | D1 (Low) | Yes |
| `[Convex][T1D1–T1D3]`, `[Concave/Convex][T2D1, T3D1]` | — | No |

The scripts process all ten available sets (five specimens × concave/convex). Only the three concave density-series sets above are reported in the manuscript; the remaining sets are processed by the same pipeline and are not part of the published results.

## Method

All lengths are normalized by the radius R of the circular field of view of the CT reconstruction, which is detected in each section (unit-disc coordinates). Parameters are defined at the top of `kink_maps.py`.

| Step | Parameter | Value |
|---|---|---|
| Chain traced elements into polylines when end points coincide | `CHAIN_TOL` | 0.02 R |
| Merge joints closer than | `MIN_SEG` | 0.01 R |
| Turning angle θ at each interior joint | — | 0° = straight, 90° = right angle, 180° = reversal |
| Resampling step along each polyline | `DS` | 0.004 R |
| Local turning intensity I(s) = Σ θ_k exp[−(s − s_k)² / (2σ²)] | `SIGMA_S` | σ = 0.06 R |
| Color scale of I(s) (HSB, S 0.8, B 1.0) | `KMAX` | 0° blue (hue 225°) → ≥60° red (hue 0°) |

Pooled statistics per set (all ten sections of one specimen), as defined in Supporting Note S1:

| Manuscript statistic | Column in `kink_metrics.csv` |
|---|---|
| (i) Sharp cusps (θ ≥ 45°) per unit trace length (1/R) | `cusps45_per_R` |
| (ii) Total turning per unit trace length (°/R) | `turning_deg_per_R` |
| (iii) Length fraction with I(s) > 30° | `frac_len_I_gt30` |

The CSV also contains auxiliary columns that are not used in the manuscript: `n_polylines`, `n_joints`, `n_cusps_ge45`, `n_cusps_ge30` (count of joints with θ ≥ 30°), `total_length_R`, and `mean_I`.

The kink-intensity surface is an adaptive-kernel local mean of I(s) over the disc. It is an interpolated visualization and is not used for quantitative comparison. Statistics (i)–(iii) are. All values depend on the granularity of the manual tracing.

## Pipeline

Run the scripts in this order from one folder. Paths are relative to the current directory; the defaults are shown and can be changed with the options.

```
python parse_traces.py   --input raw/unpacked --work work
python circle_map.py     --work work
python kink_maps.py      --work work --out out_kink
python kink_composites.py --out out_kink
```

`kink_maps.py` also accepts set names to process only those sets, e.g. `python kink_maps.py Concave_T1D1`. Rows of `kink_metrics.csv` for the processed sets are replaced and the other rows are kept, so reruns do not create duplicate rows. Each script prints its usage with `-h`.

| Script | Input | Output |
|---|---|---|
| `parse_traces.py` | Annotated section slides (see Input) | `traces.json`: traced line elements and the section image of each slide |
| `circle_map.py` | `traces.json` | `traces_circ.json`: elements in unit-disc coordinates; overlay sheets in `work/check/` (traced lines in red, detected field of view in green) for visual checking |
| `kink_maps.py` | `traces_circ.json` | `kink_metrics.csv`, `lines/`, `filled/`, `filled_faded/`, `qa/` |
| `kink_composites.py` | Outputs of `kink_maps.py` | Composite figures and `kink_metrics_chart.png` |

### Output files of `kink_maps.py`

| Folder | Files | Content |
|---|---|---|
| `lines/` | `kink_lines_<set>.png`, `kink_lines_plain_<set>.png`, `kink_lines_marked_<set>.png` | Polylines colored by I(s): annotated, plain, and with rings at cusps ≥ 60° |
| `filled/` | `kink_fill_<set>.png` | Fully interpolated disc |
| `filled_faded/` | `kink_fill_faded_<set>.png` | Interpolated disc faded to white away from traces (colors there are extrapolated) |
| `qa/` | `qa_kink_fill_<set>.png` | Interpolated disc with the polylines overlaid |

## Input

The traced sections are not included in this repository. They are available from the corresponding authors upon reasonable request, as stated in the Data Availability Statement of the manuscript.

Each set is a PowerPoint file named `[Curvature][TxDy]_raw.pptx` (e.g. `[Concave][T1D1]_raw.pptx`) with ten slides. Each slide holds one CT section image and the wrinkle contour traced over it as straight-line connector shapes. Only red (`FF0000`) straight-line connector shapes are read as traces; other shapes are reported and ignored. The first slide of each file is treated as a title slide. Unzip each file into its own folder of the same name under the input directory (default `raw/unpacked/`) before running `parse_traces.py`:

```
raw/unpacked/[Concave][T1D1]_raw/ppt/...
raw/unpacked/[Concave][T1D2]_raw/ppt/...
```

## Requirements

Python 3 with NumPy, SciPy, scikit-image, Matplotlib, and Pillow.

```
pip install numpy scipy scikit-image matplotlib pillow
```

The pipeline was checked end to end on synthetic traces; refactoring for this release (command-line paths, metrics-table handling) did not change any computed value or map.

## Notes

The scripts were written with the assistance of a generative AI tool (Claude, Anthropic) and were reviewed and verified by the authors, as disclosed in the Acknowledgements of the manuscript.

## Citation

If you use this code, please cite the article above and the archived version of this repository:

> J. Lim, *Kink-map analysis of CT contour traces*, version 1.0, Zenodo (2026). https://doi.org/10.5281/zenodo.XXXXXXX

## License

MIT License. See `LICENSE`.
