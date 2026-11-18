Analysis code for the wrinkle-contour kink statistics in:

J. Lim et al., "Porosity-Mediated Conformality of Nanomesh Electrodes to Skin Roughness and Curvature," Advanced Materials Interfaces (submitted, 2026).

This release is the version used to produce Figure 4b,c, Figures S5, S7, S8 and S9, and Table S3, as described in Supporting Notes S1 and S2.

## Changes from version 1.0.0

* Traced elements are chained into polylines only within the CT section they were traced on. Version 1.0.0 chained the elements of all ten sections of a specimen together, which joined end points of different sections and created spurious reversals.
* Statistic (i), formerly *sharp cusps*, is renamed *kinks*. The kink threshold and the turning-intensity threshold are command-line options, both 20° by default (version 1.0.0: 45° and 30°).
* New scripts for several specimens in one PowerPoint file, replicate specimens, sensitivity to the criteria, and the effect of the chaining mode.
* `kink_maps.py --chain set --kink-deg 45 --intensity-deg 30` reproduces version 1.0.0.

## Contents

* `split_pptx_sets.py` – splits one PowerPoint file with several specimens into one input folder per specimen
* `parse_traces.py` – extracts the manually traced contour elements from the annotated section slides (unchanged)
* `circle_map.py` – detects the reconstruction field of view and normalizes coordinates to the unit disc (radius R) (unchanged)
* `kink_maps.py` – chains polylines, computes the kink statistics (i)–(iii), and renders the kink-intensity maps
* `kink_composites.py` – assembles the composite figures and the summary chart
* `kink_replicates.py` – composites and group statistics of replicate specimens (mean ± SD, Welch t, exact Mann–Whitney U, Hedges g)
* `kink_sensitivity.py` – recomputes the statistics over a range of criteria and gives the turning-angle distribution
* `compare_chain_modes.py` – compares chaining within and across sections and lists the joints created between sections

## Parameters

Chaining within each section; chain tolerance 0.02 R; joint merge 0.01 R; sampling step 0.004 R; Gaussian width σ = 0.06 R; kink threshold θ ≥ 20°; turning-intensity threshold I(s) > 20°.

## Notes

* The traced CT sections are not included; they are available from the corresponding authors upon reasonable request.
* See README.md for the run order, input format, and requirements (Python 3, NumPy, SciPy, scikit-image, Matplotlib, Pillow).
