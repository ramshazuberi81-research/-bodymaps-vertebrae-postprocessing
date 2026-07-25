Vertebrae Segmentation Post-processing — BodyMaps Submission

Submitted by: Ramsha Zuberi
Files: `AbdomenAtlasDemoPredict/` (refined segmentations, 5 cases), `postprocessing_vertebrae.py`, `qc_report.txt`, this README

---

 Pipeline overview

`postprocessing_vertebrae.py` takes a raw multi-label vertebrae segmentation
(NIfTI) and runs it through four stages: three geometric cleanup steps, and
one data-driven QC layer.

![Pipeline diagram](assets/pipeline_diagram.png)

1. Morphological closing — reconnects components that got split by thin
   gaps in the raw prediction (a single vertebra should usually be one
   connected piece).
   
3. Component filtering — keeps a connected component if it's either the
   largest piece for that label, or at least 5% of the largest piece's size.
   Anything below an absolute floor (20 voxels) is dropped as near-certain
   noise regardless of relative size. This is deliberately *not* "keep only
   the single largest component" — some labels legitimately have more than
   one substantial piece (see sub-verse502 below).
   
4. Hole filling — fills small internal gaps left inside an otherwise
   solid mask.
   
6. QC checks — cross-case volume outlier detection and vertebral
   ordering consistency (details below).

---

Result: before vs. after

The visible fragmentation and noisy specks from the raw predictions are
cleaned up across all 5 cases:

![Before/after slices](assets/before_after_slices.png)

Quantified as connected-component count per label:

![Component counts before and after](assets/component_counts.png)

Most labels start out as a single clean component. Where fragmentation
existed (sub-verse501 label 75: 3 pieces → 1; sub-verse504 label 76: 2 → 1;
sub-verse505 label 76: 2 → 1), the pipeline merges the real piece back
together and drops the noise. Sub-verse502 label 75 stays at 2 components
on purpose  both pieces are substantial (258 and 157 voxels), so the
script doesn't guess which one is "real" and discards neither. I'd rather
flag genuine ambiguity than silently delete a legitimate piece.

---

Beyond geometric cleanup: data-driven plausibility checks

Geometric cleanup alone can't tell you whether a segmentation is
*anatomically* sensible — a smooth, single-component blob in the wrong place
or the wrong size is still wrong. Since I didn't have access to an external
anatomical atlas or ground truth for these 5 cases, I built two checks that
are data-driven from the submitted population itself, rather than
hardcoded assumptions:

- Volume outlier check — for each label, flags any case whose voxel
  count is more than 2 standard deviations from the mean for that label
  across the other submitted cases.
- Vertebral ordering consistency check — vertebrae should stack in a
  consistent order along the spine axis. For every pair of labels that
  co-occur in 2+ cases, the script takes a majority vote on which one sits
  more superior, then flags any case that disagrees with the majority. This
  would catch, e.g., a labeling error where two adjacent vertebrae got
  swapped.

![Volume outlier check per label](assets/volume_outlier_check.png)

Both checks ran clean on this submission (full detail in `qc_report.txt`) —
no volume outliers and no ordering violations across the 5 cases. That's a
real result from the data, not something engineered to pass.

---

 Usage

```bash
# single file
python postprocessing_vertebrae.py --input seg.nii.gz --output seg_refined.nii.gz

# batch, with QC report
python postprocessing_vertebrae.py \
  --input_dir AbdomenAtlasDemoPredict_raw \
  --output_dir AbdomenAtlasDemoPredict \
  --seg_name seg.nii.gz --out_name seg_refined.nii.gz \
  --qc_report
```

---

 Honest limitations

- No external ground truth was available for these 5 cases, so the
  plausibility checks are relative to the submitted population, not an
  absolute anatomical reference. With more cases (or a reference atlas),
  the volume-outlier check would get meaningfully more reliable.
- The 5% component-size threshold and 2-std outlier cutoff are reasonable
  defaults, not tuned on a validation set — an obvious next step if this
  becomes an ongoing project would be to tune these against labeled data.
- Sub-verse502's two-component case (above) couldn't be disambiguated as
  noise vs. real anatomy without ground truth, so both were kept rather
  than guessed at.

---

 Why this matters for BodyMaps

Scaling annotation across large datasets means manual QA can't scale
linearly with case count. Automated plausibility checks like these  even
simple, population-relative ones  are the kind of lightweight signal that
could help triage which cases need a human to look twice, without needing a
full second model or an external atlas. That connects directly to the
"Scaling Annotations" effort described in the program: catching likely
errors automatically is what makes scaling to thousands of cases feasible
with a small mentored team.
