"""
postprocessing_vertebrae.py

Post-processing for vertebrae segmentation predictions (AbdomenAtlas-style
NIfTI label maps). Cleans up raw model predictions by:

  1. Removing spurious small connected components per label
     (keeps only the components that make up a real anatomical structure).
  2. Filling small internal holes within each label's mask.
  3. Optionally applying a light morphological closing to reconnect
     components that were split by thin gaps from imperfect predictions.

Usage:
    python postprocessing_vertebrae.py --input <in.nii.gz> --output <out.nii.gz>
    python postprocessing_vertebrae.py --input_dir <dir_of_case_folders> --output_dir <out_dir>

Each case folder under --input_dir is expected to contain a segmentation
NIfTI file (default name: seg.nii.gz, override with --seg_name). Output
folders mirror the same case names and each contains seg_refined.nii.gz.
"""

import argparse
import os
import glob

import numpy as np
import nibabel as nib
from scipy import ndimage


def keep_significant_components(mask: np.ndarray, min_fraction: float = 0.05, min_voxels: int = 20):
    """
    Given a binary mask for a single label, keep only connected components
    that are either:
      - the largest component, or
      - at least `min_fraction` of the largest component's size
    and drop anything smaller than `min_voxels` outright (near-certain noise).

    This is less aggressive than "keep only the single largest component",
    which matters for structures that legitimately have multiple separate
    pieces in a single slab (e.g. a vertebra cut by the volume boundary).
    """
    labeled, n = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    if n == 0:
        return mask

    sizes = ndimage.sum(mask, labeled, range(1, n + 1))
    max_size = sizes.max()

    keep_ids = [
        i + 1
        for i, s in enumerate(sizes)
        if s >= min_voxels and s >= min_fraction * max_size
    ]

    cleaned = np.isin(labeled, keep_ids)
    return cleaned


def fill_small_holes(mask: np.ndarray):
    """Fill internal holes within a binary mask, slice-wise and in 3D."""
    filled = ndimage.binary_fill_holes(mask)
    return filled


def close_gaps(mask: np.ndarray, iterations: int = 1):
    """Light morphological closing to reconnect near-touching fragments."""
    structure = ndimage.generate_binary_structure(3, 1)
    closed = ndimage.binary_closing(mask, structure=structure, iterations=iterations)
    return closed


def refine_segmentation(data: np.ndarray, min_fraction: float = 0.05,
                         min_voxels: int = 20, closing_iterations: int = 1,
                         verbose: bool = True):
    """
    Run the full cleanup pipeline on a multi-label integer segmentation array.
    Returns a new array of the same shape/dtype.
    """
    refined = np.zeros_like(data)
    labels = [l for l in np.unique(data) if l != 0]

    for lab in labels:
        mask = data == lab

        before_voxels = int(mask.sum())
        before_n = ndimage.label(mask, structure=np.ones((3, 3, 3)))[1]

        if closing_iterations > 0:
            mask = close_gaps(mask, iterations=closing_iterations)

        mask = keep_significant_components(mask, min_fraction=min_fraction, min_voxels=min_voxels)
        mask = fill_small_holes(mask)

        after_voxels = int(mask.sum())
        after_n = ndimage.label(mask, structure=np.ones((3, 3, 3)))[1]

        if verbose:
            print(f"  label {int(lab):>3}: components {before_n} -> {after_n}, "
                  f"voxels {before_voxels} -> {after_voxels}")

        refined[mask] = lab

    return refined.astype(data.dtype)


def process_one_file(input_path: str, output_path: str, **kwargs):
    img = nib.load(input_path)
    data = img.get_fdata().astype(np.uint8)

    print(f"Processing {input_path}")
    refined = refine_segmentation(data, **kwargs)

    out_img = nib.Nifti1Image(refined, affine=img.affine, header=img.header)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    nib.save(out_img, output_path)
    print(f"  saved -> {output_path}\n")


def compute_case_stats(seg_path: str):
    """Compute per-label voxel count and z-centroid (spine axis) for one case."""
    img = nib.load(seg_path)
    data = img.get_fdata()
    stats = {}
    for lab in np.unique(data):
        if lab == 0:
            continue
        zs = np.where(data == lab)[2]
        stats[int(lab)] = {
            "voxels": int((data == lab).sum()),
            "centroid_z": float(zs.mean()),
        }
    return stats


def run_qc_report(refined_dir: str, seg_name: str, report_path: str):
    """
    Cross-case quality-control report.

    Two data-driven plausibility checks (no external anatomical atlas assumed,
    since it's derived only from the population of cases actually submitted):

      1. Volume outliers: for each label, flag any case whose voxel count is
         more than 2 standard deviations from the mean voxel count observed
         for that same label across all cases.

      2. Ordering consistency: vertebral labels should stack in a consistent
         order along the spine (z) axis across cases. For every pair of labels
         that co-occur in 2+ cases, take a majority vote on which one sits
         more superior, then flag any case that disagrees with the majority.
    """
    case_dirs = sorted(
        d for d in glob.glob(os.path.join(refined_dir, "*")) if os.path.isdir(d)
    )

    all_stats = {}
    for case_dir in case_dirs:
        case_name = os.path.basename(case_dir)
        seg_path = os.path.join(case_dir, seg_name)
        if os.path.exists(seg_path):
            all_stats[case_name] = compute_case_stats(seg_path)

    lines = []
    lines.append("QC Report: postprocessing_vertebrae.py")
    lines.append("=" * 60)
    lines.append(f"Cases analyzed: {len(all_stats)}")
    lines.append("")

    lines.append("1) Volume outlier check (per label, across cases)")
    lines.append("-" * 60)
    label_volumes = {}
    for case, stats in all_stats.items():
        for lab, s in stats.items():
            label_volumes.setdefault(lab, []).append((case, s["voxels"]))

    any_outlier = False
    for lab, entries in sorted(label_volumes.items()):
        vols = np.array([v for _, v in entries])
        mean, std = vols.mean(), vols.std()
        lines.append(f"  label {lab}: mean={mean:.0f} voxels, std={std:.0f}, n={len(vols)}")
        if std > 0:
            for case, v in entries:
                z = (v - mean) / std
                if abs(z) > 2:
                    any_outlier = True
                    lines.append(f"    FLAG: {case} label {lab} = {v} voxels (z-score {z:+.1f})")
    if not any_outlier:
        lines.append("  No volume outliers beyond 2 std across the submitted cases.")
    lines.append("")

    lines.append("2) Vertebral ordering consistency check (spine axis)")
    lines.append("-" * 60)
    pair_votes = {}
    for case, stats in all_stats.items():
        labels = list(stats.keys())
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, b = labels[i], labels[j]
                key = tuple(sorted((a, b)))
                a_above_b = stats[a]["centroid_z"] > stats[b]["centroid_z"]
                pair_votes.setdefault(key, []).append((case, a_above_b if key[0] == a else not a_above_b))

    any_violation = False
    for (a, b), votes in sorted(pair_votes.items()):
        if len(votes) < 2:
            continue
        true_count = sum(1 for _, v in votes if v)
        majority = true_count >= len(votes) / 2
        for case, v in votes:
            if v != majority:
                any_violation = True
                lines.append(f"  FLAG: {case} - label {a} vs {b} ordering disagrees with majority "
                              f"({true_count}/{len(votes)} cases agree label {a} is more superior)")
    if not any_violation:
        lines.append("  No ordering inconsistencies detected across co-occurring label pairs.")
    lines.append("")

    report = "\n".join(lines)
    with open(report_path, "w") as f:
        f.write(report)
    print(report)
    print(f"\nSaved QC report -> {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Post-process vertebrae segmentation predictions.")
    parser.add_argument("--input", type=str, default=None, help="Path to a single input .nii.gz file")
    parser.add_argument("--output", type=str, default=None, help="Path to save the refined .nii.gz file")
    parser.add_argument("--input_dir", type=str, default=None,
                         help="Directory containing per-case subfolders, each with a segmentation file")
    parser.add_argument("--output_dir", type=str, default=None,
                         help="Directory to write refined per-case subfolders into")
    parser.add_argument("--seg_name", type=str, default="seg.nii.gz",
                         help="Filename of the segmentation inside each case folder (default: seg.nii.gz)")
    parser.add_argument("--out_name", type=str, default="seg_refined.nii.gz",
                         help="Filename to save the refined segmentation as (default: seg_refined.nii.gz)")
    parser.add_argument("--min_fraction", type=float, default=0.05,
                         help="Minimum size (as a fraction of the largest component) for a component to be kept")
    parser.add_argument("--min_voxels", type=int, default=20,
                         help="Absolute minimum voxel count for a component to be kept")
    parser.add_argument("--closing_iterations", type=int, default=1,
                         help="Number of binary closing iterations to reconnect fragmented pieces (0 to disable)")
    parser.add_argument("--qc_report", action="store_true",
                         help="After batch processing, run cross-case QC checks (volume outliers + "
                              "vertebral ordering consistency) and save a report")
    parser.add_argument("--qc_report_path", type=str, default=None,
                         help="Where to save the QC report (default: <output_dir>/qc_report.txt)")

    args = parser.parse_args()

    kwargs = dict(
        min_fraction=args.min_fraction,
        min_voxels=args.min_voxels,
        closing_iterations=args.closing_iterations,
    )

    if args.input:
        output_path = args.output or args.input.replace(".nii.gz", "_refined.nii.gz")
        process_one_file(args.input, output_path, **kwargs)

    elif args.input_dir:
        if not args.output_dir:
            raise ValueError("--output_dir is required when using --input_dir")

        case_dirs = sorted(
            d for d in glob.glob(os.path.join(args.input_dir, "*"))
            if os.path.isdir(d)
        )

        for case_dir in case_dirs:
            case_name = os.path.basename(case_dir)
            seg_path = os.path.join(case_dir, args.seg_name)
            if not os.path.exists(seg_path):
                print(f"Skipping {case_name}: no {args.seg_name} found")
                continue

            output_path = os.path.join(args.output_dir, case_name, args.out_name)
            process_one_file(seg_path, output_path, **kwargs)

        if args.qc_report:
            report_path = args.qc_report_path or os.path.join(args.output_dir, "qc_report.txt")
            run_qc_report(args.output_dir, args.out_name, report_path)

    else:
        parser.error("Provide either --input/--output for a single file, "
                     "or --input_dir/--output_dir for a batch of cases.")


if __name__ == "__main__":
    main()
