#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
OLD="$ROOT/package/five_dimensional_ness_density_2026-09-15"
DEST="$ROOT/package/five_dimensional_ness_density_r11_2026-09-15"
ZIP="$ROOT/package/five_dimensional_ness_density_r11_2026-09-15_publication.zip"

if [[ -e "$DEST" || -e "$ZIP" ]]; then
  echo "refusing to overwrite an existing R11 package" >&2
  exit 2
fi

mkdir -p "$DEST"
rsync -a --exclude '.DS_Store' "$OLD/" "$DEST/"

# Preserve the previous fine two-moment endpoint before replacing the primary
# fine model and audit with the R11 three-moment endpoint.
cp -R "$DEST/models/fine" "$DEST/models/fine_two_moment"
cp -R "$DEST/audit/fine" "$DEST/audit/fine_two_moment"
mkdir -p "$DEST/models/fine" "$DEST/audit/fine" "$DEST/calibration/fine_three_moment"
rsync -a "$ROOT/recovery_r11_fine/models/" "$DEST/models/fine/"
rsync -a "$ROOT/recovery_r11_fine/calibration_raw/" "$DEST/calibration/fine_three_moment/"
for split in driven_validation driven_test equilibrium_validation equilibrium_test; do
  mkdir -p "$DEST/audit/fine/$split"
  rsync -a "$ROOT/recovery_r11_fine/$split/" "$DEST/audit/fine/$split/"
done

mkdir -p "$DEST/r11_analysis" "$DEST/final_timestep_comparison_two_vs_three"
rsync -a "$ROOT/recovery_r11_fine/analysis/" "$DEST/r11_analysis/"
rsync -a "$ROOT/final_timestep_comparison_two_vs_three_2026-09-15/" \
  "$DEST/final_timestep_comparison_two_vs_three/"
rsync -a "$ROOT/final_analysis/" "$DEST/final_analysis/"
rsync -a "$ROOT/final_figures/" "$DEST/final_figures/"
rsync -a "$ROOT/final_timestep_comparison/" "$DEST/final_timestep_comparison/"

cp "$ROOT/FINAL_VERDICT.md" "$DEST/FINAL_VERDICT.md"
cp "$ROOT/VALIDATION_REPORT.md" "$DEST/VALIDATION_REPORT.md"
cp "$ROOT/PACKAGE_README.md" "$DEST/PACKAGE_README.md"
cp "$ROOT/PROVENANCE.md" "$DEST/PROVENANCE.md"
cp "$ROOT/R11_VALIDATION_REPORT.md" "$DEST/R11_VALIDATION_REPORT.md"
cp "$ROOT/R11_COMMANDS.txt" "$DEST/R11_COMMANDS.txt"
cp "$ROOT/R11_SOURCE_HASHES.sha256" "$DEST/R11_SOURCE_HASHES.sha256"
cp "$ROOT/RECOVERY_PROTOCOL_R11_FINE_THREE_MOMENT.md" "$DEST/protocols/"
cp "$ROOT/FORMAL_BLIND_TEST_R11_FINE.md" "$DEST/protocols/"
cp "$ROOT/report/fine_three_moment_calibration_report.tex" "$DEST/report/"
cp "$ROOT/report/fine_three_moment_calibration_report.pdf" "$DEST/report/"
cp "$ROOT/report/five_dimensional_ness_density_report.tex" "$DEST/report/"
cp "$ROOT/report/five_dimensional_ness_density_report.pdf" "$DEST/report/"
cp "$ROOT/paper/section4_1_ness_density.tex" "$DEST/paper/"

for script in \
  calibrate_three_moments_fine.py create_zero_triple_equilibrium.py \
  evaluate_triple_calibrated_v6.py compare_three_moment_timesteps.py \
  summarize_r11_fine.py final_ness_density.py \
  triple_calibrated_transport_ratio_model.py transport_ratio_model.py \
  regular_ratio_model.py; do
  cp "$ROOT/scripts/$script" "$DEST/scripts/$script"
done

(cd "$DEST" && find . -type f ! -name SHA256SUMS.tsv -print0 | sort -z | \
  xargs -0 shasum -a 256 > SHA256SUMS.tsv)
(cd "$ROOT/package" && COPYFILE_DISABLE=1 /usr/bin/zip -qry \
  "$(basename "$ZIP")" "$(basename "$DEST")")
shasum -a 256 "$ZIP" > "$ZIP.sha256"
unzip -t "$ZIP" >/dev/null

echo "$ZIP"
cat "$ZIP.sha256"

