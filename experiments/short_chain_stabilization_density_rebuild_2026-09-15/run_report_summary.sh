#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
/usr/bin/python3 build_report_summary.py > analysis/report_summary.log 2>&1
cat analysis/report_summary.log
