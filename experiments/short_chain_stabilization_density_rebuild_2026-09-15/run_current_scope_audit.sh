#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 audit_current_scope.py \
  > analysis/current_scope_audit.log 2>&1
cat analysis/current_scope_audit.log
