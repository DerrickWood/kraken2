#!/usr/bin/env bash
# Regenerates tests/golden from an installed Kraken 2.
#
# Usage: tests/update_golden.sh INSTALL_DIR
#
# Run this only when a change to classification or report output is intended.
# It reruns tests/e2e.sh with E2E_UPDATE_GOLDEN=1, which copies each output to
# tests/golden instead of comparing, then reruns it normally to confirm the
# new files are stable.  Review the resulting diff before committing; the PR
# summary job in CI shows the same changes read by read.

set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 INSTALL_DIR" >&2
  exit 64
fi

here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/k2golden.XXXXXX")

E2E_UPDATE_GOLDEN=1 E2E_XPASS_OK=1 "$here/e2e.sh" "$1" "$work/update" > "$work/update.log" 2>&1 || true
"$here/e2e.sh" "$1" "$work/verify"

echo
echo "Golden files changed:"
git -C "$here/.." status --short -- tests/golden
git -C "$here/.." diff --stat -- tests/golden
