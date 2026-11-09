#!/usr/bin/env bash
# Use the package's installed environment without activating a shell environment.
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -x "$package_dir/.env/bin/python" ]]; then
  echo 'Run bash setup_linux.sh first.' >&2
  exit 2
fi
# External programs spawned by Python must see the same Conda environment.
export MPLCONFIGDIR="$package_dir/.cache/matplotlib"
export LC_ALL=C
export MAMBA_ROOT_PREFIX="$package_dir/.micromamba"
exec "$package_dir/.tools/bin/micromamba" run --prefix "$package_dir/.env" \
  python "$package_dir/run_pipeline.py" "$@"
