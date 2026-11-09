#!/usr/bin/env bash
# Install the locked Linux environment inside this package; no sudo required.
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'This installer targets Linux x86_64 (Intel/AMD), not macOS or Graviton.' >&2
  exit 2
fi
glibc_version="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')"
if [[ -z "$glibc_version" ]]; then
  echo 'A glibc-based Linux distribution is required.' >&2
  exit 2
fi
if [[ "$(printf '%s\n' 2.28 "$glibc_version" | sort -V | head -n 1)" != 2.28 ]]; then
  echo "glibc >=2.28 is required; found $glibc_version. Use a newer Linux AMI." >&2
  exit 2
fi
for command in curl tar bzip2; do
  command -v "$command" >/dev/null || { echo "Missing system command: $command" >&2; exit 2; }
done
mkdir -p "$package_dir/.tools"
micromamba_binary="$package_dir/.tools/bin/micromamba"
if [[ ! -x "$micromamba_binary" ]]; then
  curl --fail --location --retry 3 --connect-timeout 20 \
    https://micro.mamba.pm/api/micromamba/linux-64/2.9.0 \
    --output "$package_dir/.tools/micromamba.tar.bz2"
  tar -xjf "$package_dir/.tools/micromamba.tar.bz2" -C "$package_dir/.tools" bin/micromamba
fi
if [[ -e "$package_dir/.env" ]]; then
  echo 'The .env directory already exists. Checking the existing environment.'
else
  MAMBA_ROOT_PREFIX="$package_dir/.micromamba" "$micromamba_binary" create --yes \
    --prefix "$package_dir/.env" --file "$package_dir/ec2-linux-64.explicit.txt"
fi
MAMBA_ROOT_PREFIX="$package_dir/.micromamba" "$micromamba_binary" run \
  --prefix "$package_dir/.env" python "$package_dir/run_pipeline.py" --doctor
echo
echo 'Environment is ready. Run:'
printf '  %q %q\n' "$package_dir/.env/bin/python" "$package_dir/smoke_test.py"
