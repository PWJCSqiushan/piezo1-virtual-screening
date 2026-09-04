#!/usr/bin/env bash
set -euo pipefail

FPOCKET_VERSION="4.2.3"
SOURCE_DIR="/opt/fpocket-${FPOCKET_VERSION}"

if [[ "${EUID}" -ne 0 ]]; then
  echo "Run this installer as root inside WSL." >&2
  exit 1
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  ca-certificates \
  gcc \
  g++ \
  git \
  make \
  libnetcdf-dev \
  python3

if [[ ! -d "${SOURCE_DIR}/.git" ]]; then
  git clone --branch "${FPOCKET_VERSION}" --depth 1 \
    https://github.com/Discngine/fpocket.git "${SOURCE_DIR}"
fi

installed_commit="$(git -C "${SOURCE_DIR}" rev-parse HEAD)"
expected_commit="$(git -C "${SOURCE_DIR}" rev-list -n 1 "${FPOCKET_VERSION}")"
if [[ "${installed_commit}" != "${expected_commit}" ]]; then
  echo "Existing ${SOURCE_DIR} is not fpocket ${FPOCKET_VERSION}; refusing to overwrite it." >&2
  exit 1
fi

# fpocket 4.2.3's bundled Qhull rules are not parallel-safe on a clean tree:
# qvoronoi/qconvex can start before generated/copied headers are available.
make -C "${SOURCE_DIR}"
make -C "${SOURCE_DIR}" install

command -v fpocket
fpocket -h | head -n 5 || true
echo "FPOCKET_WSL_INSTALL_OK version=${FPOCKET_VERSION} commit=${installed_commit}"
