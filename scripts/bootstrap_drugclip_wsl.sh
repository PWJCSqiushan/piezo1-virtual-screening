#!/usr/bin/env bash
set -euo pipefail

ENV_DIR="/opt/drugclip-venv"
UV_DIR="/opt/uv"
UNICORE_DIR="/opt/Uni-Core-0.0.3"
UNICORE_COMMIT="44f6386f4dcd7137fc1e5d5e768117d635d64a26"
PYTORCH_INDEX="https://download.pytorch.org/whl/cu128"
PYPI_INDEX="https://pypi.tuna.tsinghua.edu.cn/simple"

if [[ "${EUID}" -ne 0 ]]; then
  echo "Run this installer as root inside WSL." >&2
  exit 1
fi

export DEBIAN_FRONTEND=noninteractive
export UV_HTTP_TIMEOUT=300
export UV_HTTP_RETRIES=10
apt-get update
apt-get install -y --no-install-recommends \
  ca-certificates \
  curl \
  git \
  build-essential \
  libxrender1 \
  libxext6

if [[ ! -x "${UV_DIR}/uv" ]]; then
  mkdir -p "${UV_DIR}"
  installer="${UV_DIR}/installer.sh"
  curl --fail --location --retry 5 --retry-delay 2 \
    https://astral.sh/uv/install.sh --output "${installer}"
  UV_INSTALL_DIR="${UV_DIR}" sh "${installer}"
fi

"${UV_DIR}/uv" python install 3.10
if [[ ! -x "${ENV_DIR}/bin/python" ]]; then
  "${UV_DIR}/uv" venv --python 3.10 "${ENV_DIR}"
fi

PYTHON="${ENV_DIR}/bin/python"
"${UV_DIR}/uv" pip install --python "${PYTHON}" --upgrade \
  --default-index "${PYPI_INDEX}" \
  pip==25.2 \
  setuptools==69.5.1 \
  wheel==0.45.1
if [[ -n "${DRUGCLIP_WHEEL_DIR:-}" && -d "${DRUGCLIP_WHEEL_DIR}" ]]; then
  "${UV_DIR}/uv" pip install --python "${PYTHON}" \
    --no-index --find-links "${DRUGCLIP_WHEEL_DIR}" \
    torch==2.7.1+cu128 \
    torchvision==0.22.1+cu128 \
    torchaudio==2.7.1+cu128
else
  "${UV_DIR}/uv" pip install --python "${PYTHON}" \
    torch==2.7.1 \
    torchvision==0.22.1 \
    torchaudio==2.7.1 \
    --index-url "${PYTORCH_INDEX}"
fi
"${UV_DIR}/uv" pip install --python "${PYTHON}" \
  --default-index "${PYPI_INDEX}" \
  numpy==1.23.5 \
  scipy==1.10.1 \
  pandas==1.5.3 \
  scikit-learn==1.2.2 \
  rdkit-pypi==2022.9.5 \
  lmdb==1.4.1 \
  tqdm==4.65.0 \
  tensorboardX==2.6.2.2 \
  tokenizers==0.13.3 \
  ml_collections==0.1.1 \
  iopath==0.1.10 \
  wandb==0.15.12 \
  biopandas==0.4.1 \
  ipython==8.18.1
"${UV_DIR}/uv" pip install --python "${PYTHON}" --default-index "${PYPI_INDEX}" \
  setuptools==69.5.1

if [[ -n "${DRUGCLIP_UNICORE_SOURCE_DIR:-}" && -d "${DRUGCLIP_UNICORE_SOURCE_DIR}/.git" ]]; then
  source_commit="$(git -C "${DRUGCLIP_UNICORE_SOURCE_DIR}" rev-parse HEAD)"
  if [[ "${source_commit}" != "${UNICORE_COMMIT}" ]]; then
    echo "Local Uni-Core source commit mismatch: ${source_commit}" >&2
    exit 1
  fi
  target_commit=""
  if [[ -d "${UNICORE_DIR}/.git" ]]; then
    target_commit="$(git -C "${UNICORE_DIR}" rev-parse HEAD 2>/dev/null || true)"
  fi
  if [[ "${target_commit}" != "${UNICORE_COMMIT}" ]]; then
    if [[ -e "${UNICORE_DIR}" ]]; then
      mv "${UNICORE_DIR}" "${UNICORE_DIR}.incomplete.$(date -u +%Y%m%dT%H%M%SZ)"
    fi
    cp -a "${DRUGCLIP_UNICORE_SOURCE_DIR}" "${UNICORE_DIR}"
  fi
else
  if [[ ! -d "${UNICORE_DIR}/.git" ]]; then
    git clone https://github.com/dptech-corp/Uni-Core.git "${UNICORE_DIR}"
  fi
  git -C "${UNICORE_DIR}" fetch --depth 1 origin "${UNICORE_COMMIT}"
  git -C "${UNICORE_DIR}" checkout --detach "${UNICORE_COMMIT}"
fi
actual_commit="$(git -C "${UNICORE_DIR}" rev-parse HEAD)"
if [[ "${actual_commit}" != "${UNICORE_COMMIT}" ]]; then
  echo "Uni-Core commit mismatch: ${actual_commit}" >&2
  exit 1
fi

(
  cd "${UNICORE_DIR}"
  "${PYTHON}" setup.py install --disable-cuda-ext
)

"${PYTHON}" - <<'PY'
import lmdb
import numpy
import rdkit
import sklearn
import torch
import unicore

assert torch.cuda.is_available(), "PyTorch cannot see the WSL CUDA GPU"
capability = torch.cuda.get_device_capability(0)
x = torch.tensor([1.0], device="cuda")
print(
    "DRUGCLIP_WSL_INSTALL_OK",
    f"python={__import__('sys').version.split()[0]}",
    f"torch={torch.__version__}",
    f"cuda={torch.version.cuda}",
    f"gpu={torch.cuda.get_device_name(0)}",
    f"capability={capability}",
    f"rdkit={rdkit.__version__}",
    f"unicore={getattr(unicore, '__version__', 'unknown')}",
    f"cuda_test={x.item()}",
)
PY
