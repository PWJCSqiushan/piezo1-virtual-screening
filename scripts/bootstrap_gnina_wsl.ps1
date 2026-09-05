param(
    [string]$DistroName = "Ubuntu-fpocket",
    [switch]$SkipDownload,
    [switch]$SkipRuntimeInstall
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$ToolConfig = Get-Content -Raw -LiteralPath (Join-Path $RepoRoot "config\tool_sources.json") | ConvertFrom-Json
$Gnina = $ToolConfig.tools.gnina
$ToolDir = Join-Path $RepoRoot "tools\gnina"
$FinalPath = Join-Path $ToolDir "gnina.cuda12.8.static"
$PartPath = "$FinalPath.part"
New-Item -ItemType Directory -Force -Path $ToolDir | Out-Null

function Test-GninaBinary([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    $Item = Get-Item -LiteralPath $Path
    if ($Item.Length -ne [int64]$Gnina.binary_size_bytes) { return $false }
    $Hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
    return $Hash -eq ([string]$Gnina.binary_sha256).ToLowerInvariant()
}

if (-not (Test-GninaBinary $FinalPath)) {
    if ($SkipDownload) { throw "Pinned GNINA binary is missing or invalid: $FinalPath" }
    if (Test-Path -LiteralPath $FinalPath) {
        throw "An invalid final GNINA file exists. Move it aside manually before retrying; this script will not delete it."
    }
    Write-Host "Downloading pinned GNINA $($Gnina.version) binary (about 2.06 GB)..." -ForegroundColor Yellow
    & curl.exe -L --retry 5 --retry-delay 2 --continue-at - --output $PartPath ([string]$Gnina.binary_url)
    if ($LASTEXITCODE -ne 0) { throw "GNINA download failed; the .part file is retained for resume." }
    if (-not (Test-GninaBinary $PartPath)) { throw "Downloaded GNINA .part file failed size or SHA256 verification." }
    Move-Item -LiteralPath $PartPath -Destination $FinalPath
}

$Drive = $FinalPath.Substring(0, 1).ToLowerInvariant()
$Tail = $FinalPath.Substring(3).Replace('\', '/')
$GninaWsl = "/mnt/$Drive/$Tail"
& wsl.exe -d $DistroName -u root -- chmod +x $GninaWsl
if ($LASTEXITCODE -ne 0) { throw "Could not mark GNINA executable inside WSL." }

$MissingLibraries = & wsl.exe -d $DistroName -u root -- bash -lc "ldd '$GninaWsl' | grep 'not found' || true; ldconfig -p | grep -q 'libnvrtc.so.12' || echo 'libnvrtc.so.12 => not found'"
if ($MissingLibraries -and -not $SkipRuntimeInstall) {
    Write-Host "Installing official NVIDIA CUDA 12.8 runtime and cuDNN 9 packages (about 2 GB download)..." -ForegroundColor Yellow
    $RuntimeInstall = @'
set -euo pipefail
if [ ! -f /usr/share/keyrings/cuda-archive-keyring.gpg ]; then
  cd /tmp
  curl -fL --retry 5 -o cuda-keyring_1.1-1_all.deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb
  dpkg -i cuda-keyring_1.1-1_all.deb
fi
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get -y install cudnn9-cuda-12 cuda-cudart-12-8 libcublas-12-8 libcusparse-12-8 libcufft-12-8 libcusolver-12-8 libnvjitlink-12-8 cuda-nvrtc-12-8
ldconfig
'@
    & wsl.exe -d $DistroName -u root -- bash -lc $RuntimeInstall
    if ($LASTEXITCODE -ne 0) { throw "GNINA runtime dependency installation failed." }
} elseif ($MissingLibraries) {
    throw "GNINA runtime libraries are missing and -SkipRuntimeInstall was supplied: $MissingLibraries"
}
& wsl.exe -d $DistroName -u root -- $GninaWsl --version
if ($LASTEXITCODE -ne 0) { throw "GNINA version check failed." }
Write-Host "GNINA_BOOTSTRAP_OK $FinalPath" -ForegroundColor Green
