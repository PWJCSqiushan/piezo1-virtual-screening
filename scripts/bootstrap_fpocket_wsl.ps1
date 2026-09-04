param(
    [string]$DistroName = "Ubuntu-fpocket",
    [string]$PdbId = "8YEZ",
    [switch]$SkipRun
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$WslDir = Join-Path $RepoRoot "tools\wsl"
$Archive = Join-Path $WslDir "ubuntu-24.04.4-wsl-amd64.wsl"
$InstallLocation = Join-Path $WslDir $DistroName
$UbuntuUrl = "https://releases.ubuntu.com/24.04.4/ubuntu-24.04.4-wsl-amd64.wsl"
$UbuntuSha256 = "9b2f7730dc68227dd04a9f3e5eab86ad85caf556b8606ad94f1f29ff5c4fd3f5"

function ConvertTo-WslPath([string]$WindowsPath) {
    $fullPath = [IO.Path]::GetFullPath($WindowsPath)
    if ($fullPath -notmatch '^([A-Za-z]):\\(.*)$') {
        throw "Expected an absolute Windows drive path, got: $fullPath"
    }
    $drive = $Matches[1].ToLowerInvariant()
    $tail = $Matches[2].Replace('\', '/')
    return "/mnt/$drive/$tail"
}

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw "WSL is not installed. Run 'wsl --install --no-distribution' as Administrator and restart Windows."
}

New-Item -ItemType Directory -Force -Path $WslDir | Out-Null
if (-not (Test-Path -LiteralPath $Archive)) {
    Write-Host "Downloading Ubuntu 24.04.4 WSL image (about 373 MB)..."
    & curl.exe -L --fail --retry 5 --retry-delay 2 -o $Archive $UbuntuUrl
    if ($LASTEXITCODE -ne 0) { throw "Ubuntu image download failed." }
}

$ActualSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $Archive).Hash.ToLowerInvariant()
if ($ActualSha256 -ne $UbuntuSha256) {
    throw "Ubuntu image SHA256 mismatch. Expected $UbuntuSha256, got $ActualSha256."
}
Write-Host "Ubuntu image checksum verified."

$InstalledDistros = @(& wsl.exe --list --quiet 2>$null) -replace "`0", "" | ForEach-Object { $_.Trim() }
if ($DistroName -notin $InstalledDistros) {
    Write-Host "Installing isolated WSL distribution '$DistroName'..."
    & wsl.exe --install --from-file $Archive --name $DistroName --location $InstallLocation --version 2 --no-launch
    if ($LASTEXITCODE -ne 0) {
        throw "WSL distribution installation failed. If Windows optional features were just enabled, restart Windows and rerun this script."
    }
}

$SetupScript = Join-Path $PSScriptRoot "bootstrap_fpocket_wsl.sh"
$SetupScriptWsl = ConvertTo-WslPath $SetupScript
$RepoRootWsl = ConvertTo-WslPath $RepoRoot

Write-Host "Installing fpocket 4.2.3 and its build dependencies inside $DistroName..."
& wsl.exe -d $DistroName -u root -- bash $SetupScriptWsl
if ($LASTEXITCODE -ne 0) { throw "fpocket compilation or installation failed." }

if (-not $SkipRun) {
    $PdbId = $PdbId.ToUpperInvariant()
    Write-Host "Running fpocket for $PdbId..."
    & wsl.exe -d $DistroName -u root --cd $RepoRootWsl -- python3 scripts/run_fpocket.py --pdb-id $PdbId
    if ($LASTEXITCODE -ne 0) { throw "fpocket run failed for $PdbId." }

    & wsl.exe -d $DistroName -u root --cd $RepoRootWsl -- python3 scripts/09_build_consensus.py --pdb-id $PdbId
    if ($LASTEXITCODE -ne 0) { throw "Three-tool consensus generation failed for $PdbId." }
}

Write-Host "FPOCKET_PIPELINE_OK distro=$DistroName pdb_id=$PdbId"
