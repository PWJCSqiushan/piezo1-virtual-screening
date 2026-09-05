param(
    [Parameter(Mandatory = $true)][string]$PdbId,
    [Parameter(Mandatory = $true)][string]$ConsensusId,
    [Parameter(Mandatory = $true)][string]$RankedCsv,
    [ValidateRange(1, 100)][int]$TopN = 20,
    [ValidateRange(1, 256)][int]$Exhaustiveness = 8,
    [ValidateRange(1, 20)][int]$NumModes = 9,
    [string]$DistroName = "Ubuntu-fpocket",
    [int]$GpuId = 0,
    [switch]$NoGpu,
    [switch]$AllowNonPiezo1Scope
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$PdbId = $PdbId.ToUpperInvariant()
$ConsensusId = $ConsensusId.ToUpperInvariant()

function ConvertTo-WslPath([string]$WindowsPath) {
    $fullPath = [IO.Path]::GetFullPath($WindowsPath)
    if ($fullPath -notmatch '^([A-Za-z]):\\(.*)$') { throw "Expected an absolute Windows drive path: $fullPath" }
    return "/mnt/$($Matches[1].ToLowerInvariant())/$($Matches[2].Replace('\', '/'))"
}

$RankedPath = (Resolve-Path -LiteralPath $RankedCsv).Path
$GninaPath = Join-Path $RepoRoot "tools\gnina\gnina.cuda12.8.static"
if (-not (Test-Path -LiteralPath $GninaPath -PathType Leaf)) {
    throw "Pinned GNINA binary is missing. Run scripts\bootstrap_gnina_wsl.ps1 first."
}
$RunId = [DateTime]::UtcNow.ToString("yyyyMMddTHHmmssZ")
$InputDir = Join-Path $RepoRoot "runs\gnina_inputs\$PdbId\$ConsensusId\$RunId"
$OutputDir = Join-Path $RepoRoot "runs\gnina\$PdbId\$ConsensusId\$RunId"
$RepoWsl = ConvertTo-WslPath $RepoRoot
$InputWsl = ConvertTo-WslPath $InputDir
$OutputWsl = ConvertTo-WslPath $OutputDir
$RankedWsl = ConvertTo-WslPath $RankedPath
$GninaWsl = ConvertTo-WslPath $GninaPath

$Prepare = @(
    "-d", $DistroName, "-u", "root", "--cd", $RepoWsl, "--",
    "/opt/drugclip-venv/bin/python", "scripts/15_prepare_gnina_inputs.py",
    "--pdb-id", $PdbId, "--consensus-id", $ConsensusId,
    "--ranked-csv", $RankedWsl, "--top-n", $TopN,
    "--output-dir", $InputWsl
)
if ($AllowNonPiezo1Scope) { $Prepare += "--allow-non-piezo1-scope" }
& wsl.exe @Prepare
if ($LASTEXITCODE -ne 0) { throw "GNINA input preparation failed." }

$Run = @(
    "-d", $DistroName, "-u", "root", "--cd", $RepoWsl, "--",
    "/opt/drugclip-venv/bin/python", "scripts/16_run_gnina.py",
    "--gnina", $GninaWsl, "--input-dir", $InputWsl, "--output-dir", $OutputWsl,
    "--exhaustiveness", $Exhaustiveness, "--num-modes", $NumModes, "--gpu-id", $GpuId
)
if ($NoGpu) { $Run += "--no-gpu" }
& wsl.exe @Run
if ($LASTEXITCODE -ne 0) { throw "GNINA docking failed. Inspect $OutputDir." }

Write-Host "GNINA_PIPELINE_OK pdb_id=$PdbId consensus_id=$ConsensusId" -ForegroundColor Green
Write-Host "Input evidence: $InputDir"
Write-Host "Docked SDF: $(Join-Path $OutputDir 'docked.sdf')"
Write-Host "Best poses CSV: $(Join-Path $OutputDir 'best_poses.csv')"
