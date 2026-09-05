param(
    [Parameter(Mandatory = $true)][string]$PdbId,
    [Parameter(Mandatory = $true)][string]$ConsensusId,
    [Parameter(Mandatory = $true)][string]$Compounds,
    [string]$Checkpoint,
    [string]$DistroName = "Ubuntu-fpocket",
    [int]$GpuId = 0
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$PdbId = $PdbId.ToUpperInvariant()
$ConsensusId = $ConsensusId.ToUpperInvariant()

function ConvertTo-WslPath([string]$WindowsPath) {
    $fullPath = [IO.Path]::GetFullPath($WindowsPath)
    if ($fullPath -notmatch '^([A-Za-z]):\\(.*)$') {
        throw "Expected an absolute Windows drive path, got: $fullPath"
    }
    $drive = $Matches[1].ToLowerInvariant()
    $tail = $Matches[2].Replace('\', '/')
    return "/mnt/$drive/$tail"
}

$CompoundsPath = (Resolve-Path -LiteralPath $Compounds).Path
if (-not $Checkpoint) {
    $Checkpoint = Join-Path $RepoRoot "tools\drugclip\official\checkpoint_best.pt"
}
$CheckpointPath = (Resolve-Path -LiteralPath $Checkpoint).Path
$RunId = [DateTime]::UtcNow.ToString("yyyyMMddTHHmmssZ")
$InputDir = Join-Path $RepoRoot "runs\drugclip_inputs\$PdbId\$ConsensusId\$RunId"
$ResultDir = Join-Path $RepoRoot "runs\drugclip\$PdbId\$ConsensusId\$RunId"

$RepoRootWsl = ConvertTo-WslPath $RepoRoot
$CompoundsWsl = ConvertTo-WslPath $CompoundsPath
$InputDirWsl = ConvertTo-WslPath $InputDir
$PrepareArguments = @(
    "-d", $DistroName, "-u", "root", "--cd", $RepoRootWsl, "--",
    "/opt/drugclip-venv/bin/python", "scripts/10_prepare_drugclip_inputs.py",
    "--pdb-id", $PdbId, "--consensus-id", $ConsensusId,
    "--compounds", $CompoundsWsl, "--output-dir", $InputDirWsl
)
& wsl.exe @PrepareArguments
if ($LASTEXITCODE -ne 0) { throw "DrugCLIP input preparation failed." }

$InputMetadata = Get-Content -LiteralPath (Join-Path $InputDir "input_run.json") -Raw | ConvertFrom-Json
$Tier = [string]$InputMetadata.tier
$RunParameters = @{
    PdbId = $PdbId
    ConsensusId = $ConsensusId
    Tier = $Tier
    MolLmdb = (Join-Path $InputDir "mols.lmdb")
    PocketLmdb = (Join-Path $InputDir "pocket.lmdb")
    Checkpoint = $CheckpointPath
    OutputDir = $ResultDir
    DistroName = $DistroName
    GpuId = $GpuId
}
& (Join-Path $PSScriptRoot "run_drugclip_wsl.ps1") @RunParameters

$RankedWsl = ConvertTo-WslPath (Join-Path $ResultDir "embeddings\ranked_compounds.txt")
$ManifestWsl = ConvertTo-WslPath (Join-Path $InputDir "molecule_manifest.csv")
$OutputCsv = Join-Path $ResultDir "ranked_compounds.csv"
$OutputCsvWsl = ConvertTo-WslPath $OutputCsv
$NormalizeArguments = @(
    "-d", $DistroName, "-u", "root", "--cd", $RepoRootWsl, "--",
    "/opt/drugclip-venv/bin/python", "scripts/12_normalize_drugclip_results.py",
    "--ranked", $RankedWsl, "--manifest", $ManifestWsl,
    "--output", $OutputCsvWsl, "--pdb-id", $PdbId,
    "--consensus-id", $ConsensusId, "--tier", $Tier
)
& wsl.exe @NormalizeArguments
if ($LASTEXITCODE -ne 0) { throw "DrugCLIP result normalization failed." }

Write-Host "DRUGCLIP_PIPELINE_OK pdb_id=$PdbId consensus_id=$ConsensusId tier=$Tier"
Write-Host "Input provenance: $InputDir"
Write-Host "Ranked CSV: $OutputCsv"

