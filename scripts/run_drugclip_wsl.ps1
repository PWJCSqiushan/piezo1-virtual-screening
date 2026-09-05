param(
    [Parameter(Mandatory = $true)][string]$PdbId,
    [Parameter(Mandatory = $true)][string]$ConsensusId,
    [Parameter(Mandatory = $true)][ValidateSet("T1", "T2")][string]$Tier,
    [Parameter(Mandatory = $true)][string]$MolLmdb,
    [Parameter(Mandatory = $true)][string]$PocketLmdb,
    [Parameter(Mandatory = $true)][string]$Checkpoint,
    [Parameter(Mandatory = $true)][string]$OutputDir,
    [string]$DistroName = "Ubuntu-fpocket",
    [int]$GpuId = 0
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

function ConvertTo-WslPath([string]$WindowsPath) {
    $fullPath = [IO.Path]::GetFullPath($WindowsPath)
    if ($fullPath -notmatch '^([A-Za-z]):\\(.*)$') {
        throw "Expected an absolute Windows drive path, got: $fullPath"
    }
    $drive = $Matches[1].ToLowerInvariant()
    $tail = $Matches[2].Replace('\', '/')
    return "/mnt/$drive/$tail"
}

foreach ($InputPath in @($MolLmdb, $PocketLmdb, $Checkpoint)) {
    if (-not (Test-Path -LiteralPath $InputPath -PathType Leaf)) {
        throw "Missing DrugCLIP input: $InputPath"
    }
}

$RepoRootWsl = ConvertTo-WslPath $RepoRoot
$MolWsl = ConvertTo-WslPath $MolLmdb
$PocketWsl = ConvertTo-WslPath $PocketLmdb
$CheckpointWsl = ConvertTo-WslPath $Checkpoint
$OutputWsl = ConvertTo-WslPath $OutputDir

$WslArguments = @(
    "-d", $DistroName, "-u", "root", "--cd", $RepoRootWsl, "--",
    "/opt/drugclip-venv/bin/python", "scripts/11_run_drugclip.py",
    "--pdb-id", $PdbId, "--consensus-id", $ConsensusId, "--tier", $Tier,
    "--mol-lmdb", $MolWsl, "--pocket-lmdb", $PocketWsl,
    "--checkpoint", $CheckpointWsl, "--output-dir", $OutputWsl,
    "--gpu-id", $GpuId
)
& wsl.exe @WslArguments
if ($LASTEXITCODE -ne 0) { throw "DrugCLIP retrieval failed. Inspect stdout.log/stderr.log in $OutputDir." }
