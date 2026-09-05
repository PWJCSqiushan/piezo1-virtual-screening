param(
    [string]$DistroName = "Ubuntu-fpocket"
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

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw "WSL is not installed. Run the fpocket WSL bootstrap first."
}
$InstalledDistros = @(& wsl.exe --list --quiet 2>$null) -replace "`0", "" | ForEach-Object { $_.Trim() }
if ($DistroName -notin $InstalledDistros) {
    throw "WSL distribution '$DistroName' is missing. Run bootstrap_fpocket_wsl.ps1 first."
}

$SetupScriptWsl = ConvertTo-WslPath (Join-Path $PSScriptRoot "bootstrap_drugclip_wsl.sh")
$WheelDir = Join-Path $RepoRoot "tools\drugclip\wheels\cu128-py310"
$WheelDirWsl = ConvertTo-WslPath $WheelDir
$UniCoreSource = Join-Path $RepoRoot "tools\Uni-Core"
$UniCoreSourceWsl = ConvertTo-WslPath $UniCoreSource
Write-Host "Installing the isolated DrugCLIP environment inside $DistroName..."
$WslArguments = @(
    "-d", $DistroName, "-u", "root", "--", "env",
    "DRUGCLIP_WHEEL_DIR=$WheelDirWsl",
    "DRUGCLIP_UNICORE_SOURCE_DIR=$UniCoreSourceWsl",
    "bash", $SetupScriptWsl
)
& wsl.exe @WslArguments
if ($LASTEXITCODE -ne 0) { throw "DrugCLIP WSL environment installation failed." }

Write-Host "DRUGCLIP_BOOTSTRAP_OK distro=$DistroName env=/opt/drugclip-venv"
