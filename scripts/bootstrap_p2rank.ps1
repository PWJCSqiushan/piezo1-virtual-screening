param(
    [string]$Version = "2.5.1",
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
& $Python (Join-Path $PSScriptRoot "bootstrap_p2rank.py") --version $Version
if ($LASTEXITCODE -ne 0) { throw "P2Rank bootstrap failed" }
