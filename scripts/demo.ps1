param(
    [ValidateSet("8YEZ", "8ZU3", "8YFC", "9VMX")]
    [string]$PdbId = "8YEZ",
    [string]$Python = "",
    [switch]$Pause
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
if (-not $Python) {
    $PythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if (-not $PythonCommand) { $PythonCommand = Get-Command py -ErrorAction SilentlyContinue }
    if (-not $PythonCommand) { throw "Python was not found on PATH. Pass -Python with an executable path." }
    $Python = $PythonCommand.Source
}

Set-Location -LiteralPath $ProjectRoot
$Host.UI.RawUI.WindowTitle = "PIEZO1 Virtual Screening Demo"

Write-Host "PIEZO1 single-conformation pocket-screening demo" -ForegroundColor Cyan
Write-Host "Project: $ProjectRoot"
$Rules = Get-Content -Raw -LiteralPath config\consensus_rules.json | ConvertFrom-Json
Write-Host "Active structures: $($Rules.active_pdb_ids -join ', ')"
Write-Host "Selected independent conformation: $PdbId" -ForegroundColor Green
Write-Host ""

Write-Host "[1/4] Verify cached structures and completed results" -ForegroundColor Yellow
& $Python scripts\08_verify_mvp.py --pdb-id $PdbId
if ($LASTEXITCODE -ne 0) { throw "MVP verification failed" }

Write-Host ""
Write-Host "[2/4] Re-run P2Rank 2.5.1 on $PdbId" -ForegroundColor Yellow
& $Python scripts\04_run_p2rank.py --pdb-id $PdbId --threads 8
if ($LASTEXITCODE -ne 0) { throw "P2Rank failed" }

Write-Host ""
Write-Host "[3/4] Rebuild the current two-tool precheck on $PdbId" -ForegroundColor Yellow
& $Python scripts\05_match_pockets.py --pdb-id $PdbId --top-p2rank 20
if ($LASTEXITCODE -ne 0) { throw "Result integration failed" }

Write-Host ""
Write-Host "[4/4] Key results" -ForegroundColor Yellow
$Matches = Import-Csv -LiteralPath "results\pocket_matching\${PdbId}_p2rank_dogsite3_matches.csv"
$Supported = @($Matches | Where-Object { $_.two_tool_candidate -eq "True" }).Count
$P2RankPredictions = Import-Csv -LiteralPath "runs\p2rank\${PdbId}\output\${PdbId}.pdb_predictions.csv"
$LatestJob = Get-Content -Raw -LiteralPath "runs\dogsite3\${PdbId}\latest_job.json" | ConvertFrom-Json
$DogSiteResult = Get-Content -Raw -LiteralPath "runs\dogsite3\${PdbId}\$($LatestJob.job_id)\result.json" | ConvertFrom-Json

Write-Host "Selected conformation: $PdbId"
Write-Host "P2Rank pockets: $($P2RankPredictions.Count)"
Write-Host "DoGSite3 pockets: $($DogSiteResult.residues.Count)"
Write-Host "Preliminary pairwise support among P2Rank top 20: $Supported/20"
$FpocketLatestPath = "runs\fpocket\${PdbId}\latest_run.json"
if (Test-Path -LiteralPath $FpocketLatestPath) {
    & $Python scripts\09_build_consensus.py --pdb-id $PdbId
    if ($LASTEXITCODE -ne 0) { throw "Three-tool consensus generation failed" }
    $FpocketLatest = Get-Content -Raw -LiteralPath $FpocketLatestPath | ConvertFrom-Json
    $Consensus = Import-Csv -LiteralPath "results\consensus\${PdbId}_consensus_pockets.csv"
    $Tier1 = @($Consensus | Where-Object { $_.support_count -eq "2" }).Count
    $Tier2 = @($Consensus | Where-Object { $_.support_count -eq "3" }).Count
    Write-Host "fpocket pockets: $($FpocketLatest.pocket_count)"
    Write-Host "Tools completed: 3/3 (P2Rank + DoGSite3 + fpocket)" -ForegroundColor Green
    Write-Host "Final consensus regions: $($Consensus.Count)"
    Write-Host "T1 (support_count=2): $Tier1"
    Write-Host "T2 (support_count=3): $Tier2" -ForegroundColor Green
} else {
    Write-Host "Tools completed: 2/3 (P2Rank + DoGSite3)" -ForegroundColor Yellow
    Write-Host "Third tool fpocket: PENDING" -ForegroundColor Red
    Write-Host "Final T1/T2 classification: NOT GENERATED" -ForegroundColor Red
}
Write-Host ""
Write-Host "Correct pipeline:" -ForegroundColor Green
Write-Host "$PdbId -> P2Rank / DoGSite3 / fpocket -> normalize pockets -> count tool support"
Write-Host "T1 = exactly $($Rules.tiers.T1.support_count) of 3 tools support the region"
Write-Host "T2 = all $($Rules.tiers.T2.support_count) tools support the region"
Write-Host ""
Write-Host "The old 8YEZ-vs-8ZU3 comparison has been removed from this demo." -ForegroundColor Red
Write-Host "Each of 8YEZ, 8ZU3, 8YFC and 9VMX must be analyzed independently." -ForegroundColor Yellow
if (Test-Path -LiteralPath $FpocketLatestPath) {
    Write-Host "This structure now has final computational T1/T2 pocket candidates." -ForegroundColor Green
    Write-Host "These are predicted pockets, not experimentally validated binding sites or drug candidates." -ForegroundColor Yellow
} else {
    Write-Host "Current two-tool matching is only a progress check; it is not yet T1." -ForegroundColor Yellow
    Write-Host "Tier labels can be assigned only after fpocket finishes." -ForegroundColor Yellow
}
Write-Host "Demo completed successfully. No drug candidate is claimed." -ForegroundColor Green
if ($Pause) {
    Write-Host "Press Enter to close this window."
    [void](Read-Host)
}
