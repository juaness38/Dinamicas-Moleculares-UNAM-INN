<#!
.SYNOPSIS
    Prepares the Conda environment and launches the umbrella demo on Windows.
#>

param(
    [string]$EnvName = "bsm-lancad-env",
    [string]$CondaPath
)

$ErrorActionPreference = "Stop"

$scriptRoot = Split-Path -Path $MyInvocation.MyCommand.Definition -Parent
$repoRoot = Resolve-Path (Join-Path $scriptRoot "..")
$runner = Join-Path $repoRoot "scripts\run_umbrella.ps1"

if (-not (Test-Path $runner)) {
    Write-Error "Unable to locate run_umbrella.ps1 at $runner"
}

$arguments = @{ EnvName = $EnvName }
if ($CondaPath) {
    $arguments["CondaPath"] = $CondaPath
}

try {
    & $runner @arguments
} catch {
    Write-Warning $_
    Write-Host "If Conda is not on PATH, provide the executable path, e.g.:" -ForegroundColor Yellow
    Write-Host "  .\scripts\bootstrap_windows.ps1 -CondaPath '$env:USERPROFILE\miniconda3\Scripts\conda.exe'" -ForegroundColor Yellow
    throw
}
