<#!
.SYNOPSIS
    Launches the WNK umbrella sampling demo in synthetic mode on Windows.
.DESCRIPTION
    Creates the Conda environment from environment.yml if it is missing and runs
    `python -m Chronosfold.umbrella_suite.run_wnk_pipeline --synthetic` inside it.
#>

param(
    [string]$EnvName = "bsm-lancad-env",
    [string]$CondaPath
)

$ErrorActionPreference = "Stop"

function Get-CondaCommand {
    param([string]$Override)

    if ($Override) {
        $resolved = Resolve-Path -LiteralPath $Override -ErrorAction SilentlyContinue
        if ($resolved) {
            return $resolved.Path
        }
        throw "Conda executable not found at $Override"
    }

    if ($env:CONDA_EXE -and (Test-Path $env:CONDA_EXE)) {
        return $env:CONDA_EXE
    }

    $command = Get-Command conda -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    $candidates = @(
        (Join-Path $env:USERPROFILE "miniconda3\Scripts\conda.exe"),
        (Join-Path $env:USERPROFILE "miniconda3\condabin\conda.bat"),
        (Join-Path $env:USERPROFILE "miniconda3\Scripts\conda.bat"),
        (Join-Path $env:USERPROFILE "Anaconda3\Scripts\conda.exe"),
        (Join-Path $env:USERPROFILE "anaconda3\Scripts\conda.exe")
    )

    foreach ($candidate in $candidates) {
        if (Test-Path $candidate) {
            return $candidate
        }
    }

    throw "Conda command not available. Install Miniconda/Anaconda or provide --CondaPath."
}

function Ensure-CondaEnv {
    param(
        [string]$Name,
        [string]$EnvironmentFile,
        [string]$CondaCommand
    )

    $envList = & $CondaCommand env list | Select-String -Pattern "^\s*$Name(\s|$)"
    if (-not $envList) {
        Write-Host "Creating Conda environment '$Name' from $EnvironmentFile" -ForegroundColor Cyan
        & $CondaCommand env create -f $EnvironmentFile
    }
}

$scriptRoot = Split-Path -Path $MyInvocation.MyCommand.Definition -Parent
$repoRoot = Resolve-Path (Join-Path $scriptRoot "..")
$environmentFile = Join-Path $repoRoot "environment.yml"

$condaCommand = Get-CondaCommand -Override $CondaPath
Write-Host "Using Conda at $condaCommand" -ForegroundColor Yellow

Ensure-CondaEnv -Name $EnvName -EnvironmentFile $environmentFile -CondaCommand $condaCommand

Write-Host "Running umbrella sampling demo (synthetic mode)" -ForegroundColor Green
& $condaCommand run -n $EnvName python -m Chronosfold.umbrella_suite.run_wnk_pipeline --synthetic
