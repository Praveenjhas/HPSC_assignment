param(
    [switch]$Clean
)

$ErrorActionPreference = "Stop"

function Step($message) {
    Write-Host ""
    Write-Host ("=" * 72) -ForegroundColor DarkCyan
    Write-Host $message -ForegroundColor Cyan
    Write-Host ("=" * 72) -ForegroundColor DarkCyan
}

function Run-Checked($command, $workingDir = $PWD.Path) {
    Push-Location $workingDir
    try {
        Write-Host ">> $command" -ForegroundColor Yellow
        Invoke-Expression $command
        if ($LASTEXITCODE -ne 0) {
            throw "Command failed with exit code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }
}

function Get-ToolPath($name) {
    $tool = Get-Command $name -ErrorAction SilentlyContinue
    if ($null -eq $tool) {
        return $null
    }
    return $tool.Source
}

function Print-MissingToolHelp($name) {
    Write-Host ""
    Write-Host "Missing required tool: $name" -ForegroundColor Red
    switch ($name) {
        "clang++" {
            Write-Host "Install LLVM/Clang, for example with winget:" -ForegroundColor Yellow
            Write-Host "  winget install LLVM.LLVM" -ForegroundColor Yellow
        }
        "python" {
            Write-Host "Install Python 3, for example with winget:" -ForegroundColor Yellow
            Write-Host "  winget install Python.Python.3.13" -ForegroundColor Yellow
        }
        "pip" {
            Write-Host "Install Python with pip enabled. On most systems pip comes with Python." -ForegroundColor Yellow
        }
        "pdflatex" {
            Write-Host "Install a LaTeX distribution with pdflatex, for example MiKTeX:" -ForegroundColor Yellow
            Write-Host "  winget install MiKTeX.MiKTeX" -ForegroundColor Yellow
        }
    }
}

function Assert-Tool($name) {
    $path = Get-ToolPath $name
    if ($null -eq $path -or $path -eq "") {
        Print-MissingToolHelp $name
        throw "Cannot continue without $name"
    }
    Write-Host ("Found {0}: {1}" -f $name, $path) -ForegroundColor Green
}

Step "Checking required tools"
Assert-Tool "clang++"
Assert-Tool "python"
Assert-Tool "pip"
Assert-Tool "pdflatex"

if ($Clean) {
    Step "Cleaning previous generated outputs"
    if (Test-Path "results") {
        Remove-Item -Recurse -Force "results"
    }
    if (Test-Path "report\\figures") {
        Remove-Item -Recurse -Force "report\\figures"
    }
    New-Item -ItemType Directory -Force "report\\figures" | Out-Null
}

Step "Building the solver"
Run-Checked "clang++ -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp src/main.cpp -o dem_solver -fopenmp"

Step "Running verification cases"
Run-Checked ".\\dem_solver verification results"

Step "Running individual base cases"
Run-Checked ".\\dem_solver free_fall results"
Run-Checked ".\\dem_solver constant_velocity results"
Run-Checked ".\\dem_solver bounce results"

Step "Running main particle-count experiment"
Run-Checked ".\\dem_solver experiment results"

Step "Running profiling and OpenMP scaling cases"
Run-Checked ".\\dem_solver scaling results"

Step "Running Section 18 neighbor-search bonus"
Run-Checked ".\\dem_solver neighbor_bonus results"

Step "Running Section 19 scientific investigation bonus"
Run-Checked ".\\dem_solver science_bonus results"

Step "Installing Python plotting dependencies if needed"
Run-Checked "python -m pip install -r requirements.txt"

Step "Generating report figures"
Run-Checked "python scripts\\plot_results.py"

Step "Compiling the report PDF"
Run-Checked "pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex" "report"
Run-Checked "pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex" "report"

Step "Pipeline completed"
Write-Host "Final report:" -ForegroundColor Green
Write-Host "report\\white_paper.pdf" -ForegroundColor Green
Write-Host ""
Write-Host "Generated data:" -ForegroundColor Green
Write-Host "results\\" -ForegroundColor Green
Write-Host ""
Write-Host "One-command workflow succeeded." -ForegroundColor Green
