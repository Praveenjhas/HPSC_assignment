# HPSC Assignment 1 Submission

This repository contains the final submission package for HPSC 2026 Assignment 1.

The project implements a compact three-dimensional discrete element method (DEM) solver for spherical particles in a rectangular box. The work includes the serial baseline solver, verification studies, runtime profiling, OpenMP parallelisation of the particle-contact loop, a neighbour-search bonus extension, and scientific bonus studies on damping and particle-cloud settling.

Included work:

- serial 3D DEM solver for spherical particles
- OpenMP parallel particle-contact loop
- verification cases
- profiling and thread-scaling runs
- Section 18 bonus: neighbour search
- Section 19 bonus: damping and timestep studies

This submission covers the assignment through Section 19 only.

## Platform note

The provided automation is written for Windows PowerShell. The project can still be built and run manually on other systems if equivalent tools are available, but the packaged one-command workflow is intended for Windows.

## Main command

On Windows, the recommended command is:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1
```

For a fresh rerun that deletes old generated outputs first:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -Clean
```

This command does the full workflow:

1. checks required tools
2. builds the solver
3. runs the individual base cases (`free_fall`, `constant_velocity`, `bounce`)
4. runs the grouped simulation cases
5. installs Python plotting packages
6. generates all figures
7. compiles the final report PDF

During long runs, the solver also prints live progress lines in the terminal showing:

- current step
- simulation time
- kinetic energy
- center-of-mass height
- active contacts
- maximum speed

Final outputs:

- `results/`
- `report/figures/`
- `report/white_paper.pdf`

The `results/` directory stores CSV diagnostics, text summaries, and selected particle snapshots for each run. The plotting script reads these outputs and writes publication figures into `report/figures/`, and the LaTeX report then includes those figures when compiling `report/white_paper.pdf`.

## Required tools

The one-command workflow expects these tools to already be installed:

- `clang++`
- `python`
- `pip`
- `pdflatex`

Python plotting packages are installed automatically from `requirements.txt`.

If one of the required tools is missing, `run_all.ps1` stops early and prints what needs to be installed.

## Optional make shortcuts

If `make` is installed, these shortcuts are also available:

```powershell
make full
make run-free-fall
make run-constant-velocity
make run-bounce
make run-verify
make run-experiment
make run-scaling
make run-neighbor
make run-science
```

`make full` is only a shortcut for:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1
```

Most users should use `run_all.ps1` directly.

## What each solver mode does

- `free_fall` runs the single-particle analytical verification under gravity.
- `constant_velocity` runs the zero-force constant-velocity verification.
- `bounce` runs the single-particle wall-bounce case.
- `verification` runs the grouped verification set, including timestep sensitivity and damping variants.
- `experiment` runs the increasing-particle-count serial study.
- `scaling` runs profiling, OpenMP strong scaling, weak scaling, and serial-versus-parallel correctness checks.
- `neighbor_bonus` runs the brute-force versus neighbour-grid comparison.
- `science_bonus` runs the damping studies and particle-cloud settling studies.

## Manual commands

If you want to run steps manually, use these commands.

Build:

```powershell
clang++ -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp src/main.cpp -o dem_solver -fopenmp
```

Run solver modes:

```powershell
.\dem_solver free_fall results
.\dem_solver constant_velocity results
.\dem_solver bounce results
.\dem_solver verification results
.\dem_solver experiment results
.\dem_solver scaling results
.\dem_solver neighbor_bonus results
.\dem_solver science_bonus results
```

Generate figures:

```powershell
python -m pip install -r requirements.txt
python scripts\plot_results.py
```

Compile report:

```powershell
cd report
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
```

## Project files

- `src/main.cpp` : main DEM solver, experiment drivers, and output logic
- `run_all.ps1` : recommended end-to-end Windows workflow
- `Makefile` : optional command shortcuts
- `requirements.txt` : Python plotting dependencies
- `scripts/plot_results.py` : figure generation from CSV outputs
- `results/` : generated diagnostics, summaries, and particle snapshots
- `report/figures/` : generated report figures
- `report/white_paper.tex` : report source
- `report/white_paper.pdf` : compiled final report

## Notes

- The code is written in C++17.
- OpenMP is used in the brute-force particle-contact loop.
- The neighbour-search implementation is used for the Section 18 bonus study.
- The scientific bonus section includes damping and particle-cloud settling studies.
- The report includes a short LLM-assistance disclosure, as required by the assignment instructions.
