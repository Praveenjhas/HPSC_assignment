# HPSC Assignment 1 Submission

This repository contains my submission for HPSC 2026 Assignment 1.

The project implements a compact three-dimensional discrete element method (DEM) solver for spherical particles in a rectangular box. The submitted work includes:

- a serial C++17 baseline solver
- verification cases
- runtime profiling
- OpenMP parallelisation of the particle-contact loop
- Section 18 bonus: neighbour search
- Section 19 bonus: damping and particle-cloud settling studies

This submission covers the assignment through the second bonus section. The larger research problem was not attempted.

## Repository layout

- `src/main.cpp`: solver, experiment drivers, output writing, profiling, OpenMP, and bonus-study logic
- `scripts/plot_results.py`: reads generated CSV/text outputs and produces report figures
- `requirements.txt`: Python packages needed for plotting
- `report/white_paper.tex`: LaTeX source for the report
- `report/figures/`: generated figures used in the report
- `results/`: generated simulation outputs, summaries, diagnostics, and selected snapshots
- `run_all.ps1`: end-to-end Windows workflow
- `Makefile`: optional shortcuts for common tasks

## What gets produced

Running the workflow generates three main outputs:

- `results/`: per-case diagnostics, summaries, and selected particle snapshots
- `report/figures/`: plots created from the generated results
- `report/white_paper.pdf`: final compiled report

The report figures are not hand-made. They are generated from the data in `results/` by `scripts/plot_results.py`.

## Solver modes

The executable supports these modes:

- `free_fall`: one-particle analytical verification under gravity
- `constant_velocity`: zero-force constant-velocity verification
- `bounce`: one-particle wall-bounce case
- `verification`: grouped verification cases, including timestep sensitivity
- `experiment`: increasing-particle-count serial experiment
- `scaling`: profiling, OpenMP strong scaling, weak scaling, and serial-versus-parallel checks
- `neighbor_bonus`: brute-force versus neighbour-search comparison
- `science_bonus`: damping and particle-cloud settling studies

## Requirements

To reproduce the runs and rebuild the report, the following tools are needed:

- a C++17 compiler with OpenMP support
- Python 3
- `pip`
- `pdflatex`

Python packages used for plotting:

- `matplotlib`
- `numpy`
- `pandas`

They can be installed with:

```bash
python -m pip install -r requirements.txt
```

## Quick start on Windows

The repository includes a Windows PowerShell script for a one-command run:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1
```

For a fresh rerun that removes old generated outputs first:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -Clean
```

This workflow:

1. checks required tools
2. builds the solver
3. runs all verification, experiment, scaling, neighbour-search, and science-bonus cases
4. installs Python plotting dependencies
5. regenerates the report figures
6. compiles the final PDF

## Manual workflow on Linux or macOS

The packaged automation script is Windows-specific, but the project can be reproduced manually on Linux or macOS with equivalent tools.

### 1. Clone the repository and enter it

```bash
git clone https://github.com/Praveenjhas/HPSC_assignment.git
cd HPSC_assignment
```

### 2. Install Python plotting dependencies

```bash
python3 -m pip install -r requirements.txt
```

### 3. Build the solver

With `clang++`:

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp src/main.cpp -o dem_solver -fopenmp
```

If your system uses `g++` with OpenMP:

```bash
g++ -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp src/main.cpp -o dem_solver
```

### 4. Run the simulation modes

```bash
./dem_solver verification results
./dem_solver free_fall results
./dem_solver constant_velocity results
./dem_solver bounce results
./dem_solver experiment results
./dem_solver scaling results
./dem_solver neighbor_bonus results
./dem_solver science_bonus results
```

### 5. Generate the figures

```bash
python3 scripts/plot_results.py
```

### 6. Compile the report

```bash
cd report
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
```

The final PDF will be:

```text
report/white_paper.pdf
```

## Minimal reproduction options

If a full rerun is not needed, the stages can be run separately.

Build only:

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp src/main.cpp -o dem_solver -fopenmp
```

Run only one mode:

```bash
./dem_solver free_fall results
```

Regenerate only plots:

```bash
python3 scripts/plot_results.py
```

Recompile only the report:

```bash
cd report
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
pdflatex -interaction=nonstopmode -halt-on-error white_paper.tex
```

## Notes on outputs

- Progress information is printed during longer runs, including timestep, simulated time, kinetic energy, center-of-mass height, active contacts, and maximum speed.
- Reported runtimes in the generated summaries are wall-clock runtimes of the simulation runs.
- Different experiments use different physical simulation times depending on purpose, so the reported runtime should not be confused with simulated physical time.

## Optional Makefile shortcuts

If `make` is available, these shortcuts can also be used:

```bash
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

On Windows, `make full` is just a shortcut for `run_all.ps1`.

## Reproducibility note

All figures in the report are generated from repository data and scripts. The intended order is:

1. compile the solver
2. run the required simulation modes
3. generate the figures from `results/`
4. compile the LaTeX report

If these steps are followed, the report can be regenerated from the repository contents.
