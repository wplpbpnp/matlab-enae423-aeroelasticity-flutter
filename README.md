# MATLAB ENAE423 Aeroelasticity and Flutter

MATLAB coursework scripts for ENAE423, including panel flutter analysis and supporting assignments.

## Course Summary (Testudo)

ENAE423 (`Aeroelasticity`) covers static and dynamic aeroelastic phenomena and methods of analysis for divergence, control reversal, flutter, and related instability/response problems.

Testudo: `https://app.testudo.umd.edu/soc/search?courseId=ENAE423&sectionId=&termId=202605`

## Highlight

- `panel_flutter.m` (primary portfolio candidate in this repo)

## Contents

- Homework scripts (`hw*.m`)
- Midterm/practice scripts (`midterm.m`, `practice_midterm.m`)
- Example/reference script (`example.m`)
- Published PDFs in `html/`

## Getting Started (MATLAB)

1. Open MATLAB in the repo root.
2. Start with `panel_flutter.m` (primary project-style script in this repo).
3. The script builds nondimensional mass/stiffness/aerodynamic matrices and sweeps `q_bar` to estimate flutter onset from modal behavior.
4. Use the `hw*.m` scripts for supporting course examples and derivation checks.
5. Compare results with the published PDFs in `html/`, especially `html/panel_flutter.pdf`.

## Dependencies / Compatibility Notes

- `panel_flutter.m` uses Symbolic Math (`syms`, symbolic differentiation/integration) before converting matrices to numeric form.
- The script generates figures for modal damping/frequency trends and prints matrices to the MATLAB console.

## Suggested Showcase Files

- `panel_flutter.m`
- `hw7_redo.m` or `hw8*.m` (if they contain stronger analysis plots)
- `html/panel_flutter.pdf`
