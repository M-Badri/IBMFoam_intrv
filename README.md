# IBMFoam

## Overview
IBMFoam is a **free, open-source CFD simulation library** built on top of [OpenFOAM](https://openfoam.org).  
It enables simulation of **flows laden with arbitrarily shaped solids** using a combination of advanced numerical techniques for fluid–particle interaction.

The solver implements the **Hybrid Fictitious Domain – Immersed Boundary Method**, coupled with the **Discrete Element Method (DEM)** for accurate particle dynamics.  
Its initial implementation is derived from the work of Federico Municchi, with additional capabilities for large-scale parallel simulation, post-processing, and flexible geometry handling.

---

## Key Features
- **Centralized simulation control** via `IBMFoamdict` (see the detailed documentation).
- **Particle geometry options**: spherical particles or arbitrary STL-defined shapes (see `IBMFoam/geom_models`).
- **Two-phase simulations** (solid–fluid) with full coupling between phases (`pimpleIBMFoam` solver).
- **Standard DEM mode** for solid-phase-only simulations.
- **Adaptive mesh refinement (AMR)** based on particle positions for higher resolution in regions of interest.
- **Spring–dashpot contact model**:
  - Elastic modulus–based stiffness
  - Damping calculated from coefficient of restitution (see `IBMFoam/contactModels`)
- **Adhesive force model** extension to contact physics.
- **Flexible solid phase initialization**:
  - Random spatial distribution
  - Uniform sizing
  - Import of pre-defined body arrangements (see `IBMFoam/add_models`)
- **Instructive tutorials**:  
  - Single particle falling through fluid
  - Interaction between a particle and complex impeller geometries (see `Tutorials`)
- **Application examples** for large-scale and high-performance runs (see `examples`).
- **Fully parallelized** for distributed-memory architectures.
- **Optional debug and profiling mode** for performance analysis (new).

---

## Compatibility
Tested with:
- OpenFOAM v8 ([release details](https://openfoam.org/version/8/))
- GCC 7.x or later
- Linux x86_64 systems (Ubuntu and CentOS tested)  
Other configurations may work but are not officially validated.

---

## Compilation Instructions
**Important:** Run all scripts from a terminal where OpenFOAM v8 has been sourced.

**To build all components:**
```bash
./compileAll.sh
```
This will compile:
- Core `IBMFoam` library
- `pimpleIBMFoam` solver
- Geometry, contact, and additional model libraries

For selective compilation, see `compile.sh` scripts in individual subdirectories.

---

## Usage Notes
- `IBMFoamdict` contains all main input parameters; see **Documentation for IBMFoamdict** for details.
- STL geometry files for complex bodies should be placed in `constant/triSurface/`.
- Tutorials provide working starting points for most applications.
- Additional **benchmark cases** are available in the `examples` directory for performance testing.

---

## Development & Contributions
IBMFoam is developed by a small research team and is provided **as-is** without warranty, in the spirit of open science.  
We welcome:
- Bug reports
- Feature requests
- Pull requests with tested improvements

If you use IBMFoam in academic work, please consider citing our repository.

---

## License
IBMFoam is licensed under the **GNU General Public License v3** (or later).  
You are free to use, modify, and redistribute the software under the GPL terms.  
See [GNU GPL v3](http://www.gnu.org/licenses/) for details.
