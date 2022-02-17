## Simulation study for Spatial Statistics paper

This repository has two branches, one for "sequential" simulation (the `main`
branch) and another one for "parallel" simulation (the `parallel`
branch). Although both of them are reproducible, the latter uses a
[SingularityCE](https://sylabs.io/singularity/) container that assures the
reproducibility and allows for running the code in _HPC_ clusters.

Both branches share a similar structure, having the following directories:

* `src` - contains the source code for:
  - Our `R` package `tpsa`;
  - A supplementary package called `tpsautils` that contains helper functions
    for data simulation
  - A zip file (`biom13115-sup-0001-gcops_codec.zip`) containing the software
    Lavancier et al. 2020 submitted to Biometrics accompanying their paper.
* `scripts` - the `R` scripts to run the run the simulations;
* `data` - directory to receive the results from the scripts.
* `lavancier` - directory where Lavancier software was unzipped and
  installed.^[Note that, if using the "sequential" you have to install this
  software by yourself at the exact same locations where we did so.]

Also, the branches contain the `spatstat.def` file, which is used to build the
container used in the "parallel" branch.

Also, we have to comment the lines 161 and 2141 in the
`programs/colocalizationTest.cpp` file from Lavancier's software. Although this
does not change the results of their functions, it was needed to run their
software from `R`. Again, the results are EXACTLY the same regardless whether
those lines were commented out or not.

---

If using the "parallel" branch, you need to install
[SingularityCE](https://sylabs.io/singularity/) on your own computer. Moreover,
if you do not use Linux, you will need to set up a virtual linux machine to run
it. The command
```
sudo singularity build --force spatstat.sif spatastat.def
```
generates the singularity image to run the simulations. The image can be used to
run the simulations on a HPC cluster.

The file `run-sim.sh` (available only in the "parallel" branch) uses
`spatstat.sif` to run the simulations using
[SLURM](https://slurm.schedmd.com/overview.html).

---
