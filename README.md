## Simulation study for Spatial Statistics paper

This repository has two options for running the simulations, a "sequential" one
(using files ending in `_sequential` branch) and a "parallel" one (using files
ending in `_parallel`). Although both of them are reproducible, the latter uses
a [SingularityCE](https://sylabs.io/singularity/) container that assures the
reproducibility and allows for running the code on _HPC_ clusters.

The structure of the repository is as follows

* `src` - contains the source code for:
  - Our `R` package `tpsa`;
  - A supplementary package called `tpsautils` that contains helper functions
    for data simulation
  - A zip file (`biom13115-sup-0001-gcops_codec.zip`) containing the software
    Lavancier et al. 2020 submitted to Biometrics accompanying their paper.
* `scripts` - the `R` scripts to run the run the simulations;
* `data` - directory to receive the results from the scripts.
* `lavancier` - directory where Lavancier software was unzipped and
  installed[^1].

Also, the `spatstat.def` file can be used to build the container used in the
"parallel" simualtions.

Moreover, we have to comment the lines 161 and 2141 in the
`programs/colocalizationTest.cpp` file from Lavancier's software. Although this
does not change the results of their functions, it was needed to run their
software from `R`. Again, the results are EXACTLY the same regardless whether
those lines were commented out or not.

---

If using the "parallel" option, you need to install
[SingularityCE](https://sylabs.io/singularity/) on your own
computer. Furthermore, if you do not use Linux, you will need to set up a
virtual linux machine to run it. The command
```
sudo singularity build --force spatstat.sif spatastat.def
```
generates the singularity image to run the simulations. The image can be used to
run the simulations on a HPC cluster.

The file `run-sim.sh` uses `spatstat.sif` to run the simulations using
[SLURM](https://slurm.schedmd.com/overview.html).

---

[^1]: Note that, if using the "sequential" option you have to Lavancier's
    software from the sources in `biom13115-sup-0001-gcops_codec.zip` (This is
    the implementation they made available along with their paper).
