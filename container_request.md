# Request: new AFNI Apptainer image with current R (separate from the existing one)

**Please build this as a new image, and leave the existing AFNI image (`/apps/afni/afni_make_build_latest.sif`) and the cluster's R module as they are.**
Each Apptainer image is a single read-only file with its own copy of R inside. A new image can't change the old one or the cluster's R.

## Why the existing image can't just be updated (option 2)
Some AFNI programs (the group statistics and reliability tools) are written in R. They need R add-on packages that the current image doesn't have.

- **The image's R is too old for those packages.** The current image is built on Ubuntu 18.04 and includes R 3.4.4 (from 2018).
  - The main package we need, `lme4`, depends on `Matrix`. Current `Matrix` requires R 4.4 or newer, so current packages can't be installed on R 3.4.4.
  - We also tried installing older 2019-era versions into the old image, to match its R. Several of them failed to compile, so that workaround doesn't work either.
- **R can't be upgraded inside that image.** CRAN's Ubuntu repository, the usual source for newer R, no longer supports Ubuntu 18.04.

So the fix is a new image on a newer Ubuntu, with a current R.

## What the new image needs
AFNI publishes official setup scripts for Ubuntu 24.04. Following them avoids guesswork:
[admin script](https://github.com/afni/afni/blob/master/src/other_builds/OS_notes.linux_ubuntu_24_64_a_admin.txt) and
[user script](https://github.com/afni/afni/blob/master/src/other_builds/OS_notes.linux_ubuntu_24_64_b_user.tcsh).

1. **Ubuntu 24.04 as the base.**
2. **The system packages from AFNI's admin script.** That list includes the brms extras (`libgdal-dev libopenblas-dev libnode-dev libudunits2-dev`).
3. **R 4.4 or newer, from CRAN's Ubuntu repository** ([instructions](https://cran.r-project.org/bin/linux/ubuntu/)). Ubuntu 24.04's built-in R is 4.3.3, which is still too old.
4. **AFNI's Ubuntu 24 binaries** (package `linux_ubuntu_24_64`), installed inside the image.
5. **These R packages, installed inside the image** (in R's system-wide library, not anyone's home folder):
   - AFNI's standard set: afex, phia, snow, nlme, lmerTest, gamm4, data.table, paran, psych, corrplot, metafor
   - Also: lme4, blme, brms, rstan, coda, ggplot2, irr
   - Optional: MuMIn, Hmisc

   Everything has to be in the image. Compute nodes have no internet, and AFNI tries to download any missing package at run time, which fails there. If a package doesn't install, the build should stop with an error; by default R only prints a warning and carries on.
6. **Keep the C++ compilers in the finished image.** brms compiles a small program every time it runs.
7. **Make the image ignore users' personal R packages.** Set `R_LIBS_USER=NULL` in the image's environment section. Otherwise R could load packages from someone's home folder instead of the image's. It must be the word `NULL`; an empty value doesn't work.
8. **Give it a new filename that includes the AFNI and R versions,** so it's clear which image is which.

## How to tell it works
After building, please confirm that:
- **R is version 4.4 or newer.**
- **All the packages above load without errors.**
- **AFNI's R-based programs start up** (3dICC, 3dLMEr, 3dMEMA and TRR). Each connects AFNI to R through a small compiled file. If that file was built for a different R version, the program crashes as soon as it starts, even just printing its help text. Starting them is therefore a quick test that AFNI and R match.
- **brms can fit a tiny test model,** which confirms the compilers work.
- **AFNI's own system check is clean** (`afni_system_check.py -check_all`).

## Appendix: test commands
These run the checks above, in the same order. Put them in the image's `%test` section, or run each one as `apptainer exec --cleanenv <new image> <command>`:
```bash
R --version | head -1                  # R 4.4 or newer
afni -ver                              # AFNI version, for the record
rPkgsInstall -pkgs ALL -check          # nothing should say "not installed"
Rscript -e 'for (p in c("lme4","blme","metafor","snow","data.table","lmerTest","phia","brms","rstan","coda","ggplot2","irr")) library(p, character.only=TRUE)'
for p in 3dICC 3dLMEr 3dMEMA TRR; do $p -help > /dev/null && echo "$p ok"; done
Rscript -e 'brms::brm(y ~ 1, data = data.frame(y = rnorm(20)), chains = 1, iter = 200, refresh = 0)'
afni_system_check.py -check_all
```
