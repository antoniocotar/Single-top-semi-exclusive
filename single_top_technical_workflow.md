# Single-Top Semi-Exclusive Technical Workflow

## Scope

This document summarizes the technical workflow used for the single-top study in proton-proton collisions at a center-of-mass energy of 13 TeV, considering semi-exclusive production via Pomeron-induced processes and photoproduction. The purpose of the workflow, up to the current stage, is to produce processed samples and inspect the reconstructed observables that can later be used to separate signal from the dominant backgrounds, especially central `pp -> t tbar`.

## 1. Signal and Background Samples

### 1.1 Pomeron-induced single-top generation

The Pomeron single-top samples were generated in:

- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/top_pomeron_pomset6_pomflux1_default`
- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/top_pomeron_pomset3_pomflux6`
- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/top_pomeron_pomset6_pomflux7`
- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/top_pomeron_pomset4_pomflux7`
- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/top_pomeron_pomset6_pomflux8`

After generation, these samples were sent to the Delphes-level storage area:

- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/delphes`

### 1.2 Photoproduction signal sample

The photoproduction single-top sample was already available as:

- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/twphoto.root`

### 1.3 Background samples

The background samples were taken from:

- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/delphes_ttbar_ap_upc_NLO_pp13TeV.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/delphes_ttbar_pa_upc_NLO_pp13TeV.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/delphes_ttbar_pom_default_pp_13TeV.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/delphes_ttbar_central_NLO_pp13TeV.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/dy.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/w1.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/w2.root`
- `/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/w3.root`

These files include the dominant central and semi-exclusive backgrounds considered in the study.

## 2. Event Processing

The generated and pre-existing Delphes samples were processed with the general analysis macro located in:

- `/lustre/home/amendieta/single-top_semi-exclusive/github/complete_chain/processed/general_sample_analysis/compute_general_sample.C`
- `/lustre/home/amendieta/single-top_semi-exclusive/github/complete_chain/processed/general_sample_analysis/run_compute_general_sample.sh`

The working directory on the cluster was:

- `/lustre/home/amendieta/single-top_semi-exclusive/github/complete_chain/processed/general_sample_analysis`

This stage was used to compute the event-level quantities needed for the subsequent analysis and variable inspection.

The macro writes:

- a `PomeronTree` with `x_pomeron`,
- an `Events` tree with reconstructed observables,
- and several histograms related to track multiplicity, rapidity-gap observables, and partonic momentum fractions.

The main reconstructed variables stored in the `Events` tree are:

- `rapidityGapPlus`
- `rapidityGapMinus`
- `hfSumPlus`
- `hfSumMinus`
- `nchPV`
- `nForwardJets`
- `nTotalJets`
- `leadLightJetPt`
- `maxHFEnergy`
- `minHFEnergy`
- `MET`
- `maxHFEnergyPlus`
- `maxHFEnergyMinus`
- `minOfMaxHF`
- `minOfHFSums`

In addition, the macro fills distributions for:

- charged-particle multiplicity (`hNTrk`),
- forward and backward rapidity gaps (`hFRG`, `hBRG`),
- Bjorken-like momentum fractions for photon, gluon, and Pomeron components (`hPhotonX`, `hGluonX`, `hPomeronX`),
- and light-quark/gluon momentum fractions from the Pomeron and proton (`hPomQuark`, `hPomGluon`, `hProtQuark`, `hProtGluon`, `hU_Pom`, `hD_Pom`, `hS_Pom`, `hU_Pro`, `hD_Pro`, `hS_Pro`).

## 3. Parallelized Processing Stage

The macro was parallelized in:

- `/lustre/home/amendieta/guille/Pp_Alessio/singletop/samples_processed`

After this step, the processed outputs were collected in:

- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado`

The integrated processed files are:

- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/delphes_ttbar_ap_pa_NLO_pp13Tev_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/delphes_ttbar_central_NLO_pp13TeV_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/delphes_ttbar_pa_upc_NLO_pp13TeV_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/delphes_ttbar_pom_default_pp_13TeV_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/dy_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/top_pomeron_pomset3_pomflux6_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/top_pomeron_pomset4_pomflux7_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/top_pomeron_pomset6_pomflux1_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/top_pomeron_pomset6_pomflux7_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/top_pomeron_pomset6_pomflux8_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/twphoto_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/w1_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/w2_processed.root`
- `/lustre/proyectos/p018/single-top_semi-exclusive/samples/integrado/w3_processed.root`

## 4. Plotting and Exploratory Analysis on lxplus

After the processed ROOT files were integrated, they were analyzed with notebooks on `lxplus`. The plotting and exploratory analysis stage was carried out in:

- `/eos/user/g/guillera/Single-top-semi-exclusive/notebooks`

At this stage, distributions were produced for the signal and background samples in order to inspect the separation power of the reconstructed observables.
