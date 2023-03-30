# HN analysis workflow


## From .lhe to EDM4HEP

The .lhe files are processed by Pythia+Delphes through a Pythia card where the lhe file is specified:
```
DelphesPythia8_EDM4HEP FCC-config/FCCee/Delphes/card_IDEA.tcl edm4hep_output_config.tcl pythia_card.cmd outputFile.root
```

+ The `FCC-config` repository is https://github.com/HEP-FCC/FCC-config. The `spring2021` branch has been used up to *now*.
+ The EMD4HEP congif output contains the collections to be saved in the output file: 
    ```
    module EDM4HepOutput EDM4HepOutput {
    add ReconstructedParticleCollections EFlowTrack EFlowPhoton EFlowNeutralHadron
    add GenParticleCollections           Particle
    add JetCollections                   Jet
    add MuonCollections                  AllMuon
    add ElectronCollections              Electron
    add PhotonCollections                Photon
    add MissingETCollections             MissingET
    add ScalarHTCollections              ScalarHT
    set RecoParticleCollectionName       ReconstructedParticles
    set MCRecoAssociationCollectionName  MCRecoAssociations
    }
    ```
+ Finally, the Pythia card is as follows:
    ```
    Main:numberOfEvents = 1000000000
    Beams:frameType             = 4
    Beams:setProductionScalesFromLHEF = off
    Beams:LHEF = XXX.lhe
    Main:timesAllowErrors = 100000
    ```

## From EDMH4HEP to n-tuples

The first analysis is done with the `fccanalysis` software. Tha analysis is a python script, to be launched with:
```
fccanalysis run analysis_tesv3.py
```
In the script, the dataset is a local EDM4HEP file or it's chosen from official collections. The spring2021 ones are [at the following link](http://fcc-physics-events.web.cern.ch/fcc-physics-events/FCCee/spring2021/Delphesevents_IDEA.php). The py code snippet linking to productions on eos is:
```python
processList = {
        'p8_ee_Zbb_ecm91':{'chunks':100},
        'p8_ee_Zcc_ecm91':{'chunks':100,'fraction':1},
        ...
}

```

## My external analysis

The n-tuples are processed by root macros independent from the fccanalysis. The main analysis is implemented in `MyExternalAnalysis/MyExternalAnalysis.C/h`. It produced an `AnalysisResults.root` with some debug information and the `eventsTree` tree.
To run the macro one shall use `runMyExternalAnalysis.C`. It requires a manual link to the fastjetlibraries. To automatize it, a script can be used:
```bash
./runrunMyExternalAnalysis.sh filelist.txt [output_suffix]
```
which created an additional root macros containing the loading of `libfastjet` end the execution of `runMyExternalAnalysis.C`. Such macro is removed at the end. *All the file lists are in `MyExternalAnalysis/`* and must be kept updated with the proper outputs of the `fccanalysis`. For the Zbb and Zcc samples of spring2021 (~1B events), the files are split into 8 files and the list of the list names is written in `LIST_of_Z*.txt`. The parallelization can be done with the following script
```
paralle -j 8 ./script_for_bkg_parallel.sh < LIST_of_Z...
```
It produces in the local directory 8 output files and 8 LOG files (no output on the screen) having as suffix the .txt file list names. They must be merged through `hadd` (see naming conventions below).

Each change in the code should be preceeded by a manual change of the macro version in `MyExternalAnalysis.C`. E.g.:
```cpp
Int_t MacroVersion = 230322;
```
The output files produced with a given version should go in a separate folder `MyExternalAnalysis/results-V<version>`. The one in use **MUST** be linked as `results/`. Furthermore the `results/` directory will contain one or more `skimmed/` directories with the output of the skimmed task. The macro version should finally correspond to the branch name of the "fccmacro" repo that can be used to analyize those files. 

### Multiple jet algorithms

The Jet 4-momenta variables in the `eventsTree` output tree (`oPxJet1, ... , oEJet2`), as well as,  `oNJet` are std vectors. in the macro version 230322, the components are as follows:
+ **jalg==0** -> Durham kt algorithm with number of jet fixed to 2. (`durham_kt_njets_mode`)
+ **jalg==1** -> Durham kt algorithm with 2 jets possibly merged if rtd < 5.0. (`durham_kt_njets_mode1`)
+ **jalg==2** -> First, it clusterizes in any number of jets according to Durham kt with rtd cut = 10.0 (`durham_kt_dcut_mode`). If the number of resulting jets is > 2, the method jalg==0 is used.
 
jalg==0 (`durham_kt_njets_mode`) is the method used since the version before 230222.

### Output file naming conventions

The naming convention is defined into the `fccmacro/lumisettings.h` file by the `AnalysisResults` function. It is as follows:
```cpp
if      (opt == "Zbb")     toret = "AnalysisResults-Zbb_highstat.root";
else if (opt == "Zcc")     toret = "AnalysisResults-Zcc_highstat.root";
else if (opt == "Zuds")    toret = "AnalysisResults-Zuds_highstat.root";
else if (opt == "Zmumu")   toret = "AnalysisResults-Zmumu.root";
else if (opt == "Ztautau") toret = "AnalysisResults-Ztautau.root";
else if (opt == "munuqq")  toret = "AnalysisResults-munuqq.root";
else if (opt == "signal" && lifetime == "n/a")  toret = Form("AnalysisResults-signal-M-%s.root",HNMass.Data());
else if (lifetime == "n/a") toret = Form("AnalysisResults-signal-M-%s.root",opt.Data());
else if (lifetime != "n/a") toret = Form("AnalysisResults-signal_10k_%s_%s.root",HNMass.Data(),lifetime.Data());
```

## Filtered trees

A skimming task is available to filter the `AnalysisResults.root` with a first event selection. It reduces the Zbb size by a factor >1000. The root macro `SkimmingTask.C` is in `MyExternalAnalysis/`. 
```cpp
void SkimmingTask(TString fin, TString outdir="./results/skimmed/", TString suffix="", TString opt="extraloose")
```
A) if fin is a .txt file, output file is `<outdir>/AnalysisResults<suffix>.root`
B) if fin is `AnalysisResult<any>.root`, outputfile is `<outdir>/AnalysisResults<any><suffix>.root`
C) if fin is a .root file not beginning with "AnalysisResults", option A. is used.
    
A script is also available to process with default option all the `AnalysisResults*.root` files present in `results/`:
```bash
./skimm_all.sh
```

The output of the macro is a `AnalysisResults.root` file with the tree `eventsTree` having the same structure of the original one, but less events.


The options currently implemented are:
+ `extraloose`
+ `loose`

The events are discarded if the conditions are not fullfilled for *at least one of the clustering alogorithms*. This allows to run on the skimmed file using any of the jet methods.

**extraloose** implements the following cuts:
+ |cos(pmiss)| < 0.94
+ cos(pmiss,mu) < 0.8
+ min(E_j) >= 3
+ cos(j1,j2) > -0.8
+ max [cos(ji,mu)] < 0.96
+ min [M(ji)] > 0.2
+ M(pmiss+pvis) > 80 

**loose** adds the following cut:
+ min [cos(ji,mu)] > -0.98q

## Analysis macros

The root macros are in the `workdir/fccmacro` directory, local copy form github
+ https://github.com/nicolovalle/fccmacro

There are different branch and the main one. The branch name should correspond the the `MyExternalAnalysis` macro version and the `results-V.../` files that can be processed with that version of the macros.

:::spoiler To created a new branch
+ Synchronize local/upstream
+ `git checkout -b <new_branch>`
+ `git push origin <new_branch>`
https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
:::