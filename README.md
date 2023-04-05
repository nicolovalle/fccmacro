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

#### Batch mode

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

To run on batch, set:
```python
outputDir = "/eos/user/n/nvalle/FCC/workdir/...."
outputDirEos = "/eos/user/n/nvalle/FCC/workdir/..."
eosType = "eosuser"
runBatch = True
batchQueue = "longlunch"
```

#### Local mode

Set the output directory, comment the batch mode or set it to `False`. To run on local EDM4HEP files, the command is
```
fccanalysis run analysis_testv3.py --files-list edm4hepfile.root
```
The output name will be `output.root`

## My external analysis

The n-tuples are processed by root macros independent from the fccanalysis. The main analysis is implemented in `MyExternalAnalysis/MyExternalAnalysis.C/h`. It produced an `AnalysisResults.root` with some debug information and the `eventsTree` tree.
To run the macro one shall use `runMyExternalAnalysis.C`. It requires a manual link to the fastjetlibraries. To automatize it, a script can be used:
```bash
./runrunMyExternalAnalysis.sh filelist.txt [output_suffix]
```
which created an additional root macros containing the loading of `libfastjet` end the execution of `runMyExternalAnalysis.C`. Such macro is removed at the end. *All the file lists are in `MyExternalAnalysis/`* and must be kept updated with the proper outputs of the `fccanalysis`. For the Zbb and Zcc samples of spring2021 (~1B events), the files are split into 8 files and the list of the list names is written in `LIST_of_Z*.txt`. The parallelization can be done with the following script
```
parallel -j 8 ./script_for_bkg_parallel.sh < LIST_of_Z...
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
void SkimmingTask(TString fin, TString outdir="./results/skimmed/", TString suffix="", TString opt="loose")
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
+ `loose`
+ `forcutvar`

The events are discarded if the conditions are not fullfilled for *any of the clustering alogorithms*. This allows to run on the skimmed file using any of the jet methods. See the table below for the cuts implemented with each option. 





## Analysis macros

The root macros are in the `workdir/fccmacro` directory, local copy form github
+ https://github.com/nicolovalle/fccmacro

There are different branch and the main one. The branch name should correspond the the `MyExternalAnalysis` macro version and the `results-V.../` files that can be processed with that version of the macros.


#### lumisettings.h

In this file the Luminosity, cross sections, weights, coupling, file names are written or computed by dedicated methods receiving as input the sample type, mass, lifetime,...

In addition, the function `MapPoint()` can be used to display the grid of available signal data points.

#### Cut.h

Here the main TChain `TREE` is defined and linked to `eventsTree`. The branch variables are declared and the addresses are set. To use it in a macro looping over the events tree, one has to:
+ **IMPORTANT** Reset the chain at the beginning: `TREE->Reset();`. If not done, files are constantly added when the same function (e.g. CutFlowOK) is called recursively by other functions.
+ Add root files to the chain: `TREE->Add(<filename>);`
+ Set the branches addresses before the loop: `TREESetBranch();`
+ Call `BUILD_DERIVATE(int jalg)` at the beginning of each interation of the loop. This initialize some derivate variables, also delcared in `Cut.h`.
+ Call `SELECTION_*(...)` in the loop, to apply cuts (a collection of boolean functions returning `true` if the event passes selections).

Finally, the function `std::vector<float> GetFloatArray(TStrning analysis_opt)` is used to intepret the substring of "analysis_opt" in between `[`,`]` (see below).
`Cut.h` is the code stiring the cut selection used in the `CutFlowOK.C` and `getvalues.C` macros. The plotting macros and other analyses use `CutFlowOK.C` and `getvalues.C` to get the cut flows or the observable spectra.

#### CutFlowOK.C

```cpp
std::map<int, std::vector<double>> CutFlowOK(TString opt="signal", Int_t mass=80, TString lifetime = "n/a", TString dir="../MyExternalAnalysis/results/", Long64_t RunOnN = -1, Bool_t CutByCutFlow = false, Int_t jalg = 2, TString analysis_opt="< d2d dsigma anymass1L2M")
```

*The implementation with `CutByCutFlow = true` is not supported*
To be read like this: `CutFlowOK[M]` where `M`(int) is the mass in GeV. CutFlowOK[M] is a vector:
+ `0`: Nnocut(opt,mass,lifetime). This is written in `lumisettings.h`, not read from the root file.
+ `1`: (tree entries && nOneMuon==1)  `2`: selection  `3`: sliding[M]. 
+ `CutFlowOK[M][dcut_id(dcut)]`, where `dcut` is an integer:
  + In units of sigmas for analyses with option `dsigma` (see below for analyses options)
  + In units of mm/100 for analysese with option `dmm`

```cpp
void makeHTMLtable(TString outfile="./summary.html", int iopt = 1, TString AnalysisResPath = "../MyExternalAnalysis/results/skimmed/", Int_t jalg = 2,  Bool_t ProcessOnlySignal = false, Int_t RunOnN = -1, Bool_t upload = false, TString comments_on_top="")
```

*The implementation with `iopt = 2` was intended to print cut by cut efficiency but it is not supported*



#### getvalues.C

```cpp
std::pair<std::vector<Double_t>, Double_t> getvalues(TString obsID, TString opt="signal", Int_t mass = 50, TString lifetime = "n/a", Long64_t RunOnN = -1, Double_t d0cut=8, Int_t jalg = 2, TString analysis_opt = "< d2d dsigma anymass1L2M", TString dir="../MyExternalAnalysis/results/")
```
See below for analyses options. See the code for possible observables IDs.
+ `getvalues.first` is non binned array of values for the observable
+ `getvalues.second` is the Scale Factor, smaller than 1 in case `RunOnN` is positive and smaller than the number of entries of the processed tree.


### Analyses options

This is string inteprepret by many of the macros reading the events tree. According to this string, different selections can be chosed. The event selection must be always done using the methods implemented in `Cut.h`

```cpp
TString analysis_opt = "< d2d dsigma anymass1L2M"
```
Analysis options available (separated by space):
+ `>` or `<` for the type of cut in impact parameter
+ `dsigma` or `dmm` to cut in number of sigma or in unit 10-5 m
+ `d2d` or `d3d` to cut on D0 or D0+Z0
+ Type of analysis:
  + `anymass1L2M` : mass independent, it uses old analysis for all cases with 2 Jets and LM-1j analysis for all cases with 1 jet
  + `anymass1L2L` : mass independent, it uses LM-1j analysis or LM-2j analyses according to number of jets
  + `cutvariation` : this shall be used together with options between brackets `[...]` which are then decoded by `Cut.h::GetFloatArray`
    + For the moment, the option between brackets must be formatted like this: `[ABCDEfghiJKLM...]`, i.e. as a list of 5-characters numbers.


### Cut list

That's the list of cuts currently implemented.

|                 | 1J LM  | 2J LM  | 2J MM   | loose      | forcutvar |
| --------------- | ------ | ------ | ------- | ---------- | --------- |
| cos(pmiss)      | < 0.94 | 1JLM   | < 0.94  | < 0.94     | < 0.94    |
| cos(pmiss,mu)   | < 0.50 | 1JLM   | < 0.80  | < 0.80     | < 0.80    |
| min [Ej]        | >= 3   | 1JLM   | >= 3    | >= 3       | >= 3      |
| cos(jj)         |        | > -0.8 | > -0.8  | > -0.8     | > -0.96   |
| cos(jj)         |        |        | < 0.98  | any value  | any value |
| max [cos(j,mu)] | < 0.96 | < 0.96 | < 0.8   | < 0.96     | < 0.98    |
| max [cos(j,mu)] | > -0.5 | > -0.5 |         | any value  | any value |
| min [cos(j,mu)] |        |        | > -0.98 | any value  | any value |
| min [Mj]        | > 0.2  | > 0.2  | > 0.2   | > 0.2      | > 0.2     |
| min [Mj2]       | > 0    | > 0    | > 0     | any value  | any value |
| M [vis+miss]    | > 80   | > 80   | > 80    | > 80       | > 80      |






:::
To created a new branch
+ Synchronize local/upstream
+ `git checkout -b <new_branch>`
+ `git push origin <new_branch>`
https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
:::