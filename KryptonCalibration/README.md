# TRD Krypton calibration documentation

## Code overview:
1. AliTrdKrClusterO2 - Alex Schmah.
2. Ana_ADC_spectra - by Alex Schmah.
3. PadCalibCCDBBuilder - by Jana Crkovska.

The outputs of step 1 and 2 from the analysis of the 2022 Kr data can be found in `/misc/alidata121/alice_u/crkovska/ALICE/2022MayKrypton/example_rootfiles`. The parent folder also contains cluster data and code for the 2022 analysis.

## 1. AliTrdKrClusterO2
Code to extract pad-by-pad ADC spectra from clusters.

This code runs on `trdkrclusters_<runNumber>_<part>.root` files. These are not produced centrally. Instead they need to be produced from digits specifically for the purpose of calibration. 

### Class contains
- AliTrdKrClusterO2.cxx
- AliTrdKrClusterO2.h
- Macro_AliTrdKrClusterO2.cc

### Compilation
Compile code using your local O2, such as:
```
alienv enter O2/latest-dev-o2
.L AliTrdKrClusterO2.cxx++
```

### Run the code
To run the code, use the steering macro `Macro_AliTrdKrClusterO2.cc`. 

#### Setting the steering macro
1. Set the date-time info for your data:
    - The date-time was set up to be used to retrieve information on pressure inside ALICE but ultimately was not deemed necessary for the two Kr runs in 2021 and 2022 (stable pressure thorought each Kr run).
    - However you must set up some date, you may use what is now set as default (all zeros).
2. (Optional) Set a graph to draw gas pressure info from:
    -  Data can be found at [https://darma.cern.ch/](https://darma.cern.ch/ ).
    -  Macro to get the graph(s) from data is included in the parent repo, in `Pressure/makePressureGraphs.C`. 
    -  Use method `AliTrdKrClusterO2::setPressureData` to set the graph (you need to first uncomment and recompile). The time-date information will then be used to get the correct pressure during a given run. **The pressure information has not been used in 2021/2022 calibration (stable pressure) but the code has been left in for future convenience.**
3. Set input and output directories:
    - Set input folder using use `AliTrdKrClusterO2::setInputDir`. The path should point to the folder containing your `trdkrclusters_*.root` files. **The path MUST end with an `/`**, e.g., `AliTrdKrClusterO2->setInputDir("/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/Data/Clusters/2021/501483_new/")"`. **This applies to ALL paths you may need to set in this macro.**
    - Set output folder with `AliTrdKrClusterO2::setOutdir`. 


#### Running the steering macro
To run the macro, you need to give it as arguments:
- `runNumber`, which is used to set the proper input path and time information. With the current setup you would typically store the data in folders corresponding to dufferent runs so the runNumber parameter can be your switch. 
- `mode`, which is used to switch whether you get "merged-pad" ADC histograms (pads merged over 2 rows and 15 columns, use option 0) or pad-by-pad ADC histograms (use option 2).
- `sector`, which controls whether you run over a specific sector (numbers 0 to 17) or over the full TRD (use option -1).

The last argument is by default `sReadFile = "none"`, leave it as it is.

So, if you want to get merged histograms from sector 1 for run 123456, you run
```
root -l 'Macro_AliTrdKrClusterO2.cc( 123456, 0, 1)' 
```
Then in your chosen output directory, you will find a file `KrHistOutput_<runNumber>_mode<mode>_sector<sector>.root`, so using the above example you should look for file `KrHistOutput_123456_mode0_sector1.root`.


**In the next step you will need:**
- **`KrHistOutput_<runNumber>_mode0_sector-1.root`,**
- **`KrHistOutput_<runNumber>_mode2_sector-1.root`.** 
 
To get them, run
```
root -l 'Macro_AliTrdKrClusterO2.cc( <runNumber>, 0, -1)'
```
and then
```
root -l 'Macro_AliTrdKrClusterO2.cc( <runNumber>, 2, -1)' 
```

## 2. Ana_ADC_spectra

Code to fit the ADC spectra obtained in the previous step. The fitting procedure is split into two steps:
1. Fit spectra in merged pads (there is 72 per detecctor). The mean and sigma from each gaussian fit is filled into a TProfile, which contains the average sigma and mean per detector. These are then used as starting parameters for the poad-by-pad fits.
2. Fit spectra in each pad using the fit results from previous step to define realistic starting fit parameters and to impose reasonable limits for these parameters. 

N.B. If correcting for pressure is deemed necessary, an additional step needs to be included between the two mentioned above, during which the pressure corrected spectra are fitted.

### Class contains
- Ana_ADC_spectra.cc
- Ana_ADC_spectra.h
- bad_chambers.h


### Run the code
The fitting procedure is executed by the function
```
void Ana_ADC_spectra(<sector>, <mode>, <det>, <input_file_with_spectra>,<output_file_for_fits>,<input_file_with_merged_fits>)
```
The function uses the following arguments:
- `sector`, which controls whether you run over a specific sector (numbers 0 to 17) or over the full TRD (use option -1).
- `mode`, which is used to switch whether you get "merged-pad" ADC histograms (pads merged over 2 rows and 15 columns, use option 0) or pad-by-pad ADC histograms (use option 2).
- `det`, which allows the user to pick a single chamber to run QA checks.
- `input_file_with_spectra` is the file with merged or pad0by-pad ADC spectra, which the user obtained in the Part 1.
- `output_file_for_fits`, which controls the name of the new file to store the fit results.
- `input_file_with_merged_fits`, which allows the user to specify the name of the file with fits over merged spectra. These are used to constrain the pad-by-pad fits so that they converge.

First you need to fit the merged spectra, which you obtained from running the code described in **Part 1** in **mode 0**:
```
root -l -b -q 'Ana_ADC_spectra.cc++(-1,0,-1,"KrHistOutput_<runNumber>_mode2_sector-1.root","Fits_mode0.root")'
```

Using the fit results stored in `"Fits_mode0.root"` and the pad ADC spectra obtained in Part 1 with **mode 2**:
```
root -l -b -q 'Ana_ADC_spectra.cc++(-1,2,-1,"KrHistOutput_<runNumber>_mode2_sector-1.root","Fits_mode2.root","Fits_mode0.root")'
```

**In the next and final step, you will need**:
- **The output in mode 2, `Fits_mode2.root`.**


### Bad chambers
The header `bad_chambers.h` allows the user to define chambers flagged as bad, e.g., the ones with rediced HV or off.

This feature is not currently used in the code. The header is anyway included in the repo as the current code requires it to run.