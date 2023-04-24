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
root -l -b -q 'Ana_ADC_spectra.cc++(-1,0,-1,"KrHistOutput_<runNumber>_mode0_sector-1.root","Fits_mode0.root")'
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

## 3. PadCalibCCDBBuilder


Class `PadCalibCCDBBuilder` is used to populate detector maps with results of the fits to the ADC pad-level spectra from the previous step, inter/extrapolating missing information, normalizing the gain to the detector average, and preparing a `LocalGainFactor` object for CCDB upload.

### Class contains

- PadCalibCCDBBuilder.cxx
- PadCalibCCDBBuilder.h
- CreateCCDBLocalGainFactor.C

The code is commited to O2 under the `Detectors/TRD/.` directory.

### Compilation


In principle, the class can be run fully locally without relying on O$^2$ (by downloading the class and macro locally and simply commenting `namespace o2` and `namespace trd` and checking dependencies point to local files). Doing this may be advantageous if one needs to modify and test the code as doing so allows one to work on the code without having to rebuild O$^2$. The only part which *requires* O$^2$ is the one populating the LocalGainFactor and uploading it to the CCDB. 

To run locally without O2, implement the following changes if you wish to test the code locally without having to rely on O2:
1. Comment out all notions of `namespace o2` and `namespace trd`.
2. Upate the list of headers in the `.cxx` file, the only relevant header not native to C++ or ROOT is `#include "PadCalibCCDBBuilder.h"`. This header also needs to be included in the steering macro. 
3. In the `.cxx` file, add the following global variables:
```
int NCOLUMN = 144;
int NROWC0 = 12;
int NROWC1 = 16;
```

Then, compile the class using your local ROOT installation:
```
root -l
.L PadCalibCCDBBuilder.cxx++
```

### Run the code
Once the class is compiled, run the steering macro `CreateCCDBLocalGainFactor.C`. As input file, feed it the output from the previous step with pad-by-pad fits (so option `(-1,2,-1)`). The class run over the `nt_Krypton` TTree, which stores all the fit results.

**However**, if you use the code curently (i.e., April 2, 2023) committed to O2, you need to implement two changes:
1.  Loop over all chambers: `idet < 1` --> `idet < constants::MAXCHAMBER` (or `idet < 540` if you run without O2).
2.  Call `PadCalibCCDBBuilder::populateEmptyNormalizedMap` once you normalize all filled maps (**before upload!**) so that **all the maps are filled with values that are not 0**. With the macro as it is now, some maps would remain empty and this will then mess up the workflow - the correction is applied in the denominator so you absolutely don’t want this empty! 
3.  The committed macro has a return statement that currently cuts the code before CCDB upload. This was used for some testing so all you need to do if you want to upload something to the CCDB is to comment this away.
4.  (Optional) You may comment the line `cout << calObject.getValue(0, 52, 13) << endl;` (see [l. 52](https://github.com/AliceO2Group/AliceO2/blob/33ecc0535ceae089c2bae81e2d751c21c85bb18d/Detectors/TRD/macros/CreateCCDBLocalGainFactor.C#L52) of the committed code). It was there just to check that the maps are indeed filled and contain some reasonable values.

Remember that in order to be able to upload stuff to the CCDB, you need a valid token.

#### Upload path and start time
The maps are set to be uploaded into `"TRD/Calib/LocalGainFactor"`.

The calibration is set to be valid from `start_time` to `end_time`. Leave the latter as is and for the former set some reasonable time, e.g., at 00:00 on the first day of the Krypton run.

### What is stored in the uploaded object
What is actually uploaded to the CCDB is a `LocalGainFactor` object, storing a local gain value for each pad as:
```
(detector ID , column ID , row ID , relative gain)
```
The relative gain has to be a number larger than 0. 

### The analysis procedure
The code takes the fitted gain values, pad-by-pad, and populates 2D representations of individual detectors. These 2D maps are defined as arrays of x=column $\times$ y=row pads. For most, this means 144 $\times$ 16 pads. For detectors in sector 2, there is  144 $\times$ 12 pads.

#### Selecting "good fits"
The root file produced in step 2 may in principle contain results of failed fits. Only pads whose fit pass a quality selection will be filled. To define appropriate cuts, the user will need to check the results of the previous step and set cuts that will work the best for a given data-taking period. In the code committed to the $O^2$, the following selection criteria have been hard-coded:
- $\chi^2>0$
- amplitude $>0$
- stand. dev $>0$ && stand. dev $<1000$

See method `PadCalibCCDBBuilder::getDetectorMap`, [link](https://github.com/AliceO2Group/AliceO2/blob/972204bd5ae56e62b182144de938e64c7f9e3203/Detectors/TRD/calibration/src/PadCalibCCDBBuilder.cxx#L386).
Only values passing these criteria will be filled into the 2D maps. These criteria should serve as a reasonable first estimate for code testing etc.

**In the here available code, a setter method allowing the user to easily set these criteria has been added.**

#### "Smoothening" the maps
The populated maps likely contain "hot" pads or regions with large difference between neigboring pad conent. From previous data-taking campaigns, we know that within one chamber, the pad-by-pad gains should create a smooth distribution, maybe with some low-gain pads at the edges.

The method `PadCalibCCDBBuilder::smoothenTheDetector` serves to get rid of isolated clusters of hot pads or pairs of pads with a too large a difference between their respective gains. In other words, all pads deemed as sus by this method will be set to 0. The "too large a distance" can be set by the user, the default value is 1000.

The result will be a detector map with empty pads in place where an inhomogoeneity was found.

**Right now, the code does not *always* handle hot pads among other filled pads correctly. The issue is being looked into (Apr 21, 2023) and hopefully will be fixed by end of May 2023. However, if you see this comment in June 2023 onwardsm the issue is still there.**

#### Filling the empty pads in populated maps
The method `PadCalibCCDBBuilder::fillTheMap` takes all maps which are **not empty** and fills the empty pads in an iterative procedure. This method does not fill the original histogram but rather creates a new one.

The filling scheme is as follows:
1. The routine checks for empty pad and stores their coordinates.
1. Each of these pads is filled with the content of a pad mirrored with respect to the y-axis. For example, an empty pad (column=0, row=0) will be filled with the content of pad (143,0).
2. If the bertically mirrored pad is empty, the pad is filled with the content of the pad mirrored with rerspect to the x-axis. So the example pad from above (0,0) will be filled with the content of the pad (0,15) (or in case of some chamebrs with (0,11)).
3. If both vertically and hotizontally mirrored pads are empty, the pad is filled with the average of contents of its closest neighbors. At least one neighboring pad must be filled. Otherwise the pad is skipped in the given round and filled later on.
4. The entire process is repeated until there are no more empty pads. So for this to work, the original map needs to contain at least one filled pad.

All the values in the filled pads are set to be negative. This is to keep the info on which values were inter/extrapolated ("*estimated*") and which come from the fitting procedure itself.

Please note that maps for chambers with no response, this method will do nothing. These need to be dealt with later.

#### Normalizing the maps
This is done using the method `PadCalibCCDBBuilder::createNormalizedMap`. This method takes a filled map, computes its mean (also including the estimated pad values), and creates a new map which will have pads filled with a ratio of the gain in the original map with respect to the map average:
```
normalized gain = gain in pad / map average gain.
```

#### Filling the empty maps
At last, the empty maps are to be filled. 
The method `PadCalibCCDBBuilder::populateEmptyNormalizedMap` will set all the pads to -1.
**This method is included it in the code available in O^2, but is not used in the steering macro that comes with the code.**

<span style="color:red">**Do not forget to fill the empty maps too. The correction is applied in the denominator, so leaving this empty will mess up everything.**</span>

The code also contains the method `PadCalibCCDBBuilder::transformMapIntoAbsoluteValues`, which takes the filled, normalized maps and replaces all the estimated values with their respective absolute values. This makes it easier to visualy inspect the maps and check if they are indeed smooth and filled properly.

#### Plotting user friendly maps
The final maps to be uploaded to the database contain a bunch of pads with negative gain. This makes reading the maps hard. For user convenience, the class contains a method `PadCalibCCDBBuilder::transformMapIntoAbsoluteValues`, which takes the filled map and creates a new one with all pads set to their absolute gain.

### Uploading the maps to the CCDB
In 2022, the maps were to be uploaded to the Test-CCDB by the user, the subsequent upload from test to actual CCDB had to be requested. Check this is still true and then call:

```
o2::ccdb::CcdbApi ccdb;
ccdb.init("http://ccdb-test.cern.ch:8080");
```
Next you can define metadata for your upload - this can be any kind of note you think may be useful or name of the user who did the upload, e.g.,:
```
std::map<std::string, std::string> metadata;
metadata.insert({"Responsible","User Name (user.email@cern.ch)"});
metadata.insert({"Note","calibration to be tested with PID"});
```
The map `metadata` MUST be initiated but you can leave it empty (as in the code example in this repo).

The validity of the calibration is handled by its start and end timestamps. These can be computed from time_t as
```
time_t start_time = 1631311200; // Fri Sep 10 2021 22:00:00 GMT+0000
time_t end_time = 2208985200;   // Sat Dec 31 2039 23:00:00 GMT+0000
auto timeStamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::from_time_t(start_time).time_since_epoch()).count();
auto timeStampEnd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::from_time_t(end_time).time_since_epoch()).count();
```

The `LocalGainFactor calObject` is then uploaded with
```
ccdb.storeAsTFileAny(&calObject, "TRD/Calib/LocalGainFactor", metadata, timeStamp, timeStampEnd);
```

You can check that your object is uploaded to the path you entered in your browser at [http://ccdb-test.cern.ch:8080/browse/TRD/Calib/LocalGainFactor](http://ccdb-test.cern.ch:8080/browse/TRD/Calib/LocalGainFactor) (update for the path you used).
