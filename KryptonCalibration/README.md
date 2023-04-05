# TRD Krypton calibration documentation

## Code overview:
1. AliTrdKrClusterO2 - Alex Schmah.
2. Ana_ADC_spectra - by Alex Schmah.
3. PadCalibCCDBBuilder - by Jana Crkovska.


## 1. AliTrdKrClusterO2
Code to extract pad-by-pad ADC spectra from clusters.

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
    -  Macro to get the graph(s) from data: 
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


**For next steps you will need `KrHistOutput_<runNumber>_mode2_sector-1.root`.** To get it, run
```
root -l 'Macro_AliTrdKrClusterO2.cc( <runNumber>, 2, -1)' 
```