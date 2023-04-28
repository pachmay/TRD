
// Alexander Schmah, August 2018

#include "Ana_ADC_spectra.h"

void Ana_ADC_spectra(Int_t sector_use = -1, Int_t mode_use = 0, Int_t det_use = -1, TString sinputfile = "KrHistOutput_515974_mode0_sector-1.root",
                    TString soutputfile = "Fits.root", TString sfitinputfile = "Fits_mode0.root")
{

    // 1. root -b -q Ana_ADC_spectra.cc++\(-1,0,-1\)
    // 2. root -b -q Ana_ADC_spectra.cc++\(-1,2,-1\)


    // sector_use: sector you want to analyze
    // det_use:    for QA, pick a single detector from the sector you are analyzing to see the ADC specta and fits
    // set det_use to -1 for not showing those histograms at all

    // For pressure correction: Check all <ADC> vs pressure plots, the fits need to be OK.
    // Check vec_c_slopes_1: slope/b needs to be relatively narrow sigma/mu ~< 0.06
    // Some chambers are not working or have some parts which are off, that is usually OK.
    // If something doesn't look right then check the ADC fits of this chamber by setting det_use to the chamber number id

    // mode_use:
    // 0: use merged histograms (72 per detector) to determine pressure correction
    // 1: load pressure corrected pad wise ADC spectra and perform fits //// --> NOT AVAILABLE IN THIS VERSION BUT THE METHODS EXIST IN .H!!!
    // 2: closure test, load fully corrected (pressure + Krypton) pad wise ADC spectra


    //------------------------------------------------------------------------
    cout << "Ana_ADC_spectra opened" << endl;
    gStyle->SetPalette(56); // 53 = black body radiation, 56 = inverted black body radiator, 103 = sunset, 87 == light temperature
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    // main_data_dir = "/misc/alidata141/alice_u/jcrkovsk/TRD/KrOutputs/";
    // main_data_dir = "/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/output/2021/Output/";
    main_data_dir = "/home/ceres/crkovska/ALICE/TRD_Krypton/test_code_alex/";
    Init_input(sector_use,mode_use,sinputfile);
    //------------------------------------------------------------------------


    Int_t sector_for_det_use = det_use/30;
    if(sector_for_det_use != sector_use) det_use = -1;

    N_columns_merge = N_TRD_columns/merge_N_columns;
    N_rows_merge    = N_TRD_rows/merge_N_rows;

    Init_functions();
    vec_merge_fit_par.resize(N_TRD);
    vec_merge_mean_pressure_fit_par.resize(N_TRD);
    
    //------------------------------------------------------------------------
    if(mode_use == 0)
    {
        init_ADC_pad_merge(sector_use);
        fit_ADC_pad_merge(sector_use);
        save_pressure_correction(sector_use, mode_use, soutputfile);

    }
    else if ( mode_use == 2)
    {
        printf("init_ADC_pad_wise \n");
        init_ADC_pad_wise(sector_use,mode_use,soutputfile);
        // printf("init_peak_positions \n");
        init_peak_positions(sfitinputfile); //("Fits_501_mode_1.root");
        printf("fit_ADC_pad_wise \n");
        fit_ADC_pad_wise(sector_use);
        printf("save_Krypton_calibration \n");
        save_Krypton_calibration(sector_use,mode_use);
        printf("done \n");
    }
    

    cout << "bye bye" << endl;


}
