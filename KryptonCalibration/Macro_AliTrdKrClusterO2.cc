        //    Kr runs timing:
        //    Kr runs 15.09 - 17.09
        //
        //    501318, 2 hrs from 17:00 to 19:30
        //    501354, ~ 7 hours, 9/16/2021, 2:25:04 AM, 103719585 triggers, 136 TB
        //    501376, 78.84 minutes, 9/16/2021, 11:01:27 AM, 18927746 triggers .
        //    501384, ~ 2 hours, 9/16/2021, 1:14:36 PM, 30114577 triggers, 39.8 TB.
        //    501483, 65 TB, 5.26 GB/s, 49578800 triggers. 10h02 to 13h46 17/09
        //
R__LOAD_LIBRARY(AliTrdKrClusterO2_cxx.so);

void Macro_AliTrdKrClusterO2(int runNumber = 0, int mode = 0, int sector = 0, TString sReadFile = "none")
{
    // First compile with
    // alienv enter O2/latest-dev-o2
    // .L AliTrdKrClusterO2.cxx++

    // mode:
    // 0 : merged histograms
    // 2 : pad-by-pad histograms

    cout << runNumber << endl;

    int setyear=0, setmonth=0, setday=0, sethour=0, setmin=0, setsec=0;

    switch ( runNumber ) {
        default:
            setyear = 1980; setmonth = 5; setday = 21; sethour = 1; setmin = 2; setsec = 3;
            // break;
        case 501318: 
            setyear = 2021; setmonth = 9; setday = 15; sethour = 15; setmin = 0; setsec = 0;
            break;
        case 501354:
            setyear = 2021; setmonth = 9; setday = 16; sethour = 0; setmin = 25; setsec = 4;
            break;
        case 501376:
            setyear = 2021; setmonth = 9; setday = 16; sethour = 9; setmin = 1; setsec = 27;
            break;
        case 501384:
            setyear = 2021; setmonth = 9; setday = 16; sethour = 11; setmin = 14; setsec = 36;
            break;
        case 501483:
            setyear = 2021; setmonth = 9; setday = 17; sethour = 8; setmin = 2; setsec = 0;
            break;
    }

    cout << "Macro_AliTrdKrClusterO2 started" << endl;

    gSystem->Load("AliTrdKrClusterO2_cxx.so");

    AliTrdKrClusterO2 *AliTrdKrClusterO2_Ana = new AliTrdKrClusterO2();

    //// set run # and initial time (of first timeframe)
    AliTrdKrClusterO2_Ana->setRunNumber( runNumber );
    AliTrdKrClusterO2_Ana->setSector( sector );
    AliTrdKrClusterO2_Ana->setInitialTimePerRun( setyear, setmonth, setday, sethour, setmin, setsec);
    if( runNumber == 501318 )
    {
        AliTrdKrClusterO2_Ana->setTreeName("o2sim","KRCLUSTER","TRGKRCLS");
    }

    // /misc/alidata121/alice_u/crkovska/ALICE/TRD_Kr_Jan2022/Krypton/501483

    //// set the graph to draw pressure info from
    TFile* file = TFile::Open("/home/ceres/crkovska/ALICE/TRD_Krypton/2021Dec_Krypton/trdkr_pressure_2021Sep5-20.root");
    TGraph* gr = (TGraph*) file->Get("grPressureIn");
    // AliTrdKrClusterO2_Ana->setPressureData(gr);

    //AliTrdKrClusterO2_Ana->setInputDir( Form("/misc/alidata121/alice_u/crkovska/ALICE/TRD_Kr_Jan2022/Krypton/%i/", runNumber) );
    //if ( runNumber == 501318 ) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/Data/Clusters/2021/");

    //------------------
    // HV scan

    //1st HV scan:
    // 0 V: 500993, 9/12/2021, 11:50:55 AM, 15 mins, 3732651 triggers, 5TB
    //+10V: 500996, 9/12/2021, 12:10:55 PM, 15 mins, 3622902 triggers, 4.85TB
    //+20V: 500998, 9/12/2021, 12:29:03 PM (real start 12:40), 15 mins, 3651026 triggers, 4.9 TB
    //-20V: 500999, 9/12/2021, 12:57:25 PM, 15 mins, 3672330 triggers, 4.9 TB
    //-10V: 501002, 9/12/2021, 1:16:58 PM, 15 mins, 3758192 triggers, 5 TB


    //2nd HV scan 15.09.21:
    // 0 V: 501266, 9/15/2021, 2:49:19 PM, 15 mins, 3797479 triggers, 5.04 TB;
    //+10V: 501276, 9/15/2021, 3:10:59 PM, 15 mins, 3768536 triggers, 5.0 TB;
    //+20V: 501283, 9/15/2021, 3:32:29 PM, 15 mins, 3712832 triggers, ~ 5.0 TB;
    //-20V: 501294, 9/15/2021, 3:55:58 PM, 15 mins, 3681439 triggers, 4.87 TB;
    //-10V: 501307, 9/15/2021, 4:22:55 PM, 15 mins, 3698432 triggers, 4.9 TB.

    // Largest Krypton run
    // 501354
    if( runNumber == 0 ) AliTrdKrClusterO2_Ana->setInputDir("/home/ceres/crkovska/ALICE/TRD_Krypton/test_code_alex/");

    if(runNumber == 515974) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata121/alice_u/crkovska/ALICE/2022MayKrypton/clusters/515974/");

    if(runNumber == 500998) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/500998/"); // +20 V
    if(runNumber == 500999) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/500999/"); // -20 V
    if(runNumber == 501002) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501002/"); // -10 V
    if(runNumber == 501276) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501276/");
    if(runNumber == 501283) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501283/");
    if(runNumber == 501294) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501294/");
    if(runNumber == 501307) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501307/");
    if(runNumber == 501266) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/HVscan/501266/");
    if(runNumber == 501354) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/KrData/501354/"); // large singel Krypton run
    if(runNumber == 501)    AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/2021Data/KrData/KrRnd2/"); // all Krypton runs with similar pressure
    //------------------

    if(runNumber == 501483) AliTrdKrClusterO2_Ana->setInputDir("/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/Data/Clusters/2021/501483_new/");

    AliTrdKrClusterO2_Ana->setOutdir("/home/ceres/crkovska/ALICE/TRD_Krypton/test_code_alex/test_output/");
    AliTrdKrClusterO2_Ana->createInputList(); //("/home/ceres/crkovska/ALICE/TRD_Krypton/QA_ClusterFinder/file.txt");
    AliTrdKrClusterO2_Ana->createChain();
    AliTrdKrClusterO2_Ana->Init_histograms(mode);
    AliTrdKrClusterO2_Ana->loopClusters(-1,mode); // -1 = use all events
    AliTrdKrClusterO2_Ana->writeOutput(mode);

}
