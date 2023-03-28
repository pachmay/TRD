
#ifndef AliTrdKrCluster_hh
#define AliTrdKrCluster_hh

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <vector>
// ROOT includes
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TProfile.h>
#include "TProfile2D.h"
#include "TRandom.h"
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTimeStamp.h>

#include <dirent.h>
#include <errno.h>
#include <string>     // std::string, std::stoi

#include "DataFormatsTRD/KrCluster.h"
#include "DataFormatsTRD/KrClusterTriggerRecord.h"
#include "TRDCalibration/KrClusterFinder.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "Math/IntegratorOptions.h"
#include <fairlogger/Logger.h>

using namespace o2::trd;
using namespace o2::trd::constants;
using BCData = o2::InteractionRecord;

class AliTrdKrClusterO2 {
public:
    AliTrdKrClusterO2();
    ~AliTrdKrClusterO2();

    void  setInputDir(TString ipinputdir) {pinputdir = ipinputdir;};
    void  setOutdir(TString ioutputdir)   {outputdir = ioutputdir;};
    void createInputList(string sListOfFiles);
    Int_t createInputList();
    void  createChain();
    void  Init_histograms(Int_t mode);
    void  loopClusters(Long64_t nEventsUse, Int_t mode);
    void  writeOutput(Int_t mode);

    void  setRunNumber(int irunnumber = 0) { runNumber = irunnumber; };
    // void  setPressureData(TGraph* graph) { grP = (TGraph*) graph->Clone("grP"); }; 
    void  setInitialTimePerRun(int iyear, int imonth, int iday, int ihour, int imin, int isec); //{ setyear = iyear; setmonth = imonth; setday = iday; sethour = ihour; setmin = imin; setsec = isec; timeStamp.Set( setyear, setmonth, setday, sethour, setmin, setsec, 0, kTRUE, 0); }; 
    // void  setInitialTimePerRun (int irunnumber = 0);
    void setTreeName(TString treename = "krClusters", TString clustername = "KrCluster", TString triggername = "TriggerRecord" );
    void setFileToRead(TString sinputfile) { sReadFile = sinputfile; };
    void setSector(int isector) { useSector = isector; };

    TTimeStamp getClusterStamp(std::vector<o2::trd::KrClusterTriggerRecord>* vTriggerRecordPtr);
    // TDatime convertTimestampToDatime(TTimeStamp timestamp);

private:

    const Int_t nTrdChambers = 540;
    const Int_t N_TRD = 540;
    const Int_t N_TRD_sectors = 18;
    const Int_t N_TRD_layers  = 6;
    const Int_t N_TRD_stacks  = 5;
    const Int_t N_TRD_rows    = 16;
    const Int_t N_TRD_columns = 144;
    const Int_t N_cuts        = 8;
    Int_t merge_N_rows    = 2;
    Int_t merge_N_columns = 16;
    Int_t N_columns_merge;
    Int_t N_rows_merge;
    Int_t N_run_ids = 1;

    TRandom ran;
    TString HistName;
    char NoP[50];
    TString pinputdir, outputdir;
    TString sReadFile = "none";
    vector<string> files_root;
    Long64_t file_entries_total;

    int useSector = -1;
    Int_t runNumber;
    Int_t setyear;
    Int_t setmonth;
    Int_t setday;
    Int_t sethour;
    Int_t setmin;
    Int_t setsec;
    TTimeStamp timeStamp;


    TString sKrClusterTreeName = "krClusters";
    TString sKrClusterBranchName = "KrCluster";
    TString sKrClusterTriggerRecordBranchName = "TriggerRecord";

    std::vector<o2::trd::KrCluster> vKrCluster;
    std::vector<o2::trd::KrCluster>* vKrClusterPtr = &vKrCluster;
    std::vector<o2::trd::KrClusterTriggerRecord> vKrClusterTriggerRecord;
    std::vector<o2::trd::KrClusterTriggerRecord>* vKrTriggerPtr = &vKrClusterTriggerRecord;

    TChain* inputChain;

    vector<TH1F*> vec_h1D_KrClustersDet;
    vector<TH2F*> vec_h2D_KrClVsTimeDet;

    vector<TH1F*> vec_h_ADC_TRD_chambers;
    vector<TH1F*> vec_h_ADC_TRD_chambers_cut;
    vector<UInt_t> vec_run_ids;
    vector< vector<TH1F*> > vec_h_ADC_TRD_chambers_runid;
    // vector< vector< vector<TH1F*> > > vec_ADC_pads;
    vector<TH2F*> vec_h2D_cls_time_vs_ADC;
    vector<TH1F*> vec_h_ADC_cut;
    TH1F* h_det;
    TH1F* hSector;
    TH2F* h2D_ADC_xy_TRD;
    TH2F* h2D_ADC_xy_TRD_merged;
    TH2F* h2D_xy_TRD_merged;
    TH2F* h2D_rmsTime_vs_ADC;
    TH2F* h2D_cls_size_vs_ADC;
    TH2F* h2D_cls_dim_vs_ADC;
    TH2F* h2D_cls_time_vs_ADC;
    TH2F* h2D_cls_col_vs_ADC;
    TH2F* h2D_cls_row_vs_ADC;
    TH2F* h2D_run_ids_vs_det;
    TH2F* h2D_ADC_pressure;
    TProfile2D* hProfile2D_ADC_xy_TRD_merged;
    vector< vector<TH1F*> > vec_h_lower_upper_cuts;
    vector<TH1F*> vec_h_peak_params;
    // vector< vector< vector< vector<TH1F*> > > > vec_ADC_pads_merge;
    // vector< vector< vector< vector<TH1F*> > > > vec_ADC_pads_merge_corr;
    vector< vector< vector<TH1F*> > > vec_ADC_pads;
    vector< vector< vector<TH1F*> > > vec_ADC_pads_merge;
    vector< vector< vector<TH1F*> > > vec_ADC_pads_merge_corr;
    vector< vector<TH1F*> > vec_h_single_pad_pressure_cat;
    TH2F* h2D_old_Krypton;
    vector<TH2F*> vec_h2D_Krypton_cat;
    vector< vector<TH1F*> > vec_ADC_spec_cat;

    vector<TH1F*> vec_h_deltaTime;
    vector<TH1F*> vec_h_deltaColumn;
    vector<TH1F*> vec_h_deltaRow;
    vector<TH1F*> vec_h_adcRMS;
    vector<TH1F*> vec_h_timeRms;
    vector<TH1F*> vec_h_adcSum;
    vector<TH1F*> vec_h_adcSumOverT;
    TH1F* h_detector_digits;

    TGraph* grP;    //!

    ClassDef(AliTrdKrClusterO2,2)

};


#endif
