#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <vector>
// ROOT includes
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <dirent.h>
#include <errno.h>
#include <string>     // std::string, std::stoi

#include "DataFormatsTRD/KrCluster.h"
#include "DataFormatsTRD/KrClusterTriggerRecord.h"
#include "TRDCalibration/KrClusterFinder.h"
#include "Math/IntegratorOptions.h"
#include <fairlogger/Logger.h>

using namespace o2::trd;
using namespace o2::trd::constants;

///////// /misc/alidata141/alice_u/jcrkovsk/TRD/Calibration/Krypton/Data/raw/2018/o2test

void KrClustersQA( TString sInputFile = "trdkrclusters_test.root", TString sOutputFileName = "output.root" )
{

    //------------------------------------------------------------------------
    // Add data to input chain
    TString sInputDirName = "";
    // TString pinputdir = "/misc/alidata141/alice_u/schmah/TRD/Calibration/Krypton/Data/raw/2018/hists/";
    uint64_t file_entries_total = 0;

    TString sKrClusterTreeName   = "krClusters";
    TString sKrClusterBranchName = "KrCluster";
    TString sKrClusterTriggerRecordBranchName = "TriggerRecord";

    std::vector<o2::trd::KrCluster> vKrCluster;
    std::vector<o2::trd::KrCluster>* vKrClusterPtr = &vKrCluster;
    std::vector<o2::trd::KrClusterTriggerRecord> vKrClusterTriggerRecord;
    std::vector<o2::trd::KrClusterTriggerRecord>* vKrTriggerPtr = &vKrClusterTriggerRecord;


    TChain* inputChain;
    inputChain  = new TChain( sKrClusterTreeName.Data(), sKrClusterTreeName.Data() );
    uint64_t entries_save;
    vector<string> files_root;

    ifstream inFileList( sInputFile );
    string filename = "";
    while( inFileList >> filename ) {
        files_root.push_back( filename );
    }

    // files_root.push_back( sInputFile.Data() );

    for(int i_file = 0; i_file < (int)files_root.size(); i_file++)
    {
        TString addfile = sInputDirName;
        addfile += files_root[i_file];
        inputChain->AddFile(addfile.Data(),-1, sKrClusterTreeName.Data() );
        uint64_t file_entries = inputChain->GetEntries();
        cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
        entries_save = file_entries;
    }
    inputChain->SetBranchAddress( sKrClusterBranchName, &vKrClusterPtr );
    inputChain->SetBranchAddress( sKrClusterTriggerRecordBranchName, &vKrTriggerPtr );


    file_entries_total = inputChain->GetEntries();
    cout << "Total number of events in tree: " << file_entries_total << endl;
    //------------------------------------------------------------------------


    TH1F* hAdcMaxA = new TH1F("hAdcMaxA", "hAdcMaxA", 400, 0, 20000 );
    hAdcMaxA->GetXaxis()->SetTitle("ADC^{max}_{A}");

    TH1F* hAdcMaxB = new TH1F("hAdcMaxB", "hAdcMaxB", 400, 0, 20000 );
    hAdcMaxB->GetXaxis()->SetTitle("ADC^{max}_{B}");

    TH1F* hAdcSum = new TH1F( "hAdcSum", "hAdcSum", 400, 0, 20000 );
    hAdcSum->GetXaxis()->SetTitle("sum ADC");

    TH1F* hAdcSumOverThr = new TH1F( "hAdcSumOverThr", "hAdcSumOverThr", 400, 0, 20000 );
    hAdcSumOverThr->GetXaxis()->SetTitle("sum ADC over threshold");

    TH1F* hAdcRms = new TH1F("hAdcRms", "hAdcRms", 500, 0, 1000);
    hAdcRms->GetXaxis()->SetTitle("RMS(ADC)");

    TH1F* hTimeRms = new TH1F("hTimeRms", "hTimeRms", 50, 0, 50);
    hTimeRms->GetXaxis()->SetTitle("RMS(t)");

    TH1F* hDeltaTime = new TH1F("hDeltaTime", "hDeltaTime", 50, 0, 50);
    hDeltaTime->GetXaxis()->SetTitle("#Deltat");

    TH1F* hDeltaRow = new TH1F("hDeltaRow", "hDeltaRow", 5, 0, 5);
    hDeltaRow->GetXaxis()->SetTitle("#Deltarow");

    TH1F* hDeltaColumn = new TH1F("hDeltaColumn", "hDeltaColumn", 5, 0, 5);
    hDeltaColumn->GetXaxis()->SetTitle("#Deltacolumn");

    TH1F* hSector = new TH1F( "hSector", "hSector", 20, 0, 20);
    hSector->GetXaxis()->SetTitle("sector");

    TH2F* hAdcMaxAVsAdcRms = new TH2F( "hAdcMaxAVsAdcRms", "hAdcMaxAVsAdcRms", 400, 0, 20000, 50, 0, 50);
    hAdcMaxAVsAdcRms->GetXaxis()->SetTitle("ADC^{max}_{A}");
    hAdcMaxAVsAdcRms->GetYaxis()->SetTitle("RMS(ADC)");

    TH2F* hTrdDigitsDistribution = new TH2F( "hTrdDigitsDistribution", "hTrdDigitsDistribution", 2600, 0, 2600, 500, 0, 500 );
    hTrdDigitsDistribution->GetXaxis()->SetTitle("x_{TRD} = col + sec*144");
    hTrdDigitsDistribution->GetYaxis()->SetTitle("y_{TRD} = row + stack*16 + layer*16*5" );

    TH2F* hAdcMaxAVsAdcMaxB = new TH2F( "hAdcMaxAVsAdcMaxB", "hAdcMaxAVsAdcMaxB", 2000, 0, 20000, 2000, 0, 20000 );
    hAdcMaxAVsAdcMaxB->GetXaxis()->SetTitle("ADC^{max}_{A}");
    hAdcMaxAVsAdcMaxB->GetYaxis()->SetTitle("ADC^{max}_{B}"); 

    TH2F* hDeltaTimeAdcMaxA = new TH2F( "hDeltaTimeAdcMaxA", "hDeltaTimeAdcMaxA", 50, 0, 50, 2000, 0, 20000 );
    hDeltaTimeAdcMaxA->GetXaxis()->SetTitle("#Deltat");
    hDeltaTimeAdcMaxA->GetYaxis()->SetTitle("ADC^{max}_{A}");

    TH2F* hDeltaTimeAdcMaxB = new TH2F( "hDeltaTimeAdcMaxB", "hDeltaTimeAdcMaxB", 50, 0, 50, 2000, 0, 20000 );
    hDeltaTimeAdcMaxB->GetXaxis()->SetTitle("#Deltat");
    hDeltaTimeAdcMaxB->GetYaxis()->SetTitle("ADC^{max}_{B}");

    TH2F* hAdcRmsTimeRms = new TH2F( "hAdcRmsTimeRms", "hAdcRmsTimeRms", 500, 0, 1000, 50, 0, 50 );
    hAdcRmsTimeRms->GetXaxis()->SetTitle("RMS(ADC)");
    hAdcRmsTimeRms->GetYaxis()->SetTitle("RMS(t)");

    // timecut

    TH1F* hAdcMaxATimecut = new TH1F("hAdcMaxATimecut", "hAdcMaxATimecut", 400, 0, 20000 );
    hAdcMaxATimecut->GetXaxis()->SetTitle("ADC^{max}_{A}");

    TH1F* hAdcMaxBTimecut = new TH1F("hAdcMaxBTimecut", "hAdcMaxBTimecut", 400, 0, 20000 );
    hAdcMaxBTimecut->GetXaxis()->SetTitle("ADC^{max}_{B}");

    TH1F* hAdcSumTimecut = new TH1F( "hAdcSumTimecut", "hAdcSumTimecut", 400, 0, 20000 );
    hAdcSumTimecut->GetXaxis()->SetTitle("sum ADC");

    TH1F* hAdcSumOverThrTimecut = new TH1F( "hAdcSumOverThrTimecut", "hAdcSumOverThrTimecut", 400, 0, 20000 );
    hAdcSumOverThrTimecut->GetXaxis()->SetTitle("sum ADC over threshold");

    TH2F* hAdcMaxAVsAdcMaxBTimecut = new TH2F( "hAdcMaxAVsAdcMaxBTimecut", "hAdcMaxAVsAdcMaxBTimecut", 2000, 0, 20000, 2000, 0, 20000 );
    hAdcMaxAVsAdcMaxBTimecut->GetXaxis()->SetTitle("ADC^{max}_{A}");
    hAdcMaxAVsAdcMaxBTimecut->GetYaxis()->SetTitle("ADC^{max}_{B}"); 

    // time + row + column cuts

    TH1F* hAdcMaxATimeColRowcut = new TH1F("hAdcMaxATimeColRowcut", "hAdcMaxATimeColRowcut", 400, 0, 20000 );
    hAdcMaxATimeColRowcut->GetXaxis()->SetTitle("ADC^{max}_{A}");

    TH1F* hAdcMaxBTimeColRowcut = new TH1F("hAdcMaxBTimeColRowcut", "hAdcMaxBTimeColRowcut", 400, 0, 20000 );
    hAdcMaxBTimeColRowcut->GetXaxis()->SetTitle("ADC^{max}_{B}");

    TH1F* hAdcSumTimeColRowcut = new TH1F( "hAdcSumTimeColRowcut", "hAdcSumTimeColRowcut", 400, 0, 20000 );
    hAdcSumTimeColRowcut->GetXaxis()->SetTitle("sum ADC");

    TH1F* hAdcSumOverThrTimeColRowcut = new TH1F( "hAdcSumOverThrTimeColRowcut", "hAdcSumOverThrTimeColRowcut", 400, 0, 20000 );
    hAdcSumOverThrTimeColRowcut->GetXaxis()->SetTitle("sum ADC over threshold");

    TH1F* hAdcRmsTimeColRowcut = new TH1F("hAdcRmsTimeColRowcut", "hAdcRmsTimeColRowcut", 500, 0, 1000);
    hAdcRmsTimeColRowcut->GetXaxis()->SetTitle("RMS(ADC)");

    TH1F* hTimeRmsTimeColRowcut = new TH1F("hTimeRmsTimeColRowcut", "hTimeRmsTimeColRowcut", 50, 0, 50);
    hTimeRmsTimeColRowcut->GetXaxis()->SetTitle("RMS(t)");

    TH2F* hAdcMaxAVsAdcMaxBTimeColRowcut = new TH2F( "hAdcMaxAVsAdcMaxBTimeColRowcut", "hAdcMaxAVsAdcMaxBTimeColRowcut", 2000, 0, 20000, 2000, 0, 20000 );
    hAdcMaxAVsAdcMaxBTimeColRowcut->GetXaxis()->SetTitle("ADC^{max}_{A}");
    hAdcMaxAVsAdcMaxBTimeColRowcut->GetYaxis()->SetTitle("ADC^{max}_{B}"); 

    TH2F* hAdcRmsTimeRmsTimeColRowcut = new TH2F( "hAdcRmsTimeRmsTimeColRowcut", "hAdcRmsTimeRmsTimeColRowcut", 500, 0, 1000, 50, 0, 50 );
    hAdcRmsTimeRmsTimeColRowcut->GetXaxis()->SetTitle("RMS(ADC)");
    hAdcRmsTimeRmsTimeColRowcut->GetYaxis()->SetTitle("RMS(t)");

    const int nDet = 540;
    int nPlotDetector = 57;

    TH1F* hAdcSum_det[nDet];
    TH1F* hAdcSumOverThr_det[nDet];
    TH1F* hAdcMaxA_det[nDet];
    TH1F* hAdcMaxB_det[nDet];
    TH2F* hAdcMaxAVsAdcMaxB_det[nDet];
    //
    TH1F* hAdcSumTimeColRowcut_det[nDet];
    TH1F* hAdcSumOverThrTimeColRowcut_det[nDet];
    TH1F* hAdcMaxATimeColRowcut_det[nDet];
    TH1F* hAdcMaxBTimeColRowcut_det[nDet];
    TH2F* hAdcMaxAVsAdcMaxBTimeColRowcut_det[nDet];

    TString sHistoName = "";

    for (int iDet = 0; iDet < nDet; iDet++ ) {

        sHistoName = Form("hAdcSum_det%i", iDet);
        hAdcSum_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcSum_det[iDet]->GetXaxis()->SetTitle("sum ADC");

        sHistoName = Form("hAdcSumOverThr_det%i", iDet);
        hAdcSumOverThr_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcSumOverThr_det[iDet]->GetXaxis()->SetTitle("sum ADC");

        sHistoName = Form("hAdcMaxA_det%i", iDet);
        hAdcMaxA_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcMaxA_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{A}");

        sHistoName = Form("hAdcMaxB_det%i", iDet);
        hAdcMaxB_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcMaxB_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{B}");

        sHistoName = Form("hAdcMaxAVsAdcMaxB_det%i", iDet);
        hAdcMaxAVsAdcMaxB_det[iDet] = new TH2F( sHistoName.Data(), sHistoName.Data(), 2000, 0, 20000, 2000, 0, 20000 );
        hAdcMaxAVsAdcMaxB_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{A}");
        hAdcMaxAVsAdcMaxB_det[iDet]->GetYaxis()->SetTitle("ADC^{max}_{B}");

        /////

        sHistoName = Form("hAdcSumTimeColRowcut_det%i", iDet);
        hAdcSumTimeColRowcut_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcSumTimeColRowcut_det[iDet]->GetXaxis()->SetTitle("sum ADC");

        sHistoName = Form("hAdcSumOverThrTimeColRowcut_det%i", iDet);
        hAdcSumOverThrTimeColRowcut_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcSumOverThrTimeColRowcut_det[iDet]->GetXaxis()->SetTitle("sum ADC");

        sHistoName = Form("hAdcMaxATimeColRowcut_det%i", iDet);
        hAdcMaxATimeColRowcut_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcMaxATimeColRowcut_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{A}");

        sHistoName = Form("hAdcMaxBTimeColRowcut_det%i", iDet);
        hAdcMaxBTimeColRowcut_det[iDet] = new TH1F( sHistoName.Data(), sHistoName.Data(), 400, 0, 20000 );
        hAdcMaxBTimeColRowcut_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{B}");

        sHistoName = Form("hAdcMaxAVsAdcMaxBTimeColRowcut_det%i", iDet);
        hAdcMaxAVsAdcMaxBTimeColRowcut_det[iDet] = new TH2F( sHistoName.Data(), sHistoName.Data(), 2000, 0, 20000, 2000, 0, 20000 );
        hAdcMaxAVsAdcMaxBTimeColRowcut_det[iDet]->GetXaxis()->SetTitle("ADC^{max}_{A}");
        hAdcMaxAVsAdcMaxBTimeColRowcut_det[iDet]->GetYaxis()->SetTitle("ADC^{max}_{B}");

    }

    //******************
    //** cuts
    int timecut = 28;   // maximum value of delta t
    int columncut = 1;  // maximum value of delta column
    int rowcut = 0;     // maximum value of delta row
    //******************

    for (int ievent=0; ievent < file_entries_total; ievent++ ) {
    
    inputChain->GetEntry( ievent );

    int nNumberKlusters = vKrClusterPtr->size();
    cout << "Number of Kr clusters: " << nNumberKlusters << endl;

    // int totalClusters = vKrTriggerPtr->size();
    // o2::trd::KrClusterTriggerRecord record = vKrTriggerPtr->at( ievent ); 
    // cout << "Total nr of clusters: " <<  record.getNumberOfClusters() << " written in vector of size " << totalClusters <<  endl;

    // return;

    for (int i=0; i<nNumberKlusters; i++ ) {

        int detector = (vKrClusterPtr->at(i)).getDetector();
        int sector   = (vKrClusterPtr->at(i)).getSector();
        int column   = (vKrClusterPtr->at(i)).getColumn();
        int row      = (vKrClusterPtr->at(i)).getRow();
        int layer    = (vKrClusterPtr->at(i)).getLayer();
        int stack    = (vKrClusterPtr->at(i)).getStack();

        hSector->Fill( sector );

        int xtrd = column + sector*144;
        int ytrd = row + stack*16 + layer*16*5;

        hTrdDigitsDistribution->Fill( xtrd, ytrd );

        int adcRms = (vKrClusterPtr->at(i)).getAdcRms();
        hAdcRms->Fill( adcRms );

        int timeRms = (vKrClusterPtr->at(i)).getTimeRms();
        hTimeRms->Fill( timeRms );

        int deltaTime = (vKrClusterPtr->at(i)).getClSizeTime();
        hDeltaTime->Fill( deltaTime );

        int deltaRow = (vKrClusterPtr->at(i)).getClSizeRow();
        hDeltaRow->Fill( deltaRow );

        int deltaColumn = (vKrClusterPtr->at(i)).getClSizeCol();
        hDeltaColumn->Fill( deltaColumn );

        int adcMaxA = (vKrClusterPtr->at(i)).getAdcMaxA();
        // if ( adcMaxA < 20 ) continue;
        hAdcMaxA->Fill( adcMaxA );

        int adcMaxB = (vKrClusterPtr->at(i)).getAdcMaxB();
        hAdcMaxB->Fill( adcMaxB );

        int adcSum = (vKrClusterPtr->at(i)).getAdcSum();
        hAdcSum->Fill( adcSum );

        int adcSumOverT = (vKrClusterPtr->at(i)).getAdcSumEoverT();
        hAdcSumOverThr->Fill( adcSumOverT );

        hAdcMaxAVsAdcMaxB->Fill( adcMaxA, adcMaxB );
        hDeltaTimeAdcMaxA->Fill( deltaTime, adcMaxA );
        hDeltaTimeAdcMaxB->Fill( deltaTime, adcMaxB );

        hAdcMaxAVsAdcRms->Fill( adcMaxA, adcRms );

        hAdcRmsTimeRms->Fill( adcRms, timeRms );

        //// fill per detector histos
        hAdcSum_det[detector]->Fill( adcSum );
        hAdcSumOverThr_det[detector]->Fill( adcSumOverT );
        hAdcMaxA_det[detector]->Fill( adcMaxA );
        hAdcMaxB_det[detector]->Fill( adcMaxB );
        hAdcMaxAVsAdcMaxB_det[detector]->Fill( adcMaxA, adcMaxB );

        ////*************************************************
        ////*************************************************
        //// NOW ALL SUBSEQUENT WILL HAVE TIMECUT APPLIED
        ////*************************************************
        ////*************************************************
        if ( deltaTime > timecut ) { 
            continue; 
        }

        hAdcMaxATimecut->Fill( adcMaxA );
        hAdcMaxBTimecut->Fill( adcMaxB );
        hAdcSumTimecut->Fill( adcSum );
        hAdcSumOverThrTimecut->Fill( adcSumOverT );
        hAdcMaxAVsAdcMaxBTimecut->Fill( adcMaxA, adcMaxB );


        ////*************************************************
        ////*************************************************
        //// NOW ALL SUBSEQUENT WILL HAVE ALL CUTS APPLIED
        ////*************************************************
        ////*************************************************
        if ( deltaTime > timecut || deltaColumn > columncut || deltaRow > rowcut ) { 
            continue; 
        }
            
        hAdcSumTimeColRowcut_det[detector]->Fill( adcSum );
        hAdcSumOverThrTimeColRowcut_det[detector]->Fill( adcSumOverT );
        hAdcMaxATimeColRowcut_det[detector]->Fill( adcMaxA );
        hAdcMaxBTimeColRowcut_det[detector]->Fill( adcMaxB );
        hAdcMaxAVsAdcMaxBTimeColRowcut_det[detector]->Fill( adcMaxA, adcMaxB );

        hAdcMaxATimeColRowcut->Fill( adcMaxA );
        hAdcMaxBTimeColRowcut->Fill( adcMaxB );
        hAdcSumTimeColRowcut->Fill( adcSum );
        hAdcSumOverThrTimeColRowcut->Fill( adcSumOverT );
        hAdcMaxAVsAdcMaxBTimeColRowcut->Fill( adcMaxA, adcMaxB );

        hTimeRmsTimeColRowcut->Fill( timeRms );
        hAdcRmsTimeColRowcut->Fill( adcRms );
        hAdcRmsTimeRmsTimeColRowcut->Fill( adcRms, timeRms );

        

    }

    }//ievent loop

    TFile* fout = new TFile( sOutputFileName, "RECREATE");
    hSector->Write();
    hTrdDigitsDistribution->Write();
    hTimeRms->Write();
    hDeltaTime->Write();
    hDeltaRow->Write();
    hDeltaColumn->Write();
    hAdcMaxA->Write();
    hAdcMaxB->Write();
    hAdcSum->Write();
    hAdcSumOverThr->Write();
    hAdcRms->Write();
    hAdcRmsTimeRms->Write();
    hAdcMaxAVsAdcMaxB->Write();
    hDeltaTimeAdcMaxA->Write();
    hDeltaTimeAdcMaxB->Write();
    hAdcMaxAVsAdcRms->Write();
    hAdcMaxATimecut->Write();
    hAdcMaxBTimecut->Write();
    hAdcSumTimecut->Write();
    hAdcSumOverThrTimecut->Write();
    hAdcMaxAVsAdcMaxBTimecut->Write();
    hAdcMaxATimeColRowcut->Write();
    hAdcMaxBTimeColRowcut->Write();
    hAdcSumTimeColRowcut->Write();
    hAdcSumOverThrTimeColRowcut->Write();
    hAdcMaxAVsAdcMaxBTimeColRowcut->Write();
    hTimeRmsTimeColRowcut->Write();
    hAdcRmsTimeColRowcut->Write();
    hAdcRmsTimeRmsTimeColRowcut->Write();
    //// per detector
    fout->mkdir("AdcSumPerDet");
    fout->cd("AdcSumPerDet");
    for( int i=0; i<nDet; i++) {
        hAdcSum_det[i]->Write();
        hAdcSumTimeColRowcut_det[i]->Write();
    }
    fout->cd();
    //
    fout->mkdir("AdcSumOverThrPerDet");
    fout->cd("AdcSumOverThrPerDet");
    for( int i=0; i<nDet; i++ ) {
        hAdcSumOverThr_det[i]->Write();
        hAdcSumOverThrTimeColRowcut_det[i]->Write();
    }
    fout->cd();
    //
    fout->mkdir("AdcMaxAPerDet");
    fout->cd("AdcMaxAPerDet");
    for( int i=0; i<nDet; i++) {
        hAdcMaxA_det[i]->Write();
        hAdcMaxATimeColRowcut_det[i]->Write();
    }
    fout->cd();
    //
    fout->mkdir("AdcMaxBPerDet");
    fout->cd("AdcMaxBPerDet");
    for( int i=0; i<nDet; i++) {
        hAdcMaxB_det[i]->Write();
        hAdcMaxBTimeColRowcut_det[i]->Write();
    }
    fout->cd();
    //
    fout->mkdir("AdcMaxAVsAdcMaxBPerDet");
    fout->cd("AdcMaxAVsAdcMaxBPerDet");
    for( int i=0; i<nDet; i++) {
        hAdcMaxAVsAdcMaxB_det[i]->Write();
        hAdcMaxAVsAdcMaxBTimeColRowcut_det[i]->Write();
    }
    fout->cd();
    //
    fout->Write();

    fout->Close();

}
