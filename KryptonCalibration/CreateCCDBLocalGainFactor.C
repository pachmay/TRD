// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CreateCCDBLocalGainFactor.C
/// \brief x
/// \author Jana Crkovska

#include "PadCalibCCDBBuilder.h"

#include <ctime>
#include <sstream>
#include <string>

// using namespace o2::trd; /* UNCOMMENT to run with upload */
// using namespace o2::trd::constants; /* UNCOMMENT to run in O2, this enables the use of constants::MAXCHAMBER */

void CreateCCDBLocalGainFactor(TString sOpenFile)
{
  PadCalibCCDBBuilder* krCalibToCCDB = new PadCalibCCDBBuilder();

  TFile* file = TFile::Open(sOpenFile);
  TTree* tree = (TTree*)file->Get("nt_Krypton");
  // LocalGainFactor calObject; /* UNCOMMENT to run with upload */

  const int MAXCHAMBER = 540; /* COMMENT to run in O2 because this variable is already defined in o2::trd::constants */
  TH2F* hDet[MAXCHAMBER];
  TH2F* hDetPlot[MAXCHAMBER];

  for (int idet = 0; idet < MAXCHAMBER ; idet++) {
    krCalibToCCDB->setSelection(0,0,0,500);
    TH2F* hDetDef = krCalibToCCDB->getDetectorMap(tree, idet);

    if ( hDetDef->GetEntries() == 0 ) { /* this code will act only on EMPTY maps */
      krCalibToCCDB->populateEmptyNormalizedMap(hDetDef);
      hDet[idet] = (TH2F*) hDetDef->Clone( Form("%s_normalized", hDetDef->GetName()));
    } else { /* the following code will populate all the maps with at least 1 hit */
      krCalibToCCDB->smoothenTheDetector(hDetDef, 500);
      TH2F* hDetFilled = krCalibToCCDB->fillTheMap(hDetDef);
      hDet[idet] = krCalibToCCDB->createNormalizedMap(hDetFilled, Form("%s_normalized", hDetDef->GetName()) );
      delete hDetDef;
      delete hDetFilled;
    }

    for (int irow = 0; irow < hDet[idet]->GetNbinsY(); irow++) {
      for (int icol = 0; icol < hDet[idet]->GetNbinsX(); icol++) {
        float relativeGain = hDet[idet]->GetBinContent(hDet[idet]->GetXaxis()->FindBin(icol), hDet[idet]->GetYaxis()->FindBin(irow));
        // calObject.setPadValue(idet, icol, irow, relativeGain); /* UNCOMMENT to run with upload */
      }
    }

    hDetPlot[idet] = krCalibToCCDB->transformMapIntoAbsoluteValues(hDet[idet], Form("hDetSwapped_%i",idet) ); /* for when you want to check the results visually */
    // delete hDet; /* UNCOMMENT to run with upload */
    cout << idet << endl;
  }
  cout << "done" << endl;

  TFile* fout = new TFile("krypton_maps.root", "recreate");
  fout->mkdir("normalized");
  fout->mkdir("swapped");
  for (int idet = 0; idet < MAXCHAMBER ; idet++) {
    fout->cd("normalized");
    hDet[idet]->Write();
    fout->cd();
    fout->cd("swapped");
    hDetPlot[idet]->Write();
    fout->cd();
  }
  fout->Write();
  fout->Close();

  // cout << calObject.getValue(0, 52, 13) << endl; /* if you want to check that calObject is indeed filled */

  delete tree; 
  delete file; 

  return;

  /* 
  * UNCOMMENT all of the following code to upload to the CCDB 
  */

  // o2::ccdb::CcdbApi ccdb;
  // ccdb.init("http://ccdb-test.cern.ch:8080");
  // std::map<std::string, std::string> metadata;

  // time_t start_time = 1631311200; // Fri Sep 10 2021 22:00:00 GMT+0000
  // time_t end_time = 2208985200;   // Sat Dec 31 2039 23:00:00 GMT+0000
  // auto timeStamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::from_time_t(start_time).time_since_epoch()).count();
  // auto timeStampEnd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::from_time_t(end_time).time_since_epoch()).count();
  // ccdb.storeAsTFileAny(&calObject, "TRD/Calib/LocalGainFactor", metadata, timeStamp, timeStampEnd);
}
