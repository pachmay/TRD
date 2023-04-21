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

// using namespace o2::trd; /* uncomment to run with upload */

void CreateCCDBLocalGainFactor(TString sOpenFile)
{
  PadCalibCCDBBuilder* krCalibToCCDB = new PadCalibCCDBBuilder();

  TFile* file = TFile::Open(sOpenFile);
  TTree* tree = (TTree*)file->Get("nt_Krypton");
  // LocalGainFactor calObject; /* uncomment to run with upload */

  TH2F* hDetPlot;

  for (int idet = 0; idet < 1 ; idet++) {
    krCalibToCCDB->setSelection(0,0,0,500);
    TH2F* hDetDef = krCalibToCCDB->getDetectorMap(tree, idet);
    krCalibToCCDB->smoothenTheDetector(hDetDef, 500);
    TH2F* hDetFilled = krCalibToCCDB->fillTheMap(hDetDef);
    delete hDetDef;
    TH2F* hDet = krCalibToCCDB->createNormalizedMap(hDetFilled);
    delete hDetFilled;
    for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
      for (int icol = 0; icol < hDet->GetNbinsX(); icol++) {
        float relativeGain = hDet->GetBinContent(hDet->GetXaxis()->FindBin(icol), hDet->GetYaxis()->FindBin(irow));
        // calObject.setPadValue(idet, icol, irow, relativeGain); /* uncomment to run with upload */
      }
    }

    hDetPlot = krCalibToCCDB->transformMapIntoAbsoluteValues(hDet, "hDetPlot"); /* for when you want to check the results visually */
    // delete hDet; /* uncomment to run with upload */
    cout << idet << endl;
  }
  hDetPlot->Draw("colz");

  // delete tree; /* uncomment to run with upload */
  // delete file; /* uncomment to run with upload */

  // cout << calObject.getValue(0, 52, 13) << endl; /* if you want to check that calObject is indeed filled */

  return;

  /* 
  * uncomment all of the following code to upload to the CCDB 
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
