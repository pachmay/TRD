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

/// \file PadCalibCCDBBuilder.cxx
/// \brief Krypton calibration - class to smoothen and inter/extrapolate gain in chambers and calculate normalized ADC gain per chamber
/// \author Jana Crkovska

#include "PadCalibCCDBBuilder.h"
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

// using namespace o2::trd::constants;

// namespace o2::trd
// {

/**
 * COMMENT the next 3 lines to run in O2!
 * */
int NCOLUMN = 144;
int NROWC0 = 12;
int NROWC1 = 16;

/**
 * Check if the pad is a part of a small cluster that can be contained within an area of NxN pads
 * where N < areaMustBeContainedWithin.
 * 
 * @param hDet A pointer to the TH2F object representing the detector.
 * @param coordinates A vector of integers representing the coordinates of the two pads being compared.
 *                     Must have size 4 with column and row of pad1 and pad2.
 * @param upperLimit A float representing the upper limit for the pad gain. Pads with a gain above this limit
 *                   are considered hot pads and will be checked for isolation.
 * @param areaMustBeContainedWithin An integer representing the maximum size of the hot pad area to be considered.
 *                                  Hot pads with an area greater than or equal to this will not be replaced.
 * 
 * !!! TO DO !!!
 * now the code in solatedHotPadsContainmentSize does not allow to set the size of the hot area
 * maximum is 3x3, from 4 up it will not consider anything as a hot pad area
 * so need to add a parameter to isolatedHotPadsContainmentSize, which would control this
 */

void PadCalibCCDBBuilder::checkIfIsolatedHotPadCandidate(TH2F* hDet, std::vector<int> coordinates, float upperLimit, int areaMustBeContainedWithin)
{ 
  auto numberOfCoordinates = coordinates.size();  /* must be a vector of size 4 with x and y of pad1 and pad2 */
  if (numberOfCoordinates != 4) {
    std::cerr << "Invalid size of the coordinates vector!" << std::endl;
    return;
  }

  float averageGain = computeDetectorAverage(hDet); /* computed the average gain in the detector to compare gain in each pad against */

  // Check the first pad for hotness and isolation
  if (hDet->GetBinContent(coordinates[0], coordinates[1]) > upperLimit * averageGain) {
    int sizeOfHotArea = isolatedHotPadsContainmentSize(hDet, coordinates[0], coordinates[1]);
    if (sizeOfHotArea > 0 && sizeOfHotArea < areaMustBeContainedWithin) {
      replaceIsolatedHotPads(hDet, coordinates[0], coordinates[1], sizeOfHotArea);
    }
  }
  // Check the second pad for hotness and isolation 
  else if (hDet->GetBinContent(coordinates[2], coordinates[3]) > upperLimit * averageGain) {
    int sizeOfHotArea = isolatedHotPadsContainmentSize(hDet, coordinates[2], coordinates[3]);
    if (sizeOfHotArea > 0 && sizeOfHotArea < areaMustBeContainedWithin) {
      replaceIsolatedHotPads(hDet, coordinates[2], coordinates[3], sizeOfHotArea);
    }
  }
  // Check for the case where the first pad has zero gain and the second pad does not
  /* N.B. 
   * I only added check for when the first pad is empty,
   * because the code identifies these sus pads by sweeping from left to right or bottom to top. 
   * So if a potential hot area is detected, this will happen on its left or bottom edge
   * and the coordinates of the empty bin will always be stored first.
   * If the code decides that indeed this is a hot pad area, it will remove it entirely, 
   * i.e., no need to check again on its right/top edge.
   */ 
  else if (hDet->GetBinContent(coordinates[0], coordinates[1]) == 0 && hDet->GetBinContent(coordinates[2], coordinates[3]) != 0) {
    int sizeOfHotArea = isolatedHotPadsContainmentSize(hDet, coordinates[2], coordinates[3]);
    if (sizeOfHotArea > 0 && sizeOfHotArea < areaMustBeContainedWithin) {
      replaceIsolatedHotPads(hDet, coordinates[2], coordinates[3], sizeOfHotArea);
    }
  }
}

/**
 * Takes a pair of suspect pads
 * and computes the distance between each pad and the center of the detector'
 * If the pad at a closer distance has lower gain than the one further,
 * this closer pad will be set to 0.
 * 
 * 
 * DISCALIMER
 * This method is no longer used in the code!
 * This premise that all maps should have a gain smoothly increasing towards the center is a priori false
 * so it makes no sense to keep using the method that calls this one to replace the bin closer to the center.
 * The code was however left if for some reason this could prove useful in the future.
 */
void PadCalibCCDBBuilder::checkIfSmallerCloserToCenter(TH2F* hDet, std::vector<int> coordinates, float allowedDifference)
{

  if (hDet->GetBinContent(coordinates[0], coordinates[1]) == 0 || hDet->GetBinContent(coordinates[2], coordinates[3]) == 0) {
    return;
  }

  float xCenter = hDet->GetNbinsX() / 2.;
  float yCenter = hDet->GetNbinsY() / 2.;

  std::vector<float> vCenter, vPad1, vPad2;
  vCenter.push_back(hDet->GetNbinsX() / 2.);
  vCenter.push_back(hDet->GetNbinsY() / 2.);

  vPad1.push_back(hDet->GetXaxis()->GetBinCenter(coordinates[0]));
  vPad1.push_back(hDet->GetYaxis()->GetBinCenter(coordinates[1]));

  vPad2.push_back(hDet->GetXaxis()->GetBinCenter(coordinates[2]));
  vPad2.push_back(hDet->GetYaxis()->GetBinCenter(coordinates[3]));

  float dist1 = computeDistance(vPad1, vCenter);
  float dist2 = computeDistance(vPad2, vCenter);

  if ((dist1 < dist2) && ((hDet->GetBinContent(coordinates[2], coordinates[3]) - hDet->GetBinContent(coordinates[0], coordinates[1])) > allowedDifference)) {
    replacePadCloserToCenter(hDet, coordinates[0], coordinates[1]);
  }

  if ((dist1 > dist2) && ((hDet->GetBinContent(coordinates[0], coordinates[1]) - hDet->GetBinContent(coordinates[2], coordinates[3])) > allowedDifference)) {
    replacePadCloserToCenter(hDet, coordinates[2], coordinates[3]);
  }
}

/**
 * Compare the gain of two neighboring pads in a TH2F histogram, and return
 * the coordinates of the pads if the gain difference is greater than the
 * specified threshold.
 * 
 * N.B. Here the coordinates do not give column/row as a TRD coordinate
 * but rather a x-/y-index of the corresponding bin within the TH2F.
 *
 * @param hDet The TH2F histogram to analyze.
 * @param column The X index of the first pad. 
 * @param row The Y index of the first pad.
 * @param shiftcolumn The X shift of the second pad (should be 0 or 1).
 * @param shiftrow The Y shift of the second pad (should be 0 or 1).
 * @param allowedDifference The maximum allowed difference between the pad gains.
 * @return A vector containing the coordinates (x+y-indices) of the pads if their gains differ
 * by more than the allowed difference, or {-1, -1, -1, -1} if the pads are invalid (e.g., at edges)
 * or their gains do not differ enough.
 */
std::vector<int> PadCalibCCDBBuilder::compareGain(TH2F* hDet, int column, int row, int shiftcolumn, int shiftrow, float allowedDifference)
{
  // Initialize a vector to hold the coordinates of the two neigboring pads
  std::vector<int> coordinates = {-1, -1, -1, -1};

  // Get the maximum number of columns and rows in the histogram
  int colMax = hDet->GetNbinsX();
  int rowMax = hDet->GetNbinsY();

  // Check that the shift is along only one axis and only by one pad
  if (!(shiftcolumn == 1 && shiftrow == 0) && !(shiftcolumn == 0 && shiftrow == 1)) {
    return coordinates;
  }

  // Check that the pad is valid
  if ((column >= colMax) || (row >= rowMax)) {
    return coordinates;
  }

  // Do not compare with overflow
  if ((column == colMax && shiftcolumn == 1) || (row == rowMax && shiftrow == 1)) {
    return coordinates; // do not compare with overflow
  }

  // Get the gain values of the two pads to be compared
  float gain1 = hDet->GetBinContent(column, row);
  float gain2 = hDet->GetBinContent(column + shiftcolumn, row + shiftrow);

  // If the gain values differ by more than the allowed difference, or one gain is 0 while the other is not,
  // save the coordinates of the two pads in the coordinates vector
  if ((gain1 == 0 && gain2 > 0) || (gain1 > 0 && gain2 == 0) || (abs(gain1 - gain2) > allowedDifference)) {
    coordinates[0] = column;
    coordinates[1] = row;
    coordinates[2] = column + shiftcolumn;
    coordinates[3] = row + shiftrow;
  }

  // Return the coordinates vector, which may be {-1, -1, -1, -1} if the two pads have similar enough gain
  return coordinates;
}

/**
 * Computes the average gain value over filled cells of a 2D histogram representing a detector.
 * Cells are accessed through their TRD column and row coordinates, not bin indices.
 * The average is computed over absolute values as interpolated/extrapolated cells may have negative values.
 * 
 * @param hDet The 2D histogram representing the detector.
 * @return The average gain value over filled cells of the detector.
 */
float PadCalibCCDBBuilder::computeDetectorAverage(TH2F* hDet)
{ 
  float average = 0.;
  int nBinsUsed = 0;

  // loop over all bins in the detector map
  for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
    for (int icol = 0; icol < hDet->GetNbinsX(); icol++) {
      // get the current gain and skip the bin if it's not filled
      float currentGain = TMath::Abs(hDet->GetBinContent(icol + 1, irow + 1));
      if (currentGain == 0) {
        continue;
      }
      // add the current gain to the average and count the number of filled bins
      average += currentGain;
      nBinsUsed++;
    }
  }
  // compute the average gain over filled bins
  average /= nBinsUsed;

  return average;
}

/**
 * Computes the distance between two pads 
 * 
 * @param pad1 The coordinates of the first pad as a vector of floats
 * @param pad2 The coordinates of the second pad as a vector of floats
 * @return float The distance between the two pads
 */
float PadCalibCCDBBuilder::computeDistance(std::vector<float> pad1, std::vector<float> pad2)
{
  float distance = -1.;

  // Check that the input vectors have the same dimension
  auto numberOfCoordinates = pad1.size();
  if (numberOfCoordinates != pad2.size()) {
    std::cerr << "Something fishy with the pad coordinates!" << std::endl;
    return distance;
  }

  // Compute the Euclidean distance
  for (int i = 0; i < numberOfCoordinates; i++) {
    distance += TMath::Power(pad1[i] - pad2[i], 2);
  }
  distance = TMath::Sqrt(distance);

  return distance;
}

/**
 * Clones the filled 2D map and computes the average gain.
 * Replaces each bin content with the ratio of the gain over the average.
 *
 * @param hDet  a pointer to a TH2F histogram containing the detector information.
 * @param sNewName  a TString with the new name for the normalized map. If an empty string is given,
 *                  a default name is generated by appending "_normalized" to the original histogram's name.
 * @return  a pointer to a newly created TH2F histogram with normalized values.
 */
TH2F* PadCalibCCDBBuilder::createNormalizedMap(TH2F* hDet, TString sNewName)
{ 
  // If no name is provided, create a default name
  if (sNewName == "") { 
    sNewName = hDet->GetName();
    sNewName += "_normalized";
  }

  // Clone the input histogram
  TH2F* hDetNormalized = (TH2F*)hDet->Clone(sNewName.Data());

  // Compute the average gain over all filled bins
  float average = computeDetectorAverage(hDet);

  // Loop over all bins in the histogram
  for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
    for (int icol = 0; icol < hDet->GetNbinsX(); icol++) {
      // Get the current gain for this bin
      float currentGain = hDet->GetBinContent(hDet->GetXaxis()->FindBin(icol), hDet->GetYaxis()->FindBin(irow));
      // Calculate the new gain as the ratio of current gain over average gain
      float newGain = currentGain / average;
      // Set the new gain as the content of the bin in the normalized histogram
      hDetNormalized->SetBinContent(hDetNormalized->GetXaxis()->FindBin(icol), hDetNormalized->GetYaxis()->FindBin(irow), newGain);
    }
  }

  // Return the normalized histogram
  return hDetNormalized;
}

/**
 * Replaces a gain in an empty pad with a new value.
 * This function changes the input map! 
 * 
 * @param hDet pointer to the map 
 * @param column X index of the bin in the map 
 * @param row Y index of the bin in the map
 * @param newGain the new value to set in the bin  
 */
void PadCalibCCDBBuilder::fillInTheGap(TH2F* hDet, int column, int row, float newGain)
{ 
  // get the current value in our pad of interest
  float currentGain = hDet->GetBinContent(column + 1, row + 1);

  if (currentGain != 0) {
    return; // make sure we don't replace existing non-zero gain, just fill in gaps
  }
  // set the pad to the new value
  hDet->SetBinContent(column + 1, row + 1, newGain);
}

/**
 * Clones the map and fills the clone. Fills any map with at least 1 hit. 
 * 
 * @param hDet pointer to the input map 
 * @param sNewName name of the output map 
 * @param nbuffer size of the neighborhood to consider when calculating the average gain 
 * @return TH2F* pointer to the filled map 
 */
TH2F* PadCalibCCDBBuilder::fillTheMap(TH2F* hDet, TString sNewName, int nbuffer)
{ 
  // Clone the input map 
  if (sNewName == "") { // create default name
    sNewName = hDet->GetName();
    sNewName += "_filled";
  }
  TH2F* hDetFilled = (TH2F*)hDet->Clone(sNewName);

  // Create a temporary copy of the input map 
  TH2F* hDetTemp = (TH2F*)hDet->Clone("hDetTemp");

  // Find empty bins in the temporary map 
  std::vector<std::vector<int>> emptyBinsColRow = findEmpty(hDetTemp);
  
  int firstFilledX = hDetFilled->FindFirstBinAbove();
  int lastFilledX = hDetFilled->FindLastBinAbove();
  int nBinsX = hDetFilled->GetNbinsX();

  int firstFilledY = hDetFilled->FindFirstBinAbove(0, 2);
  int lastFilledY = hDetFilled->FindLastBinAbove(0, 2);
  int nBinsY = hDetFilled->GetNbinsY();

  // Loop over empty bins of the clone and fill them with the average gain
  // calculated from the gain in neighboring bins in the "Temp" map
  while (!emptyBinsColRow.empty()) {

    // Loop over each empty bin 
    // for (const auto& ibin : emptyBinsColRow) {
    for (int ibin = 0; ibin < emptyBinsColRow.size(); ibin++) {

      int flippedCoordinateX = nBinsX - hDetFilled->GetXaxis()->FindBin(emptyBinsColRow[ibin][0]) + 1;
      int flippedCoordinateY = nBinsY - hDetFilled->GetYaxis()->FindBin(emptyBinsColRow[ibin][1]) + 1;
      // Get mirrored gain from the filled map
      float mirroredGain = hDetFilled->GetBinContent(flippedCoordinateX, hDetFilled->GetYaxis()->FindBin(emptyBinsColRow[ibin][1]));
      
      // If mirrored gain is not found mirroring horizontally, check vertical mirror
      if (mirroredGain == 0) {
        mirroredGain = hDetFilled->GetBinContent(hDetFilled->GetXaxis()->FindBin(emptyBinsColRow[ibin][0]), flippedCoordinateY);
      }
      
      // If mirrored gain is still not found, calculate the average gain from neighboring bins
      if (mirroredGain == 0) {
        mirroredGain = getAverageFromNeighbors(hDetTemp, emptyBinsColRow[ibin][0], emptyBinsColRow[ibin][1], nbuffer);
      }
      
      // Determine factor based on the sign of mirroredGain
      // This is necessary because already "artificially" filled bins will have a gain < 0
      // and you do not want to flip this back to >0
      float factor = mirroredGain > 0 ? -1 : 1;

      fillInTheGap(hDetFilled, emptyBinsColRow[ibin][0], emptyBinsColRow[ibin][1], factor * mirroredGain);
    }
    //
    auto nEmptyPrevious = emptyBinsColRow.size();
    emptyBinsColRow.clear();

    hDetTemp = (TH2F*)hDetFilled->Clone("hDetTemp");
    emptyBinsColRow = findEmpty(hDetTemp);

    if (emptyBinsColRow.size() == nEmptyPrevious) {
      break; // will break out of the loop if no more empty pads can be filled
    }
  } // will continue the loop till all bins are filled
    // or until the algo is not able to fill any of the remaining empty pads
    //
    // Note: Likely you could replace the break condition to
    // "if (emptyBinsColRow.size() == 0) break;"
    // When I was writing the code, I needed (for reasons) a safeguard
    // that would make sure I always can leave the loop at some point.
    // Feel free to test and improve.

  delete hDetTemp;

  return hDetFilled;
}

/**
 * Finds the bin indices of all empty bins in the given 2D histogram
 * 
 * @param hDetectorMap a pointer to the input 2D histogram
 * @return a vector of vectors containing the bin indices of all empty bins
 */ 
std::vector<std::vector<int>> PadCalibCCDBBuilder::findEmpty(TH2F* hDetectorMap)
{ 
  std::vector<std::vector<int>> emptyBins;

  // Loop over all bins in the 2D histogram
  for (int irow = 0; irow < hDetectorMap->GetNbinsY(); irow++) {
    for (int icolumn = 0; icolumn < hDetectorMap->GetNbinsX(); icolumn++) {
      // Check if the bin is empty (i.e. has a gain of 0)
      float gain = hDetectorMap->GetBinContent(icolumn + 1, irow + 1);
      if (gain == 0) {
        // If the bin is empty, store its bin indices in a vector and add it to the list of empty bins
        std::vector<int> coordinates; 
        coordinates.push_back(icolumn);
        coordinates.push_back(irow);
        emptyBins.push_back(coordinates);
      }
    } // loop over columns
  }   // loop over rows
  return emptyBins;
}

/**
 * Finds bins that have 1+ neighbor with significantly different content
 * @param hDet Pointer to the TH2F histogram representing the detector map.
 * @param allowedDifference The maximum allowed difference between the pad gains.
 * @return A vector of vectors with bin indices of suspicious pad pairs,
 * the vector of indices is in the form 
 * (x-index-of-first-pad, y-index-of-first-pad, x-index-of-second-pad, y-index-of-second-pad).
 */
std::vector<std::vector<int>> PadCalibCCDBBuilder::findInhomogeneities(TH2F* hDet, float allowedDifference)
{ 
  std::vector<std::vector<int>> suspiciousBinPairs;

  // check for edges - this operation goes before filling
  // edges along X:
  int xFirst = hDet->FindFirstBinAbove();
  int xLast = hDet->FindLastBinAbove(); // stop the process at penultimate vs ultimate row/column
  // // edges along Y:
  int yFirst = hDet->FindFirstBinAbove(0, 2);
  int yLast = hDet->FindLastBinAbove(0, 2); // stop the process at penultimate vs ultimate row/column


  for (int irow = 0; irow < hDet->GetNbinsY(); irow++) { // loop over all columns in the chamber

    int thisRow = hDet->GetXaxis()->FindBin(irow); //This must be here because the bin labels 
    // (which correspond to the row/column as it comes from the TRD) go from 0 till NCOLUMN-1/NROW-1
    // so this can just be accessed via looping from 0 up to whatever the bin number is.
    // But we want bin indices, which start from 1 and go to NCOLUMN/NROW.

    for (int icolumn = 0; icolumn < hDet->GetNbinsX(); icolumn++) { // loop over all rows in the chamber

      int thisColumn = hDet->GetXaxis()->FindBin(icolumn); // see above

      std::vector<int> pair = compareGain(hDet, thisColumn, thisRow, 1, 0, allowedDifference); // compare the pad with its right-hand-side neighbor (shift aliong X)
      if (!(std::any_of(pair.begin(), pair.end(), [](int i) { return i == -1; }))) {
        suspiciousBinPairs.push_back(pair); // any vector of indices, which DOES NOT have all 4 elements == -1, is stored
      }
      pair.clear();

      pair = compareGain(hDet, thisColumn, thisRow, 0, 1, allowedDifference); // compare the pad with its upper neighbor (shift along Y)
      if (!(std::any_of(pair.begin(), pair.end(), [](int i) { return i == -1; }))) {
        suspiciousBinPairs.push_back(pair);
      }

    } // loop over columns

  } // loop over rows

  return suspiciousBinPairs;
}

/**
 * Takes bin (column,row) from hDet map and averages the gain in its (up to) 8 neighbours
 *
 * @param hDet: pointer to TH2F histogram of detector
 * @param column: column index of bin to average neighbors of
 * @param row: row index of bin to average neighbors of
 * @param nbuffer: size of buffer around bin to average neighbors of
 * @return average gain of neighbors
 */
float PadCalibCCDBBuilder::getAverageFromNeighbors(TH2F* hDet, int column, int row, int nbuffer)
{ 
  float average = 0.;
  std::vector<float> gainNeighbor;  // vector to store the values of neighboring bins
  int offset = nbuffer / 2;

  for (int irow = 0; irow < nbuffer; irow++) {

    int rowNeighbor = (row - offset) + irow;
    if (rowNeighbor < 0 || rowNeighbor >= hDet->GetNbinsY()) {
      continue; // skip out-of-bounds neighbors (avoids under and overflow)
    }
    for (int icol = 0; icol < nbuffer; icol++) {
      if (icol == 1 && irow == 1) {
        continue; // exclude self
      }
      int colNeighbor = (column - offset) + icol;
      if (colNeighbor < 0 || colNeighbor >= hDet->GetNbinsX()) {
        continue; // skip out-of-bounds neighbors (avoids under and overflow)
      }
      float tempGain = TMath::Abs(hDet->GetBinContent(hDet->GetXaxis()->FindBin(colNeighbor), hDet->GetYaxis()->FindBin(rowNeighbor)));
      if (tempGain <= 0) {
        continue; // skip empty or negative bins
      }
      gainNeighbor.push_back(tempGain);
    }
  }
  auto numberOfNeighbors = gainNeighbor.size();
  if (numberOfNeighbors < 1) {
    return 0; // if there are no neighboring bins, return 0
  }

  average = std::accumulate(gainNeighbor.begin(), gainNeighbor.end(), decltype(gainNeighbor)::value_type(0.0)) / numberOfNeighbors;

  return average; // return the average gain of the neighboring bins
}

/**
 * Creates a TH2F map of a detector nDet, fills it by looping over a TTree and
 * setting ADC gain to its corresponding bin in the TH2F. Limits the accepted
 * ADC gain by setting mingain and maxgain (default range 0-10k).
 *
 * @param tree Pointer to TTree with detector data.
 * @param nDet Integer specifying the detector number to create a map for.
 * @param mingain Float specifying the minimum ADC gain to accept (default 0).
 * @param maxgain Float specifying the maximum ADC gain to accept (default 10000).
 * @param sDetName TString specifying the name of the TH2F to create (default "hDetN").
 * @return Pointer to the created TH2F.
 */
TH2F* PadCalibCCDBBuilder::getDetectorMap(TTree* tree, int nDet, float mingain, float maxgain, TString sDetName)
{ 
  // checks if a name has been provided for the TH2F and sets it accordingly
  if (sDetName == "") {
    sDetName = Form("hDet%i", nDet);
  }

  // sets the number of rows in the TRD for the specified detector
  // remember there is 12 rows in stack 2
  int nTrdRows = NROWC1;
  int detInSupermodule = nDet % 30; 
  if (detInSupermodule >= 12 && detInSupermodule <= 17) {
    nTrdRows = NROWC0;
  }
  
  // creates the 2D histogram to be filled
  TH2F* hDetector = new TH2F(sDetName.Data(), sDetName.Data(), NCOLUMN, 0, NCOLUMN, nTrdRows, 0, nTrdRows);
  
  // sets the branches of the TTree
  setTreeBranches(tree);
  
  // loops over the TTree and fills the histogram
  int nentries = tree->GetEntries();
  for (int ientry = 0; ientry < nentries; ientry++) {
    tree->GetEntry(ientry);
    if ((int)mDet != nDet) {
      continue;
    }
    if (mChi <= mChiMin || mAmp <= mAmpMin || mSgm <= mSgmMin || mSgm > mSgmMax) {
      continue;
    }
    if (mAdc < mingain || mAdc > maxgain) {
      continue;
    }
    hDetector->SetBinContent(hDetector->GetXaxis()->FindBin(mCol), hDetector->GetYaxis()->FindBin(mRow), mAdc);
  }
  return hDetector;
}

/**
 *  Checks if an area can be considered isolated, i.e., if the local gain
 *  is much different in the area wrt its surroundings.
 * 
 * @param hDet The map to be investigated.
 * @param column x-index of the first pad of a potential isolated are (bottom left corner)
 * @param row y-index of the first pad of a potential isolated area
 * @param matrixSize the size of the potential isolated area (must be within a square of an area matrixSize*matrixSize)
 * @return a bool isIsolated:
 *  - kTRUE if the cluster will be considered as isolated area of hot pads
 *  - kFALSE if not
 * 
 *  isolated = EITHER no filled pads around
 *  OR it means that the surrounding pads have a mean gain close to the average gain of the detector
 *  so in principle this should also find hot pads inside a filled map
 * 
 *  Takes into account the pads close to edge, where the number of surrounding pads is modified. 
 */
bool PadCalibCCDBBuilder::isHotAreaIsolated(TH2F* hDet, int column, int row, int matrixSize)
{
  bool isIsolated = kFALSE;
  float averageGain = computeDetectorAverage(hDet);
  int nSurroundingNotHot = 0; /* this variable will hold the number of EMPTY surrounding pads */
  int nMaxHot = TMath::Power(matrixSize + 2, 2) - TMath::Power(matrixSize, 2); /* computes the number of pads that have to surround an isolated cluster of size matrixSize*matrixSize */
  float averageAround = 0.; /*  variable to hold the mean gain of surrounding pads that are filled */
  int nUsedBins = 0;  /* variable to keep how many surrounding pads are NOT empty*/

  int nHotOffset = 0; // offsets the nMaxHot criterion for cells at edges
  if ((column == 1 || column == hDet->GetNbinsX()) && row > 1 && row < hDet->GetNbinsY()) {
    nHotOffset = matrixSize + 2;
  } else if ((row == 1 || row == hDet->GetNbinsY()) && column > 1 && column < hDet->GetNbinsX()) {
    nHotOffset = matrixSize + 2;
  } else if ((column == 1 || column == hDet->GetNbinsX()) && (row == 1 || row == hDet->GetNbinsY())) {
    nHotOffset = nMaxHot / 2 + 1;
  }

  for (int i = 0; i < matrixSize + 2; i++) {
    int icol = i - 1;
    for (int j = 0; j < matrixSize + 2; j++) {
      int jrow = j - 1;
      if ((-1 < icol) && (icol < matrixSize) && (-1 < jrow) && (jrow < matrixSize)) {
        continue;
      }
      float temp = TMath::Abs(hDet->GetBinContent(column + icol, row + jrow));
      if (temp != 0) { /* if surrounding pad is not empty than its gain is considered for the mean gain over surrounding pads */
        nUsedBins++;
        averageAround += temp;
      } else if (temp == 0) { /* counts how many surrounding pads are empty */
        nSurroundingNotHot++;
      }
    }
  }

  if (nUsedBins > 0) {
    averageAround /= nUsedBins; /* computes the mean gain over surrounding pads that are FILLED (empties not counted!) */
  }
  if (averageAround < 2 * averageGain && nSurroundingNotHot <= nHotOffset) { /* checks that the surrounding filled pads are not hot, i.e., their mean gain does not exceed double the average gain per the whole detector */
  //
  /* TO DO : It would make a lots of sense to make this criterion more flexible so that the user can impose how much above the average the pads can be. Didn't have time to implement, feel free to do so if you like. */
  //
    isIsolated = kTRUE;
  } else if (nSurroundingNotHot == nMaxHot) { /* If all surrounding are empty than yes this is a hot cluster. */
    isIsolated = kTRUE;
  }

  return isIsolated;
}

/**
  * Returns the size of the area within which the cluster of hot pads can be contained
  * the area is NxN
  * 
  * @param hDet The map to be investigated.
  * @param column x-index of the first pad of a potential isolated are (bottom left corner)
  * @param row y-index of the first pad of a potential isolated area
  * @param areaMustBeContainedWithin not used in this version (smh), was intended to control the size of the maximum checked hot area
  * 
  * @return N if N <= 3
  * @return -1 if N > 4 (we do not consider those contained clusters)
  */
int PadCalibCCDBBuilder::isolatedHotPadsContainmentSize(TH2F* hDet, int column, int row)
{ 

 // TO DO : add a parameter to control the size of the hot area

  int nsize = 0;
  bool isContained = kFALSE;
  while (!isContained) {
    nsize++;
    isContained = isHotAreaIsolated(hDet, column, row, nsize); // returns kTRIE if hot pads within a square of nsize*nsize
    if (nsize == 4) {
      return -1;
    }
  }
  return nsize;
}
/**
 * Sets all cells in an empty normalized map to a given value.
 * By default all cells are set to -1.
 * 
 * @param hDet The map to be filled.
 * @param valueToSet The value you want to set the gain in the map to.
 */ 
void PadCalibCCDBBuilder::populateEmptyNormalizedMap(TH2F* hDet, float valueToSet)
{ 
  if (!hDet || hDet->GetEntries() != 0) {
    std::cout << "Histogram does not exist or is not empty!";
    return;
  }
  for (int icol = 0; icol < hDet->GetNbinsX(); icol++) {
    for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
      hDet->SetBinContent(icol + 1, irow + 1, valueToSet);
    }
  }
}

/**
 * Sets all cells in edges (along Y axis = rows) to 0
 * by default, the two pads closest to the edge are removed (i.e., the entire columns 0, 1 and 142, 143 are set to 0).
 * The number of pads to be removed can be set to another value.
 * @param hDet a pointer to the 2D histogram to modify
 * @param nsize the number of pads to remove from each edge (default 2)
 */
void PadCalibCCDBBuilder::removeEdges(TH2F* hDet, int nsize)
{ 
  if (!hDet || hDet->GetEntries() == 0) {
    return;
  }
  for (int icol = 0; icol < nsize; icol++) {
    // Loop over all the rows in the histogram
    for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
      // Set the bin content to zero for the first nsize columns
      hDet->SetBinContent(icol + 1, irow + 1, 0); 
      // Set the bin content to zero for the last nsize columns
      hDet->SetBinContent(hDet->GetNbinsX() - icol, irow, 0); 
    }
  }
}

/**
 * Removes cells with gain values that are outside of a given range.
 *
 * @param hDet Pointer to the TH2F histogram representing the detector map.
 * @param upperLimit The upper limit for gain values (default is 2, in multiples of average gain).
 * @param lowerLimit The lower limit for gain values (default is 0.5, in multiples of average gain).
 */
void PadCalibCCDBBuilder::removeExtremePads(TH2F* hDet, float upperLimit, float lowerLimit)
{ 
  // If the input histogram is null or empty, return immediately
  if (!hDet || hDet->GetEntries() == 0) {
    return;
  }

  // Calculate the average gain value for the detector map
  float average = computeDetectorAverage(hDet);

  // Loop over all cells in the detector map
  for (int irow = 0; irow < hDet->GetNbinsY(); irow++) {
    for (int icol = 0; icol < hDet->GetNbinsX(); icol++) {
      float value = hDet->GetBinContent(icol + 1, irow + 1);
      // If the cell's gain value is above the upper limit or below the lower limit, set it to 0
      if (value > upperLimit * average || value < lowerLimit * average) {
        hDet->SetBinContent(icol + 1, irow + 1, 0);
      }
    }
  }
}

/**
 * this method has been originally conceived to set to 0 those pads
 * in the pairs of pads label as problematic,
 * which are closer to the detector center as compared to the other pad in the pair
 * 
 * what it really does is to set the pad with the given BIN corrdinates to 0
 * (so its coordinates within X:Y of the TH2*, not its row and column corrdinates)
 *
 * 
 * DISCALIMER
 * This method is no longer used in the code!
 * This premise that all maps should have a gain smoothly increasing towards the center is false
 * so it makes no sense to keep using the method that calls this one to replace the bin closer to the center.
 * The code was, however left, if for some reason this could prove useful in the future.
 */
void PadCalibCCDBBuilder::replacePadCloserToCenter(TH2F* hDet, int xcloser, int ycloser)
{

  if (!hDet || hDet->GetEntries() == 0) {
    std::cerr << "invalid histogram!" << std::endl;
    return;
  }

  float newGain = 0.;
  hDet->SetBinContent(xcloser, ycloser, newGain);
}

/**
 * Replaces a square block of pads in a detector histogram with zero gain values.
 *
 * @param hDet - Pointer to the 2D histogram representing the detector
 * @param column - Column index of the top-left pad in the block to be replaced
 * @param row - Row index of the top-left pad in the block to be replaced
 * @param nsize - The size of the square block to be replaced (number of rows and columns)
 */
void PadCalibCCDBBuilder::replaceIsolatedHotPads(TH2F* hDet, int column, int row, int nsize)
{
  // Loop over the NxN area of hot pads and set their bin contents to zero
  // where N == nsize
  for (int jrow = 0; jrow < nsize; jrow++) {
    for (int icol = 0; icol < nsize; icol++) {
      hDet->SetBinContent(column + icol, row + jrow, 0);
    }
  }
}

/**
 * Sets the branch addresses for a given TTree object, so that data from the tree can be read into the PadCalibCCDBBuilder object.
 * The branch addresses are set for the following data members of the PadCalibCCDBBuilder object:
 *     mDet: an integer representing the detector ID                // this is the chamber number as it comes from the detector!
 *     mCol: an integer representing the column of the pad          // this is the column number as it comes from the detector!
 *     mRow: an integer representing the row of the pad             // this is the row number as it comes from the detector!
 *     mAdc: a float representing the mean ADC value in the pad     // from the fit to the ADC spectrum
 *     mChi: a float representing the chi2 value in the pad         // from the fit to the ADC spectrum
 *     mSgm: a float representing the sigma value in the pad        // from the fit to the ADC spectrum
 *     mAmp: a float representing the amplitude value in the pad    // from the fit to the ADC spectrum
 */
void PadCalibCCDBBuilder::setTreeBranches(TTree* tree)
{
  tree->SetBranchAddress("det", &mDet);
  tree->SetBranchAddress("col", &mCol);
  tree->SetBranchAddress("row", &mRow);
  tree->SetBranchAddress("mean", &mAdc);
  tree->SetBranchAddress("chi2", &mChi);
  tree->SetBranchAddress("sigma", &mSgm);
  tree->SetBranchAddress("amplitude", &mAmp);
}

/**
 * Cleans up the map, i.e., removes found inhomogeneities between pads.
 * The map is homogeneous if neighboring pads have gain within @param allowedDifference.
 */ 
void PadCalibCCDBBuilder::smoothenTheDetector(TH2F* hDet, float allowedDifference)
{
  std::vector<std::vector<int>> SuspiciousBinPairs = findInhomogeneities(hDet, allowedDifference); // find pairs of pad with too much of a difference in gain
  auto numberOfPairs = SuspiciousBinPairs.size();
  if (numberOfPairs == 0) {
    return;
  }

  for (auto susPair : SuspiciousBinPairs) { 
    checkIfIsolatedHotPadCandidate(hDet, susPair);
  }
  for (auto element : SuspiciousBinPairs ) {
    if ( hDet->GetBinContent(element[0], element[1]) == 0 || hDet->GetBinContent(element[2], element[3]) == 0 ) continue;
        hDet->SetBinContent(element[0], element[1], 0);
        hDet->SetBinContent(element[2], element[3], 0);
  }
}

/**
 * Transforms the input histogram 'hDet' into a new histogram with absolute values and returns it.
 * If 'sName' is provided, sets it as the name of the new histogram, otherwise the default name is used.
 * 
 * Returns a new TH2F!
 * 
 * @param hDet - the input histogram to be transformed
 * @param sName - the name of the new histogram to be created (optional)
 * @return a pointer to the new histogram with absolute values
 */
TH2F* PadCalibCCDBBuilder::transformMapIntoAbsoluteValues(TH2F* hDet, TString sName)
{
  // if input histogram is null or empty, return null
  if (!hDet || hDet->GetEntries() == 0) {
    return nullptr;
  }
  // set a name for the new histo
  if (sName == "") {  // if no name provided, set the default name
    sName = hDet->GetName();
    sName += "_transformed";
  }

  int nCols = hDet->GetNbinsX();
  int nRows = hDet->GetNbinsY();

  // create a new histogram with the specified name and number of bins
  TH2F* hDetCopy = new TH2F(sName.Data(), sName.Data(), nCols, 0, nCols, nRows, 0, nRows);

  for (int icol = 0; icol < nCols; icol++) {
    for (int irow = 0; irow < nRows; irow++) {
      // set the bin content of the new histogram to the absolute value of the corresponding bin in the input histogram
      hDetCopy->SetBinContent(icol + 1, irow + 1, TMath::Abs(hDet->GetBinContent(icol + 1, irow + 1)));
    }
  }

  return hDetCopy;
}

// } // namespace o2::trd