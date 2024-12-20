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


#include "TMath.h"


struct IndexTableUtils {
  int getZetaBin(float zeta);
  int getMultBin(float mult);
  int getBinIndex(float zeta, float mult);
  int mZetaBins = 20, mMultBins = 20;
  float minZeta = -10, maxZeta = 10;
  float minMult = 0., maxMult = 100;
};

inline int IndexTableUtils::getZetaBin(float zeta)
{
  float deltaZeta = (maxZeta - minZeta) / (mZetaBins);
  int bZeta = (zeta - minZeta) / deltaZeta; // bins recentered to 0
  return bZeta;
};

inline int IndexTableUtils::getMultBin(float mult)
{
  float deltaMult = (maxMult - minMult) / (mMultBins);
  int bMult = (mult - minMult) / deltaMult; // bin recentered to 0
  return bMult;
}

inline int IndexTableUtils::getBinIndex(float zeta, float mult)
{
  float deltaMult = (maxMult - minMult) / (mMultBins);
  float deltaZeta = (maxZeta - minZeta) / (mZetaBins);
  int bZeta = getZetaBin(zeta);
  int bMult = getMultBin(mult);
  // std::cout << "Zeta bin: " << bZeta << ", Mult bin: " << bMult << std::endl;
  return (bZeta >= mZetaBins || bMult >= mMultBins || bZeta < 0 || bMult < 0) ? mZetaBins * mMultBins : bZeta + mZetaBins * bMult;
}