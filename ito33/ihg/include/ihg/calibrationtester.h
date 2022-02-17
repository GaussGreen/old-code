/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/calibrationtester.h
// Purpose:     Test class for calibration
// Author:      David
// Created:     2004/01/14
// RCS-ID:      $Id: calibrationtester.h,v 1.5 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/calibrationtester.h
    @brief Test class for calibration

    Test class for calibration
*/

#ifndef _IHG_CALIBRATIONTESTER_H_
#define _IHG_CALIBRATIONTESTER_H_

#include "ito33/array.h"

#include "ihg/basetester.h"

namespace ito33 
{
  
namespace finance
{
  enum OptionType;
  enum ExerciseType;
}

namespace ihg 
{

/**
  @brief Test class for calibration

  Helper functions for testing calibration in IHG

*/
class CalibrationTester : public ihg::BaseTester
{

public:
  CalibrationTester() 
  {
    m_nNbRuns = 0;
  }

  virtual ~CalibrationTester() { }

  /// Read the entire input file
  void ReadInputFile(std::string& sFilename);

  /// Read in call data specified by prices
  void ReadPriceCallData(std::ifstream& sIn);

  /// Read in call data specified by prices
  void ReadImpVolCallData(std::ifstream& sIn);

  /// Convert implied vols to prices
  void ConvertImpVolsToPrices();

  /// Do the calibration
  void Calibrate(size_t nNbTests);

  /// Base class virtual function. Doesn't make sense here
  void RunPricer()
  {
    ASSERT_MSG(false, "RunPricer function has no meaning for calibration");    
  }

  /// Read calibration data
  void ReadContractParams(std::ifstream& /* sIn */)
  {
  }

  void ReadContractParams(std::string& /* sInstrumentFile */)
  {
  }

protected:

  /// Call option data
  int m_iCallDataType;
  size_t m_nNbCalls;
  Array<int> m_piCallMaturities;
  Array<double> m_pdCallStrikes;
  Array<double> m_pdCallPrices;

  /// If the call data is provided via implied vols
  Array<double> m_pdImpliedVols;



}; // class CalibrationTester


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CALIBRATIONTESTER_H_
