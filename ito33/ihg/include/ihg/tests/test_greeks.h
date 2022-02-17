/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/test_greeks.h
// Purpose:     Compare computed greeks  wih financial level
//                greeks done by shifting
// Created:     March 31, 2005
// RCS-ID:      $Id: test_greeks.h,v 1.9 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_TESTS_GREEKS_H_
#define _IHG_TESTS_GREEKS_H_

#include "ito33/beforestd.h"
#include <list>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/ihg/modeloutput.h"

#include "ito33/xml/write.h"

namespace ito33
{
 
namespace finance
{
  class QualityControl;
  class SessionData;
  class DerivativeVisitorGoodType;
  class Derivative;
  class ModelOutput;
}


namespace ihg
{

  class TheoreticalModel;

  
class TestGreeks
{
public:

  TestGreeks();

  /**
      Runs test for greeks based on a list of input filname

      @param fileList list of file name
   */
  void Run(const std::list<std::string>& fileList);

private:
  
  /**
      Gets the spot stored in sessiondata

      @param sFileName file to read
      @return spot spot value in sessiondata
   */
  double GetSpot(const std::string& sFileName);

  /**
      Stores the output data

      @param output Model output
   */
  void SetGreeksData(const finance::ModelOutput& output);
  
  /**
      The parameters for quality control

      @param pQualityControl parameters for quality control
   */
  void 
  SetQualityControl(const shared_ptr<finance::QualityControl>& pQualityControl)
  {
    m_pQualityControl = pQualityControl;
  }

  /**
      Checks that the finance level greek matches the computed greek.

      @param fileName file to test
   */
  void CheckGreekVega(const std::string& sFileName);
   
  /**
      Checks that the finance level rho greek matches the computed greek.

      @param fileName file to test
   */
  void CheckGreekRho(const std::string& sFileName);
   
  /**
      Checks that the finance level fx delta greek matches the computed greek.

      @param fileName file to test
   */
  void CheckGreekFXDelta(const std::string& sFileName);
   
  /**
      Checks the convergence of the FX delta greek.

      @param fileName file to test
   */
  void CheckGreekFXDeltaConv(const std::string& sFileName);

  /**
      Checks if the error is less than 5%

      @param dError error
      @param dOriginalGreek original greek value
      @param dGreek greek value computed by shifting
      @param sTestName name of test being run
      @param sFileName file being tested
   */
  void CheckResult(double dError, double dOriginalGreek, 
        double dGreek, std::string sTestName, std::string sFileName);
  
  /**
      Checks if the error is less than 5% for the greek surface
      at analysis date.

      @param pdOriginalGreeks Original greek values at analysis date
      @param pdOriginalPrices Original price values at analysis date
      @param pdOriginalSpots Original spot values at analysis date
      @param pdShiftedPrices Shifted price values at analysis date
      @param dInverseShift The inverse of the shift
      @param sTestName name of test being run
      @param sFileName file being tested
   */
  void CheckGreekSurfaceAtAnalysisDateResults
       ( const finance::Values& pdOriginalGreeks, 
         const finance::Values& pdOriginalPrices,
         const finance::Values& pdOriginalSpots,
         const finance::Values& pdShiftedPrices, 
         const finance::Values& pdShiftedSpots, 
         double dInverseShift,
         std::string sTestName, std::string sFileName );

  /// output file for test report
  std::ofstream m_OutputFile;

  /// xml tag
  ito33::XML::RootTag m_RootTag;
  
  /// Quality controle
  shared_ptr<finance::QualityControl> m_pQualityControl;

};

void ReadXMLFile(std::string filename,
                 shared_ptr<finance::SessionData>& sessionData,
                 shared_ptr<finance::Derivative>& pDerivative,
                 shared_ptr<ihg::TheoreticalModel>& pModel);

void PriceContract(finance::Derivative& pDerivative,
                   const shared_ptr<ito33::ihg::TheoreticalModel>& pModel,
                   finance::ModelOutput& output);

} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_TESTS_GREEKS_H_
