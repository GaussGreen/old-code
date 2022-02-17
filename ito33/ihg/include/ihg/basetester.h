/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/basetester.h
// Purpose:     Base class for testing IHG projects
// Author:      David
// Created:     2003/12/10
// RCS-ID:      $Id: basetester.h,v 1.28 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/basetester.h
    @brief Base class for testing IHG projects

    Base class for testing IHG projects
*/

#ifndef _IHG_BASETESTER_H_
#define _IHG_BASETESTER_H_

#include "ito33/beforestd.h"
#include <fstream>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"

namespace ito33 
{

namespace finance
{
  class Dividends;
  class YieldCurve;
}

namespace numeric
{
  class MeshParams;
  class NumParams;
}

namespace ihg
{

class Volatility;
class HazardRate;


/**
  @brief Base class for testing IHG projects

  This base class provides helper functions for testing IHG projects.  In 
  particular, it provides functions for reading common parameters
  (yield curve, foreign curve, dividends, numerical params, etc) from test
  text files, and provides some functions for reporting test results.

  Functions are provided for each parameter type because:
  - different contracts may order data in different ways
  - if the data representation changes, only one simple function needs
    to be changed

  Derived classes are expected to implement functions for reading contract
  specific data, as well as functions for contract specific tests. Derived
  classes must also implement the pure virtual function ReadContractParams.
  Derived classes may need to reimplement the virtual functions for reading
  input files.
*/

class BaseTester
{

public:

  /// Constructor 
  BaseTester() 
  {
    std::cout.precision(10);

    // By default, compute nothing
    m_bComputeDelta = false;
    m_bComputeGamma = false;
    m_bComputeVega = false;
    m_bComputeRho = false;
    m_bComputeTheta = false;
    m_bComputeArrays = false;
    m_bComputeSurfaces = false;
  }

  virtual ~BaseTester() { }

  /** 
    Functions to read and set pricing data
  */

  /// Set the current spot. Used when reading from multiple files.
  void SetSpotSharePrice(double dSpot) { m_dSpot = dSpot; }

  /// Read and store the pricing date
  void ReadValuationDate(std::ifstream& sIn);
  void ReadValuationDate(std::string& sValuationDateFile);
  void SetValuationDate(Date ValuationDate) { m_ValuationDate = ValuationDate; }

  /// Read and construct the yield curve
  void ReadYieldCurve(std::ifstream& sIn);
  void ReadYieldCurve(std::string& sYieldFile);

  /// Read and construct the foreign curve
  void ReadForeignCurve(std::ifstream& sIn);
  void ReadForeignCurve(std::string& sForeignFile);

  /// Read and construct the dividends
  void ReadDividends(std::ifstream& sIn);
  void ReadDividends(std::string& sDividendFile);

  /// Read in the volatility
  void ReadVolatility(std::ifstream& sIn);
  void ReadVolatility(std::string& sVolatilityFile);
  void SetVolatility(const shared_ptr<Volatility>& volatility)
  {
    m_pVolatility = volatility;
  }

  /// Read in the hazard rate
  void ReadHazardRate(std::ifstream& sIn);
  void ReadHazardRate(std::string& sHazardRateFile);
  void SetHazardRate(const shared_ptr<HazardRate>& hazardRate)
  {
    m_pHazardRate = hazardRate;
  }

  /// Read and store the mesh parameters
  void ReadMeshParams(std::ifstream& sIn);
  void ReadMeshParams(std::string& sMeshParamFile);

  /// Read and store the numerical parameters
  void ReadNumParams(std::ifstream& sIn);
  void ReadNumParams(std::string& sNumParamFile);

  /// All testers must be able to read an input file
  virtual void ReadInputFile(std::string& sFilename);

  /// Read from multiple input files
  virtual void ReadFromMultipleFiles(std::string& sYieldCurveFile,
    std::string& sForeignCurveFile,
    std::string& sDividendFile,
    std::string& sVolatilityFile,
    std::string& sHazardRateFile,
    std::string& sMeshParamFile,
    std::string& sNumParamFile,
    std::string& sValuationDateFile,
    std::string& sContractFile);

  /// Read contract specific data. Pure virtual. 
  virtual void ReadContractParams(std::ifstream& sIn) = 0;
  virtual void ReadContractParams(std::string& sContractFile) = 0;

  /// Get the final price
  double GetPrice();

  /// Compute the final convergence rate for the price
  double GetConvergenceRate();

  /// Construct and run the appropriate pricer
  virtual void RunPricer() = 0;

  /// Output a surface
  void OutputSurface(const finance::SharedSurface& surface, 
    std::string sFileName, double dSLeft, double dSRight); 

  void OutputSurface(const shared_ptr<finance::SurfaceFlag>& surface, 
    std::string sFileName, double dSLeft, double dSRight); 

  /// Ouptut a single curve. 
  void OutputArray(const std::vector<double>& pdSpots, const std::vector<double>& pdValues,
    std::string sFileName);

  /** 
    Automatic testing and reporting functions
  */

  /// Run convergence testing
  void RunConvergenceTest(size_t nNbRuns);

  /// Given an array of values, compute and report convergence values
  void ReportConvergence(size_t nNbValues, double* pdValues, double* pdTimes,
                         std::string& strHeading);

  /// Run tests for the price
  bool RunPriceTests(double dExpected, double dTol = 1.e-4);

  /// Run convergence tests for the price
  bool RunPriceConvergenceTests(double dExpected, size_t nNbTests = 4);

  /// Set the output stream
  void SetOutputStream(std::ostream* outputStream)
  {
    m_outputStream = outputStream; 
  }

  /// Report if the two numbers are within the tolerance
  bool ReportPass(double dExpected, double dComputed, double dTol);


protected:

  /// Helper class to read yield/foreign curves
  void ReadCurve(std::ifstream&, int iWhich);

  /// The current spot price
  double m_dSpot;

  /// The pricing date
  Date m_ValuationDate;

  /// The yield curve
  shared_ptr<finance::YieldCurve> m_pYieldCurve;

  /// The foreign curve
  shared_ptr<finance::YieldCurve> m_pForeignCurve;

  /// The dividends
  shared_ptr<finance::Dividends> m_pDividends;

  /// The volatility
  shared_ptr<Volatility> m_pVolatility;

  /// The hazard rate
  shared_ptr<HazardRate> m_pHazardRate;

  /// Numerical parameters
  shared_ptr<numeric::NumParams> m_pNumParams;

  /// Mesh parameters
  shared_ptr<numeric::MeshParams> m_pMeshParams;

  /// Output stream used for testing
  std::ostream* m_outputStream;

  /// Output from a single pricing
  shared_ptr<finance::ModelOutput> m_pOutput;

  /// The number of runs for convergence testing
  size_t m_nNbRuns;

  /// These arrays are used for convergence testing
  Array<double> m_pdTimes;
  Array<double> m_pdPrices;
  Array<double> m_pdDeltas;
  Array<double> m_pdGammas;
  Array<double> m_pdVegas;
  Array<double> m_pdRhos;
  Array<double> m_pdThetas;

  /// Flags of what to compute
  bool m_bComputeDelta;
  bool m_bComputeGamma;
  bool m_bComputeVega;
  bool m_bComputeRho;
  bool m_bComputeTheta;
  bool m_bComputeArrays;
  bool m_bComputeSurfaces;


  NO_COPY_CLASS(BaseTester);
};


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_BASETESTER_H_
