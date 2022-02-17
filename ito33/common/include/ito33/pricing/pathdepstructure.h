/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/pathdepstructure.h
// Purpose:     path dependent structure class
// Author:      Yann and David
// Created:     2004/07/07
// RCS-ID:      $Id: pathdepstructure.h,v 1.15 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/pathdepstructure.h
    @brief path dependent structure class.
 */

#ifndef _ITO33_PRICING_PATHDEPSTRUCTURE_H_
#define _ITO33_PRICING_PATHDEPSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/instdata.h"

namespace ito33
{

namespace pricing
{
  class Params;
  class BackwardMeshManagerWithSpaceMesh;
  class StepperTriDiag;
  class PathDepEvent;
  

/// Path dependent structure class.
class PathDepStructure
{

  /**
      Keep this structure, or strip everything out into individual arrays?
   */
  struct Path
  {
    shared_ptr<Params> params;
    shared_ptr<BackwardMeshManagerWithSpaceMesh> meshes;
    shared_ptr<InstData> instdata;
    shared_ptr<StepperTriDiag> stepper;
  };

public:

  /// Empty destructor.
  virtual ~PathDepStructure() {}


  /** 
      Initializes the path to save.

      Called by the pricer after the last path dep event
   */
  virtual void InitPathToSave() = 0;

  /**
      Updates the numerical output.  
      
      Must be done in derived class since type information is needed.

      @param nIdx the path to update
   */
  virtual void UpdateOutput(size_t nIdx) = 0;

  /**
      Updates the numerical output at the end of grid. 

      Must be done in derived class since type information is needed.

      @param nIdx the path to update
   */
  virtual void UpdateOutputEndOfGrid(size_t nIdx) = 0;

  /**
      Finalizes the numerical output.

      Must be done in derived class since type information is needed.
   */
  virtual void Finalize() = 0;  

  /**
      Sets the path to save.

      @param nIdx the path to save
   */
  virtual void SetPathToSave(size_t nIdx);

  /**
      Gets the path to save.

      @return path to save
   */
  size_t GetPathToSave() const
  {
    return m_nPathToSave;
  }

  /**
      Gets the data for the specified path.

      Allows different models to return model appropriate data. Default 
      to returning the price array in the path instdata.

      @param nIdxPath The path to get the data from

      @return pointer to the price data
   */
  virtual double* GetPriceData(size_t nIdxPath)
  {
    return m_path[nIdxPath].instdata->m_pdPrices.Get();
  }

  /**
      Linearly interpolates in the x direction.

      @param nPath path number
      @param dS value to interpolate at

      @return interpolated value
   */
  double LinearInterpolate(size_t nPath, double dS);

  /**
      Quadratically interpolates in the x direction.

      @param nPath path number
      @param dS value to interpolate at

      @return quadratically interpolated value
   */
  double QuadraticInterpolate(size_t nPath, double dS);

  /**
      Similarity inerpolation assuming two paths.

      The 2nd path is used for extrapolation.

      @param nIdPathLow path number of lower value
      @param nIdPathHigh path number of higher value
      @param dY interpolated value in the y direction
      @param dS value of the spot
      @param pdGridY grid in the y direction
      @param nHomogeneousDegree degree of the homogeneity updating rule
      @param bIsLinearInterpolation indicate if linear or quadratic
                 interpolation is used, by default set to true
   */
  double SimilarityInterpolation(size_t nIdPathLow, 
    size_t nIdPathHigh, double dY, double dS, 
    const std::vector<double>& pdGridY, size_t nHomogeneousDegree,
    bool bIsLinearInterpolation = true);
   
  /**
      Similarity interpolation.

      @param nIdPath path number
      @param dY interpolated value in the y direction
      @param dS value of the spot
      @param dYStar value of the path for similarity reduction
      @param nHomogeneousDegree degree of the homogeneity updating rule
      @param bIsLinearInterpolation indicate if linear or quadratic
                 interpolation is used, by default set to true
   */
  double SimilarityInterpolation(size_t nIdPath, double dY, 
    double dS, double dYStar, size_t nHomogeneousDegree,
    bool bIsLinearInterpolation = true);
  

  /// Vector of all the underlying classes for the 1D solves
  std::vector<Path> m_path; 

  /// Data structure for the grids (excluding the 1D 'S' or 'X' grids)
  std::vector< std::vector< std::vector<double> > > m_pppdGrids; 

  /// State which paths are currently active
  std::vector<bool> m_pbIsActive;

protected:

  /**
      Ctor.

      The base class does not construct any objects, since it does not know 
      what structure is required.
   */
  PathDepStructure() { }

  /// Path to save for output
  size_t m_nPathToSave;
  
  /// Copy constructor is not allowed
  NO_COPY_CLASS(PathDepStructure);
};


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_PATHDEPSTRUCTURE_H_
