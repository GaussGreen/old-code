/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/backwardnumoutput.h
// Purpose:   implementation of base backward numoutput class 
// Author:    ZHANG Yunzhi
// RCS-ID:    $Id: backwardnumoutput.h,v 1.25 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/backwardnumoutput.h
    @brief implementation of base backward numoutput class 

    Note that time only numoutput classes should not be derived from
    this class.

    The flags inside BackwardNumOutput come from two parts: general flags as
    ComputeSurface and analysis date come directly from pricer and flag for 
    greeks computed by PDE come from InstData.
 */

#ifndef _IHG_BACKWARDNUMOUTPUT_H_
#define _IHG_BACKWARDNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/surfacedouble.h"
#include "ihg/numoutput.h"


namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace pricing
{
  class Params;
}

namespace ihg
{

class BackwardInstData;

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class BackwardNumOutput : public NumOutput
{
public:
  
  /// default ctor
  BackwardNumOutput() : NumOutput(), m_bHasConstraintFlags(false)
  { 
    m_nNbS = 0;

    m_dAnalysisTime = - 100;
  }

  /// virtual dtor
  virtual ~BackwardNumOutput() { }

  /**
     If required, store the pricing data at the current timestep
   */
  virtual void UpdateMe(BackwardInstData& instdata, double dTime);
 
  /**
     If required, store the pricing data at the current timestep which
     is end of grid
   */
  virtual void UpdateMeAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /**
     Get the price values at the analysis date

     @return The price values at the analysis date
   */
  const finance::Values& GetPricesAtAnalysisDate() const 
  { 
    return m_pdAnalysisValues; 
  }

  /**
     Get the price surface

     @return The price surface
   */
  shared_ptr<numeric::SurfaceGeneral> GetPriceSurface()
  {
    return m_pPriceSurface; 
  }

  /**
     Get the fugit surface

     @return The fugit surface
   */
  shared_ptr<numeric::SurfaceGeneral> GetFugitSurface()
  {
    return m_pFugitSurface; 
  }

  /**
     Get the flag surface

     @return The flag surface
   */
  shared_ptr<numeric::SurfaceFlag> GetFlagSurface()
  {
    return m_pFlagSurface; 
  }

  /**
     Saves the Final results such as prices, vegas for other use. For example,
     the initial condition of another problem may depend on them.

     Normally we call this function when the flag bFinalSave is set as true
    
     @param instdata instdata at the final time where we get general information
   */
  virtual void Finalize(BackwardInstData& instdata);


protected:
  
  /**
     Calculates the scalar data at the final time

     @param instdata instdata at the final time where we get general information
   */
  virtual void CalculateFinalScalarResult(BackwardInstData& instdata);

  virtual void SaveAnalysisData(BackwardInstData& instdata);
  
  virtual void SaveSurface(BackwardInstData& instdata, double dTime) = 0;

  virtual void SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /**
     Called by UpdateMe function and do other work than Analysis Data Saving,
     Surface Saving.

     @param instdata instdata
   */
  virtual void DoSpecialWorkWith(BackwardInstData& instdata, double dTime);

  void GetOutput(finance::ModelOutput* pOutput);

  void SaveSurfaceWithFixMesh(BackwardInstData& instdata,
                              numeric::DomainFixedSpaceMesh& domain,
                              double dTime);

  virtual void SaveSurfaceDataFrom
               (
                BackwardInstData& instdata,
                bool bEndOfGrid = false
               );
  
  /**
     This function is to be called by the Init() function of specific NumOutput classes
     for the problems using fix space mesh.

     @param instdata backward inst data
     @param params base params to get Analysis time
     @param bHasConstraintFlags whether allocate Flag Surface
   */
  void InitWithFixMesh
        (
          BackwardInstData& instdata,
          const pricing::Params& params,
          bool bHasConstraintFlags
        );

  /**
     Copy the space mesh from instdata, adding extra boundary points if needed.

     @param instdata backward instdata reference

     @return A copy of the current space mesh
   */
  std::vector<double> CopySpaceMesh(const BackwardInstData& instdata);

  /**
     Copy data from the source array to the target array, adding extra boundary
     points if needed.

     Templated so that it works for double arrays (price, most Greeks)
     and int arrays (constraint flags).

     @param instdata backward instdata reference
     @param pdSource where to copy the data from
     @param pdTarget the vector to copy the data to
   */
  template <typename T>
  void CopyData(const BackwardInstData& instdata, 
                const T* pdSource, T* pdTarget);

  /**
     Determine the number of extra points to add at boundaries.

     Extra points are only added at Dirichlet nodes.

     @param instdata backward instdata reference
     @param nNbLeft (output) number of points to add at left boundary
     @param nNbRight (output) number of points to add at right boundary
   */
  void GetNbExtraPoints(const BackwardInstData& instdata, 
                        size_t& nNbLeft, size_t& nNbRight);

protected:

  /// The size of the space mesh
  size_t m_nNbS;

  /// pointer to space mesh
  const double *m_pdS;

  /// vega on spot at final time
  double m_dVega;

  /// fugit at valuation date
  double m_dFugit;

  /// The pricing domain. Only set if the surface is computed
  shared_ptr<numeric::Domain> m_pDomain;

  /// The price surface. Only set if surfaces are computed
  shared_ptr<numeric::SurfaceGeneral> m_pPriceSurface;

  // The delta and gamma surface are calculated by price gamma.
  // they are not saved.

  /// The vega surface. Only set if surfaces and vega are computed
  shared_ptr<numeric::SurfaceGeneral> m_pVegaSurface;

  /// The Fugit surface. Only set if surfaces and Fugit are computed
  shared_ptr<numeric::SurfaceGeneral> m_pFugitSurface;

  /// The American constraint surface.  Only set for American options
  shared_ptr<numeric::SurfaceFlag> m_pFlagSurface;

  /// Boolean to indicate if there has constraint flags
  bool m_bHasConstraintFlags;

  /// Analysis time
  double m_dAnalysisTime;

  /// The prices at the analysis date
  finance::Values m_pdAnalysisSpots;

  /// The prices at the analysis date
  finance::Values m_pdAnalysisValues;

  /// The vega at the analysis date
  finance::Values m_pdAnalysisVegas;

  /// The fugit at the analysis date
  finance::Values m_pdAnalysisFugits;

  /// The theta at the analysis date
  finance::Values m_pdAnalysisThetas;

  /// The constraint values at the analysis date
  finance::Flags m_pAnalysisFlags;

  /// The final vega array (presumably at valuation date)
  std::vector<double> m_pdFinalVegas;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(BackwardNumOutput);

}; // class BackwardNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_BACKWARDNUMOUTPUT_H_
