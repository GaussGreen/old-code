/////////////////////////////////////////////////////////////////////////////
// Name:      hg/backwardnumoutput.h
// Purpose:   implementation of HG base backward numoutput class 
// Created:   2005/01/13 
// RCS-ID:    $Id: backwardnumoutput.h,v 1.25 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/backwardnumoutput.h
    @brief implementation of HG base backward numoutput class 

    Note that time only numoutput classes should not be derived from
    this class.

    The flags inside BackwardNumOutput come from three parts: general flags as
    ComputeSurface and analysis date come directly from pricer; flag for greeks
    computed by PDE come from InstData; flags for greeks computed by 
    perturbation come from the pricer.
 */

#ifndef _HG_BACKWARDNUMOUTPUT_H_
#define _HG_BACKWARDNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/surfacedouble.h"

#include "hg/numoutput.h"
#include "hg/sensitivitybyadjointdata.h"

namespace ito33
{

namespace numeric
{
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace pricing
{
  class Params;
}

namespace finance
{
  class ModelOutput;
}

namespace hg
{

class BackwardInstData;
class ModelOutput;
class MultiOutput;

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class BackwardNumOutput : public NumOutput
{
public:

  /**
     Constructor by params
 
     @param params reference to pricing params
   */
  BackwardNumOutput(pricing::Params& params) 
                  : NumOutput(), 
                    m_params(params),
                    m_bHasConstraintFlags(false),
                    m_bSensitivityOnObjectif(false)
  {
  }

  /// virtual dtor
  virtual ~BackwardNumOutput() { }

  /**
     Initialize the class variables, the default implementations considers
     that the space mesh is fixed.

     instdata::SetInitialState must have been called
  */
  void Init(BackwardInstData& instdata);

  /**
     If required, store the pricing data at the current timestep
   */
  virtual void UpdateMe(BackwardInstData& instdata, double dTime);
 
  /**
     If required, store the pricing data at the current timestep which
     is end of grid
   */
  void UpdateMeAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /**
     Get the spots at the anaysis date

     @return The spots at the analysis date
   */
  const finance::Values& GetSpotsAtAnalysisDate() const 
  { 
    return m_pdAnalysisSpots; 
  }

  /**
     Get the price values at the anaysis date

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
  shared_ptr<numeric::SurfaceGeneral> GetPriceSurface(size_t nRegime = 0)
  {
    ASSERT_MSG(nRegime < m_nNbRegimes, 
               "Invalid regime in call to GetPriceSurface(nRegime)");

    if ( m_ppPriceSurfaces.empty() )
      return shared_ptr<numeric::SurfaceGeneral>();

    return m_ppPriceSurfaces[nRegime]; 
  }

  /**
     Saves the Final results such as prices, fugits for other use. For example,
     the initial condition of another problem may depend on them.

     Normally we call this function when the flag bFinalSave is set as true
    
     @param instdata instdata at the final time where we get general information
   */
  virtual void Finalize(BackwardInstData& instdata);

  /**
      Sets the recovery value at valuation date.
   
      @param dValueAfterDefault the recovery value at valuation date
   */
  void SetValueAfterDefault(double dValueAfterDefault)
  {
    m_dValueAfterDefault = dValueAfterDefault;
  }

  /**
     Gets the recovery values.

     @return The recovery values
   */
  const std::vector<double>& GetRecoveryValues()
  {
    return m_pdRecoveryValues;
  }

  /**
     Gets the pricing times.

     @return The pricing times
   */
  const std::vector<double>& GetTimes()
  {
    return m_pdTimes;
  }

  /**
     Subtract the specified numoutput data from this data.

     Used when a payoff is A - B, and the A and B parts are priced separately
     and need to be combined.  An example is conditional variance swaps.

   */
  void SubractNumOutput(const BackwardNumOutput& otherNumOutput);

  virtual shared_ptr<ModelOutput> GetModelOutput();

  virtual shared_ptr<MultiOutput> GetMultiOutput();

protected:
  
  /**
     Calculates the scalar data at the final time

     @param instdata instdata at the final time where we get general information
   */
  virtual void CalculateFinalScalarResult(BackwardInstData& instdata);

  virtual void SaveAnalysisData(BackwardInstData& instdata);
  
  /**
     Saves the surface, default implementation considers that the space mesh
     is fixed, do doesn't change with the time.
     
     @param instdata The instdata holds the data to be saved.
   */
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  virtual void SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /**
     Called by UpdateMe function and do other work than Analysis Data Saving,
     Surface Saving.

     @param instdata instdata
   */
  virtual void DoSpecialWorkWith(BackwardInstData& instdata, double dTime);

  void SaveSurfaceDataFrom
       ( BackwardInstData& instdata, bool bEndOfGrid = false );
  
  /**
     This function is to be called by the Init() function of specific NumOutput classes
     for the problems using fix space mesh.

     @param instdata backward inst data
     @param meshes fix backward mesh manager
     @param params base params to get Analysis time
     @param bHasConstraintsFlags whether allocate Flag Surface
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
     @param nNbRegimes the number of regimes to copy
   */
  template <typename T>
  void CopyData(const BackwardInstData& instdata, 
                const T* pdSource, T* pdTarget, size_t nNbRegimes);

  /**
     Determine the number of extra points to add at boundaries.

     Extra points are only added at Dirichlet nodes.

     @param instdata backward instdata reference
     @param nNbLeft (output) number of points to add at left boundary
     @param nNbRight (output) number of points to add at right boundary
   */
  void GetNbExtraPoints(const BackwardInstData& instdata, 
                        size_t& nNbLeft, size_t& nNbRight);

  void SensitivityByAdjoint(BackwardInstData& instdata);

  void SolveTransposeSystem
       (BackwardInstData& instdata, const double* pdRHS, double* pdValues);

  /**
     Fills the given output.

     @param output the output to be filled by this numoutput
   */
  void GetOutput(finance::ModelOutput& output);

protected:

  void SensitivityByAdjointFD(BackwardInstData& instdata);

  void SolveTransposeSystemFD
       (BackwardInstData& instdata, const double* pdRHS, double* pdValues);

  void SensitivityByAdjointFE(BackwardInstData& instdata);

  void SolveTransposeSystemFE
       (BackwardInstData& instdata, const double* pdRHS, double* pdValues);

  /// The size of the values
  size_t m_nNbX;

  /// The size of the space mesh
  size_t m_nNbS;

  /// The number of regimes
  size_t m_nNbRegimes;

  /// pointer to space mesh
  const double *m_pdS;

  /// Price at default regime
  double m_dValueAfterDefault;

  /// fugit at valuation date
  double m_dFugit;

  std::vector<SensitivityByAdjointData> m_pAdjointDatas;

  /// The pricing domain. Only set if the surface is computed
  shared_ptr<numeric::Domain> m_pDomain;

  /// The price surfaces (one per regime). Only set if surfaces are computed.
  std::vector< shared_ptr<numeric::SurfaceGeneral> > m_ppPriceSurfaces;

  // The delta and gamma surface are calculated from the price surface.
  // They are not saved.

  /// The Fugit surface. Only set if surfaces and Fugit are computed
  shared_ptr<numeric::SurfaceGeneral> m_pFugitSurface;

  /// The American constraint surfaces.  Only set for American options
  shared_ptr<numeric::SurfaceFlag> m_pFlagSurface;

  /// Boolean to indicate if there has constraint flags
  bool m_bHasConstraintFlags;

  /// Analysis time
  double m_dAnalysisTime;

  /// The prices at the analysis date
  finance::Values m_pdAnalysisSpots;

  /// The prices at the analysis date
  finance::Values m_pdAnalysisValues;

  /// The fugit at the analysis date
  finance::Values m_pdAnalysisFugits;

  /// The theta at the analysis date
  finance::Values m_pdAnalysisThetas;

  /// The constraint values at the analysis date
  finance::Flags m_pAnalysisFlags;

  /// The recovery values for each time
  std::vector<double> m_pdRecoveryValues;

  /// The (non-duplicated due to events) pricing times
  std::vector<double> m_pdTimes;

  pricing::Params& m_params;

  bool m_bSensitivityOnObjectif;

  /// The prices at each of the observation points
  std::vector<double> m_pdPricesAtObservations;

  double m_dObjectif;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(BackwardNumOutput);

}; // class BackwardNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_BACKWARDNUMOUTPUT_H_
