/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/hedgeoutput.h
// Purpose:     hg hedging output class
// Created:     2005/09/23
// RCS-ID:      $Id: hedgeoutput.h,v 1.9 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/hedgeoutput.h

    Hedging output class for the homogeneous (HG) model.

    HERO is only computed if requested.
 */

#ifndef _ITO33_HG_HEDGEOUTPUT_H_
#define _ITO33_HG_HEDGEOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/hg/common.h"

namespace ito33
{

namespace hg
{

class ITO33_HG_DLLDECL HedgeRatioData;

/**
    Hedging output class for the homogeneous (HG) model.

    @iterator GetHedgeRatioData
    @nocreate
 */
class ITO33_HG_DLLDECL HedgeOutput
{
public:

  /**
      Default constructor.

      By default, HERO is not computed.
   */
  HedgeOutput()
    : m_bHasHERO(false)
  { }

  /// Dummy virtual dtor
  virtual ~HedgeOutput() { }


  /**
      @internal
      @brief The hedge ratio of the underlying.

      @noexport
   */
  void SetUnderlyingHedgeRatio(double dUnderlyingHedgeRatio) 
  { 
    m_dUnderlyingHedgeRatio = dUnderlyingHedgeRatio; 
  }

  /**
      @internal
      @brief The model output of the target contract.

      @noexport
   */
  void 
  SetTargetModelOutput(const shared_ptr<finance::ModelOutput>& pTargetOutput)
  { 
    m_pTargetModelOutput = pTargetOutput; 
  }

  /**
      @internal     
      @brief The HERO model output.

      @noexport
   */
  void SetHEROModelOutput(const shared_ptr<finance::ModelOutput>& pHEROOutput)
  { 
    m_bHasHERO = true;
    m_pHEROModelOutput = pHEROOutput; 
  }

  /**
      @internal
      @brief The hedging ratio data for all hedging contracts.

      Does not include the information for the target or underlying.

      @noexport
   */
  void SetHedgeRatioData(
    const std::vector< shared_ptr<HedgeRatioData> >& ppHedgeRatioData) 
  { 
    m_ppHedgeRatioData = ppHedgeRatioData; 
  }


  /**
      The hedging ratio data for all hedging contracts.

      Does not include the target or underlying information.

      @return the hedging ratio data for all hedge contracts
      @noexport COM
   */
  const std::vector< shared_ptr<HedgeRatioData> >& GetHedgeRatioData() const
  { 
    return m_ppHedgeRatioData; 
  }

  /**
      The hedge ratio of the underlying.

      @return the hedge ratio of the underlying
   */
  double GetUnderlyingHedgeRatio() const
  { 
    return m_dUnderlyingHedgeRatio; 
  }

  /**
      The model output of the target contract.

      @return the model output of the target contract
   */
  const shared_ptr<finance::ModelOutput>& GetTargetModelOutput() const
  { 
    return m_pTargetModelOutput; 
  }

  /**
      Check if the HERO was computed.

      @return true if HERO was computed, false otherwise
   */
  bool HasHERO() const 
  { 
    return m_bHasHERO; 
  }

  /**
      Check if the HERO surface was computed.

      @return true if the HERO surface was computed, false otherwise
   */
  bool HasHEROSurface() const 
  { 
    if ( m_bHasHERO == false )
      return false;

    return m_pHEROModelOutput->HasPriceSurface(); 
  }

  /**
      Check if the HERO was computed at the analysis date.

      @return true if HERO was computed at the analysis date, false otherwise
   */
  bool HasHEROAtAnalysisDate() const 
  { 
    if ( m_bHasHERO == false )
      return false;

    return m_pHEROModelOutput->HasPriceAtAnalysisDate(); 
  }

  /**
      The HERO value.

      @return the HERO value, if computed
   */
  double GetHERO() const
  { 
    if ( !HasHERO() )
      ThrowHERONotAvailable();

    return m_pHEROModelOutput->GetPrice(); 
  }

  /**
      The HERO surface.

      @return the HERO surface, if computed
   */
  shared_ptr<finance::SurfaceDouble> GetHEROSurface() const
  { 
    if ( !HasHERO() )
      ThrowHERONotAvailable();

    return m_pHEROModelOutput->GetPriceSurface(); 
  }

  /**
      The HERO analysis date spots.

      @return the HERO analysis date spots, if computed
   */
  const std::vector<double>& GetHEROSpotsAtAnalysisDate() const
  {
    if ( !HasHERO() )
      ThrowHERONotAvailable();

    return m_pHEROModelOutput->GetSpotsAtAnalysisDate(); 
  }

  /**
      The HERO analysis date values.

      @return the HERO analysis date values, if computed
   */
  const std::vector<double>& GetHEROValuesAtAnalysisDate() const
  { 
    if ( !HasHERO() )
      ThrowHERONotAvailable();

    return m_pHEROModelOutput->GetPricesAtAnalysisDate(); 
  }


protected:

  /// Throw error if hero is requested but not computed
  static void ThrowHERONotAvailable();

  /// Is the HERO computed/stored or not?
  bool m_bHasHERO;

  /// The pricing model output of the target contract
  shared_ptr<finance::ModelOutput> m_pTargetModelOutput;

  /// The hedging ratio for the underlying
  double m_dUnderlyingHedgeRatio;

  /// The hedging data for the hedge contracts
  std::vector< shared_ptr<HedgeRatioData> > m_ppHedgeRatioData;


private:

  /// The HERO output (accessed by wrapper functions)
  shared_ptr<finance::ModelOutput> m_pHEROModelOutput;

}; // class HedgeOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_HEDGEOUTPUT_H_
