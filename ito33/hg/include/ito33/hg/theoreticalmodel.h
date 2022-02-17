/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/theoreticalmodel.h
// Purpose:     hg dynamic class
// Created:     2005/01/13
// RCS-ID:      $Id: theoreticalmodel.h,v 1.45 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/theoreticalmodel.h
    @brief HG theoretical model class
 */

#ifndef _ITO33_HG_THEORETICALMODEL_H_
#define _ITO33_HG_THEORETICALMODEL_H_

#include "ito33/vector.h"
#include "ito33/typeinfo.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/hg/underlyingprocess.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Option;
  class ITO33_DLLDECL OneTouch;
  class ITO33_DLLDECL FXOneTouch;
  class ITO33_DLLDECL CDSLike;
  class ITO33_DLLDECL CDS;
  class ITO33_DLLDECL ReferenceCDS;
  class ITO33_DLLDECL EDS;
  class ITO33_DLLDECL ParBond;
  class ITO33_DLLDECL Bond;
  class ITO33_DLLDECL ConvertibleBond;
  class ITO33_DLLDECL CBOption;
  class ITO33_DLLDECL Reset;
  class ITO33_DLLDECL PERCSLike;
  class ITO33_DLLDECL PEPSLike;
  class ITO33_DLLDECL GeneralizedPEPSLike;
  class ITO33_DLLDECL AsianOption;
  class ITO33_DLLDECL LogContract;
  class ITO33_DLLDECL VarianceSwapLike;
  class ITO33_DLLDECL VarianceSwap;
  class ITO33_DLLDECL GammaVarianceSwap;
  class ITO33_DLLDECL ConditionalVarianceSwap;
  class ITO33_DLLDECL OptionVarianceSwap;
  class ITO33_DLLDECL VarianceSwaption;

  class ITO33_DLLDECL Derivatives;
}

namespace pricing
{
  class VarianceSwap;
}

namespace hg
{

class ITO33_HG_DLLDECL HedgeOutput;
class ITO33_HG_DLLDECL HEROFlags;

class TransitionProbabilityOutput;

/// HG theoretical model implementation class.
class ITO33_HG_DLLDECL TheoreticalModel : public finance::TheoreticalModel
{
public:

  /// The type of RegisterPriceFunction() last parameter.
  typedef shared_ptr<finance::ModelOutput>
    (TheoreticalModel::*ModelPriceFunction)(const finance::Derivative&) const;

  /**
     ctor using a pre-defined underlying process. 
     
     @param pUnderlyingProcess The underlying process associated to this model
   */
  TheoreticalModel(const shared_ptr<UnderlyingProcess>& pUnderlyingProcess);

  /**
     Gets the risk neutral underlying process associated to this model.

     @return The underlying process associated to this model.
   */
  const shared_ptr<UnderlyingProcess>& GetUnderlyingProcess() const
  {
    return m_pUnderlyingProcess;
  }

  /** 
     The sharpe ratio of the real underlying process.

     @return The sharpe ratio of the real underlying process.
   */
  double GetSharpeRatio() const { return m_dSharpeRatio; }

  /** 
     Sets the risk neutral underlying process specifying the model parameters.

     @param pUnderlyingProcess a shared pointer to an underlying process for HG
   */
  void SetUnderlyingProcess
       (const shared_ptr<UnderlyingProcess>& pUnderlyingProcess)
  {
    m_pUnderlyingProcess = pUnderlyingProcess;
  }

  /**
     @internal
     @brief The underlying process used for mesh construction.

     This function should only be called for internally cloned TMs.

     @param pUnderlyingProcessForMesh a shared pointer to an underlying process
    
     @noexport
   */
  void SetUnderlyingProcessForMesh
    (const shared_ptr<UnderlyingProcess>& pUnderlyingProcessForMesh)
  {
    ASSERT_MSG( pUnderlyingProcessForMesh, 
                "Invalid process for mesh in TheoreticalModel class" );
    
    m_bHasExternalProcessForMesh = true;

    m_pUnderlyingProcessForMesh = pUnderlyingProcessForMesh;
  }

  /**
     @internal
     @brief The underlying process used for mesh construction.

     If the process for mesh has not been set, return the usual process.

     @return A shared pointer to the underlying process for mesh construction
    
     @noexport
   */
  const shared_ptr<UnderlyingProcess>& GetUnderlyingProcessForMesh() const
  {
    // If the mesh process was not defined, the normal process is used.
    // Always return a valid process due to the forward pricer requirements.
    if (m_bHasExternalProcessForMesh)
      return m_pUnderlyingProcessForMesh;
    else
      return m_pUnderlyingProcess;
  }

  /**
     The sharpe ratio of the real underlying process.

     @param dSharpeRatio The sharpe ratio of the real underlying process.
   */
  void SetSharpeRatio(double dSharpeRatio);

  /**
     Computes the real(historic) underlying process corresponding to the
     risk neutral underlying process.

     @return The real underlying process
   */
  shared_ptr<UnderlyingProcess> GetRealUnderlyingProcess() const
  {
    return m_pUnderlyingProcess->ComputeUnderlyingProcess(m_dSharpeRatio);
  }

  /**
     Computes the price and other characteristics of the given derivative.

     Throws if it is impossible to price the given kind of instruments with
     the hg model.

     Note that, to access the specifique data of a certain derivative, the
     user has to cast the output pointer to shared_ptr<Real_Output_Type>. 

     @param derivative the instrument to price
     @return ModelOutput shared pointer to the output object.
   */
  shared_ptr<finance::ModelOutput>
  Compute(const finance::Derivative& derivative) const;

#ifndef __CPP2ANY__

  /**
      Computes the derivative with the specified price function of the model.
     
      @param derivative the instrument to price
      @param pf function to be used for pricing of the given instrument
      @return ModelOutput shared pointer to the output object.

      @todo It is possible to make this function a template so the
            non registered method can be called without any explicit cast.

      @noexport
   */    
  shared_ptr<finance::ModelOutput>
  Compute(const finance::Derivative& derivative, ModelPriceFunction pf) const;

#endif

  /**
     @internal
     @brief Computes the (homogeneous) transition probability from time 0 to
            a given time.
             
     @param dPeriod time

     @return shared pointer to transition probability output class.

     @noexport
   */
  shared_ptr<TransitionProbabilityOutput>
  ComputeTransitionProbability(double dPeriod);

  /**
     Computes the hedge ratios by using a list of instruments against
     a target instrument.
  
     @param target The target instrument to be hedged
     @param hedgeInstruments The list of instruments used to hedge the target

     @return shared pointer to hedge output class, containing the hedge ratios,
             pricing outputs.
   */
  shared_ptr<HedgeOutput>
  Hedge(const finance::Derivative& target, 
        const finance::Derivatives& hedgeInstruments);

  /**
     Computes the hedge ratios and the hero by using a list of instruments
     against a target instrument.
  
     @param target The target instrument to be hedged
     @param hedgeInstruments The list of instruments used to hedge the target
     @param flags flags for additional output of HERO

     @return shared pointer to hedge output class, containing the hedge ratios,
             pricing outputs, and hero.
   */
  shared_ptr<HedgeOutput>
  ComputeHERO(const finance::Derivative& target, 
              const finance::Derivatives& hedgeInstruments,
              shared_ptr<HEROFlags> flags = shared_ptr<HEROFlags>() );

#ifndef __CPP2ANY__
  /**
     Register the given method for pricing of the given kind of derivatives.

     It is ok to call this method multiple times.

     Instead of calling it directly it is preferable to use
     RegisterPriceFunction helper template class below.

     @param tiDerivative the typeinfo of the instrument class
     @param pf the price function to be used for this model and instrument
   */
  static void DoRegisterPriceFunction(const TypeInfo& tiDerivative,
                                      ModelPriceFunction pf);

  /**
     @name pricing functions
  */
  //@{

  /**
     @internal
     @brief Prices an option contract.

     @param option Option to be priced.
    
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceOption(const finance::Option& option);

  /**
     @internal
     @brief Prices a one touch contract.

     @param oneTouch one touch to be priced.
    
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceOneTouch(const finance::OneTouch& oneTouch);

  /**
     @internal
     @brief Prices a FX one touch contract.

     @param oneTouch one touch to be priced.
    
     @noexport
   */
  shared_ptr<finance::ModelOutput>
  PriceFXOneTouch(const finance::FXOneTouch& oneTouch);

  /**
     @internal
     @brief Prices a CDS contract.

     @param cds CDS to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceCDS(const finance::CDS& cds);

  /**
     @internal
     @brief Prices a reference CDS contract.

     @param refCDS reference CDS to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceReferenceCDS(const finance::ReferenceCDS& refCDS);

  /**
     @internal
     @brief Prices a EDS contract.

     @param eds EDS to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceEDS(const finance::EDS& eds);

  /**
     @internal
     @brief Prices a ParBond contract.

     @param parBond ParBond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceParBond(const finance::ParBond& parBond);

  /**
     @internal
     @brief Prices a bond contract.

     @param bond bond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceBond(const finance::Bond& bond);

  /**
     @internal
     @brief  Prices a convertible bond contract.

     @param cb convertible bond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceConvertibleBond(const finance::ConvertibleBond& cb);

  /**
     @internal
     @brief Prices a convertible bond option contract.

     @param cboption convertible bond option to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput>
  PriceCBOption(const finance::CBOption& cboption);

  /**
     @internal
     @brief Prices a PERCS-like contract.

     @param percsLike PERCS-like to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput>   
  PricePERCSLike(const finance::PERCSLike& percsLike);

  /**
     @internal
     @brief Prices a PEPS-like contract.

     @param pepsLike PEPS-like to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput>
  PricePEPSLike(const finance::PEPSLike& pepsLike);

  /**
     @internal
     @brief Prices a Generalized PEPS-like contract.

     @param pepsLike Generalized PEPS-like to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceGeneralizedPEPSLike(const finance::GeneralizedPEPSLike& pepsLike);

  /**
     @internal
     @brief Prices a resettable convertible bond contract.

     @param reset resettable convertible bond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceReset(const finance::Reset& reset);

  /**
     @internal
     @brief Prices an Asian option contract.

     @param asianOption Asian option contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceAsianOption(const finance::AsianOption& asianOption);

  /**
     @internal
     @brief Prices a log contract.

     @param logContract Log contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceLogContract(const finance::LogContract& logContract);

  /**
     @internal
     @brief Prices a standard variance swap contract.

     @param varianceSwap standard variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceVarianceSwap(const finance::VarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices a standard variance swap contract by using the log contract

     Note: this is an unregistered function and must be called in the following
     way:

     Compute( varianceSwap, 
              reinterpret_cast<ModelPriceFunction>(PriceVarianceSwapByLog) ).

     Same mechanism can be applied to other (non) registered functions.

     @param varianceSwap standard variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceVarianceSwapByLog(const finance::VarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices a gamma variance swap contract.

     @param varianceSwap gamma variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceGammaVarianceSwap(const finance::GammaVarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices a conditional variance swap contract.

     @param varianceSwap conditional variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceConditionalVarianceSwap(
    const finance::ConditionalVarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices an option variance swap contract.

     @param varianceSwap option variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceOptionVarianceSwap(const finance::OptionVarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices a variance swaption contract.

     @param varianceSwaption variance swaption contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceVarianceSwaption(const finance::VarianceSwaption& varianceSwaption);

  //@} // @name pricing functions 

  TheoreticalModel* Clone() const;

  void Dump(ito33::XML::Tag& tagParent) const;

protected:

  /// The homogeneous underlying process
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcess;

  /// The homogeneous underlying process for mesh construction
  mutable shared_ptr<UnderlyingProcess> m_pUnderlyingProcessForMesh;

  /// Whether or not the process for mesh was set externally
  bool m_bHasExternalProcessForMesh;

  /// The sharpe ratio of the real underlying process.
  double m_dSharpeRatio;

private:

  /**
      Computes the price and the greeks for the given derivative.

      Note that, to access the specific data of a certain derivative, the
      user has to cast the output pointer to shared_ptr<Real_Output_Type>. 

      @param derivative the instrument to price
      @param pf the appropriate function to price with
      @return shared pointer to the model output object.
   */
  shared_ptr<finance::ModelOutput>
  ComputeAll(const finance::Derivative& derivative, ModelPriceFunction pf) const;

  shared_ptr<finance::ModelOutput> 
  DoCompute(const finance::Derivative& derivative) const
  {
    return Compute(derivative);
  }

  /**
     Implementation of Hedge() and ComputeHERO().
  
     @param target The target instrument to be hedged
     @param hedgeInstruments The list of instruments used to hedge the target
     @param bComputeHERO whether hero is requested
     @param flags flags for additional output of HERO

     @return shared pointer to hedge output class, containing the hedge ratios,
             pricing outputs, and hero.
   */
  shared_ptr<HedgeOutput>
  DoComputeHERO(const finance::Derivative& target, 
                const finance::Derivatives& hedgeInstruments,
                bool bcomputeHERO,
                shared_ptr<HEROFlags> flags);

  /**
     Computes the price and the greeks for a CDS-like contract.

     Called by the specific (derived) CDS-like price functions, such as
     CDS and ReferenceCDS.

     @param cdsLike CDS-like contract to be priced.

     @return ModelOutput shared pointer to the output object.
   */
  shared_ptr<finance::ModelOutput> 
  DoPriceCDSLike(const finance::CDSLike& cdsLike);

  /**
     Computes the price and the greeks for a VarianceSwap-like contract.

     Called by the specific (derived) VarianceSwap-like price functions, such
     as for GammaVarianceSwap or ConditionaVarianceSwap.

     @param varianceSwapLike VarianceSwap-like contract to be priced
     @param varSwap contains specific information from the derived variance
                    swap classes

     @return ModelOutput shared pointer to the output object.
   */
  shared_ptr<finance::ModelOutput> 
  DoPriceVarianceSwapLike(const finance::VarianceSwapLike& varianceSwapLike,
                          pricing::VarianceSwap& varSwap);

  // find the function to be used for pricing of the given derivative
  // (common part of CanPrice() and ComputePrice())
  ModelPriceFunction 
  FindPriceFunction(const finance::Derivative& derivative) const;
  
}; // class TheoreticalModel

/**
   This is just a wrapper around TheoreticalModel::DoRegisterPriceFunction().

   It is a convenient wrapper nevertheless as you don't have to call typeid()
   manually when using it and you must also specify the function with correct
   signature unlike when using DoRegisterPriceFunction directly as this class
   will do the necessary casts.

   Notice that it is impossible to make DoRegisterPriceFunction() a template
   function itself because MSVC6 can't compile it then.

   Template parameters are:
     - D the instrument class deriving from Derivative
 */
template <class D>
class RegisterPriceFunction
{
public:
  /**
     Ctor registers the specified pricing function for pricing our instrument
     class with our model.

     A static RegisterPriceFunction is usually created globally to ensure that
     the price function is registered before it might be used. But, of course,
     it is also possible to do it later.

     @param pf the price function to be used for this model and instrument
   */
  RegisterPriceFunction( 
    shared_ptr<finance::ModelOutput> (TheoreticalModel::*pf)(const D&) )
  {
    // ensure that D can be converted to Derivative, otherwise reinterpret_cast
    // below is not going to work
    static_cast<finance::Derivative *>(static_cast<D *>(0));

    TheoreticalModel::ModelPriceFunction 
      mpf = reinterpret_cast<TheoreticalModel::ModelPriceFunction>(pf);

    TheoreticalModel::DoRegisterPriceFunction
                      (
                        typeid(D),
                        static_cast<TheoreticalModel::ModelPriceFunction>(mpf)
                      );
  }
#endif // __CPP2ANY__
};


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_THEORETICALMODEL_H_
