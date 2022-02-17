/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/theoreticalmodel.h
// Purpose:     ihg dynamic class
// Author:      ZHANG Yunzhi
// Created:     Feb 24, 2004
// RCS-ID:      $Id: theoreticalmodel.h,v 1.97 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/theoreticalmodel.h
    @brief ihg theoretical model class
 */

#ifndef _ITO33_IHG_THEORETICALMODEL_H_
#define _ITO33_IHG_THEORETICALMODEL_H_

#include "ito33/vector.h"
#include "ito33/typeinfo.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/common.h"
#include "ito33/ihg/underlyingprocess.h"

namespace ito33
{

class ITO33_DLLDECL Date;

namespace finance
{
  class ITO33_DLLDECL AsianOption;
  class ITO33_DLLDECL Option;
  class ITO33_DLLDECL CDSLike;
  class ITO33_DLLDECL CDS;
  class ITO33_DLLDECL ReferenceCDS;
  class ITO33_DLLDECL EDS;
  class ITO33_DLLDECL OneTouch;
  class ITO33_DLLDECL FXOneTouch;
  class ITO33_DLLDECL VarianceSwapLike;
  class ITO33_DLLDECL VarianceSwap;
  class ITO33_DLLDECL GammaVarianceSwap;
  class ITO33_DLLDECL OptionVarianceSwap;
  class ITO33_DLLDECL Bond;
  class ITO33_DLLDECL AttachedWarrantConvertibleBond;
  class ITO33_DLLDECL ConvertibleBond;
  class ITO33_DLLDECL CBOption;
  class ITO33_DLLDECL Reset;
  class ITO33_DLLDECL PERCSLike;
  class ITO33_DLLDECL PEPSLike;
  class ITO33_DLLDECL GeneralizedPEPSLike;
  class ITO33_DLLDECL ParBond;
  class ModelParametersConsumer;
}

namespace pricing
{
  class VarianceSwap;
}

namespace ihg
{

  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL HazardRate;
  class ITO33_IHG_DLLDECL PerfectHedgeRatios;

/**
    IHG theoretical model implementation class.
 */
class ITO33_IHG_DLLDECL TheoreticalModel : public finance::TheoreticalModel
{
public:

  /// The type of RegisterPriceFunction() last parameter.
  typedef shared_ptr<finance::ModelOutput>
          (ihg::TheoreticalModel::*PriceFunction)
          (
            const finance::Derivative&
          ) const;

  /// ctor
  TheoreticalModel();

  /**
     Override default constructor so that a tm can be constructed
     directly by passing underlying process instead of calling the set
     function.

     @param pUnderlyingProcess shared pointer to underlying process

     @noexport
  */
  TheoreticalModel(const shared_ptr<UnderlyingProcess>& pUnderlyingProcess);

  /**
     Override default constructor so that a tm can be constructed
     directly by passing vol and hr instead of calling the set
     functions.

     @param pVolatility shared pointer to volatility
     @param pHazardRate shared pointer to hazard rate

     @noexport
  */
  TheoreticalModel(const shared_ptr<Volatility>& pVolatility, 
    const shared_ptr<HazardRate>& pHazardRate);

  /**
      @name Initialization functions.
   */
  //@{

  /** 
     Sets the underlying process specifying the model parameters.

     @param pUnderlyingProcess a shared pointer to an underlying process for 
                               IHG.
   */
  void 
  SetUnderlyingProcess(const shared_ptr<UnderlyingProcess>& pUnderlyingProcess);

  /**
     @internal
     @brief The underlying process for the mesh.

     @param pUnderlyingProcessForMesh shared pointer to an underlying process
    
     @noexport
   */
  void SetUnderlyingProcessForMesh(const shared_ptr<UnderlyingProcess>& 
                                   pUnderlyingProcessForMesh)
  {
    ASSERT_MSG( pUnderlyingProcessForMesh, 
                "Invalid underlying process for mesh in TheoreticalModel "
                "class." );
    
    m_pUnderlyingProcessForMesh = pUnderlyingProcessForMesh;
  }

  /**
     The volatility of the model.

     @param pVolatility shared pointer to a Volatility
   */
  void SetVolatility(const shared_ptr<Volatility>& pVolatility);

  /**
     The hazard rate of the model.

     @param pHazardRate shared pointer to a HazardRate
   */
  void SetHazardRate(const shared_ptr<HazardRate>& pHazardRate);
 

  //@} // name Initialization functions
  
  // see base class
  virtual const shared_ptr<ihg::UnderlyingProcess>& 
  GetUnderlyingProcess() const;

  /**
     @internal
     @brief The underlying process for mesh of the model.

     @return shared pointer to the internal underlying process for mesh
    
     @noexport
   */
  const shared_ptr<UnderlyingProcess>& GetUnderlyingProcessForMesh() const;

  /**
     The volatility of the model.

     @return shared pointer to the internal volatility
   */
  const shared_ptr<Volatility>& GetVolatility() const;

  /**  
     @internal

     The volatility for mesh of the model.

     @return shared pointer to the internal volatility for mesh

     @noexport
   */
  const shared_ptr<Volatility>& GetVolatilityForMesh() const
  {
    return GetUnderlyingProcessForMesh()->GetVolatility();
  }

  /**
     Gets the hazard rate of the model.

     @return shared pointer to the internal hazard rate
   */
  const shared_ptr<HazardRate>& GetHazardRate() const;

  /**
      Computes the price and other characteristics of the given derivative.

      Throws if it is impossible to price the given kind of instruments with
      the ihg model.

      Note that, to access the specific data of a certain derivative, the
      user has to cast the output pointer to shared_ptr<Real_Output_Type>. 

      @param derivative the instrument to price
      @return ModelOutput shared pointer to the output object.
   */
  shared_ptr<finance::ModelOutput>
  Compute(const finance::Derivative& derivative) const;

  /**
     Computes the implied brownian volatility of the given derivative.

     Throws if it can't find the implied brownian vol.

     @param derivative the instrument to price
     @return the implied brownian volatility, with known hazard rate
   */
  double
  ComputeImpliedBrownianVol(const finance::Derivative& derivative) const;

  /**
     Computes the cumulative probability of non default of the given model.

     @param sessionData the session data
     @param maturityDates the array of maturities at which the probability will 
                          be computed.
     @return the array of the cumulative default probability
     
   */
  std::vector<double> 
  ComputeCumulativeDefaultProbability(const finance::SessionData& sessionData,
                                      const std::vector<Date>& maturityDates);

  /**
     Computes for the derivative instrument to be hedged, the hedge ratios so 
     as to match the jump to default and the delta of the considered derivative 
     to be hedged.

     @param target Derivative to be hedged. 
     @param hedge Derivative used to hedge the target derivative. 
     @return A PerfectHedgeRatios object.
    
     @todo The local characteristic of the underlying must be modified after
           the modification of the output delta.(Now, delta has a dimension
           for the cross currency case but we will change this as soon as 
           possible).

     Note: If the derivative is an exchangeable bond, an exception is thrown.
   */
  PerfectHedgeRatios ComputePerfectHedgeRatios
                            (
                              const finance::Derivative& target,
                              const finance::Derivative& hedge
                            ) const;

#ifndef __CPP2ANY__
  /**
      @internal
      @brief Register the given method for pricing of the given kind of
             derivatives.

      It is ok to call this method multiple times.

      Instead of calling it directly it is preferable to use
      RegisterPriceFunction helper template class below.

      @param tiDerivative the typeinfo of the instrument class
      @param pf the price function to be used for this model and instrument
   */
  static void DoRegisterPriceFunction(const TypeInfo& tiDerivative,
                                      PriceFunction pf);

  /// @name pricing functions
  //@{

  /**
     @internal
     @brief Prices an asian option contract.

     @param asianOption asian option to be priced
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceAsianOption(const finance::AsianOption& asianOption);

  /**
     @internal
     @brief Prices an option contract.

     @param option Option to be priced.
    
     @noexport
   */
  shared_ptr<finance::ModelOutput> PriceOption(const finance::Option& option);

  /**
     @internal
     @brief Prices a CDS contract.

     @param cds CDS to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> PriceCDS(const finance::CDS& cds);

  /**
     @internal
     @brief Prices a ReferenceCDS contract.

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
  shared_ptr<finance::ModelOutput> PriceEDS(const finance::EDS& eds);

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
     @brief Prices a bond contract.

     @param bond bond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> PriceBond(const finance::Bond& bond);

  /**
     @internal
     @brief Prices attached warrant convertible bond contract.

     @param warrant an attached warrant convertible bond to be priced

     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceAttachedWarrantConvertibleBond
  ( const finance::AttachedWarrantConvertibleBond& warrant );

  /**
     @internal
     @brief Prices a convertible bond contract.

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
     @brief Prices a par bond.

     @param parbond par bond to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceParBond(const finance::ParBond& parbond);

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
     @brief Prices a gamma variance swap contract.

     @param varianceSwap gamma variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceGammaVarianceSwap(const finance::GammaVarianceSwap& varianceSwap);

  /**
     @internal
     @brief Prices an option variance swap contract.

     @param varianceSwap option variance swap contract to be priced.
     
     @noexport
   */
  shared_ptr<finance::ModelOutput> 
  PriceOptionVarianceSwap(const finance::OptionVarianceSwap& varianceSwap);

  //@} // @name pricing functions 

  virtual TheoreticalModel* Clone() const;

  // defined in copytheoreticalmodel.cpp
  virtual TheoreticalModel* DeepCopy() const;

  virtual void Dump(ito33::XML::Tag& tagParent) const;
  
  virtual void GetModelParameters
               (finance::ModelParametersConsumer& visitor) const;

protected:

  /// checks if the underlying process is valid.
  void CheckUnderlyingProcess() const;  

  /// checks if model parameters are valid.
  void CheckAll() const;

private:
  
  /// @name IHG specific parameters.
  //@{

  /// The underlying process of the ihg model 
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcess;

  /// The underlying process for mesh of the model 
  mutable shared_ptr<UnderlyingProcess> m_pUnderlyingProcessForMesh;

  //@}
  
private:
  
  /**
      Computes the price and the greeks for the given derivative.

      Note that, to access the specific data of a certain derivative, the
      user has to cast the output pointer to shared_ptr<Real_Output_Type>. 

      @param derivative the instrument to price
      @param PriceFunction pf the appropriate function to price with
      @return ModelOutput shared pointer to the output object.
   */
  shared_ptr<finance::ModelOutput>
  ComputeAll(const finance::Derivative& derivative, PriceFunction pf) const;

  virtual shared_ptr<finance::ModelOutput> 
    DoCompute(const finance::Derivative& derivative) const
  {
    return Compute(derivative);
  }

  // find the function to be used for pricing of the given derivative
  // (common part of CanPrice() and ComputePrice())
  PriceFunction FindPriceFunction(const finance::Derivative& derivative) const;
  
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
  RegisterPriceFunction
                (
                shared_ptr<finance::ModelOutput> (TheoreticalModel::*pf)(const D&)
                )
  {
    // ensure that D can be converted to Derivative, otherwise reinterpret_cast
    // below is not going to work
    static_cast<finance::Derivative *>(static_cast<D *>(0));

    typedef shared_ptr<finance::ModelOutput> (TheoreticalModel::*ModelPriceFunction)
      (const finance::Derivative&) const;

    ModelPriceFunction mpf = reinterpret_cast<ModelPriceFunction>(pf);

    TheoreticalModel::DoRegisterPriceFunction
                      (
                        typeid(D),
                        static_cast<TheoreticalModel::PriceFunction>(mpf)
                      );
  }
#endif // __CPP2ANY__
};

} // namespace ihg


} // namespace ito33

#endif // _ITO33_IHG_THEORETICALMODEL_H_
