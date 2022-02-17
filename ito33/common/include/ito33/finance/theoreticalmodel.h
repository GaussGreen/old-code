/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/theoreticalmodel.h
// Purpose:     Abstract base class for calibrating, pricing and hedging.
// Author:      Vadim Zeitlin
// Created:     Dec 6, 2003
// RCS-ID:      $Id: theoreticalmodel.h,v 1.49 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/theoreticalmodel.h
    @brief Abstract base class for calibrating, pricing and hedging.
 */

#ifndef _ITO33_FINANCE_THEORETICALMODEL_H_
#define _ITO33_FINANCE_THEORETICALMODEL_H_

#include "ito33/debug.h"
#include "ito33/sharedptr.h"

namespace ito33
{

struct DebugParameters;

namespace XML
{
  class Tag;
}

namespace numeric
{
  class NumParams;
}

namespace finance
{

class ITO33_DLLDECL QualityControl;
class ITO33_DLLDECL ComputationalFlags;
class ITO33_DLLDECL Derivative;
class ITO33_DLLDECL ModelOutput;
class ModelParametersConsumer;

/**
    This class abstracts the model used to price derivative instruments.

    This class is somewhat too general because different models can be used to
    do different things (i.e not all models can compute prices of all
    derivatives) and in a different manner (different models produce different
    output). However we still want to present a somewhat unified API to them.

    To solve the problem of not being able to price arbitrary derivative in any
    model, we use the double dispatch mechanism. Basically, each model must
    register itself to price the derivatives it does handle.

    @rename TheoreticalModelBase
    @nocreate
 */
class ITO33_DLLDECL TheoreticalModel
{
public:

  /// Default ctor
  TheoreticalModel();

  // default copy ctor and assignment operator are ok

  /// Virtual dtor as for any base class
  virtual ~TheoreticalModel();

  /**
      Computes the price and other output quantities for the given instrument.

      @param derivative the instrument to price
      @return the generic ModelOutput object containing all outputs
      @noimpl
   */
  shared_ptr<ModelOutput> Compute(const Derivative& derivative) const
  {
    return DoCompute(derivative);
  }

  /**
      The parameters for quality control.

      @param pQualityControl parameters for quality control
   */
  void SetQualityControl(const shared_ptr<QualityControl>& pQualityControl)
  {
    m_pQualityControl = pQualityControl;
  }
  
  /**
      The parameters for quality control.

      @return parameters for quality control
   */
  const shared_ptr<QualityControl>& GetQualityControl() const
  {
    return m_pQualityControl;
  }

  /**
      @internal
      @brief Defines the external computational flags that will mainly be
             used for calibration.

      @noexport
   */
  void SetExternalFlags
       ( const shared_ptr<ComputationalFlags>& pComputFlags )
  {
    ASSERT( pComputFlags );

    m_bHasExternalFlags = true;
    m_pComputFlags = pComputFlags;
  }
  
  /**
      @internal
      @brief Defines the external computational flags as the default flags

      @noexport
   */
  void SetExternalFlagsToDefaults();

  /** 
      @internal
      @brief The computational flags determining which Greeks to
             compute, how much data to store, etc. 

      @return shared pointer to the internal compute flags for model

      @noexport
   */
  const shared_ptr<ComputationalFlags>& GetComputationalFlags() const
  {    
    return m_pComputFlags;    
  }

  /**
      @name Debugging output.

      The model can produce debug output containing the value of all input
      parameters and all outputs calculated by ComputePrice() if
      EnableDebugOutput(true) is called. The output is dumped to a default
      location (system-dependent) but can be explicitly changed by calling
      SetDebugOutputFile().
   */
  //@{

  /**
      Flag telling us whether or not to produce debugging output.

      @property DebugOutput
   */
  void EnableDebugOutput(bool bEnable = true);

  /**
      Returns the current value of "debug output" flag.

      @property DebugOutput
   */
  bool IsDebugOutputEnabled() const;

  /**
      The location of the debug output file.

      If the file name is empty, it reverts to default. If the file name is
      just the basename (i.e without path components), it is created in the
      default directory, otherwise it is used as is.

      Also calls EnableDebugOutput(true), so it is not necessary to call it in
      addition to this method.
   */
  void SetDebugOutputFile(const std::string& filename);

  /**
      The location of the debug output file.

      If EnableDebugOutput() had been called before, the value set by it is
      returned. Otherwise, the default location of the debug file for this
      platform is returned.

      This method may be called when debug output is enabled or not, i.e
      independently of the value returne by IsDebugOutputEnabled().

      @return the full path to the which is or would be used for debug output
   */
  std::string GetDebugOutputFile() const;

  //@}

  /**
      @internal 
      @brief Clones the essential part of the model, and disable all the
             flags, including debug output and computationalflags. Note that
             the new TheoreticalModel object still shares pointers (for
             example, to volatility) with original TheoreticalModel object.

      This is an ugly implementation, but this is really because model is
      intended to be used for all kinds of pricing, but pricing might have
      different needs of the flags. The better soluton is to move these
      flags outside of the model IMO.

      @noexport
   */
  virtual TheoreticalModel* Clone() const = 0;

  /**
      @internal 
      @brief Clones the essential part of the model, and disable all the flags, 
             including debug output and computationalflags. Unlike Clone(), 
             this function creates a really new copy of original 
             TheoreticalModel object.

      This is an ugly implementation, but this is really because model is
      intended to be used for all kinds of pricing, but pricing might have
      different needs of the flags. The better soluton is to move these
      flags outside of the model IMO.

      @noexport
   */
  virtual TheoreticalModel* DeepCopy() const;

  /**
      @internal
      @brief Dumps all data of this model in XML format.

      @param tagParent the parent tag under which our tag(s) should be created
     
      @noexport
   */
  virtual void Dump(ito33::XML::Tag& tagParent) const = 0;

  /**
      @internal 
      @brief Serializes the model parameters

      @noexport
   */
  virtual
  void GetModelParameters(ModelParametersConsumer& visitor) const;
  

protected:
  /**
      Dumps all data common to finance::TM.

      @param tagParent the parent tag under which our tag(s) should be created
   */
  void DumpMe(ito33::XML::Tag& tagParent) const;

  /** 
      Gets the computational flags from the derivative. 

      @param derivative the derivative
   */
  void SetFlagsFrom(const Derivative& derivative) const;

  /**
      Translates quality control to num params.

      @param derivative The instrument to price
   */
  numeric::NumParams* 
  GetNumParams(const Derivative& derivative) const;

  /// the quality control parameters
  shared_ptr<QualityControl> m_pQualityControl; 

  /**
      Boolean indicates if we should use the flags from derivative or
      the flags set externally by call to SetExternalFlags.
   */
  bool m_bHasExternalFlags;

  /**
      The flags for compute greeks etc, holding temporarily flags from
      derivative(pricing) or external(calibration etc).
   */
  mutable shared_ptr<ComputationalFlags> m_pComputFlags;

private:

  virtual shared_ptr<finance::ModelOutput> 
  DoCompute(const Derivative& derivative) const = 0;
  
  shared_ptr<DebugParameters> m_pDebug;
};

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_THEORETICALMODEL_H_
