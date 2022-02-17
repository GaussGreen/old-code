/////////////////////////////////////////////////////////////////////////////
// Name:        hg/instdatatimeonly.h
// Purpose:     instdata class for time only problem with HG model
// Created:     2005/06/09
// RCS-ID:      $Id: instdatatimeonly.h,v 1.3 2006/03/30 15:48:15 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/instdatatimeonly.h

   @remark As a time only equation, the method used here to solve the system 
           is different from others. The recovery value at the previuous time 
           step is used and so needs to be setup correctly.
 */

#ifndef _HG_INSTDATATIMEONLY_H_
#define _HG_INSTDATATIMEONLY_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{
  class Params;
  class MeshManager;
}

namespace hg
{

  class Model;

/// Struct to help computing sensitivity
struct TimeOnlySensitivityData
{
  bool m_bIsJumpToDefault;
  size_t m_nIdxR1;
  size_t m_nIdxR2;
  double m_dIntensity;
};

/// Inst data class for time only problem in the HG model
class InstDataTimeOnly : public pricing::InstDataTimeOnly
{
public:

  InstDataTimeOnly(pricing::Params& params,
                   Model& model,
                   pricing::MeshManager& meshes);

  // Default dtor is ok

  void Init();

  void SetInitialValue();
  
  virtual void UpdateBeforeStep();

  virtual void UpdateRecoveryValue() = 0;

  /** 
     Normally it's implemented in Stepper, but it's just more convient
     to do it here
   */
  void Run();

  std::vector<double> GetSensitivities() const;
 
  /// The number of regimes
  size_t m_nNbRegimes;

  /// the price at this time
  Array<double> m_pdPrices;

  /// the price at the old time step
  Array<double> m_pdOldPrices;
  
  /// Current default value
  double m_dRecoveryValue;

  /// Flags indicate if we compute sensitivities
  std::vector<bool> m_pbComputeSensitivities;

  size_t m_nNbSensitivities;

  Array< Array<double> > m_ppdSensitivities;

  Array< Array<double> > m_ppdOldSensitivities;


protected:
    
  void ApplyEvent(const pricing::Event* pEvent); 

  Model& m_model;

  /// the jumps to default
  std::vector<double> m_pdJumpsToDefault;
  
  /// The helper matrix
  double m_ppdA[3][3];
  
  /// The old recovery value
  double m_dOldRecoveryValue;
  
  /// old interest rate (yield)
  double m_dOldRate;


private:

  std::vector<TimeOnlySensitivityData> m_sensitivityDatas;
  
  NO_COPY_CLASS(InstDataTimeOnly);

}; // class InstDataTimeOnly


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CDSINSTDATA_H_
