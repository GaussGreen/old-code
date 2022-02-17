/////////////////////////////////////////////////////////////////////////////
// Name:        hg/heroinstdata.h
// Purpose:     HG HERO instdata class
// Created:     2005/09/26
// RCS-ID:      $Id: heroinstdata.h,v 1.2 2005/12/01 18:43:00 dave Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/heroinstdata.h
   @brief HG HERO instdata class
 */

#ifndef _HG_HEROINSTDATA_H_
#define _HG_HEROINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace hg
{

class HeroParams;
class HeroMeshManager;

class HeroInstData : public BackwardInstData
{

public:

  HeroInstData(HeroParams& params,
               Model& model,
               HeroMeshManager& meshes);

  virtual ~HeroInstData() { }

  virtual void Init();

  void UpdateBeforeStep();

  virtual void SetInitialValue();

  /// The hero PDE terms at this timestep
  std::vector<double> m_pdHeroTerms;

protected:

  /// The HERO params
  HeroParams& m_heroParams;
  
  /// The HERO mesh manager
  HeroMeshManager& m_heroMeshes;  

private:

  NO_COPY_CLASS(HeroInstData);

}; // class HeroInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_HEROINSTDATA_H_

