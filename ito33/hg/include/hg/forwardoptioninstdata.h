/////////////////////////////////////////////////////////////////////////////
// Name:        hg/forwardoptioninstdata.h
// Purpose:     HG forward option instdata class
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptioninstdata.h,v 1.3 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/forwardoptioninstdata.h
   @brief HG forward option instdata class
 */

#ifndef _HG_FORWARDOPTIONINSTDATA_H_
#define _HG_FORWARDOPTIONINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class ForwardOptionParams;
  class ForwardOptionMeshManager;
}

namespace hg
{


class ForwardOptionInstData : public BackwardInstData
{

public:

  ForwardOptionInstData(pricing::ForwardOptionParams& params,
                        Model& model,
                        pricing::ForwardOptionMeshManager& meshes);

  virtual ~ForwardOptionInstData() { }

  virtual void Init();

  void UpdateBeforeStep();

  virtual void SetInitialValue();

  void SetupFlags(const finance::ComputationalFlags& flags);


protected:

  /// The params related to option
  pricing::ForwardOptionParams& m_forwardOptionParams;
  
  /// the option mesh manager
  pricing::ForwardOptionMeshManager& m_forwardOptionMeshes;


private:

  NO_COPY_CLASS(ForwardOptionInstData);

}; // class ForwardOptionInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_FORWARDOPTIONINSTDATA_H_
