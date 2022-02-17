/////////////////////////////////////////////////////////////////////////////
// Name:        hg/optioninstdata.h
// Purpose:     HG Option instdata class
// Created:     2005/01/13
// RCS-ID:      $Id: optioninstdata.h,v 1.6 2005/08/26 19:15:09 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/optioninstdata.h
   @brief HG Option instdata class
 */

#ifndef _HG_OPTIONINSTDATA_H_
#define _HG_OPTIONINSTDATA_H_

#include "ito33/pricing/minconstraint.h"
#include "hg/instdatawithconstraints.h"

namespace ito33
{

namespace pricing
{
  class OptionParams;
  class OptionMeshManager;
}

namespace hg
{


class OptionInstData : public InstDataWithConstraints
{

public:

  OptionInstData(pricing::OptionParams& params,
                 Model& model,
                 pricing::OptionMeshManager& meshes);

  virtual ~OptionInstData() { }

  virtual void Init();

  void UpdateBeforeStep();

  virtual void SetInitialValue();

protected:

  /// The params related to option
  pricing::OptionParams& m_optionParams;
  
  /// The option mesh manager
  pricing::OptionMeshManager& m_optionMeshes;

  /// The American constraint
  pricing::MinConstraint m_constraint;

private:

  NO_COPY_CLASS(OptionInstData);

}; // class OptionInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_OPTIONINSTDATA_H_

