/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/optioninstdata.h
// Purpose:     Option instdata class, including Greek data
// Author:      David Pooley, WANG Xuewen
// Created:     2003/12/17
// RCS-ID:      $Id: optioninstdata.h,v 1.23 2005/02/01 12:16:24 nabil Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/optioninstdata.h
    @brief Option instdata class, including Greek data
*/

#ifndef _IHG_OPTIONINSTDATA_H_
#define _IHG_OPTIONINSTDATA_H_

#include "ito33/pricing/minconstraint.h"
#include "ihg/instdatawithconstraints.h"

namespace ito33
{

namespace pricing
{
  class OptionParams;
  class OptionMeshManager;
}

namespace ihg
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

  /// the constraints
  pricing::MinConstraint m_constraint;


protected:

  /// The params related to option
  pricing::OptionParams& m_optionParams;
  
  /// the option mesh manager
  pricing::OptionMeshManager& m_optionMeshes;


private:

  NO_COPY_CLASS(OptionInstData);

}; // class OptionInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_OPTIONINSTDATA_H_

