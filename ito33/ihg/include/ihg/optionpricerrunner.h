/****************************************************************************
 * Name:        ihg/optionpricerrunner.h
 * Purpose:     option pricer runner class
 * Author:      ZHANG Yunzhi
 * Created:     2004/07/26
 * RCS-ID:      $Id: optionpricerrunner.h,v 1.2 2004/10/04 18:04:04 pedro Exp $
 * Copyright:   (c) 2003-2003 Trilemma LLP
 ****************************************************************************/

/**
    @file ihg/optionpricerrunner.h
    @brief Option pricer runner class

    Implementation of the pricer class for options.
*/

#ifndef _IHG_OPTIONPRICERRUNNER_H_ 
#define _IHG_OPTIONPRICERRUNNER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

#include "ihg/optioninstdata.h"
#include "ihg/optionstepper.h"
#include "ihg/optionnumoutput.h"

namespace ito33
{

namespace ihg
{

class OptionPricerRunner
{
public:
  OptionPricerRunner(
                      pricing::OptionParams &params,
                      pricing::OptionMeshManager &meshes,
                      OptionInstData &instdata,
                      OptionStepper &stepper,
                      OptionNumOutput &numoutput)
                      : m_engine(params, meshes, instdata, stepper, numoutput)
  {
  }

  void Run()
  {
    m_engine.Run();
  }

private:
  OptionEngine m_engine;

private:

  NO_COPY_CLASS(OptionPricerRunner);

}; // class OptionPricerRunner


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_OPTIONPRICERRUNNER_H_

