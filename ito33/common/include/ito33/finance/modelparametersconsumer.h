/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/modelparametersconsumer.h
// Purpose:     Class helps to access model parameters
// Created:     2005/08/11
// RCS-ID:      $Id: modelparametersconsumer.h,v 1.8 2006/01/10 16:53:16 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/modelparametersconsumer.h
    @brief ModelParametersConsumer class declaration
 */

#ifndef _ITO33_FINANCE_MODELPARAMETERSCONSUMER_H_
#define _ITO33_FINANCE_MODELPARAMETERSCONSUMER_H_

#include "ito33/vector.h"
#include "ito33/string.h"

namespace ito33
{
  class ITO33_DLLDECL Date;

namespace finance
{

#define MODEL_PARAM_NAME_VOL_FLAT          "Flat volatility"
#define MODEL_PARAM_NAME_VOL_TIMEONLY      "Time-dependent only volatility"
#define MODEL_PARAM_NAME_VOL_TANH          "Tanh volatility"
#define MODEL_PARAM_NAME_VOL_POWER         "Power volatility"

#define MODEL_PARAM_NAME_HR_FLAT           "Flat hazard rate"
#define MODEL_PARAM_NAME_HR_TIMEONLY       "Time-dependent only hazard rate"
#define MODEL_PARAM_NAME_HR_POWER          "Power hazard rate"
#define MODEL_PARAM_NAME_HR_TIMESPOT       "Time-and-spot-dependent hazard rate"

#define SCALAR_MODEL_PARAM_NAME_VALUE        "Value"
#define SCALAR_MODEL_PARAM_NAME_BETA         "Beta"
#define SCALAR_MODEL_PARAM_NAME_ALPHA        "Alpha"
#define SCALAR_MODEL_PARAM_NAME_LEFT_LIMIT   "Left limit"
#define SCALAR_MODEL_PARAM_NAME_RIGHT_LIMIT  "Right limit"
#define SCALAR_MODEL_PARAM_NAME_SCALE        "Scale"
#define SCALAR_MODEL_PARAM_NAME_S0           "S0"

/**
  struct holding all information of a scalar parameter

  @noexport
  */
struct ScalarModelParameter
{
  ScalarModelParameter() : expressedInPercents(true) {}

  /// Name of the parameter, like "beta"
  std::string name;

  /// Value of the parameter
  double value;

  /// Whether the value should be displayed as percentage
  bool expressedInPercents;
};


/**    
     Class helps to get calibration results.

     @noexport
 */
class ModelParametersConsumer 
{
public:


  virtual ~ModelParametersConsumer() {}

  /**
    Handles a collection of ScalarModelParameters belong to a "subset"
    of model parameters qualified by categoryName.

    @param categoryName like "volatility flat" or "harzard rate time only"
    @param parameters collection of scalar model parameters
    */
  virtual void OnScalarValues
               (
                 const std::string& categoryName,
                 const std::vector<ScalarModelParameter>& parameters
               ) = 0;

  virtual void OnTimeComponentValues
               (
                 const std::string& categoryName,
                 const std::vector<Date>& dates,
                 const std::vector<double>& values
               ) = 0;
};
  

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_MODELPARAMETERSCONSUMER_H_
