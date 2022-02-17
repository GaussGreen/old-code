/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testparametrization.h
// Purpose:     testing for parametrization visitors
// Author:      Ito33 Canada
// Created:     08/08/2005
// RCS-ID:      $Id: testparametrization.h,v 1.2 2005/08/15 20:04:11 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <string.h>
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/date.h"

#include "ito33/finance/modelparametersconsumer.h"

namespace ito33
{
namespace finance
{
 
  
class ModelParametersConsumerTest : public ModelParametersConsumer 
{
 
public:

  void OnScalarValues
       (
         const std::string& categoryName,
         const std::vector<ScalarModelParameter>& parameters
       )      
  {     
    m_sCategoryName = categoryName;
    m_pParameters   = parameters;    
  }
     
  void OnTimeComponentValues
    (
      const std::string& categoryName,
      const std::vector<Date>& dates,
      const std::vector<double>& values
    )      
  {
    m_sCategoryName = categoryName;
    m_pDates = dates;
    m_pdVal  = values;     
  }

  std::vector<double> m_pdVal;
  std::vector<Date> m_pDates;
  std::string m_sCategoryName;
  std::vector<ScalarModelParameter> m_pParameters;


}; 


} //namespace finance
} //namespace ito33

