/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/arraycheckers.cpp
// Purpose:     Some useful checkers for array
// Author:      Wang
// Created:     2004/06/04
// RCS-ID:      $Id: arraycheckers.cpp,v 1.4 2004/10/05 09:13:42 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/arraycheckers.h"

#include "ito33/useexception.h"
#include "ito33/error.h"

extern const ito33::Error ITO33_ARRAY_NONPOSITIVE;
extern const ito33::Error ITO33_ARRAY_NEGATIVE;

namespace ito33
{


const double* 
CheckNonNegativity(const double *pdValues,
                   size_t nNumber,
                   const std::string& pcArrayName)
{
  CheckArraySize(nNumber, pcArrayName);

  for (size_t nIdx = 0; nIdx < nNumber; nIdx++)
    if(pdValues[nIdx] < 0)
    {
      if(pcArrayName == "")
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NEGATIVE,
                ITO33_ARRAY_NEGATIVE.GetMessage()
                );
      else
      {
        std::string msg(TRANS(pcArrayName.c_str()));
        msg += ": ";
        msg += ITO33_ARRAY_NEGATIVE.GetMessage();
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NEGATIVE,
                msg.c_str()
                );
      }
    }

  return pdValues;
}

const std::vector<double>& 
CheckNonNegativity(const std::vector<double>& pdValues,
                   const std::string& pcArrayName)
{
  CheckNonNegativity(&pdValues[0], pdValues.size(), pcArrayName);
  
  return pdValues;
}

const std::vector<double>& 
CheckPositivity(const std::vector<double>& pdValues,
                const std::string& pcArrayName)
{
  CheckPositivity(&pdValues[0], pdValues.size(), pcArrayName);
  
  return pdValues;
}

const double* 
CheckPositivity(const double *pdValues,
                size_t nNumber,
                const std::string& pcArrayName)
{
  CheckArraySize(nNumber, pcArrayName);

  for (size_t nIdx = 0; nIdx < nNumber; nIdx++)
    if(pdValues[nIdx] <= 0)
    {
      if(pcArrayName == "")
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NONPOSITIVE,
                ITO33_ARRAY_NONPOSITIVE.GetMessage()
                );
      else
      {
        std::string msg(TRANS(pcArrayName.c_str()));
        msg += ": ";
        msg += ITO33_ARRAY_NONPOSITIVE.GetMessage();
        throw EXCEPTION_MSG
                (
                ITO33_ARRAY_NONPOSITIVE,
                msg.c_str()
                );
      }
    }

  return pdValues;
}



} // namespace ito33
