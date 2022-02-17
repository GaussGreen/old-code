/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/termstructurederivative.cpp
// Purpose:     A general derivative term structure class declaration
// Author:      Ito33
// Created:     2004/12/14
// RCS-ID:      $Id: termstructurederivative.cpp,v 1.8 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file common/src/finance/termstructurederivative.h
    @brief Declaration of  a general derivative term structure class.
 */

#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/termstructurederivative.h"

#include "ito33/xml/write.h"

extern const ito33::finance::Error
  ITO33_TERMSTRUCTURE_ADD_NULL,
  ITO33_TERMSTRUCTURE_ADD_DUPLICATED,
  ITO33_INCONSISTENT_SESSION_DATA,
  ITO33_TERMSTRUCTURE_INCOMPLETE;

namespace ito33
{

namespace finance
{

void ThrowExceptionWhenAddingNull()
{
  throw EXCEPTION(ITO33_TERMSTRUCTURE_ADD_NULL);
}

void ThrowExceptionWhenAddingDuplicated()
{
  throw EXCEPTION(ITO33_TERMSTRUCTURE_ADD_DUPLICATED);
}

void ThrowInconsistentSessionDatas()
{
  throw EXCEPTION(ITO33_INCONSISTENT_SESSION_DATA);
}

void CheckSize(const TermStructure<Derivative>& tsDerivative,
                size_t nSize) 
{
  CHECK_COND( tsDerivative.GetAll().size() >= nSize,
              ITO33_TERMSTRUCTURE_INCOMPLETE );
}

void 
SetTermStructureSessionData(const TermStructure<Derivative>& tsDerivative, 
 const shared_ptr<SessionData>& pSessionData)
{
  const TermStructure<Derivative>::Elements& 
    derivatives = tsDerivative.GetAll();

  TermStructure<Derivative>::Elements::const_iterator i;
  for ( i = derivatives.begin(); i != derivatives.end(); ++i)
    (*i)->SetSessionData(pSessionData);
}


void DumpTermStructureDerivative(const TermStructure<Derivative>& tsDerivative,
                                 XML::Tag& tagParent, std::string sXmlTagName)
{
  XML::Tag tagTS(sXmlTagName, tagParent);

  const TermStructureDerivative::Elements& elements(tsDerivative.GetAll());
  TermStructureDerivative::Elements::const_iterator i;
  for (i = elements.begin(); i != elements.end(); ++i)
    (*i)->Dump(tagTS);
}

} // namespace finance

} // namespace ito33
