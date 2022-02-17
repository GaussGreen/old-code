/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructurederivative.h
// Purpose:     A general derivative term structure class declaration
// Created:     2004/06/14
// RCS-ID:      $Id: termstructurederivative.h,v 1.8 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructurederivative.h
    @brief Declaration of a general derivative term structure class.
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTUREDERIVATIVE_H_
#define _ITO33_FINANCE_TERMSTRUCTUREDERIVATIVE_H_

#include "ito33/finance/termstructure.h"
#include "ito33/finance/derivative.h"

namespace ito33
{

namespace XML
{
  class Tag;
}


namespace finance
{


/**
    a general derivative term structure
 */
typedef TermStructure<Derivative> TermStructureDerivative;


/** 
   @internal
   @brief Check that the number of elements corresponds to the one actually 
          requested.

   @param tsDerivative The term structure whose size will be checked
   @param nSize number of element requested
 */
void CheckSize(const TermStructure<Derivative>& tsDerivative,
               size_t nSize);

/**
   @internal
   @brief Sets the session data for each derivative inside the term structure.

   @param termStructureDerivative The term structure holding the derivatives
   @param pSessionData the session data to be set to each derivative.

   @noexport
 */
void 
SetTermStructureSessionData
(const TermStructure<Derivative>& tsDerivative,
 const shared_ptr<SessionData>& pSessionData);

/**
    @internal
    @brief Dumps a term structure of derivatives in XML format.

    This method is usually called by the function doing the pricing,
    calibration &c but can also be called "manually" if needed.

    @param tsDerivative The term structure holding the derivatives
    @param tagParent the parent tag under which our tag(s) should be created
    @param sXmlTagName the tag name for the term structure
    
    @noexport
  */
void DumpTermStructureDerivative(const TermStructure<Derivative>& tsDerivative,
                                 XML::Tag& tagParent, std::string sXmlTagName);
} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTUREDERIVATIVE_H_
