/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructure_enumerator.h
// Purpose:     Visitor for termstructures
// Created:     2004-05-11
// RCS-ID:      $Id: termstructure_enumerator.h,v 1.7 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructure_enumerator.h
    @brief Visitor for termstructures
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTRUE_ENUMERATOR_
#define _ITO33_FINANCE_TERMSTRUCTRUE_ENUMERATOR_

#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/termstructureoption.h"

namespace ito33
{

namespace finance
{

/**
    TermStructureEnumerator get shared ptr to the right type of the object.
 */
class TermStructureEnumerator
{
public:
  
  /**
      Gets the term structure as a termstructure of derivative.

      We are implictly supposing that there is only one term structure.
      
      @return The term structure of the derivative, might be empty if there
              isn't any term structure
   */
  TermStructureDerivative GetTermStructure() const
  {
    if ( m_pTermStructureCDS )
      return *m_pTermStructureCDS;
    else if ( m_pTermStructureParBond )
      return *m_pTermStructureParBond;
    else if ( m_pTermStructureEDS )
      return *m_pTermStructureEDS;
    else if ( m_pTermStructureOption )
      return *m_pTermStructureOption;

    return TermStructureDerivative();
  }

  /// Gets visited CDS term structure, otherwise returns 0
  shared_ptr<TermStructureCDS> GetTermStructureCDS() const
  { 
    return m_pTermStructureCDS;
  } 

  /// Gets visited ParBond term structure, otherwise returns 0
  shared_ptr<TermStructureParBond> GetTermStructureParBond() const
  { 
    return m_pTermStructureParBond;
  } 

  /// Gets visited EDS term structure, otherwise returns 0
  shared_ptr<TermStructureEDS> GetTermStructureEDS() const
  { 
    return m_pTermStructureEDS;
  } 

  /// Gets visited option term structure, otherwise returns 0
  shared_ptr<TermStructureOption> GetTermStructureOption() const 
  { 
    return m_pTermStructureOption;
  } 
 
public:
  
  /// Sets the session data of the derivatives inside the term structures 
  void SetSessionData(const shared_ptr<SessionData>& pSessionData)
  {
    if ( m_pTermStructureCDS )
      SetTermStructureSessionData(*m_pTermStructureCDS, pSessionData);
       
    if ( m_pTermStructureParBond )
      SetTermStructureSessionData(*m_pTermStructureParBond, pSessionData);

    if ( m_pTermStructureEDS )
      SetTermStructureSessionData(*m_pTermStructureEDS, pSessionData);
                   
    if ( m_pTermStructureOption )
      SetTermStructureSessionData(*m_pTermStructureOption, pSessionData);
  }

  /// Called for a cds term structure
  void SetTermStructureCDS(const shared_ptr<TermStructureCDS>& tsCDS)  
  {   
    m_pTermStructureCDS = tsCDS;
  } 

  /// Called for a cds term structure
  void SetTermStructureParBond(const shared_ptr<TermStructureParBond>& tsParBond)  
  {   
    m_pTermStructureParBond = tsParBond;
  } 

  /// Called for a eds term structure
  void SetTermStructureEDS(const shared_ptr<TermStructureEDS>& tsEDS)  
  {   
    m_pTermStructureEDS = tsEDS;
  } 

  /// Called for an option term structure
  void SetTermStructureOption(const shared_ptr<TermStructureOption>& tsOption)  
  {   
    m_pTermStructureOption = tsOption;
  } 


private:

  shared_ptr<TermStructureCDS> m_pTermStructureCDS;

  shared_ptr<TermStructureParBond> m_pTermStructureParBond;

  shared_ptr<TermStructureEDS> m_pTermStructureEDS;

  shared_ptr<TermStructureOption> m_pTermStructureOption;

};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTRUE_ENUMERATOR_
