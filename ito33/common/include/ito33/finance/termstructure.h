/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructure.h
// Purpose:     TermStructure template class declaration
// Author:      Wang
// Created:     2004/06/10
// RCS-ID:      $Id: termstructure.h,v 1.23 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructure.h
    @brief class represents a term structure of an instrument 
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTURE_H_
#define _ITO33_FINANCE_TERMSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <list>
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{

void ThrowExceptionWhenAddingDuplicated();
void ThrowExceptionWhenAddingNull();
void ThrowInconsistentSessionDatas();

template<class T>
inline bool 
HasSmallerMaturity(const shared_ptr<T>& t1, const shared_ptr<T>& t2)
{
  return t1->GetMaturityDate() >= t2->GetMaturityDate();
}

/**
    A term structure is a set of instruments (same type) with increasing and 
    distinct maturities.  It often appears as input to a calibration problem.
 */
template<class T>
class TermStructure
{
public:

  /// type of the container we use for storing the instruments
  typedef std::list< ito33::shared_ptr<T> > Elements;

  /// Default empty ctor
  TermStructure() { }

  /// Dummy virtual dtor for base class
  virtual ~TermStructure() { }

  /**
      @internal
      @brief Ctor constructs a term structure from a compatible structure

      @param termStructure a term structure of a derived class
     
      @noexport
   */
  template <class U>
  TermStructure(const TermStructure<U>& termStructure)
  {
    const typename TermStructure<U>::Elements& 
      elements = termStructure.GetAll();

    typename TermStructure<U>::Elements::const_iterator i;
    
    for (i = elements.begin(); i != elements.end(); ++i)
      m_elements.push_back(*i);
  }

  /**
      Adds an instrument to the term structure list.

      The function keeps the list sorted.

      The function throw an exception when a list holds already an
      element whose maturity date is same as that of the given instrument. 

      @param pElement an instrument to be added to the list
   */
  void Add(const shared_ptr<T>& pElement)
  {  
    if ( !pElement )
      ThrowExceptionWhenAddingNull();

    typename Elements::iterator
      iterResult = std::search_n( m_elements.begin(), 
                                  m_elements.end(),
                                  1,
                                  pElement,
                                  HasSmallerMaturity<T>);

    if ( iterResult != m_elements.end() &&
         (*iterResult)->GetMaturityDate() == pElement->GetMaturityDate() )
    {
      ThrowExceptionWhenAddingDuplicated();
    }

    m_elements.insert(iterResult, pElement);
  }

  /**
      @internal
      @brief Gets the list of instruments.

      @return the list of instruments
     
      @noexport
   */
  const Elements& GetAll() const { return m_elements; }

  /**
      @internal
      @brief Checks that all elements of the term structure are of the same
      type. 
    
      @return true/false 
      @noexport
   */
  bool Validate() const
  { 
    typename Elements::const_iterator iter;

    for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
    {
      if ( !dynamic_pointer_cast<T>(*iter) )
        return false;
    }

    return true;
  }
 
  /**
      @internal
      @brief Checks that all elements of the term structure have the
             same session data. Note that this function can only be called once
             the session data have been initialized otherwise an exception will
             be thrown.
 
      @noexport
   */
  void ValidateSessionDatas()
  {
    //list empty or contains only one element
    if ( m_elements.size() <= 1 )
      return;

    typename Elements::const_iterator iter = m_elements.begin();
    shared_ptr<SessionData> pSessionData( (*iter)->GetSessionData() );
    ++iter; //increment pointer to start past first element

    for (iter; iter != m_elements.end(); ++iter)
      if ( pSessionData != (*iter)->GetSessionData() )
        ThrowInconsistentSessionDatas();
  }

protected:

  /// The list of instruments
  Elements m_elements;

}; // class TermStructure


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTURE_H_
