/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructureoption.h
// Purpose:     An Option term structure class 
// Created:     2005/03/04
// RCS-ID:      $Id: termstructureoption.h,v 1.8 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructureoption.h
    @brief An Option term structure class.
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTUREOPTION_H_
#define _ITO33_FINANCE_TERMSTRUCTUREOPTION_H_

#include "ito33/finance/termstructure.h"
#include "ito33/finance/option.h"

namespace ito33
{

namespace finance
{

#ifndef __CPP2ANY__

/**
    Class defining a term structure of option.
 */
class ITO33_DLLDECL TermStructureOption : public TermStructure<Option>
{
public:

  TermStructureOption() : TermStructure<Option>() { }

  /**
      Dumps all data stored in this object in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c but can also be called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

}; // class TermStructureOption

#else // __CPP2ANY__

/**
    Class for a term structure of option.  
    It often appears as input to a calibration problem.
 */
class ITO33_DLLDECL TermStructureOption
{
public:

  /// type of the container used for storing the instruments
  typedef std::list< ito33::shared_ptr<Option> > Elements;

  /// Default empty ctor
  TermStructureOption() { }

  /// Default dtor is ok

  /**
     @internal
     @brief Ctor constructs an option term structure from a compatible structure.

     @param termStructure a term structure of a derived class
     
     @noexport
   */
  template <class U>
  TermStructureOption(const TermStructure<U>& termStructure)
  {
    const typename TermStructure<U>::Elements& 
      elements = termStructure.GetAll();

    typename TermStructure<U>::Elements::const_iterator i;
    
    for (i = elements.begin(); i != elements.end(); ++i)
      m_elements.push_back(*i);
  }

  /**
     Adds an option to the list.

     The function keeps the list sorted.

     The function throw an exception when a list holds already an
     element whose maturity date is same as that of the given instrument. 

     @param pElement an instrument to be added to the list
   */
  void Add(const shared_ptr<Option>& pElement)
  {  
    if ( !pElement ) 
      ThrowExceptionWhenAddingNull();

    Elements::iterator
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
     @brief Gets the list of instruments, it doesn't throw, it's up to the caller to 
     ensure the required size.

     @return the list of instruments
     
     @noexport
   */
  const Elements& GetAll() const { return m_elements; }

  /**
      @internal
      @brief Dumps all data stored in this object in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c but can also be called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

  /// @noexport
  void Visit(DerivativeVisitor& visitor) const;


protected:

  /// The list of instruments
  Elements m_elements;

}; // class TermStructureOption

#endif  // __CPP2ANY__


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTUREOPTION_H_

