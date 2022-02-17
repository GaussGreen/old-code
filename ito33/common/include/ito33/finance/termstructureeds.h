/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructureeds.h
// Purpose:     An EDS term structure class declaration
// Author:      ITO33
// Created:     2005/02/10
// RCS-ID:      $Id: termstructureeds.h,v 1.8 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructureeds.h
    @brief Declaration of the EDS term structure class.
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTUREEDS_H_
#define _ITO33_FINANCE_TERMSTRUCTUREEDS_H_

#include "ito33/finance/eds.h"
#include "ito33/finance/termstructure.h" 

#include "ito33/dlldecl.h"

#ifndef __CPP2ANY__

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
    Class defining a term structure of equity default swaps.

    @todo Add constructor creating special structure for forward calibration
 */
class ITO33_DLLDECL TermStructureEDS : public TermStructure<EDS>
{
public:

  TermStructureEDS() : TermStructure<EDS>() { }

  /**
     Adds an EDS shared pointer to the term structure.

     Throws an exception if there is already an eds in the term structure 
     that has the same maturity

     @param eds EDS to be added to the term structure
   */
  void Add(const shared_ptr<EDS>& eds);


  /**
      Dumps all data stored in this object in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c but can also be called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

}; // class TermStructureEDS


} // namespace finance

} // namespace ito33

#else  // __CPP2ANY__ =================== Temporary workaround ==========

#include "ito33/beforestd.h"
#include <list>
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/sharedptr.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{


namespace XML
{
  class Tag;
}

namespace finance
{


template<class T>
inline bool 
HasSmallerMaturity(const shared_ptr<T>& t1, const shared_ptr<T>& t2)
{
  return t1->GetMaturityDate() >= t2->GetMaturityDate();
}

/**
    Class defining a term structure of credit default swaps.

    @noexport
 */
class ITO33_DLLDECL TermStructure_EDS
{
public:

  /// type of the container we use for storing the instruments
  typedef std::list< ito33::shared_ptr<EDS> > Elements;

  /// Default empty ctor
  TermStructure_EDS() { }

  /// Dummy virtual dtor for base class
  virtual ~TermStructure_EDS() { }

  /**
     @internal
     @brief Ctor constructs a term structure from a compatible structure.

     @param termStructure a term structure of a derived class
     
     @noexport
   */
  template <class U>
  TermStructure_EDS(const TermStructure<U>& termStructure)
  {
    const typename TermStructure<U>::Elements& 
      elements = termStructure.GetAll();

    typename TermStructure<U>::Elements::const_iterator i;
    
    for (i = elements.begin(); i != elements.end(); ++i)
      m_elements.push_back(*i);
  }

  /**
     Adds an EDS to the list.

     The function keeps the list sorted.

     The function throw an exception when a list holds already an
     element whose maturity date is same as that of the given instrument. 

     @param pElement an instrument to be added to the list
   */
  void Add(const shared_ptr<EDS>& pElement)
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

  
protected:

  /// The list of instruments
  Elements m_elements;

}; // class TermStructure_CDS

/**
    an equity default swap term structure

    @todo Add constructor creating special structure for forward calibration
 */
class ITO33_DLLDECL TermStructureEDS : public TermStructure_EDS
{
public:

  TermStructureEDS() : TermStructure_EDS() { }

  /**
     Adds an eds shared pointer to the term structure

     Throws an exception if there is already an eds in the term structure 
     that has the same maturity

     @param eds eds to be added to the term structure
   */
  void Add(const shared_ptr<EDS>& eds);

  /// @noexport
  void Visit(DerivativeVisitor& visitor) const
  {
    visitor.OnTermStructureEDS(*this);
  }


}; // class TermStructureEDS


} // namespace finance

} // namespace ito33

#endif  // __CPP2ANY__  ================================================

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTUREEDS_H_

