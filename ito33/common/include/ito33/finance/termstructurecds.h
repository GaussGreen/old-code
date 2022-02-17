/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/termstructurecds.h
// Purpose:     A cds term structure class declaration
// Author:      ZHANG Yunzhi
// Created:     05.06.2004
// RCS-ID:      $Id: termstructurecds.h,v 1.27 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/termstructurecds.h
    @brief Declaration of the cds term structure class.
 */

#ifndef _ITO33_FINANCE_TERMSTRUCTURECDS_H_
#define _ITO33_FINANCE_TERMSTRUCTURECDS_H_

#include "ito33/finance/cdslike.h"
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
    Class defining a term structure of credit default swaps.

    @todo Add constructor creating special structure for forward calibration
 */
class ITO33_DLLDECL TermStructureCDS : public TermStructure<CDSLike>
{
public:

  TermStructureCDS() : TermStructure<CDSLike>() { }

  /**
     Adds a CDSLike shared pointer to the term structure.

     Throws an exception if there is already a cds in the term structure 
     that has the same maturity

     @param cds CDSLike to be added to the term structure
   */
  void Add(const shared_ptr<CDSLike>& cds);


  /**
      Dumps all data stored in this object in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c but can also be called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

}; // class TermStructureCDS


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
class ITO33_DLLDECL TermStructure_CDS
{
public:

  /// type of the container we use for storing the instruments
  typedef std::list< ito33::shared_ptr<CDSLike> > Elements;

  /// Default empty ctor
  TermStructure_CDS() { }

  /// Dummy virtual dtor for base class
  virtual ~TermStructure_CDS() { }

  /**
     @internal
     @brief Ctor constructs a term structure from a compatible structure.

     @param termStructure a term structure of a derived class
     
     @noexport
   */
  template <class U>
  TermStructure_CDS(const TermStructure<U>& termStructure)
  {
    const typename TermStructure<U>::Elements& 
      elements = termStructure.GetAll();

    typename TermStructure<U>::Elements::const_iterator i;
    
    for (i = elements.begin(); i != elements.end(); ++i)
      m_elements.push_back(*i);
  }

  /**
     Adds a CDSLike to the list.

     The function keeps the list sorted.

     The function throw an exception when a list holds already an
     element whose maturity date is same as that of the given instrument. 

     @param pElement an instrument to be added to the list
   */
  void Add(const shared_ptr<CDSLike>& pElement)
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
    a credit default swap term structure

    @todo Add constructor creating special structure for forward calibration
 */
class ITO33_DLLDECL TermStructureCDS : public TermStructure_CDS
{
public:

  TermStructureCDS() : TermStructure_CDS() { }

  /**
     Adds a CDSLike shared pointer to the term structure.

     Throws an exception if there is already a cds in the term structure 
     that has the same maturity

     @param cds cds to be added to the term structure
   */
  void Add(const shared_ptr<CDSLike>& cds);

  /// @noexport
  void Visit(DerivativeVisitor& visitor) const
  {
    visitor.OnTermStructureCDS(*this);
  }

}; // class TermStructureCDS


} // namespace finance

} // namespace ito33

#endif  // __CPP2ANY__  ================================================

#endif // #ifndef _ITO33_FINANCE_TERMSTRUCTURECDS_H_

