/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/jumps.h
// Purpose:     Collection of jumps from one regime to another
// Created:     2005/01/13
// RCS-ID:      $Id: jumps.h,v 1.6 2005/06/08 14:44:13 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/hg/jumps.h
   @brief Collection of jumps from one regime to another.
 */
#ifndef _ITO33_HG_JUMPS_H_
#define _ITO33_HG_JUMPS_H_

#include "ito33/vector.h"

#include "ito33/hg/jump.h"

namespace ito33
{

namespace hg
{


/**
   Collection of jumps from one regime to another.

   @indexed_iterator GetAll
 */
class ITO33_HG_DLLDECL Jumps
{
public:

  /// std container used as storage
  typedef std::vector<Jump> Elements;
  
  /// Embedded type for iteration on the container
  typedef Elements::const_iterator const_iterator;

  /// Default ctor contains no jump
  Jumps() { }

  /**
     Default ctor from an existing container.

     @jumps The collection of jumps as a std container
   */
  Jumps(const Elements& jumps) : m_elements(jumps) { }

  /**
     Adds a jump.

     @param dIntensity The instensity of the jump
     @param dAmplitude The amplitude(size) of the jump
   */
  void Add(double dIntensity, double dAmplitude)
  {
    m_elements.push_back( Jump(dIntensity, dAmplitude) );
  }

  /**
     Adds a jump.

     @param jump A jump to be added to the collection
     @noexport 
   */
  void Add(const Jump& jump) { m_elements.push_back(jump); }

  /**
     Sets at once all the jumps.
   
     @param jumps A vector of jumps

     @noexport
   */
  void SetAll(const Elements& jumps) { m_elements = jumps; }
   
  /** 
     Gets all the jumps as a vector.

     @return all the jumps in the collection

     @noexport COM
   */
  const Elements& GetAll() const { return m_elements; } 
  
  /** 
     Removes all the exiting jumps.
   */
  void Clear() { clear(); }  
  
  /// @name std container method
  //@{

  /**
     Adds a jump.

     @noexport
   */
  void push_back(const Jump& jump) { m_elements.push_back(jump); }

  /** 
     Removes all the existing jumps.

     @noexport
   */
  void clear() { m_elements.clear(); }

  /**
     @internal

     @brief Gets the iterator points to the first cash flow

     @noexport
   */
  const_iterator begin() const { return m_elements.begin(); }

  /**
     @internal

     @brief Returns an iterator that addresses the location succeeding
            the last cash flow.

     @noexport
   */
   const_iterator end() const { return m_elements.end(); }

  /**
     @internal

     @brief Returns the number of cash flows.

     @noexport
   */
  size_t size() const { return m_elements.size(); }

  /**
     @internal

     @brief Check if the collection is empty.

     @noexport
   */
  bool empty() const { return m_elements.empty(); }

  //@}


private:

  /// The jumps
  Elements m_elements;

}; // class Jumps


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_JUMPS_H_
