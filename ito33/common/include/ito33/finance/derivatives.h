/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivatives.h
// Purpose:     A collection of derivatives
// Created:     2005/05/17
// RCS-ID:      $Id: derivatives.h,v 1.19 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivatives.h
    @brief A collection of derivatives
 */

#ifndef _ITO33_FINANCE_DERIVATIVES_H_
#define _ITO33_FINANCE_DERIVATIVES_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/// forward declaration
class ITO33_DLLDECL Derivative;
class ITO33_DLLDECL SessionData;


/**
    Derivatives is a collection of derivatives.
 */
class ITO33_DLLDECL Derivatives
{
public:
  
  /// Typedef for one element
  typedef std::pair< shared_ptr<Derivative>, double > Element;

  /// Typedef for the container
  typedef std::vector<Element> Elements;

  /// Typedef for the const_iterator
  typedef Elements::const_iterator const_iterator;

  /// Constructor contains no derivative.
  Derivatives() { }

  /// default dtor is ok

  /**
      Adds a derivative to the collection.

      @param pDerivative A shared pointer to derivative to be added
   */
  void Add(const shared_ptr<Derivative>& pDerivative) 
  { 
    CheckDerivative(pDerivative);

    m_elements.push_back( std::make_pair(pDerivative, 1.0) );
  }

  /**
      Adds a derivative to the collection with an associated weight.

      @param pDerivative A shared pointer to derivative to be added.
      @param dWeight Weighting factor.
   */
  void AddWithWeight(const shared_ptr<Derivative>& pDerivative, double dWeight) 
  { 
    CheckDerivative(pDerivative);

    CheckWeight(dWeight);

    m_elements.push_back( std::make_pair(pDerivative, dWeight) );
  }

  /**
      Iterator to the begin of derivatives.

      @return iterator to the begin of the container

      @noexport
   */
  const_iterator begin() const { return m_elements.begin(); }

  /**
      Iterator to the end of derivatives.

      @return iterator to the end of the container

      @noexport
   */
  const_iterator end() const { return m_elements.end(); }

  /**
      @internal
      @brief Gets all the derivatives in the collection.

      @return The container of the derivatives

      @noexport
   */
  const Elements& GetAll() const { return m_elements; }

  /**
      @internal
      @brief Gets all the derivatives in the collection, non const version

      @return The container of the derivatives

      @noexport
   */
  Elements& GetAll() { return m_elements; }

  /**
      @internal
      @brief Checks that all the derivatives have the same session data. Note
             that this function can only be called once the session data have
             been initialized otherwise an exception will be thrown.
 
      @noexport
   */
  void ValidateSessionDatas();

  /**
      Sets session data for each derivative inside Derivatives object

      @param pSessionData given session data
   */
  void SetSessionData(const shared_ptr<SessionData>& pSessionData);

  /**
      @internal 
      @brief Gets the session data of the derivative.

      @return The session data for all the derivatives in this collection

      @noexport
   */
  shared_ptr<SessionData> GetSessionData() const;

protected:

  /// check the derivative pointer, throw if null
  void CheckDerivative(const shared_ptr<finance::Derivative>& pDerivative);

  /// check the derivative weight, throw is non-positive
  void CheckWeight(double dWeight);

  /// The container for the derivatives
  Elements m_elements;

}; // class Derivatives


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVES_H_
