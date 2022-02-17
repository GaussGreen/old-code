/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivative.h
// Purpose:     Base derivative class declaration
// Created:     16.12.2003
// RCS-ID:      $Id: derivative.h,v 1.56 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivative.h
    @brief Declaration of the base derivative class.

    This class includes 
      - a shared pointer to the session data associated to the derivative
      - the market price of the derivative.
 */

#ifndef _ITO33_FINANCE_DERIVATIVE_H_
#define _ITO33_FINANCE_DERIVATIVE_H_

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL SessionData;
class ITO33_DLLDECL Numeraire;
class ITO33_DLLDECL Issuer;
class ITO33_DLLDECL ComputationalFlags;

class DerivativeVisitor;
class DerivativeModifyingVisitor;

/**
    Derivative is the base class for derivative instruments.

    @nocreate
 */
class ITO33_DLLDECL Derivative
{
public:
  /**
      Constructor for a derivative.

      Currency might be null which means that this derivative is
      in the same currency as its underlying.
   */
  Derivative();

  /// virtual destructor since Derivative can be used polimorphically
  virtual ~Derivative() { }

  /**
      @name Methods for initializing the derivative.
   */
  //@{

  /**
      The session data associated with this derivative.
 
      Derivatives with another derivative as underlying (CBOption for example)
      should reimplement this function to define as well the session data of
      the underlying derivative.

      @param pSessionData session data of this derivative
   */
  virtual void SetSessionData(const shared_ptr<SessionData>& pSessionData);
  
  /**
      The currency of this derivative.

      @param pCurrency the currency
   */
  void SetNumeraire(const shared_ptr<Numeraire>& pCurrency);

  /**
      The issuer of this derivative, might be null.

      @param pIssuer the issuer of this derivative
   */
  void SetIssuer(const shared_ptr<Issuer>& pIssuer)
  {
    m_pIssuer = pIssuer;
  }

  /*
    Note: CDS price can be positive and negative as well. 
          RefCDS price can be set by the spread. 
          In these cases SetPrice should be implemented differently. 
          That is why this function is virtual.
  */
  /**
      The market price of this derivative.

      @param dPrice market price
   */
  virtual void SetMarketPrice(double dPrice);

  /** 
      The computational flags determining which Greeks to
      compute, how much data to store, etc. 

      @param pFlags a shared pointer to a computational flags
   */
  void SetComputationalFlags(const shared_ptr<ComputationalFlags>& pFlags);

  //@}  // end @name Methods for initializing the Derivative.

  /**
      @name Methods for accessing the Derivative.
   */
  //@{

  /**
      The session data associated with this derivative.

      @return The session data associated with the derivative
   */
  const shared_ptr<SessionData>& GetSessionData() const;

  /**
      The issuer of this derivative, might be null.

      @return The issuer of this derivative, if any
   */
  const shared_ptr<Issuer>& GetIssuer() const
  {
    return m_pIssuer;
  }

  /*
    This function is virtual because market price of some instruments can be
    quoted differently than market price. 
  */
  /**
      Checks if the market price is set.
      
      The user should check this function before
      calling GetMarketPrice(). @sa GetMarketPrice().

      @return true if the market price is set, false otherwise.
   */
  virtual bool HasMarketPrice() const;

  /**
      Returns the market price of this derivative or throws if it
      is not available.

      The market price is either previously set by SetMarketPrice() or 
      indirectly defined by other functions, such as SetImpliedVol() 
      for an option.

      The usage of this function should look like this:
      @code
        ...
        if ( derivative.HasMarketPrice() )
          price = derivative.GetMarketPrice();
        ...
      @endcode 

      @return The market price of the derivative.
   */
  double GetMarketPrice() const 
  {
    CheckMarketPrice();
    return DoGetMarketPrice(); 
  }

  /** 
      Returns the maturity date of this derivative.

      @return The maturity date of this derivative
   */
  virtual Date GetMaturityDate() const = 0;

  /** 
      The computational flags determining which Greeks to
      compute, how much data to store, etc. 

      Might be null which means that default values will be computed.

      @return shared pointer to the internal computational flags
   */
  const shared_ptr<ComputationalFlags>& GetComputationalFlags() const
  {
    return m_pFlags;    
  }

  //@} 
  
  /**
      @name Functions useful for cross-currency instruments.
   */
  //@{

  /**
      Checks if the derivative and the underlying share are denominated 
      in different currencies. 
      
      Throws exception if the session data is not set.

      @return true if it is a cross-currency instrument.
   */
  bool IsCrossCurrency() const;

  /**
      The currency of this derivative.

      A null object is returned if the currency has not been previously set.

      @return The currency of this derivative
   */
  const shared_ptr<Numeraire>& GetNumeraire() const
  {
    return m_pNumeraire;
  }

  //@}

  /**
      Validates the contents of this derivative. The consistency between the
      derivative and its session data is not verified.

      Default implementation assumes that the derivative doesn't need 
      additional validation.
   */
  virtual void Validate() const { }

  /**
      Validates both the derivative and the consistency between the derivative
      and its session data.

      The session data of this derivative must be set.
   */
  virtual void ValidateAll() const;

  /**
      Validates both the derivative and the consistency between the derivative
      and the given session data.

      Useful only when the session data is not set in the derivative.

      @param sessionData The session data of this derivative
   */
  virtual void ValidateWith(const SessionData& sessionData) const;

  /**
      @internal      
      @brief Perturbs the FX rate. The derivative must have been already
             validated with its session data.

      This function is used to compute the FXDelta.

      @param dFXRateShift the perturbation parameter, its meaning varies for
                          different derived classes.
      
      @noexport
   */
  virtual void PerturbFXRate(double dFXRateShift) const;

  /**
      @internal
      @brief This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      derivatives without adding virtual functions to do all of them in the
      base Derivative class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      DerivativeVisitor and do whatever is required in its methods.

      @noexport
   */
  virtual void Visit(DerivativeVisitor& visitor) const = 0;

  /**
      @internal
      @brief Allows visitor to modify "this" derivative.

      @noexport
   */
  virtual void Visit(DerivativeModifyingVisitor& visitor);

  /**
      @internal
      @brief Dumps all data of this instrument in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c for all instruments involved at once but can also be
      called "manually" if needed.

      Note that this method does @b not dump the contents of the associated
      SessionData because when we dump several instruments using the same 
      session data we don't have to duplicate the common session data many
      times.

      @param tagParent the parent tag under which our tag(s) should be created
      @return the tag we created (so that caller could add more stuff to it
              before closing it)

      @noexport
   */
  virtual XML::Tag Dump(XML::Tag& tagParent) const = 0;


protected:

  /// Checks if the session data is set. Throw exception if not.
  void CheckSessionData() const;

  /**
      Check if a derivative is a cross currency with the given session data.

      @param sessionData The session data of the derivative
      @return true if the derivative is cross currency with the session data
              false otherwise
   */
  bool IsCrossCurrency(const SessionData& sessionData) const;

  /**
      Verifies if the market price is available. Actually it verifies
      whether HasMarketPrice() returns true and throws if it returns false.
   */
  void CheckMarketPrice() const;
  
  /// Sets the value of m_dPrice.
  void DoSetMarketPrice(double dPrice);

  /**
      Gets the market price.

      A specific derivative needs to overload this function if its market price
      can be implicitly defined by other methods. This is the case for option,
      whose price can be given as implied volatility.

      This function is called by GetMarketPrice() and CheckMarketPrice() has
      already been called. Therefore, the function should neither throw nor
      return invalid value.
   */
  virtual double DoGetMarketPrice() const
  {
    return m_dPrice;
  }

  /// Sets m_dPrice to an invalid value, make it unavailable
  void UnsetMarketPrice();

  /// Returns true if m_dPrice is correctly set
  bool MarketPriceIsSet() const;

  /**
      Dumps entries common to all derivative.

      @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;

  /**
      Dumps value of m_dPrice if it is valid. 
      
      @param tagParent the parent tag under which our tag(s) should be created
   */
  void DumpMarketPrice(XML::Tag& tagParent) const;

  /// the market price or INVALIDPRICE if not set
  double m_dPrice;

  /// Shared pointer to the session data object
  shared_ptr<SessionData> m_pSessionData;

  /// The currency of the derivative
  shared_ptr<Numeraire> m_pNumeraire;

  /// The issuer of the derivative, might be null
  shared_ptr<Issuer> m_pIssuer;

  /// The flags determining which Greeks to compute, how much data to store, etc
  shared_ptr<ComputationalFlags> m_pFlags;

}; // class Derivative


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVE_H_
