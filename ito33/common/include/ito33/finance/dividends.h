/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/dividends.h
// Purpose:     Dividends is a collection of Dividend objects
// Author:      Vadim Zeitlin
// Created:     03.06.03
// RCS-ID:      $Id: dividends.h,v 1.25 2006/06/18 20:26:07 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/dividends.h
    @brief Declaration of Dividend and Dividends classes.

    Dividend is simply an amount, expressed as either the percentage of the
    current price or as an absolute cash amount, associated with a date. 
    Dividends is just a collection of all dividends associated with something.
 */

#ifndef _ITO33_FINANCE_DIVIDENDS_H_
#define _ITO33_FINANCE_DIVIDENDS_H_

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/dlldecl.h"


namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
    An object of this class represents a single dividend.

    The only subtlety is that dividends may be expressed either in percents of
    current stock price (yield) or in absolute value (cash). Take care to 
    specify the type correctly when creating them.

    @noexport
 */
class ITO33_DLLDECL Dividend
{
public:
  /**
      Type of dividend.

      Note that the values of the first two enum elements are chosen to be
      identical to what FreeBound2kSetDividende() uses so don't change them!
   */
  enum Type
  {
    /// value contains the absolute cash amount of the dividend,
    Cash,

    /// value is the amount in percents of the current stock price
    Yield,

    /**
        Pseudo cash dividend.

        A dividend of this type doesn't really exist but is used
        internally to make sure that the dividend paid never exceeds the 
        stock price. It is a mix of the cash and yield cases: the value
        of the dividend is a percentage of the current stock price until 
        some cut off value starting from which the dividend is constant.

        We need two parameters to describe the dividends of this type:
        value is used as in the Cash case and is the dividend value for stock 
        prices larger than the cut off. For smaller prices, the dividend is 
        given by a percentage equal to pseudoYield with 0 < pseudoYield < 1. 
        The cut off stock price is then value/pseudoYield.
     */
    PseudoCash
  };

  /**
      Constructor: as the object is immutable, it must be used to initialize
      all members.

      @param type_ the type of the dividend, element of Type enum above
      @param date_ the date when the dividend is to be paid
      @param value_ the value of the dividend, as absolute cash amount or in
                    percents depending on the value of type_
      @param theta_ the pseudoYield parameter for pseudo cash case, unused by
                    default
   */
  Dividend(Type type_, Date date_, double value_, double theta_ = 0.)
    : type(type_), date(date_), value(value_), pseudoYield(theta_)
  {
  }

  /// the type of the dividend; immutable
  const Type type;

  /// the date the dividend is to be paid on; immutable
  const Date date;

  /// the value of the dividend in type-dependent units; immutable
  const double value;

  /// the parameter of the pseudo cash dividends, unused for the other ones
  const double pseudoYield;


  /// we have to define operator<() to allow Dividends sorting
  bool operator<(const Dividend& div) const
  {
   return date < div.date;
  }

  /**
      Dump this dividend in XML format.

      @param tagParent parent tag
      @return tag containing the description of this dividend
   */
  XML::Tag Dump(XML::Tag& tagParent) const;

private:
  // no assignment operator (the standard one can't be generated and we
  // want to suppress the warning about being unable to generate it)
  Dividend& operator=(const Dividend&);
};

/**
    Dividends is a simple collection of Dividend objects.

    This class provides special functions for adding dividends,but for accessing
    them you can simply use the container directly. Do @b not rely on its exact
    type because it may change later. Always use Elements typedef.
 */
class ITO33_DLLDECL Dividends
{
public:
  /// the type of the container we use for storing the individual Dividend's
  typedef std::list<Dividend> Elements;

  /// Default constructor
  Dividends()
  {
    m_bValidated = true;
  }

  // default dtor, copy ctor and assignment operators are ok

  /**
      Adds a dividend to the list.

      @param div the Dividend object.

      @noexport
   */
  void Add(const Dividend& div);

  /**
      Adds a dividend to the list.

      @param type the kind of the dividend (absolute or percents)
      @param date the dividend date
      @param value the value of the dividend, either in dollars or percents
      @param pseudoYield the pseudoYield parameter for pseudo cash case, unused
                         otherwise

      @noexport
   */
  void 
  Add(Dividend::Type type, Date date, double value, double pseudoYield = 0.)
  {
    Add(Dividend(type, date, value, pseudoYield));
  }

  /**
      Adds a new dividend with value expressed as absolute amount.

      @param date the dividend date
      @param value its value in dollars
   */
  void AddCash(Date date, double value)
  {
    Add(Dividend::Cash, date, value);
  }

  /**
      Adds a new dividend with value expressed in percents of current stock 
      price.

      @param date the dividend date
      @param percent its value in percents of current stock price
   */
  void AddYield(Date date, double percent)
  {
    Add(Dividend::Yield, date, percent);
  }

  /**
      Adds a new pseudo cash dividend.

      @param date the dividend date
      @param value the value of dividend in dollars after cut off
      @param pseudoYield the slope of the dividend curve before the cut off
                         (must be > 0 and <= 1)

      @noexport
   */
  void AddPseudoCash(Date date, double value, double pseudoYield)
  {
    Add(Dividend::PseudoCash, date, value, pseudoYield);
  }

  /// Remove all dividends we store.
  void Clear()
  {
    m_elements.clear();

    m_bValidated = true;
  }

  /**
      This method must be called after adding all the dividends.

      It checks that all the dividends which were added were valid and sorts
      them by date. It does nothing if are already validated.

      If there are any invalid dividends, an exception is thrown.

      Note that it is not necessary to call this method, it is called
      automatically by the functions accessing the dividends if needed. It
      may be useful to call it "manually" if you want to report a possible
      error with the dividends as soon as possible.
   */
  void Validate() const
  {
    if ( !m_bValidated )
    {
      const_cast<Dividends *>(this)->DoValidate();
    }
  }

#ifndef __CPP2ANY__

  /**
      Sets all dividends at once.

      This is equivalent to calling Clear() and then AddAbs() in a loop.

      @sa SetAbs
      @sa SetYield

      @param type the type of the dividend (cash or yield)
      @param pDates pointer to the first date
      @param pEndDates pointer one beyond the last date
      @param pValues pointer to the first value
   */
  template <typename T, typename U>
  void Set(Dividend::Type type, T pDates, T pEndDates, U pValues)
  {
    for ( Clear(); pDates != pEndDates; ++pDates, ++pValues )
    {
      Add(type, *pDates, *pValues);
    }

    Validate();
  }

  /**
      Sets all dividends at once assuming they are all of absolute kind.

      This is equivalent to calling Clear() and then AddAbs() in a loop.

      The parameters may be either the raw pointers or iterators.

      @param pDates pointer to the first date
      @param pEndDates pointer one beyond the last date
      @param pValues pointer to the first value
   */
  template <typename T, typename U>
  void SetCash(T pDates, T pEndDates, U pValues)
  {
    Set(Dividend::Cash, pDates, pEndDates, pValues);
  }

  /**
      Sets all dividends at once assuming they are all of percents kind.

      This is exactly equivalent to SetCash() except that the values are
      expressed in percents of the current stock price.

      @param pDates pointer to the first date
      @param pEndDates pointer one beyond the last date
      @param pYields pointer to the first value (in percents)
   */
  template <typename T, typename U>
  void SetYield(T pDates, T pEndDates, U pYields)
  {
    Set(Dividend::Yield, pDates, pEndDates, pYields);
  }

#endif // #ifndef __CPP2ANY__

  /**
      Sets all dividends at once assuming they are all of absolute kind.

      This is equivalent to calling Clear() and then AddAbs() in a loop.

      The parameters may be either the raw pointers or iterators.

      @param pDates dividend dates
      @param pdValues dividend values as an absolute amount
   */
  void SetCash(const std::vector<Date>& pDates,
               const std::vector<double>& pdValues)
  {
    SetCash( pDates.begin(), pDates.end(), pdValues.begin() );
  }

  /**
      Sets all dividends at once assuming they are all of percents kind.

      This is exactly equivalent to SetCash() except that the values are
      expressed in percents of the current stock price.

      @param pDates dividend dates
      @param pdYields dividend values as percentage of current stock price 
   */
  void SetYield(const std::vector<Date>& pDates,
                const std::vector<double>& pdYields)
  {
    SetYield( pDates.begin(), pDates.end(), pdYields.begin() );
  }

  /**
      Check if there is cash dividends between two given dates(inclusive).

      Cash dividend can render the problem inhomegeneous.

      @param dateBegin The first date
      @param dateEnd The second date, must not be before the first date

      @noexport
   */
  bool HasCashBetween(Date dateBegin, Date dateEnd);

  /**
      Get the container containing all dividends.

      This function throws if we haven't been validated.

      @noexport
   */
  const Elements& GetAll() const { Validate(); return m_elements; }

  /**
      Get the container containing all dividends (non const version).

      This function throws if we haven't been validated.

      @noexport
   */
  Elements& GetAll() { Validate(); return m_elements; }

  /**
      @internal
      @brief Dumps all dividends in this collection in XML format.

      We should have some dividends, otherwise we're going to create an empty
      tag which will be strange (even though not catastrophic).

      @param tagParent parent tag

      @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:
  /// Do unconditionally validate our elements
  void DoValidate();

  /// all the dividends we store
  Elements m_elements;

  /// true if Validate() has been called and Add() hasn't been called since
  bool m_bValidated;
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DIVIDENDS_H_

