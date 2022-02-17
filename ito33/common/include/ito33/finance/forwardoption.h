/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/forwardoption.h
// Purpose:     financial forward option class
// Author:      ITO33
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoption.h,v 1.9 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/forwardoption.h
    @brief declaration of the financial forward option class
 */

#ifndef _ITO33_FINANCE_FORWARDOPTION_H_
#define _ITO33_FINANCE_FORWARDOPTION_H_

#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL Option;
class ITO33_DLLDECL SessionData;

/**
    ForwardOption can be used internally to price at the same time a list of
    European options.
 */
class ForwardOption
{
public:

  struct Element
  {
    Element(const shared_ptr<Option>& pOption, double dWeight)
          : m_pOption(pOption), m_dWeight(dWeight) { }

    shared_ptr<Option> m_pOption;
    double m_dWeight;
  };

  typedef std::list<Element> Elements;

  /// Default ctor
  ForwardOption() { }

  /**
      ctor from a list of options.

      @param options The European options
   */
  ForwardOption(const std::list< shared_ptr<Option> >& options)
  {
    std::list< shared_ptr<Option> >::const_iterator iter;

    for (iter = options.begin(); iter != options.end(); iter++)
      Add(*iter);
  }

  // Default ctor and dtor are ok.


  /**
      Adds an European option.

      @param pOption The European option to be added
      @param dWeight The weight of the option in the calibration
   */
  void Add(const shared_ptr<Option>& pOption, double dWeight = 1.);

  /**
      @name Methods for accessing the ForwardOption.
   */
  //@{

  /**
      The session data associated with the options.

      @return The session data associated with the options
   */
  const shared_ptr<SessionData>& GetSessionData() const;

  /**
      Gets the biggest maturity date of the options.

      @return The biggest maturity date of the options.
   */
  Date GetMaturityDate() const;

  /**
      Checks if there are options.

      @return false if there are options, true otherwise
   */
  bool IsEmpty() const { return m_elements.empty(); }

  /**
      Gets all the options.

      @return All the European options
   */
  Elements& GetAll() { return m_elements; }

  /**
      Gets all the options, const version.

      @return All the European options
   */
  const Elements& GetAll() const { return m_elements; }

  //@}

private:

  Elements m_elements;

}; // class ForwardOption


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_FORWARDOPTION_H_
