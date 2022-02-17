/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/numeraire.h
// Purpose:     Declaration of the Numeraire class
// Created:     2004/05/19
// RCS-ID:      $Id: numeraire.h,v 1.17 2006/04/04 16:29:46 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/numeraire.h
    @brief Declaration of Numeraire class.

    Numeraire represents the currency (dollar, euro, ...) of a given money 
    market.

    Use Currency as type name would be better but there has problem with .net. 
 */

#ifndef _ITO33_FINANCE_NUMERAIRE_H_
#define _ITO33_FINANCE_NUMERAIRE_H_

#include "ito33/string.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Numeraire represents the currency (dollar, euro, etc) of a given money 
    market.
 */
class ITO33_DLLDECL Numeraire
{

public:
  /**
      Creates a currency by its code.

      @param strCode code name
   */
  Numeraire(const std::string& strCode) : m_strCode(strCode) {}

  /**
      @internal
      @brief Creates a currency by its code.

      @param strCode code name

      @noexport
   */
  Numeraire(const char* strCode) : m_strCode(strCode) {}

  /**
      Gets code name.

      @return code name
   */
  const std::string& GetCode() const { return m_strCode; }

private:

  std::string m_strCode;
};


/// compares two Numeraires
inline bool operator==(const Numeraire& first, const Numeraire& second)
{
  return first.GetCode() == second.GetCode();
}


/// compares two Numeraires
inline bool operator!=(const Numeraire& first, const Numeraire& second)
{
  return first.GetCode() != second.GetCode();
}


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_NUMERAIRE_H_
