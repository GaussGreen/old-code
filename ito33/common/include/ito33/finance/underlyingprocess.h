/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/underlyingprocess.h
// Purpose:     Base class for underlying process
// Created:     2006/06/01
// RCS-ID:      $Id: underlyingprocess.h,v 1.3 2006/06/24 14:14:12 wang Exp $
// Copyright:   (c) 2004-2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/underlyingprocess.h
    @brief Base class for underlying process
 */

#ifndef _ITO33_FINANCE_UNDERLYINGPROCESS_H_
#define _ITO33_FINANCE_UNDERLYINGPROCESS_H_

#include "ito33/common.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
    Base class for underlying process.

    @rename UnderlyingProcessBase
    @nocreate
    @noexport COM
 */
class ITO33_DLLDECL UnderlyingProcess
{
  
public:

  /// Virtual dtor for any base class
  virtual ~UnderlyingProcess() {}
  
  /**
      The post default volatility. 
      By default it is set to zero.

      @param dPostDefaultVolatility post default volatility
   */
  void SetPostDefaultVolatility(double dPostDefaultVolatility);

  /**
      The post default volatility.

      @return the post default volatility.
   */
  double GetPostDefaultVolatility() const
  {
    return m_dPostDefaultVolatility;
  }


protected:

  /// Ctor
  UnderlyingProcess(); 

  /**
      Writes myself to tag parent.

      @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;


private:  

  /// Post default Volatility 
  double m_dPostDefaultVolatility;
  
}; // class UnderlyingProcess


} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_UNDERLYINGPROCESS_H_
