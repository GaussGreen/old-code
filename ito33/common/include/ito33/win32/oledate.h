/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/oledate.h
// Purpose:     Dummy wrapper class for win32 DATE
// Author:      Wang
// Created:     2004/08/04
// RCS-ID:      $Id: oledate.h,v 1.2 2004/10/05 09:13:40 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/win32/oledate.h
    @brief Dummy wrapper class for win32 DATE

    The DATE in win32 is a typedef of double. And the type double is mapped 
    already to VT_R8. We define a dummy classto map it to VT_DATE.
 */

#ifndef _ITO33_WIN32_OLEDATE_H_
#define _ITO33_WIN32_OLEDATE_H_

#include "ito33/win32/winwrap.h"

namespace ito33
{

namespace Win32
{


/**
   Dummy wrapper class for win32 DATE
 */
class OleDate
{
public:
 
  OleDate(DATE date)
  {
    m_date = date;
  }

  operator DATE() { return m_date; }


private:

  DATE m_date;

}; // class OleDate


} // namespace Win32

} // namespace ito33

#endif // _ITO33_WIN32_OLEDATE_H_

