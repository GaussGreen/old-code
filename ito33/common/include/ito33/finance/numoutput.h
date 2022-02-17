/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/numoutput.h
// Purpose:     finance level numerical output
// Created:     September 27, 2005
// RCS-ID:      $Id: numoutput.h,v 1.2 2006/02/28 16:13:48 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/numoutput.h
    @brief base numerical output class
 */

#ifndef _ITO33_FINANCE_NUMOUTPUT_H_
#define _ITO33_FINANCE_NUMOUTPUT_H_

namespace ito33
{

namespace finance
{

/// Base output class for all numerical output.
class NumOutput
{
public:
  
  NumOutput() { }

  /// Dummy virtual dtor for base class
  virtual ~NumOutput() { }

}; // class NumOutput

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_NUMOUTPUT_H_
