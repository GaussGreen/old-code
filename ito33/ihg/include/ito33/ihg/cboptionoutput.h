/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/cboptionoutput.h
// Purpose:     Output class for cb option contracts
// Author:      Nabil
// Created:     2005/06/22
// RCS-ID:      $Id: cboptionoutput.h,v 1.4 2006/05/22 10:15:43 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/cboptionoutput.h
    @brief Output class for cb option contracts.
*/

#ifndef _ITO33_IHG_CBOPTIONOUTPUT_H_
#define _ITO33_IHG_CBOPTIONOUTPUT_H_

#include "ito33/finance/bondlike/cboptionoutput.h"
#include "ito33/ihg/dlldecl.h"

namespace ito33
{

namespace ihg
{

/**
   ModelOutput class for cb option contracts.

   @nocreate
 */
class CBOptionOutput : public finance::CBOptionOutput
{
};  // class CBOptionOutput


} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_CBOPTIONOUTPUT_H_
