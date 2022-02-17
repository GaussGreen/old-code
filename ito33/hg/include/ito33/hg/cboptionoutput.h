/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/cboptionoutput.h
// Purpose:     Output class for cb option contracts
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionoutput.h,v 1.4 2006/05/22 10:16:59 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/cboptionoutput.h
    @brief Output class for cb option contracts.
 */

#ifndef _ITO33_HG_CBOPTIONOUTPUT_H_
#define _ITO33_HG_CBOPTIONOUTPUT_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/cboptionoutput.h"

#include "ito33/hg/dlldecl.h"

namespace ito33
{

namespace hg
{

/**
    ModelOutput-derived class for cb option contracts.

    @nocreate
 */
class CBOptionOutput : public finance::CBOptionOutput
{
};  // class CBOptionOutput

} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_CBOPTIONOUTPUT_H_
