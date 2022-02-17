/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/modeloutput.h
// Purpose:     model output class for regular options
// Author:      Based on ICARE version
// Created:     2003/12/10
// RCS-ID:      $Id: modeloutput.h,v 1.34 2006/05/22 10:15:43 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/modeloutput.h

    Model output class of simple instruments using IHG model.
 */

#ifndef _ITO33_IHG_MODELOUTPUT_H_
#define _ITO33_IHG_MODELOUTPUT_H_

#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/dlldecl.h"

namespace ito33
{

namespace ihg
{

/**
   Base output class for the inhomogeneous (IHG) model

   @nocreate
 */
class ModelOutput : public finance::ModelOutput
{
};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_MODELOUTPUT_H_
