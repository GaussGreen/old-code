/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/modeloutput.h
// Purpose:     hg base model output class
// Created:     2005/01/13
// RCS-ID:      $Id: modeloutput.h,v 1.10 2006/05/22 10:16:59 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/modeloutput.h

    Base model output class for the homogeneous (HG) model

    Only price is computed for now.
 */

#ifndef _ITO33_HG_MODELOUTPUT_H_
#define _ITO33_HG_MODELOUTPUT_H_

#include "ito33/finance/modeloutput.h"

#include "ito33/hg/dlldecl.h"

namespace ito33
{

namespace hg
{
/**
   Base output class for the homogeneous (HG) model.

   @nocreate
 */
class ModelOutput : public finance::ModelOutput
{
};

} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_MODELOUTPUT_H_
