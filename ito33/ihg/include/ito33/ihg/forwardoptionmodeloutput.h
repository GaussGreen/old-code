/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/forwardoptionmodeloutput.h
// Purpose:     model output class for option using forward PDE
// Author:      Wang
// Created:     2004/03/10
// RCS-ID:      $Id: forwardoptionmodeloutput.h,v 1.3 2006/02/28 17:33:21 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/forwardoptionmodeloutput.h
    @brief model output class for option using forward PDE
 */

#ifndef _ITO33_IHG_FORWARDOPTIONMODELOUTPUT_H_
#define _ITO33_IHG_FORWARDOPTIONMODELOUTPUT_H_

#include "ito33/finance/modeloutput.h"

namespace ito33
{

namespace ihg
{

class ForwardOptionModelOutput : public finance::ModelOutput
{
public:

  ForwardOptionModelOutput() : finance::ModelOutput() { }

  // Default dtor is ok


private:

  NO_COPY_CLASS(ForwardOptionModelOutput);

};  // class ForwardOptionModelOutput 


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_FORWARDOPTIONMODELOUTPUT_H_

