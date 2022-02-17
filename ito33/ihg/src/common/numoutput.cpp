/////////////////////////////////////////////////////////////////////////////
// Name:        common/numoutput.cpp
// Purpose:     implementation of NumOutput class 
// Author:      October 5, 2005
// RCS-ID:      $Id: numoutput.cpp,v 1.4 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"

#include "ihg/numoutput.h"

extern const ito33::finance::Error ITO33_OUTPUT_NOT_AVAILABLE;

// implement the SharedPtrDeleter for NumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::NumOutput);
}

namespace ito33
{

namespace ihg
{

double NumOutput::GetImpliedBrownianVol()
{
  if ( !m_bHasImpliedBrownianVol ) 
       throw EXCEPTION_MSG
         (
           ITO33_OUTPUT_NOT_AVAILABLE,
           TRANS("Brownian implied volatility not available")
         );

  return m_dImpliedBrownianVol;
}

} //end namespace ihg

} //end namepace ito33
