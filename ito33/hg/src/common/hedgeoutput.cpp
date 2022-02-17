/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/hedgeoutput.cpp
// Purpose:     implementation of HG hedge output class 
// RCS-ID:      $Id: hedgeoutput.cpp,v 1.6 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/error.h"

#include "ito33/hg/hedgeoutput.h"

extern const ito33::finance::Error ITO33_OUTPUT_NOT_AVAILABLE;

using ito33::hg::HedgeOutput;

/* static */
void HedgeOutput::ThrowHERONotAvailable()
{
   throw EXCEPTION_MSG
         (
           ITO33_OUTPUT_NOT_AVAILABLE,
           TRANS("HERO is not available!")
         );
}
