///////////////////////////////////////////////////////////////////////////////
// Name:        ihg/error.cpp
// Purpose:     implementation of ito33::Error for ihg
// Author:      ZHANG Yunzhi
// Created:     06.30.04
// RCS-ID:      $Id: error.cpp,v 1.20 2006/07/26 14:51:51 nabil Exp $
// Copyright:   (c) 2004- Trilemma LLP
///////////////////////////////////////////////////////////////////////////////


#include "ito33/ihg/error.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{

namespace ihg
{

  ITO33_IMPLEMENT_ERROR_CLASS;

}

}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------
using ito33::ihg::Error;

extern const Error ITO33_IHG_LICENSE
("No license for using 'ihg_base' feature.");

extern const Error ITO33_IHG_LICENSE_CB
("No license for using 'ihg_cb' feature.");

extern const Error ITO33_IHG_NULL_VOLATILITY
("Invalid volatility in inhomogeneous underlying process.");

extern const Error ITO33_IHG_NULL_HAZARDRATE
("Invalid hazard rate in inhomogeneous underlying process.");

extern const Error ITO33_IHG_CANT_COMPUTE
("Inhomogeneous model can't compute price for given derivative.");

extern const Error ITO33_IHG_IMPLIEDVOL
("Cannot find the implied volatility for this derivative: " 
"the given price can't be matched");

extern const Error ITO33_IHG_PERFECT_HEDGE_FAILED
("Perfect hedge is impossible for given instruments. " 
"Simple Delta hedge is suggested in this case.");

extern const Error ITO33_IHG_PERFECT_HEDGE_EXCHANGEABLE
("Perfect hedge is not available for (or with) an exchangeable.");

extern const Error ITO33_IHG_PERFECT_HEDGE_CBOPTION_EXCHANGEABLE
("Perfect hedge is not available for (or with) a CB option with "
"an exchangeable CB.");

extern const Error ITO33_IHG_PERFECT_HEDGE_BAD_CROSSCURRENCY_HEDGER
("Perfect hedge: The cross-currency derivative used to hedge the target " 
"derivative has a currency different from the currency of the target " 
"derivative.");

extern const Error ITO33_IHG_PERFECT_HEDGE_DIFFERENT_UNDERLYING
("Can't hedge with an instrument which has a session " 
"different from the one of the instrument to hedge.");

//----------------------hazard rate----------------------------------------
extern const Error ITO33_IHG_HRTO_CDSCURVE_EMPTYCURVE
("Empty cds curve.");

//---------------------- calibration --------------------------------------
extern const Error ITO33_IHG_CALIBRATION_DEGENERATE
("Calibration failed, with possibly degenerate data. Please try "
"a time-only parametrization.");

extern const Error ITO33_IHG_INCORRECT_DERIVATIVELIST 
("Incorrect derivative list: an even number of derivatives is needed, "
"and 3 (or more) derivatives cannot have the same maturity date.");

