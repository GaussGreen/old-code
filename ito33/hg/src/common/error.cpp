///////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/error.cpp
// Purpose:     HG error codes implementation
// Created:     2005/01/13
// RCS-ID:      $Id: error.cpp,v 1.9 2006/06/12 17:03:49 zhang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

#include "ito33/hg/error.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{

namespace hg
{

  ITO33_IMPLEMENT_ERROR_CLASS;

}

}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------
using ito33::hg::Error;

extern const Error ITO33_HG_LICENSE
("No license for using 'hg' feature.");

extern const Error ITO33_HG_UNDERLYINGPROCESS
("Invalid underlying process.");

extern const Error ITO33_HG_SHARPERATIO
("Invalid Sharpe ratio. The Sharpe ratio must not be negative.");

extern const Error ITO33_HG_Intensity
("Negative jump intensity."); 

extern const Error ITO33_HG_REGIME
("At least one regime must be specified.");

extern const Error ITO33_HG_REGIMEINDEX
("Required regime does not exist.");

extern const Error ITO33_HG_REGIMESIZE
("Dimension of the model data is not the same as the number of regimes.");

extern const Error ITO33_HG_AMPLITUDE
("Jump amplitude cannot be less than -100%.");   

extern const Error ITO33_HG_COMPUTEFUNCTION
("Homogeneous model cannot compute price for this derivative.");

extern const Error ITO33_HG_PATHDEPENDENT
("Path dependent features not yet supported.");

extern const Error ITO33_HG_EXCHANGEABLE
("Exchangeable convertible not yet supported.");

extern const Error ITO33_HG_CALIBRATIONFLAGS
("Not enough flags are specified for calibration parameters.");

//_________________________HERO________________________________________

extern const Error ITO33_HG_HERO_MISSINGDATA
("Missing data: HERO calculations require full surface data for the " 
"target contract, and surface data from the valuation date to the " 
"lesser of the target maturity date and hedge contract maturity date " 
"for all hedge contracts (is there a path-dependent contract?).");

