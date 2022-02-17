///////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/optionerror.cpp
// Purpose:     implementation of ito33::Error for option
// Author:      ITO 33 Canada
// Created:     March 30, 2005
// RCS-ID:      $Id: optionerror.cpp,v 1.11 2006/07/21 12:41:43 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////


#include "ito33/finance/optionerror.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{
namespace finance
{
  ITO33_IMPLEMENT_GIVEN_ERROR_CLASS(OptionError);
}
}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------

using ito33::finance::OptionError;

//_____________Common______________________________________________________
extern const OptionError ITO33_OPTION_NEGATIVE_STRIKE
("Option definition: Setting negative strike.");

extern const OptionError ITO33_OPTION_INVALID_MATURITY_DATE
("Option definition: Setting invalid expiry date.");

extern const OptionError ITO33_OPTION_INVALID_EXERCISETYPE
("Option definition: Setting invalid exercise type.");

extern const OptionError ITO33_OPTION_NEGATIVE_IMPLIED_VOLATILITY
("Option definition: Setting invalid implied volatility.");

extern const OptionError ITO33_OPTION_NO_IMPLIED_VOLATILITY
("Option definition: no implied volatility nor market price is set.");

//______________ options____________________________________________
extern const OptionError ITO33_OPTION_INVALID_OPTIONTYPE
("Option definition: Setting invalid option type.");

//_____________Asian options_______________________________________________
extern const OptionError ITO33_ASIANOPTION_NEGATIVE_CURRENTAVERAGE
("Asian option definition: The entered current average is negative.");

extern const OptionError ITO33_ASIANOPTION_AVG_START_DATE_AFTER_MATURITY_DATE
("Asian option definition: average start date is after the maturity date.");

extern const OptionError ITO33_ASIANOPTION_AVG_END_DATE_BEFORE_AVG_START_DATE
("Asian option definition: average end date before average start date.");
  
extern const OptionError ITO33_ASIANOPTION_AVG_END_DATE_AFTER_MATURITY_DATE
("Asian option definition: average end date after maturity date.");

extern const OptionError ITO33_ASIANOPTION_AVG_END_DATE_NOT_SET
("Asian option definition: average end date not set.");

extern const OptionError ITO33_ASIANOPTION_INVALID_NBSAMPLINAVERAGES
("Asian option definition: there must be at least one sampling average.");

extern const OptionError ITO33_ASIANOPTION_INVALID_NB_SAMPLES_USED
("Asian option definition: the number of samples used to compute the " 
"current average must be less than the total number of sampling " 
"averages.");

extern const OptionError ITO33_ASIANOPTION_ZERO_NB_SAMPLES_USED
("Asian option definition: the number of samples used to compute the "
"current average must be at least one.");

//_____________One touches_______________________________________________
extern const OptionError ITO33_ONETOUCH_NEGATIVE_BSBARRIER
("One touch definition: the Black-Scholes barrier value must not be "
"negative.");

extern const OptionError ITO33_ONETOUCH_BARRIER_NOTFOUND
("One touch definition: cannot determine the barrier level from the "
"Black-Scholes barrier (the market quote might be too high).");

extern const OptionError ITO33_ONETOUCH_USE_MARKET_QUOTE
("One touch definition: the market price can only be set by calling "
"SetMarketQuote.");

extern const OptionError ITO33_ONETOUCH_NO_MARKET_QUOTE
("One touch definition: the market quote has not been set.");

extern const OptionError ITO33_ONETOUCH_QUOTE_TOO_LARGE
("One touch definition: problem encountered while setting the market quote "
"(either the market quote or Black-Scholes barrier price is too large).");

