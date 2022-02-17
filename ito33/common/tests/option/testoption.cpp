/////////////////////////////////////////////////////////////////////////////
// Name:        testoption.cpp
// Purpose:     Acceptance test for option 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testoption.cpp,v 1.4 2006/07/28 21:02:49 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/sharedptr.h"

#include "ito33/finance/option.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/exercisetype.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testoption.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Conversion schedule reset
// ----------------------------------------------------------------------------

void OptionTest::NegativeStrike()
{
  double dStrike = -1;
  Date maturityDate(2004,Date::Jun,1);
  ExerciseType exercisetype = ExerciseType_American;
  OptionType optiontype     = Option_Call;


  Option opt(dStrike,maturityDate,optiontype,exercisetype); 
                                 
}//OptionTest::NegativeStrike()

  
void OptionTest::Dump()
{
  double dStrike = 100;
  Date maturityDate(2004,Date::Jun,1);
  ExerciseType exercisetype = ExerciseType_American;
  OptionType optiontype     = Option_Call;


  Option opt(dStrike,maturityDate,optiontype,exercisetype); 

  std::ostringstream oss;

    ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<option>\n"
    "<type>call</type>\n"
    "<exercise>american</exercise>\n"
    "<maturity>2004-06-01</maturity>\n" 
    "<strike>100</strike>\n"
    "</option>\n"
    "</root>\n");

  ito33::XML::RootTag root("root",oss);

  opt.Dump(root);

}
