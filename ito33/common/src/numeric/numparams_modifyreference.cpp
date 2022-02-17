/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams_reference.cpp
// Purpose:     Definition of static functions of NumParamsReference
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: numparams_modifyreference.cpp,v 1.5 2006/05/26 13:34:31 nabil Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/numparams_modifyreference.h"


#ifdef ITO33_TEST_MODENV

namespace ito33
{

namespace numeric
{

enum SchemeType;

namespace NumParamsReference
{

extern size_t g_nNbTimeStepsFor5Years;

extern size_t g_nMinNbSpaceSteps;

extern double g_dMaxDeltaLogS;

extern SchemeType g_iSchemeType;


void SetNbTimeStepsFor5Years(size_t nNumber)
{
  g_nNbTimeStepsFor5Years = nNumber;
}

void SetMinNbSpaceSteps(size_t nNumber)
{
  g_nMinNbSpaceSteps = nNumber;
}

void SetMaxDeltaLogS(double dValue)
{
  g_dMaxDeltaLogS = dValue;
}

void SetSchemeType(SchemeType schemeType)
{ 
  g_iSchemeType = schemeType; 
}

}

}

}

#endif
