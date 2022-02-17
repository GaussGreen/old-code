/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams_reference.cpp
// Purpose:     Definition of static functions of NumParamsReference
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: numparams_reference.cpp,v 1.6 2006/05/26 13:34:31 nabil Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/numparams_reference.h"

#include "ito33/numeric/schemetype.h"

namespace ito33
{

namespace numeric
{

namespace NumParamsReference
{

size_t g_nNbTimeStepsFor5Years = 150;

size_t g_nMinNbSpaceSteps = 100;

double g_dMaxDeltaLogS = 0.15;

SchemeType g_iSchemeType = SchemeType_ThreeLevel;

size_t GetMinNbSpaceSteps()
{
  return g_nMinNbSpaceSteps;
}

size_t GetNbTimeStepsFor5Years()
{
  return g_nNbTimeStepsFor5Years;
}

double GetMaxDeltaLogS()
{
  return g_dMaxDeltaLogS;
}


SchemeType GetSchemeType()
{
  return g_iSchemeType; 
}

}

}

}


