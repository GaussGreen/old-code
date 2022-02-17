/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams.cpp
// Purpose:     do necessary for numerical parameters class
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: numparams.cpp,v 1.8 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/debug.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/numparams_reference.h"

#include "ito33/finance/qualitycontrol.h"

namespace ito33
{
  using namespace finance;

  namespace numeric
  {
    NumParams::NumParams(const QualityControl & control, double dYears)
      // get default values from NumParamsReference
      : m_dMaxDeltaLogS(NumParamsReference::GetMaxDeltaLogS()),
        m_SchemeType(NumParamsReference::GetSchemeType())
    {
      // please check numparams_reference.h for detail of the implementation.
      double dNbTimeSteps;

      double dNb = NumParamsReference::GetNbTimeStepsFor5Years();
      if(dYears <= 1)
        dNbTimeSteps = dNb / 5.;
      else if(dYears <= 5)
        dNbTimeSteps = dNb / 5. * dYears;
      else if(dYears <= 15)
        dNbTimeSteps = dNb;
      else
        dNbTimeSteps = dYears * dNb / 15;

      switch (control.GetComputationQuality())
      {
      case ComputationQuality_Standard:
        m_nNbSpaceSteps = NumParamsReference::GetMinNbSpaceSteps();
        break;
      case ComputationQuality_High: // double the precision
        dNbTimeSteps += dNbTimeSteps;
        m_nNbSpaceSteps = NumParamsReference::GetMinNbSpaceSteps() * 2;
        break;
      case ComputationQuality_VeryHigh_But_MuchSlower: 
        // 4 times the precision
        dNbTimeSteps *= 4;
        m_nNbSpaceSteps = NumParamsReference::GetMinNbSpaceSteps() * 4;
        break;
        
      default:
        ASSERT_MSG(false, "This ComputationQuality hasn't been implemented.");
      }

      // set Number of time steps and Number of space steps
      m_nNbTimeSteps = (int)dNbTimeSteps;
      if(m_nNbSpaceSteps < m_nNbTimeSteps)
        m_nNbSpaceSteps =  m_nNbTimeSteps;

    }
  }
}

