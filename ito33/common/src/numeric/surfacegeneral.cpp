/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/surfacegeneral.cpp
// Purpose:     implementation of the general double value surface class
// Author:      Wang
// Created:     2004/05/07
// RCS-ID:      $Id: surfacegeneral.cpp,v 1.12 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 -   Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/numeric/surfacegeneral.cpp
   @brief implementation of the general double value surface class
 */

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/debugparameters.h"
#include "ito33/dateutils.h"
#include "ito33/constants.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/domain.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/deltagamma.h"

namespace ito33
{

namespace numeric
{


void SurfaceGeneral::GetValuesAt
     (size_t nIdxT, const Spots& pdS, Doubles& values) const
{
  values.clear();

  values.resize( pdS.size() );
  
  const Doubles& valuesTmp = GetValuesAt(nIdxT);

  const Spots& spotsTmp = m_pDomain->GetOutputSpaceMeshAt(nIdxT);

  Interpolate(&spotsTmp[0], &valuesTmp[0], valuesTmp.size(),
              &pdS[0], &values[0], pdS.size(),
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void SurfaceGeneral::GetFirstValuesAt
     (size_t nIdxT, const Spots &pdS, Doubles &values) const
{
  values.clear();

  values.resize( pdS.size() );
  
  const Doubles& valuesTmp = GetFirstValuesAt(nIdxT);

  const Spots& spotsTmp = m_pDomain->GetFirstSpaceMeshAt(nIdxT);

  Interpolate(&spotsTmp[0], &valuesTmp[0], valuesTmp.size(),
              &pdS[0], &values[0], pdS.size(),
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void SurfaceGeneral::GetLastValuesAt
     (size_t nIdxT, const Spots &pdS, Doubles &values) const
{
  values.clear();

  values.resize( pdS.size() );
  
  const Doubles& valuesTmp = GetLastValuesAt(nIdxT);

  const Spots& spotsTmp = m_pDomain->GetLastSpaceMeshAt(nIdxT);

  Interpolate(&spotsTmp[0], &valuesTmp[0], valuesTmp.size(),
              &pdS[0], &values[0], pdS.size(),
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void SurfaceGeneral::GetFirstValuesAt
     (double dTime, const Spots& pdS, Doubles& values) const
{

  // Find the time indices bounding the specified time.
  // Want to err on the side of larger times, since we assume backward 
  // pricing, and want to avoid issues with events. 
  double dTimeTmp = dTime + 0.5 * TIMETOLERANCE;
  const std::vector<double>& pdTimes = m_pDomain->GetTimes();
  size_t nNbTimes = pdTimes.size();

  // Special case: If past maturity, return zero.
  // Backward pricers have maturity at time zero.
  if ( IsAfter(dTime, pdTimes[0]) )
  {
    for (size_t nIdx = 0; nIdx < pdS.size(); nIdx++)
      values[nIdx] = 0.0;

    return;
  }

  // Usual case
  size_t nIdxTime = 1;
  while ( pdTimes[nIdxTime] > dTimeTmp && nIdxTime < nNbTimes - 1 )
    nIdxTime++;

  // Check for exact match
  if ( AreTimesEqual(dTime, pdTimes[nIdxTime-1]) )
  {
    const Doubles& valuesTmp = GetFirstValuesAt(nIdxTime-1);
    const Spots& spotsTmp = m_pDomain->GetFirstSpaceMeshAt(nIdxTime-1);

    Interpolate(&spotsTmp[0], &valuesTmp[0], valuesTmp.size(),
                &pdS[0], &values[0], pdS.size(),
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);

    return;
  }

  if ( AreTimesEqual(dTime, pdTimes[nIdxTime]) )
  {
    const Doubles& valuesTmp = GetFirstValuesAt(nIdxTime);
    const Spots& spotsTmp = m_pDomain->GetFirstSpaceMeshAt(nIdxTime);

    Interpolate(&spotsTmp[0], &valuesTmp[0], valuesTmp.size(),
                &pdS[0], &values[0], pdS.size(),
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);

    return;
  }
   

  // Values at dTime should be between nIdxTime-1 and nIdxTime
  // Do the space interpolation
  const Doubles& valuesTmp1 = GetLastValuesAt(nIdxTime-1);
  const Spots& spotsTmp1 = m_pDomain->GetLastSpaceMeshAt(nIdxTime-1);

  Doubles pdValues1(pdS.size());
  Interpolate(&spotsTmp1[0], &valuesTmp1[0], valuesTmp1.size(),
              &pdS[0], &pdValues1[0], pdS.size(),
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  const Doubles& valuesTmp2 = GetFirstValuesAt(nIdxTime);
  const Spots& spotsTmp2 = m_pDomain->GetFirstSpaceMeshAt(nIdxTime);

  Doubles pdValues2(pdS.size());
  Interpolate(&spotsTmp2[0], &valuesTmp2[0], valuesTmp2.size(),
              &pdS[0], &pdValues2[0], pdS.size(),
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  // Now do time interpolation
  double dTime1 = pdTimes[nIdxTime-1];
  double dTime2 = pdTimes[nIdxTime];
  double dWeight = (dTime - dTime1) / ( dTime2 - dTime1);
  
  for (size_t nIdx = 0; nIdx < pdS.size(); nIdx++)
    values[nIdx] = dWeight * pdValues2[nIdx] + (1.-dWeight) * pdValues1[nIdx];
  
}

void SurfaceGeneral::GetLastValuesAt
     (double /*dTime*/, const Spots& /*pdS*/, Doubles& /*values*/) const
{
    ASSERT_MSG(false, "Not yet implemented");
}

void SurfaceGeneral::GetDeltaAndGamma
                      (
                        shared_ptr<SurfaceDouble>& pDeltaSurfaceIn,
                        shared_ptr<SurfaceDouble>& pGammaSurfaceIn
                      ) const
{
  // The new surfaces replace whatever was passed in.  Assigned at the bottom,
  // to avoid type problems
  shared_ptr<SurfaceGeneral> pDeltaSurface( new SurfaceGeneral( m_pDomain ) );
  shared_ptr<SurfaceGeneral> pGammaSurface( new SurfaceGeneral( m_pDomain ) );

  // vector to store the derivs. 
  Doubles pdDeltaTmp;
  Doubles pdGammaTmp;

  size_t nNbS;
  SurfaceDataDouble::const_iterator iter;
  size_t nIdxTime; // also save the count

  // loop over the times in the domain/surface
  for(iter = m_surface.begin(), nIdxTime = 0;
      iter != m_surface.end();
      iter++, nIdxTime++)
  {
//    const SurfaceDataAtSameTime& currentData = *iter;

    size_t nCurrentDataSize = iter->SizeOfStandardValues();

    // always start by first
    {
      const Spots& spots = m_pDomain->GetFirstSpaceMeshAt(nIdxTime);
      const Doubles& prices = m_surface.GetFirstValuesAt(nIdxTime);

      nNbS = spots.size();
      pdDeltaTmp.resize( nNbS );
      numeric::ComputeDelta(&spots[0], &prices[0], nNbS, &pdDeltaTmp[0]);
      pdGammaTmp.resize( nNbS );
      numeric::ComputeGammaFD(&spots[0], &prices[0], nNbS, &pdGammaTmp[0]);

      pDeltaSurface->Insert(pdDeltaTmp, nIdxTime);
      pGammaSurface->Insert(pdGammaTmp, nIdxTime);
    }

    // work on last standard values if size is 2 here
    if(nCurrentDataSize == 2)
    {
      const Spots& spots = m_pDomain->GetOutputSpaceMeshAt(nIdxTime);
      const Doubles& prices = m_surface.GetValuesAt(nIdxTime);

      nNbS = spots.size();
      pdDeltaTmp.resize( nNbS );
      numeric::ComputeDelta(&spots[0], &prices[0], nNbS, &pdDeltaTmp[0]);
      pdGammaTmp.resize( nNbS );
      numeric::ComputeGammaFD(&spots[0], &prices[0], nNbS, &pdGammaTmp[0]);

      pDeltaSurface->Insert(pdDeltaTmp, nIdxTime);
      pGammaSurface->Insert(pdGammaTmp, nIdxTime);
    }

    // deal with the last if any
    if(iter->HasInitValuesForNextStep())
    {
      const Spots& spots = m_pDomain->GetLastSpaceMeshAt(nIdxTime);
      const Doubles& prices = m_surface.GetLastValuesAt(nIdxTime);

      nNbS = spots.size();
      pdDeltaTmp.resize( nNbS );
      numeric::ComputeDelta(&spots[0], &prices[0], nNbS, &pdDeltaTmp[0]);
      pdGammaTmp.resize( nNbS );
      numeric::ComputeGammaFD(&spots[0], &prices[0], nNbS, &pdGammaTmp[0]);

      pDeltaSurface->InsertInitConditionForNextStep(pdDeltaTmp, nIdxTime);
      pGammaSurface->InsertInitConditionForNextStep(pdGammaTmp, nIdxTime);
    }
  }

  ASSERT_MSG(nIdxTime == m_pDomain->GetTimes().size(), 
          "Bad input, domain size is not coherrent with the surface size.");

  // Set the output. This will possibly delete whatever the user passed in,
  // and the user will get a surface general whether they wanted to or not
  pDeltaSurfaceIn = pDeltaSurface;
  pGammaSurfaceIn = pGammaSurface;

}


void SurfaceGeneral::GetThetaBackwardOnly
                          (shared_ptr<SurfaceDouble>& pThetaSurfaceIn) const
{
  ASSERT_MSG( m_pDomain->GetTimes()[0] > m_pDomain->GetTimes()[1],
              "GetThetaBackwardOnly() is only available for backward problem."
            );

  // The new surface replaces whatever was passed in.  Assigned at the bottom,
  // to avoid type problems
  shared_ptr<SurfaceGeneral> pThetaSurface( new SurfaceGeneral( m_pDomain ) );

  // vector to store the derivs. 
  Doubles pdThetaTmp;

  SurfaceDataDouble::const_iterator iterRight, iterHere;

  iterRight = m_surface.begin();
  iterHere = m_surface.begin();
  iterHere++;

  ASSERT_MSG( iterHere != m_surface.end(),
              "Not enough points to construct theta surface");
  
  size_t nNbTimes = m_pDomain->GetTimes().size();

  // Do the last time by backward differencing. These values will
  // be the same as forward differencing at the previous timestep,
  // assuming the same grid.  If the grid is different, interpolating
  // should be good enough
  {
    // Begin of ComputationOfThetaCurve
    //=========================================================================
    //

    // Get the right values at this time, and the left values at the next,
    // thereby avoiding event effects (ie. there can be no intervening
    // events before this pair of data)
    const Spots&
      spotsTheta      = m_pDomain->GetFirstSpaceMeshAt(1),
      spotsRight      = m_pDomain->GetLastSpaceMeshAt(0);
    const Doubles& 
      valuesLeft      = iterHere->GetFirstValues(),
      valuesRightInit = iterRight->GetLastValues();
    size_t 
      nNbLeft         = spotsTheta.size(),
      nNbRight        = spotsRight.size();

    // The "right" values may don't use the same space mesh, do interpolation 
    Doubles valuesRight(nNbLeft);
    numeric::QuadraticInterpolate
                (
                  &spotsRight[0], &valuesRightInit[0], nNbRight, 
                  &spotsTheta[0], &valuesRight[0], nNbLeft, 
                  numeric::ExtrapolationMode_Linear, 
                  numeric::ExtrapolationMode_Linear
                );
    
    // Get the time difference, which becomes the "epsilon" shift.  Probably
    // a rather large epsilon, but as the time grid is refined the 
    // approximation does get better. Reverse the order of the difference
    // since we want deriv wrt current time, and not the maturity
    double dInvTimeShift = 1.0 /
      (m_pDomain->GetTimeAt(0) - m_pDomain->GetTimeAt(1));

    // Simply do the forward differencing.
    pdThetaTmp.resize(nNbLeft);
    for (size_t nIdxS = 0; nIdxS < nNbLeft; nIdxS++)
    {
      pdThetaTmp[nIdxS] = dInvTimeShift
                        * (valuesRight[nIdxS] - valuesLeft[nIdxS]);
    }
    
    //        now the theta curve is defined by 
    //                  spotsTheta and pdThetaTmp
    //=========================================================================
    // End of ComputationOfThetaCurve

    // The size of output at this time. There must be at least one data entry
    size_t nNbDoubles = iterRight->SizeOfStandardValues();

    // Now check how many data points we have.
    // Since Theta is only for financial use, we should set theta curve in
    // right place where GetValueAt() works. As for other values, we just don't
    // care.
    // note that we ignore the effect of event.
    if (nNbDoubles == 1)
    {
      pThetaSurface->Insert(pdThetaTmp, 0);
    }
    else if (nNbDoubles == 2)
    {
      const Spots& spotsTmp = m_pDomain->GetOutputSpaceMeshAt(0);
      size_t nSizeTmp = spotsTmp.size();
      Doubles pdTmp1(nSizeTmp);
      numeric::QuadraticInterpolate
                  (
                    &spotsTheta[0], &pdThetaTmp[0], spotsTheta.size(), 
                    &spotsTmp[0], &pdTmp1[0], nSizeTmp, 
                    numeric::ExtrapolationMode_Linear, 
                    numeric::ExtrapolationMode_Linear
                  );

      pThetaSurface->Insert(pdThetaTmp, pdTmp1, 0);
    }

    // we just do nothing for theta at end of grid if we are there
    if(iterRight->HasInitValuesForNextStep())
    {
      const Spots& spotsTmp1 = m_pDomain->GetLastSpaceMeshAt(0);
      Doubles pdTmp1(spotsTmp1.size(), 0);

      pThetaSurface->InsertInitConditionForNextStep(pdTmp1, 0);
    }
  }

  for (size_t nIdxTime = 1;
       nIdxTime < nNbTimes;
       nIdxTime++, iterRight++, iterHere++)
  {

    // Begin of ComputationOfThetaCurve
    //=========================================================================
    //

    // Get the right values at this time, and the left values at the next,
    // thereby avoiding event effects (ie. there can be no intervening
    // events before this pair of data)
    const Spots&
      spotsTheta      = m_pDomain->GetFirstSpaceMeshAt(nIdxTime),
      spotsRight      = m_pDomain->GetLastSpaceMeshAt(nIdxTime - 1);
    const Doubles& 
      valuesLeft      = iterHere->GetFirstValues(),
      valuesRightInit = iterRight->GetLastValues();
    size_t 
      nNbLeft         = spotsTheta.size(),
      nNbRight        = spotsRight.size();

    // The "right" values may don't use the same space mesh, do interpolation 
    Doubles valuesRight(nNbLeft);
    numeric::QuadraticInterpolate
                (
                  &spotsRight[0], &valuesRightInit[0], nNbRight, 
                  &spotsTheta[0], &valuesRight[0], nNbLeft, 
                  numeric::ExtrapolationMode_Linear, 
                  numeric::ExtrapolationMode_Linear
                );
    
    // Get the time difference, which becomes the "epsilon" shift.  Probably
    // a rather large epsilon, but as the time grid is refined the 
    // approximation does get better. Reverse the order of the difference
    // since we want deriv wrt current time, and not the maturity
    double dInvTimeShift = 1.0 /
      (m_pDomain->GetTimeAt(nIdxTime - 1) - m_pDomain->GetTimeAt(nIdxTime));

    // Simply do the forward differencing.
    pdThetaTmp.resize(nNbLeft);
    for (size_t nIdxS = 0; nIdxS < nNbLeft; nIdxS++)
    {
      pdThetaTmp[nIdxS] = dInvTimeShift
                        * (valuesRight[nIdxS] - valuesLeft[nIdxS]);
    }
    
    //        now the theta curve is defined by 
    //                  spotsTheta and pdThetaTmp
    //=========================================================================
    // End of ComputationOfThetaCurve

    // The size of output at this time. There must be at least one data entry
    size_t nNbDoubles = iterRight->SizeOfStandardValues();

    // Now check how many data points we have.
    // Since Theta is only for financial use, we should set theta curve in
    // right place where GetValueAt() works. As for other values, we just don't
    // care.
    // note that we ignore the effect of event.
    if (nNbDoubles == 1)
    {
      pThetaSurface->Insert(pdThetaTmp, nIdxTime);
    }
    else if (nNbDoubles == 2)
    {
      const Spots& spotsTmp = m_pDomain->GetOutputSpaceMeshAt(nIdxTime);
      size_t nSizeTmp = spotsTmp.size();
      Doubles pdTmp1(nSizeTmp);
      numeric::QuadraticInterpolate
                  (
                    &spotsTheta[0], &pdThetaTmp[0], spotsTheta.size(), 
                    &spotsTmp[0], &pdTmp1[0], nSizeTmp, 
                    numeric::ExtrapolationMode_Linear, 
                    numeric::ExtrapolationMode_Linear
                  );

      pThetaSurface->Insert(pdThetaTmp, pdTmp1, nIdxTime);
    }

    // we just do nothing for theta at end of grid if we are there
    if(iterRight->HasInitValuesForNextStep())
    {
      const Spots& spotsTmp1 = m_pDomain->GetLastSpaceMeshAt(nIdxTime);
      Doubles pdTmp1(spotsTmp1.size(), 0);

      pThetaSurface->InsertInitConditionForNextStep(pdTmp1, nIdxTime);
    }

  } // loop over the times

  // Set the output. This will possibly delete whatever the user passed in,
  // and the user will get a surface general whether they wanted to or not
  pThetaSurfaceIn = pThetaSurface;

}

shared_ptr<SurfaceDouble> SurfaceGeneral::ComputeFiniteDifference
  (const SurfaceDouble& shiftedSurface, double dInverseShift) const
{
  ASSERT_MSG( dynamic_cast<const SurfaceGeneral*>(&shiftedSurface), 
              "Surfaces do not have the same type.");

  const SurfaceGeneral&
    shiftedGeneralSurface = static_cast<const SurfaceGeneral&>(shiftedSurface);

  // Create the return surface
  shared_ptr<SurfaceGeneral> pNewSurface( new SurfaceGeneral( GetDomain() ) );

  SurfaceDataDouble::const_iterator
    iterSelf = m_surface.begin(),
    iterShifted = shiftedGeneralSurface.m_surface.begin();

  SurfaceDataDouble::iterator
    iterNew = pNewSurface->m_surface.begin();

  size_t nIdxTime;

  size_t nNbTimes = m_pDomain->GetTimes().size();

  for(nIdxTime = 0;
      nIdxTime < nNbTimes;
      nIdxTime++, iterSelf++, iterShifted++, iterNew++)
  {
    ASSERT_MSG
      (
           iterSelf != m_surface.end()
        || iterShifted != shiftedGeneralSurface.m_surface.end()
        || iterNew != pNewSurface->m_surface.end(),
        "Surfaces do not have the same structure."
      );

    iterSelf->ComputeFiniteDifference(*iterShifted, dInverseShift, *iterNew);
    
  } // outermost loop over the times

  return pNewSurface;
}

void SurfaceGeneral::DumpToFiles() const
{
  std::string debugDir(DebugParameters::GetDebugDir());
  std::string suffix = ".txt";
 
  size_t nIdxS;

  SurfaceDataDouble::const_iterator iter;
  size_t nIdxT = 0;
  for (iter = m_surface.begin(); iter != m_surface.end(); ++iter, ++nIdxT)
  {
    const SurfaceDataDouble::DataAtSameTime& dataAtSameTime = *iter;

    double dTime = m_pDomain->GetTimeAt(nIdxT);
    Date date = GetDateFrom(dTime);

    std::string extra;

    // if the time corresponding to an exact date, add the date to file name
    if ( AreTimesEqual(dTime, GetDoubleFrom(date)) )
      extra = "_" + date.Format("%Y_%m_%d");

    // for backward problem, first computed data(last index) first added,
    // so inverse the data index
    std::string fileName = debugDir 
        + String::Printf("%03d", m_pDomain->GetTimes().size() - nIdxT - 1) + extra;
    
    std::ofstream of((fileName + suffix).c_str());
    of.precision(12);
    const Spots& spots = m_pDomain->GetFirstSpaceMeshAt(nIdxT);
    for (nIdxS = 0; nIdxS < spots.size(); nIdxS++)
      of << spots[nIdxS] << '\t' << dataAtSameTime.GetFirstValues()[nIdxS] << '\n';
    of.close();

    if (dataAtSameTime.SizeOfStandardValues() == 2) // has event at current time, 
    {
      std::ofstream of((fileName + "_event" + suffix).c_str());
      of.precision(12);
      const Spots& spots = m_pDomain->GetOutputSpaceMeshAt(nIdxT);
      for (nIdxS = 0; nIdxS < spots.size(); nIdxS++)
        of << spots[nIdxS] << '\t' << dataAtSameTime.GetValues()[nIdxS] << '\n';
      of.close();
    }

    if ( dataAtSameTime.HasInitValuesForNextStep() ) // end of grid
    {
      std::ofstream of((fileName + "_constraint" + suffix).c_str());
      of.precision(12);
      const Spots& spots = m_pDomain->GetLastSpaceMeshAt(nIdxT);
      for (nIdxS = 0; nIdxS < spots.size(); nIdxS++)
        of << spots[nIdxS] << '\t' << dataAtSameTime.GetLastValues()[nIdxS] << '\n';
      of.close();
    }
  }
}


} // namespace numeric

} // namespace ito33

