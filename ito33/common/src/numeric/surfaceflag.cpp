/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/surfaceflag.cpp
// Purpose:     implementation of the flag value surface class
// Created:     2004/05/06
// RCS-ID:      $Id: surfaceflag.cpp,v 1.8 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 - 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/numeric/surfaceflag.cpp
   @brief implementation of the flag value surface class

   @todo For the moment, the smaller value is chosen for a point between 
         two spots, which may need to be changed.
 */

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/dateutils.h"
#include "ito33/debugparameters.h"

#include "ito33/numeric/domain.h"
#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/surfaceflag.h" 

namespace ito33
{

namespace numeric
{


void SurfaceFlag::GetValuesAt
                  (size_t nIdxT, const Spots& pdS, Flags& values) const
{
  const Flags& flagsTmp = GetValuesAt(nIdxT);
  const Spots& spotsTmp = m_pDomain->GetOutputSpaceMeshAt(nIdxT);
  size_t nNbSTmp = flagsTmp.size();

  values.clear();
  values.reserve( pdS.size() );
  
  size_t nIdxSTmp = 0;

  for (size_t nIdxS = 0; nIdxS < pdS.size(); nIdxS++)
  {
    while (nIdxSTmp < nNbSTmp && spotsTmp[nIdxSTmp] < pdS[nIdxS])
      nIdxSTmp++;
    
    if (nIdxSTmp == 0)
      values.push_back(flagsTmp[0]);
    else if ( nIdxSTmp == nNbSTmp )
      values.push_back(flagsTmp[nNbSTmp - 1]);
    else
    {
      Flag flag = (flagsTmp[nIdxS] > flagsTmp[nIdxS - 1]) 
                ? flagsTmp[nIdxS - 1] : flagsTmp[nIdxS];

      values.push_back(flag);
    }
  }
}

void SurfaceFlag::DumpToFiles() const
{
  std::string debugDir(DebugParameters::GetDebugDir());
  std::string suffix = ".txt";
 
  size_t nIdxS;

  SurfaceDataInt::const_iterator iter;
  size_t nIdxT = 0;
  for (iter = m_surface.begin(); iter != m_surface.end(); ++iter, ++nIdxT)
  {
    const SurfaceDataInt::DataAtSameTime& dataAtSameTime = *iter;

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


