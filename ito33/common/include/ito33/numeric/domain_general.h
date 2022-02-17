/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/domain_general.h
// Purpose:     Domain class with seperated space mesh for each time
// Author:      ZHANG Yunzhi
// Created:     2004/09/07
// RCS-ID:      $Id: domain_general.h,v 1.3 2005/09/14 14:39:32 nabil Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/domain_general.h
    @brief Numeric Domain class defines a fixed space mesh. 
 */

#ifndef _ITO33_NUMERIC_DOMAIN_GENERAL_H_
#define _ITO33_NUMERIC_DOMAIN_GENERAL_H_

#include <math.h>

#include "ito33/constants.h"
#include "ito33/numeric/domain.h"
#include "ito33/numeric/surfacedata.h"

namespace ito33
{

namespace numeric
{

/**
   Domain class with same space mesh for all times. 
 */
class DomainGeneral : public Domain
{
 
public:

  /**
     Ctor constructs a void domain
   */
  DomainGeneral()
  {
  }

  /**
     Gets the first space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetFirstSpaceMeshAt(size_t nIdxT) const;
 
  /**
     Gets the last space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetLastSpaceMeshAt(size_t nIdxT) const;

  /**
     Gets the output space mesh at the given time index

     @param nIdxT the index of the time at which the space mesh is required

     @return the space mesh at current time
   */
  const Spots& GetOutputSpaceMeshAt(size_t nIdxT) const;

  /**
     Add a time point and spot points correspond to
     "initial values for next step" to the numerical domain 

     @param dTime A time point to be added to the domain
     @param pdSpots spots to be added to the domain
     @param bIsInitValuesForNextStep is true if the given spots is 
          special and is correspond to"initial values for next step"
   */
  void AddSpotsAtTime
            (
              const Domain::Spots& pdSpots,
              double dTime,
              bool bIsInitValuesForNextStep = false
            )
  { 
    numeric::Domain::AddTime(dTime);

    if(bIsInitValuesForNextStep)
      m_spots.AddInitValuesForNextStep(pdSpots);
    else
      m_spots.Add(pdSpots, GetCurrentIndex());
  }

private:

  /// the space mesh points
  SurfaceDataDouble m_spots;
   
}; // class DomainGeneral


inline const Domain::Spots& 
DomainGeneral::GetFirstSpaceMeshAt(size_t nIdxT) const
{
  return m_spots.GetFirstValuesAt(nIdxT);
}

inline const Domain::Spots& 
DomainGeneral::GetLastSpaceMeshAt(size_t nIdxT) const
{
  return m_spots.GetLastValuesAt(nIdxT);
}

inline const Domain::Spots& 
DomainGeneral::GetOutputSpaceMeshAt(size_t nIdxT) const
{
  return m_spots.GetValuesAt(nIdxT);
}


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_DOMAIN_GENERAL_H_
