/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/specialtimes.h
// Purpose:     special times help for time mesh construction
// Author:      Wang
// Created:     2004/10/08
// RCS-ID:      $Id: specialtimes.h,v 1.2 2004/11/08 12:13:55 wang Exp $
// Copyright:   (c) 2004  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/numeric/mesh/specialtimes.h
   @brief special times helps for time mesh construction

   We are not doing the same thing for space mesh because space mesh is 
   constructed very differently from time mesh.
 */

#ifndef _ITO33_NUMERIC_MESH_SPECIALTIMES_H_
#define _ITO33_NUMERIC_MESH_SPECIALTIMES_H_

#include "ito33/list.h"

#include "ito33/numeric/predicatetime.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{


/// Refine level for a given time point
enum RefineLevel
{
  /// Standard refine level
  RefineLevel_Standard,

  // High refine level
  RefineLevel_High,

  // Very High refine level
  RefineLevel_VeryHigh,

  RefineLevel_Max

}; // enum RefineLevel


/// special time point having a refine level
class SpecialTime
{
  
public:

  SpecialTime(double dTime, RefineLevel refineLevel = RefineLevel_Standard)
    :m_dTime(dTime), m_refineLevel(refineLevel) { }

  double GetTime() const { return m_dTime; }

  RefineLevel GetRefineLevel() const { return m_refineLevel; }


private:

  /// the time 
  double m_dTime;

  /// the refine level at this time point
  RefineLevel m_refineLevel;

}; // class SpecialTime


inline bool operator < (const SpecialTime& time1, const SpecialTime& time2)
{
  return IsBefore( time1.GetTime(), time2.GetTime() ); 
}

inline bool operator > (const SpecialTime& time1, const SpecialTime& time2)
{
  return IsAfter( time1.GetTime(), time2.GetTime() ); 
}

class SpecialTimes : public std::list<SpecialTime>
{
};


} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_SPECIALTIMES_H_
