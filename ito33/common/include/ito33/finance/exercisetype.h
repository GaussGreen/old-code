/******************************************************************************
 * File.........: ito33/finance/exercisetype.h
 * Purpose......: Enum of exercise type for option contract
 * Author.......: Wang
 * Created......: 2004/03/02
 * RCS-ID.......: $Id: exercisetype.h,v 1.6 2006/03/23 09:27:14 yann Exp $
 * Copyright....: (c) 2003 Trilemma LLP
 ******************************************************************************/

/**
   @file ito33/finance/exercisetype.h
   @brief Enumerations of the exercise type of option contract.
 */

#ifndef _ITO33_FINANCE_EXERCISETYPE_H_
#define _ITO33_FINANCE_EXERCISETYPE_H_

namespace ito33
{

namespace finance
{

/// Exercise type of the option: American or European
enum ExerciseType
{
  /// European exercise type
  ExerciseType_European = 0,

  /// American exercise type
  ExerciseType_American = 1
  
  #ifndef __CPP2ANY__
  ,
  /// noexport
  ExerciseType_Max
  #endif

};

/**
   Checks if given exercise type is valid

   @param exerciseType given value
   @return true if exercise type is valid , false if not 
   @noexport
 */
inline bool IsValidExerciseType(ExerciseType exerciseType)
{
  return exerciseType == ExerciseType_European ||
         exerciseType == ExerciseType_American;
  
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXERCISETYPE_H_

