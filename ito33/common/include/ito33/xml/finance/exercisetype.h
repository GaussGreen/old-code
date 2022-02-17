/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/exercisetype.h
// Purpose:     Names of elements and attributes used in XML for exercise type.
// Created:     November 7, 2005
// RCS-ID:      $Id: exercisetype.h,v 1.3 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/exercisetype.h
    @brief 
 */

#ifndef _ITO33_XML_FINANCE_EXERCISE_TYPE_H_
#define _ITO33_XML_FINANCE_EXERCISE_TYPE_H_

#include "ito33/finance/exercisetype.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_OPTION_EXERCISE_TYPE                     "exercise"

#define XML_VALUE_OPTION_EXERCISE_TYPE_EUROPEAN          "european"

#define XML_VALUE_OPTION_EXERCISE_TYPE_AMERICAN          "american"

//@}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<finance::ExerciseType> g_exerciseTypes[] =
  {
    { XML_VALUE_OPTION_EXERCISE_TYPE_EUROPEAN, finance::ExerciseType_European},
    { XML_VALUE_OPTION_EXERCISE_TYPE_AMERICAN, finance::ExerciseType_American},
  };

}

#endif // #define _ITO33_XML_FINANCE_EXERCISE_TYPE_H_
