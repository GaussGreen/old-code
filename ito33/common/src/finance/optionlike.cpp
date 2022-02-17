/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/optionlike.cpp
// Purpose:     Implementation of shared ptr for Option class
// Author:      ZHANG Yunzhi
// Created:     April 28, 2003
// RCS-ID:      $Id: optionlike.cpp,v 1.7 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/optionlike.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/optionerror.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/optiontype.h"
#include "ito33/xml/finance/exercisetype.h"

extern const ito33::finance::OptionError
  ITO33_OPTION_INVALID_MATURITY_DATE,
  ITO33_OPTION_INVALID_OPTIONTYPE,
  ITO33_OPTION_INVALID_EXERCISETYPE;

namespace ito33
{

namespace finance
{


OptionLike::OptionLike(Date maturityDate,
                       OptionType optionType,
                       ExerciseType exerciseType)
                       :Derivative()
{

  m_maturityDate = maturityDate;
  m_optionType   = optionType;
  m_exerciseType = exerciseType;

  CHECK_COND( m_maturityDate.IsValid(), ITO33_OPTION_INVALID_MATURITY_DATE);

  CHECK_COND( IsValidExerciseType(m_exerciseType), 
              ITO33_OPTION_INVALID_EXERCISETYPE );

  CHECK_COND( IsValidOptionType(m_optionType) , 
              ITO33_OPTION_INVALID_OPTIONTYPE );
}

void OptionLike::DumpMe(XML::Tag& tagParent) const
{
  Derivative::DumpMe(tagParent);

  tagParent.Element(XML_TAG_OPTION_TYPE)
                     (
                      GetNameFromEnumValue(
                        m_optionType,
                        SIZEOF(g_optionTypes),
                        g_optionTypes)
                     );

  tagParent.Element(XML_TAG_OPTION_EXERCISE_TYPE)
               (
                 GetNameFromEnumValue(
                        m_exerciseType,
                        SIZEOF(g_exerciseTypes),
                        g_exerciseTypes)
                );

  tagParent.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);
 
  DumpMarketPrice(tagParent); 
}


} // namespace finance

} // namespace ito33
