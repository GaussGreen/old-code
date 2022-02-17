/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/issuer.cpp
// Purpose:     do the necessary for issuer class
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: issuer.cpp,v 1.14 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/date.h"
#include "ito33/arraycheckers.h"

#include "ito33/finance/error.h"
#include "ito33/finance/issuer.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"
#include "ito33/xml/finance/issuer.h"

extern const ito33::finance::Error ITO33_INVALID_FISCAL_YEAR,
                          ITO33_ISSUER_INVALID_DEFAULT_INTENSITY,
                          ITO33_ISSUER_NO_DEFAULT_INTENSITY;

namespace ito33
{

namespace finance
{
  
void Issuer::SetFiscalYearStartDate(Date fiscalYearStartDate)
{ 
  // make sure the date is valid before setting
  CHECK_COND( fiscalYearStartDate.IsValid(), ITO33_INVALID_FISCAL_YEAR );

  m_fiscalYearStartDate = fiscalYearStartDate; 
} 

void Issuer::SetDefaultIntensity
             ( const std::vector<Date>& pDates,
               const std::vector<double>& pdDefaultIntensities )
{
  CHECK_COND(   !pdDefaultIntensities.empty() 
             && pDates.size() == pdDefaultIntensities.size(), 
             ITO33_ISSUER_INVALID_DEFAULT_INTENSITY);

  CheckNonNegativity( pdDefaultIntensities, 
                      ITO33_ISSUER_INVALID_DEFAULT_INTENSITY.GetMessage() );

  m_pDates = pDates;
  m_pdDefaultIntensities = pdDefaultIntensities;
}

bool Issuer::HasDefaultIntensity() const
{
  return !m_pdDefaultIntensities.empty();
}

const std::vector<Date>& Issuer::GetDatesOfDefaultIntensity() const
{
  CHECK_COND( HasDefaultIntensity(), ITO33_ISSUER_NO_DEFAULT_INTENSITY);

  return m_pDates;
}

const std::vector<double>& Issuer::GetValuesOfDefaultIntensity() const
{
  CHECK_COND( HasDefaultIntensity(), ITO33_ISSUER_NO_DEFAULT_INTENSITY);

  return m_pdDefaultIntensities;
}

void Issuer::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagIssuer(XML_TAG_ISSUER_ROOT, tagParent);

  tagIssuer.Element(XML_TAG_ISSUER_FISCALYEARSTART_DATE)
                   (m_fiscalYearStartDate);

  if ( HasDefaultIntensity() )
  {
    XML::Tag tag(XML_TAG_ISSUER_DEFAULTINTENSITY, tagIssuer);
 
    XML::DumpVector( tag, m_pDates,
                     XML_TAG_FINANCE_DATES, XML_TAG_FINANCE_DATE);

    XML::DumpVector( tag, m_pdDefaultIntensities,
                     XML_TAG_FINANCE_VALUES, XML_TAG_FINANCE_VALUE);
  }
}

} // namespace finance

} // namespace ito33
