/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/numeric/domain.cpp
// Purpose:     implementation of the base domain class
// Author:      Wang
// Created:     2004/05/03
// RCS-ID:      $Id: domain.cpp,v 1.11 2006/08/19 23:10:11 wang Exp $
// Copyright:   (c) 2004 -   Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/numeric/domain.cpp
   @brief implementation of the base domain class
 */

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/debug.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/domain.h"

extern const ito33::Error ITO33_OUTOFBOUND;

using ito33::numeric::Domain;

// the temporary vector datesTmp may not be needed
void Domain::GenerateOutputDates()
{
  size_t nIdxT;

  size_t nNbTimes = m_pdTimes.size();

  // if the time vector is empty, there is nothing to do. Actually, this
  // likely indicates an error in the code
  ASSERT_MSG(nNbTimes > 0, 
             "No time mesh available when generating output domain.");

  std::vector<Date> datesTmp( nNbTimes );

  for (nIdxT = 0; nIdxT < nNbTimes; nIdxT++)
    datesTmp[nIdxT] = ito33::GetDateFrom(m_pdTimes[nIdxT]);

  Date date;
  Date dateOld = datesTmp[0];
 
  m_dates.clear();
  m_pnIdxDates.clear();

  std::unique_copy( datesTmp.begin(), datesTmp.end(), 
                    std::insert_iterator< std::vector<Date> >
                    (
                      m_dates, m_dates.begin()
                    ) 
                  );

  nIdxT = 0;
  for (size_t nIdxDate = 0; nIdxDate < m_dates.size(); nIdxDate++)
  {
    while ( nIdxT < datesTmp.size() && datesTmp[nIdxT] == m_dates[nIdxDate] )
      nIdxT++;

    m_pnIdxDates.push_back(nIdxT - 1);
  }

}

/* static */
void Domain::ThrowInvalidDateIndex()
{
  throw EXCEPTION_MSG
        (
          ITO33_OUTOFBOUND,
          TRANS("Invalid date index")
        );
}

