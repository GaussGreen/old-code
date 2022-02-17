/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/impliedcdsspreads.cpp
// Purpose:     Class for computing implied CDS spreads
// Created:     2005/02/08
// RCS-ID:      $Id: impliedcdsspreads.cpp,v 1.11 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/useexception.h"
#include "ito33/arraycheckers.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/error.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/finance/impliedcdsspreads.h"

extern const ito33::Error 
   ITO33_BAD_DATA,
   ITO33_DIV0;
extern const ito33::finance::Error
   ITO33_INVALID_FREQUENCY,
   ITO33_INVALID_RECOVERYRATE_1,
   ITO33_VALUATION_BEFORE_ISSUE,
   ITO33_MATURITYBEFOREVALUATION;

namespace ito33
{

namespace finance 
{


// Implementation is based on the fact that the price is a linear function
// of the spread.
ImpliedCDSSpreads::ImpliedCDSSpreads(Date contractingDate,
                                     Date firstPaymentDate,
                                     Date::DayCountConvention dcc,
                                     Frequency freq,
                                     double dRecoveryRate)                         
                                   : m_contractingDate(contractingDate),
                                     m_firstPaymentDate(firstPaymentDate),
                                     m_dcc(dcc),
                                     m_freq(freq),
                                     m_dRecoveryRate(dRecoveryRate)
{ 
  // Verify the data
  CHECK_COND_1(dRecoveryRate >= 0. && dRecoveryRate <= 1.,
               ITO33_INVALID_RECOVERYRATE_1,
               dRecoveryRate);

	CHECK_COND_MSG
  ( 
    m_firstPaymentDate > contractingDate,
    ITO33_BAD_DATA,
    TRANS("Implied CDS spread calculator: First payment date equal or "
          "before issue date.")
  );

  CHECK_COND( IsValid(m_freq), ITO33_INVALID_FREQUENCY );
}

std::vector<double> 
ImpliedCDSSpreads::Compute(const shared_ptr<TheoreticalModel>& model, 
                           const shared_ptr<SessionData>& pSessionData, 
                           const std::vector<Date>& pMaturityDates)
{
  // Clone the model, don't use the flags inside the original model
  shared_ptr<TheoreticalModel> myModel( model->Clone() );

  // Don't use flags inside Derivative either, use default flags
  myModel->SetExternalFlagsToDefaults();

  // Verify the data
  CheckIncreasingOrder(pMaturityDates);

  if ( pMaturityDates[0] < m_firstPaymentDate ) 
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Implied CDS spread calculator: Maturity dates must be " \
                  "after the first payment date.")
          );

  CHECK_COND( pSessionData->GetValuationDate() >= m_contractingDate,
              ITO33_VALUATION_BEFORE_ISSUE );

  CHECK_COND( pMaturityDates[0] > pSessionData->GetValuationDate(),
              ITO33_MATURITYBEFOREVALUATION );

  // Compute using two different spreads.
  double dSpread1 = 0.1;
  double dSpread2 = 0.2;

  // create the return vector of implied spreads
  size_t nNbDates = pMaturityDates.size();
  std::vector<double> pdImpliedSpreads(nNbDates);

  // If nothing is passed in, return an empty vector.  Could also throw
  // an exception
  if (nNbDates == 0)
    return pdImpliedSpreads;

  // The temporary output objet
  shared_ptr<ModelOutput> pMO;
  
  for (size_t nIdx = 0; nIdx < nNbDates; nIdx++)
  {
    shared_ptr<CashFlowStreamUniform> 
      pSpreadStream( new CashFlowStreamUniform
                         ( m_contractingDate,
                           m_firstPaymentDate, 
                           pMaturityDates[nIdx],
                           dSpread1,
                           m_dcc,
                           m_freq 
                         )
                   );

    shared_ptr<CDS> 
      pCDS( new CDS(m_dRecoveryRate, pSpreadStream) );

    pCDS->SetSessionData(pSessionData);

    pMO = myModel->Compute(*pCDS);
    
    double dPrice1 = pMO->GetPrice();

    pSpreadStream = make_ptr( new CashFlowStreamUniform
                                  ( m_contractingDate,
                                    m_firstPaymentDate, 
                                    pMaturityDates[nIdx],
                                    dSpread2,
                                    m_dcc,
                                    m_freq 
                                  ) );

    pCDS = make_ptr( new CDS(m_dRecoveryRate, pSpreadStream) );

    pCDS->SetSessionData(pSessionData);

    pMO = myModel->Compute(*pCDS); 

    double dPrice2 = pMO->GetPrice();

    double dSlope = (dPrice2 - dPrice1) / (dSpread2 - dSpread1);
    double dC = dPrice1 - dSlope * dSpread1;
    
    CHECK_COND( fabs(dSlope) > 1.e-20, ITO33_DIV0 );

    pdImpliedSpreads[nIdx] = - dC / dSlope;
  }

  return pdImpliedSpreads;
}


} // namespace finance

} // namespace ito33
