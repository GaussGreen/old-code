/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/com/terms_impl_for_inclusion.cpp
// Purpose:     Implementation for functions creating pure financial object
// Created:     2006/06/12
// RCS-ID:      $Id: terms_impl_for_inclusion.cpp,v 1.11 2006/08/22 08:29:48 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

STDMETHODIMP
ClassImpl::NewTimeFunction(SAFEARRAY(DATE) *dates,
                           SAFEARRAY(double) *values,
                           ITimeFunction **ppTimeFunction)
{
  try
  {
    *ppTimeFunction = new TimeFunctionImpl(make_ptr(
                         new finance::TimeFunction
                             (C2A::COM::ToVector<Date, DATE>(dates),
                              C2A::COM::ToVector<double, double>(values))
                       ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}


STDMETHODIMP
ClassImpl::NewNumeraire(BSTR currencyCode,
                        INumeraire **ppNumeraire)
{
  try
  {
    *ppNumeraire = new NumeraireImpl(make_ptr(
                         new finance::Numeraire
                             (C2A::COM::Translate<std::string>::From(currencyCode))
                       ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewYieldCurveFlat(double dRate, IYieldCurveFlat **ppYC)
{
  try
  {
    *ppYC = new YieldCurveFlatImpl(make_ptr(new finance::YieldCurveFlat(dRate)));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewYieldCurveAnnuallyCompounded(DATE dateOrigin, 
                                           IYieldCurveAnnuallyCompounded **ppYC)
{
  try
  {
    *ppYC = new YieldCurveAnnuallyCompoundedImpl(make_ptr(
                  new finance::YieldCurveAnnuallyCompounded
                      ( C2A::COM::Translate<ito33::Date>::From(dateOrigin) )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewYieldCurveSwap(DATE dateOrigin, IYieldCurveSwap **ppYC)
{
  try
  {
    *ppYC = new YieldCurveSwapImpl(make_ptr(
                  new finance::YieldCurveSwap
                      ( C2A::COM::Translate<ito33::Date>::From(dateOrigin) )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewYieldCurveZeroCoupon(DATE dateOrigin, 
                                   IYieldCurveZeroCoupon **ppYC)
{
  try
  {
    *ppYC = new YieldCurveZeroCouponImpl(make_ptr(
                  new finance::YieldCurveZeroCoupon
                      ( C2A::COM::Translate<ito33::Date>::From(dateOrigin) )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewMoneyMarket(INumeraire *pNumeraire,
                          IYieldCurve *pYieldCurve,                         
                          IMoneyMarket **ppMM)
{
  try
  {
    *ppMM = new MoneyMarketImpl(make_ptr(
                  new finance::MoneyMarket
                      (
                        C2A::COM::Translate<INumeraire *>::From(pNumeraire),
                        C2A::COM::Translate<IYieldCurve *>::From(pYieldCurve)     
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP ClassImpl::NewEquity(double dSpotSharePrice,
                                  INumeraire *pNumeraire,
                                  IEquity **ppEquity)
{
  try
  {
      *ppEquity = new EquityImpl(make_ptr(
                        new finance::Equity
                            (
                              dSpotSharePrice,
                              C2A::COM::Translate<INumeraire *>::From(pNumeraire)
                            )  
                      ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP ClassImpl::NewSessionData(IRateData *pRateData,
                                       IEquity *pEquity,
                                       DATE valuationDate,
                                       ISessionData **ppSessionData)
{
  try
  {
    *ppSessionData = new SessionDataImpl(make_ptr(
                           new finance::SessionData
                               (
                                 C2A::COM::Translate<IRateData*>::From(pRateData),
                                 C2A::COM::Translate<IEquity*>::From(pEquity),
                                 C2A::COM::Translate<ito33::Date>::From(valuationDate)
                               )
                         ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewCashFlowStreamUniformWithLastPaymentIrregular(DATE contractingDate,
                                                            DATE firstPaymentDate,
                                                            DATE lastPaymentDate,
                                                            double annualAmt,
                                                            DayCountConvention dcc,
                                                            Frequency freq,
                                                            LastPaymentType lastPaymentType,
                                                            ICashFlowStreamUniform **ppUFS)
{
  try
  {
    *ppUFS = new CashFlowStreamUniformImpl(make_ptr(
                   new finance::CashFlowStreamUniform
                       (
                         C2A::COM::Translate<ito33::Date>::From(contractingDate),
                         C2A::COM::Translate<ito33::Date>::From(firstPaymentDate),
                         C2A::COM::Translate<ito33::Date>::From(lastPaymentDate),
                         annualAmt,
                         static_cast<Date::DayCountConvention>(dcc),
                         static_cast<finance::Frequency>(freq),
                         static_cast<finance::LastPaymentType>(lastPaymentType)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewCashFlowStreamUniform(DATE contractingDate,
                                    DATE firstPaymentDate,
                                    DATE lastPaymentDate,
                                    double annualAmt,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    ICashFlowStreamUniform **ppUFS)
{
  try
  {
    *ppUFS = new CashFlowStreamUniformImpl(make_ptr(
                   new finance::CashFlowStreamUniform
                       (
                         C2A::COM::Translate<ito33::Date>::From(contractingDate),
                         C2A::COM::Translate<ito33::Date>::From(firstPaymentDate),
                         C2A::COM::Translate<ito33::Date>::From(lastPaymentDate),
                         annualAmt,
                         static_cast<Date::DayCountConvention>(dcc),
                         static_cast<finance::Frequency>(freq)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewCashFlowStreamUniformAdjusted(DATE contractingDate,
                                            SAFEARRAY(DATE) *paymentDates,
                                            SAFEARRAY(double) *paymentAmounts,
                                            double annualAmt,
                                            DayCountConvention dcc,
                                            Frequency freq,
                                            ICashFlowStreamUniform **ppUFS)
{
  try
  {
    *ppUFS = new CashFlowStreamUniformImpl(make_ptr(
                   new finance::CashFlowStreamUniform
                       (
                         C2A::COM::Translate<ito33::Date>::From(contractingDate),
                         C2A::COM::ToVector<Date, DATE>(paymentDates),
                         C2A::COM::ToVector<double, double>(paymentAmounts),
                         annualAmt,
                         static_cast<Date::DayCountConvention>(dcc),
                         static_cast<finance::Frequency>(freq)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewCashFlowStreamUniformSimple(DATE contractingDate,
                                          long nMonths,
                                          Frequency freq,
                                          double annualAmt,
                                          DayCountConvention dcc,                                     
                                          ICashFlowStreamUniform **ppUFS)
{
  try
  {
    *ppUFS = new CashFlowStreamUniformImpl(make_ptr(
                   new finance::CashFlowStreamUniform
                       (
                         C2A::COM::Translate<ito33::Date>::From(contractingDate),
                         nMonths,
                         static_cast<finance::Frequency>(freq),
                         annualAmt,
                         static_cast<Date::DayCountConvention>(dcc)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
STDMETHODIMP 
ClassImpl::NewCashFlowStreamGeneral(DATE contractingDate,
                                    SAFEARRAY(DATE) *paymentDates,
                                    SAFEARRAY(double) *paymentAmounts,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    ICashFlowStreamGeneral **ppCFG)
{
  try
  {
    *ppCFG = new CashFlowStreamGeneralImpl(make_ptr(
                   new finance::CashFlowStreamGeneral
                       (
                         C2A::COM::Translate<ito33::Date>::From(contractingDate),
                         C2A::COM::ToVector<Date, DATE>(paymentDates),
                         C2A::COM::ToVector<double, double>(paymentAmounts),
                         static_cast<Date::DayCountConvention>(dcc),
                         static_cast<finance::Frequency>(freq)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewFloatingRates(double dMargin,
                            DATE startOfAccruedDate, 
                            DATE firstUnknownPaymentDate,
                            DATE lastUnknownPaymentDate,
                            Frequency freq,
                            IFloatingRates **ppFR)
{
  try
  {
    *ppFR = new FloatingRatesImpl(make_ptr(
                  new finance::FloatingRates
                      (
                        dMargin,
                        C2A::COM::Translate<ito33::Date>::From(startOfAccruedDate),
                        C2A::COM::Translate<ito33::Date>::From(firstUnknownPaymentDate),
                        C2A::COM::Translate<ito33::Date>::From(lastUnknownPaymentDate),
                        static_cast<finance::Frequency>(freq)
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewFloatingRatesWithNoUnknownPayment(DATE startOfAccruedDate,
                                                Frequency freq,
                                                IFloatingRates **ppFR)
{
  try
  {
    *ppFR = new FloatingRatesImpl(make_ptr(
                  new finance::FloatingRates
                      (
                        C2A::COM::Translate<ito33::Date>::From(startOfAccruedDate),
                        static_cast<finance::Frequency>(freq)
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewFloatingRatesWithLastPaymentIrregular(double dMargin,
                                                    DATE startOfAccruedDate,  
                                                    DATE firstUnknownPaymentDate,
                                                    DATE lastUnknownPaymentDate,
                                                    Frequency freq,
                                                    LastPaymentType lastPaymentType,
                                                    IFloatingRates **ppFR)
{
  try
  {
    *ppFR = new FloatingRatesImpl(make_ptr(
                  new finance::FloatingRates
                      (  
                        dMargin,
                        C2A::COM::Translate<ito33::Date>::From(startOfAccruedDate),
                        C2A::COM::Translate<ito33::Date>::From(firstUnknownPaymentDate),
                        C2A::COM::Translate<ito33::Date>::From(lastUnknownPaymentDate),
                        static_cast<finance::Frequency>(freq),
                        static_cast<finance::LastPaymentType>(lastPaymentType)
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewParBond(DATE contractingDate,
                      long nMaturity,
                      double dYTM,
                      double dSpread,
                      Frequency freq,
                      DayCountConvention dcc,
                      double dRecoveryRate,
                      IParBond **ppParBond)
{
  try
  {
    *ppParBond = new ParBondImpl(make_ptr(
                       new finance::ParBond
                           (
                             C2A::COM::Translate<ito33::Date>::From(contractingDate),
                             nMaturity,
                             dYTM,
                             dSpread,
                             static_cast<finance::Frequency>(freq),
                             static_cast<Date::DayCountConvention>(dcc),
                             dRecoveryRate
                           )
                     ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewReferenceCDS(long maturity,
                           Frequency freq,
                           DayCountConvention dcc,
                           double dRecoveryRate,
                           IReferenceCDS **ppRefCDS)
{
  try
  {
    *ppRefCDS = new ReferenceCDSImpl(make_ptr(
                      new finance::ReferenceCDS
                          (
                            maturity,
                            static_cast<finance::Frequency>(freq),
                            static_cast<Date::DayCountConvention>(dcc),                         
                            dRecoveryRate
                          )
                    ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewCDS(double recoveryRate,
                  ICashFlowStreamUniform * pSpreadStream,
                  ICDS **ppCDS)
{
  try
  {
    *ppCDS = new CDSImpl(make_ptr(
                   new finance::CDS
                       (
                         recoveryRate,
                         C2A::COM::Translate<ICashFlowStreamUniform *>::From(pSpreadStream)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewEDS(double recoveryRate,
                  ICashFlowStreamUniform * pSpreadStream,
                  double barrier,
                  IEDS **ppEDS)
{
  try
  {
    *ppEDS = new EDSImpl(make_ptr(
                   new finance::EDS
                       (
                         recoveryRate,
                         C2A::COM::Translate<ICashFlowStreamUniform *>::From(pSpreadStream),
                         barrier
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewOption(double dStrike,
                     DATE maturityDate,
                     OptionType optionType,
                     ExerciseType exerciseType,
                     IOption **ppOption)
{
  try
  {
    *ppOption = new OptionImpl(make_ptr(
                      new finance::Option
                          (
                            dStrike,
                            C2A::COM::Translate<ito33::Date>::From(maturityDate),
                            static_cast<finance::OptionType>(optionType),
                            static_cast<finance::ExerciseType>(exerciseType)
                          )
                    ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewVarianceSwapTerms(DATE maturityDate,
                                SwapType swapType,
                                DATE startOfSamplingPeriod,
                                long nNbSamplingReturns,
                                IVarianceSwapTerms **ppTerms)
{
  try
  {
    *ppTerms = new VarianceSwapTermsImpl(make_ptr(
                     new finance::VarianceSwapTerms
                         (
                           C2A::COM::Translate<ito33::Date>::From(maturityDate),
                           static_cast<finance::SwapType>(swapType),
                           C2A::COM::Translate<ito33::Date>::From(startOfSamplingPeriod),
                           nNbSamplingReturns
                         )
                   ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewVarianceSwap(IVarianceSwapTerms* pTerms,
                           double dVolatilityStrike,
                           IVarianceSwap **ppVarianceSwap)
{
  try
  {
    *ppVarianceSwap = new VarianceSwapImpl(make_ptr(
                            new finance::VarianceSwap
                                ( 
                                  C2A::COM::Translate<IVarianceSwapTerms*>::From(pTerms),
                                  dVolatilityStrike
                                )
                          ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewGammaVarianceSwap(IVarianceSwapTerms* pTerms,
                                double dVolatilityStrike,
                                IGammaVarianceSwap **ppVarianceSwap)
{
  try
  {
    *ppVarianceSwap = new GammaVarianceSwapImpl(make_ptr(
                            new finance::GammaVarianceSwap
                                ( 
                                  C2A::COM::Translate<IVarianceSwapTerms*>::From(pTerms),
                                  dVolatilityStrike
                                )
                          ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewConditionalVarianceSwap(IVarianceSwapTerms* pTerms,
                                      double dVolatilityStrike,
                                      IConditionalVarianceSwap **ppVarianceSwap)
{
  try
  {
    *ppVarianceSwap = new ConditionalVarianceSwapImpl(make_ptr(
                            new finance::ConditionalVarianceSwap
                                ( 
                                  C2A::COM::Translate<IVarianceSwapTerms*>::From(pTerms),
                                  dVolatilityStrike
                                )
                          ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewOptionVarianceSwap(IVarianceSwapTerms* pTerms,
                                 double dVolatilityStrike,
                                 OptionType optionType,
                                 IOptionVarianceSwap **ppOptionVarianceSwap)
{
  try
  {
    *ppOptionVarianceSwap = new OptionVarianceSwapImpl(make_ptr(
                                  new finance::OptionVarianceSwap
                                      ( 
                                        C2A::COM::Translate<IVarianceSwapTerms*>::From(pTerms),
                                        dVolatilityStrike,
                                        static_cast<finance::OptionType>(optionType)
                                      )
                                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewVarianceSwaption(IVarianceSwapTerms* pTerms,
                               OptionType optionType,
                               double dStrike,
                               DATE maturityDate,
                               IVarianceSwaption **ppVarianceSwaption)
{
  try
  {
    *ppVarianceSwaption = new VarianceSwaptionImpl(make_ptr(
                                new finance::VarianceSwaption
                                    ( 
                                      C2A::COM::Translate<IVarianceSwapTerms*>::From(pTerms),
                                      static_cast<finance::OptionType>(optionType),
                                      dStrike,
                                      C2A::COM::Translate<ito33::Date>::From(maturityDate)
                                    )
                              ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewOneTouch(DATE maturityDate,
                       double dBarrier,
                       BarrierType barrierType,
                       RebateType rebateType,
                       IOneTouch **ppOneTouch)
{
  try
  {
    *ppOneTouch = new OneTouchImpl(make_ptr(
                        new finance::OneTouch
                            (
                              C2A::COM::Translate<ito33::Date>::From(maturityDate),
                              dBarrier,
                              static_cast<finance::BarrierType>(barrierType),
                              static_cast<finance::RebateType>(rebateType)
                            )
                      ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewFXOneTouch(DATE maturityDate,
                         double BSBarrier,                          
                         BarrierType barrierType,
                         double vol,
                         IFXOneTouch **ppFXOneTouch)
{
  try
  {
    *ppFXOneTouch = new FXOneTouchImpl(make_ptr(
                          new finance::FXOneTouch
                              (
                                C2A::COM::Translate<ito33::Date>::From(maturityDate),
                                BSBarrier,
                                static_cast<finance::BarrierType>(barrierType),
                                vol
                              )
                        ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewAveragePriceOption(double dFixedStrike,
                                 DATE maturityDate,
                                 OptionType optionType,
                                 ExerciseType exerciseType,
                                 DATE averageStartDate,
                                 long numberOfSamplingAverages,
                                 IAsianOption** ppAO)
{
  try
  {
    *ppAO = new AsianOptionImpl(make_ptr(
                  new finance::AsianOption
                      (
                        dFixedStrike,
                        C2A::COM::Translate<ito33::Date>::From(maturityDate),
                        static_cast<finance::OptionType>(optionType),
                        static_cast<finance::ExerciseType>(exerciseType),
                        C2A::COM::Translate<ito33::Date>::From(averageStartDate),
                        numberOfSamplingAverages
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewAverageStrikeOption(DATE maturityDate,
                                  OptionType optionType,
                                  ExerciseType exerciseType,
                                  DATE averageStartDate,
                                  long numberOfSamplingAverages,
                                  IAsianOption** ppAO)
{
  try
  {
    *ppAO = new AsianOptionImpl(make_ptr(
                  new finance::AsianOption
                      (
                        C2A::COM::Translate<ito33::Date>::From(maturityDate),
                        static_cast<finance::OptionType>(optionType),
                        static_cast<finance::ExerciseType>(exerciseType),
                        C2A::COM::Translate<ito33::Date>::From(averageStartDate),
                        numberOfSamplingAverages
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewBondLikeTerms(DATE issueDate, 
                            double dIssuePrice, 
                            DATE maturityDate, 
                            double dNominal,
                            double dRecoveryRate,
                            IBondLikeTerms **ppBLT)
{
  try
  {
    *ppBLT = new BondLikeTermsImpl(make_ptr(
                   new finance::BondLikeTerms
                       (
                         C2A::COM::Translate<ito33::Date>::From(issueDate),
                         dIssuePrice,
                         C2A::COM::Translate<ito33::Date>::From(maturityDate),
                         dNominal,
                         dRecoveryRate
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewBondTerms(DATE issueDate, 
                        double dIssuePrice,
                        DATE maturityDate,
                        double dNominal,                       
                        double dRedemptionPrice,
                        double dRecoveryRate,
                        IBondTerms **ppBT)
{
  try
  {
    *ppBT = new BondTermsImpl(make_ptr(
                  new finance::BondTerms
                      (
                        C2A::COM::Translate<ito33::Date>::From(issueDate),
                        dIssuePrice,
                        C2A::COM::Translate<ito33::Date>::From(maturityDate),
                        dNominal,
                        dRedemptionPrice,
                        dRecoveryRate
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewCallPeriodWithStrike(DATE startDate, 
                                   DATE endDate, 
                                   double dStrike, 
                                   ICallPeriod **ppCP)
{
  try
  {
    *ppCP = new CallPeriodImpl(
                  finance::CallPeriod::CreateWithStrike
                  ( 
                    C2A::COM::Translate<ito33::Date>::From(startDate),
                    C2A::COM::Translate<ito33::Date>::From(endDate),
                    dStrike
                  )
                );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewCallPeriodWithYield(DATE startDate, 
                                  DATE endDate, 
                                  double dYield, 
                                  ICallPeriod **ppCP)
{
  try
  {
    *ppCP = new CallPeriodImpl(
                  finance::CallPeriod::CreateWithYield
                  ( 
                    C2A::COM::Translate<ito33::Date>::From(startDate),
                    C2A::COM::Translate<ito33::Date>::From(endDate),
                    dYield
                  )
                );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewBond(IBondTerms *pBT, IBond **ppBond)
{
  try
  {
    *ppBond = new BondImpl(make_ptr(
                    new finance::Bond
                        ( C2A::COM::Translate<IBondTerms *>::From(pBT) )
                  ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewConversionPeriod(DATE startDate, 
                               DATE endDate, 
                               double dRatio, 
                               IConversionPeriod **ppCP)
{
  try
  {
    *ppCP = new ConversionPeriodImpl(make_ptr(
                  new finance::ConversionPeriod
                      ( 
                        C2A::COM::Translate<ito33::Date>::From(startDate),
                        C2A::COM::Translate<ito33::Date>::From(endDate),
                        dRatio
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewConvertibleBond(IBondTerms *pBT, 
                              IConversionSchedule *pConversionSchedule,
                              IConvertibleBond **ppCB)
{
  try
  {
    *ppCB = new ConvertibleBondImpl(make_ptr(
                  new finance::ConvertibleBond
                      ( 
                        C2A::COM::Translate<IBondTerms *>::From(pBT),
                        C2A::COM::Translate<IConversionSchedule *>::From
                        (pConversionSchedule)
                      )
                ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewCBOption(IConvertibleBond *pConvertibleBond,
                       IFloatingRates *pFloatingRates, 
                       DATE maturityDate,
                       ICBOption **ppCBOption)
{
  try
  {
    *ppCBOption = new CBOptionImpl(make_ptr(
                        new finance::CBOption
                            (   
                              C2A::COM::Translate<IConvertibleBond *>::From(pConvertibleBond),
                              C2A::COM::Translate<IFloatingRates *>::From(pFloatingRates),
                              C2A::COM::Translate<ito33::Date>::From(maturityDate)
                            )
                      ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewShareDependentConversion(DATE startDate, 
                                       DATE endDate, 
                                       double dBaseRatio, 
                                       double dIncrementalShareFactor, 
                                       IShareDependentConversion **ppSDC)
{
  try
  {
    *ppSDC = new ShareDependentConversionImpl(make_ptr(
                   new finance::ShareDependentConversion
                       (                            
                         C2A::COM::Translate<ito33::Date>::From(startDate),
                         C2A::COM::Translate<ito33::Date>::From(endDate),
                         dBaseRatio, 
                         dIncrementalShareFactor
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewAttachedWarrantConvertibleBond(IBondTerms *pBondTerms,
                                             IShareDependentConversion *pSDC,
                                             IAttachedWarrantConvertibleBond **ppAWCB)
{
  try
  {
    *ppAWCB = new AttachedWarrantConvertibleBondImpl(make_ptr(
                    new finance::AttachedWarrantConvertibleBond
                        ( 
                          C2A::COM::Translate<IBondTerms *>::From(pBondTerms),
                          C2A::COM::Translate<IShareDependentConversion *>::From(pSDC)
                        )
                  ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewConversionPriceReset(DATE resetDate,
                                   double dFloorRate,
                                   IConversionPriceReset **ppCPR)
{
  try
  {
    *ppCPR = new ConversionPriceResetImpl(make_ptr(
                   new finance::ConversionPriceReset
                       ( 
                         C2A::COM::Translate<ito33::Date>::From(resetDate),
                         dFloorRate
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}


STDMETHODIMP 
ClassImpl::NewResetConversionSchedule(DATE startDate,
                                      DATE endDate, 
                                      double dInitialConversionPrice,
                                      double dCurrentConversionPrice,
                                      ResetFlooredBy flooredBy,
                                      IResetConversionSchedule **ppRCS)
{
  try
  {
    *ppRCS = new ResetConversionScheduleImpl(make_ptr(
                   new finance::ResetConversionSchedule
                       ( 
                         C2A::COM::Translate<ito33::Date>::From(startDate),
                         C2A::COM::Translate<ito33::Date>::From(endDate),
                         dInitialConversionPrice,
                         dCurrentConversionPrice,
                         static_cast<finance::ResetFlooredBy>(flooredBy)
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewReset(IBondTerms *pBT, 
                    IResetConversionSchedule *pRCS,
                    IReset **ppReset)
{
  try
  {
    *ppReset = new ResetImpl(make_ptr(
                     new finance::Reset
                         ( 
                           C2A::COM::Translate<IBondTerms *>::From(pBT),
                           C2A::COM::Translate<IResetConversionSchedule *>::From(pRCS)
                         )
                   ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewCallfixedShare(DATE startDate, 
                             DATE endDate, 
                             double dRatio,
                             ICallFixedShare **ppCFS)
{
  try
  {
    *ppCFS = new CallFixedShareImpl(make_ptr(
                   new finance::CallFixedShare
                       ( 
                         C2A::COM::Translate<ito33::Date>::From(startDate),
                         C2A::COM::Translate<ito33::Date>::From(endDate),
                         dRatio
                       )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewPERCSLike(IBondLikeTerms *pBondLikeTerms,
                        double dCapPrice, 
                        double dConversionRatio,
                        IPERCSLike **ppPERCS)
{
  try
  {
    *ppPERCS = new PERCSLikeImpl(make_ptr(
                     new finance::PERCSLike
                         ( 
                           C2A::COM::Translate<IBondLikeTerms *>::From(pBondLikeTerms),
                           dCapPrice,
                           dConversionRatio
                         )
                   ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewPEPSLike(IBondLikeTerms *pBondLikeTerms,
                       double dMaxRatio, 
                       double dMinRatio,
                       IPEPSLike **ppPEPS)
{
  try
  {
    *ppPEPS = new PEPSLikeImpl(make_ptr(
                    new finance::PEPSLike
                        ( 
                          C2A::COM::Translate<IBondLikeTerms *>::From(pBondLikeTerms),
                          dMaxRatio,
                          dMinRatio
                        )
                  ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP
ClassImpl::NewPEPSAveragingPeriodWithStock(DATE startDate, 
                                           DATE endDate,
                                           long numberOfSamplingAverages,
                                           IPEPSAveragingPeriod** ppPAP)
{
  try
  {
    *ppPAP = new PEPSAveragingPeriodImpl(
                   finance::PEPSAveragingPeriod::CreateWithStock
                   (
                     C2A::COM::Translate<ito33::Date>::From(startDate),
                     C2A::COM::Translate<ito33::Date>::From(endDate),
                     numberOfSamplingAverages
                   )
                 );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewPEPSAveragingPeriodWithConversionRatio(DATE startDate, 
                                                     DATE endDate,
                                                     long numberOfSamplingAverages,
                                                     IPEPSAveragingPeriod** ppPAP)
{
  try
  {
    *ppPAP = new PEPSAveragingPeriodImpl(
                   finance::PEPSAveragingPeriod::CreateWithConversionRatio
                   (
                     C2A::COM::Translate<ito33::Date>::From(startDate),
                     C2A::COM::Translate<ito33::Date>::From(endDate),
                     numberOfSamplingAverages
                   )
                 );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewGeneralizedPEPSLikeCall(DATE startDate, 
                                      DATE endDate, 
                                      double dTriggerRate,
                                      GeneralizedPEPSLikeCallType type,
                                      IGeneralizedPEPSLikeCall **ppCall)
{
  try
  {
    *ppCall = new GeneralizedPEPSLikeCallImpl(make_ptr(
                    new finance::GeneralizedPEPSLikeCall
                        ( 
                          C2A::COM::Translate<ito33::Date>::From(startDate),
                          C2A::COM::Translate<ito33::Date>::From(endDate),
                          dTriggerRate,
                          static_cast<finance::GeneralizedPEPSLikeCallType>(type)
                        )
                 ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ClassImpl::NewGeneralizedPEPSLike(IBondLikeTerms *pBondLikeTerms,
                                  double dDownsideConversionRatio, 
                                  double dLowerStrike,
                                  double dUpsideBaseConversionRatio, 
                                  double dHigherStrike,
                                  IGeneralizedPEPSLike **ppPEPS)
{
  try
  {
    *ppPEPS = new GeneralizedPEPSLikeImpl(make_ptr(
                    new finance::GeneralizedPEPSLike
                        ( 
                          C2A::COM::Translate<IBondLikeTerms *>::From(pBondLikeTerms),
                          dDownsideConversionRatio,
                          dLowerStrike,
                          dUpsideBaseConversionRatio,
                          dHigherStrike
                        )
                  ));

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
