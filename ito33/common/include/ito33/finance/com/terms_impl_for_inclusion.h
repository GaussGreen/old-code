/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/com/terms_impl_for_inclusion.h
// Purpose:     declaration of class imlpementing ITerms interface
// Created:     2006/06/08
// RCS-ID:      $Id: terms_impl_for_inclusion.h,v 1.10 2006/08/22 08:29:48 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_FINANCE_COM_TERMS_IMPL_FOR_INCLUSION_H_
#define _ITO33_FINANCE_COM_TERMS_IMPL_FOR_INCLUSION_H_

STDMETHODIMP NewTimeFunction(SAFEARRAY(DATE) *dates,
                             SAFEARRAY(double) *values,
                             ITimeFunction **ppTimeFunction);

STDMETHODIMP NewNumeraire(BSTR currencyCode,
                          INumeraire **ppNumeraire);

STDMETHODIMP NewYieldCurveFlat(double dRate,
                               IYieldCurveFlat **ppYC);

STDMETHODIMP NewYieldCurveAnnuallyCompounded(DATE dateOrigin,
                                             IYieldCurveAnnuallyCompounded **ppYC);

STDMETHODIMP NewYieldCurveSwap(DATE dateOrigin, IYieldCurveSwap **ppYC);

STDMETHODIMP NewYieldCurveZeroCoupon(DATE dateOrigin,
                                     IYieldCurveZeroCoupon **ppYC);

STDMETHODIMP NewMoneyMarket(INumeraire *pNumeraire,
                            IYieldCurve *pYC,
                            IMoneyMarket **ppMM);

STDMETHODIMP NewEquity(double dSpotSharePrice,
                       INumeraire *pNumeraire,
                       IEquity **ppEquity);


STDMETHODIMP NewSessionData(IRateData *pRateData,
                            IEquity *pEquity,
                            DATE valuationDate,
                            ISessionData **ppSessionData);

STDMETHODIMP NewCashFlowStreamUniformWithLastPaymentIrregular(DATE contractingDate,
                                                              DATE firstPaymentDate,
                                                              DATE lastPaymentDate,
                                                              double annualAmt,
                                                              DayCountConvention dcc,
                                                              Frequency freq,
                                                              LastPaymentType lastPaymentType,
                                                              ICashFlowStreamUniform **ppUFS);

STDMETHODIMP NewCashFlowStreamUniform(DATE contractingDate,
                                      DATE firstPaymentDate,
                                      DATE lastPaymentDate,
                                      double annualAmt,
                                      DayCountConvention dcc,
                                      Frequency freq,
                                      ICashFlowStreamUniform **ppUFS);

STDMETHODIMP NewCashFlowStreamUniformSimple(DATE contractingDate,
                                            long nMonths,
                                            Frequency freq,
                                            double annualAmt,
                                            DayCountConvention dcc,                                     
                                            ICashFlowStreamUniform **ppUFS);

STDMETHODIMP NewCashFlowStreamUniformAdjusted(DATE contractingDate,
                                              SAFEARRAY(DATE) *paymentDates,
                                              SAFEARRAY(double) *paymentAmounts,
                                              double annualAmt,
                                              DayCountConvention dcc,
                                              Frequency freq,
                                              ICashFlowStreamUniform **ppUFS);

STDMETHODIMP NewCashFlowStreamGeneral(DATE contractingDate,
                                      SAFEARRAY(DATE) *paymentDates,
                                      SAFEARRAY(double) *paymentAmounts,
                                      DayCountConvention dcc,
                                      Frequency freq,
                                      ICashFlowStreamGeneral **ppCFG);

STDMETHODIMP NewFloatingRates(double dMargin,
                              DATE startOfAccruedDate, 
                              DATE firstUnknownPaymentDate,
                              DATE lastUnknownPaymentDate,
                              Frequency freq,
                              IFloatingRates **ppFR);

STDMETHODIMP NewFloatingRatesWithNoUnknownPayment(DATE startOfAccruedDate,
                                                  Frequency freq,
                                                  IFloatingRates **ppFR);

STDMETHODIMP NewFloatingRatesWithLastPaymentIrregular(double dMargin,
                                                      DATE startOfAccruedDate, 
                                                      DATE firstUnknownPaymentDate,
                                                      DATE lastUnknownPaymentDate,
                                                      Frequency freq,
                                                      LastPaymentType lastPaymentType,
                                                      IFloatingRates **ppFR);

STDMETHODIMP NewParBond(DATE contractingDate,
                        long nMaturity,
                        double dYTM,
                        double dSpread,
                        Frequency freq,
                        DayCountConvention dcc,
                        double dRecoveryRate,
                        IParBond **ppParBond);

STDMETHODIMP NewReferenceCDS(long maturity,
                             Frequency freq,
                             DayCountConvention dcc,
                             double dRecoveryRate,
                             IReferenceCDS **ppRefCDS);

STDMETHODIMP NewCDS(double recoveryRate,
                    ICashFlowStreamUniform * pSpreadStream,
                    ICDS **ppCDS);

STDMETHODIMP NewEDS(double recoveryRate,
                    ICashFlowStreamUniform *pSpreadStream,
                    double barrier,
                    IEDS **ppEDS);

STDMETHODIMP NewOption(double dStrike,
                       DATE maturityDate,
                       OptionType optionType,
                       ExerciseType exerciseType,
                       IOption **ppOption);

STDMETHODIMP NewVarianceSwapTerms(DATE maturityDate,
                                  SwapType swapType,
                                  DATE startOfSamplingPeriod,
                                  long nbSamplingReturns,
                                  IVarianceSwapTerms **ppVarianceSwapTerms);

STDMETHODIMP NewVarianceSwap(IVarianceSwapTerms* pTerms,
                             double dVolatilityStrike,
                             IVarianceSwap **ppVarianceSwap);

STDMETHODIMP NewGammaVarianceSwap(IVarianceSwapTerms* pTerms,
                                  double dVolatilityStrike,
                                  IGammaVarianceSwap **ppVarianceSwap);

STDMETHODIMP NewConditionalVarianceSwap(IVarianceSwapTerms* pTerms,
                                        double dVolatilityStrike,
                                        IConditionalVarianceSwap **ppVarianceSwap);

STDMETHODIMP NewOptionVarianceSwap(IVarianceSwapTerms* pTerms,
                                   double dVolatilityStrike,
                                   OptionType optionType,
                                   IOptionVarianceSwap **ppOptionVarianceSwap);

STDMETHODIMP NewVarianceSwaption(IVarianceSwapTerms* pTerms,
                                 OptionType optionType,
                                 double dStrike,
                                 DATE maturityDate,
                                 IVarianceSwaption **ppVarianceSwaption);

STDMETHODIMP NewOneTouch(DATE maturityDate,
                         double barrier,                          
                         BarrierType barrierType,
                         RebateType rebateType,
                         IOneTouch **ppOneTouch);

STDMETHODIMP NewFXOneTouch(DATE maturityDate,
                           double BSBarrier,                          
                           BarrierType barrierType,
                           double vol,
                           IFXOneTouch **ppFXOneTouch);

STDMETHODIMP NewAveragePriceOption(double dFixedStrike,
                                   DATE maturityDate,
                                   OptionType optionType,
                                   ExerciseType exerciseType,
                                   DATE averageStartDate,
                                   long numberOfSamplingAverages,
                                   IAsianOption** ppAO);

STDMETHODIMP NewAverageStrikeOption(DATE maturityDate,
                                    OptionType optionType,
                                    ExerciseType exerciseType,
                                    DATE averageStartDate,
                                    long numberOfSamplingAverages,
                                    IAsianOption** ppAO);

STDMETHODIMP NewBondLikeTerms(DATE issueDate, 
                              double dIssuePrice, 
                              DATE maturityDate, 
                              double dNominal,
                              double dRecoveryRate,
                              IBondLikeTerms **ppBLT);

STDMETHODIMP NewBondTerms(DATE issueDate, 
                          double dIssuePrice,
                          DATE maturityDate,
                          double dNominal,                       
                          double dRedemptionPrice,
                          double dRecoveryRate,
                          IBondTerms **ppBT);

STDMETHODIMP NewBond(IBondTerms *pBT, IBond **ppBond);

STDMETHODIMP NewCallPeriodWithStrike(DATE startDate, 
                                     DATE endDate, 
                                     double dStrike, 
                                     ICallPeriod **ppCP);

STDMETHODIMP NewCallPeriodWithYield(DATE startDate, 
                                    DATE endDate, 
                                    double dYield, 
                                    ICallPeriod **ppCP);

STDMETHODIMP NewConversionPeriod(DATE startDate, 
                                 DATE endDate, 
                                 double dRatio,
                                 IConversionPeriod **ppCP);

STDMETHODIMP NewConvertibleBond(IBondTerms *pBT, 
                                IConversionSchedule *pConversionSchedule,
                                IConvertibleBond **ppCB);

STDMETHODIMP NewCBOption(IConvertibleBond *pConvertibleBond,
                         IFloatingRates *pFloatingRates, 
                         DATE maturityDate,
                         ICBOption **ppCBOption);

STDMETHODIMP NewShareDependentConversion(DATE startDate,
                                         DATE endDate,  
                                         double dBaseRatio,
                                         double dIncrementalShareFactor,
                                         IShareDependentConversion **ppSDC);

STDMETHODIMP NewAttachedWarrantConvertibleBond(IBondTerms *pBondTerms,
                                               IShareDependentConversion *pSDC,
                                               IAttachedWarrantConvertibleBond **ppAWCB);

STDMETHODIMP NewConversionPriceReset(DATE resetDate,
                                     double dFloorRate,
                                     IConversionPriceReset **ppCPR);

STDMETHODIMP NewResetConversionSchedule(DATE startDate,
                                        DATE endDate, 
                                        double dInitialConversionPrice,
                                        double dCurrentConversionPrice,
                                        ResetFlooredBy flooredBy,
                                        IResetConversionSchedule **ppRCS);

STDMETHODIMP NewReset(IBondTerms *pBT, 
                      IResetConversionSchedule *pRCS,
                      IReset **ppReset);

STDMETHODIMP NewCallfixedShare(DATE startDate, 
                               DATE endDate, 
                               double dRatio,
                               ICallFixedShare **ppCFS);

STDMETHODIMP NewPERCSLike(IBondLikeTerms *pBondLikeTerms,
                          double dCapPrice, 
                          double dConversionRatio,
                          IPERCSLike **ppPERCS);

STDMETHODIMP NewPEPSLike(IBondLikeTerms *pBondLikeTerms,
                         double dMaxRatio, 
                         double dMinRatio,
                         IPEPSLike **ppPEPS);

STDMETHODIMP NewPEPSAveragingPeriodWithStock(DATE startDate, 
                                             DATE endDate,
                                             long numberOfSamplingAverages,
                                             IPEPSAveragingPeriod** ppPAP);

STDMETHODIMP NewPEPSAveragingPeriodWithConversionRatio(DATE startDate, 
                                                       DATE endDate,
                                                       long numberOfSamplingAverages,
                                                       IPEPSAveragingPeriod** ppPAP);

STDMETHODIMP NewGeneralizedPEPSLikeCall(DATE startDate, 
                                        DATE endDate, 
                                        double dTriggerRate,
                                        GeneralizedPEPSLikeCallType type,
                                        IGeneralizedPEPSLikeCall **ppCall);

STDMETHODIMP NewGeneralizedPEPSLike(IBondLikeTerms *pBondLikeTerms,
                                    double dDownsideConversionRatio, 
                                    double dLowerStrike,
                                    double dUpsideBaseConversionRatio, 
                                    double dHigherStrike,
                                    IGeneralizedPEPSLike **ppPEPS);

#endif // _ITO33_FINANCE_COM_TERMS_IMPL_FOR_INCLUSION_H_
