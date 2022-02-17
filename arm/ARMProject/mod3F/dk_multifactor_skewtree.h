#include "DKMaille.h"
#include "DKMaille2D.h"

void dk_1factor_skewtree(int no_of_steps,DKMaille<double> &t,DKMaille<double> &df,
                         double dVolatilityChoice,double dBetaA,double dSigmaA1,double dSigmaA2,
                         double dSmile,int GreenIndex,DKMaille<double> *GreensFunction,
                         DKMaille<double> *localvol,DKMaille<double> *HJM,
                         DKMaille<double> *LowestG);

double get_rate_dk1factor_skewtree(DKMaille<double> &rates,DKMaille<double> &time,
                                   DKMaille<double> &localvol, int i, int j,double smile0,
                                   double forate);

void dk_2factor_skewtree(int no_of_steps,DKMaille<double> &t,DKMaille<double> &df,
                         double dVolatilityChoice,double dBetaA,double dSigmaA1,double dSigmaA2,
                         double dBetaB,double dSigmaB1,double dSigmaB2,double dCorrelation,
                         double dSmile,int greenindex,DKMaille2D<double> *Green,
                         DKMaille<double> *localvolA,DKMaille<double> *localvolB,DKMaille<double> *HJM,
                         DKMaille<double> *LowestG);


void dk_2factor_skewtree_quanto(int no_of_steps,DKMaille<double> &t,DKMaille<double> &dfBase,DKMaille<double> &dfForeign,
                                double dVolatilityChoice,double dBetaBase,double dSigma1Base,double dSigma2Base,
                                double dBetaForeign,double dSigma1Foreign,double dSigma2Foreign,double dCorrelationBaseForeign,
                                double dSmileBase,double dSmileForeign,int greenindex,DKMaille2D<double> *GreenBase,
                                DKMaille2D<double> *GreenForeign,DKMaille<double> *localvolBase,DKMaille<double> *localvolForeign,
                                DKMaille<double> *HJMBase,DKMaille<double> *HJMForeign,DKMaille<double> *LowestGBase,
                                DKMaille<double> *LowestGForeign,double dFXForeignCorrelation,DKMaille<double> *UniGreenBase,
                                DKMaille<double> *UniGreenForeign,double dFXVolatility);

double get_rate_dk2factor_skewtree(DKMaille<double> &rates,DKMaille<double> &time, 
                                DKMaille<double> &localvolA,DKMaille<double> &localvolB,
                                int i,int j,int k,double smile0,double forate);

double rebuild_dk1factor_gtree(DKMaille<double> &rates, DKMaille<double> &time, DKMaille<double> &vol, int i, int j);

double TreeDK1F_SWAPTION_SYSTEM(int frequency,int swapmaturity,int optionexpiration,int today,double strike,
                                char type[300],
                                DKMaille<double> times,DKMaille<double> rates,DKMaille<double> vols,
                                DKMaille2D<double> GreensFunction,int steps,double smile1,
                                double smile2,double smile3,double smile4,int smoothing_flag,
                                DKMaille<double> forates,DKMaille<double> payment_dates);

double TreeDK1F_SWAPTION_BERMUDA_SYSTEM(int frequency,int strikesize,DKMaille<double> &swapmaturityset,
                                        DKMaille<double> &optionexpirationset,int today,DKMaille<double> &strikeset,
                                        char type[300],double notional,DKMaille<double> &penalty,DKMaille<double> &t,
                                        DKMaille<double> &rate,DKMaille<double> &volatility,double smile0,int no_of_steps,
                                        DKMaille<double> &forates,DKMaille<double> &payment_dates);

double TreeDK1F_SWAPTION_BCIF_SYSTEM(int output,int strikesize,DKMaille<double> &swapmaturityset,
                                     DKMaille<double> &optionexpirationset,int today,DKMaille<double> &strikeset,
                                     char type[300],int capfloor,double notional,DKMaille<double> &spread_on_libor,
                                     DKMaille<double> &t,DKMaille<double> &rate,DKMaille<double> &volatility,
                                     double smile0,double smile0Cap,int no_of_steps,DKMaille<double> &forates,
                                     DKMaille<double> &payment_dates,DKMaille<double> &rateCap,
                                     DKMaille<double> &volatilityCap,DKMaille<double> &foratesCap);

double TreeDK2F_SWAPTION_SYSTEM(int frequency,double swapmaturity,double optionexpiration,double today,double strike,char type[300],
                                DKMaille<double> &times,DKMaille<double> &rates,DKMaille<double> &vols,
						                    DKMaille<double> &volsB,DKMaille2D<double> &GreensFunction,int steps,double smile1,
						                    double smile2,double smile3,double smile4,double dCorrelation,
                                DKMaille<double> &forates,DKMaille<double> &payment_dates,int smooching_flag);

double TreeDK2F_SWAPTION_BERMUDA_SYSTEM(int frequency,int strikesize,DKMaille<double> &swapmaturityset,
                                        DKMaille<double> &optionexpirationset,int today,DKMaille<double> &strikeset,
                                        char type[300],double notional,DKMaille<double> &penalty,DKMaille<double> &t,
                                        DKMaille<double> &rate,DKMaille<double> &volatility,DKMaille<double> &volatilityB,
                                        double smile0,int no_of_steps,double dCorrelation,DKMaille<double> &forates,
                                        DKMaille<double> &payment_dates);

double TreeDK2F_SWAPTION_CMS_SYSTEM(int frequency,int strikesize,DKMaille<double> &strikeset,
                                        DKMaille<double> &optionexpirationset,int today,DKMaille<double> &funding,
                                        char type[300],double notional,DKMaille<double> &t,
                                        DKMaille<double> &rate,DKMaille<double> &volatility,DKMaille<double> &volatilityB,
                                        double smile0,int no_of_steps,double dCorrelation,DKMaille<double> &forates,
                                        DKMaille<double> &payment_dates,int firstoption,int lastoption);

double TreeDK1F_SWAPTION_CMS_SYSTEM(int frequency,int strikesize,DKMaille<double> &strikeset,
                                        DKMaille<double> &optionexpirationset,int today,DKMaille<double> &funding,
                                        char type[300],double notional,DKMaille<double> &t,
                                        DKMaille<double> &rate,DKMaille<double> &volatility,
                                        double smile0,int no_of_steps,DKMaille<double> &forates,
                                        DKMaille<double> &payment_dates,int firstoption,int lastoption);

double TreeDK1F_SWAPTION_SYSTEM_ANALYTICS(int frequency,double today,double strike,char type[300],
                                DKMaille<double> times,DKMaille<double> rates,DKMaille<double> vols,
                                DKMaille<double> GreensFunction,int steps,double smile1,
                                double smile2,double smile3,double smile4,int smoothing_flag,DKMaille<double> forates,
                                DKMaille<double> payment_dates,int no_of_paymentdates,DKMaille<double> short_rate,
                                double dt, DKMaille<double> modelcoeffs,
                                DKMaille<double> T,DKMaille<double> R,int no_of_rates);

void Bondcalc_Analytics(DKMaille<double> T,DKMaille<double> R,int number,int j1, 
							 double paymentdate,DKMaille<double> short_rate,double dtshort,
							 double expirydate,DKMaille<double> modelcoeffs,DKMaille<double> green,
               DKMaille<double> *bond,double smile0);

void Bondcalc(int j2, int j1,DKMaille<double> rates,DKMaille<double> vols,DKMaille<double> times,
              DKMaille2D<double> *bond,double smile0,DKMaille<double> forates);

double TreeDK2F_SWAPTION_BCIF_SYSTEM(int strikesize,
                                     DKMaille<double> &swapmaturityset,
                                     DKMaille<double> &optionexpirationset,
                                     int today,
                                     DKMaille<double> &strikeset,
                                     char type[300],
                                     int capfloor,
                                     double notional,
                                     DKMaille<double> &funding,
                                     DKMaille<double> &t,
                                     DKMaille<double> &rate,
                                     DKMaille<double> &volatility,
                                     DKMaille<double> &volatilityB,
                                     double dCorrelation,
                                     double smile0,
                                     double smile0Cap,
                                     int no_of_steps,
                                     DKMaille<double> &forates,
                                     DKMaille<double> &payment_dates,
                                     DKMaille<double> &rateCap,
                                     DKMaille<double> &volatilityCap,
                                     DKMaille<double> &volatilityBCap,
                                     double dCorrelationCap,
                                     DKMaille<double> &foratesCap);

void CalcSwapRate(int freq,int i,int *paymentnode,int payments,DKMaille<double> rate,DKMaille<double> volatility,DKMaille<double> t, 
                  DKMaille<double> *swapRate,DKMaille<double> *annuity,int no_of_steps,double smile0,DKMaille<double> forates);

double HJMDrift(double short_g,DKMaille<double> &t,DKMaille<double> &df,double dVolatilityChoice,
                         double dBetaA,double dSigmaA1,double dSigmaA2,double smile0);

void CalcSwapRateAnalytics(int freq,int i,DKMaille<double> payment_dates,int payments,DKMaille<double> t, 
                           DKMaille<double> *swapRate,DKMaille<double> *annuity,int no_of_steps,double smile0,
                           DKMaille<double> short_rate,double dt,DKMaille<double> modelcoeffs,DKMaille<double> T,
                           DKMaille<double> R,int no_of_rates,DKMaille<double> green);


void CreateStaticFXLattice(int no_of_steps, int iActualIndex,DKMaille<double> &t_sys,
                          double dSpotFX,double dFXVolatility,
                          DKMaille<double> *SpotFXTable);

void CreateStaticFXLatticeNew(int no_of_steps, int iActualIndex,
                         DKMaille<double> &t_sys,double dSpotFX,double dFXVolatility,
                         int iAddToTop,int iAddToBottom,DKMaille<double> *SpotFXTableNew);

void FindEdgesFXLattice(int no_of_steps,int iSize,DKMaille<double> &t_sys,
                        double dSpotFX,double dFXVolatility,
                        DKMaille<double> &dMaximumDrift,DKMaille<double> &dMinimumDrift,
                        DKMaille<int> *iAddToTop,DKMaille<int> *iAddToBottom);

void CalculateFKFunction(int iActualIndex,
                         DKMaille<double> &SpotFXTable,
                         DKMaille<double> *FeynmanKacFunction,
                         double dStrike,
                         double dCallPut);

void CalculateForwardFKFunction(int iActualIndex,
                                DKMaille<double> &SpotFXTable,
                                DKMaille<double> *FeynmanKacFunction,
                                double dStrike,
                                double dLongShort);

void CreateThreeFactorFXLattice(int iTimeLine,
                                DKMaille<double> &t_sys,DKMaille<double> &LowestGBase,
                                DKMaille<double> &LowestGForeign,double forateBase,double forateForeign,
                                DKMaille<double> &LocalVolBase,DKMaille<double> &LocalVolForeign,
                                double smileBase,double smileForeign,DKMaille<double> &CurrentSpotFXTable,
                                DKMaille<double> &ForwardSpotFXTable,
                                DKMaille<double> *Current3FSpotFXTable,
                                DKMaille<unsigned int> *CurrentSpotFXConnections,
                                DKMaille<double> *CurrentSpotFXUpProbs,
                                double dFXVolatility);

void CreateThreeFactorFXLatticeNew(int iTimeLine,int iFXSize,int iFXSizeNext,DKMaille<double> &t_sys,DKMaille<double> &LowestGBase,
                                  DKMaille<double> &LowestGForeign,double forateBase,double forateForeign,
                                  DKMaille<double> &LocalVolBase,DKMaille<double> &LocalVolForeign,
                                  double smileBase,double smileForeign,DKMaille<double> &CurrentSpotFXTableNew,
                                  DKMaille<double> &ForwardFXTableNew,
                                  DKMaille<double> *Current3FSpotFXTableNew,
                                  DKMaille<unsigned int> *CurrentSpotFXConnectionsNew,
                                  DKMaille<double> *CurrentSpotFXUpProbsNew,
                                  double dFXVolatility);

void CreateStaticFXLattice_TD(int no_of_steps, int iActualIndex,
                              DKMaille<double> &t_sys,double dSpotFX,
                              DKMaille<double> dFXVolatility,DKMaille<double> *SpotFXTable);

void CalculateOptionFKFunction(int iActualIndex,
                               DKMaille<double> &SpotFXTable,
                               DKMaille<double> *FeynmanKacFunction,
                               double dCallPut,
                               double dStrike,
                               double dCapStrike,
                               double dSmoothing);


int iLinearFromCube(int j1,int j2,int j3,int size);
int iLinearFromCubeNew(int j1,int j2,int j3,int size1,int size2,int size3);

void CreateThreeFactorCashFlowExchange(int iTimeLine,
                                       double dFundingAccrualPeriod,
																			 double dFXCouponAccrualPeriod,
                                       DKMaille<double> &LIBOR,
                                       DKMaille<double> &DomesticDiscount,
                                       double dFunding,
                                       double dLongShort,
                                       double dFXFixing,
                                       DKMaille<double> &DigitalPayoff,
                                       DKMaille<double> *Structure);

void CalculateBinaryFKFunction(int iActualIndex,
                               DKMaille<double> &SpotFXTable,
                               DKMaille<double> *FeynmanKacFunction,
                               double dStrike,
                               double dPayoffAbove,
                               double dPayoffBelow,
                               double dSmoothing);

void CreateThreeFactorFundingFunction(int iLastResetIndex,
                                      int iPaymentIndex,
                                      DKMaille<double> &t_sys,
                                      DKMaille<double> &LowestGBase,
                                      DKMaille<double> &foratesBase,
                                      DKMaille<double> &LocalVolBase,
                                      double smileBase,
                                      DKMaille<double> *LIBOR,
                                      DKMaille<double> *DomesticDiscount);

void CreateThreeFactorFXCallFunction(int iTimeLine,DKMaille<double> &Payoff,DKMaille<double> *Call);

void CreateThreeFactorIntermediateFXCallFunction(int iTimeLine,
                                                 DKMaille<double> &BackwardsPayoff,
                                                 DKMaille<double> &BackwardsCallTree,
                                                 DKMaille<double> *IntermediateCallTree);


void CreateThreeFactorFXFunction(int iTimeLine,DKMaille<double> &FeynmanKac,DKMaille<double> *Payoff);
void CreateThreeFactorFXFunctionNew(int iTimeLine,int iFXAdditional,DKMaille<double> &FeynmanKac,DKMaille<double> *Payoff);

void RollBackwardsFXFunction(int iTimeLine,
                             DKMaille<double> &Payoff,
                             DKMaille<double> &dDiscountsBase,
                             DKMaille<unsigned int> &PreviousSpotFXConnections,
                             DKMaille<double> &PreviousSpotFXUpProbs,
                             DKMaille<double> *BackwardsPayoff,
                             double dCorrelation12,
                             double dCorrelation23,
                             double dCorrelation13);

void RollBackwardsFXFunctionNew(int iTimeLine,
                                int iFXSize,
                                int iFXSizeNext,
                                DKMaille<double> &Payoff,
                                DKMaille<double> &dDiscountsBase,
                                DKMaille<unsigned int> &SpotFXConnections,
                                DKMaille<double> &SpotFXUpProbs,
                                DKMaille<double> *BackwardsPayoff,
                                double dCorrelation12,
                                double dCorrelation23,
                                double dCorrelation13);



void CopyOver(DKMaille<double> &from, DKMaille<double> &to);
void Multiply(double dCoupon,DKMaille<double> *Payoff);
void AddTo(DKMaille<double> &NewPayoff, DKMaille<double> &BackwardsPayoff);

DKMaille<double> TreeDK2F_SWAPTION_BCIF_QUANTO(int strikesize,
                                               DKMaille<double> &swapmaturityset,
                                               DKMaille<double> &optionexpirationset,
                                               int today,
                                               DKMaille<double> &strikeset,
                                               double dBasis,
                                               double notional,
                                               DKMaille<double> &funding,
                                               DKMaille<double> &t,
                                               int no_of_steps,
                                               DKMaille<double> &payment_dates,
                                               DKMaille<double> &rateBase,
                                               DKMaille<double> &volatilityBase,
                                               double smile0Base,
                                               DKMaille<double> &foratesBase,
                                               DKMaille<double> &rateForeign,
                                               DKMaille<double> &volatilityForeign,
                                               double smile0Foreign,
                                               DKMaille<double> &foratesForeign,
                                               double dCorrelationBaseForeign,
                                               double dNoticePeriod);


int orderingintree(DKMaille<double> time_array, int array_size, double target_time, int *nearest_node);

int calage(DKMaille<double> time_array,int array_size,int *nearest_node2,int *nearest_node1,double target_time1,double target_time2);



void AmericanFXFunction(int iTimeLine,
                        DKMaille<double> &BackwardsPayoff,
                        DKMaille<double> &SpotFXTable,
											  double dStrike,
											  double dCallPut,
                        DKMaille<double> *NewPayoff);










