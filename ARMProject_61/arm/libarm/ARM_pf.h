#ifndef ARM_PF_H
#define ARM_PF_H

#include "ARM_result.h"



extern long ARM_PF  (const VECTOR<long>& insts,
                     const VECTOR<double>& coeffs,
                     const VECTOR<double>& marketPrices,
                     const VECTOR<double>& precisions,
					 ARM_result& result,
					 long objId = -1);

extern long ARM_MKTPRICE (long idPF,
						  ARM_result& result);

extern long ARM_PFSPLESTIMATE (long pfId,
                               long splineId,
                               double settlement,
                               ARM_result& result,
							   long objId = -1);

extern long ARM_PFVSKESTIMATE (long pfId,
                               double settlement,
                               double a,
                               ARM_result& result,
							   long objId = -1);

extern long ARM_PFGYCFIT (long pfId,
                          double curveId,
                          double settlement,
                          ARM_result& result,
						  long objId = -1);

extern long ARM_PFTHEOPRICE (long pfId,
                             long modelId,
                             double settlement,
                             ARM_result& result);

extern long ARM_PFGYCSIGVARFIT (long pfId,
								long curveId,
                                const VECTOR<double>& matCurve,
                                double in_min_meanRev,
                                double in_max_meanRev,
                                double in_min_vol,
                                double in_max_vol,
                                double in_precision_meanRev,
                                double in_precision_vol,
                                long nbMaxIter,
                                ARM_result& result, long objId = -1);

extern long ARM_PFMODFIT (const CCString& modName,
                          long pfId,
                          double settlement,
                          long zcId,
                          VECTOR<double>& vList,
                          VECTOR<double>& fList,
                          long step,
                          double horizon,
                          long nag_algo,
                          ARM_result& result,
                          long objId = -1);



extern long ARM_PFGYCSIGVARGLOBFIT(long pfId, 
                                   long curveId,
                                   const VECTOR<double>& matCurve,
                                   double start_meanRev,
                                   const VECTOR<double>& start_vol,
                                   ARM_result& result,
                                   long objId = -1);


extern long ARM_PFGYCTWOFACTSIGVARFIT(long pfId,
                                      long curveId,
                                      const VECTOR<double>& matCurve,
                                      double in_a,
                                      double in_b,
                                      double in_sigmaRatio,
                                      double in_rho,
                                      double in_min_vol,
                                      double in_max_vol,
                                      long   in_precision_vol,
                                      long   nbMaxIter,
                                      ARM_result& result, long objId = -1);

extern long ARM_PFHW2FPENALFIT(long pfId,
							   long curveId,
                               const VECTOR<double>& initVect,
                               const VECTOR<double>& penal_vect,
                               const VECTOR<double>& coeff_vect,
                               double theAccuracy,
                               ARM_result& result,
                               long objId = -1);

extern long ARM_PFGYCSIGVARPENALFIT(long pfId,
									long curveId,
                                    const VECTOR<double>& matCurve,
                                    double start_meanRev,
                                    const VECTOR<double>& start_vol_vect,
                                    const VECTOR<double>& penal_vect,
                                    const VECTOR<double>& coeff_vect,
                                    double theAccuracy,
                                    ARM_result& result,
                                    long objId = -1);


extern long ARM_PFLOGDECVOLFIT(long pfId, 
                               long curveId,
                               const VECTOR<double>& proba,
                               double theAccuracy,
                               long shapeType,
                               double decay,
                               double slope,
                               double asymptote,
                               const VECTOR<double>& matCurve,
                               const VECTOR<double>& volinit_vect,
                               const VECTOR<double>& coeff_vect,
                               ARM_result& result,
                               long objId = -1);

extern long ARM_INSTLOGDECVOLFIT(long secId,
								 long curveId,
								 long volBSId,
								 const VECTOR<double>& proba,
								 long LiborType,
								 double strike,
								 double theAccuracy,
								 long shapeType,
								 double decay,
								 double slope,
								 double asymptote,
								 const VECTOR<double>& matCurve,
								 const VECTOR<double>& volinit_vect,
								 const VECTOR<double>& coeff_vect,
								 ARM_result& result,
								 long objId = -1);

extern long ARM_PFINSTLOGDECVOLFIT(long pfId,
                                   long secId, 
                                   long curveId,
                                   const VECTOR<double>& proba,
                                   long UsePFResetDates,
                                   double theAccuracy,
                                   long shapeType,
                                   double decay,
                                   double slope,
                                   double asymptote,
                                   long VolBSId,
                                   const VECTOR<double>& matCurve,
                                   const VECTOR<double>& volinit_vect,
                                   const VECTOR<double>& coeff_vect,
                                   ARM_result& result,
                                   long objId = -1);

#endif    // ARM_PF_H

/*---- End Of File ----*/
// EOF %M%