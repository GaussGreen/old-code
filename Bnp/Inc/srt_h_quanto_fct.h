#ifndef QUANTO_H_
#define QUANTO_H_

Err srt_f_quanto_fra(long FraStartDate, long FraEndDate, String refratecode,
                     char *szFXUndName, double *adjfra);

Err srt_f_quanto_swap(long SwapStartDate, long SwapEndDate,
                      SrtCompounding FixedFreq, SrtBasisCode FixedBasis,
                      String FloatRefRate, char *FxUndName, double *AdjSwap);

Err srt_SwapFraWeigths(Date SwapStartDate, Date SwapEndDate,
                       SrtCompounding FixedFreq, SrtBasisCode FixedBasis,
                       SrtCompounding FloatFreq, SrtBasisCode FloatBasis,
                       String YieldCurveName, double **WeightsVector,
                       long *WeightsVectorSize);

Err srt_Weigths(Date SwapStartDate, Date SwapEndDate, String FixedFreqStr,
                String FixedBasisStr, String FloatFreqStr, String FloatBasisStr,
                String YieldCurveName, double **WeightsVector,
                long *WeightsVectorSize);

Err srt_FraSwapWeights(long StartDate, long EndDate, SrtCompounding FixedFreq,
                       SrtBasisCode FixedBasis, SrtCompounding FloatFreq,
                       SrtBasisCode FloatBasis, String YieldCurveName,
                       /*OUTPUT*/
                       double ***FraSwapWeightsMatrix,
                       long *FraSwapWeightsSize);

#endif