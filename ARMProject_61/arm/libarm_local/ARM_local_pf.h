#ifndef ARM_LOCAL_PF_H
#define ARM_LOCAL_PF_H

#include "ARM_result.h"



extern long ARMLOCAL_PF  (const VECTOR<long>& insts,
						  const VECTOR<double>& coeffs,
						  const VECTOR<double>& marketPrices,
						  const VECTOR<double>& precisions,
						  ARM_result& result,
						  long objId = -1);

extern long ARMLOCAL_PFFILL_Create(const VECTOR<long >& assetsIdVec,
		 VECTOR<double >& weightsVec,
         VECTOR<double >& mktpricesVec,
         VECTOR<double >& vegasVec,
		long portfolioId,
         ARM_result&	result, 
         long        objId );

extern long ARMLOCAL_PF_Merge(
		 VECTOR<long> pfIds,
         ARM_result&	result,
         long        objId );

extern long ARMLOCAL_PFGYCSIGVARFIT (long pfId,
									 long curveId,
									 const VECTOR<double>& matCurve,
									 double in_min_meanRev,
									 double in_max_meanRev,
									 double in_min_vol,
									 double in_max_vol,
									 double in_precision_meanRev,
									 double in_precision_vol,
									 long nbMaxIter,
									 ARM_result& result,
									 long objId = -1);

extern long ARMLOCAL_PFLOGDECVOLFIT(long pfId,
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

extern long ARMLOCAL_PFGYCSIGVARPENALFIT(long pfId,
										 long curveId,
										 const VECTOR<double>& matCurve,
										 double start_meanRev,
										 const VECTOR<double>& start_vol_vect,
										 const VECTOR<double>& penal_vect,
										 const VECTOR<double>& coeff_vect,
										 double theAccuracy,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_PFINSTLOGDECVOLFIT(long pfId,
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

extern long ARMLOCAL_PFMODFIT (const CCString& modName,
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


extern long CC_LOCAL_PFModEstimate(char* mod_name, int pfId, char* settl, int zcId,
                            int size, double* initVect, int* flagEstim,
                            int step, char* horizon, int NAG_ALGO, int MOD_ID, long & retour);

extern long CC_LOCAL_SetPFModEstimated(int modId, char* mod_name, int pfId, 
                                  char* settl, int zcId,
                                  int size, double* initVect, int* flagEstim,
                                  int step, char* horizon, int NAG_ALGO, long & retour);

extern long CC_LOCAL_PFModEstimateAlgo1(char* mod_name, int pfId, char* settl,
                                   int zcId, 
                                   int size, double* initVect, int* flagEstim,
                                   int step, char* horizon, long & retour);

extern long CC_LOCAL_SetPFModEstimatedAlgo1(int modId, char* mod_name, int pfId, 
        char* settl, int zcId, int size, double* initVect, int* flagEstim,
        int step, char* horizon, long & retour);

extern long ARMLOCAL_CalibrationPF(long zcId,
							long volId,
							long secId,
							long modeType,
							long portfolioType,
							double shift,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_MRSCalibrationPF(long zcId,
							long volId,
							long secId,
							long portfolioId,
							int freq,
							bool ATMStrikes,
							ARM_result& result,
							long objId = -1);

extern string SecurityTypeToClass(long portfolioId, long assetIdIndex);

extern long ARMLOCAL_GetAssetFromPF(long portfolioId, 
                                    long assetIdIndex,
                                    ARM_result& result,
                                    long    objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_GetSizeOfPF(long portfolioId, 
                                    long& sizeofPf,
                                    ARM_result& result);

extern long ARMLOCAL_ComputeAll(long pfId, 
								long modelId,
								ARM_result& result);


#endif /* ARM_LOCAL_PF_H */