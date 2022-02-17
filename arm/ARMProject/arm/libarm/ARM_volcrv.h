#ifndef ARM_VOLCRV_H
#define ARM_VOLCRV_H



#include "ARM_result.h"



extern long ARM_volcurv(const VECTOR<double>& matu,
						const VECTOR<double>& strikes,
				        const VECTOR<double>& vols,
						double date,
				        long strikeType,
						long volType,
				        ARM_result& result,
						long objId = -1);

extern long ARM_GetVolFromSummit(const CCString& index,
								 const CCString& currency, 
								 const CCString& cvName,
								 double date, 
								 const CCString& vtype,
								 const CCString& matuIndex,
								 ARM_result& result,
								 long objId = -1);

extern long ARM_GetVolCubeFromSummit(const CCString& index, 
									 const CCString& currency, 
 						             const CCString& cvName,
									 double date, 
						             const CCString& vtype,
							         VECTOR<CCString>& tenors,
							         const CCString& smileOrNot,
						             ARM_result& result, 
									 long objId = -1);

extern long ARM_ComputeVolatility(long idCurve,
								  double matu, 
								  double strike,
								  double tenor,
								  ARM_result& result);

extern long ARM_volflat(double vol,
						double date,
						ARM_result& result,
						long objId = -1);


extern long ARM_ARM_GetSummitVolMatrix(const CCString& index,
									   const CCString& currency, 
								       const CCString& cvName,
									   double date,
								       const CCString& vtype,
									   long summitFormat,
								       ARM_result& result);

extern long ARM_VolCube(long ATMVolId,
						const VECTOR<long>& volCurveIds,
						const VECTOR<double>& tenors,
						ARM_result& result,
						long objId = -1);

extern long ARM_GetFXVolFromSummit(const CCString& ccy1,
								   const CCString& ccy2,
								   double date,
								   const CCString& cvName,
								   ARM_result& result,
								   long objId = -1);

extern long ARM_ARM_BumpVolatility(long VolId,
								   double valueToBump,
								   long nthLine,
								   long nthCol,
								   long cumulId,
								   ARM_result& result,
								   long objId = -1);
#endif	// ARM_VOLCRV_H

/*---- End Of File ----*/
// EOF %M%