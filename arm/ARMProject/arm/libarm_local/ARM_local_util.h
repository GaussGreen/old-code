#ifndef ARM_LOCAL_UTIL_H
#define ARM_LOCAL_UTIL_H

#ifndef MATU_LINES_SIZE
	#define MATU_LINES_SIZE 13
	#define MATU_COLS_SIZE 13
#endif

struct lCurve
{
	CCString Currency;
	CCString Index;
	CCString CvName;
	long id;
};

struct lForex
{
	CCString domCcy;
	CCString forCcy;

	double fxSpot;
};


extern long transformInMonth(const CCString& pStr);

extern double calcMoy(ARM_Vector* Data1);

extern long ARMLOCAL_GetCorrelInst(double date1,
								   double date2,
								   const CCString& ccy1,
								   const CCString& index1,
								   const CCString& expiry1,
								   const CCString& tenor1,
								   const CCString& ccy2,
								   const CCString& index2,
								   const CCString& expiry2,
								   const CCString& tenor2,
								   const CCString& curve1_ccy1,
								   const CCString& curve2_ccy1,	
								   const CCString& nbmonths_curve1_ccy1,
								   const CCString& curve1_ccy2,
								   const CCString& curve2_ccy2,	
								   const CCString& nbmonths_curve1_ccy2,
								   long typeId,
								   double lambda,
								   double precision,
								   const CCString& ccy,
								   ARM_result& result);

extern long ARMLOCAL_GetMoyCorrel(double date1,
								  double date2,
								  const CCString& ccy1,
								  const CCString& index1,
								  const CCString& expiry1,
								  const CCString& tenor1,
								  const CCString& ccy2,
								  const CCString& index2,
								  const CCString& expiry2,
								  const CCString& tenor2,
								  const CCString& curve1_ccy1,
								  const CCString& curve2_ccy1,	
								  const CCString& nbmonths_curve1_ccy1,
								  const CCString& curve1_ccy2,
								  const CCString& curve2_ccy2,	
								  const CCString& nbmonths_curve1_ccy2,
								  long typeId,
								  double lambda,
								  double precision,
								  const CCString& ccy,
								  ARM_result& result);

extern long ARMLOCAL_GetCorrelQuanto(double date1,
									 double date2,
									 const CCString& ccy,
									 const CCString& index,
									 const CCString& expiry,
									 const CCString& tenor,
									 const CCString& cvname1,
									 const CCString& cvname2,
									 const CCString& switchinmonth,
									 const CCString& domccy,
									 const CCString& domindex,
									 const CCString& forccy,
									 const CCString& forindex,
									 long typeId,
									 double lambda,
									 double precision,
									 const CCString& calccy,
									 long fwdOrNot,
									 ARM_result& result);

extern long ARMLOCAL_GetMoyCorrelQuanto(double date1,
										double date2,
										const CCString& ccy,
										const CCString& index,
										const CCString& expiry,
										const CCString& tenor,
										const CCString& cvname1,
										const CCString& cvname2,
										const CCString& switchinmonth,
										const CCString& domccy,
										const CCString& domindex,
										const CCString& forccy,
										const CCString& forindex,
										long typeId,
										double lambda,
										double precision,
										const CCString& calccy,
										long fwdOrNot,
										ARM_result& result);

extern long ARMLOCAL_GetDealsFromSummitFilter (const CCString& filter,
											   VECTOR<CCString>& listDeals,
											   VECTOR<CCString>& listTypes,
											   ARM_result& result);

extern long ARMLOCAL_GetAsOfVolOrRate(double date1,
									  const CCString& ccy,
									  const CCString& index,
									  const CCString& cvName,
									  const CCString& expiry,
									  const CCString& matu,
									  long yieldOrValId,
									  long calcModId,
									  const CCString& volType,
									  ARM_result& result);

extern long ARMLOCAL_GetFwdRatesMatrix(double asOf,
									   const CCString& ccy,
									   const CCString& index,
									   const CCString& cvName,
									   ARM_result& result);

extern long ARMLOCAL_ARM_GetInfo (long secId,
								  const CCString& secClass,
								  const CCString& type,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_ARM_GetModelFactor(double AsOfDate,
										const CCString& model,
										const CCString& type,
										const CCString& factorName,
										const CCString& ccy,
										const CCString& index,
										const CCString& cvName,
										long calcModId,
										ARM_result& result,
										long objId = -1);

#endif