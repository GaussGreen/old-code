#ifndef ICM_LOCAL_LEG_H
#define ICM_LOCAL_LEG_H

#include <ICMKernel\glob\icm_enums.h>
class ARM_result;


extern long ICMLOCAL_FIXEDLEG (double		startDateIn, 
					    double		endDateIn,
						double		fixedRateIn,
						qPAYMENT_PREMIUM_LEG			AccruedOnDefaultIn,
						int			AccruedDayCountIn,
						double		LastIndexFixingIn,
						int			rcvOrPayIn,
						int			freqIn,
						int			dayCountIn,
						int			decompFreqIn,
						int			payTimingIn,
						int			intRuleIn,
						int			stubRuleIn,
						long		discountCcyIn,
						CCString	payCalNameIn,
						int			nxChangeIn,
						double		refDateIn,
						ARM_result& result,
						long		objId = -1);

extern long ICMLOCAL_FIXEDLEGGEN (	double		startDateIn, 
									double		endDateIn,
									double		fixedRateIn,
									int			AccruedDayCountIn,
									int			freqIn,
									int			dayCountIn,
									int			payTimingIn,
									int			intRuleIn,
									int			stubRuleIn,
									long		discountCcyIn,
									CCString	payCalNameIn,
									double		refDateIn,
									ARM_result& result,
									long		objId = -1);


extern long ICMLOCAL_GenLeg (  const double&	startDateIn, 
							   const double&	endDateIn,
							   const double&	fixedRate,
							   const double&	fixedNotional,
							   const long&		VarNotId,
							   const long&		VarRateId,
							   const long&		ExchangeNotId,
							   const int&		frequency,
							   const int&		daycount,
							   const int&		payTiming,
							   const int&		intRule,
							   const int&		stubRule,
							   const long&		CcyId,
							   const CCString&	payCal,
							   const double&	refDate,
							   const bool&		includematurity,
							   const int&		adjuststartDate,
							   const int&		legtype,
							   const int&		IndexId,
							   const int&		creditlag,
							   const double&	binary,
							   const CCString&	name,
							   const int&		nxchange,
							   qPAYMENT_PREMIUM_LEG iaccruedOnDef,
							   ARM_result&	result,
							   long			objId = -1);

extern long ICMLOCAL_SetVariableSpread (long secId,
									long rId,
									ARM_result& result);

extern long ICMLOCAL_SetCoupons(long secId,
								 long couponId,
								 long styleId,
								 long partId,
								 ARM_result& result);

extern long ICMLOCAL_SetRiskyProfile(long secId,
							  long typeId,
							  long dateId,
							  ARM_result& result);


extern long ICMLOCAL_SetLeg (long secId,
							 long refValId,
							 long option,	
							 ARM_result& C_result);


extern long ICMLOCAL_Index (const std::string& IndexName,
							const std::vector<std::string>& labels,
							int Basis,
							int ResetFreq,
							int PayFreq,
							ARM_Vector yearterm,
							ARM_Vector Spread,
							const std::string& ccy,
							long Method,
							int DefaultCurveId,
							int fwdRule,
							int resetTiming,
							int resetGap,
							int payTiming,
							int payGap,
							int intRule,
							qCDS_ADJ AdjCalType,
							int cm_resetWeekDay,
							int cm_resetOccur,
							long objId = -1);


extern long ICMLOCAL_CORRIDORLEG(const CCString& name,
								 double startDate,
								 double endDate,
								 double refdate,
								 double fstcpneffdate,
								 long RefValueSpreads,
								 //double notional,
								 long floatingIdx,
								 double leverageFloatIdx,
								 long creditIdx,
								 long refvalueKUP,
								 long refvalueKDW,
								 qPAYMENT_PREMIUM_LEG accondef,
								 long accdaycount,
								 long payfreq,
								 long resetfreq,
								 long daycount,
								 long paytiming,
								 long intrule,
								 long stubrule,
								 const CCString& disc_ccy,
								 const CCString& paycalname,
								 ARM_result& result,
								 long objId = -1);

extern long ICMLOCAL_CORRIDORLEG_SCHE(const string& name,
								 double Notional,
								 long recieveOrPay,
								 long RefValueSpreads,
								 long floatingIdx,
								 double leverageFloatIdx,
								 long creditIdx,
								 long refvalueKUP,
								 long refvalueKDW,
								 long scheduleId,
								 qPAYMENT_PREMIUM_LEG accondef,
								 const string& disc_ccy,
								 ARM_result& result,
								 long objId = -1);

extern long ICMLOCAL_IRLEGTOCREDITLEG(int SwapLegId,
							   int LegType,
							   int creditindexId,
							   int PricerId,
							   ARM_result& result,
							   long	objId=-1);

#endif	// ICM_LOCAL_LEG_H

/*----End Of File ----*/
// EOF %M%