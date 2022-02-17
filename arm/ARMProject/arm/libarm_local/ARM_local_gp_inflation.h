/*! \file ARM_local_infcurv.h
 *
 *  \brief file for the local interface of inflation objects
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2003
 */


#ifndef ARMLOCAL_INFCURV_H
#define ARMLOCAL_INFCURV_H

//// forward declaration
class ARM_Date;

#include "firstToBeIncluded.h"
#include "ARM_result.h"
#include "ARM_local_gp_genericaddin.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>

#include "gpbase/curvetypedef.h"
#include <string>
#include <ARM\local_xlarm\ARM_local_interglob.h>

enum instrument {CPI, ZCRate};

typedef struct _CPI_REF
{
    double CPIDate;
    double ZCDate;
    double CPIValue;
} 
ARM_CPI_REF;

/// function to convert gap from date or gap
int GetGapOrJulianDate( int gapInput );

/// create a curve
extern long ARMLOCAL_InfCurv_Create( 
	double asOfDate,
	const CCString&	indexName,
	double CPIIndexValue,
	double CPIIndexDate,
	const VECTOR<CCString>&	maturities,
	const VECTOR<double>&	values,
	long MonthlyInterpType,
	long DailyInterpType,
	long DCFMonthly,
	long DCFDaily,
	long ExtrapolType,
	long resetManagerId,
	long seasonManagerId,
	ARM_result&	result,		// here to allow default values after
	long objId				= ARM_NULL_OBJECT_ID );


/// central call for both CPI and ZCRate
extern long ARMLOCAL_InfCurv_Interp(
	long					idCurve,
	double					CPIDate, 
	ARM_result&				result,
	const CCString&			DCFLag,
	long					DailyInterpType,
	const CCString&			ResetLag,
	double					weight,
	instrument				what );


/// to create an inflation index
extern long ARMLOCAL_InfIdx_Create(
	const CCString& indexName,
	const CCString& ResetLag,
	const CCString& DCFLag,
	long ccyId,long infCurveId,
	ARM_result&	result,
	long objId	= ARM_NULL_OBJECT_ID );


/// to create an inflation leg generic
extern long ARMLOCAL_InfLegwDateStrip_Create(
	const CCString& indexName,
	long rcvOrPay,
	int interpType,
	double multiple,
    double constant,
	long finalNotionalType,
	long numStripDateId,
	long denomStripDateId,
	ARM_result&	result,
	long objId 	= ARM_NULL_OBJECT_ID );


/// to create an inflation leg with all inputs
extern long ARMLOCAL_InfLegAllInputs_Create(
	double startDate, 
	double endDate,
	const CCString& indexName,
	int swapType,
	int rcvOrPay,
	int interpType,
	double multiple,
	double CoMultiple,
	double constant,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetNumGap,
	int resetDenomGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	int adjFirstDate,
	int finalNotionalType,
	double firstReset,
	ARM_result&	result,
	long objId = - 1);


/*****************************************************************************************************/
class ARM_Corridor_InfIr_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_Corridor_InfIr_CreateFunctor() {}

	virtual char* ClassName() const { return "LCINF"; }

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_HybridInfIrLeg_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_HybridInfIrLeg_CreateFunctor() {}

	virtual char* ClassName() const { return "LHINF"; }

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_HybridInfIrPayOff_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_HybridInfIrPayOff_CreateFunctor() {}

	virtual char* ClassName() const { return "LHIPO"; }

	virtual	long operator()( ARM_result& result, long objId );
};


class ARM_HybridInfIrModel_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_HybridInfIrModel_CreateFunctor() {}

	virtual char* ClassName() const { return "LHIMO"; }

	virtual	long operator()( ARM_result& result, long objId );
};

extern long ARMLOCAL_HybridInfIrMkt_Create(		const ARM_Date&			asOf,
												const vector<string>&	keys,
												const vector<long>&		mods,
												ARM_result&				result,	
												long					objId = -1);

extern long ARMLOCAL_HybridInfIr_Load	(		const long & ins, 
												const long & mkt,
												const long & mod,
												const long & pay,
												ARM_result & result,	
												long	objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_HybridInfIrPayOff_Create(	const long   & cstCpnCoef,
												const long   & cstOptCoef,
												
												const string & mainCpnName,
												const long	 & mainCpnCoef,
												const string & mainOptName, 
												const long	 & mainOptCoef, 
												
												const string & subCpnName,
												const long	 & subCpnCoef,
												const string & subOptName, 
												const long   & subOptCoef,
												
												const string & supCpnName,
												const long	 & supCpnCoef, 
												const string & supOptName, 
												const long   & supOptCoef,
												ARM_result	 & result,	
												long		   objId= -1);

extern long ARMLOCAL_InfSpreadCap_Create(		const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strike, 
												const long	 & notional, 
												ARM_result	 & result,
												long		   objId= -1);

extern long ARMLOCAL_InfSpreadDigital_Create(	const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strike, 
												const long	 & notional, 
												ARM_result	 & result,	
												long		   objId= -1);

extern long ARMLOCAL_InfDoubleDigital_Create(	const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const long	 & mainStrike,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & subStrike, 
												const long	 & notional, 
												const long	 & spread,
												ARM_result	 & result,	
												long		   objId =-1);

extern long ARMLOCAL_InfCorridor_Create(		const string & mainIndex,
												const string & mainType,
												const long	 & mainLeverage,
												const string & subIndex,
												const string & subType,
												const long	 & subLeverage,
												const long	 & strikeInf, 
												const long	 & strikeSup,
												const long	 & notional, 
												ARM_result	 & result,	
												long		   objId);


extern long ARMLOCAL_Inf_GetPrice(				const long				& pricerId,
												const string			& key,
												ARM::ARM_GramFctorArg	& argResult,
												ARM_result				& result);

extern long ARMLOCAL_Inf_GetSchedule(			const long				& legId,
												const string			& key,
												ARM::ARM_GramFctorArg	& argResult,
												ARM_result				& result);


/*****************************************************************************************************/



extern long ARMLOCAL_infYCmod(
	long idZeroCurve, 
	long idDiscCurve,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_fixZC(
	double startDate, 
	double endDate,
	double fixRate,
	int rcvOrPay,
	int dayCount,
	int intRule,
	int stubRule,
	int payGap,
	const CCString& payCalendar,
	long ccyId,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );
    

extern long ARMLOCAL_GetData(
	const VECTOR<CCString>& vec,
	int nbcolumns,
	int nbrows,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_SparseVolCube_CreateNFill(
	double asOfDate,
	double lastKnownDate,
	const CCString& indexName,
	long dim1Type,
	const VECTOR<CCString>& dim1Value,
	long dim2Type,
	double dim2Value,
	const VECTOR<double>& strikes,
	const VECTOR<double>& vols,
	long strikeType,
	long volType,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_SparseVolCube_Fill(
	long dim1Type,
	const VECTOR<CCString>& dim1,
	long dim2Type,
	double dim2Value,
	const VECTOR<double>& strikes,
	const VECTOR<double>& vols,
	long sparseVolCubeId,
	ARM_result& result,
	long objId );


extern long ARMLOCAL_VolCubeFromSparseVolCube(
	long sparseVolCubeId,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_infBSMod(	double asOfDateDble,
	long zeroCurveId, 
	long fwdInfCurveId, 
	long infCapVolCurveId,
	long correlManagerId,
	long IRBSModelId,
	long infSwoptCurveId,
	long IRSwoptCurveId,
	ARM_result& result,
	long objId				= ARM_NULL_OBJECT_ID);



extern long ARMLOCAL_infBSSmiledModel(	double		C_AsOfDate,
										long		C_DiscountCurvId,
										long		C_InfFwdCurvId,
										long		C_VolSigmaId,
										long		C_VolNuId,
										long		C_VolRhoId,
										long		C_VolBetaId,
										long		C_VolAtmIrId,
										long		C_CorrelId,
										long		C_CorrelAdjId,
										ARM_result& result, 
										long		objId = ARM_NULL_OBJECT_ID);




extern long ARMLOCAL_InfCapFloor_Create(
	double startDate,
	double endDate,
	const CCString& indexName,
	int capOrFloor,
	double strike,
	double leverage,
	double spread,
	int swapType,
	int rcvOrPay,
	int interpType,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	int adjFirstDate,
	double firstReset,
	long ccyId,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_CorrelMat_Create(
	double asOfDate,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>&	Z,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_CorrelManager_Create(
	const CCString& mktTag, 
	const CCString& intraMktTag,
	double	asOfDate,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>& Z,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_CorrelManager_Fill(
	const CCString& mktTag,
	const CCString& intraMktTag,
	const VECTOR<CCString>& X,
	const VECTOR<CCString>& Y,
	const VECTOR<double>& Z,
	long correlManagerId,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_CorrelManager_FillFromMat(
	const CCString& mktTag,
	const CCString& intraMktTag,
	long CorrelMatId,
	long correlManagerId,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );



extern long ARMLOCAL_CorrelManager_CreateFromMat(
	const CCString& mktTag, 
	const CCString& intraMktTag,
	long CorrelMatId,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


extern long	ARMLOCAL_CreateGenCorrelManager (VECTOR<CCString>& mktTags,
											 VECTOR<CCString>& intraMktTags,
											 vector<long>& correlVolIds, 
											 ARM_result& result, 
											 long objId = ARM_NULL_OBJECT_ID);



extern long ARMLOCAL_ComputeCorrelFromCorrelManager(
	const CCString& tmpMktTag,
	const CCString& tmpIntraMktTag,
	double x,
	double y,
	long correlManagerId,
	ARM_result& result );


extern long ARMLOCAL_Vol_to_Cor(
		const VECTOR<double>& ZCVolV, 
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& MaturityV,
		VECTOR<double>& CorV,
		ARM_result& result );

extern long ARMLOCAL_YtYCor_to_ZC(
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& CorV, 
		const VECTOR<double>& MaturityV,
		VECTOR<double>& ZCVolV,
		ARM_result& result );

extern long ARMLOCAL_ZCCor_to_YtY(
		const VECTOR<double>& ZCVolV, 
		const VECTOR<double>& CorV, 
		const VECTOR<double>& MaturityV,
		VECTOR<double>& YtYVolV,
		ARM_result& result );

extern long ARMLOCAL_Bounds(
		const VECTOR<double>& ZCVolV,
		const VECTOR<double>& YtYVolV,
		const VECTOR<double>& CorV,
		const VECTOR<double>& MaturityV,
		VECTOR<double>& UBoundV,
		VECTOR<double>& LBoundV,
		const string& choice, 
		string& TBound, 
		ARM_result& result );

extern long ARMLOCAL_HmgVol_to_Cor(
		const VECTOR<double>& ZCVolV, 
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& MaturityV,
		VECTOR<double>& CorV,
		const int length,
		ARM_result& result );

extern long ARMLOCAL_HmgZCCor_to_YtY( 
		const VECTOR<double>& pZCVol, 
		const VECTOR<double>& pCor, 
		const VECTOR<double>& pMaturity,
		const int length,
		VECTOR<double>& YtYVol,
		ARM_result& result );

extern long ARMLOCAL_HmgYtYCor_to_ZC(
		const VECTOR<double>& YtYVolV, 
		const VECTOR<double>& CorV, 
		const VECTOR<double>& MaturityV,
		const int length,
		VECTOR<double>& ZCVolV,
		ARM_result& result );


extern long ARMLOCAL_VolYoY_to_VolSwp(
		const VECTOR<double>&  DFactor, 
		const VECTOR<double>&  FwdCPI, 
		const VECTOR<double>&  Vol_DF,
		const VECTOR<double>&  Vol_YoY, 
		const VECTOR<double>&  AvgCor, 
		const VECTOR<double>&  Dates, 
		const VECTOR<double>&  Tenors, 
		const double SwpRate, 
		double &Vol_Swp,
		ARM_result& result );

extern long ARMLOCAL_SeasonalityManager_Create(
			const VECTOR<CCString>&	monthsList,
			const VECTOR<double>&	seasonSpreadList,
			const VECTOR<double>&	seasonHorizonList,
			long SeasonalityCorrectionType,
			ARM_result&				result,
			long					objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_GetInfZcFromSummit_Create(
			const CCString&	index,
			const CCString&	ccy,
			const CCString&	cvname,
			double date,
			long seasonAdjId,
			long seasonAdjModeId,
			ARM_result&				result,
			long					objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_LIVRETACURVE (double asOfDateDble,
								   long infCurvId,
								   long euribCurvId,
								   long flagArrondi,
								   long infresetManagerId,
								   long fixingEuribId,
								   long fixingLivretAId,
								   long monthForAugustId,
								   long monthForFebruaryId,
								   ARM_result& result,
								   long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_LIVRETACURVEGETRATEDATE (long livretACurvId, 
											  double dateDble, 
											  ARM_result& result);


extern long ARMLOCAL_InfSwoVolCurveFromModel_Create(
			const ARM_Date& asOfDate,
			long InfIRModelId,
			const VECTOR<double>& tenors,
			const VECTOR<double>& expiries,
			long computationMethod,
			ARM_result&				result,
			long					objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_InfSwoVolCubeFromModel_Create(
			const ARM_Date& asOfDate,
			long InfIRModelId,
			const VECTOR<double>& tenors,
			const VECTOR<double>& expiries,
			const VECTOR<double>& smiledTenors,
			const VECTOR<double>& strikes,
			long computationMethod,
			ARM_result&				result,
			long					objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_InfOATSwoVolCurveFromModel_Create(
			const ARM_Date&			asOfDate,
			long					InfIRModelId,
			const VECTOR<double>&	tenors,
			const VECTOR<double>&	expiries,
			double coupon,
			long					choice,
			long					computationMethod,
			ARM_result&				result,
			long					objId = ARM_NULL_OBJECT_ID );


long ARMLOCAL_infMultiBSMod_Create(	
	const VECTOR<long>& infMultiBSModIdVec,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


long ARMLOCAL_InfCurv_SetResetManager(
	const long& infCurvId,
	const long& resetManagerId,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


long ARMLOCAL_InfCurv_SetSeasonalityManager(
	const long&	infCurvId,
	const long& seasonalityManagerId,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID );


extern long ARMLOCAL_GetSeasonMgrFromSummit (const CCString& index,
											 const CCString& ccy,
											 const CCString& cvname,
											 double asOf,
											 long modeId,
											 ARM_result& result,
											 long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_GetResetMgrFromSummit (double C_asof,
											const CCString& C_index,
											const CCString& C_source,
											const CCString& C_ccy,
											long isInflatIndexId,
											const CCString& C_term,
											ARM_result& result,
											long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_GetReset (long resetMgrId,
							   double date,
							   ARM_result& result);

extern long ARMLOCAL_GP_INFCAPFLOOR (long swapId,
									 long CF,
									 double strike,
									 long strikeId,
									 ARM_result&	result, 
									 long        objId);

extern long ARMLOCAL_GP_INFDIGITAL (long payLegId,
									long digitLegId,
									long payOffType,
									double barrier,
									long barrierId,
									long CFId,
									long RecOrPayId,
									ARM_result&	result, 
									long        objId);

extern long ARMLOCAL_GP_INFCALLSPREAD (	long			swapId,
										long			CF,
										double			strike,
										long			strikeId,
										ARM_result&		result, 
										long			objId);

extern long ARMLOCAL_InfEqHwSV_Laplace(		const long&				modelId,
											const double&			evalTime,
											const double&			startTime,
											const double&			endTime,
											const double&			xt,
											const double&			vt,
											const double&			k_real,
											const double&			XL_k_imag,
											const bool&				XL_isReal,
											ARM_result&				result );

extern long ARMLOCAL_InfEqHwSV_Density(		const long&				modelId,
											const double&			evalTime,
											const double&			startTime,
											const double&			endTime,
											const double&			xt,
											const double&			vt,
											const double&			x,
											const double&			period,
											const double&			frequency,
											ARM_result&				result );

extern long ARMLOCAL_Inf_GetAdjCorrel(		const long		& infBsSmiledId,
											ARM_result		& result,	
											long			objId =-1);




#endif	// ARMLOCAL_INFCURV_H


