/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_IRPATHDEP_CALCULATORS_H
#define ARMLOCAL_GP_IRPATHDEP_CALCULATORS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;


/****************************************************************************
						TARN CALCULATOR + TARNSB
*****************************************************************************/

long ARMLOCAL_TARNCalculator_Create(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        long cpnResetGap,
		long intRule,
        double leverage,
        long leverageId,
		double cpnMin,
        long cpnMinId,
		double cpnMax,
        long cpnMaxId,
		double lifeTimeCap,
		bool globalCapFlag,
		double lifeTimeFloor,
        double fundSpread,
        long fundSpreadId,
        long fundFreq,
        long fundDayCount,
        double nominal,
        long nominalId,
		double fees,
		long feesId,
        long fundNomId,
		const vector<double>& nbIterations,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
        const vector< string >& calibFlags,
        const vector< string >& outputFlags,
        long mktDataManagerId,
        const vector< string >& keys,
        bool isCustomResetFlag,
		long resetDatesId,
        double AsOfDate,
		ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_TARNSBCalculator_Create(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
		double coupon0,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        long cpnResetGap,
		long intRule,
		string intSubRule,
        double leverage,
        long leverageId,
		double cpnMin,
		long cpnMinId,
		double cpnMax,
		long cpnMaxId,
		double levPrev,
        long levPrevId,
		long reverse,
		double lifeTimeCap,
		bool globalCapFlag,
		double lifeTimeFloor,
        double fundSpread,
        long fundSpreadId,
        long fundFreq,
        long fundDayCount,
        double nominal,
        long nominalId,
		double fees,
        long feesId,
		double fundNom,
		long fundNomId,
		const vector< double >& nbIterations,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
        const vector< string >& calibFlags,
        const vector< string >& outputFlags,
        long mktDataManagerId,
        const vector< string >& keys,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_TARN_Get(
        long tarnId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_TARN_Set(
        long tarnId,
        long dataId,
        const string& setType,
		const vector<string>& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

long ARMLOCAL_TARNCalculator_GetPricingData(
		long tarnCalculatorId,
		const string& key,
		ARM_GramFctorArg& argResult,
		ARM_result& result );

extern long ARMLOCAL_TARN_SetProductToPrice(
	long tarnId,
	vector<string>& prodToPrice,
	ARM_result&	result,
	long objId);


/****************************************************************************
						MATURITY CAP CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_MaturityCapCalculator_Create(
        double startDate,
        double endDate,
		double underlyingEndDate,
        long longShort,
        long capFloor,
        long resetFreq,
		long payFreq,
        const string& indexTerm,
        long dayCount,
		long intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		long maturityCapMode,
		double coeff,
		double amortizing,
		long amortizingId,
		long resetGap,
        const string& resetCal,
        const string& payCal,
		long calibrationMode,
		long nbIterations,
		const vector< string >& flags,
        long mktDataManagerId,
        const vector< string >& keys,
        ARM_result&	result, 
        long        objId );


extern long ARMLOCAL_MaturityCap_Get(
        long maturityCapId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_MaturityCap_Set(
        long maturityCapId,
        long dataId,
        const string& setType,
		const vector<string>& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

long ARMLOCAL_MaturityCap_GetPricingData(
		long maturityCapId,
		const string& key,
		ARM_GramFctorArg& argResult,
		ARM_result& result );

/****************************************************************************
						CALLABLE SB CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_CallableSBCalculator_Create(
        double startDate,
        double endDate,
		long cpnFreq,
		long cpnDaycount,
		long cpnIndexTiming,
		const string& cpnIdxTerm,
		long cpnIndexDaycount,
		const string& cpnResetCal,
		const string& cpnPayCal,
		long cpnIntRule,
		long cpnIndexResetLag,
		long payRec,
		long CF,
		double notional,
		long notionalId,
		double fundnotional,
		long fundnotionalId,
		double dConst,
		long constId,
		double lPrevCpn,
		long lPrevCpnId,
		double lNewOpt,
		long lNewOptId,
		double strikeOpt,
		long strikeOptId,
		double minCpn,
		long minCpnId,
		double maxCpn,
		long maxCpnId,
		long fundFreq,
		long fundDaycount,
		double fundCoeff,
		long fundCoeffId,
		double fundMargin,
		long fundMarginId,
		long NotifDays,
		long NonCall,
		const string& exerciseCal,
		bool CallSBOrStacky,
		long feesId,
		long NbPathBounding,
		long NbPathPricing,
		long NbMaxBucket,
		long FixSteps,
		const string& USMethod,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
		const string& ModelType,
		const string& TriggerMode,
		int calibSwopFreq,
		const string& regressors,
		const vector< double >& hkDatas,
        ARM_result&	result, 
        long        objId );


long ARMLOCAL_CALLABLESB_Get(
        long callableSBId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CALLABLESB_Set(
        const long& csbId,
        const long& dataId,
        const string& setType,
		const vector< string >& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId );


/****************************************************************************
						GLOBAL CAP CALCULATOR
*****************************************************************************/

long ARMLOCAL_GlobalCapCalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		double notional,
		long notionalId,
		string fundIndexType,
		int fundDayCount,
		int fundFreq,
		int fundResetGap,
		int fundPayGap,
		int fundResetTiming,
		int fundPayTiming,
		string fundResetCal,
		string fundPayCal,
		int fundAdjRule,
		int fundIntRule,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		ARM_Vector* globalCapParams,
		long pastFixingsId,
		double capFixed,
		long capFixedId,
		double capStrike,
		long capStrikeId,
		double capSpread,
		long capSpreadId,
		int	nbSteps,
		vector<string> randGenerator,
		int samplerType,
		ARM_Vector* calibParams,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId);

long ARMLOCAL_GlobalCapCalculator_Create(
		long globalCapId,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		int	nbSteps,
		vector<string> randGenerator,
		int samplerType,
		ARM_Vector* calibParams,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId);

long ARMLOCAL_GlobalCapCalculator_Create_WithoutMktData(
		double asOfDate,
		long globalCapId,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		ARM_result&	result, 
        long objId);

extern long ARMLOCAL_GlobalCap_Set(
		long globalCapId,
		long dataId,
		const string& setType,
		const vector<string>& keys,
		bool isUpdated,
		ARM_result& result, 
		long objId );

/****************************************************************************
						SNOW RANGE CALCULATOR
*****************************************************************************/

long ARMLOCAL_SnowRangeCalculator_Create(
		ARM_Currency ccy,
		double endDate,
		double startDate,
		int payReceive,
		double notional,
		long notionalId,
		string fundingIndexTerm,
		int fundingDayCount,
		string couponIndexTerm,
		int couponDayCount,
		int resetFreq,
		int payFreq,
		int resetTiming,
		int payTiming,
		int resetGap,
		string resetCal,
		string payCal,
		int adjRule,
		int intRule,
		double spread,
		long spreadId,
		double strike,
		long strikeId,
		double ratchet,
		long ratchetId,
		double cashFlow,
		long cashFlowId,
		double fixedRate,
		long fixedRateId,
		double leverage,
		long leverageId,
		ARM_Vector* snowRangeParams,
		ARM_Vector* calibParams,
		string modelName,
		ARM_Vector* modelParams,
		int	nbSteps,
		string generatorType,
		string inversionMethod,
		bool antithetic,
		int samplerType,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId);

//===========================================================//
// VolBondCalculator										 //
//===========================================================//

extern long ARMLOCAL_VolBondCalculator_Create(
		const double nominal,
		const double startDate,
		const double endDate,
		const string payFreq,
		const string resetFreq,
		const string dayCount,
		const string tenor,
		const string intrule,
		const string stubrule,
		const double resetgap,
		const string paycalendar,
		const string resetcalendar,
		const string type_odesolver,
		const vector<double>& RK4parameters,
		const vector<double>& RK5parameters,
		const double mcnbsteps,
		const double mcbucket,
		const double mcnbstepsbyyear,
		const vector<string>& randomgenerator,
		const long marketdatamanager_id,
		const vector<string>& marketdatamanagerkeys,
		const string payofftype,
		const vector<string>& productstoprice,
		ARM_result& result,
        long objId = -1);


#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
