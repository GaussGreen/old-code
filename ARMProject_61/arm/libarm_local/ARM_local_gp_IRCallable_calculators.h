/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_IRCALLABLE_CALCULATORS_H
#define ARMLOCAL_GP_IRCALLABLE_CALCULATORS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;


/****************************************************************************
						CRF CALCULATOR
*****************************************************************************/

bool CRFCalculatorCheckAutoCalFlags(const vector< string >& flags,size_t nbFlagsToCheck);

extern long ARMLOCAL_CRFCalculator_Create(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        long payRec,
        double fixEndDate,
        long fixDayCount,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long  cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
		long stubRule,
        long cpnResetGap,
        double leverage,
        long leverageId,
        double cpnMin,
        long cpnMinId,
        double cpnMax,
        long cpnMaxId,
        double fundSpread,
        long fundSpreadId,
        long fundFreq,
        long fundDayCount,
        double nominal,
        long nominalId,
        long exerGap,
        long nbNonCall,
        double exerFee,
        long exerFeeId,
        const vector< string >& flags,
        long mktDataManagerId,
        const vector< string >& keys,
        double fundnominal,
        long fundnominalId,
        ARM_result&	result, 
        long        objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_Crude_CRFCalculator_Create(
        const double & startDate,
        const double & endDate,
        const double & strike,
        const long & strikeId,
        const long & payRec,
        const double & fixEndDate,
        const long & fixDayCount,
        const long & cpnDayCount,
        const long & cpnFreq,
        const long & cpnTiming,
        const string & cpnIndexTerm,
        const long &  cpnIndexDayCount,
        const string & cpnResetCal,
        const string & cpnPayCal,
		const long & stubRule,
        const long & cpnResetGap,
        const double & leverage,
        const long & leverageId,
        const double & cpnMin,
        const long & cpnMinId,
        const double & cpnMax,
        const long & cpnMaxId,
        const double & fundSpread,
        const long & fundSpreadId,
        const long & fundFreq,
        const long & fundDayCount,
        const double & nominal,
        const long & nominalId,
        const long & exerGap,
        const long & nbNonCall,
        const double & exerFee,
        const long & exerFeeId,
        const double & fundnominal,
        const long & fundnominalId,
        ARM_result&	result, 
        long        objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_CRF_Get(
        long crfId,
        const string& getType,
        ARM_result&	result, 
        long        objId );
		
extern long ARMLOCAL_CRF_GetMRS(
        long crfId,
        ARM_result&	result,
		long objId = -1);

extern long ARMLOCAL_CRF_Set(
        long crfId,
        long dataId,
        const string& setType,
		const vector<string>& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CRF_SetAutoCalFlags(
	long crfId,
    vector< string >& flags,
    ARM_result&	result);

extern long ARMLOCAL_CRF_SetOneCalFlag(
	const long& crfId,
    const string& flagStr,
	const string& typeStr,
    ARM_result&	result ,
	long objId);

extern long ARMLOCAL_CRF_SetAutoCalFlagsAndClone(
	const long& crfId,
    const vector< string >& flags,
    ARM_result&	result,
	long        objId);


bool CRF_SetProductToPrice(ARM_Object* crf, const string& prodToPrice);

extern long ARMLOCAL_CRF_SetProductToPrice(
	long crfId,
	string prodToPrice,
	ARM_result&	result );

/****************************************************************************
						CRA CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_CRA_GetModel(
		long craId,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CRA_GetCalibData(
		long craId,
		string calibOrPortfolio,
		string dataType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CRA_GetUnderlying(
		long craId,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CRALocalCalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		int fundFreq,
		int fundDayCount,
		int cpnPayFreq,
		string cpnResetCal,
		string cpnPayCal,
		int boostedIndexType,
		string boostedVarTerm,
		int boostedResetGap,
		int boostedResetTiming,
		int boostedDayCount,
		int boostedAdjRule,
		int boostedRule,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnResetGap,
		int refIndexType,
		string  refTerm,
		int refDayCount,
		double refCoeff,
		double notional,
		long notionalId,
		long callFeesId,
		double fundSpread,
		long fundSpreadId,
		double cpnSpread,
		long cpnSpreadId,
		double boostedFix,
		long boostedFixId,
		double bDown,
		long bDownId,
		double bUp,
		long bUpId,
		bool localModel,
		int localModelType,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		int localResetFreq,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool isStdCalib,
		ARM_result&	result, 
        long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_LocalCRA_Set(long craId,
							 	  long dataId,
								  const string& setType,
								  const vector<string>& keys,
								  bool isUpdated,
								  ARM_result& result, 
								  long objId );

extern long ARMLOCAL_LocalCRA_Get(long craId,
						 		  const string& getType,
								  ARM_result& result, 
								  long objId );


extern long ARMLOCAL_CRACalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		int fundFreq,
		int fundDayCount,
		int cpnPayFreq,
		string cpnResetCal,
		string cpnPayCal,
		int boostedIndexType,
		string boostedVarTerm,
		int boostedResetGap,
		int boostedResetTiming,
		int boostedDayCount,
		int boostedAdjRule,
		int boostedRule,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnResetGap,
		int refIndexType,
		string  refTerm,
		int refDayCount,
		double refCoeff,
		double notional,
		long notionalId,
		long callFeesId,
		double fundSpread,
		long fundSpreadId,
		double cpnSpread,
		long cpnSpreadId,
		double boostedFix,
		long boostedFixId,
		double bDown,
		long bDownId,
		double bUp,
		long bUpId,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool isStdCalib,
		ARM_result&	result, 
        long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_CRACalculator_Create(
		int OptionPfId,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool localModel,
		int localModelType,
		int localResetFreq,
		bool isStdCalib,
		ARM_result&	result, 
        long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_CRA_Set(long craId,
							 long dataId,
							 const string& setType,
							 const vector<string>& keys,
							 bool isUpdated,
							 ARM_result& result, 
							 long objId );

extern long ARMLOCAL_CRA_Get(long craId,
							 const string& getType,
							 ARM_result& result, 
							 long objId );

extern long ARMLOCAL_CRA_SetRecalibFlags(long craId,
									     int mrsFlag,
									     int betaFlag,
									     ARM_result& result,
										 long objId = -1);

/****************************************************************************
					BERMUDA SWAPTION CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_BERMUDASWAPTION_Get(
		long bsId,
        const string& getType,
        ARM_result&	result, 
        long        objId );
		
extern long ARMLOCAL_BERMUDASWAPTION_Set(
        long bsId,
        long dataId,
        const string& setType,
		const vector< string >& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_BERMUDASWAPTIONCalculator_Create(
		ARM_Currency ccy,
	    double startDate,
		double endDate,
		double notional,
		long notionalId,
		double strike,
		long strikeId,
		double fees,
		long feesId,
		double spread,
		long spreadId,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		double firstCallDate,
		double lastCallDate,
		int fixFreq,
		int fixBasis,
		string fixPayCal,
		int fixAdjRule,
		int fixRule,
		int fixPayGap,
		bool isZc,
		int varFreq,
		int varBasis,
		string varResetCal,
		string varPayCal,
		string varIndexTerm,
		int varAdjRule,
		int varRule,
		int varResetGap,
		int varPayGap,
		int stubRule,
		int genSecType,
		vector<int>* controlVariates,
		vector<double>* controlPrices,
		vector<string> mdmKeys,
		long mktDataManagerId,
		int modelType,
        vector<string> modelParams,
		vector<string> calibFlags,
		int numMethodType,
		int amcIter,
		int mcIter,
		int maxBucketSize,
		string genType1,
		string genType2,
		string pathOrder,
		string pathScheme,
		int firstNbDims,
		int treeSteps,
		vector<int> portfolioMode,
		bool boundaryFlag,
		bool approxMarginFlag,
		bool freezeBetasFlag,
		bool calculateProbaFlag,
		ARM_result&	result, 
        long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_ARM_BermudaRootMrs(
        long		bsId,
        double		targetPrice,
		double		fTolerance,
		int			maxIter,
        ARM_result&	result);


/****************************************************************************
						CAPTION CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_CaptionCalculator_Create(
        double startDate,
        double endDate,
		const string& cpnIdxTerm,
		long payRec,
		long CF,
		double coupon,
        long couponId,
		const string& FundIdxTerm,
		long NotifDays,
		long NonCall,
		double exercise,
		long exerStyleId,
		double notional,
        long notionalId,
		long cpnDayCount,
		long cpnResetTiming,
		const string& cpnResetCal,
        const string& cpnPayCal,
		double cpnSpread,
        long cpnSpreadId,
		long fundDayCount,
		const string& fundResetCal,
        const string& fundPayCal,
		double fundSpread,
        long fundSpreadId,
		long factorNb,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys,
        ARM_result&	result, 
        long        objId);


extern long ARMLOCAL_Caption_Get(
        long captionId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_Caption_Set(
        long captionId,
        long dataId,
        const string& setType,
		const vector<string>& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

long ARMLOCAL_CaptionCalculator_GetPricingData(
		long captionCalculatorId,
		const string& key,
		ARM_GramFctorArg& argResult,
		ARM_result& result );


/****************************************************************************
						CSO CALCULATOR
*****************************************************************************/

/// CSO calculator: standard version
extern long ARMLOCAL_CSOCalculator_Create(
						double startDate,
						double endDate,
						long CMS1Type,
						long CMS2Type,
						long cpnFreq,
						long cpnDaycount,
						long fundFreq,
						long fundDaycount,
						long exerFreq,
						long NotifDays,
						double notional,
						long notionalId,
						double minCpn,
						long minCpnId,
						double maxCpn,
						long maxCpnId,
						double leverage,
						long leverageId,
						double fundMargin,
						long fundMarginId,
						long fixCpnId,
						long feesId,
						const vector< string >& CalibMod,
						const vector< string >& ProductFlags,
						const vector< double >& ModelDatas,
						long mktDataManagerId,
						const vector< string >& keys,
						ARM_result&	result, 
						long        objId );

/// CSO calculator: extended version
long ARMLOCAL_ExtendedCSOCalculator_Create (
        const double&		startDate,
        const double&		endDate,
        const long &		CMS1Type,
        const long &		CMS2Type,
        const long &		cpnFreq,
        const long &		cpnDaycount,
        const long &		cpnResetTiming,
        const long &		fundFreq,
        const long &		fundDaycount,
        const long &		fundResetTiming,
        const long &		exerFreq,
        const long &		NotifDays,
        const long &        payRec,
        const double &		cpnnotional,
        const long &		cpnnotionalId,
        const double &		minCpn,
        const long &		minCpnId,
        const double &		maxCpn,
        const long &		maxCpnId,
        const double &		leverageLong,
        const long &		leverageLongId,
        const double &		leverageShort,
        const long &		leverageShortId,
        const double &		strike,
        const long &		strikeId,
        const double &		fundMargin,
        const long &		fundMarginId,
		const double &		fundLeverage,
        const long &		fundLeverageId,
        const long &		fixCpnId,
        const long &		feesId,
		const bool &		switchFlag,
		const long &		fundType,
        const vector< string >&	CalibMod,
        const vector< string >&	ProductFlags,
        const vector< double >&	ModelDatas,
        const long &		mktDataManagerId,
        const vector< string >&	keys,
        ARM_result&			result, 
        long				objId );


extern long ARMLOCAL_BasicCSOCalculator_Create(const double&	asOfDate,
                                               const double&	startDate,
                                               const double&	inFixEndDate,
                                               const double&	endDate,
                                               const long&		CMS1Type,
                                               const long&		CMS2Type,
                                               const long&		cpnFreq,
                                               const long&		cpnDaycount,
                                               const long&		cpnResetTiming,
                                               const long&		fundFreq,
                                               const long&		fundDaycount,
                                               const long&		fundResetTiming,
                                               const long&		exerFreq,
                                               const long&		NotifDays,
                                               const long&		payRec,
                                               const double&	cpnnotional,
                                               const long&		cpnnotionalId,
                                               const double&	minCpn,
                                               const long&		minCpnId,
                                               const double&	maxCpn,
                                               const long&		maxCpnId,
                                               const double&	leverageLong,
                                               const long&		leverageLongId,
                                               const double&	leverageShort,
                                               const long&		leverageShortId,
                                               const double&	strike,
                                               const long&		strikeId,
                                               const double&	fundNotional,
                                               const long&		fundNotionalId,
                                               const double&	fundMargin,
                                               const long&		fundMarginId,
		                                       const double&	fundLeverage,
                                               const long&		fundLeverageId,
                                               const long&		fixCpnId,
                                               const long&		feesId,
                                               const CCString&	CpnCcy,
                                               const CCString&	FundCcy,
											   const long&		nbNoCall,
                                               ARM_result&		result, 
                                               long				objId );


extern long ARMLOCAL_BasisCSOCalculator_Create(const double&	startDate,
		                                       const double&	fixEndDate,
                                               const double&	endDate,
		                                       const string&     cpnCcy,
		                                       const string&     fundCcy,
                                               const long &		CMS1Type,
                                               const long &		CMS2Type,
                                               const long &		cpnFreq,
                                               const long &		cpnDaycount,
                                               const long &		cpnResetTiming,
		                                       const long &     stubRule,
                                               const string&    cpnResetCal,
                                               const string&    cpnPayCal,
		                                       const long &		cpnnotionalId,
		                                       const long &		strikeId,
		                                       const long &		leverageLongId,
                                               const long &		leverageShortId,
                                               const long &		minCpnId,
                                               const long &		maxCpnId,
                                               const long &		fundFreq,
                                               const long &		fundDaycount,
                                               const long &		fundResetTiming,
		                                       const long &		fundingnotionalId,
		                                       const long &		fundMarginId,
                                               const long &		exerFreq,
                                               const long &		NotifDays,
                                               const long &		payRec,
		                                       const long &     nbNCall,
                                               const long &		feesId,
                                               const vector< string >&	CalibMod,
                                               const vector< string >&	ProductFlags,
                                               const vector< double >&	ModelDatas,
                                               const long &		mktDataManagerId,
                                               ARM_result&			result, 
                                               long				objId);

extern long ARMLOCAL_CSO_Get(
        const long& CSOId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CSO_Set(
        const long& calculatorId,
        const long& dataId,
        const string& setType,
		const vector< string >& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId );


/****************************************************************************
						CRA Spread CALCULATOR
*****************************************************************************/

long ARMLOCAL_CRASpreadCalculator_Create(
		const ARM_Currency& ccy,
		const ARM_Currency& fundCcy,
		const double& startDate,
		const double& endDate,
		const int& payReceive,
		const int& callFreq,
		const int& callNotice,
		const string& callCal,
		const int& fundFreq,
		const int& fundDayCount,
		const int& cpnDayCount,
		const int& cpnPayFreq,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const int& cpnResetFreq,
		const int& cpnResetTiming,
		const int& cpnResetGap,
		const int& refIndex1,
		const int& refIndex2,
		const int& payIndex,
		const int& payIndexResetTiming,
		const double& notional,
		const long& notionalId,
		const double& fundNotionaldble,
		const long& callFeesId,
		const double& fundSpread,
		const long& fundSpreadId,
		const double& boostedFix,
		const long& boostedFixId,
		const double& payIndexMult,
		const long& payIndexMultId,
		const double& bDown,
		const long& bDownId,
		const double& bUp,
		const long& bUpId,
		const double& coeff1,
		const double& coeff2,
		const int& refCond1Index,
		const double& bCond1Down,
		const long& bCond1DownId,
		const double& bCond1Up,
		const long& bCond1UpId,
		const long& boostedFix2Id,
		const long& bDown2Id,
		const long& bUp2Id,
		const long& boostedFix3Id,
		const long& bDown3Id,
		const long& bUp3Id,
		const vector<string>& tenorVec,
		const vector<string>& calibMod,
		const vector<string>& mdmKeys,
		const long& mktDataManagerId,
		const vector<string>& pricingFlags,
		const vector<string>& localCalibFlags,
		const vector<double>& optimResetData,
		const vector<double>& exerProbas,
		ARM_result&	result,
        long objId);

long ARMLOCAL_CRASpreadCalculator_Create(int OptionPfId,
										int refResetFreq,
										vector<string> mdmKeys,
										vector<string> calibMod,
										double payIndexMultValue,
										long payIndexMultId,
										long mktDataManagerId,
										vector<string> pricingFlags,
										ARM_result&	result, 
										long objId = ARM_NULL_OBJECT_ID);

long ARMLOCAL_CRASpreadCalculator_Create(double asOfDate,
										int securityId,
										int refResetFreq,
										double payIndexMultValue,
										long payIndexMultId,
										ARM_result&	result, 
										long objId = ARM_NULL_OBJECT_ID);

long ARMLOCAL_CRASpreadCalculator_Create(const ARM_Currency& ccy,
										 const double& startDate,
										 const double& endDate,
										 const int& payReceive,
										 const int& callFreq,
										 const int& callNotice,
										 const string& callCal,
										 const int& fundFreq,
										 const int& fundDayCount,
										 const int& cpnDayCount,
										 const int& cpnPayFreq,
										 const string& cpnResetCal,
										 const string& cpnPayCal,
										 const int& cpnResetFreq,
										 const int& cpnResetTiming,
										 const int& cpnResetGap,
										 const int& refIndex1,
										 const int& refIndex2,
										 const int& payIndex,
										 const int& payIndexResetTiming,
										 const double& notional,
										 const long& notionalId,
										 const long& callFeesId,
										 const double& fundSpread,
										 const long& fundSpreadId,
										 const double& boostedFix,
										 const long& boostedFixId,
										 const double& payIndexMult,
										 const long& payIndexMultId,
										 const double& bDown,
										 const long& bDownId,
										 const double& bUp,
										 const long& bUpId,
										 const double& coeff1,
										 const double& coeff2,
										 ARM_result&	result,
										 long objId);

extern long ARMLOCAL_CRASpread_Get(
        const long& CRASpreadId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_CRASpread_Set(
        const long& CRASpreadId,
        const long& dataId,
        const string& setType,
		const vector<string>& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId );

#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
