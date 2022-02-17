/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file HybridIRFXCalculator.h
 *
 *  \brief file for the hybrid IR FX Calculator
 *	\author P. Lam
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INGPCALCULATORS_HYBRIDIRFXCALCULATOR_H
#define _INGPCALCULATORS_HYBRIDIRFXCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gencalculator.h"
#include "gpmodels/Mixture_FX.h"
//#include "crv/fixingsched.h"

#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_2IRFXModel;
class ARM_1IRFXModel;
class ARM_NP1IRNFXModel;

/// H&W vol range [10bp,500bp]
const double HWVOL_LOWER_BOUND      = 0.00001;
const double HWVOL_UPPER_BOUND      = 0.05;

/// Q vol range [2.5%,100%]
const double QVOL_LOWER_BOUND       = 0.025;
const double QVOL_UPPER_BOUND       = 1.0;

/// 0.001bp of vega to be selected in portfolio for volatility bootstrapping 
const double IR_VEGA_MIN_TO_SELECT=0.0000001;
const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

/// 10-3 bp of vega to be selected in portfolio for FX volatility bootstrapping 
const double FX_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double FX_DEFAULT_WEIGHT      = 1.0;
const double FX_DEFAULT_PRICE       = 1.0e+30;

/// Reference schedules for TARN date structure
const unsigned int CPNFX_SCHED = 0;
const unsigned int NB_HybridIRFX_SCHED = 1;

///-----------------------------------------------------------------------------
/// \class ARM_HybridIRFXCalculator
/// \brief
///  Class that implements a TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_HybridIRFXCalculator : public ARM_GenCalculator 
{
	public:
		
		static const string ColNamesTable [];
		static const int ProductToPriceColumns[];

		enum productToPriceAlias
		{
			FundingPrice,
			CouponPrice,
			DFValue,
			PaidCouponPrice,
			RealCouponPrice,
			IsAliveValue,
			CapPrice,
			FloorPrice,
			RedemptionPrice,
			RealRedemptionPrice,
			RealFundingPrice,
			DurationValue,
			ProbaValue,
			TARNPrice,
			NbProductsToPrice
		};

		enum ColAlias
		{
			ResetDate = 0,
			NextResetDate,
			StartDate,
			EndDate,
			FundingStartDate,
			FundingEndDate,
			IT,
			RedemptionResetDate,
			MinCpn,
			MaxCpn,
			Notional,
			Target,
			DF,
			Offset,
			Width,
			FixedFunding,
			Funding,
			InitFX,
			DomCpn,
			FgnCpn,
			FxSpot1,
			FxSpot2,
			FxSpot3,
			Coupon,
			PaidCoupon,
			SumCoupons,
			Fees,
			IsAlive,
			RealCoupon,
			DiscountFunding,
			RealFunding,
			Cap,
			Floor,
			RedemptionStrike,
			Redemption,
			RealRedemption,	
			TARN,	
			Duration,
			Proba,
			InitFX2,
			DomCpn2,
			FgnCpn2,
			RedemptionStrike2,
			InitFX3,
			DomCpn3,
			FgnCpn3,
			RedemptionStrike3
		};

		enum ModelType
		{
			Model1IRFX,
			Model2IRFX,
			ModelNP1IRNFX
		};

		ARM_HybridIRFXCalculator(const ARM_Date& asOfDate,
								 const ARM_Date& startDate,
								 const ARM_Date& endDate,
								 const ARM_Currency& DomCcy,
								 const vector<ARM_Currency>& ForCcy,
								 const ARM_Currency& CpnCcy,
								 const ARM_Currency* FundCcy,
								 int payRec,
								 int cpnDayCount,
								 int cpnFreq,
								 int cpnResetGap,
								 const string& cpnResetCal,
								 const string& cpnPayCal,
								 int stubRule,
								 int cpnTiming,
								 int cpnIntRule,
								 const ARM_Curve* cpnNominal,
								 const ARM_Curve* domesticCpn,
								 const ARM_Curve* foreignCpn,
								 const ARM_Curve* MinCpn,
								 const ARM_Curve* MaxCpn,
								 const ARM_Curve* InitialFX,
								 int fundFreq,
								 int fundDayCount,
								 const ARM_Curve* fundNominal,
								 const ARM_Curve* fundSpread,
								 int redemptionType,
								 int redemptionGap,
								 double redemptionStrike,
								 const ARM_Curve* target,
								 const ARM_Curve* fees,
								 const string& FXChoice,
								 int intermediatePrices,
								 const ARM_StringVector& columnsToPrice,
								 const ARM_FXTARNPayoffType& payOffName,
								 ARM_FixingSched* pastFixings = NULL,
								 char* refDate = GETDEFAULTVALUESTR,
								 ARM_Date& effDate = ARM_Date(ARM_DEFAULT_DATE));

		ARM_HybridIRFXCalculator(const ARM_Date& asOfDate,
								 const ARM_Date& startDate,
								 const ARM_Date& endDate,
								 const ARM_Currency& DomCcy,
								 const vector<ARM_Currency>& ForCcy,
								 const ARM_Currency& FundCcy,
								 int payRec,
								 int cpnDayCount,
								 int cpnFreq,
								 int cpnResetGap,
								 const string& cpnResetCal,
								 const string& cpnPayCal,
								 int stubRule,
								 int cpnTiming,
								 int cpnIntRule,
								 const ARM_Curve& cpnNominal,
								 const vector<ARM_Curve>& domesticCpn,
								 const vector<ARM_Curve>& foreignCpn,
								 const ARM_Curve& MinCpn,
								 const ARM_Curve& MaxCpn,
								 const vector<ARM_Curve>& InitialFX,
								 int fundFreq,
								 int fundDayCount,
								 const ARM_Curve& fundNominal,
								 const ARM_Curve& fundSpread,
								 int redemptionType,
								 int redemptionGap,
								 std::vector<double>& redemptionStrike,
								 const ARM_Curve& target,
								 const ARM_Curve& fees,
								 const string& FXChoice,
								 int intermediatePrices,
								 const ARM_StringVector& columnsToPrice,
								 const ARM_FXTARNPayoffType& payOffName,
								 ARM_FixingSched* pastFixings = NULL,
								 char* refDate = GETDEFAULTVALUESTR,
								 ARM_Date& effDate = ARM_Date(ARM_DEFAULT_DATE));

		void Init(
			int nbSimul,
			int bucketSize,
			int randGendType1,
			int randGenAlpgo1,
			int randGenType2,
			int randGenAlgo2,
			int firstNbDims,
			int firstNbTimes,
			int factorNb,
			int timeStepNb,
			int spaceStepNb,
			double stdDevNb,
			int skipPDE,
			int rescalling,
			int modelType,
			int smileFlag,
			int mixCalib,
			int oneFactorFlag,
			int correlType,
			const ARM_MarketData_ManagerRep& mktDataManager);

		void Init(
			int nbSimul,
			int bucketSize,
			int randGendType1,
			int randGenAlpgo1,
			int randGenType2,
			int randGenAlgo2,
			int firstNbDims,
			int firstNbTimes,
			int factorNb,
			int timeStepNb,
			int spaceStepNb,
			double stdDevNb,
			int skipPDE,
			int rescalling,
			int modelType,
			int smileFlag,
			int mixCalib,
			int oneFactorFlag,
			int correlType,
			vector<ARM_ZeroCurve*> zeroCurves,
			vector<ARM_ZeroCurve*> basisCurves,
			vector<ARM_Forex*> forex,
			vector<ARM_VolCurve*> ATMVol, //for swopt BSGen
			vector<ARM_VolCurve*> fxVol, //for BS fx models
			vector<ARM_ParamsMixture_Fx*> mixtureParams, //for mixture fx models
			vector<ARM_CurveModelParam*> mrsParams,
			vector<ARM_CurveModelParam*> QParams,
			ARM_GP_Matrix* correlMatrix);

		///copy constructor, assignment constructor, destructor
		ARM_HybridIRFXCalculator( const ARM_HybridIRFXCalculator& rhs );
		ASSIGN_OPERATOR(ARM_HybridIRFXCalculator)
		~ARM_HybridIRFXCalculator();

		virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const {};

		// Basis
		void ComputeDomesticBasis();
		
		/// Core of the calculator
		ARM_2IRFXModel* CreateModel2IRFX();
		ARM_1IRFXModel* CreateModel1IRFX();
		ARM_NP1IRNFXModel* CreateModelNP1IRNFX();
		virtual void CreateAndSetModel();
		vector<ARM_DensityFunctor*> CreateDensityFunctor();
		ARM_StdPortfolioPtr CreateDiagonalSwaption(ARM_Currency& ccy);
		void ComputeSwaptionPrice(ARM_Currency& ccy, ARM_CalibMethod* IRVolCalibMethod);
		ARM_CalibMethod* CreateIRCalibMethod(const ARM_StdPortfolioPtr& diagonalSwaptionPF,int modelIdx);
		ARM_StdPortfolioPtr CreateFxOption(ARM_Currency& foreignCcy);
		void ComputeFxOptionPrices(ARM_CalibMethod* FXVolCalibMethod, int Ccyi);
		ARM_CalibMethod* CreateNumFxCalibMethod();
		virtual ARM_CalibMethod* CreateCalibration(ModelType modelType);
		ARM_CalibMethod* CreateFxCalibMethod(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);
		virtual void CreateAndSetCalibration();
		virtual void UpdateModel();
		virtual void UpdateCalibration(bool isUpdateStrike=true);
		virtual void Calibrate();
		virtual double Price();
		virtual ARM_Vector* ComputeAll() {return NULL;};
		virtual void CheckData(){};
		virtual void CheckMktData(){};
		virtual void ComputePricingData() const;

		/// TARN FX deal description creation functions
		virtual ARM_RowInfo ColumnNames() const;
		virtual ARM_DateStripCombiner DatesStructure() const;
		virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;

		virtual double PastLiborValue(ARM_DateStripPtr dateStrip, double& IT, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec){return 0.0;};
		virtual string FXValue(ARM_DateStripPtr dateStrip, size_t eventIdx, int FXNb=0);
		virtual void DoPastReset(ARM_DateStripPtr dateStrip, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec){};

		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_HybridIRFXCalculator(*this); }

		virtual string ExportShortName() const { return "LIRFX";}

		/// Dates Strip
		inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
		inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
		inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

		/// Vector
		inline virtual std::vector<double> GetvCpnNominal() const  { return std::vector<double>(0);};
		inline virtual std::vector<double> GetvFundNominal() const { return std::vector<double>(0);};
		inline virtual std::vector<double> GetvFundSpread() const { return std::vector<double>(0);};

		/// the discount curve
		inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  { return NULL;};
		inline virtual ARM_Forex* GetForex() const  { return NULL;};
		inline bool IsNormal() const { return itsIsNormal; };

	protected :

		ARM_Date				itsStartDate;
		ARM_Date				itsEndDate;
		ARM_Date				itsRedemptionResetDate;
		char*					itsReferenceDate;
		ARM_Date				itsEffectiveDate;
		int						itsNbFX;
		ARM_Currency			itsDomesticCcy;
		vector<ARM_Currency>	itsForeignCcy;
		ARM_Currency			itsFundingCcy;
		ARM_Currency			itsCouponCcy;
		int						itsPayRec;
		int						itsCpnDayCount;
		int						itsCpnFreq;
		int						itsCpnResetGap;
		string					itsCpnResetCal;
		string					itsCpnPayCal;
		int						itsStubRule;
		int						itsCpnTiming;
		int						itsCpnIntRule;
		ARM_Curve				itsCpnNominalCv;
		ARM_Curve				itsFundNominalCv;
		ARM_Curve				itsFundSpreadCv;
		vector<ARM_Curve>		itsDomesticCpnCv;
		vector<ARM_Curve>		itsForeignCpnCv;
		ARM_Curve				itsMinCpnCv;
		ARM_Curve				itsMaxCpnCv;
		vector<ARM_Curve>		itsInitialFXCv;
		int						itsFundFreq;
		int						itsFundDayCount;
		int						itsRedemptionType;
		int						itsRedemptionGap;
		std::vector<double>			itsRedemptionStrike;
		ARM_Curve				itsFees;
		int						itsNbPastFixings;
		ARM_Curve				itsTarget;
		int						itsIntermediatePrices;
		
		ARM_IntVector			itsvFundIndex;
		string					itsFXChoice;
		ARM_FXTARNPayoffType	itsPayOffName;
		string					itsPayModelName;
		bool					itsIsNormal;

		int itsNbSimul;
		int itsBucketSize;
		int itsRandGenType1;
		int itsRandGenAlgo1;
		int itsRandGenType2;
		int itsRandGenAlgo2;
		int itsFirstNbDims;
		int itsFirstNbTimes;
		int itsFactorsNb;
		int itsTimeStepNb;
		int itsSpaceStepNb;
		int itsStdDevNb;
		int itsSkipPDE;
		int itsRescalling;
		int itsModelType;
		int itsSmileFlag;
		int itsMixCalib;
		int itsOneFactorFlag;
		int itsCorrelType;
		
		int	itsFirstEventIdx;
		mutable ARM_DateStripPtr itsCpnDateStrip;
		mutable ARM_DateStripPtr itsDomDateStrip;
		mutable ARM_DateStripPtr itsFundDateStrip;
		ARM_StringVector itsColumnsToPrice;
		ARM_StdPortfolioPtr itsDomSwaptionPF;
		vector<ARM_StdPortfolioPtr> itsForSwaptionPF;
		vector<ARM_StdPortfolioPtr> itsFxOptionPF;
		bool itsHasBeenComputed;

		ARM_CalibMethodPtr itsSmiledFxCalibMethod;

		vector<string> itsRandGenParams;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

