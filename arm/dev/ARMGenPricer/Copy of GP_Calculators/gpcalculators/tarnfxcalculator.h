/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.h
 *
 *  \brief file for the TARN FX Calculator
 *	\author  A. Lekrafi
 *	\version 1.0
 *	\date August 2006
 */


#ifndef _INGPCALCULATORS_TARNFXCALCULATOR_H
#define _INGPCALCULATORS_TARNFXCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "hybridirfxcalculator.h"
#include "gpmodels/Mixture_FX.h"

#include "typedef.h"
#include "refvalue.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_2IRFXModel;
class ARM_1IRFXModel;
class ARM_NP1IRNFXModel;


///-----------------------------------------------------------------------------
/// \class ARM_TARNFXCalculator
/// \brief
///  Class that implements a TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_TARNFXCalculator : public ARM_HybridIRFXCalculator 
{
	public:
		
		static const string TARNFXColNamesTable [];
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

	public:

		enum TARNFXColAlias
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
			Funding,
			InitFX,
			DomCpn,
			FgnCpn,
			FxSpot,
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

		/// constructor for deal from term sheet
	/*	ARM_TARNFXCalculator(const ARM_Date& asOfDate,
							 const ARM_Date& startDate,
							 const ARM_Date& endDate,
							 const ARM_Currency& DomCcy,
							 const ARM_Currency& ForCcy,
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
							 double TargetRedemption,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 double redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_IntVector& productsToPrice);*/

		ARM_TARNFXCalculator(const ARM_Date& asOfDate,
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
							 double TargetRedemption,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 ARM_GP_Vector& redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_IntVector& productsToPrice,
							 const ARM_FXTARNPayoffType& payOffName = ARM_TARNFXPayoffType::TARNFX);

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
		ARM_TARNFXCalculator( const ARM_TARNFXCalculator& rhs );
		ASSIGN_OPERATOR(ARM_TARNFXCalculator)
		~ARM_TARNFXCalculator();

		// Create the constant manager
		ARM_CstManagerPtr CreateCstManager();

		/// Initialisation to 0 of all columns of the deal description that could be priced
		void InitPriceableColumns(vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

		// Basis
//		void ComputeDomesticBasis();
		
		/// Core of the calculator
//		ARM_2IRFXModel* CreateModel2IRFX();
//		ARM_1IRFXModel* CreateModel1IRFX();
		ARM_NP1IRNFXModel* CreateModelNP1IRNFX();
		virtual void CreateAndSetModel();
//		ARM_StdPortfolioPtr CreateDiagonalSwaption(ARM_Currency& ccy);
//		void ComputeSwaptionPrice(ARM_Currency& ccy, ARM_CalibMethod* IRVolCalibMethod);
//		ARM_CalibMethod* CreateIRCalibMethod(const ARM_StdPortfolioPtr& diagonalSwaptionPF,int modelIdx);
//		ARM_StdPortfolioPtr CreateFxOption(ARM_Currency& foreignCcy);
//		void ComputeFxOptionPrices(ARM_CalibMethod* FXVolCalibMethod, int Ccyi);
//		ARM_CalibMethod* CreateFxCalibMethod(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);
		vector<ARM_DensityFunctor*> CreateDensityFunctor();
		ARM_CalibMethod* CreateNumFxCalibMethod();
		virtual ARM_CalibMethod* CreateCalibration(ModelType modelType);
//		virtual void CreateAndSetCalibration();
//		virtual void UpdateModel();
//		virtual void UpdateCalibration(bool isUpdateStrike=true);
		virtual void Calibrate();
		virtual double Price();
		virtual ARM_Vector* ComputeAll();
		virtual void CheckData();
		virtual void CheckMktData();
		virtual void ComputePricingData() const;

		/// TARN FX deal description creation functions
		virtual ARM_RowInfo ColumnNames() const;
		virtual ARM_StringVector ProductsToPriceColumnNames();
//		virtual ARM_DateStripCombiner DatesStructure() const;
//		virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_TARNFXCalculator(*this); }

		virtual string ExportShortName() const { return "LTRFX";}

		/// Dates Strip
		inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
		inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
		inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

		/// Vector
		inline virtual ARM_GP_Vector GetvCpnNominal() const  { return ARM_GP_Vector(0);};
		inline virtual ARM_GP_Vector GetvFundNominal() const { return ARM_GP_Vector(0);};
		inline virtual ARM_GP_Vector GetvFundSpread() const { return ARM_GP_Vector(0);};

		/// the discount curve
		inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;};
		inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  { return NULL;};
		inline virtual ARM_Forex* GetForex() const  { return NULL;};

	protected :
		ARM_IntVector			itsvFundIndex;
		int						itsCpnIntRule;
		double					itsTargetRedemption;
		string					itsFXChoice;
		ARM_FXTARNPayoffType     itsPayOffName;
		int						itsIntermediatePrices;
		ARM_IntVector			itsProductsToPrice;

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

