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


CC_BEGIN_NAMESPACE( ARM )

class ARM_2IRFXModel;
class ARM_1IRFXModel;
class ARM_NP1IRNFXModel;


///-----------------------------------------------------------------------------
/// \class ARM_HybridIRFXCalculator
/// \brief
///  Class that implements a TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_HybridIRFXCalculator : public ARM_GenCalculator 
{
	public:
		
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
								 int redemptionType,
								 int redemptionGap,
								 double redemptionStrike,
								 const ARM_Curve* fees);

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
								 ARM_GP_Vector& redemptionStrike,
								 const ARM_Curve& fees);

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
		virtual void CreateAndSetModel();
		ARM_StdPortfolioPtr CreateDiagonalSwaption(ARM_Currency& ccy);
		void ComputeSwaptionPrice(ARM_Currency& ccy, ARM_CalibMethod* IRVolCalibMethod);
		ARM_CalibMethod* CreateIRCalibMethod(const ARM_StdPortfolioPtr& diagonalSwaptionPF,int modelIdx);
		ARM_StdPortfolioPtr CreateFxOption(ARM_Currency& foreignCcy);
		void ComputeFxOptionPrices(ARM_CalibMethod* FXVolCalibMethod, int Ccyi);
		virtual ARM_CalibMethod* CreateCalibration(ModelType modelType);
		ARM_CalibMethod* CreateFxCalibMethod(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);
		virtual void CreateAndSetCalibration();
		virtual void UpdateModel();
		virtual void UpdateCalibration(bool isUpdateStrike=true);
		virtual void Calibrate();
		virtual double Price() {return 0.0;};
		virtual ARM_Vector* ComputeAll() {return NULL;};
		virtual void CheckData(){};
		virtual void CheckMktData(){};
		virtual void ComputePricingData() const {};

		/// TARN FX deal description creation functions
		virtual ARM_RowInfo ColumnNames() const;
		virtual ARM_StringVector ProductsToPriceColumnNames();
		virtual ARM_DateStripCombiner DatesStructure() const;
		virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_HybridIRFXCalculator(*this); }

		virtual string ExportShortName() const { return "LIRFX";}

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

		ARM_Date		itsStartDate;
		ARM_Date		itsEndDate;
		ARM_Date		itsRedemptionResetDate;
		int				itsNbFX;
		ARM_Currency	itsDomesticCcy;
		vector<ARM_Currency>	itsForeignCcy;
		ARM_Currency	itsFundingCcy;
		int				itsPayRec;
		int				itsCpnDayCount;
		int				itsCpnFreq;
		int				itsCpnResetGap;
		string			itsCpnResetCal;
		string			itsCpnPayCal;
		int				itsStubRule;
		int				itsCpnTiming;
		int				itsCpnIntRule;
		ARM_Curve	itsCpnNominalCv;
		ARM_Curve	itsFundNominalCv;
		ARM_Curve	itsFundSpreadCv;
		vector<ARM_Curve>	itsDomesticCpnCv;
		vector<ARM_Curve>	itsForeignCpnCv;
		ARM_Curve	itsMinCpnCv;
		ARM_Curve	itsMaxCpnCv;
		vector<ARM_Curve>	itsInitialFXCv;
		int				itsFundFreq;
		int				itsFundDayCount;
		int itsRedemptionType;
		int				itsRedemptionGap;
		ARM_GP_Vector	itsRedemptionStrike;
		ARM_Curve	itsFees;

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

