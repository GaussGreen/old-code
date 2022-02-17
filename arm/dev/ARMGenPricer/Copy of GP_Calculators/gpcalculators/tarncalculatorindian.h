/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarncalculatorindian.h
 *
 *  \brief file for the TARN FX Calculator
 *	\author  A. Lekrafi
 *	\version 1.0
 *	\date August 2006
 */


#ifndef _INGPCALCULATORS_TARNCALCULATORINDIAN_H
#define _INGPCALCULATORS_TARNCALCULATORINDIAN_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
//#include "gencalculator.h"
#include "tarnfxcalculator.h"
#include "gpmodels/Mixture_FX.h"

#include "typedef.h"
#include "refvalue.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_2IRFXModel;
class ARM_1IRFXModel;
class ARM_NP1IRNFXModel;


///-----------------------------------------------------------------------------
/// \class ARM_TARNCalculatorIndian
/// \brief
///  Class that implements a "indian" TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_TARNCalculatorIndian : public ARM_TARNFXCalculator 
{
	public:
		
		static const string TARNIndianColNamesTable[];
		static const int ProductToPriceColumns[];

		enum productToPriceAlias
		{
			SpotFwdPrice,
			TwiceCallPutPrice,
			CallMinusCallPrice,
			PutMinusPutPrice,
			CallDigitalPrice,
			PutDigitalPrice,
			IsAliveValue,
			CouponPrice,
			PaidCouponPrice,
			RealCouponPrice,
			ProbaValue,
			TARNPrice,
			NbProductsToPrice
		};

	public:

		enum TARNIndianColAlias
		{
			ResetDate=0,
			PayDate,
			BarrierUp,
			BarrierDown,
			Strike,
			Notional,
			Target,
			Fees,
			SpotFwd,
			TwiceCallPut,
			CallMinusCall,
			PutMinusPut,
			Epsilon,
			CallDigital,
			PutDigital,
			IsAlive,
			Coupon,
			SumCoupon,
			PaidCoupon,
			RealCoupon,
			Proba,
			TARN
		};

		/// constructor for deal from term sheet
		ARM_TARNCalculatorIndian(const ARM_Date& asOfDate,
								 const ARM_Date& startDate,
								 const ARM_Date& endDate,
								 const ARM_Currency& DomCcy,
								 const vector<ARM_Currency>& ForCcy,
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
								 const ARM_Curve& barrierUp,
								 const ARM_Curve& barrierDown,
								 double target,
								 double epsilon,
								 const ARM_Curve& fees,
								 int intermediatePrices,
								 const ARM_IntVector& productsToPrice);

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
		ARM_TARNCalculatorIndian( const ARM_TARNCalculatorIndian& rhs );
		ASSIGN_OPERATOR(ARM_TARNCalculatorIndian)
		~ARM_TARNCalculatorIndian();

		void CheckData();
		void CheckMktData();

		/// Initialisation to 0 of all columns of the deal description that could be priced
		void InitPriceableColumns(vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

		/// TARN indian deal description creation functions
		virtual ARM_RowInfo ColumnNames() const;
		virtual ARM_StringVector ProductsToPriceColumnNames();
		virtual ARM_DateStripCombiner DatesStructure() const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_TARNCalculatorIndian(*this); }

		virtual string ExportShortName() const { return "LTRIN";}

	protected :
		ARM_Curve		itsBarrierUp;
		ARM_Curve       itsCpnNominalCv;
		ARM_Curve		itsBarrierDown;
		double			itsEpsilon;

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

