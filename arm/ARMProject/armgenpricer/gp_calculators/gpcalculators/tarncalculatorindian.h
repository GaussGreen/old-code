/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarncalculatorindian.h
 *
 *  \brief file for the TARN Calculator Indian
 *	\author  P. Lam
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INGPCALCULATORS_TARNCALCULATORINDIAN_H
#define _INGPCALCULATORS_TARNCALCULATORINDIAN_H

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
/// \class ARM_TARNCalculatorIndian
/// \brief
///  Class that implements a "indian" TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_TARNCalculatorIndian : public ARM_HybridIRFXCalculator 
{
	public:
		
		static const string IndianColNamesTable[];
		static const int ProductToPriceColumns[];

		enum productToPriceAlias
		{
			SpotFwdPrice,
			ShortPrice,
			LongPrice,
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

		enum IndianColAlias
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
			Short,
			Long,
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

		enum IndianType
		{
			DownUp,
			Trigger,
			DigitalTrigger
		};

		/// constructor for deal from term sheet
		ARM_TARNCalculatorIndian(const ARM_Date& asOfDate,
								 const ARM_Date& startDate,
								 const ARM_Date& endDate,
								 const ARM_Currency& DomCcy,
								 const vector<ARM_Currency>& ForCcy,
								 const ARM_Currency& CpnCcy,
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
								 const ARM_Curve& strike,
								 const ARM_Curve& target,
								 double epsilon,
								 const ARM_Curve& fees,
								 int intermediatePrices,
								 const ARM_StringVector& columnsToPrice,
								 const int optionType,
								 const int indianType,
								 ARM_FixingSched* pastFixings = NULL);

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
		virtual ARM_DateStripCombiner DatesStructure() const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_TARNCalculatorIndian(*this); }

		virtual string ExportShortName() const { return "LTRIN";}

	protected :

		ARM_Curve		itsBarrierUp;
		ARM_Curve		itsBarrierDown;
		double			itsEpsilon;
		ARM_Curve		itsStrike;
		int				itsOptionType;
		int				itsIndianType;
		
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

