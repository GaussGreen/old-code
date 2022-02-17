/*!
 *
 * Copyright (c) NATIXIS July 2007 Paris
 *
 *	\file FxSnowBallKoCalculator.h
 *
 *  \brief file for the IR FX Snwo Ball Ko Calculator
 *	\author N. Belgrade
 *	\version 1.0
 *	\date July 2007
 */

#ifndef _INGPCALCULATORS_FXSNOWBALLKOCALCULATOR_H
#define _INGPCALCULATORS_FXSNOWBALLKOCALCULATOR_H

//#include "hybridirfxcalculator.h"


#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_2IRFXModel;
class ARM_1IRFXModel;


///-----------------------------------------------------------------------------
/// \class ARM_FxSnowBallKoCalculator
/// \brief
///  Class that implements a "KO" FX Snow Ball calculator
///-----------------------------------------------------------------------------

class ARM_FxSnowBallKoCalculator : public ARM_HybridIRFXCalculator
{
	public :
		static const string fxSbKoColNamesTable[];
		static const int ProductToPriceColumns[];
	
		enum productToPriceAlias
		{
			SpotFwdPrice,
			OptionPrice,
			IsAliveValue,
			ProbaValue,
			CouponPrice,
			RealCouponPrice,
			FundingPrice,
			RealFundingPrice,
			FXSBKOPrice,
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
			Notional,
			Strike,
			Barrier,
			Leverage,
			Coeff,
			FixedFunding,
			Offset,
			Width,
			FxSpot,
			IsAlive,
			Option,
			Coupon,
			RealCoupon,
			Funding,
			DiscountFunding,
			RealFunding,
			DF,
			Proba,
			FXSBKO
		};

		enum ModelType
		{
			Model1IRFX,
			Model2IRFX,
		};

		// Constructor from term sheet
		ARM_FxSnowBallKoCalculator(const ARM_Date& asOfDate,
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
								 const ARM_Curve& cpnNominal,
								 int fundFreq,
								 int fundDayCount,
								 const ARM_Curve& fundNominal,
								 const ARM_Curve& fundSpread,
								 const ARM_Curve& strike,
								 const ARM_Curve& barrier,
								 const ARM_Curve& leverage,
								 const ARM_Curve& coeff,
								 int intermediatePrices,
								 const ARM_StringVector& columnsToPrice,
								 int optionType = K_CALL,
								 bool isPerform = false,
								 ARM_FixingSched* pastFixings = NULL);

		// Copy constructor, assignment constructor, destructor
		ARM_FxSnowBallKoCalculator (const ARM_FxSnowBallKoCalculator& rhs );
		ASSIGN_OPERATOR(ARM_FxSnowBallKoCalculator);
		~ARM_FxSnowBallKoCalculator() {};

		void CheckData();
		void CheckMktData();

		// Initialisation to 0 of all columns of the deal description that could be priced
		void InitPriceableColumns(vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

		// Fx Snow Ball Ko deal description creation functions
		virtual ARM_RowInfo ColumnNames() const;
		virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
		virtual string toString(const string& indent="",const string& nextIndent="") const;

		/// Standard ARM support
		virtual ARM_Object* Clone() const { return new ARM_FxSnowBallKoCalculator(*this); }

		virtual string ExportShortName() const { return "LSBKO"; }

	protected :
		int			itsOptionType;
		bool		itsIsPerform;
		ARM_Curve	itsStrikeCv;
		ARM_Curve	itsBarrierCv;
		ARM_Curve	itsLeverageCv;
		ARM_Curve	itsCoeffCv;
};

CC_END_NAMESPACE()

#endif