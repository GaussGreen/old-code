/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarncalculatorindian.h
 *
 *  \brief file for the TARN Calculator Indian
 *	\author  P. Lam
 *	\version 1.0
 *	\date June 2007
 */


#ifndef _INGPCALCULATORS_PRDKOSWITCHERCALCULATOR_H
#define _INGPCALCULATORS_PRDKOSWITCHERCALCULATOR_H

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
/// \class ARM_PRDKOSwitcherCalculator
/// \brief
///  Class that implements a "indian" TARN FX calculator
///-----------------------------------------------------------------------------
class ARM_PRDKOSwitcherCalculator : public ARM_HybridIRFXCalculator 
{
	public:
		
		static const string SwitcherColNamesTable[];
		static const int ProductToPriceColumns[];

		enum productToPriceAlias
		{
			SpotFwdPrice,
			IsAliveValue,
			CouponPrice,
			RealCouponPrice,
			ProbaValue,
			OptionPrice,
			NbProductsToPrice
		};

	public:

		enum SwitcherColAlias
		{
			ResetDate=0,
			PayDate,
			IT,
			Barrier,
			CallBarrier,
			Leverage,
			Notional,
			SpotFwd,
			IsAlive,
			Coupon,
			RealCoupon,
			Proba,
			Option
		};

		/// constructor for deal from term sheet
		ARM_PRDKOSwitcherCalculator( const ARM_Date& asOfDate,
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
									 const ARM_Curve& barrier,
									 const ARM_Curve& callBarrier,
									 const ARM_Curve& leverage,
									 int intermediatePrices,
									 const ARM_StringVector& columnsToPrice,
									 ARM_FixingSched* pastFixings = NULL);

		///copy constructor, assignment constructor, destructor
		ARM_PRDKOSwitcherCalculator( const ARM_PRDKOSwitcherCalculator& rhs );
		ASSIGN_OPERATOR(ARM_PRDKOSwitcherCalculator)
		~ARM_PRDKOSwitcherCalculator();

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
		virtual ARM_Object* Clone() const { return new ARM_PRDKOSwitcherCalculator(*this); }

		virtual string ExportShortName() const { return "LSWCH";}

	protected :

		ARM_Curve		itsBarrier;
		ARM_Curve		itsCallBarrier;
		ARM_Curve		itsLeverage;
		
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

