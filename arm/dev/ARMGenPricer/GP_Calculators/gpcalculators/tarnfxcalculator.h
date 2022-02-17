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
							 const ARM_Curve& target,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 ARM_GP_Vector& redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_StringVector& columnsToPrice,
							 const ARM_FXTARNPayoffType& payOffName = ARM_TARNFXPayoffType::TARNFX,
							 ARM_FixingSched* pastFixings = NULL);

		///copy constructor, assignment constructor, destructor
		ARM_TARNFXCalculator( const ARM_TARNFXCalculator& rhs );
		ASSIGN_OPERATOR(ARM_TARNFXCalculator)
		~ARM_TARNFXCalculator();

		// Create the constant manager
		ARM_CstManagerPtr CreateCstManager();

		/// Initialisation to 0 of all columns of the deal description that could be priced
		void InitPriceableColumns(vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

		virtual void CheckData();
		virtual void CheckMktData();

		virtual double PastLiborValue(ARM_DateStripPtr dateStrip, double& IT, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec);
		virtual void DoPastReset(ARM_DateStripPtr dateStrip, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec);

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
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

