/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prcscalculator.h
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date June 2006
 */


#ifndef _INGPCALCULATORS_PRDKOCALCULATOR_H
#define _INGPCALCULATORS_PRDKOCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "prdcalculator.h"


CC_BEGIN_NAMESPACE( ARM )

///-----------------------------------------------------------------------------
/// \class ARM_PRDKOCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_PRDKOCalculator : public ARM_PRDCalculator 
{
    static const string PRDKOColNamesTable [];

public:

    enum PRDKOColAlias
    {
		Barrier,
        FwdFX,
		PRDKO,
        PRDKOption,
    };
	/// constructor for deal from term sheet
	ARM_PRDKOCalculator(const ARM_Date& asofDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& switchDate,
		const ARM_Date& endDate,
		const ARM_Currency& DomCcy,
		const ARM_Currency& ForCcy,
		const ARM_Currency& FundCcy,
		int cpnDayCount,
		int cpnFreq,
		int FxResetGap,
		int stubRule,
		int resetTiming,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& domesticCpn,
		const ARM_Curve& foreignCpn,
		const ARM_Curve& MinCpn,
		const ARM_Curve& MaxCpn,
		const ARM_Curve& initialFx,
		const ARM_Curve& barrier,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,		
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		int redemptionGap,
		double redemptionStrike,		
		const ARM_RedemptionType& redemptionType    = ARM_PRCSRedemptionType::standard,
		const ARM_StringVector& productsToPrice		= ARM_StringVector(1,PRDKOColNamesTable[PRDKOption]),
		bool fxLocalModelFlag						= true,
		bool basisIRCalibFlag					    = true,
		const ARM_Curve& fundlevrage                = ARM_FlatCurve(1));

	///copy constructor, assignment constructor, destructor
	ARM_PRDKOCalculator( const ARM_PRDKOCalculator& rhs );
	ASSIGN_OPERATOR(ARM_PRDKOCalculator)
	~ARM_PRDKOCalculator(){};

    /// Initialisation to 0 of all columns of the deal description that could be priced
    void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	virtual ARM_StdPortfolioPtr CreateFxOption();
	virtual void UpdateFxOption(size_t eventIdx, ARM_Option*& fxOption);

    /// PRDC deal description creation functions
	virtual ARM_RowInfo ColumnNames() const;
	void CheckData();
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_PRDKOCalculator(*this); }

	virtual string ExportShortName() const { return "LPRDK";}

private:

	/// PRDC financial datas
	ARM_GP_Vector	itsvBarrier;
	size_t          itsNbCallable;
	ARM_BoolVector	itsvCpnIsCallable;

	/// Convert ARM_Curves into relevant vectors
	void ComputeProductVectorsFromCurves(const ARM_Curve& barrier );

    //virtual void Calibrate();
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

