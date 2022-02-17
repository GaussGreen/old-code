/*!
 *
 * Copyright (c) CDC IXIS CM February 2005 Paris
 *
 *	\file vanillaargnumeric.h
 *
 *  \brief vanilla spreadoption
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date February 2005
 */


#ifndef _INGPCALIB_VANILLAARGNUMERIC_H
#define _INGPCALIB_VANILLAARGNUMERIC_H

#include "gpbase/env.h"
#include "gpbase/valuetype.h"
#include "vanillaarg.h"
#include "gpcalib/typedef.h"
#include "gpinfra/nummethod.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStripCombiner;


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaEqOption
/// \brief common vanilla arg for fx option
///////////////////////////////////////////////////////////////
struct ARM_VanillaArgNumeric: public ARM_VanillaArg
{
	ARM_VanillaArgNumeric(const string& curveName,
		double evalTime,
		int CallPut,
		double expiryTime);

	ARM_VanillaArgNumeric(const ARM_VanillaArgNumeric& arg);
    ARM_VanillaArgNumeric& operator=(const ARM_VanillaArgNumeric& rhs);

	virtual ~ARM_VanillaArgNumeric(){};
	virtual double Price(ARM_PricingModel* model ) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;

	/// standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	void CopyNoCleanUp(const ARM_VanillaArgNumeric& rhs);
	bool							isGenSecurityCreated;
	ARM_GenSecurityPtr				itsGenSecurity;
	ARM_NumMethodPtr				itsNumMethod;
	ARM_CstManagerPtr				itsCstManager;

public:
	/// Set accessors
	void SetGenSecurity( const ARM_GenSecurityPtr& genSecurity );
	void SetNumMethod( const ARM_NumMethodPtr& nummethod );
	
	/// Get accessors
	inline const ARM_GenSecurityPtr GetGenSecurity() const {return itsGenSecurity;}
	inline ARM_GenSecurityPtr GetGenSecurity() {return itsGenSecurity;}   // for non const functions
	inline const ARM_NumMethodPtr GetNumMethod() const {return itsNumMethod;}
	inline ARM_NumMethodPtr GetNumMethod() {return itsNumMethod;} // for non const functions
	inline void setConstantManager( const ARM_CstManagerPtr cstManager ) { itsCstManager = cstManager; }
	inline const ARM_CstManagerPtr& getConstantManager() const { return itsCstManager; }

	virtual ARM_NumMethod::GP_PricingDirection GetItsPricingDirection() const = 0;

	virtual ARM_RowInfo NumArgColumnNames() const = 0;
	virtual ARM_RowInfo NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates ) const = 0;
	virtual ARM_GP_VectorPtr DatesStructure() const = 0;
	virtual const ARM_CstManagerPtr securityCstManager() const = 0;
	virtual void CreateAndSetNumMethod(const ARM_Date& asOfDate);	
	
	void ValidateGenSecurity( const ARM_GenSecurityPtr& genSecurity );
	void FillRowInfo( const ARM_RowInfo& rowInfo,
		CC_NS(std,vector)<string>::iterator& textIter, 
		CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator& formatIter );

    /// Fct to build the generic security of the calculator
	void CreateAndSetGenSec(const string& curveName, const ARM_Date& asOfDate);	
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
