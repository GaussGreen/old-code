/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file discretisationscheme.cpp
 *  \brief computes the inflation swaption volatility from year to year volatility
 *	\author  E Benhamou, N. Belgrade
 *	\version 1.0
 *	\date September 2004
 */

    
/*----------------------------------------------------------------------------*/

/*! \file infswopvol.h
 *
 *  \brief computation of model inmplied inflation swaption volatility
 *
 *	\author  E. Benhamou, N. Belgrade
 *	\version 1.0
 *	\date September 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_INFSWOPVOL_H
#define _INGPINFLATION_INFSWOPVOL_H

#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

/// kernel
#include <glob/armglob.h>

#include <string>
CC_USING_NS(std,string)

/// forward declaration
class ARM_Matrix;
class ARM_VolCurve;
class ARM_VolCube;
class ARM_CorrelManager;
class ARM_BSModel;
class ARM_Date;

/// ARM namespace
CC_BEGIN_NAMESPACE( ARM )

class ARM_InfCurv;
class ARM_InfBSModel;

/////////////////////////////////////////////////////////
/// common class
/// function are pure virtual
/// to force redefinition
/// static are protected for easy access by derived classes
/////////////////////////////////////////////////////////

class ARM_InfVolComputation_Producer : public ARM_Object
{
private:
	/// part to have some Default for the tenors and expiries
	static double DefaultTenors[];
	static size_t DefaultTenorsSize;
	static double DefaultExpiries[];
public: // FIXMEFRED: mig.vc8 (31/05/2007 10:30:25): someone access to DefaultExpiriesSize
	static size_t DefaultExpiriesSize;
	
/// protected for easy access
protected:
	ARM_InfBSModel*	itsInfBSModel;

	void GetDefaultDataIfMissing(
		ARM_GP_Vector*& tenors,
		ARM_GP_Vector*& expiries );

public:
	/// computes the inflation swaption volatility from 
	///		- inflation year on year inflation vol
	///		- inflation CPI fwd curve
	///		- IR vol
	///		- IR df
	/// computation of the whole vol curve
	/// there is no concept of strike!
/*	virtual ARM_VolCurve* GenerateSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries )=0;
*/

	/// constructors
	ARM_InfVolComputation_Producer( ARM_InfBSModel*	infIRBSModel );
	ARM_InfVolComputation_Producer( const ARM_InfVolComputation_Producer& rhs );
	ARM_InfVolComputation_Producer& operator=( const ARM_InfVolComputation_Producer& );
	virtual ~ARM_InfVolComputation_Producer();

	/// standard ARM_Object support
	virtual ARM_Object* Clone() = 0;
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	virtual string toString( const string& indent = "" ) = 0;

};


/////////////////////////////////////////////////////////
/// class ARM_InfVolComputation_Producer_Std
///		computes the inflation swaption volatility
///		according to standard vol swap vol FRA relationship
/////////////////////////////////////////////////////////

class ARM_InfVolComputation_Producer_Std : public  ARM_InfVolComputation_Producer
{
private:

	/// Standard relation as a analogous of vol swp-vol FRA relationship between the 
	/// year on year vols and the periodic annual inflation swap

	double VolYoY_to_VolSwp(
		/// fix leg
		ARM_GP_Vector* pFixLegDF, 
		ARM_GP_Vector* pFixLegIRYearFrac, 
		ARM_GP_Vector* pFixLegDFVol,
		//// inflation
		ARM_GP_Vector* pFwdCPIRatio,		/// this is supposed to be convexity corrected
		ARM_GP_Vector* pVol_YoY,
		ARM_GP_Vector* pInfLegDF,
		ARM_GP_Vector* pInfLegDFVol,
		ARM_Matrix* AvgFloatCPIIRCorrel,
		ARM_Matrix* AvgFloatFixCPIIRCorrel,
		ARM_Matrix* AvgFloatCPICPICorrel,
		ARM_Matrix* AvgFixIRIRCorrel,
		ARM_Matrix* AvgFloatIRIRCorrel,
		ARM_Matrix* AvgFloatFixIRIRCorrel );

	/// Standard relation as a analogous of vol swp-vol FRA relationship between the 
	/// zero coupon vols and the periodic OAT inflation swap

		double VolZC_to_VolOATSwp(
		/// fix leg
		ARM_GP_Vector* pFixLegDF,
		ARM_GP_Vector* pFixLegIRYearFrac,
		ARM_GP_Vector* pFixLegDFVol,
		ARM_GP_Vector* pFwdCPI,
		ARM_GP_Vector* pVol_ZC,
		ARM_GP_Vector* pInfLegDF,
		ARM_GP_Vector* pInfLegDFVol,
		ARM_Matrix* AvgFloatCPIIRCorrel,
		ARM_Matrix* AvgFloatFixCPIIRCorrel,
		ARM_Matrix* AvgFloatCPICPICorrel,
		ARM_Matrix* AvgFixIRIRCorrel,
		ARM_Matrix* AvgFloatIRIRCorrel,
		ARM_Matrix* AvgFloatFixIRIRCorrel,
		double coupon,
		bool mode);

		
public:
	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries );


	virtual ARM_VolCube* GenerateSwopVolCube( 
		const ARM_Date&			asOfDate,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ARM_GP_Vector*			smiledTenors,
		ARM_GP_Vector*			strikes);

	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateOATSwopVolCurve(
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries,
		double				coupon,
		bool				mode);

	virtual ARM_VolCurve* CompleteOATSwopVolCurve(
		const ARM_Date&		asOfDate,
		ARM_VolCurve* vol,
		double coupon);

	/// constructor
	ARM_InfVolComputation_Producer_Std(	ARM_InfBSModel*	infIRBSModel );
	ARM_InfVolComputation_Producer_Std( const ARM_InfVolComputation_Producer_Std& rhs );
	ARM_InfVolComputation_Producer_Std& operator=( const ARM_InfVolComputation_Producer_Std& );
	virtual ~ARM_InfVolComputation_Producer_Std();

	/// standard ARM_Object support
	virtual ARM_Object* Clone();
	virtual string toString( const string& indent = "" );
};



/////////////////////////////////////////////////////////
/// ARM_InfVolComputation_Producer_EqualWeight 
/// just computation using same weights
/////////////////////////////////////////////////////////

class ARM_InfVolComputation_Producer_EqualWeight : public  ARM_InfVolComputation_Producer
{
public:
	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries );

	virtual ARM_VolCube* GenerateSwopVolCube( 
		const ARM_Date&			asOfDate,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ARM_GP_Vector*			smiledTenors,
		ARM_GP_Vector*			strikes);

	/// constructor
	ARM_InfVolComputation_Producer_EqualWeight( ARM_InfBSModel*	infIRBSModel );
	ARM_InfVolComputation_Producer_EqualWeight( const ARM_InfVolComputation_Producer_EqualWeight& rhs );
	ARM_InfVolComputation_Producer_EqualWeight& operator=( const ARM_InfVolComputation_Producer_EqualWeight& rhs );
	virtual ~ARM_InfVolComputation_Producer_EqualWeight();

	/// standard ARM_Object support
	virtual ARM_Object* Clone();
	virtual string toString( const string& indent = "" );
};



/////////////////////////////////////////////////////////
/// ARM_InfVolComputation_Producer_EqualWeight 
/// just computation using same weights
/////////////////////////////////////////////////////////

class ARM_InfVolComputation_Producer_SimpleWeight : public  ARM_InfVolComputation_Producer
{
private:
	bool itsUsedSquareVersion;
public:
	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateSwopVolCurve( 
		const ARM_Date&			asOfDate,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries );

	virtual ARM_VolCube* GenerateSwopVolCube( 
		const ARM_Date&			asOfDate,
		ARM_GP_Vector*			tenors,
		ARM_GP_Vector*			expiries,
		ARM_GP_Vector*			smiledTenors,
		ARM_GP_Vector*			strikes);

	/// constructor
	ARM_InfVolComputation_Producer_SimpleWeight( ARM_InfBSModel* infIRBSModel, bool usedSquareVersion );
	ARM_InfVolComputation_Producer_SimpleWeight( const ARM_InfVolComputation_Producer_SimpleWeight& rhs );
	ARM_InfVolComputation_Producer_SimpleWeight& operator=( const ARM_InfVolComputation_Producer_SimpleWeight& );
	virtual ~ARM_InfVolComputation_Producer_SimpleWeight();

	/// standard ARM_Object support
	virtual ARM_Object* Clone();
	virtual string toString( const string& indent = "" );
};



class ARM_InfOATVolComputation_Producer_SimpleWeight : public  ARM_InfVolComputation_Producer
{
private:
	bool itsUsedSquareVersion;
public:
	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateOATSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries );

		
	/// constructor
	ARM_InfOATVolComputation_Producer_SimpleWeight( ARM_InfBSModel* infIRBSModel, bool usedSquareVersion );
	ARM_InfOATVolComputation_Producer_SimpleWeight( const ARM_InfOATVolComputation_Producer_SimpleWeight& rhs );
	ARM_InfOATVolComputation_Producer_SimpleWeight& operator=( const ARM_InfOATVolComputation_Producer_SimpleWeight& );
	virtual ~ARM_InfOATVolComputation_Producer_SimpleWeight();

	/// standard ARM_Object support
	virtual ARM_Object* Clone();
	virtual string toString( const string& indent = "" );
};


/////////////////////////////////////////////////////////
/// ARM_InfVolComputation_Producer_EqualWeight 
/// just computation using same weights
/////////////////////////////////////////////////////////

class ARM_InfOATVolComputation_Producer_EqualWeight : public  ARM_InfVolComputation_Producer
{
public:
	/// computation of the whole vol curve
	/// there is no concept of strike!
	virtual ARM_VolCurve* GenerateOATSwopVolCurve( 
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries );


	/// constructor
	ARM_InfOATVolComputation_Producer_EqualWeight( ARM_InfBSModel*	infIRBSModel );
	ARM_InfOATVolComputation_Producer_EqualWeight( const ARM_InfOATVolComputation_Producer_EqualWeight& rhs );
	ARM_InfOATVolComputation_Producer_EqualWeight& operator=( const ARM_InfOATVolComputation_Producer_EqualWeight& rhs );
	virtual ~ARM_InfOATVolComputation_Producer_EqualWeight();

	/// standard ARM_Object support
	virtual ARM_Object* Clone();
	virtual string toString( const string& indent = "" );
};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
