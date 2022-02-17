/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: infswopfactory.cpp,v $
 * Revision 1.1  2004/09/09 16:39:43  ebenhamou
 * Initial revision
 *
 */


/*! \file infswopfactory.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou N. Belgrade
 *	\version 1.0
 *	\date September 2004
 *
 *	\version 2.0
 *	\date March 2005
 */


#include "gpinflation/infswopvol.h"
#include "gpinflation/infswopfactory.h"
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////
/// 

ARM_VolCurve* ARM_InfVolComputation_Factory::GenerateSwopVolCurve(
	const ARM_Date&			asOfDate,
	ARM_InfBSModel*			infIRBSModel,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ComputationMethodType	method)
{
	switch( method )
	{
	case INF_SWOPT_STD:
		{
			ARM_InfVolComputation_Producer_Std volComputer(infIRBSModel);
			return volComputer.GenerateSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	case INF_SWOPT_EQUAL_WEIGHT:
		{
			ARM_InfVolComputation_Producer_EqualWeight volComputer(infIRBSModel);
			return volComputer.GenerateSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	case INF_SWOPT_DF_WEIGHT:
		{
			ARM_InfVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,false);
			return volComputer.GenerateSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	case INF_SWOPT_DF_WEIGHT_SQUARE:
		{
			ARM_InfVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,true);
			return volComputer.GenerateSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": unknown type." );
	}
}


ARM_VolCurve* ARM_InfVolComputation_Factory::GenerateOATSwopVolCurve(
	const ARM_Date&			asOfDate,
	ARM_InfBSModel*			infIRBSModel,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ComputationMethodType	method,
	double					coupon,
	long					choice)
{

	bool mode = true;
		switch( choice )
		{
		case 1:
			{
				mode = true;
			}
			break;
		case 0:
			{
				mode = false;
			}
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": unknown type." );
	}

	
	switch( method )
	{
	case INF_SWOPT_STD:
		{
			ARM_InfVolComputation_Producer_Std volComputer(infIRBSModel);
			return volComputer.GenerateOATSwopVolCurve( asOfDate, tenors, expiries, coupon, mode );
		}
		break;
	case INF_SWOPT_EQUAL_WEIGHT:
		{
			ARM_InfOATVolComputation_Producer_EqualWeight volComputer(infIRBSModel);
			return volComputer.GenerateOATSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	case INF_SWOPT_DF_WEIGHT:
		{
			ARM_InfOATVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,false);
			return volComputer.GenerateOATSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	case INF_SWOPT_DF_WEIGHT_SQUARE:
		{
			ARM_InfOATVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,true);
			return volComputer.GenerateOATSwopVolCurve( asOfDate, tenors, expiries );
		}
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": unknown type swaption vol curve." );
	}

}

ARM_VolCube* ARM_InfVolCubeComputation_Factory::GenerateSwopVolCube(
	const ARM_Date&			asOfDate,
	ARM_InfBSModel*			infIRBSModel,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ARM_GP_Vector*			smiledTenors,
	ARM_GP_Vector*			strikes,
	ARM_InfVolComputation_Factory::ComputationMethodType	method)
{
	switch( method )
	{
	case ARM_InfVolComputation_Factory::INF_SWOPT_STD:
		{
			ARM_InfVolComputation_Producer_Std volComputer(infIRBSModel);
			return volComputer.GenerateSwopVolCube( asOfDate, tenors, expiries, smiledTenors, strikes );
		}
		break;
	case ARM_InfVolComputation_Factory::INF_SWOPT_EQUAL_WEIGHT:
		{
			ARM_InfVolComputation_Producer_EqualWeight volComputer(infIRBSModel);
			return volComputer.GenerateSwopVolCube( asOfDate, tenors, expiries, smiledTenors, strikes );
		}
		break;
	case ARM_InfVolComputation_Factory::INF_SWOPT_DF_WEIGHT:
		{
			ARM_InfVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,false);
			return volComputer.GenerateSwopVolCube( asOfDate, tenors, expiries, smiledTenors, strikes );
		}
		break;
	case ARM_InfVolComputation_Factory::INF_SWOPT_DF_WEIGHT_SQUARE:
		{
			ARM_InfVolComputation_Producer_SimpleWeight volComputer(infIRBSModel,true);
			return volComputer.GenerateSwopVolCube( asOfDate, tenors, expiries, smiledTenors, strikes);
		}
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": unknown type for swaption vol cube." );
	}
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

