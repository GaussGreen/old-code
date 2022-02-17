/*----------------------------------------------------------------------------
 *
 *
/*! \file inflcapfloor.h
 *
 *  \brief Vanilla inflation cap and floor
 *	\Copyright (c) NATIXIS July 2007 Paris
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date July 2007
 */

/*----------------------------------------------------------------------------*/

#ifndef _INGPINFLATION_INFCAPFLOOR__H
#define _INGPINFLATION_INFCAPFLOOR__H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

/// inflation
#include "infleg_.h"

/// kernel
#include <inst/capfloor.h>
#include <glob/dates.h>

CC_USING_NS(std,map)



CC_BEGIN_NAMESPACE( ARM )



class InfCapFloor : public InfLeg/*public Instrument, public ARM_CapFloor*/
{
	public:
		InfCapFloor(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow, 
						int capFloor, double strike,
						//for backward compatibility : to remove
						int interpType,	double firstReset );
		InfCapFloor(	const InfCapFloor& rhs );
		InfCapFloor&	operator=( const InfCapFloor& rhs );
		virtual ~InfCapFloor(){}

		//virtual void View(char* id = NULL, FILE* ficOut = NULL);

		virtual void CptCashFlowValues();
		virtual void performCalculation() ;

		void StoreVol( double vol, double volLookupStrike, double strike, double timeToStart, 
			double tenor, int k );

	protected :
		int itsCapFloor;
		double itsStrike ;
		ARM_Matrix* itsVolInfo;
		ARM_GP_Vector itsCFValues;
	private:
		virtual void SetModel(ARM_CountedPtr<ARM_Model> model)  ;
		virtual void ComputePrice() ;


};
/*
	TODO : gicler ces merdes !!
*/


class InfCapFloorContext : public ARM_Object
{
private:
	ARM_Date itsNumDate;
 	ARM_Date itsDenomDate;
public:
	InfCapFloorContext(
		const ARM_Date& numDate,
 		const ARM_Date& denomDate)
	:
		itsNumDate( numDate ),
 		itsDenomDate( denomDate )
	{}
	
	inline const ARM_Date& GetNumDate() const {	return itsNumDate; }
	inline const ARM_Date& GetDenomDate() const { return itsDenomDate; }
};

class StoreVolInfoWithNewCF 
{
public:
	StoreVolInfoWithNewCF( double strike, double timeToStart, double tenor, 
		int capletNumber, ARM_CountedPtr<InfCapFloor> infCapFloor )
	: itsStrike( strike ), itsTimeToStart( timeToStart ), 
	itsTenor( tenor ), itsCapletNumber( capletNumber ),	itsInfCapFloor( infCapFloor )
	{}
	
	virtual void Store( double* data )
	{
		double vol				= data[0];
		double volLookupStrike	= data[1];
		itsInfCapFloor->StoreVol( vol, volLookupStrike, itsStrike, itsTimeToStart, itsTenor, itsCapletNumber );
	}
	ARM_CountedPtr<InfCapFloor> GetInfCapFloor() const { return itsInfCapFloor; }
private:
	double	itsStrike;
	double	itsTimeToStart;
	double	itsTenor;
	int		itsCapletNumber;
	ARM_CountedPtr<InfCapFloor>itsInfCapFloor;
};
/*
	fin des merdes
*/

CC_END_NAMESPACE()

#endif
