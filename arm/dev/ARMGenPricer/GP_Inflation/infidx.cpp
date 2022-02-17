/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 */

#include "gpinflation/infidx.h"
#include "gpinflation/infdata.h"
#include "gpbase/argconvdefault.h"

#include <string>
CC_USING_NS(std,string)

///	flag for interp
#include <armdef.h>

CC_BEGIN_NAMESPACE( ARM )



ARM_InfIdx::ARM_InfIdx( ARM_CountedPtr<ARM_Currency> ccy,
						ARM_CountedPtr<ARM_InfCurv> infCurve,
						const string& IndexName,
						const string& ResetLag,
						const string& DCFLag): 
	ARM_IRIndex(ARM_ArgConv_LgNameDayCount.GetNumber("30/360")),itsIndexName( IndexName ),itsInfFwdCurve(infCurve),
	itsResetLag( ResetLag == GETDEFAULTVALUESTR?InfData::GetResetLag( itsIndexName.c_str() ):ResetLag ), 
	itsDCFLag( DCFLag == GETDEFAULTVALUESTR?InfData::GetDCFLag( itsIndexName.c_str() ):DCFLag )
{
	ARM_Currency* tmpCcy ;

	if( ccy.IsNull() )
		tmpCcy = new ARM_Currency( InfData::GetCurrency( itsIndexName.c_str() ) );
	else
		tmpCcy = (ARM_Currency*)( ccy->Clone() );
	SetCurrencyUnit( tmpCcy ); 
	delete tmpCcy;
	SetName( ARM_INFIDX );
}

ARM_InfIdx::ARM_InfIdx(	const string& IndexName )
: ARM_IRIndex(ARM_ArgConv_LgNameDayCount.GetNumber("30/360"))
{
	itsIndexName	= IndexName ;
	SetCurrencyUnit( new ARM_Currency( InfData::GetCurrency( IndexName.c_str() ) ) );	
	SetName( ARM_INFIDX );
}

ARM_InfIdx::ARM_InfIdx( const ARM_InfIdx& rhs )
:
	ARM_IRIndex( rhs ), 
	itsIndexName(rhs.itsIndexName) ,
	itsInfFwdCurve( rhs.itsInfFwdCurve ),
	itsResetLag( rhs.itsResetLag ),
	itsDCFLag( rhs.itsDCFLag )
{}

ARM_InfIdx& ARM_InfIdx::operator =(const ARM_InfIdx& rhs )
{
	if( this != &rhs )
	{
		ARM_IRIndex::operator =(  rhs );
		itsInfFwdCurve		= rhs.itsInfFwdCurve, 
		itsResetLag		= rhs.itsResetLag;
		itsDCFLag		= rhs.itsDCFLag;
		itsIndexName	= rhs.itsIndexName ;
	}
	return *this;
}


ARM_Object* ARM_InfIdx::Clone(void)	
{ 
	return new ARM_InfIdx(*this);
}


void ARM_InfIdx::View(char* id, FILE* ficOut)
{
	/// to be consistent with other interface
	/// need to use fprintf
    FILE* fOut;
    char fOutName[40];
	
	/// do we have already a file opened?
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// printing of the field
    fprintf(fOut, "\n\n\t =====> Inflation Index\n\n");
	fprintf(fOut, "\t Index     : %s \t\t Currency    : %s\n", 
		itsIndexName.c_str(), GetCurrencyUnit()->GetCcyName() );
    if ( ficOut == NULL )
       fclose(fOut);
}


double ARM_InfIdx::FwdCPI( const ARM_Date& fixingDate )
{
	return itsInfFwdCurve->CPIInterpolate( fixingDate, fixingDate, itsInfFwdCurve->GetDailyInterpType() );

}


CC_END_NAMESPACE()

