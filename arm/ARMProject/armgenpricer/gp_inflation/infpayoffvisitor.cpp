
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf PayOff Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


#include "gpinflation\infpayoffvisitor.h"
#include "gpinflation\infoptionspreadvisitor.h"
#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation\infdoubledigitalvisitor.h"
#include "gpinflation\infcorridorvisitor.h"

CC_BEGIN_NAMESPACE( ARM )


void ARM_InfPayOffVisitor::Visit(ARM_InfHybridPayOff & payOff){
		ARM_InfHybridPayOffValue* tmp = new ARM_InfHybridPayOffValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
}; 

void ARM_InfPayOffVisitor::Visit(ARM_InfHybridCap & payOff){
		ARM_InfHybridCapValue* tmp = new ARM_InfHybridCapValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
};

void ARM_InfPayOffVisitor::Visit(ARM_InfHybridDigit & payOff){
		ARM_InfHybridDigitValue* tmp = new ARM_InfHybridDigitValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
};

void ARM_InfPayOffVisitor::Visit(ARM_InfDoubleDigit & payOff){
		ARM_InfDoubleDigitValue* tmp = new ARM_InfDoubleDigitValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
};

void ARM_InfPayOffVisitor::Visit(ARM_InfCorridor & payOff){
		ARM_InfCorridorValue* tmp = new ARM_InfCorridorValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
};

CC_END_NAMESPACE()



/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















