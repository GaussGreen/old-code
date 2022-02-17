/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/

#include "vpcashfl.h"

#include "kutilios.h"		// ios utilities

extern	"C" {
#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             // MAX 
#include "ldate.h"              // GtoDayCountFraction 
#include "stub.h"               // GTO_STUB_NONE 
#include "datelist.h"           // GtoNewDateList 
#include "date_sup.h"           // GtoFreq2DateInterval 
#include "tcurve.h"             // GtoDiscountDate 
#include "cerror.h"             // PgPrintErrMsg 
#include "date_sup.h"
#include "modlcons.h"		// 
#include "modltype.h"		// 
#include "swapton.h"		// GtoSwaptionModel() 
#include "swopt.h"		// GtoTrinomialCashFlow() 
#include "trinomv.h"		// GtoSwaptionSens() 
#include "optprop.h"

#include "zr2coup.h"		// Analytics C Library 
#include "zr2simp.h"		// Analytics C Library 

#include "fxddate.h"
#include "cashflow.h"
#include "tcurve.h"
#include "yield.h"

#include "swopblk.h"		// GtoSwaptionPVBlackVanilla() 

#include "drlio.h"		// FScanStruct 
#include "drlstr.h"		// StringLineVarScan 
#include "drloptio.h"		// Black 
#include "drltime.h"		// DrlTDateIntervalNil() 
#include "drlroot.h"

};



//---------------------------------------------------------------


KVPCashFlows::KVPCashFlows(
	const char *name,
	const KVector(TDate) payDates,
	const KVector(double) amounts,
	const char* discZcName)
	: KVPInstr(name)
{

	if (payDates.size() != amounts.size())
	    throw KFailure("KVPCashFlows::KVPCashFlows: "
		"payDates.size() != amounts.size().\n");

	mPayDates = payDates;
	mPayAmounts = amounts;
	mDiscZcName = String(discZcName);

	SetName(name);
}



//---------------------------------------------------------------

istream&
KVPCashFlows::Get(istream& is, int drw)
{

    try {
	if (drw) {
	    int	idx, numDates;

	    numDates = getInt(is, "KVPCashFlows::Get: number of dates.");

	    for (idx=0; idx<numDates; idx++) {
		mPayDates.insert(mPayDates.end(),
			getTDate(is, "KVPCashFlows::Get: pay date."));
		mPayAmounts.insert(mPayAmounts.end(),
			getDouble(is, "KVPCashFlows::Get: pay amounts."));
	    }

	} else {
		throw KFailure("KVPCashFlows::Get: format N/A.\n");
	}


	return(is);
    }
    catch(...) {
	throw KFailure("KVPCashFlows::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPCashFlows::Put(ostream& os, int indent) const
{
	int	idx, numDates;

    try {
	numDates = mPayDates.size();

	ASSERT_OR_THROW(mPayDates.size() == numDates);
	ASSERT_OR_THROW(mPayAmounts.size() == numDates);


	os << "NAME: `" << GetName() << "'" << endl;
	os << "NUMDATES: " << numDates << endl;

	for (idx=0; idx<numDates; idx++) {
		os << format(" %10s %18.6f\n",
			DrlTDatePrint(NULL, mPayDates[idx]),
			mPayAmounts[idx]);
	}
	os << "DISCOUNT CURVE: " << mDiscZcName << endl;

	// list dependencies
	this->KVPAtom::Put(os);

	return(os);
    }
    catch(...) {
	throw KFailure("KVPCashFlows::Get: failed.\n");
    }
}




//---------------------------------------------------------------

ostream&
KVPCashFlows::YacctionWrite(ostream& os, int indent)
{
	int	idx, numDates;

    try {

	if (GetWriteFlag())
	{
	    numDates = mPayDates.size();

	    ASSERT_OR_THROW(mPayDates.size() == numDates);
	    ASSERT_OR_THROW(mPayAmounts.size() == numDates);


	    os << GetName() << "=CASHFLOWS({" << endl;

	    for (idx=0; idx<numDates; idx++) {
		os << format(" %10s %18.6f\n",
			GtoFormatDate(mPayDates[idx]),
			mPayAmounts[idx]);
	    }
	    os << "\t}," << endl 
	       << "\t\"" << mDiscZcName << "\");"
	       << endl << endl;
	
	    WriteDone();
	}

	return(os);
    }
    catch(...) {
	throw KFailure("KVPCashFlows::Get: failed.\n");
    }
}





