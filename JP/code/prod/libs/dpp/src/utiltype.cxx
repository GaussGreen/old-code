/************************************************************************
 * Module:      
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/crxdpp/src/utiltype.cxx,v 1.3 2005/06/30 13:06:52 dliu Exp $
 ************************************************************************/
#define	_dpp_SRC
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>		// getenv

#include "ktypes.h"
#include "kutilios.h"

extern	"C" {
#include "cgeneral.h"
#include "cerror.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "ratelist.h"		/* TRateList */
#include "tcurve.h"
#include "date_sup.h"
#include "stub.h"

#include "drlio.h"		/* FScanStruct */
#include "drlstr.h"		/* StringLineVarScan */
#include "drltime.h"		/* TDatePrint */
#include "drloptio.h"		/* Black */
};


//---------------------------------------------------------------

KDateInterval::KDateInterval(const char *str)
{
static	char	routine[] = "KDateInterval(const char *str)";

	TDateInterval	interval;
	int		dateAdjType;
	char		holidayFile[64];
	long		badDayConv;
	char		*p, *q, buf[128];

    try {

	strncpy(buf, str, sizeof(buf));


	if ((p = strstr(buf, "-B")) != NULL) {
		dateAdjType = GTO_DATE_ADJ_TYPE_BUSINESS;

		*p = '\0';
		p += 2;
		switch (*p) {
		case 'F': 
			badDayConv = GTO_BAD_DAY_FOLLOW;
			break;
		case 'P': 
			badDayConv = GTO_BAD_DAY_PREVIOUS;
			break;
		case 'M': 
			badDayConv = GTO_BAD_DAY_MODIFIED;
			break;
		case 'N': 
			badDayConv = GTO_BAD_DAY_NONE;
			break;
		default:
			throw KFailure("Failed scanning "
				"bus day conv in `%s'.\n", p);
		}

		// Calendar
		q = holidayFile;
		++p;
		if (*p == '-') {
			++p;
			strncpy(q, p, 3);
			p += 3;
			q += 3;
			*q = '\0';
		} else {
			strcpy(holidayFile, "NONE");
		}
	} else {
		dateAdjType = GTO_DATE_ADJ_TYPE_CALENDAR;
		badDayConv = GTO_BAD_DAY_NONE;
		strcpy(holidayFile, "NONE");
	}


	IF_FAILED_THROW( DrlTDateIntervalScan(buf, &interval));


	mDateAdjIntvl.interval = interval;
	mDateAdjIntvl.isBusDays = dateAdjType;
	mDateAdjIntvl.holidayFile = NULL;
	SetHoliday(holidayFile);
	mDateAdjIntvl.badDayConv = badDayConv;


    }
    catch (KFailure) {
	throw KFailure("%s: failed scanning KDateInterval in `%s'.\n",
		routine, buf);
    }
}


//---------------------------------------------------------------

bool KDateInterval::operator==(const KDateInterval& intvl) const
{
	TDateInterval	iv1, iv2;

	iv1 = mDateAdjIntvl.interval;
	iv2 = intvl.mDateAdjIntvl.interval;


	if (mDateAdjIntvl.isBusDays != intvl.mDateAdjIntvl.isBusDays)
		return (false);

	if (((mDateAdjIntvl.holidayFile == NULL) &&
	     (intvl.mDateAdjIntvl.holidayFile != NULL)) ||
	    ((mDateAdjIntvl.holidayFile != NULL) &&
	     (intvl.mDateAdjIntvl.holidayFile == NULL)))
		return (false);

	if ((mDateAdjIntvl.holidayFile != NULL) &&
	    (intvl.mDateAdjIntvl.holidayFile != NULL)) {
		if (strcmp(mDateAdjIntvl.holidayFile,
			   intvl.mDateAdjIntvl.holidayFile))
				return (false);
	}

	if (!GtoIntervalEqualsInterval(&iv1, &iv2))
		return (false);

	return (true);
}



//---------------------------------------------------------------

bool KDateInterval::operator< (const KDateInterval& intvl) const
{
	double	x1, x2;
	TDateInterval	iv1, iv2;

	iv1 = mDateAdjIntvl.interval;
	IF_FAILED_THROW( GtoDateIntervalToYears(
		&iv1, &x1));

	iv2 = intvl.mDateAdjIntvl.interval;
	IF_FAILED_THROW( GtoDateIntervalToYears(
		&iv2, &x2));

	return(x1 < x2);
}


//---------------------------------------------------------------

istream& operator>> (istream& is, KDateInterval &intvl)
{

	char		buf[128];

	strcpy(buf, getString(is, "KDateInterval"));
	intvl = KDateInterval(buf);

	return(is);
}


//---------------------------------------------------------------

ostream& operator<< (ostream& os, const KDateInterval &intvl)
{
	os << DrlTDateIntervalPrint(NULL,
		intvl.mDateAdjIntvl.interval);

	if (intvl.mDateAdjIntvl.isBusDays) {
	    switch (intvl.mDateAdjIntvl.badDayConv) {
	    case GTO_BAD_DAY_FOLLOW:
		os << "-BF";
		break;
	    case GTO_BAD_DAY_PREVIOUS:
		os << "-BP";
		break;
	    case GTO_BAD_DAY_MODIFIED:
		os << "-BM";
		break;
	    case GTO_BAD_DAY_NONE:
		os << "-BN";
		break;
	    default:
		os << "-ERROR";
		break;
	    }
	}

	return (os);
}


//---------------------------------------------------------------


const char*
KDateInterval::Str() const
{
#undef	MAX_IDX
#define	MAX_IDX	8
static	char	tmp[MAX_IDX][64] ;
static	int	tmpIdx=0;
	char	*s;

	s = tmp[tmpIdx];
	tmpIdx++;
	if (tmpIdx > MAX_IDX-1) tmpIdx=0;

	strcpy(s, DrlTDateIntervalPrint(NULL, mDateAdjIntvl.interval));
	if (mDateAdjIntvl.isBusDays) {
		switch (mDateAdjIntvl.badDayConv) {
		case GTO_BAD_DAY_FOLLOW:
			strcat(s, "-BF");
			break;
		case GTO_BAD_DAY_PREVIOUS:
			strcat(s, "-BP");
			break;
		case GTO_BAD_DAY_MODIFIED:
			strcat(s, "-BM");
			break;
		case GTO_BAD_DAY_NONE:
			strcat(s, "-BN");
			break;
		default:
			strcat(s, "-ERROR");
			break;
		}
		strcat(s, "-");
		if (mDateAdjIntvl.holidayFile) {
			strcat(s, mDateAdjIntvl.holidayFile);
		} else {
			strcat(s, "<null>");
		}
		strcat(s, "");
	}
	return (s);
}





//---------------------------------------------------------------


KDateInterval
getKDateInterval(istream& is, const char *debugName)
{
	const char	*p = getString(is, debugName);
	KDateInterval	iv;
	DppReadFromString(p, iv);
	return(iv);
}




//---------------------------------------------------------------


void
DppHolidayLoadFromDisk(
	const char *currencies,	// (I) list of cur separated by ":"
	const char *dir)	// (I) directory (or NULL)
{
static	char	routine[] = "DppLoadCalendars";
static	char	envvar[] = "DR_HOLIDAY_DIR";

    try {
	char	pathname[256],
		curname[32],
		buf[256], *p, *q, *r;

	if (currencies == NULL) {
		/* Nothing to do */
		if ((p = getenv(envvar)) == NULL)
			return;
		strncpy(buf, p, sizeof(buf));

		if ((q = strchr(buf, ':')) == NULL) {
		    throw KFailure("%s: bad format for "
			"environment variable %s (`%s').\n",
			routine, envvar, buf);
		}
		*q = '\0';
		q++;

		DppHolidayLoadFromDisk(q, buf);

	} else {

		strncpy(buf, currencies, sizeof(buf));
		p = &buf[0];

		while ((q = strtok(p, ":")) != NULL) {
			strncpy(curname, q, sizeof(curname));

			for (r = curname; *r; tolower(*r++));
			sprintf(pathname, "%s%sholiday.%s",
				(dir ? dir : ""),
				(dir ? "/" : ""),
				curname);

			for (r = curname; *r; toupper(*r++));

	
			IF_FAILED_THROW( GtoHolidayLoadFromDisk(
				curname,
				pathname));

			if (p != NULL) p = NULL;
		}
	}
    }
    catch (KFailure) {
	throw KFailure("%s: failed (cur=`%s',dir=`%s')\n",
		routine, (currencies ? currencies : "(NULL)"),
		(dir ? dir : "(NULL)"));
    }
}






//---------------------------------------------------------------
//
//

KKnockIO	
getKnockIO(istream& is, const char *debugName)
{
	const char	*s = getString(is);
	KKnockIO x;
    try {
	switch (toupper(s[0])) {
	case 'I':
		x = CRX_KNOCK_IN;
		break;
	case 'O':
		x = CRX_KNOCK_OUT;
		break;
	case 'N':
		x = CRX_NONE;
		break;
	default:
		throw KFailure();
	}
	return (x);
    }
    catch (...) {
	throw KFailure("Failed scanning %s.\n", debugName);
    }
}




istream& operator>>(istream& is, KKnockIO &x)
{
static	char	routine[] = "operator>>(istream&, KKnockIO&)";
	const char	*s = getString(is);
	switch (toupper(s[0])) {
	case 'I':
		x = CRX_KNOCK_IN;
		break;
	case 'O':
		x = CRX_KNOCK_OUT;
		break;
	case 'N':
		x = CRX_NONE;
		break;
	default:
		throw KFailure("%s: can't read `%s'\n", 
			routine, s);
	}
	return (is);
}


ostream& operator<<(ostream& os, const KKnockIO &x)
{
	switch (x) {
	case CRX_KNOCK_IN: os << "KNOCK IN "; break;
	case CRX_KNOCK_OUT:  os << "KNOCK OUT "; break;
	case CRX_NONE:  os << "NONE "; break;
	}
	return(os);
}




//---------------------------------------------------------------
//
//

KObsType	
getKObsType(istream& is, const char *debugName)
{
	const char	*s = getString(is);
	KObsType x;
    try {
	switch (toupper(s[0])) {
	case 'D':
		x = DATES;
		break;
	case 'F':
		x = FREQ;
		break;
	case 'N':
		x = NODES;
		break;
	case 'C':
		x = CONTINUOUS;
		break;
	default:
		throw KFailure();
	}
	return (x);
    }
    catch (...) {
	throw KFailure("Failed scanning %s.\n", debugName);
    }
}




istream& operator>>(istream& is, KObsType &x)
{
static	char	routine[] = "operator>>(istream&, KObsType&)";
	const char	*s = getString(is);
	switch (toupper(s[0])) {
	case 'D':
		x = DATES;
		break;
	case 'F':
		x = FREQ;
		break;
	case 'N':
		x = NODES;
		break;
	case 'C':
		x = CONTINUOUS;
		break;
	default:
		throw KFailure("%s: can't read `%s'\n", 
			routine, s);
	}
	return (is);
}


ostream& operator<<(ostream& os, const KObsType &x)
{
	switch (x) {
	case DATES: os << "DATES "; break;
	case FREQ:  os << "FREQ "; break;
	case NODES:  os << "NODES "; break;
	case CONTINUOUS:  os << "CONTINUOUS "; break;
	}
	return(os);
}



//---------------------------------------------------------------
//
//

KSmooth	
getSmooth(istream& is, const char *debugName)
{
	const char	*s = getString(is);
	KSmooth x;
    try {
	switch (toupper(s[0])) {
	case 'D':
		x = DOUBLE_SMOOTH;
		break;
	case 'S':
		x = SINGLE_SMOOTH;
		break;
	case 'N':
		x = NO_SMOOTH;
		break;
	default:
		throw KFailure();
	}
	return (x);
    }
    catch (...) {
	throw KFailure("Failed scanning %s.\n", debugName);
    }
}




istream& operator>>(istream& is, KSmooth &x)
{
static	char	routine[] = "operator>>(istream&, KSmooth&)";
	const char	*s = getString(is);
	switch (toupper(s[0])) {
	case 'D':
		x = DOUBLE_SMOOTH;
		break;
	case 'S':
		x = SINGLE_SMOOTH;
		break;
	case 'N':
		x = NO_SMOOTH;
		break;
	default:
		throw KFailure("%s: can't read `%s'\n", 
			routine, s);
	}
	return (is);
}


ostream& operator<<(ostream& os, const KSmooth &x)
{
	switch (x) {
	case DOUBLE_SMOOTH: os << "DOUBLE "; break;
	case SINGLE_SMOOTH:  os << "SINGLE "; break;
	case NO_SMOOTH:  os << "NO "; break;
	}
	return(os);
}



//---------------------------------------------------------------



KDayCc::KDayCc(const char *s)
{
	char buf[512];
	strncpy(buf, s, sizeof(buf));

    char *pos = NULL;

    // Strip out the 'D' if default accrual, e.g. something like "ACT/360D"
    //
    pos = strpbrk(buf, "D");

    if (pos != NULL)  /* check for default accrual */
        *pos = '\0';    /* remove the 'D'            */


	IF_FAILED_THROW( GtoStringToDayCountConv(
		buf,
		&mDayCc));
    
    // Flip the sign for indication of default accrual
    //
    if (pos != NULL)
        mDayCc = -mDayCc;

}



istream&
KDayCc::Get(istream& is, const char *debugName)
{
	char buf[512];
	strcpy(buf, getString(is, debugName));

	IF_FAILED_THROW( GtoStringToDayCountConv(
		buf,
		&mDayCc));

	return(is);
}



ostream& operator<<(ostream& os, const KDayCc &dcc)
{
        if (!dcc.isNegative())
           os << GtoFormatDayCountConv(dcc.mDayCc);
        else    // default accrual
           os << GtoFormatDayCountConv(-dcc.mDayCc) << "D";
    
	return(os);
}



double	DayCountFraction(TDate startDate, TDate endDate, KDayCc &dcc)
{
	double	dayCnt;

	IF_FAILED_THROW( GtoDayCountFraction(
		startDate,
		endDate,
		labs((long) dcc),    // default accrual
		&dayCnt));


	return(dayCnt);
}



//---------------------------------------------------------------


KStubConv::KStubConv(const char *s)
{
	char buf[512];
	strncpy(buf, s, sizeof(buf));

	if (toupper(buf[0]) == 'P')
		mStubConv = KV_STUB_PAR;
	else if (toupper(buf[0]) == 'O')
		mStubConv = KV_STUB_SPOT_SIMPLE;
	else
		IF_FAILED_THROW( GtoStringToStubType(
					buf,
					&mStubConv));
}


istream&
KStubConv::Get(istream& is, const char *debugName)
{
	char buf[512];
	strcpy(buf, getString(is, debugName));

	if (toupper(buf[0]) == 'P')
		mStubConv = KV_STUB_PAR;
	else if (toupper(buf[0]) == 'O')
		mStubConv = KV_STUB_SPOT_SIMPLE;
	else
		IF_FAILED_THROW( GtoStringToStubType(
					buf,
					&mStubConv));
	return(is);
}


ostream& operator<<(ostream& os, const KStubConv &stubConv)
{
	if (stubConv.mStubConv == KV_STUB_PAR)
		os << "PAR_STUB" << endl;
	else if (stubConv.mStubConv == KV_STUB_SPOT_SIMPLE)
		os << "ON_SPOT_SIMPLE" << endl;
	else
		os << GtoFormatStubType(stubConv.mStubConv);
	return(os);
}















//---------------------------------------------------------------




#ifdef	_SKIP




#define GTO_ACT_365       1L               /* Actual/365 */
#define GTO_ACT_365F      2L               /* Actual/365 Fixed */
#define GTO_ACT_360       3L               /* Actual/360 */
#define GTO_B30_360       4L               /* 30/360 */
#define GTO_B30E_360      5L               /* 30E/360 */
#define GTO_ACT_365FJ     6L               /* GTO_ACT_365FJ is maintained for
                                            * for backward compatibility even though
                                            * it is identical to the NL/365 convention.
                                            * Otherwise, object wrapper functions will
                                            * not work correctly if GTO_ACT_365FJ=NL_365.
                                            */
#define GTO_B30E_360I     7L
#define GTO_ACT_ACT      GTO_ACT_365
#define GTO_B30_360_FIXED 8L               /* For bond coupon payments */
#define GTO_B30EP_360     9L               /* 30E+/360 */
#define GTO_ACT_ACT_FRF   10L              /* For French TAM  */
#define GTO_NL_365        11L               /* Actual/365 Fixed w/No leap day */






istream& operator>>(istream& os, KDayCounyConv &dcc);


ostream& operator<<(ostream& os, const KDayCountConv &dcc);



#endif

