/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include "kemmatr.h"

#include "kutilios.h"		// I/O on stream utilities


extern	"C" {
#include "date_sup.h"
#include "ldate.h"
#include "cerror.h"
#include "convert.h"
#include "interp.h"		/* GtoInterpRate() */

#include "drlvtype.h"
#include "drlstr.h"
#include "drlio.h"
#include "drlinter.h"
#include "drltime.h"
#include "drlmem.h"		/* DoubleMatrAlloc */
#include "drllineq.h"		/* DrlRealLinSysSvd */

#include "drlsmat.h"		/* Drl routines */
};


// Because we may frequency of -1

static	TSwaptionMatrix2D*
copyMatrixSafe(TSwaptionMatrix2D *mat)
{
	TSwaptionMatrix2D *newMat;
	long	freqSave = mat->swapPayFreq;

	mat->swapPayFreq = 2;
	newMat = DrlTSwaptionMatrix2DNewCopy(mat);
	if (newMat == NULL)
		throw KFailure();

	newMat->swapPayFreq = freqSave;
	mat->swapPayFreq = freqSave;
	return(newMat);
}



//---------------------------------------------------------------

KEmMatrix::KEmMatrix()
{
	mFormat = DEF;
	mMatrix = NULL;
}

//---------------------------------------------------------------

KEmMatrix::KEmMatrix(const KEmMatrix& argument)
{
	mFormat = DEF;
	mMatrix = NULL;
	if (argument.mMatrix)
		mMatrix = copyMatrixSafe(argument.mMatrix);
}

//---------------------------------------------------------------

KEmMatrix::KEmMatrix(TSwaptionMatrix2D* argument)
{
	mFormat = DEF;
	mMatrix = NULL;
	if (argument)
		mMatrix = copyMatrixSafe(argument);
}

//---------------------------------------------------------------

KEmMatrix::~KEmMatrix()
{
	cleanup();
}

//---------------------------------------------------------------

void
KEmMatrix::CheckValid()
{
static	char	routine[] = "KEmMatrix::CheckValid()";
	int	n;

	// check frequency 
	switch (Freq()) {
	case -1:
	case 0:
	case 1:
	case 2:
	case 4:
	case 12:
		break;
	default:
		throw KFailure("%s: bad frequency %d\n", routine, Freq());
	}

	if (TExp(0) < 0e0) {
		throw KFailure("%s: exp # %d (%lf) < 0\n",
				routine, 1, TExp(0));
	}
	for (n=1; n<=NExp()-1; n++) {
	    if (TExp(n-1) >= TExp(n)) {
		throw KFailure("%s: exp #%d (%lf) > exp #%d (%lf).\n",
		   routine, n, TExp(n-1), n+1, TExp(n));
	    }
	}

	if (TMat(0) < 0e0) {
	    throw KFailure("%s: mat # %d (%lf) < 0\n",
		routine, 1, TMat(0));
	}
	for (n=1; n<=NMat()-1; n++) {
	    if (TMat(n-1) >= TMat(n)) {
		throw KFailure("%s: mat #%d (%lf) > mat #%d (%lf).\n",
		   routine, n, TMat(n-1), n+1, TMat(n));
	    }
	}

	return;
}





//---------------------------------------------------------------
// Private function for memory cleanup.

int
KEmMatrix::cleanup()
{
	DrlTSwaptionMatrix2DFree(mMatrix);
	mMatrix = NULL;
	return(SUCCESS);
}


//---------------------------------------------------------------
// Assignement operator.


KEmMatrix&
KEmMatrix::operator=(const KEmMatrix& argument)
{
	cleanup();
	if (!argument.isempty()) {
		mMatrix = copyMatrixSafe(argument.mMatrix);
	}
	return(*this);
}


//---------------------------------------------------------------
// Assignement operator.


KEmMatrix&
KEmMatrix::operator=(TSwaptionMatrix2D* argument)
{
	cleanup();
	mMatrix = copyMatrixSafe(argument);
	CheckValid();
	return(*this);
}



//---------------------------------------------------------------
// Assignement operator.


KEmMatrix&
KEmMatrix::operator=(TCurve* curve)
{
	int	i;

	// cleanup 
	cleanup();

	// Creat new matrix 
	mMatrix  = DrlTSwaptionMatrix2DNew(
			FALSE, 	/* vertical */
			2,
			curve->fNumItems,
			NULL,
			1,
			NULL);
	ASSERT_OR_THROW(mMatrix != NULL);

	mMatrix->swapPayFreq = (long) curve->fBasis;
	mMatrix->diagonal = FALSE;

	TMat(0) = 1e0 / curve->fBasis;
	for (i=0; i<=NExp()-1; i++) {
		TExp(i) = (curve->fArray[i].fDate -
				curve->fBaseDate) / 365e0;
		Val(i,0) = curve->fArray[i].fRate;
	}

	return(*this);
}




//---------------------------------------------------------------

istream&
KEmMatrix::Get(istream& is)
{
static	char	routine[] = "operator>>(istream&, KEmMatrix&)";

const	char	*p;
	int	i, j, nMat, nExp, c;



    try {
	switch (mFormat) {
	case KEmMatrix::DEF:
	case KEmMatrix::CURM:
		p = getNextLine(is);
		ASSERT_OR_THROW(sscanf(p, "NUMBER_OF_MAT: %d", &nMat) == 1);

		p = getNextLine(is);
		ASSERT_OR_THROW(sscanf(p, "NUMBER_OF_EXP: %d", &nExp) == 1);

		// Creat new matrix 
		mMatrix  = DrlTSwaptionMatrix2DNew(
			FALSE, 	/* vertical */
			2,
			nExp,
			NULL,
			nMat,
			NULL);
		ASSERT_OR_THROW(mMatrix != NULL);


		p = getNextLine(is);
		ASSERT_OR_THROW(sscanf(p, "FREQ: %d", &c) == 1);
		mMatrix->swapPayFreq = (long) c;

		p = getNextLine(is);
		ASSERT_OR_THROW(sscanf(p, "DIAG: %d", &c) == 1);
		mMatrix->diagonal = (c == 0 ? FALSE : TRUE);

		if (mMatrix->diagonal) {
			throw KFailure("%s: not support diag matrix.\n",
				routine);
		}


		for (i=0; i<=nMat-1;i++) {
			mMatrix->table->dim2Values[i] =
				getDoubleInterval(is);
		}

		for (j=0; j<=nExp-1; j++) {
			mMatrix->table->dim1Values[j] =
				getDoubleInterval(is);
	        	for (i=0; i<=nMat-1;i++) {
		    	mMatrix->table->matrix->data[j][i] = getDouble(is);
			}
		}   
		break;

	case KEmMatrix::BVMAT: {
		// Base volatility matrix format (possibly IMM dates)
		// of the form
		// \begin{verbatim}
		// "Maturity"   "Exp\Und"  0.08330 0.25000 0.50000 ...
		// "SN"      "20-May-1998" 0.43000 0.45000 0.41000 ...
		// "IMM1"    "17-Jun-1998" 0.45000 0.42000 0.42000 ...
		// "IMM2"    "16-Sep-1998" 0.42000 0.43000 0.46000 ...
		// "IMM3"    "16-Dec-1998" 0.43000 0.44000 0.48000 ...
		// ...        ...          ...     ...     ...     ...
		// \end{verbatim}
		// The date {\tt baseDate} is used to convert intervals 
		// to years (assuming a 30/360 day count).
		// If {\tt baseDate} is 0, the reference date is assumed
		// to be the first date (SN).
		// Returns 0 iff successful. 


		TDate	baseDate;		// reference date
		int	freq = 2;		// frequency (overwritten!)

		char	nextLine[1024];
		TDate	expDate;
		int	nExp, nMat, idxE, idxM;
#define		nExpMax	128
#define		nMatMax	128
		double	tExp[nExpMax],
			tMat[nMatMax],
			values[nExpMax][nMatMax];

        	TDayCount dayCount = GTO_B30_360; // dcc for year fraction 

		// Read first line
		strcpy(nextLine, getNextLine(is, "First line"));
		istrstream  fls(nextLine, strlen(nextLine));

		getString(fls, "Maturity token");
		getString(fls, "Exp/Und token");

		for (nMat=0; !fls.eof(); nMat++) {
			tMat[nMat] = getDouble(fls, "time to mat");
			if (nMat >= nMatMax) {
				throw KFailure("%s: too many maturities.\n",
					routine);
			}
		}
		if (nMat <= 0)
			throw KFailure("%s: no maturities.\n", routine);


		// Read lines
		for(nExp=0; !is.eof(); nExp++) {
			strcpy(nextLine, getNextLine(is));
			istrstream  ls(nextLine, strlen(nextLine));

			// We must have 1st date to be SN 
			//
			const char*p = getString(ls, "exp interval");	// "interval"
			if (nExp == 0) {
			    if (strcmp(p, "SN")) 
				throw KFailure("%s: first interval is not"						"SN (but %s).\n", routine, p);
			}

			expDate = getTDate(ls, "exp date");
			if (nExp == 0) baseDate = expDate-1L;


			IF_FAILED_THROW( GtoDayCountFraction(
				baseDate,
				expDate,
				dayCount,
				&tExp[nExp]));


			// read value
			for (idxM=0; idxM<=nMat-1; idxM++) {
				values[nExp][idxM] = getDouble(ls, "vol");
			}
			if (nExp >= nExpMax)
				throw KFailure("%s: too many maturities.\n",
					routine);
		}

		// Create Matrix 
		mMatrix = DrlTSwaptionMatrix2DNew(
			FALSE,		// TRUE=diagonal, FALSE=vertical 
			(long) freq,
			nExp,
			tExp,
			nMat,
			tMat);
		if (mMatrix == NULL)
			throw KFailure();

		for (idxM=0; idxM<=NMat()-1; idxM++)
		for (idxE=0; idxE<=NExp()-1; idxE++) {
			Val(idxE, idxM) = values[idxE][idxM];
		}

		// We set the frequency to -1 so that it must be set outside
		// the routine !
		Freq() = -1;

		} break;

	default:
		throw KFailure("%s: format not supported.\n",
			routine, mFormat);
	}

	CheckValid();

	return (is);
    }
    catch (...) {
	throw KFailure("%s: failed.\n", routine);
    }
}

// Convenoence routine

static	char*
printInterval(double dblValue)
{
static	char		buf[] = "ERR";
	TDateInterval	interval;
	if (GtoYearsToDateInterval(dblValue, &interval) != SUCCESS) {
		return buf;
	} else {
		return GtoFormatDateInterval(&interval);
	}
}




//---------------------------------------------------------------

ostream&
KEmMatrix::Put(ostream& os, int oneLine) const
{
static	char	routine[] = "ostream& operator<<(ostream&, KEmMatrix&)";
	int	i,j;
	int	idxE, idxM;


	if (isempty()) {
		os << "EMPTY\n";
		return(os);
	}

	switch (mFormat) {
	case KEmMatrix::CUR:
	    // Currency format
	    os << format( "NUMBER_OF_MAT: %d\n", NMat());
	    os << format( "NUMBER_OF_EXP: %d\n", NExp());
	    os << format( "FREQ: %d\n", (int) Freq());
	    os << format( "DIAG: %d\n", Diag());

	    os << format( "        ");
	    for (i=0; i<=NMat()-1; i++) {
		os << format( " %15s", printInterval(TMat(i)));
	    }
	    os << format( "\n");

	    for (j=0; j<=NExp()-1;j++) {
	    	os << format( "%8s", printInterval(TExp(j)));
	        for (i=0; i<=NMat()-1; i++) {
			os << format( " %15s",
				DrlCurPrint(NULL, Val(j,i), 0));
	        }
	        os << format( "\n");
	    }
	    break;

	case KEmMatrix::TXT:
	case KEmMatrix::PRN:
	    /*
	     * Format compatible with import is in spreadsheets
	     */
	    os << format( "\"NUMBER_OF_MAT\"\t%d\n", NMat());
	    os << format( "\"NUMBER_OF_EXP\"\t%d\n", NExp());
	    os << format( "\"FREQ\"\t%d\n", (int) Freq());
	    os << format( "\"DIAG\"\t%d\n", Diag());

	    os << format( "\"VOLS\"");
	    for (i=0; i<=NMat()-1; i++) {
		os << format( "%8.4f", TMat(i));
	    }
	    os << format( "\n");

	    for (j=0; j<=NExp()-1;j++) {
		os << format( "\t%8.4f", TExp(j));
	        for (i=0; i<=NMat()-1; i++) {
		    os << format( "\t%12.8f",
			Val(j,i)) ;
	        }
	        os << format( "\n");
	    }
	    break;
	case KEmMatrix::DRW2SWMAT:
	    /*
	     * Export for Dr Wrapper Type 2
	     */
            os << format( "# Number of expiries\n %d\n", NExp());
            os << format( "# No of forward maturities\n %d \n", NMat());
            os << format( "# Forward maturities (years), "
                "then expiries (months) and volatilities in %%\n");

	    os << format( "        ");
	    for (i=0; i<=NMat()-1; i++) {
		os << format( "\t%4.0f", TMat(i));
	    }
	    os << format( "\n");

	    for (j=0; j<=NExp()-1;j++) {
	    	os << format( "  %4.0f  ", TExp(j)*12e0);

	        for (i=0; i<=NMat()-1; i++) {
			os << format( " %7.3f", Val(j,i) * 1e2);
		}
	        os << format( "\n");
	    }
	    break;

	case KEmMatrix::LIST:
	    /*
	     * List format '<exp> <mat> <value>'
	     */

	    for (idxM=0; idxM<=NMat()-1; idxM++)
	    for (idxE=0; idxE<=NExp()-1; idxE++) {
	    	os << format(" %8.4f", TExp(idxE));
	    	os << format(" %4s", printInterval(TMat(idxM)));
		os << format(" %15s", DrlCurPrint(NULL,
				Val(idxE, idxM), 0));
		os << format("\n");
	    }

	    break;

	case KEmMatrix::LIVERATE:
	    os << format( "\"EXP/UND\" ");
	    for (i=0; i<=NMat()-1; i++) {
		os << format( "\t%4.0f", TMat(i));
	    }
	    os << format( "\n");

	    for (j=0; j<=NExp()-1;j++) {
	    	os << format( "  %8.4f  ", TExp(j));

	        for (i=0; i<=NMat()-1; i++) {
			os << format( "\t%-.4f%%", Val(j,i) * 1e2);
		}
	        os << format( "\n");
	    }
	    break;

	default:
	case KEmMatrix::DEF:
	case KEmMatrix::PERCENT:
	case KEmMatrix::CURK:
	case KEmMatrix::CURM:

	    os << format("NUMBER_OF_MAT: %d\n", NMat());
	    os << format("NUMBER_OF_EXP: %d\n", NExp());
	    os << format("FREQ: %d\n", Freq());
	    os << format("DIAG: %d\n", Diag());

	    os << format( "        ");
	    for (i=0; i<=NMat()-1; i++)
		os << format( "     %4s", printInterval(TMat(i)));
	    os << format( "\n");

	    for (j=0; j<=NExp()-1;j++) {
		os << format( "%8s", printInterval(TExp(j)));

	        for (i=0; i<=NMat()-1; i++) {
		    switch (mFormat) {
		    case KEmMatrix::PERCENT:
			os << format("\t%-.2f",
				Val(j,i) * 1e2);
			break;
		    case KEmMatrix::CUR:
			os << format(" %15s",
				DrlCurPrint(NULL, Val(j,i), 0));
			break;
		    case KEmMatrix::CURK:
			os << format( " %11s",
				DrlCurPrint(NULL, Val(j,i)*1e-3, 0));
			break;
		    case KEmMatrix::CURM:
			os << format( " %7s",
				DrlCurPrint(NULL, Val(j,i)*1e-6, 0));
			break;
		    case KEmMatrix::DEF:
		    default:
			os << format(" %s",
				DrlFloatPrint(NULL, Val(j,i), 8));
			break;
		    }
	        }
	        os << format("\n");
	    }
	    break;
	}

	return(os.flush());
}



//---------------------------------------------------------------

int
CheckSameType(const KEmMatrix& mat1, const KEmMatrix& mat2)
{
    try {
	if (mat1.isempty() && mat2.isempty())
		return (TRUE);

	if ((!mat1.isempty() &&  mat2.isempty()) ||
	    ( mat1.isempty() && !mat2.isempty())) {
		throw KFailure();
	}

	if (!DrlTSwaptionMatrix2DIsSameType(
			mat1.mMatrix,
			mat2.mMatrix)) {
		throw KFailure();
	}
	return(TRUE);
    }
    catch (KFailure) {
	throw KFailure("CheckSameType: inconsistent.\n");
    }
}


//---------------------------------------------------------------

KEmMatrix&
KEmMatrix::operator=(double argument)
{
	char equal[] = "=";

	if (isempty()) return(*this);
	IF_FAILED_THROW( DrlTSwaptionMatrix2DOperScalar(
		mMatrix,
		equal,
		argument));
	return(*this);
}

KEmMatrix&
KEmMatrix::operator+=(double argument)
{
	char plus[] = "+";

	if (isempty()) return(*this);
	IF_FAILED_THROW( DrlTSwaptionMatrix2DOperScalar(
		mMatrix,
		plus,
		argument));
	return(*this);
}

KEmMatrix&
KEmMatrix::operator*=(double argument)
{
	char mult[] = "*";

	if (isempty()) return(*this);
	IF_FAILED_THROW( DrlTSwaptionMatrix2DOperScalar(
		mMatrix,
		mult,
		argument));
	return(*this);
}

KEmMatrix&
KEmMatrix::operator+=(const KEmMatrix& mat)
{
	int	i, j;

	if (isempty()) {
		*this = mat;
	} else {
		CheckSameType(*this, mat);
		for (i=0; i<=mat.NExp()-1; i++)
		for (j=0; j<=mat.NMat()-1; j++)
			Val(i,j) += mat.Val(i,j);
	}
	return(*this);
}

KEmMatrix&
KEmMatrix::operator-=(const KEmMatrix& mat)
{
	int	i, j;

	if (isempty()) {
		*this = mat;
		*this *= (-1e0);
	} else {
		CheckSameType(*this, mat);
		for (i=0; i<=mat.NExp()-1; i++)
		for (j=0; j<=mat.NMat()-1; j++)
			Val(i,j) -= mat.Val(i,j);
	}
	return(*this);
}


//---------------------------------------------------------------

double
KEmMatrix::Sum() const
{
	double	sum = 0e0;
	int	i, j;

	for (i=0; i<=NExp()-1; i++)
	for (j=0; j<=NMat()-1; j++)
		sum += Val(i,j);
	return(sum);
}


//---------------------------------------------------------------


double
KEmMatrix::SumProd(const KEmMatrix& mat) const
{
	int	i, j;
	double	retVal = 0e0;

	CheckSameType(*this, mat);

	for (i=0; i<=mat.NExp()-1; i++)
	for (j=0; j<=mat.NMat()-1; j++)
		retVal += Val(i,j) * mat.Val(i,j);

	return(retVal);
}


//---------------------------------------------------------------


KEmMatrix&
KEmMatrix::Interp(KEmMatrix& mat, TDate baseDate, int interpFlag)
{
static	char	routine[] = "KEmMatrix::Interp";


    try {
	if (interpFlag == FALSE) {
		IF_FAILED_THROW(DrlTSwaptionMatrix2DInterpMatrix(
			mat.mMatrix,
			mMatrix,
			baseDate,
			interpFlag));
	} else {
		IF_FAILED_THROW(DrlTSwaptionMatrix2DInterpMatrix(
			mMatrix,
			mat.mMatrix,
			baseDate,
			interpFlag));
	}
    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
    return(*this);
}






//---------------------------------------------------------------


double
KEmMatrix::Interp(double expYrs, double matYrs)
{
	double	value;

	// Interp volatility volatility
	IF_FAILED_THROW( DrlTSwaptionMatrix2DInterpExpMat(
		(TSwaptionMatrix2D*) *this,
		&value,
		expYrs,
		matYrs,
		FALSE));	 // TRUE=adjoint, FALSE=direct*/

	return(value);
}





//---------------------------------------------------------------


void
KEmMatrix::InterpCurve(
	TDate baseDate,		// (I) base date
	int finalFlag,		// (I) TRUE=final, FALSE=cms 
	int finalDate,		// (I) TRUE=use date, FALSE=use interval
	TDate matDate,		// (I) final maturity date (used if final) 
	KDateInterval matInt,	// (I) fwd maturity interval (used if cms) 

	KDateInterval minMat,	// (I) minimum maturity

	KVector(TDate) &vDates,	// (O) array of vol exp dates
	KVector(double) &vExp,	// (O) array of vol exp time
	KVector(double) &vMat,	// (O) array of vol fwd mat
	KVector(int) &vFreq,	// (O) array of vol freq
	KVector(double) &vRates)	// (B) array of vol 

{
static	char	routine[] = "KEmMatrix::InterpVolCurve";

	int	idxE;
	double	expYrs, matYrs, minMatYrs;
	TDate	expDate;
	double	value;
	int	freq;
	int	spotDateAdded = FALSE;



    try {


	// empty vectors
	//
	//vDates.clear();
	//vExp.clear();
	//vMat.clear();
	//vFreq.clear();
	//vRates.clear();

	vDates.erase(vDates.begin(), vDates.end());
	vExp.erase(vExp.begin(), vExp.end());
	vMat.erase(vMat.begin(), vMat.end());
	vFreq.erase(vFreq.begin(), vFreq.end());
	vRates.erase(vRates.begin(), vRates.end());




	PG_IF_LOGGING(PG_LOG_PRI_D, dppLog << format( \
		"%s: BASEDATE: %10s\n"\
		" DATE        EXP       MAT       FR  INTERP\n", \
		routine, DrlTDatePrint(NULL, baseDate)));


	minMatYrs = minMat.Years();

	for (idxE=0; idxE < NExp(); idxE++) {

		// Check if 1st expiration is less than 1 week
		// if yes, discard it and use spot date instead
		// if no, add spot date as first date 
		//
		expYrs = TExp(idxE);
		if (spotDateAdded == FALSE) {
			if (TExp(0) > 0.010958904e0) {
				idxE--;
			}
			expYrs = 0e0;
			spotDateAdded = TRUE;
		}

		// Now add interp info for expYrs.
		//

		// Get exp date by ACT/365F advance
		expDate = baseDate + (long)(expYrs*365e0);

		/*IF_FAILED_THROW( GtoTDateAdvanceYears(
			baseDate,
			expYrs,
			&expDate));*/

		// Extract interp in years
		if (finalFlag == FALSE) {
			//
			// Vertical calibration
			//
			matYrs = matInt.Years();
			matYrs = MAX(matYrs, minMatYrs);
		} else {
			//
			// Diagonal calibration
			//
			if (finalDate) {
			    IF_FAILED_THROW( GtoDayCountFraction(
				expDate,
				matDate,
				GTO_B30_360,
				&matYrs));
			} else {
				matYrs = matInt.Years() - expYrs;
			}
			matYrs = MAX(matYrs, minMatYrs);
		}

		// Interp volatility volatility
		IF_FAILED_THROW( DrlTSwaptionMatrix2DInterpExpMat(
			(TSwaptionMatrix2D*) *this,
			&value,
			expYrs,
			matYrs,
			FALSE));	 // TRUE=adjoint, FALSE=direct*/


		// Frequency
		freq = Freq();



		// Add to arrays
		vDates.insert(vDates.end(),  expDate);
		vExp.insert(vExp.end(),  expYrs);
		vMat.insert(vMat.end(),  matYrs);
		vFreq.insert(vFreq.end(),  freq);
		vRates.insert(vRates.end(),  value);


		PG_IF_LOGGING(PG_LOG_PRI_D, dppLog << format( \
			" %10s  %8.4f  %8.4f  %2d  %10.6f\n", \
			DrlTDatePrint(NULL, expDate), expYrs, matYrs, \
			freq, value));




	}

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//


void
KEmMatrix::InterpMatCoeffs(
	int ie,			// (I) exp idx
	double tMat,		// (I) mat to interp
	int *imlo,		// (O) lo mat idx
	int *imhi,		// (O) hi mat idx
	double *wmlo,		// (O) lo mat weight
	double *wmhi)		// (O) hi mat weight
{

	IF_FAILED_THROW( DrlGetMatInterpCoeffs(
		this->mMatrix,
		ie,
		tMat,
		imlo,
		imhi,
		wmlo,
		wmhi));
}








//---------------------------------------------------------------
// This routine regresses a matrix move on a vectors
// of matrix components.
// 

void
KEmMatrixRegress(
	int numSwfact,			// (I) number of vol scen 
	KEmMatrix *swfact,		// (I) swaption scenarios 
	KEmMatrix &weiMat,		// (I) weights matr (or NULL) 
	KEmMatrix &beforeMat,		// (I) start matrix (or NULL) 
	KEmMatrix &afterMat,		// (I) end matrix 
	double *amplitudes,		// (O) amplitudes [0..numSwfact-1] 
	double *residual)		// (O) residual 
{
static	char	routine[] = "KEmMatrixRegress";

	int	idxE, idxM, i, j, nExp, nMat;
        int	m;		// (I) number of equations 
        int	n;		// (I) number of unknowns  
        double	**a = NULL;	// (I) input matrix [0..m-1][0..n-1] 
        double	*x = NULL;	// (O) solution [0..n-1] 
        double	*b = NULL;	// (I) inhomogeneous term [0..m-1] 
        double	*r = NULL;	// (I) residual [0..m-1] 
        int	regType = 0;	// (I) see DrlMatrixSvdRegularizeEV 
        double	alpha = 1e-9;	// (I) see DrlMatrixSvdRegularizeEV 
        double	param = 0e0;	// (I) see DrlMatrixSvdRegularizeEV 
	double	res0, rfct;

    try {
	// 

	nExp = afterMat.NExp();
	nMat = afterMat.NMat();

	n = numSwfact;
	m = nExp*nMat;

	// 
	ASSERT_OR_THROW((a = DrlDoubleMatrAlloc(0, m-1, 0, n-1)) != NULL);
	ASSERT_OR_THROW((x = DrlDoubleVectAlloc(0, n-1)) != NULL);
	ASSERT_OR_THROW((b = DrlDoubleVectAlloc(0, m-1)) != NULL);
	ASSERT_OR_THROW((r = DrlDoubleVectAlloc(0, m-1)) != NULL);

	// Check consistency
	//
	for (j=0; j<=n-1; j++) 
	    CheckSameType(beforeMat, swfact[j]);

	// Copy matrix as vector for SVD call
	//
	for (idxE=0; idxE<=nExp-1; idxE++)
	for (idxM=0; idxM<=nMat-1; idxM++) {
	    i = idxM + idxE * nMat;

	    if (!weiMat.IsEmpty()) {
		for (j=0; j<=n-1; j++) {
			a[i][j] = swfact[j].Val(idxE, idxM) *
				sqrt(weiMat.Val(idxE, idxM));
		}
		b[i] = afterMat.Val(idxE, idxM) *
				sqrt(weiMat.Val(idxE, idxM));
		if (!beforeMat.IsEmpty()) {
			b[i] = beforeMat.Val(idxE, idxM) *
				sqrt(weiMat.Val(idxE, idxM));
		}
	    } else {
		for (j=0; j<=n-1; j++) {
			a[i][j] = swfact[j].Val(idxE, idxM);
		}
		b[i] = afterMat.Val(idxE, idxM);
		if (!beforeMat.IsEmpty()) {
			b[i] -= beforeMat.Val(idxE, idxM);
		}
	    }
	}


	// Perform SVD
	//
	IF_FAILED_THROW (DrlRealLinSysSvd(
		m,
		n,
		a,
		x,
		b,
		regType,
		alpha,
		param));

	for (j=0; j<=n-1; j++) {
		amplitudes[j] = x[j];
	}


	// Compute residual vect.
	//
	for (idxE=0; idxE<=nExp-1; idxE++)
	for (idxM=0; idxM<=nMat-1; idxM++) {
	    i = idxM + idxE * nMat;
	    r[i] = 0e0;
	    for (j=0; j<=n-1; j++) {
		r[i] += a[i][j]*x[j];
	    }
	    r[i] -= b[i];
	}


	{
		double	mv1, mv2, mvmin, mvmax;
		double	rv1, rv2, rvmin, rvmax;

		dppLog << format("%s:", routine);

		IF_FAILED_THROW( DrlVectNorm(b, m, "L1", &mv1)); mv1 /= m;
		IF_FAILED_THROW( DrlVectNorm(b, m, "L2", &mv2)); mv2 /= sqrt((double)m);
		IF_FAILED_THROW( DrlVectNorm(b, m, "LMIN", &mvmin));
		IF_FAILED_THROW( DrlVectNorm(b, m, "LMAX", &mvmax));

		IF_FAILED_THROW( DrlVectNorm(r, m, "L1", &rv1)); rv1 /= m;
		IF_FAILED_THROW( DrlVectNorm(r, m, "L2", &rv2)); rv2 /= sqrt((double)m);
		IF_FAILED_THROW( DrlVectNorm(r, m, "LMIN", &rvmin));
		IF_FAILED_THROW( DrlVectNorm(r, m, "LMAX", &rvmax));

		dppLog << format(" L1=%5.2f/%5.2f %% ", mv1*1e2, rv1*1e2);
		dppLog << format(" L2=%5.2f/%5.2f %% ", mv2*1e2, rv2*1e2);
		dppLog << format(" MIN=%5.2f/%5.2f %% ", mvmin*1e2, rvmin*1e2);
		dppLog << format(" MAX=%5.2f/%5.2f %% ", mvmax*1e2, rvmax*1e2);


		dppLog << endl;
	}



	// Compute residual : absolute everage deviation 
	//
	*residual = 0e0;
	res0 = 0e0;
	for (idxE=0; idxE<=nExp-1; idxE++)
	for (idxM=0; idxM<=nMat-1; idxM++) {
	    i = idxM + idxE * nMat;
	    rfct = 0e0;
	    for (j=0; j<=n-1; j++) {
		rfct += a[i][j]*x[j];
	    }
	    *residual += fabs(b[i]-rfct);
	    res0 += fabs(b[i]);
	}
	*residual /= nExp*nMat;
	res0 /= nExp*nMat;

	/*dppLog << format("%s: RES0=%8.4f %%   RES1=%8.4f %%\n",
			routine, res0*1e2, *residual*1e2);*/



#ifdef	__DEBUG__

	DrlFPrintf(stdout, "MOV ");
	for (idxM=0; idxM<=nMat-1; idxM++)
		DrlFPrintf(stdout, "   %3s", DrlDoubleIntervalPrint(NULL, 
			TMat(beforeMat, idxM)));
	DrlFPrintf(stdout, "\n");
	for (idxE=0; idxE<=nExp-1; idxE++) {
		DrlFPrintf(stdout, "%3s ", DrlDoubleIntervalPrint(NULL, 
			TExp(beforeMat, idxE)));
		for (idxM=0; idxM<=nMat-1; idxM++) {
	    		i = idxM + idxE * nMat;
			DrlFPrintf(stdout, "%5.2f ", b[i]*1e2);
		}
		DrlFPrintf(stdout, "\n");
	}

	DrlFPrintf(stdout, "RESD");
	for (idxM=0; idxM<=nMat-1; idxM++)
		DrlFPrintf(stdout, "   %3s", DrlDoubleIntervalPrint(NULL, 
			TMat(beforeMat, idxM)));
	DrlFPrintf(stdout, "\n");
	for (idxE=0; idxE<=nExp-1; idxE++) {
		DrlFPrintf(stdout, "%3s ", DrlDoubleIntervalPrint(NULL, 
			TExp(beforeMat, idxE)));
		for (idxM=0; idxM<=nMat-1; idxM++) {
	    		i = idxM + idxE * nMat;
	    		for (rfct=0e0, j=0; j<=n-1; j++) rfct += a[i][j]*x[j];
			DrlFPrintf(stdout, "%5.2f ", rfct*1e2);
		}
		DrlFPrintf(stdout, "\n");
	}



#endif

	



	DrlDoubleMatrFree(a, 0, n-1, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);
	DrlDoubleVectFree(b, 0, m-1);
	DrlDoubleVectFree(r, 0, m-1);
	return;

    }
    catch (KFailure) {
	DrlDoubleMatrFree(a, 0, n-1, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);
	DrlDoubleVectFree(b, 0, m-1);
	DrlDoubleVectFree(r, 0, m-1);
	throw KFailure("%s: failed.\n");
    }
}















