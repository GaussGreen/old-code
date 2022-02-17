/************************************************************************
 * Module:	DRL
 * Submodule:	FPARS - Formula Parsing
 * File:	
 * Function:	Formula Parsing
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "kvtree.h"
#include "kvtspar.h"		// prototype consistency
#include "vtsparh.h"		// Symbol table definition
#include "ktypes.h"         // SmoothType, Ko/Ki type

#include "kutilios.h"

extern "C" {
#include "cerror.h"		// DppErrMsg
};



extern	"C" {
					// yysetinput: store input string
extern	int	vtsparsetinput(const char *s);
					// yyparse:
extern	int	vtsparparse();
};



static	KTSlice			*vtsParRetVal;
static	KVTree			*vtsParVt;
static	int			vtsParNXArgs;
static	KTSlice			*vtsParXArgs[32];
static	KMap(String,double)	*vtsConstTable;	// constants table
static	String			vtsCurveName;


/*----------------------------------------------------------------------
 * Virtual Tree dynamic time slice formula parser tool.
 * \begin{arglist}
 * \item[vt] the virtual tree,
 * \item[tpIdx] timepoint index at which evaluated,
 * \item[nx] the number of X-variables in the formula.
 * \item[x] array of length "nx" of the timeslices to be substituted
 *          to x-variables.
 * \item[ny] the number of Y-variables in the formula.
 * \item[y] array of length "ny" of the timeslices to be substituted.
 *          to y-variables.
 * \item[formula] a formula of variables labelled
 *   "x0", "x1",\dots,"xNX-1","y0",\dots,"yNY-2" that contains
 *    \begin{itemize}
 *    \item  numerical constants, such as {\tt 1.0} or {\tt 4.56e-2},
 *    \item  basic arithmetic operations $+$, $-$, $*$ and $/$
 *    \item  {\tt MIN(E1,E2)} and {\tt MAX(E1,E2)}
 *    \item  the {\tt IFPOS(E0,E1,E2)} statement which returns
 *           the expression E1 (resp. E2) if expression E0 is
 *           positive (resp.negative).
 *    \end{itemize}
 * \item[retVal] on successful exit, contains the timeslice obtained
 *    by evaluating the formula by replacing the formal arguments
 *    "x0", "x1",\dots,"xNX-1", "y0", "y1",\dots,"yNY-1",
 *     by the values given in the input array
 *   "x".
 * \end{arglist}
 * {\bf WARNING: since the routine uses static variables, the code
 * is NOT reentrant.}
 */

void
KTSliceParEval(
	KVTree &vt,			// (I) virtual tree 
	int nx,				// (I) # of X args
	KTSlice *x,			// (I) arrays of X args timeslices
	const char *formula,		// (I) formula
	KMap(String,double)& constTable,	// (I) constants table
	KTSlice *retVal)		// (O) return value timeslice
{
static	char	routine[] = "KTSliceParEval";

	int	i;

    try {

	// copy formula in buffer 
	//
	if (vtsparsetinput(formula) != SUCCESS) {
		throw KFailure("%s: yysetinput failed.\n", routine);
	}

	// copy arguments buffer
	//
	for(i=0; i<=nx-1;i++)
		vtsParXArgs[i] = &x[i];
	vtsParNXArgs = nx;

	vtsParVt = &vt;
	vtsParRetVal = retVal;
	vtsConstTable = &constTable;

	// We use the curve name of the first slice
	// to initialize temporay slices needed by the parser.
	vtsCurveName = x[0].GetCurveName();

	// call the parser: yyparse returns 0 if success,
	// a non zero code otherwise
	//
	if (vtsparparse() != 0) 
		throw KFailure("%s: yyparse failed.\n", routine);

	// Clean the Flex buffer
	vtsparLexClean();

    }
    catch (KFailure) {
	    throw KFailure("%s: Parsing `%s' failed.\n", routine, formula);
    }

}



//----------------------------------------------------------------------
//
//

void
KTSliceParEvalV(
	KVTree &vt,			// (I) virtual tree 
	const char *formula,		// (I) formula
	KMap(String,double)& constTable,// (I) constants table
	KTSlice *retVal,		// (O) return value timeslice
	int nx,				// (I) # of X args
	// KTSlice *x1			// (I) X1  arg timeslice
	// ..
	// KTSlice *xNx			// (I) Xnx arg timeslice
	...)
{
static	char	routine[] = "KTSliceParEvalV";

	int	i;
	va_list	ap;

    try {

	// copy formula in buffer 
	//
	if (vtsparsetinput(formula) != SUCCESS) {
		throw KFailure("%s: yysetinput failed.\n", routine);
	}

	// pointers to arguments buffer
	//
	va_start(ap, nx);

	for(i=0; i<=nx-1;i++) {
		vtsParXArgs[i] = (KTSlice*) va_arg(ap, KTSlice*);
	}
	va_end(ap);
	vtsParNXArgs = nx;


	vtsParVt = &vt;
	vtsParRetVal = retVal;
	vtsConstTable = &constTable;

	// call the parser: yyparse returns 0 if success,
	// a non zero code otherwise
	//
	if (vtsparparse() != 0) 
		throw KFailure("%s: yyparse failed.\n", routine);

    }
    catch (KFailure) {
	    throw KFailure("%s: Parsing `%s' failed.\n", routine, formula);
    }

}


//----------------------------------------------------------------------

static	int
vtsGetConstTable(const char *name, double *value)
{
	String	sName(name);


/*
	dppLog << "Defined variables: " << endl;
	for(KMap(String,double)::iterator q=vtsConstTable->begin();
	    q != vtsConstTable->end();
	    ++q) {
		dppLog << format("%10s  %18.10f",
			(const char*) ((*q).first).c_str(),
			(*vtsConstTable)[(*q).first]) << endl;
	}
	pLog << "Defined variables: DONE:" << endl;
*/


	KMap(String,double)::iterator p = vtsConstTable->find(sName);

	if (p == vtsConstTable->end()) {
		dppErr << "Undefined variable `" << sName << "'" << endl;
		dppErr << "Defined variables: " << endl;
		for(KMap(String,double)::iterator q=vtsConstTable->begin();
		    q != vtsConstTable->end();
		    ++q) {
			dppErr << format("%10s  %18.10f",
				(const char*) ((*q).first).c_str(),
				(*vtsConstTable)[(*q).first]) << endl;
		}


		return(FAILURE);
	}
	*value = (*vtsConstTable)[(*p).first];
	return(SUCCESS);

}





/************************************************************************
 *
 *	A C T I O N S   S E C T I O N
 *
 ************************************************************************/


extern "C" {

//----------------------------------------------------------------------
// The following routines are the action routines for the yacc parser.
// WARNING: since they use static variables, the code
// is NOT reentrant.
//


//----------------------------------------------------------------------
// Copy the slice value to the static variable return value.

int
vtsParCopyReturn(VTSSym *symb)
{
	KTSlice	*ts = (KTSlice*) symb->ptr;

	*vtsParRetVal = *ts;
	delete ts;

	return(SUCCESS);
}


//----------------------------------------------------------------------
// Copy double value to the static variable return value.

int
vtsParCopyReturnConstant(double value)
{
	*vtsParRetVal = value;
	return(SUCCESS);
}



//----------------------------------------------------------------------
// Retreives the argument # idx and puts it in the symbol.

int
vtsParGetArg(VTSSym *symbO, int idx)
{
static	char	routine[] = "vtsParGetArg";

	if ((idx < 0) ||  (idx >= vtsParNXArgs)) {
		DppErrMsg("%s: bad X argument index %d (max %d).\n",
		    routine, idx, vtsParNXArgs);
		return (FAILURE);
	}

	KTSlice	*ts = new KTSlice(
		*vtsParVt,
		"arg",
		vtsParXArgs[idx]->GetCurveName());
	*ts = *(vtsParXArgs[idx]);
	symbO->ptr = ts;

	return(SUCCESS);
}


//----------------------------------------------------------------------
// Copies a double in a time slice

int
vtsParLetDouble(VTSSym *symbO, double value)
{
static	char	routine[] = "vtsParLetDouble";

	KTSlice	*ts = new KTSlice(*vtsParVt, "arg_dbl", vtsCurveName);
	*ts = value;
	symbO->ptr = ts;

	return(SUCCESS);
}

//----------------------------------------------------------------------
// Retreives a double vrariable form the table and
// copies a double variable in a time slice

int
vtsParLetDoubleVariable(VTSSym *symbO, const char* dvarname)
{
static	char	routine[] = "vtsParLetDouble";
	double	value;

	//
	if (vtsGetConstTable(dvarname, &value) != SUCCESS)
		return(FAILURE);
	

	//
	KTSlice	*ts = new KTSlice(*vtsParVt, "arg_dbl", vtsCurveName);
	*ts = value;
	symbO->ptr = ts;

	return(SUCCESS);
}


//----------------------------------------------------------------------
//
// Basic algebracic operations +, -, *, /, log, exp, -(minus)

int
vtsParUnaryOper(VTSSym *symbO, VTSSym *symb1, char *oper)
{
	KTSlice	*ts0 = (KTSlice*) symbO->ptr;
	KTSlice	*ts1 = (KTSlice*) symb1->ptr;

	switch (oper[0]) {
	case '-':
		*ts1 *= (-1);
		break;
	case '=':
		symbO->ptr = symb1->ptr;
		symb1->ptr = NULL;
		break;
	default:
		DppErrMsg("vtsParUnaryOper: invalid operation %s\n", oper);
		return(FAILURE);
	}

	symbO->ptr = ts1;
	return(SUCCESS);
}


int
vtsParScalarOper(VTSSym *symbO, VTSSym *symb1, double value, char *oper)
{
	KTSlice	*ts1 = (KTSlice*) symb1->ptr;

	switch (oper[0]) {
	case '=':
		*ts1 = value;
		break;
	case '+':
		*ts1 += value;
		break;
	case '-':
		*ts1 -= value;
		break;
	case '*':
		*ts1 *= value;
		break;
	case '/':
		*ts1 /= value;
		break;
	case 'M':
		ts1->max(value);
		break;
	case 'm':
		ts1->min(value);
		break;
	case '^':
	default:
		DppErrMsg("vtsParScalarOper: invalid operation `%s'\n", oper);
		return(FAILURE);
	}


	symbO->ptr = ts1;
	return(SUCCESS);
}


int
vtsParBinaryOper(VTSSym *symbO, VTSSym *symb1, VTSSym *symb2, char *oper)
{
	KTSlice	*ts1 = (KTSlice*) symb1->ptr;
	KTSlice	*ts2 = (KTSlice*) symb2->ptr;

	switch (oper[0]) {
	case '=':
		*ts1 = *ts2;
		break;
	case '+':
		*ts1 += *ts2;
		break;
	case '-':
		*ts1 -= *ts2;
		break;
	case '*':
		*ts1 *= *ts2;
		break;
	case '/':
		*ts1 /= *ts2;
		break;
	case 'M':
		ts1->max(*ts2);
		break;
	case 'm':
		ts1->min(*ts2);
		break;
	default:
		DppErrMsg("vtsParBinaryOper: invalid operation %s\n", oper);
		return(FAILURE);
	}


	delete ts2;
	symbO->ptr = ts1;
	return(SUCCESS);
}



int
vtsParTernaryOper(VTSSym *symbO, VTSSym *symb1, VTSSym *symb2, VTSSym *symb3, char *oper)
{
	KTSlice	*ts1 = (KTSlice*) symb1->ptr;
	KTSlice	*ts2 = (KTSlice*) symb2->ptr;
	KTSlice	*ts3 = (KTSlice*) symb3->ptr;
	KTSlice	*ts = new KTSlice(*ts1);

	switch (oper[0]) {
	case '?':
		// Returns ts2 if ts1 is >=0, ts3 otherwise
		vtsParVt->TSliceSpecialOper(*ts, "IFPOS",
			(void*) ts1,
			(void*) ts2,
			(void*) ts3);

		break;
	default:
		DppErrMsg("vtsParTernaryOper: invalid operation %s\n", oper);
		return(FAILURE);
	}


	delete ts1;
	delete ts2;
	delete ts3;
	symbO->ptr = ts;
	return(SUCCESS);
}

int 
vtsParStepOper(VTSSym *symbO, VTSSym *symb1, VTSSym *lb, VTSSym *ub, 
               char *ioWindow, char *kiSmooth)
{
    KTSlice *ts1  = (KTSlice*) symb1->ptr;   
    KTSlice *tsLB = (KTSlice*) lb->ptr;
    KTSlice *tsUB = (KTSlice*) ub->ptr;
    KTSlice *ts   = new KTSlice(*ts1);

    KKnockIO    mIOWindow;
    KSmooth     mSmooth;

//    DppErrMsg("iowindow = %s\t kismooth = %s\n", ioWindow, kiSmooth);

    switch (kiSmooth[0])    
    {
    case 'N':
//        cout << "nosmooth" << endl;
        mSmooth = NO_SMOOTH;     
        break;

    case 'S':
        mSmooth = SINGLE_SMOOTH;     
//        cout << "single smooth" << endl;
        break;

    case 'D':
        mSmooth = DOUBLE_SMOOTH;     
//        cout << "double smooth" << endl;
        break;

    default:
		DppErrMsg("vtsParStepOper: invalid operation %s\n", kiSmooth);
		return(FAILURE);
    }

    switch (ioWindow[0])
    {
    case 'I':
//        cout << "inside" << endl;
        mIOWindow = CRX_KNOCK_IN;
        break;
    case 'O':
//        cout << "outside" << endl;
        mIOWindow = CRX_KNOCK_OUT;
        break;
    default:
        DppErrMsg("vtsParStepOper: invalid operation %s\n", ioWindow);
        return(FAILURE);
    }

    vtsParVt->TSliceSpecialOper(*ts, "STEP",
            (void*) ts1,
            (void*) tsLB,
            (void*) tsUB,
            (void*) &mIOWindow,
            (void*) &mSmooth);

    delete ts1;
    delete tsLB;
    delete tsUB;
    symbO->ptr = ts;
    return(SUCCESS);
}


};	// extern "C"

