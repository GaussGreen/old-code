/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	C. Daher, D. Liu
 * Revision:	$Header$
 ***************************************************************/
#include "kutilios.h"    /* Standard definitions & error hadling */

#define	_kmrntree_SRC
#include "kmrntree.h"

#ifndef BARRIER_TOL
#define BARRIER_TOL (1e-7)
#endif


static	double DrSmoothStep(
	double  UpValue,  /* (I) used if Index > Barrier + step  */
	double  DownValue,/* (I) used if Index < Barrier - step  */
	double  Index,    /* (I) index value                     */
	double  Barrier,  /* (I) barrier value                   */
	double  Step);    /* (I) smooth between Barrier +/- step */

static	double  DrSmoothMax(
	double  a,     /* (I) Argument 1                 */
	double  b,     /* (I) Argument 2                 */
	double  Step); /* (I) smooth in [-step, step]    */



//--------------------------------------------------------------
// Calculate array offsets. There are three components:
//
//  1)  node indices run into [-halfwidth, halfwidth] in each dimension 
//          whereas memory is allocated as [0, 2 * halfwidth]. This gives an 
//          offset of halfwidth.
//  2)  memory is allocated linearly in dimension 2 and 3. This gives an
//          offset of i * width[1] in dimension 2 and i * width[1] * width[2]
//          + j * width[1] in dimension 2.
//  3)  the axis of the ellipse has a slope. The coordinate of the axis
//          is given by (bottom2[i]+top2[i])/2 at index [i] in dimension 2 and
//          (bottom3[i][j]+top3[i][j])/2 at index [i][j] in dimension 3.
//
//

int
KMrNTree::NodeOffset(
	int Dim,			// (I) Dimension
	int j0,				// (I) Node indices
	int j1,				// (I) 
	int t)				// (I) Current time point
{
	int	Offset;       /* Returned offset */

	switch (Dim) {
	case 3:
	    	//NodeCheck(2, j0, j1, 0, t, 1);

	    	Offset = (j0 + this->mHalfWidth[0]) 
		     	* mWidth[1] * mWidth[2]
		     	+(j1-((mOutTop2[t][j0]+mOutBottom2[t][j0])>>1)
		     	+ mHalfWidth[1]) * mWidth[2]
		     	- ((mOutTop3[t][j0][j1]+mOutBottom3[t][j0][j1])>>1)
			+ mHalfWidth[2];
		break;
	case 2:
		//NodeCheck(1, j0, 0, 0, t, 1);

		Offset = (j0 + mHalfWidth[0]) * mWidth[1]
			-((mOutTop2[t][j0]+mOutBottom2[t][j0])>>1)
			+ mHalfWidth[1];
		break;
	case 1:

		Offset = mHalfWidth[0];
		break;
	}

	return (Offset);
}

//--------------------------------------------------------------
//

void
KMrNTree::NodeCheck(
	int Dim,			// (I) Dimension
	int j0,				// (I) Node indices
	int j1,				// (I) 
	int j2,				// (I) 
	int t,				// (I) Current time point
	int inOut)			// (I) 0=in, 1=out
{
static	char	routine[] = "KMrNTree::NodeCheck";

	if (Dim >= 1) {
	    if ((j0 < mOutBottom1[t]) || (j0 > mOutTop1[t]))
		throw KFailure("%s: invalid j0 (%d) at "
			"t=%d, OutBottom=%d, OutTop=%d.\n",
			routine, j0,
			t, mOutBottom1[t], mOutTop1[t]);
	 if (Dim >= 2) {
	    if ((j1 < mOutBottom2[t][j0]) || (j1 > mOutTop2[t][j0]))
		throw KFailure("%s: invalid j1 (%d) at "
			"t=%d, j0=%d, OutBottom=%d, OutTop=%d.\n",
			routine, j1,
			t, j0, mOutBottom2[t][j0], mOutTop2[t][j0]);
	  if (Dim >= 3) {
	    if ((j2 < mOutBottom3[t][j0][j1]) || (j2 > mOutTop3[t][j0][j1]))
		throw KFailure("%s: invalid j2 (%d) at "
			"t=%d, j0=%d, j1=%d, OutBottom=%d, OutTop=%d.\n",
			routine, j2,
			t, j0, j1, mOutBottom3[t][j0][j1], mOutTop3[t][j0][j1]);
	  }
	 }
	}
}


//--------------------------------------------------------------
//

double*
KMrNTree::sliceNew(int sliceDim)
{
static	char		routine[] = "KMrNTree::sliceNew";
	register int	idx;
	int		size;
	double		*ptr;

	// Allocate slice
	if (sliceDim > mNbFactor) {
		throw KFailure("%s: invalid dimension %d (nb fact is %d).\n",
			routine, sliceDim, mNbFactor);
	}

	switch (sliceDim) {
	case 1:
		size = this->mWidth[0];
		break;
	case 2:
		size = this->mWidth[0]
		     * this->mWidth[1];
		break;
	case 3:
		size = this->mWidth[0]
		     * this->mWidth[1]
		     * this->mWidth[2];
		break;
	default:
		throw KFailure("%s: invalid dimension %d.\n",
			routine, sliceDim);
	}

	ptr = new double [size];
	if (ptr == NULL)
		throw KFailure("%s: memory allocation failure.\n", routine);

	// Clear to 0
	//
	for (idx=0; idx<size; idx++) ptr[idx] = 0e0;


	return(ptr);
}


//--------------------------------------------------------------
//

void
KMrNTree::sliceDelete(double *ptr)
{
static	char	routine[] = "KMrNTree::sliceDelete";

	if (ptr == NULL) return;
	delete [] ptr;
	return;
}

//--------------------------------------------------------------
//

void
KMrNTree::slicePrint(
	double *ptr,			// (I) data
	int sliceDim,			// (I) slice dimension
	int tpIdx,			// (I) current time point index
	int outLim,			// (I) TRUE=out limits.
	ostream &os)			// (I) printed to
{
	int	j0, j1, j2;
	int	offset;
	double	*sliceL;

	if (outLim != FALSE) {
	    os << format("TP=%3d OUTER LIMITS.\n", tpIdx);

	    switch (sliceDim) {
            case 1: 
		offset = this->NodeOffset(1, 0, 0, tpIdx);
		sliceL = (double*) ptr + offset;

		for (j0 = mOutBottom1[tpIdx]; j0 <= mOutTop1[tpIdx]; j0++) {
			os << format("(%4d)          ", j0);
			os << format("    %12.8f\n", sliceL[j0]);
		}
		break;

	    case 2:
		for (j0 = mOutBottom1[tpIdx]; j0 <= mOutTop1[tpIdx]; j0++) {
			offset = this->NodeOffset(2, j0, 0, tpIdx);
			sliceL = (double*) ptr + offset;
            	for (j1 = mOutBottom2[tpIdx][j0];
					j1 <= mOutTop2[tpIdx][j0]; j1++) {
			os << format("(%4d,%4d)     ", j0, j1);
			os << format("    %12.8f\n", sliceL[j1]);
		}
		}
		break;

	    case 3:
		for (j0 = mOutBottom1[tpIdx]; j0 <= mOutTop1[tpIdx]; j0++)
            	for (j1 = mOutBottom2[tpIdx][j0]; j1 <= mOutTop2[tpIdx][j0]; 
		     j1++)
		{
               		offset = this->NodeOffset(3, j0, j1, tpIdx);
			sliceL = (double*) ptr + offset;

	    	for (j2 = mOutBottom3[tpIdx][j0][j1];
				j2 <= mOutTop3[tpIdx][j0][j1]; j2++) {
			os << format("(%4d,%4d,%4d)", j0, j1, j2);
			os << format("    %12.8f\n", sliceL[j2]);
		}
		}
		break;

		}
	} else {
	    os << format("TP=%3d INNER LIMITS.\n", tpIdx);

	    switch (sliceDim) {
            case 1: 
		offset = this->NodeOffset(1, 0, 0, tpIdx);
		sliceL = (double*) ptr + offset;

		for (j0 = mBottom1[tpIdx]; j0 <= mTop1[tpIdx]; j0++) {
			os << format("(%4d)          ", j0);
			os << format("    %12.8f\n", sliceL[j0]);
		}
		break;

	    case 2:
		for (j0 = mBottom1[tpIdx]; j0 <= mTop1[tpIdx]; j0++) {
			offset = this->NodeOffset(2, j0, 0, tpIdx);
			sliceL = (double*) ptr + offset;
            	for (j1 = mBottom2[tpIdx][j0]; j1 <= mTop2[tpIdx][j0]; j1++) {
			os << format("(%4d,%4d)     ", j0, j1);
			os << format("    %12.8f\n", sliceL[j1]);
		}
		}
		break;

	    case 3:
		for (j0 = mBottom1[tpIdx];     j0 <= mTop1[tpIdx];     j0++)
            	for (j1 = mBottom2[tpIdx][j0]; j1 <= mTop2[tpIdx][j0]; j1++) {
               		offset = this->NodeOffset(3, j0, j1, tpIdx);
			sliceL = (double*) ptr + offset;

	    	for (j2 = mBottom3[tpIdx][j0][j1];
					j2 <= mTop3[tpIdx][j0][j1]; j2++) {
			os << format("(%4d,%4d,%4d)", j0, j1, j2);
			os << format("    %12.8f\n", sliceL[j2]);
		}
		}
		break;

	    }
	} 
}


//--------------------------------------------------------------
// Performs a scalar operation on a timeslice.
//


void
KMrNTree::sliceScalarOper(
	double *data,			// (B) slice data 
	const int sliceDim,		// (I) slice dimension
	double argument,		// (I) argument of operation
	KOper oper)			// (I) operation type
{
static	char	routine[] = "KMrNTree::sliceScalarOper";

register int	i, j, k;		// Node indices
register int	j0, j1, j2;		// Node indices
	int	nFact = mNbFactor;	// used by some macros
	int     Top1, Bottom1;		// Tree limits (1rst dim)
 	int     *Top2, *Bottom2;	// Tree limits (2nd dim) 
	int     **Top3, **Bottom3;	// Tree limits (3rd dim)
	int	offset;			// Node offset


	int	t = TPIdxCurrent();	// Current time point
	double	*Slice = data;		// data slice
	double  *SliceL;

	double	maxDiff;
	double	*auxSlice = this->NewPrice;	// Auxiliary slice
	double	*auxSliceL;

	Top1    = this->mTop1[t];
	Top2    = this->mTop2[t];
	Top3    = this->mTop3[t];
	Bottom1 = this->mBottom1[t];
	Bottom2 = this->mBottom2[t];
	Bottom3 = this->mBottom3[t];


	// The loop is the same for all operations, so make it a macro
	// It is a little ugly, but efficient 

#define	DO_LOOP \
	switch (sliceDim) { \
        case 1: \
		offset = NodeOffset(1, 0, 0, t); \
		SliceL = Slice + offset; \
		for (i = Bottom1; i <= Top1; i ++) { \
			OPERATION(SliceL[i]); \
		} \
		break; \
        case 2: \
		for (i = Bottom1; i <= Top1; i ++) {\
                	offset = NodeOffset(2, i, 0, t); \
			SliceL = Slice + offset; \
			for (j = Bottom2[i]; j <= Top2[i]; j++) { \
				OPERATION(SliceL[j]); \
			} \
		} \
		break; \
        case 3: \
		for (i = Bottom1; i <= Top1; i ++) \
                for (j = Bottom2[i]; j <= Top2[i]; j++) { \
		    offset = NodeOffset(3, i, j, t); \
		    SliceL = Slice + offset; \
		    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++) { \
				OPERATION(SliceL[k]); \
		    } \
		} \
		break; \
	default: \
		KFailure("%s: invalid slice dimension (%d).\n", routine, \
			sliceDim); \
	} 


	switch (oper) {

	//----------------------------------------------
	// Standard arithmetic operations
	//----------------------------------------------

	case COPY:
		// Assignement
#undef	OPERATION
#define	OPERATION(x)	x = argument
		DO_LOOP
#undef	OPERATION
		break;
	case ADD:
		// Addition
#undef	OPERATION
#define	OPERATION(x)	x += argument
		DO_LOOP
#undef	OPERATION
		break;
	case SUB:
		// Subtraction
#undef	OPERATION
#define	OPERATION(x)	x -= argument
		DO_LOOP
#undef	OPERATION
		break;
	case MULT:
		// Multiplication
#undef	OPERATION
#define	OPERATION(x)	x *= argument
		DO_LOOP
#undef	OPERATION
		break;
	case DIV:
		// Division
#undef	OPERATION
#define	OPERATION(x)	x /= argument
		DO_LOOP
#undef	OPERATION
		break;
	case POW:
		// Power
#undef	OPERATION
#define	OPERATION(x)	x = pow(x, argument)
		DO_LOOP
#undef	OPERATION
		break;
	case LOG:
		// Logarithm
#undef	OPERATION
#define	OPERATION(x)	x = log(x)
		DO_LOOP
#undef	OPERATION
		break;

	//----------------------------------------------
	// Smooth operation: MIN
	//----------------------------------------------

	case MIN:

	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x)	x = (x <= argument ? x : argument)
			DO_LOOP
#undef	OPERATION
	} else {
	    // Copy into auxiliary slice
	    //
	    sliceUnaryOper(auxSlice, sliceDim, data, COPY);

	    // Perform smooth max
	    //
	    switch (sliceDim) {
            case 1: 
		offset    = this->NodeOffset(1, 0, 0, t);
		SliceL    = (double*) data + offset;
		auxSliceL = (double*) auxSlice + offset;

		for (j0 = Bottom1; j0 <= Top1; j0++) {
			maxDiff = GetIndexStep(data, 1, j0, 0, 0, t);
			SliceL[j0] += argument
				- DrSmoothMax(auxSliceL[j0], argument,
					maxDiff);
		}
		break;

	    case 2:
		for (j0 = Bottom1; j0 <= Top1; j0++) {
			offset    = this->NodeOffset(2, j0, 0, t);
			SliceL    = (double*) data + offset;
			auxSliceL = (double*) auxSlice + offset;
            	for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
			maxDiff = GetIndexStep(data, 2, j0, j1, 0, t);
			SliceL[j1] += argument
				- DrSmoothMax(auxSliceL[j1], argument,
					maxDiff);
		}
		}
		break;

	    case 3:
		for (j0 = Bottom1; j0 <= Top1; j0++)
            	for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
               		offset    = this->NodeOffset(3, j0, j1, t);
			SliceL    = (double*) data + offset;
			auxSliceL = (double*) auxSlice + offset;

	    	for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
			maxDiff = GetIndexStep(data, 3, j0, j1, j2, t);
			SliceL[j2] = argument
				- DrSmoothMax(auxSliceL[j2], argument,
					maxDiff);
		}
		}
		break;

	    }
	}

	break;

	//----------------------------------------------
	// Smooth operation: MAX
	//----------------------------------------------

	case MAX:
	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x)	x = (x >= argument ? x : argument)
		DO_LOOP
#undef	OPERATION
	} else {
	    // Copy into auxiliary slice
	    //
	    sliceUnaryOper(auxSlice, sliceDim, data, COPY);

	    // Perform smooth max
	    //
	    switch (sliceDim) {
            case 1: 
		offset    = this->NodeOffset(1, 0, 0, t);
		SliceL    = (double*) data + offset;
		auxSliceL = (double*) auxSlice + offset;

		for (j0 = Bottom1; j0 <= Top1; j0++) {
			maxDiff = GetIndexStep(auxSlice, 1, j0, 0, 0, t);
			SliceL[j0] = DrSmoothMax(auxSliceL[j0], argument,
					maxDiff);
		}
		break;

	    case 2:
		for (j0 = Bottom1; j0 <= Top1; j0++) {
			offset    = this->NodeOffset(2, j0, 0, t);
			SliceL    = (double*) data + offset;
			auxSliceL = (double*) auxSlice + offset;
            	for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
			maxDiff = GetIndexStep(auxSlice, 2, j0, j1, 0, t);
			SliceL[j1] = DrSmoothMax(auxSliceL[j1], argument,
					maxDiff);
		}
		}
		break;

	    case 3:
		for (j0 = Bottom1; j0 <= Top1; j0++)
            	for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
               		offset    = this->NodeOffset(3, j0, j1, t);
			SliceL    = (double*) data + offset;
			auxSliceL = (double*) auxSlice + offset;

	    	for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
			maxDiff = GetIndexStep(auxSlice, 3, j0, j1, j2, t);
			SliceL[j2] = DrSmoothMax(auxSliceL[j2], argument,
					maxDiff);
		}
		}
		break;

	    }
	}

	break;

	//----------------------------------------------
	// Smooth operation: Up binary
	//----------------------------------------------

	case GEQ:
	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x)	x = (x >= argument ? 1e0 : 0e0)
		DO_LOOP
#undef	OPERATION
	} else {
		THROW_NA;
	}

	break;

	//----------------------------------------------
	// Smooth operation: Down binary
	//----------------------------------------------

	case LEQ:
	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x)	x = (x >= argument ? 0e0 : 1e0)
		DO_LOOP
#undef	OPERATION
	} else {
		THROW_NA;
	}
	break;

	//----------------------------------------------
	// Specific operation: extract state variable
	//
	//----------------------------------------------

	case STVAR:
	{
		int	j0, j1, j2, l;
		double	jump[NFMAX][NFMAX],	// Jump coeff [0<i<=j<nbFactor]
			du;

		// Set the slice to the stochasic variable with
		// index given by the argument.
		l = (int) argument;
		if ((l < 0) || (l > sliceDim)) {
			throw KFailure("%s: cannot set stochastic variable %d "
				"on slice (sliceDim=%d, nFact=%d).\n",
				routine, l, sliceDim, nFact);
		}

		// Jump size
		du = sqrt (JUMPCOEFF * LengthJ[t-1]);
		for (i=0; i<nFact; i++)
		for (j=0; j<=i   ; j++) {
			jump[i][j] = mAweight[t-1][AIDX(i, j)] * du;
		}

		switch (this->mNbFactor) {
        	case 1: 
			offset = NodeOffset(1, 0, 0, t);
			SliceL = Slice + offset;
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				SliceL[j0] = j0 * jump[0][0];
			}
			break;
        	case 2:
			switch (l) {
			case 0:
			    for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = NodeOffset(2, j0, 0, t);
				SliceL = Slice + offset;
				for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
					SliceL[j1] = j0*jump[0][0];
				}
			    }
			    break;
			case 1:
			    for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = NodeOffset(2, j0, 0, t);
				SliceL = Slice + offset;
				for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
					SliceL[j1] = j0*jump[1][0] + j1*jump[1][1];
				}
			    }
			    break;
			}
			break;
        	case 3:
			switch (l) {
			case 0:
			    for (j0 = Bottom1;     j0 <= Top1;     j0++)
			    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				offset = NodeOffset(3, j0, j1, t);
				SliceL = Slice + offset;
				for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
					SliceL[j2] = j0*jump[0][0];
				}
			    }
			    break;

			case 1:
			    for (j0 = Bottom1;     j0 <= Top1;     j0++)
			    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				offset = NodeOffset(3, j0, j1, t);
				SliceL = Slice + offset;
				for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
					SliceL[j2] = j0*jump[1][0] + j1*jump[1][1];
				}
			    }
			    break;

			case 2:
			    for (j0 = Bottom1;     j0 <= Top1;     j0++)
			    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				offset = NodeOffset(3, j0, j1, t);
				SliceL = Slice + offset;
				for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
					SliceL[j2] = j0*jump[2][0] + j1*jump[2][1]
							+ j2*jump[2][2];
				}
			    }
			    break;
			}
			break;

		default:
			throw KFailure("%s: cannot set stochastic variable %d "
				"on slice %s (sliceDim=%d, nFact=%d).\n",
				routine, l, sliceDim, nFact);
		} 
	
	

	}
	break;
	default:
		throw KFailure("%s: undefined instruction `%s'.\n",
			routine, oper);
	}

#undef	DO_LOOP
}

//--------------------------------------------------------------
// Performs a unary operation on a timeslice.
//

void
KMrNTree::sliceUnaryOper(
	double *data,			// (B) slice data
	int sliceDim,			// (I) slice dimension
	const double *data1,		// (I) slice oper argument (same dim)
	KOper oper)			// (I) operation type
{
static	char	routine[] = "KMrNTree::sliceUnaryOper";
register int	i, j, k;		// Node indices
	int     Top1, Bottom1;		// Tree limits (1rst dim)
 	int     *Top2, *Bottom2;	// Tree limits (2nd dim) 
	int     **Top3, **Bottom3;	// Tree limits (3rd dim)
	int	offset;			// Node offset

	int	t = TPIdxCurrent();	// Current time point
	double	*Slice1, *Slice1L,
		*Slice2, *Slice2L;

	double	maxDiff;
	double	*auxSlice = this->NewPrice;	// Auxiliary slice
	double	*auxSliceL;

	Top1    = this->mTop1[t];		// Limits
	Top2    = this->mTop2[t];
	Top3    = this->mTop3[t];
	Bottom1 = this->mBottom1[t];
	Bottom2 = this->mBottom2[t];
	Bottom3 = this->mBottom3[t];

	Slice1 = (double*) data;		// More convenient
	Slice2 = (double*) data1;


	// The loop is the same for all operations, so make it a macro
	// It is s little horrible, but efficient 

#define	DO_LOOP \
	switch (sliceDim) { \
        case 1: \
		offset = NodeOffset(1, 0, 0, t); \
		Slice1L = Slice1 + offset; \
		Slice2L = Slice2 + offset; \
		for (i = Bottom1; i <= Top1; i ++) { \
			OPERATION(Slice1L[i], Slice2L[i]); \
		} \
		break; \
        case 2: \
		for (i = Bottom1; i <= Top1; i ++) {\
                	offset = NodeOffset(2, i, 0, t); \
			Slice1L = Slice1 + offset; \
			Slice2L = Slice2 + offset; \
			for (j = Bottom2[i]; j <= Top2[i]; j++) { \
				OPERATION(Slice1L[j], Slice2L[j]); \
			} \
		} \
		break; \
        case 3: \
		for (i = Bottom1; i <= Top1; i ++) \
                for (j = Bottom2[i]; j <= Top2[i]; j++) { \
		    offset = NodeOffset(3, i, j, t); \
		    Slice1L = Slice1 + offset; \
		    Slice2L = Slice2 + offset; \
		    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++) { \
				OPERATION(Slice1L[k], Slice2L[k]); \
		    } \
		} \
		break; \
	default: \
		KFailure("%s: invalid slice dimension.\n", routine); \
	} 


	// Perform actual operation

	switch (oper) {
	//----------------------------------------------
	// Standard arithmetic operations
	//----------------------------------------------

	case COPY:
#undef	OPERATION
#define	OPERATION(x, y)	x = y
		DO_LOOP
#undef	OPERATION
		break;
	case ADD:
#undef	OPERATION
#define	OPERATION(x, y)	x += y
		DO_LOOP
#undef	OPERATION
		break;
	case SUB:
#undef	OPERATION
#define	OPERATION(x, y)	x -= y
		DO_LOOP
#undef	OPERATION
		break;
	case MULT:
#undef	OPERATION
#define	OPERATION(x, y)	x *= y
		DO_LOOP
#undef	OPERATION
		break;
	case DIV:
#undef	OPERATION
#define	OPERATION(x, y)	x /= y
		DO_LOOP
#undef	OPERATION
		break;

	//----------------------------------------------
	// Smooth operation: MIN
	//----------------------------------------------

	case MIN:
	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x, y)	x = (x < y ? x : y)
		DO_LOOP
#undef	OPERATION
	} else {
	    // X1=MIN(X1,X2) calc by X1 = X2 - MAX(X2-X1,0)

	    // Copy X2-X1 into auxiliary slice
	    sliceUnaryOper(auxSlice, sliceDim, Slice2, COPY);
	    sliceUnaryOper(auxSlice, sliceDim, Slice1, SUB);

	    // Perform smooth max
	    switch (sliceDim) {
            case 1:
		offset = NodeOffset(1, 0, 0, t);
		Slice1L = Slice1 + offset;
		Slice2L = Slice2 + offset;
		auxSliceL = auxSlice + offset;
		for (i = Bottom1; i <= Top1; i ++) {
			maxDiff = GetIndexStep(auxSlice, 1, i, 0, 0, t);
			Slice1L[i] = Slice2L[i]
				- DrSmoothMax(auxSliceL[i], 0e0, maxDiff);
		}
		break;
            case 2:
		for (i = Bottom1; i <= Top1; i ++) {
                	offset = NodeOffset(2, i, 0, t);
			Slice1L = Slice1 + offset;
			Slice2L = Slice2 + offset;
			auxSliceL = auxSlice + offset;
		for (j = Bottom2[i]; j <= Top2[i]; j++) {
			maxDiff = GetIndexStep(auxSlice, 2, i, j, 0, t);
			Slice1L[j] = Slice2L[j]
				- DrSmoothMax(auxSliceL[j], 0e0, maxDiff);
		}
		}
		break;
            case 3:
		for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++) {
			offset = NodeOffset(3, i, j, t);
			Slice1L = Slice1 + offset;
			Slice2L = Slice2 + offset;
			auxSliceL = auxSlice + offset;
		for (k = Bottom3[i][j]; k <= Top3[i][j]; k++) {
			maxDiff = GetIndexStep(auxSlice, 3, i, j, k, t);
			Slice1L[k] = Slice2L[k]
				- DrSmoothMax(auxSliceL[k], 0e0, maxDiff);

		}
		}
		break;
	    default:
		KFailure("%s: invalid slice dimension.\n", routine);
	    } 
	}

	break;

	//----------------------------------------------
	// Smooth operation: MAX
	//----------------------------------------------

	case MAX:
	if (mSmoothFact == 0e0) {
#undef	OPERATION
#define	OPERATION(x, y)	x = (x > y ? x : y)
		DO_LOOP
#undef	OPERATION
	} else {
	    // X1=Max(X1,X2) calc by X1 += MAX(X2-X1,0)


	    // Copy X2-X1 into auxiliary slice
	    sliceUnaryOper(auxSlice, sliceDim, Slice2, COPY);
	    sliceUnaryOper(auxSlice, sliceDim, Slice1, SUB);

	    // Perform smooth max
	    switch (sliceDim) {
            case 1:
		offset = NodeOffset(1, 0, 0, t);
		Slice1L = Slice1 + offset;
		Slice2L = Slice2 + offset;
		auxSliceL = auxSlice + offset;
		for (i = Bottom1; i <= Top1; i ++) {
			maxDiff = GetIndexStep(auxSlice, 1, i, 0, 0, t);
			Slice1L[i] += DrSmoothMax(auxSliceL[i], 0e0,
					maxDiff);
		}
		break;
            case 2:
		for (i = Bottom1; i <= Top1; i ++) {
                	offset = NodeOffset(2, i, 0, t);
			Slice1L = Slice1 + offset;
			Slice2L = Slice2 + offset;
			auxSliceL = auxSlice + offset;
		for (j = Bottom2[i]; j <= Top2[i]; j++) {
			maxDiff = GetIndexStep(auxSlice, 2, i, j, 0, t);
			Slice1L[j] += DrSmoothMax(auxSliceL[j], 0e0,
					maxDiff);
		}
		}
		break;
            case 3:
		for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++) {
			offset = NodeOffset(3, i, j, t);
			Slice1L = Slice1 + offset;
			Slice2L = Slice2 + offset;
			auxSliceL = auxSlice + offset;
		for (k = Bottom3[i][j]; k <= Top3[i][j]; k++) {
			maxDiff = GetIndexStep(auxSlice, 3, i, j, k, t);
			Slice1L[k] += DrSmoothMax(auxSliceL[k], 0e0,
					maxDiff);

		}
		}
		break;
	    default:
		KFailure("%s: invalid slice dimension.\n", routine);
	    } 
	}

	break;

	//----------------------------------------------
	//
	//----------------------------------------------
	default:
		throw KFailure("%s: undefined instruction `%s'.\n",
			routine, oper);
	}
	return;
#undef	DO_LOOP
}

//--------------------------------------------------------------
// 
//

void
KMrNTree::sliceSpecialOper(
	double *data,			// (I) slice data
	int sliceDim,			// (I) slice dimension
	const char *what,		// (I) operaton type
	...)				// (I) arguments ...
{
static	char	routine[] = "KMrNTree::sliceSpecialOper";

	int	t = TPIdxCurrent();	// Current time point
	int	j0, j1, j2, offset;
	int	Top1      = this->mTop1[t],
		Bottom1   = this->mBottom1[t],
 		*Top2     = this->mTop2[t],
		*Bottom2  = this->mBottom2[t],
		**Top3    = this->mTop3[t],
		**Bottom3 = this->mBottom3[t];
	double	*dataL,		// slice pointers
		*sliceL;

	va_list	ap;


	if (!strcmp(what, "sum")) {
		// Sums all value in a slice (inner limits): prototype is
		// sliceSpecialOper(double *data, "sum", double *value)
		//
		va_start(ap, what);
		double	*value = va_arg(ap, double*),
			val = 0e0;
		va_end(ap);

		switch (sliceDim) {
        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);
			sliceL = (double*) data + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				val += sliceL[j0];
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);
				sliceL = (double*) data + offset;
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				val += sliceL[j1];
			}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);
				sliceL = (double*) data + offset;

	    		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
				val += sliceL[j2];
			}
			}
			break;

		}

		*value = val;

	} else if (!strcmp(what, "c*s")) {
		//------------------------------------------------------
		// Multiplies a slice by a scalar and puts the result
		// in slice "data".
		// sliceSpecialOper(double *data, "c*s",
		//		double scalar,
		//		double *slice);
		//------------------------------------------------------

		va_start(ap, what);
		double	scalar = va_arg(ap, double);
		double	*slice = va_arg(ap, double*);
		va_end(ap);

		switch (sliceDim) {
        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);
			dataL  = (double*) data  + offset;
			sliceL = (double*) slice + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				dataL[j0] = scalar * sliceL[j0];
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);
				dataL  = (double*) data  + offset;
				sliceL = (double*) slice + offset;
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				dataL[j1] = scalar * sliceL[j1];
			}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);
				dataL  = (double*) data  + offset;
				sliceL = (double*) slice + offset;

	    		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
				dataL[j2] = scalar * sliceL[j2];
			}
			}
			break;

		}

	} else if (!strcmp(what, "minmax")) {
		//------------------------------------------------------
		// Sums all value in a slice (inner limits): prototype is
		// sliceSpecialOper(double *data, "sum", double *value)
		//------------------------------------------------------

		va_start(ap, what);
		double	*valueMin = va_arg(ap, double*);
		double	*valueMax = va_arg(ap, double*);
		double	val;
		va_end(ap);

		*valueMin = DBL_MAX;
		*valueMax = -DBL_MAX;

		switch (sliceDim) {
        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);
			sliceL = (double*) data + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				val = sliceL[j0];
				*valueMin = MIN(val, *valueMin);
				*valueMax = MAX(val, *valueMax);
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);
				sliceL = (double*) data + offset;
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				val = sliceL[j1];
				*valueMin = MIN(val, *valueMin);
				*valueMax = MAX(val, *valueMax);
			}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);
				sliceL = (double*) data + offset;

	    		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++) {
				val = sliceL[j2];
				*valueMin = MIN(val, *valueMin);
				*valueMax = MAX(val, *valueMax);
			}
			}
			break;

		}

	} else {
		throw KFailure("%s: cannot process operation `%s' on slice.\n",
			routine, what);
	}

}

//==============================================================
//
//	Smoothing operations
//
//==============================================================

//--------------------------------------------------------------
// GetIndexStep:
// Obtains the difference between values of an index at adjacent nodes.
// This value is then used in the smoothing algorithm (Smooth_Step).
// This function hides the dimension generality from  the  caller, 
// but it has no means to check that the adequate amount of space has 
// been allocated under the void * being passed.
//

double
KMrNTree::GetIndexStep(
	double     *Index,      /* (I) Index pointer       */
	int         Dim,        /* (I) Index dimension     */
	int         i,          /* (I) Node indices        */
	int         j,
	int         k,
	int         t)          /* (I) Current time point  */
{

    double  *IndexL;            /* Local pointer */

    double  IndexStep;          /* Output index step */
    double  IndexVal;           /* Index value at mid node */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */



    Top1    = mTop1[t];
    Top2    = mTop2[t];
    Top3    = mTop3[t];
    Bottom1 = mBottom1[t];
    Bottom2 = mBottom2[t];
    Bottom3 = mBottom3[t];


    IndexStep = ERROR;  /* To avoid division by 0 */
                
    switch (Dim)
    {
        case 1:
        {
            IndexL = Index + NodeOffset(1, 0, 0, t);

            if (i > Bottom1)
                IndexStep = MAX (IndexStep, fabs (IndexL[i-1] - IndexL[i]));
            if (i < Top1)                                                           
                IndexStep = MAX (IndexStep, fabs (IndexL[i+1] - IndexL[i]));

            break;
        }                    
        case 2:
        {
            IndexL = Index + NodeOffset(2, i, 0, t);

            IndexVal = IndexL[j];

            if (j > Bottom2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j-1] - IndexVal));
            }
            if (j < Top2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                IndexL = Index + NodeOffset(2, i-1, 0, t);

                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }
            if (i < Top1)
            {         
                IndexL = Index + NodeOffset(2, i+1, 0, t);

                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {                                                  
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }

            break;
        }                    
        case 3:
        {
            IndexL = Index + NodeOffset(3, i, j, t);

            IndexVal = IndexL[k];

            if (k > Bottom3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k-1] - IndexVal));
            }
            if (k < Top3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexL = Index + NodeOffset(3, i-1, j, t);

                    if ((k>=Bottom3[i-1][j])&&(k<=Top3[i-1][j]))
                    {
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (i < Top1)    
            {
                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {
                    IndexL = Index + NodeOffset(3, i+1, j, t);

                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {                                                       
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (j > Bottom2[i])
            {
                IndexL = Index + NodeOffset(3, i, j-1, t);

                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
            if (j < Top2[i])
            {
                IndexL = Index + NodeOffset(3, i, j+1, t);

                if ((k>=Bottom3[i][j+1])&&(k<=Top3[i][j+1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
        }                    
        default:
        {
            break;
        }                    
    }  /* switch */


    return (IndexStep);

}


/*****  DrSmoothStep  ********************************************************/
/*
 *      Smooth version of the step function.
 *      Returns DownValue           if Index <= barrier - step
 *              'smooth curve'      if barrier - step < Index < barrier + step
 *              UpValue             if Index >= barrier + step
 *
 */
static	double DrSmoothStep(
	double  UpValue,  /* (I) used if Index > Barrier + step  */
	double  DownValue,/* (I) used if Index < Barrier - step  */
	double  Index,    /* (I) index value                     */
	double  Barrier,  /* (I) barrier value                   */
	double  Step)     /* (I) smooth between Barrier +/- step */
{
    double smoothValue, xSquare, x;

    Step = ABS(Step);

    /* do the limiting case when step = 0.0 */
    if (Step < DBL_EPSILON)
    {
        smoothValue = (Index < Barrier) ? DownValue : UpValue;
        return (smoothValue);
    }

    /* do the smooth case */
    x = (Index - Barrier)/Step;

    if (x <= -1.0)
    {
        smoothValue = DownValue;
    }
    else if (x >= 1.0)
    {
        smoothValue = UpValue;
    }
    else
    {
        xSquare = x * x;
        smoothValue = ((3.*xSquare-10.) * xSquare + 15.) * x *.0625 +.5;
        smoothValue = DownValue + (UpValue - DownValue) * smoothValue;
    } 
    return (smoothValue);

}


/*****  DrSmoothMax  ********************************************************/
/*
 *      Smooth version of the max(a,b)
 *      This is the based on the integral of the DrSmoothStep function
 *      Returns b                   if a <= b - step
 *              'smooth curve'      if b - step < a < b + step
 *              a                   if a >= b + step
 *
 */
static	double  DrSmoothMax(
	double  a,     /* (I) Argument 1                 */
	double  b,     /* (I) Argument 2                 */
	double  Step)  /* (I) smooth in [-step, step]    */
{
    double smoothValue, x, xSquare;

    Step = ABS(Step);

    /* do the limiting case when step = 0.0 */
    if (Step < DBL_EPSILON)
    {
        smoothValue = (a < b) ? b : a;
        return (smoothValue);
    }

    /* do the smooth case */
    x = (a - b)/Step;

    if (x <= -1.0)
    {
        smoothValue = b;
    }
    else if (x >= 1.0)
    {
        smoothValue = a;
    }
    else
    {
        xSquare = x * x;
        smoothValue = ((xSquare-5.) * xSquare + 15.) * xSquare * 0.03125 + 
                      0.5 * x + 0.15625;
	    smoothValue = Step * smoothValue + b;
    } 
    return (smoothValue);
}




//==============================================================
//
//	TSlice virtual functions.
//
//==============================================================


//--------------------------------------------------------------
// Allocate a slice of proper dimension.
// Name and curve index are NOT specified.

KTSlice&
KMrNTree::TSliceCreate(KTSlice &ts)
{
static	char	routine[] = "KMrNTree::TSliceCreate";

	return TSliceDimCreate(ts, mNbFactor);

}


//--------------------------------------------------------------
// Allocate a slice of given dimension.
// Name and curve index are NOT specified.

KTSlice&
KMrNTree::TSliceDimCreate(KTSlice &ts, int nDim)
{
static	char	routine[] = "KMrNTree::TSliceDimCreate";

	int	sliceDim;

	// Check that we have anything to do 
	if (!ts.IsEmpty()) return(ts);

	if (nDim < 1 || nDim > 3 )
		throw KFailure("%s: invalid slice dimension(%d).\n",
				routine, nDim);

	sliceDim = nDim;

	ts.SetSliceDim(sliceDim);

	// Slice allocation
	if (this->mWidth != NULL)	// tree limit has been set up
		ts.mData = sliceNew(sliceDim);
	else
		throw KFailure("%s: can NOT allocate slice before the "
			       "tree is set up.\n", routine);

	// Set timepoint
	ts.SetTpIdx(this->tpIdxCurrent);

	if (ts.mData == NULL)
		throw KFailure("%s: memory allocation failure.\n", routine);

	return(ts);
}


//--------------------------------------------------------------
//

void
KMrNTree::TSliceDestroy(KTSlice &ts)
{
static	char	routine[] = "KMrNTree::sliceDestroy";

	// Check that we have anything to do 
	if (ts.mData == NULL) return;

	sliceDelete((double*)ts.mData);
	ts.mData = NULL;

	return;
}


//--------------------------------------------------------------
//

double
KMrNTree::TSliceGetCenter(KTSlice &ts)
{
static	char	routine[] = "KMrNTree::TSliceGetCenter";
	int	offset;

	int	sliceDim;

	if (ts.IsEmpty()) {
		throw KFailure("%s: slice `%s' is empty.\n",
			routine, ts.GetSliceName().c_str());
	}


	if (ts.GetTpIdx() != 0) {
		throw KFailure("%s: slice `%s' has nonzero "
			"timepoint index (%d).\n",
			routine, ts.GetSliceName().c_str(), ts.GetTpIdx());
	}


	// Record prices
	sliceDim = ts.GetSliceDim();

	offset = NodeOffset(sliceDim, 0, 0, 0);

	return (((double*) ts.mData) + offset)[0];
}

//--------------------------------------------------------------
//

void
KMrNTree::TSliceSetCenter(KTSlice &ts, double value)
{
static	char	routine[] = "KMrNTree::TSliceSetCenter";
	int	offset;

	int	sliceDim;

	if (ts.IsEmpty()) {
		throw KFailure("%s: slice `%s' is empty.\n",
			routine, ts.GetSliceName().c_str());
	}


	ts.SetTpIdx(0);


	// Record prices
	sliceDim = ts.GetSliceDim();

	offset = NodeOffset(sliceDim, 0, 0, 0);

	(((double*) ts.mData) + offset)[0] = value;
	return;
}


//--------------------------------------------------------------
// Check for dimension and time point consistency

bool
KMrNTree::TSliceCompare(KTSlice &ts1, KTSlice &ts2, 
			KTSComp compType)
{
static	char	routine[] = "KMrNTree::TSliceCompare";

	// check dimensions are the same
	if (ts1.GetSliceDim() != ts2.GetSliceDim())
		return false;

	// check time points consistency
	switch (compType) {
	case SAME:
		if(ts1.GetTpIdx() != ts2.GetTpIdx())
			return false;
		break;

	case NEXT:
		if(ts1.GetTpIdx() != ts2.GetTpIdx()+1)
			return false;
		break;

	default:
		throw KFailure("%s: unknown comparison type (%d).\n",
				routine, compType);
	}

	return true;
}


//--------------------------------------------------------------
// Performs a scalar operation on a timeslice.
//


KTSlice&
KMrNTree::TSliceScalarOper(
	KTSlice &ts,
	double argument,
	KOper oper)
{
static	char	routine[] = "KMrNTree::TSliceScalarOper";

	TSliceCreate(ts);


	// Check timept consistency
	int	tpIdx = TPIdxCurrent();
	if ((oper == COPY) || (oper == STVAR)) {
		ts.SetTpIdx(tpIdx);
	} else {
	   if (ts.GetTpIdx() != tpIdx) {
		throw KFailure("%s: slice `%s' tpIdx %d != tree tpIdx %d.\n",
			routine,
			ts.GetSliceName().c_str(),
			ts.GetTpIdx(), tpIdx);
	   }
	}

	// Perform the operation
	sliceScalarOper((double*)ts.mData, ts.GetSliceDim(), argument, oper);

	return(ts);
}


//--------------------------------------------------------------
// Performs a unary operation on a timeslice.
//

KTSlice&
KMrNTree::TSliceUnaryOper(
	KTSlice &ts,
	const KTSlice &ts1,
	KOper oper)
{
static	char	routine[] = "KMrNTree::TSliceUnaryOper";
	int	tpIdx = TPIdxCurrent();

	// Ensure slice allocated
	if (ts1.IsEmpty()) 
		throw KFailure("%s: while performing operation `%s' "
			" on slice `%s', argument slice `%s' is empty.\n",
			routine, printOper(oper), 
			ts.GetSliceName().c_str(), 
			ts1.GetSliceName().c_str());


	// Ensure target ts is allocated
	TSliceCreate(ts);

	// Check consistency of two slices
	if(ts.GetSliceDim() != ts1.GetSliceDim()) {
		throw KFailure("%s: while performing operation `%s',"
			" slice `%s' dim %d !=  slice `%s' dim %d.\n",
			routine, printOper(oper),
			ts.GetSliceName().c_str(),  ts.GetSliceDim(),
			ts1.GetSliceName().c_str(), ts1.GetSliceDim());
	}


	// Check timept consistency
	if (ts1.GetTpIdx() != tpIdx) {
		throw KFailure("%s: slice `%s' tpIdx %d != tree tpIdx %d.\n",
			routine,
			ts1.GetSliceName().c_str(),
			ts1.GetTpIdx(), tpIdx);
	}
	if (oper == COPY) {
		ts.SetTpIdx(tpIdx);
		ts.SetCurveIdx(ts1.GetCurveIdx());
	} else {
	   if (ts.GetTpIdx() != tpIdx) {
		throw KFailure("%s: slice `%s' tpIdx %d != tree tpIdx %d.\n",
			routine,
			ts.GetSliceName().c_str(),
			ts.GetTpIdx(), tpIdx);
	   }
	}

	// Perform the operation
	sliceUnaryOper(
		(double*) ts.mData,
		(long) ts.GetSliceDim(),
		(double*) ts1.mData,
		oper);


	return(ts);
}



//--------------------------------------------------------------
// Performs a binary operation on two timeslices.
//

KTSlice TSliceBinaryOper(const KTSlice &ts1, 
			 const KTSlice &ts2, 
			 KOper oper)
{
static	char	routine[] = "TSliceBinayOper";

	KTSlice newSlice(*(ts1.mVTree),
			 format("%s %s %s",
				ts1.GetSliceName().c_str(),
				printOper(oper),
				ts2.GetSliceName().c_str()),
			 ts1.GetCurveName());

	(ts1.mVTree)->TSliceUnaryOper(newSlice,
				      ts1,
				      COPY);

	(ts1.mVTree)->TSliceUnaryOper(newSlice,
				      ts2,
				      oper);

	return(newSlice);

}



//--------------------------------------------------------------
// Performs a slice expansion on a timeslice.
//

KTSlice&
KMrNTree::TSliceExpand(
	int	tpIdx,		// (I) time point index
	const KTSlice &ts1,	// (I) slice in sLoDim dim
	KTSlice &ts)		// (I/0) slice in sHiDim dim
{
static	char	routine[] = "KMrNTree::TSliceExpand";

	// Ensure slice allocated
	if (ts1.IsEmpty()) 
		throw KFailure("%s: while performing operation expansion "
			" on slice `%s', argument slice `%s' is empty.\n",
			routine, 
			ts.GetSliceName().c_str(), 
			ts1.GetSliceName().c_str());

	// Ensure target ts is allocated
	TSliceCreate(ts);

	// Check slice dimensions
	if (ts.GetSliceDim() < ts1.GetSliceDim()) {
		throw KFailure("%s: while performing operation expansion,"
			" slice `%s' dim %d <  slice `%s' dim %d.\n",
			routine, 
			ts.GetSliceName().c_str(),  ts.GetSliceDim(),
			ts1.GetSliceName().c_str(), ts1.GetSliceDim());
	}


	sliceExpand(
		tpIdx,
		ts.GetSliceDim(),
		ts1.GetSliceDim(),
		(double*) ts1.mData,
		(double*) ts.mData);

	return(ts);

}



//--------------------------------------------------------------
// Performs a slice degeneration on a timeslice.
//

KTSlice&
KMrNTree::TSliceDegenerate(
	int	tpIdx,		// (I) time point index
	const KTSlice &ts1,	// (I) slice in sHiDim dim
	KTSlice &ts)		// (I/0) slice in sLoDim dim
{
static	char	routine[] = "KMrNTree::TSliceDegenerate";

	// Ensure slice allocated
	if (ts1.IsEmpty()) 
		throw KFailure("%s: while performing degeneration operation"
			" on slice `%s', argument slice `%s' is empty.\n",
			routine, 
			ts.GetSliceName().c_str(), 
			ts1.GetSliceName().c_str());

	// Ensure target ts is allocated
	TSliceCreate(ts);

	// Check slice dimensions
	if (ts.GetSliceDim() > ts1.GetSliceDim()) {
		throw KFailure("%s: while performing degeneration operation,"
			" slice `%s' dim %d >  slice `%s' dim %d.\n",
			routine,
			ts.GetSliceName().c_str(),  ts.GetSliceDim(),
			ts1.GetSliceName().c_str(), ts1.GetSliceDim());
	}


	sliceDegenerate(
		tpIdx,
		ts1.GetSliceDim(),
		ts.GetSliceDim(),
		(double*) ts1.mData,
		(double*) ts.mData);

	return(ts);

}



//--------------------------------------------------------------
// 
//

void	
KMrNTree::TSlicePut(KTSlice &ts, ostream &os, int minMaxOnly)
{
	int	tpIdx = TPIdxCurrent();	// Current time point
	
	if (!ts.IsEmpty())
	{
	    if (minMaxOnly) {
		double	valueMin, valueMax;

		sliceSpecialOper(
			(double*) ts.mData,
			(int) ts.GetSliceDim(),
			"minmax",
			&valueMin,
			&valueMax);
		os << format("[%12.8f,%12.8f]", valueMin, valueMax);

	    } else {
		os << ts.GetSliceName() 
	   	<< format(", TPDate = %s:", 
		      	  GtoFormatDate(TPDateCurrent())) << endl;

		slicePrint(
			(double*) ts.mData,
			(int) ts.GetSliceDim(),
			tpIdx,
			FALSE, 			// TRUE=out limits.
			os);			// (I) printed to
	    }

	}
	else 
		os << "Empty Slice";

	os << endl;
}


//--------------------------------------------------------------
// 
//

void
KMrNTree::TSliceSpecialOper(KTSlice &ts, char *what, ...)
{
static	char	routine[] = "KMrNTree::TSliceSpecialOper";

	int	t = TPIdxCurrent();	// Current time point
	int	sliceDim;

	int	j0, j1, j2, offset;
	int	Top1      = this->mTop1[t],
		Bottom1   = this->mBottom1[t],
 		*Top2     = this->mTop2[t],
		*Bottom2  = this->mBottom2[t],
		**Top3    = this->mTop3[t],
		**Bottom3 = this->mBottom3[t];

	double	*Slice = (double*) ts.mData;
	double	*SliceL;

	double	maxDiff;

	va_list ap;
	va_start(ap, what);

	sliceDim = (int) ts.mSliceDim;

    try {

	if (!strcmp(what, "sum")) {	
		double	*value = va_arg(ap, double*);

		sliceSpecialOper(
				(double*)ts.mData,
				sliceDim,
				"sum", 
				value);
	}
	else if (!strcmp(what, "EXTIMESET")) {
		//
		// Hack ! This should not be here !
		// Ask me for explanations !
		//
 
 
		KTSlice *undKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *optKs = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt0Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt1Ks = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *xt2Ks = (KTSlice*) va_arg(ap, KTSlice*);
		double  te     = (double)   va_arg(ap, double);
 
		double  *undTs = (double*) undKs->mData;
		double  *optTs = (double*) optKs->mData;
		double  *xt0Ts = (double*) xt0Ks->mData;
		double  *xt1Ts = (double*) xt1Ks->mData;
		double  *xt2Ts = (double*) xt2Ks->mData;

		double  *undTsL = NULL;
		double  *optTsL = NULL;
		double  *xt0TsL = NULL;
		double  *xt1TsL = NULL;
		double  *xt2TsL = NULL;


		switch (sliceDim) {

        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);

			undTsL = undTs + offset;
			optTsL = optTs + offset;
			xt0TsL = xt0Ts + offset;
			xt1TsL = xt1Ts + offset;
			xt2TsL = xt2Ts + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
			    if (undTsL[j0] > optTsL[j0])
			    {
				xt0TsL[j0] = 1e0;
				xt1TsL[j0] = te;
				xt2TsL[j0] = te*te;
			    }
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);

				undTsL = undTs + offset;
				optTsL = optTs + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

            			for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
			    	    if (undTsL[j1] > optTsL[j1])
			    	    {
					xt0TsL[j1] = 1e0;
					xt1TsL[j1] = te;
					xt2TsL[j1] = te*te;
			    	    }
				}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);

				undTsL = undTs + offset;
				optTsL = optTs + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

	    			for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; 
				     j2++) 
				{
			    	    if (undTsL[j2] > optTsL[j2])
				    {
					xt0TsL[j2] = 1e0;
					xt1TsL[j2] = te;
					xt2TsL[j2] = te*te;
			    	    }
				}
			}
			break;

		}

	}

	else if (!strcmp(what, "IFPOS")) {
	//---------------------------------------------------------
	// Calculates (X1 >= 0 ? X2 : X3) where X0, X1, X2 are
	// time slices
	// TSliceSpecialOper(KTSlice &ts, "IFPOS",
	//	(void*) KTSlice* X0,
	//	(void*) KTSlice* X1,
	//	(void*) KTSlice* X2);
	//
 
	KTSlice *xt0Ks = (KTSlice*) va_arg(ap, KTSlice*);
	KTSlice *xt1Ks = (KTSlice*) va_arg(ap, KTSlice*);
	KTSlice *xt2Ks = (KTSlice*) va_arg(ap, KTSlice*);

	double  *xt0Ts = (double*) xt0Ks->mData;
	double  *xt1Ts = (double*) xt1Ks->mData;
	double  *xt2Ts = (double*) xt2Ks->mData;

	double  *xt0TsL = NULL;
	double  *xt1TsL = NULL;
	double  *xt2TsL = NULL;

	ASSERT_OR_THROW(xt0Ks->mSliceDim == sliceDim);
	ASSERT_OR_THROW(xt1Ks->mSliceDim == sliceDim);
	ASSERT_OR_THROW(xt2Ks->mSliceDim == sliceDim);


	if (mSmoothFact == 0e0) {
		switch (sliceDim) {
       		case 1: 
			offset = this->NodeOffset(1, 0, 0, t);

			SliceL = Slice + offset;
			xt0TsL = xt0Ts + offset;
			xt1TsL = xt1Ts + offset;
			xt2TsL = xt2Ts + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				SliceL[j0] = (xt0TsL[j0] > 0e0 ?
					xt1TsL[j0] : xt2TsL[j0]);
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);

				SliceL = Slice + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
				SliceL[j1] = (xt0TsL[j1] > 0e0 ?
					xt1TsL[j1] : xt2TsL[j1]);
			}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);

				SliceL = Slice + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

	    		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++){
				SliceL[j2] = (xt0TsL[j2] > 0e0 ?
					xt1TsL[j2] : xt2TsL[j2]);
			}
			}
			break;

		}
	} else {
		switch (sliceDim) {
       		case 1: 
			offset = this->NodeOffset(1, 0, 0, t);

			SliceL = Slice + offset;
			xt0TsL = xt0Ts + offset;
			xt1TsL = xt1Ts + offset;
			xt2TsL = xt2Ts + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
			    maxDiff = GetIndexStep(xt0Ts, 1, j0, 0, 0, t);
			    SliceL[j0] = DrSmoothStep(
					xt1TsL[j0], xt2TsL[j0],
					xt0TsL[j0],
					0e0,
					maxDiff);
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);

				SliceL = Slice + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {

			    maxDiff = GetIndexStep(xt0Ts, 2, j0, j1, 0, t);
			    SliceL[j1] = DrSmoothStep(
					xt1TsL[j1], xt2TsL[j1],
					xt0TsL[j1],
					0e0,
					maxDiff);
			}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);

				SliceL = Slice + offset;
				xt0TsL = xt0Ts + offset;
				xt1TsL = xt1Ts + offset;
				xt2TsL = xt2Ts + offset;

	    		for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++){

			    maxDiff = GetIndexStep(xt0Ts, 3, j0, j1, j2, t);
			    SliceL[j2] = DrSmoothStep(
					xt1TsL[j2], xt2TsL[j2],
					xt0TsL[j2],
					0e0,
					maxDiff);
			}
			}
			break;

		}

	}

 
	}	
    else if (!strcmp(what, "STEP")) {
	//---------------------------------------------------------
	// Calculates STEP(X1, X1, X2) where X0, X1, X2 are
	// time slices
	// TSliceSpecialOper(KTSlice &ts, "STEP",
	//	(void*) KTSlice* X0,
	//	(void*) KTSlice* X1,
	//	(void*) KTSlice* X2,
    //  (void*) KKnockIo* mIOWindow,
    //  (void*) KSmooth*  mSmooth);
	//
    
	KTSlice *xt0Ks = (KTSlice*) va_arg(ap, KTSlice*);
	KTSlice *xlbKs = (KTSlice*) va_arg(ap, KTSlice*);
	KTSlice *xubKs = (KTSlice*) va_arg(ap, KTSlice*);

    KKnockIO mIOWindow = *((KKnockIO*) va_arg(ap,KKnockIO*));
    KSmooth  mKOSmooth = *((KSmooth*) va_arg(ap,KSmooth*));
    
//    cout << "mrtslice.STEP " << mIOWindow << '\t' << mKOSmooth << endl;

	double  *xt0Ts = (double*) xt0Ks->mData;
	double  *xlbTs = (double*) xlbKs->mData;
	double  *xubTs = (double*) xubKs->mData;

	double  *xt0TsL = NULL;
	double  *xlbTsL = NULL;
	double  *xubTsL = NULL;

    double  lowerStep, upperStep;

	ASSERT_OR_THROW(xt0Ks->mSliceDim == sliceDim);
	ASSERT_OR_THROW(xlbKs->mSliceDim == sliceDim);
	ASSERT_OR_THROW(xubKs->mSliceDim == sliceDim);


    if (mKOSmooth == NO_SMOOTH) 
    {
        if(mIOWindow == CRX_KNOCK_IN)
        {
            switch (sliceDim) 
            {
            case 1: 
                offset = this->NodeOffset(1, 0, 0, t);

                SliceL = Slice + offset;
                xt0TsL = xt0Ts + offset;
                xlbTsL = xlbTs + offset;
                xubTsL = xubTs + offset;

                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    if (xt0TsL[j0] > xlbTsL[j0] * (1. - BARRIER_TOL) &&
                        xt0TsL[j0] < xubTsL[j0] * (1. + BARRIER_TOL))
                            SliceL[j0] = 1.;
                    else
                            SliceL[j0] = 0.;
                }
                break;

            case 2:
                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    offset = this->NodeOffset(2, j0, 0, t);

                    SliceL = Slice + offset;
                    xt0TsL = xt0Ts + offset;
                    xlbTsL = xlbTs + offset;
                    xubTsL = xubTs + offset;

                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        if (xt0TsL[j1] > xlbTsL[j1] * (1. - BARRIER_TOL) &&
                            xt0TsL[j1] < xubTsL[j1] * (1. + BARRIER_TOL))
                                SliceL[j1] = 1.;
                        else
                                SliceL[j1] = 0.;
                    }
                }
                break;

            case 3:
                for (j0 = Bottom1; j0 <= Top1; j0++)
                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        offset = this->NodeOffset(3, j0, j1, t);

                        SliceL = Slice + offset;
                        xt0TsL = xt0Ts + offset;
                        xlbTsL = xlbTs + offset;
                        xubTsL = xubTs + offset;

                        for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++)
                        {
                            if (xt0TsL[j2] > xlbTsL[j2] * (1. - BARRIER_TOL) &&
                                xt0TsL[j2] < xubTsL[j2] * (1. + BARRIER_TOL))
                                    SliceL[j2] = 1.;
                            else
                                    SliceL[j2] = 0.;
                        }
                    }
                break;
            }// switch
        } // if (mIOWindow == CRX_KNOCK_OUT)
        else if(mIOWindow == CRX_KNOCK_OUT)
        {
            switch (sliceDim) 
            {
            case 1: 
                offset = this->NodeOffset(1, 0, 0, t);

                SliceL = Slice + offset;
                xt0TsL = xt0Ts + offset;
                xlbTsL = xlbTs + offset;
                xubTsL = xubTs + offset;

                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    if (xt0TsL[j0] < xlbTsL[j0] * (1. + BARRIER_TOL) ||
                        xt0TsL[j0] > xubTsL[j0] * (1. - BARRIER_TOL))
                            SliceL[j0] = 1.;
                    else
                            SliceL[j0] = 0.;
                }
                break;

            case 2:
                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    offset = this->NodeOffset(2, j0, 0, t);

                    SliceL = Slice + offset;
                    xt0TsL = xt0Ts + offset;
                    xlbTsL = xlbTs + offset;
                    xubTsL = xubTs + offset;

                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        if (xt0TsL[j1] < xlbTsL[j1] * (1. + BARRIER_TOL) ||
                            xt0TsL[j1] > xubTsL[j1] * (1. - BARRIER_TOL))
                                SliceL[j1] = 1.;
                        else
                                SliceL[j1] = 0.;
                    }
                }
                break;

            case 3:
                for (j0 = Bottom1; j0 <= Top1; j0++)
                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        offset = this->NodeOffset(3, j0, j1, t);

                        SliceL = Slice + offset;
                        xt0TsL = xt0Ts + offset;
                        xlbTsL = xlbTs + offset;
                        xubTsL = xubTs + offset;

                        for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++)
                        {
                            if (xt0TsL[j2] < xlbTsL[j2] * (1. + BARRIER_TOL) ||
                                xt0TsL[j2] > xubTsL[j2] * (1. - BARRIER_TOL))
                                    SliceL[j2] = 1.;
                            else
                                    SliceL[j2] = 0.;
                        }
                    }
                break;
            }// switch
        } // else if (mIOWindow == CRX_KNOCK_OUT)
        else
            throw KFailure("%s: invalid knock window type (%d). \n",
                            routine,
                            mIOWindow);
    } // if (mKOSmooth == NO_SMOOTH)
    else if (mKOSmooth == SINGLE_SMOOTH ||
             mKOSmooth == DOUBLE_SMOOTH)
    {
			   ////////////////////////////////////////////
			   // Given two step functions:
			   // 
			   //               f1
			   // 
			   //            ------------------------
			   //            |         
			   //            |         
			   //    --------         
			   //
			   //
			   //               f2
			   // 
			   //                      --------------
			   //                      |         
			   //                      |         
			   //    ------------------         
			   //
			   //
			   // Knock-in range payoff function
			   //
			   //             (1-f2)*f1
			   //
			   //            -----------
			   //            |         |
			   //            |         |
			   //    --------          -------------- 
			   //
			   //
			   // Knock-out range payoff function
			   //
			   //            1 - (1-f2)*f1
			   //
			   //    --------          -------------- 
			   //            |         |
			   //            |         |
			   //            -----------
			   //
			   //////////////////////////////////////////////
        if (mIOWindow == CRX_KNOCK_IN)
        {        
            switch (sliceDim) 
            {
            case 1: 
                offset = this->NodeOffset(1, 0, 0, t);

                SliceL = Slice + offset;
                xt0TsL = xt0Ts + offset;
                xlbTsL = xlbTs + offset;
                xubTsL = xubTs + offset;

                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    maxDiff = GetIndexStep(xt0Ts, 1, j0, 0, 0, t);
                    lowerStep = DrSmoothStep(
                        1., 
                        0.,
                        xt0TsL[j0],
                        xlbTsL[j0],
                        maxDiff);
                    upperStep = DrSmoothStep(
                        1., 
                        0.,
                        xt0TsL[j0],
                        xubTsL[j0],
                        maxDiff);

                    SliceL[j0] = (1. - upperStep) * lowerStep;
                }
                break;

            case 2:
                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    offset = this->NodeOffset(2, j0, 0, t);

                    SliceL = Slice + offset;
                    xt0TsL = xt0Ts + offset;
                    xlbTsL = xlbTs + offset;
                    xubTsL = xubTs + offset;

                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        maxDiff = GetIndexStep(xt0Ts, 2, j0, j1, 0, t);
                        lowerStep = DrSmoothStep(
                            1., 
                            0.,
                            xt0TsL[j1],
                            xlbTsL[j1],
                            maxDiff);
                        upperStep = DrSmoothStep(
                            1., 
                            0.,
                            xt0TsL[j1],
                            xubTsL[j1],
                            maxDiff);
                        SliceL[j1] = (1. - upperStep) * lowerStep;
                    }
                }
                break;

            case 3:
                for (j0 = Bottom1; j0 <= Top1; j0++)
                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        offset = this->NodeOffset(3, j0, j1, t);

                        SliceL = Slice + offset;
                        xt0TsL = xt0Ts + offset;
                        xlbTsL = xlbTs + offset;
                        xubTsL = xubTs + offset;

                        for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++)
                        {
                            maxDiff = GetIndexStep(xt0Ts, 3, j0, j1, j2, t);
                            lowerStep = DrSmoothStep(
                                1., 
                                0.,
                                xt0TsL[j2],
                                xlbTsL[j2],
                                maxDiff);
                            upperStep = DrSmoothStep(
                                1., 
                                0.,
                                xt0TsL[j2],
                                xubTsL[j2],
                                maxDiff);
                            SliceL[j2] = (1. - upperStep) * lowerStep;
                        }
                    }
                break;
            }// switch
        } // if (mIOWindow == CRX_KNOCK_IN)
        else if (mIOWindow == CRX_KNOCK_OUT)
        {        
            switch (sliceDim) 
            {
            case 1: 
                offset = this->NodeOffset(1, 0, 0, t);

                SliceL = Slice + offset;
                xt0TsL = xt0Ts + offset;
                xlbTsL = xlbTs + offset;
                xubTsL = xubTs + offset;

                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    maxDiff = GetIndexStep(xt0Ts, 1, j0, 0, 0, t);
                    lowerStep = DrSmoothStep(
                        1., 
                        0.,
                        xt0TsL[j0],
                        xlbTsL[j0],
                        maxDiff);
                    upperStep = DrSmoothStep(
                        1., 
                        0.,
                        xt0TsL[j0],
                        xubTsL[j0],
                        maxDiff);

                    SliceL[j0] = 1. - (1. - upperStep) * lowerStep;
                }
                break;

            case 2:
                for (j0 = Bottom1; j0 <= Top1; j0++) 
                {
                    offset = this->NodeOffset(2, j0, 0, t);

                    SliceL = Slice + offset;
                    xt0TsL = xt0Ts + offset;
                    xlbTsL = xlbTs + offset;
                    xubTsL = xubTs + offset;

                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        maxDiff = GetIndexStep(xt0Ts, 2, j0, j1, 0, t);
                        lowerStep = DrSmoothStep(
                            1., 
                            0.,
                            xt0TsL[j1],
                            xlbTsL[j1],
                            maxDiff);
                        upperStep = DrSmoothStep(
                            1., 
                            0.,
                            xt0TsL[j1],
                            xubTsL[j1],
                            maxDiff);
                        SliceL[j1] = 1. - (1. - upperStep) * lowerStep;
                    }
                }
                break;

            case 3:
                for (j0 = Bottom1; j0 <= Top1; j0++)
                    for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) 
                    {
                        offset = this->NodeOffset(3, j0, j1, t);

                        SliceL = Slice + offset;
                        xt0TsL = xt0Ts + offset;
                        xlbTsL = xlbTs + offset;
                        xubTsL = xubTs + offset;

                        for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; j2++)
                        {
                            maxDiff = GetIndexStep(xt0Ts, 3, j0, j1, j2, t);
                            lowerStep = DrSmoothStep(
                                1., 
                                0.,
                                xt0TsL[j2],
                                xlbTsL[j2],
                                maxDiff);
                            upperStep = DrSmoothStep(
                                1., 
                                0.,
                                xt0TsL[j2],
                                xubTsL[j2],
                                maxDiff);
                            SliceL[j2] = 1. - (1. - upperStep) * lowerStep;
                        }
                    }
                break;
            }// switch
        } //else if (mIOWindow == CRX_KNOCK_OUT)
        else
            throw KFailure("%s: invalid knock window type (%d). \n",
                            routine,
                            mIOWindow);
    }// else if (mKOSmooth == SINGLE_SMOOTH || mKOSmooth == DOUBLE_SMOOTH)
    else
            throw KFailure("%s: invalid smooth type (%d). \n",
                            routine,
                            mKOSmooth);
    }else if (!strcmp(what, "COPRODUCT")) {
		//
		// (coslice, slice, val)
		//
 
		KTSlice *sliceKs   = (KTSlice*) va_arg(ap, KTSlice*);
		KTSlice *cosliceKs = (KTSlice*) va_arg(ap, KTSlice*);
		double  *value     = (double*)  va_arg(ap, double*);
 
		double  *sliceTs   = (double*) sliceKs->mData;
		double  *cosliceTs = (double*) cosliceKs->mData;

		double  *sliceTsL   = NULL;
		double  *cosliceTsL = NULL;

		*value = 0e0;

		switch (sliceDim) {

        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);

			sliceTsL   = sliceTs   + offset;
			cosliceTsL = cosliceTs + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				*value += sliceTsL[j0] * cosliceTsL[j0];
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);

				sliceTsL   = sliceTs   + offset;
				cosliceTsL = cosliceTs + offset;

            			for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
					*value += sliceTsL[j1] * cosliceTsL[j1];
				}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);

				sliceTsL   = sliceTs   + offset;
				cosliceTsL = cosliceTs + offset;

	    			for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; 
				     j2++) 
				{
					*value += sliceTsL[j2] * cosliceTsL[j2];
				}
			}
			break;

		}

 
	} else if (!strcmp(what, "TSPRT_MINMAX")) {
		//
		//
 
		double  minVal =  1e30,
			maxVal = -1e30;
 
		double	*valL = NULL;

		ostream *os = (ostream*) va_arg(ap, ostream*);
 
		switch (sliceDim) {

        	case 1: 
			offset = this->NodeOffset(1, 0, 0, t);

			valL = (double*)ts.mData + offset;

			for (j0 = Bottom1; j0 <= Top1; j0++) {
				minVal = MIN(minVal, valL[j0]);
				maxVal = MAX(maxVal, valL[j0]);
			}
			break;

		case 2:
			for (j0 = Bottom1; j0 <= Top1; j0++) {
				offset = this->NodeOffset(2, j0, 0, t);

				valL = (double*)ts.mData + offset;

            			for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
					minVal = MIN(minVal, valL[j1]);
					maxVal = MAX(maxVal, valL[j1]);
				}
			}
			break;

		case 3:
			for (j0 = Bottom1; j0 <= Top1; j0++)
            		for (j1 = Bottom2[j0]; j1 <= Top2[j0]; j1++) {
                		offset = this->NodeOffset(3, j0, j1, t);

				valL = (double*)ts.mData + offset;

	    			for (j2 = Bottom3[j0][j1]; j2 <= Top3[j0][j1]; 
				     j2++) 
				{
					minVal = MIN(minVal, valL[j2]);
					maxVal = MAX(maxVal, valL[j2]);
				}
			}
			break;

		}
 
		(*os) << format("min: %lf  max:%lf",
			 minVal, maxVal);
 
 
        } else {
                throw KFailure("%s: operation `%s' not supported.\n",
				routine, what);
        }

	va_end(ap);
 
    }
    catch (KFailure) {
        va_end(ap);
        throw KFailure("%s: failed", routine);
    }
}




//--------------------------------------------------------------
// 
//

KTSlice&
KMrNTree::TSliceDev(KTSlice &ts, const String &disCurveName)
{
static	char	routine[] = "KMrNTree::TSliceDev";

 try {

	// Dev empty slices does nothing 
	if (ts.mData == NULL) return ts;

	// Don't do anything at last timept
	if (this->tpIdxCurrent == TPNum()) return ts;

	// Check timepoint consistency
	if ((ts.GetTpIdx() != -1) && (ts.GetTpIdx() != tpIdxCurrent+1))
		throw KFailure("%s: attempting to Ev slice `%s' "
			"with TP %d on tree at TP %d.\n",
			routine, ts.GetSliceName().c_str(), 
			ts.GetTpIdx(), tpIdxCurrent);

	// Check index
	///$$$ ... to be done ...

	// Perform EV
	sliceEv(
		(double*) ts.mData,
		(int) ts.mSliceDim,
		(double*) NULL,
		tpIdxCurrent);


	// Set new TP
	ts.SetTpIdx(tpIdxCurrent);

	return (ts);

    }
    catch (KFailure) {
        throw KFailure("%s: failed", routine);
    }
}

//--------------------------------------------------------------
// 
//

KTSlice&
KMrNTree::TSliceEv(KTSlice &ts)
{
static	char	routine[] = "KMrNTree::TSliceEv";

 try {

	// Dev empty slices does nothing 
	if (ts.mData == NULL) return ts;

	// Don't do anything at last timept
	if (this->tpIdxCurrent == TPNum()) return ts;

	// Check timepoint consistency
	if ((ts.GetTpIdx() != -1) && (ts.GetTpIdx() != tpIdxCurrent+1))
		throw KFailure("%s: attempting to Ev slice `%s' "
			"with TP %d on tree at TP %d.\n",
			routine, ts.GetSliceName().c_str(), 
			ts.GetTpIdx(), tpIdxCurrent);

	// Check index
	///$$$ ... to be done ...

	// Perform EV
	sliceEv(
		(double*) ts.mData,
		(int) ts.mSliceDim,
		(double*) NULL,
		tpIdxCurrent);


	// Set new TP
	ts.SetTpIdx(tpIdxCurrent);

	return (ts);

    }
    catch (KFailure) {
        throw KFailure("%s: failed", routine);
    }
}


//--------------------------------------------------------------
// 
//

KTSlice&
KMrNTree::TSliceFwd(KTSlice &ts)
{
static	char	routine[] = "KMrNTree::TSliceFwd";

 try {

	// Fwd empty slices does nothing 
	if (ts.mData == NULL) return ts;

	// Don't do anything at last timept
	if (this->tpIdxCurrent == TPNum()) return ts;

	// Check timepoint consistency
	if ((ts.GetTpIdx() != -1) && (ts.GetTpIdx() != tpIdxCurrent))
		throw KFailure("%s: attempting to forward slice `%s' "
			       "with TP %d on tree at TP %d.\n",
			       routine, ts.GetSliceName().c_str(), 
			       ts.GetTpIdx(), tpIdxCurrent);

	// Perform EV
	sliceFw(
		(double*) ts.mData,
		(int) ts.mSliceDim,
		(double*) NULL,
		tpIdxCurrent);


	// Set new TP
	ts.SetTpIdx(tpIdxCurrent+1);

	return (ts);

    }
    catch (KFailure) {
        throw KFailure("%s: failed", routine);
    }
}




//--------------------------------------------------------------
// Degenerate the slice from higher dimension sHiDim to lower 
// dimension sLoDim. The degenerated value is the sum over 
// the higher sHiDim-sLoDim dimensions  

void
KMrNTree::sliceDegenerate(
	int 	tpIdx,		// (I) time point index
	int    	sHiDim,		// (I) higher slice dimensions
	int    	sLoDim,		// (I) lower degenerated slice dimensions
	double 	*sHiVal,	// (I) slice values in sHiDim dimensions.
	double 	*sLoVal)	// (O) slice values in sLoDim degenrated dim.
{
static	char	routine[] = "KMrNTree::sliceDegenerate";

	int	t = tpIdx;		// more convenient ...
	double	*sHiValL,		//
		*sLoValL;		//


	int	i, j0, j1, j2;	// factor indices

    try {

	// Check that sHiDim >= sLoDim
	if (sHiDim < sLoDim)
		throw KFailure("%s: invalid factor dimensions. "
			       "sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);

	// initialize slice values to zero.
	sliceScalarOper(sLoVal, sLoDim, 0.0e0, COPY);

	switch (sHiDim) {
	case 1:

	   switch (sLoDim) {
	   case 1:
		// sHiDim = sLoDim.  Just reassign the values
		sHiValL = sHiVal + NodeOffset(1, 0, 0, t);
		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);

		for (i = mBottom1[t]; i <= mTop1[t]; i++) { 
			sLoValL[i] = sHiValL[i];
		}

		break;
	   default:
		throw KFailure("%s: invalid factor dimensions. "
			       "sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;

	case 2:  
	   switch (sLoDim) {
	   case 2:
		// sHiDim = sLoDim.  Just reassign the values
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
			sHiValL = sHiVal + NodeOffset(2, j0, 0, t);
			sLoValL = sLoVal + NodeOffset(2, j0, 0, t);

			for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
				sLoValL[i] = sHiValL[i];
			}           
		}			

		break;

	   case 1:
		// sHiDim > sLoDim. Need to degenerate from sHiDim 
		// to sLoDim dimensions, sum over the sHiDim-sLoDim 
		// dimensions.
		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
			sHiValL = sHiVal + NodeOffset(2, j0, 0, t);
			
			for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
				sLoValL[j0] += sHiValL[i];
			}
		}			

		break;
			
	   default:
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;
	
	case 3:
	   switch (sLoDim) {
	   case 3:
		// sHiDim = sLoDim.  Just reassign the values
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		{
			sHiValL = sHiVal + NodeOffset(3, j0, j1, t);
			sLoValL = sLoVal + NodeOffset(3, j0, j1, t);

			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sLoValL[j2] = sHiValL[j2];
			}
		}  // for j0, j1

		break;
	   case 2:
		// sHiDim > sLoDim. Need to degenerate from sHiDim 
		// to sLoDim dimensions, sum over the sHiDim-sLoDim 
		// dimensions.
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
		    sLoValL = sLoVal + NodeOffset(2, j0, 0, t);

		    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		    {
			sHiValL    = sHiVal    + NodeOffset(3, j0, j1, t);
			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sLoValL[j1] += sHiValL[j2];
			}
		    } // for j1
		}  // for j0

		break;
	   case 1:
		// sHiDim > sLoDim. Need to degenerate from sHiDim 
		// to sLoDim dimensions, sum over the sHiDim-sLoDim 
		// dimensions.
		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);

		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		{
			sHiValL = sHiVal + NodeOffset(3, j0, j1, t);
			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sLoValL[j0] += sHiValL[j2];
			}
		}  // for j0

		break;

	   default:
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;

        default:
		throw KFailure("%s: N/A for sHiDim=%d.\n", routine, sHiDim);
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Expand the slice at tpIdx from sLoDim to sHiDim dimensions, 
// in which they are constant in the higer sHiDim-sLoDim dimensions.

void
KMrNTree::sliceExpand(
	int 	tpIdx,		// (I) time point index
	int    	sHiDim,		// (I) higher slice dimensions
	int    	sLoDim,		// (I) lower degenerated slice dimensions
	double 	*sLoVal,	// (I) slice values in sLoDim dimensions.
	double 	*sHiVal)	// (O) slice values in sHiDim dimensions
{
static	char	routine[] = "KMrNTree::sliceExpand";

	int	t = tpIdx;		// more convenient ...
	double	*sHiValL,		//
		*sLoValL;		//


	int	i, j0, j1, j2;	// factor indices


    try {

	// Check that sHiDim >= sLoDim
	if (sHiDim < sLoDim)
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);

	switch (sHiDim) {
	case 1:

	   switch (sLoDim) {
	   case 1:
		// sHiDim = sLoDim.  Just reassign the values
		sHiValL = sHiVal + NodeOffset(1, 0, 0, t);
		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);

		for (i = mBottom1[t]; i <= mTop1[t]; i++) { 
			sHiValL[i] = sLoValL[i];
		}

		break;
	   default:
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;

	case 2:  
	   switch (sLoDim) {
	   case 2:
		// sHiDim = sLoDim.  Just reassign the values
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
			sHiValL = sHiVal + NodeOffset(2, j0, 0, t);
			sLoValL = sLoVal + NodeOffset(2, j0, 0, t);

			for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
				sHiValL[i] = sLoValL[i];
			}           
		}			

		break;

	   case 1:
		// sHiDim > sLoDim. Need to expand slice from sLoDim 
		// to sHiDim dimensions, with constant in the sHiDim-sLoDim dim.

		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);

		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
			sHiValL = sHiVal + NodeOffset(2, j0, 0, t);

			for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
				sHiValL[i] = sLoValL[j0];
			}           
		}			

		break;
			
	   default:
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;
	
	case 3:
	   switch (sLoDim) {
	   case 3:
		// sHiDim = sLoDim.  Just reassign the values
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		{
			sHiValL = sHiVal + NodeOffset(3, j0, j1, t);
			sLoValL = sLoVal + NodeOffset(3, j0, j1, t);

			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sHiValL[j2] = sLoValL[j2];
			}
		}  // for j0, j1

		break;
	   case 2:
		// sHiDim > sLoDim. Need to expand slice from sLoDim 
		// to sHiDim dimensions, with constant in the sHiDim-sLoDim dim.
		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		{
		    sLoValL = sLoVal + NodeOffset(2, j0, 0, t);

		    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		    {
			sHiValL   = sHiVal   + NodeOffset(3, j0, j1, t);
			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sHiValL[j2] = sLoValL[j1];
			}
		    } // for j1
		}  // for j0

		break;
	   case 1:
		// sHiDim > sLoDim. Need to expand slice from sLoDim 
		// to sHiDim dimensions, with constant in the sHiDim-sLoDim dim.
		sLoValL = sLoVal + NodeOffset(1, 0, 0, t);

		for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
		for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
		{
			sHiValL   = sHiVal   + NodeOffset(3, j0, j1, t);
			for (j2 = mBottom3[t][j0][j1]; j2 <= mTop3[t][j0][j1]; 
			     j2++) {
				sHiValL[j2] = sLoValL[j0];
			}
		}  // for j0

		break;

	   default:
		throw KFailure("%s: invalid factor dimensions."
			       " sHiDim=%d, sLoDim=%d.\n", 
				routine, sHiDim, sLoDim);
           } 

	   break;

        default:
		throw KFailure("%s: N/A for sHiDim=%d.\n", routine, sHiDim);
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}








//--------------------------------------------------------------
// 
//

KTSliceNode&
KMrNTree::TSliceNodeBegin(KTSlice &ts, KTSliceNode &tsNode)
{
static	char	routine[] = "KMrNTree::TSliceNodeBegin";

	int	t = TPIdxCurrent();	// Current time point

	int	sliceDim;

	int	j1, j2, j3;
	int	Bottom1   = this->mBottom1[t],
		*Bottom2  = this->mBottom2[t],
		**Bottom3 = this->mBottom3[t];


 try{

	// Check time point consistency
	//
	if(ts.GetTpIdx() != t)
		throw KFailure("%s: time point of the slice (%d) != "
			       "time point of the tree (%d).\n",
				routine, 
				ts.GetTpIdx(),
				t);

	// Slice dimension
	//
	sliceDim = (int) ts.GetSliceDim();

	//
	// Initialize the tsNode
	// Use the first index for slice node end flag (off at initialization),
	// and next "sliceDim" indices for slice indices
	//
	// Slice node index is stored in mIndex, which
	// is initialzed when tsNode.begin is called.
	// Have to clean it up before used with new slice.
	//
	tsNode.mSlice  = &ts;

	if (tsNode.mIndex != NULL)
		delete [] tsNode.mIndex;

	tsNode.mIndex  = new int [sliceDim+1];

	if (tsNode.mIndex == NULL)
		throw KFailure("%s: memory allocation failure.\n", routine);

	((int *)tsNode.mIndex)[0] =  0;


	switch (sliceDim) {

        case 1: 
		j1 = Bottom1;

		((int *)tsNode.mIndex)[1] =  j1;

		break;

	case 2:
		j1 = Bottom1;
		j2 = Bottom2[j1];

		((int *)tsNode.mIndex)[1] =  j1;
		((int *)tsNode.mIndex)[2] =  j2;

		break;

	case 3:
		j1 = Bottom1;
		j2 = Bottom2[j1];
		j3 = Bottom3[j1][j2];

		((int *)tsNode.mIndex)[1] =  j1;
		((int *)tsNode.mIndex)[2] =  j2;
		((int *)tsNode.mIndex)[3] =  j3;

		break;
	}

 
	return(tsNode);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// 
//

KTSliceNode&
KMrNTree::TSliceNodeNext(KTSliceNode &tsNode)
{
static	char	routine[] = "KMrNTree::TSliceNodeNext";

	int	t = TPIdxCurrent();	// Current time point

	int	tsNodeDim;
	int	*tsNodeIdx = (int *)tsNode.mIndex;

	int	j1Current, j2Current, j3Current;
	int	j1Next,    j2Next,    j3Next;

	int	Top1      = this->mTop1[t],
		Bottom1   = this->mBottom1[t],
 		*Top2     = this->mTop2[t],
		*Bottom2  = this->mBottom2[t],
		**Top3    = this->mTop3[t],
		**Bottom3 = this->mBottom3[t];


 try{

	// Check time point consistency
	//
	if(tsNode.mSlice->GetTpIdx() != t)
		throw KFailure("%s: time point of the slice (%d) != "
			       "time point of the tree (%d).\n",
				routine, 
				tsNode.mSlice->GetTpIdx(),
				t);

	// Check for end of slice node
	//
	if(tsNodeIdx[0])
		throw KFailure("%s: end of slice indices. "
			       "Can not get next slice indices.\n",
				routine); 


	tsNodeDim = (int) tsNode.mSlice->GetSliceDim();

	switch (tsNodeDim) {

        case 1: 
		j1Current = tsNodeIdx[1];

		// Set the flag if reaches the end point
		//
		if(j1Current < Top1){
			j1Next = j1Current + 1;	
			tsNodeIdx[1] = j1Next;
		}
		else	// end of slice node
			tsNodeIdx[0] = 1;	// Set the end of slice flag on

		break;

	case 2:
		j1Current = tsNodeIdx[1];
		j2Current = tsNodeIdx[2];

		if (j2Current < Top2[j1Current])  // inside the ellipse
		{
			j2Next = j2Current + 1;

			tsNodeIdx[2] = j2Next;
		}
		else if(j1Current <  Top1 && 
			j2Current == Top2[j1Current]) // at the boundary of 
		{				      // 2nd dim
			j1Next = j1Current + 1;
			j2Next = Bottom2[j1Next];

			tsNodeIdx[1] = j1Next;
			tsNodeIdx[2] = j2Next;
		}
		else if(j1Current == Top1 && 
			j2Current == Top2[j1Current]) // at the end of ellipse 
		{
			tsNodeIdx[0] = 1;	// Set the end of slice flag on
		}
		else
			throw KFailure("%s: invalid slice node indices on"
				       "%d-dimensional slice.",
				       "j1=%d, j2=%d.\n",
					routine, 
					tsNodeDim,
					j1Current,
					j2Current);
			
		break;

	case 3:
		j1Current = tsNodeIdx[1];
		j2Current = tsNodeIdx[2];
		j3Current = tsNodeIdx[3];

		if (j3Current < Top3[j1Current][j2Current]) // inside ellipoide
		{
			j3Next = j3Current + 1;
			tsNodeIdx[3] = j3Next;
		}
		else if(j2Current <  Top2[j1Current] &&   // at the boundary
			j3Current == Top3[j1Current][j2Current]) // of 3rd dim
		{
			j2Next = j2Current + 1;
			j3Next = Bottom3[j1Current][j2Next];

			tsNodeIdx[2] = j2Next;
			tsNodeIdx[3] = j3Next;
		}
		else if(j1Current <  Top1 &&    // at the boundary of ellipse
			j2Current == Top2[j1Current] &&   // axis of 2 & 3
			j3Current == Top3[j1Current][j2Current]) 
		{
			j1Next = j1Current + 1;
			j2Next = Bottom2[j1Next];
			j3Next = Bottom3[j1Next][j2Next];

			tsNodeIdx[1] = j1Next;
			tsNodeIdx[2] = j2Next;
			tsNodeIdx[3] = j3Next;
		}
		else if(j1Current == Top1 &&   // at the end of ellipsoid
			j2Current == Top2[j1Current] &&
			j3Current == Top3[j1Current][j2Current]) 
		{
			tsNodeIdx[0] = 1;
		}
		else
			throw KFailure("%s: invalid slice node indices on "
				       "%d-dimensional slice. "
				       "j1=%d, j2=%d, j3=%d.\n",
					routine, 
					tsNodeDim,
					j1Current,
					j2Current,
					j3Current);
			
		break;
	}
 
	return(tsNode);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// Allow different dimensions for KTSlice and KTSliceNode subject to
// tsDim <= tsNodeDim
//

double&
KMrNTree::TSliceAccessNodeValue(
	KTSlice 	&ts,			// (I) Slice 
	KTSliceNode 	&tsNode)		// (I) Slice node
{
static	char	routine[] = "KMrNTree::TSliceAccessNodeValue";

	int	t = TPIdxCurrent();	// Current time point
	int	j1, j2, j3, offset;

	int     tsDim;
	int     tsNodeDim;

	int	*tsNodeIdx = (int *)tsNode.mIndex;

	double	*dataL;

	double	*value;

 try{
	// Check time point consistency
	//
	if(tsNode.mSlice->GetTpIdx() != t || 
	   ts.GetTpIdx()     != t)
		throw KFailure("%s: inconsistent time points of the slice "
			       "node (%d), slice (%d), and tree (%d).\n",
				routine, 
				tsNode.mSlice->GetTpIdx(),
				ts.GetTpIdx(),
				t);

	// Check NOT end of slice
	//
	if(tsNodeIdx[0])
		throw KFailure("%s: end of slice indices. "
			       "Can not get slice value.\n",
				routine); 

	tsDim = (int)ts.GetSliceDim();
	tsNodeDim = (int)(tsNode.mSlice)->GetSliceDim();

	//
	// Check dimension validity
	//
	if(tsNodeDim < tsDim)
		throw KFailure("%s: dimension of slice node (%d) < "
			       "dimension of slice (%d).\n",
				routine, 
				tsNodeDim,
				tsDim);

	switch (tsDim) {
        case 1: 
		j1 = tsNodeIdx[1];

		offset = this->NodeOffset(1, 0, 0, t);
		dataL = (double*) ts.mData + offset;

		value = &dataL[j1];

		break;

	case 2:
		j1 = tsNodeIdx[1];
		j2 = tsNodeIdx[2];

		offset = this->NodeOffset(2, j1, 0, t);
		dataL = (double*) ts.mData + offset;

		value = &dataL[j2];

		break;

	case 3:
		j1 = tsNodeIdx[1];
		j2 = tsNodeIdx[2];
		j3 = tsNodeIdx[3];

               	offset = this->NodeOffset(3, j1, j2, t);
		dataL = (double*) ts.mData + offset;

		value = &dataL[j3];

		break;
	}

	return *value;

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}




//----------------------------------------------------------------
// Allow different dimensions for KTSlice and KTSliceNode subject to
// tsDim <= tsNodeDim
//

double
KMrNTree::TSliceGetNodeStepMax(
	KTSlice 	&ts,			// (I) Slice 
	KTSliceNode 	&tsNode)		// (I) Slice node
{
static	char	routine[] = "KMrNTree::TSliceGetNodeStepMax";

	int	t = TPIdxCurrent();	// Current time point
	int	j1, j2, j3, offset;

	int     tsDim;
	int     tsNodeDim;

	int	*tsNodeIdx = (int *)tsNode.mIndex;

	int	Top1      = this->mTop1[t],
		Bottom1   = this->mBottom1[t],
 		*Top2     = this->mTop2[t],
		*Bottom2  = this->mBottom2[t],
		**Top3    = this->mTop3[t],
		**Bottom3 = this->mBottom3[t];

	double	*dataL;

	double	indexValue;
	double	indexStepMax;

 try{
	// Check time point consistency
	//
	if(tsNode.mSlice->GetTpIdx() != t || 
	   ts.GetTpIdx()     != t)
		throw KFailure("%s: inconsistent time points of the slice "
			       "node (%d), slice (%d), and tree (%d).\n",
				routine, 
				tsNode.mSlice->GetTpIdx(),
				ts.GetTpIdx(),
				t);

	// Check NOT end of slice
	//
	if(tsNodeIdx[0])
		throw KFailure("%s: end of slice indices. "
			       "Can not get slice value.\n",
				routine); 


	tsDim = (int)ts.GetSliceDim();
	tsNodeDim = (int)(tsNode.mSlice)->GetSliceDim();

	//
	// Check dimension validity
	//
	if(tsNodeDim < tsDim)
		throw KFailure("%s: dimension of slice node (%d) < "
			       "dimension of slice (%d).\n",
				routine, 
				tsNodeDim,
				tsDim);


	// Get the maximum difference of node values in nearest neighbors
	// on each dimension, MAX(value[i-1], value[i], value[i+1]).  
	// At the boundary, only do one side, MAX(value[n], value[n-1]),
	// or MAX(value[-n], value[-n+1]).
	//

	indexStepMax = ERROR;	// ERROR = (1e-7), avoid to be divided by 0

	switch (tsDim) {
        case 1: 
		j1 = tsNodeIdx[1];

		offset = this->NodeOffset(1, 0, 0, t);
		dataL = (double*) ts.mData + offset;

		if(j1 > Bottom1)  // include j1 = Top1
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j1-1]-dataL[j1]));

		if(j1 < Top1)  	  // include j1 = Bottom1
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j1+1]-dataL[j1]));

		break;

	case 2:
		j1 = tsNodeIdx[1];
		j2 = tsNodeIdx[2];

		offset = this->NodeOffset(2, j1, 0, t);
		dataL = (double*) ts.mData + offset;

		indexValue = dataL[j2];

		if(j2 > Bottom2[j1])  	// include j2 = Top2[j1]
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j2-1]-indexValue));

		if(j2 < Top2[j1])  	// include j2 = Bottom2[j1]
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j2+1]-indexValue));

		if(j1 > Bottom1)
		{
			offset = this->NodeOffset(2, j1-1, 0, t);
			dataL = (double*) ts.mData + offset;

			if ((j2>=Bottom2[j1-1])&&(j2<=Top2[j1-1]))
				indexStepMax = MAX(indexStepMax,
					          fabs (dataL[j2]-indexValue));	
		}

		if(j1 < Top1)
		{
			offset = this->NodeOffset(2, j1+1, 0, t);
			dataL = (double*) ts.mData + offset;

			if ((j2>=Bottom2[j1+1])&&(j2<=Top2[j1+1]))
				indexStepMax = MAX(indexStepMax,
					          fabs (dataL[j2]-indexValue));	
		}

		break;

	case 3:
		j1 = tsNodeIdx[1];
		j2 = tsNodeIdx[2];
		j3 = tsNodeIdx[3];

               	offset = this->NodeOffset(3, j1, j2, t);
		dataL = (double*) ts.mData + offset;

		indexValue = dataL[j3];

		// Along the 3rd-axis, both up and down 
		//
		if(j3 > Bottom3[j1][j2]) // include j3 = Top3[j1][j2]
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j3-1]-indexValue));

		if(j3 < Top3[j1][j2])  	 // include j3 = Bottom3[j1][j2]
			indexStepMax = MAX(indexStepMax, 
					   fabs(dataL[j3+1]-indexValue));


		// Along the 1st-axis, both up and down 
		//
		if(j1 > Bottom1)
		{
		    if ((j2>=Bottom2[j1-1])&&(j2<=Top2[j1-1]))
		    {
			offset = this->NodeOffset(3, j1-1, j2, t);
			dataL = (double*) ts.mData + offset;

			if ((j3>=Bottom3[j1-1][j2])&&(j3<=Top3[j1-1][j2]))
				indexStepMax = MAX(indexStepMax,
					           fabs (dataL[j3]-indexValue));	
		    }
		}

		if(j1 < Top1)
		{
		    if ((j2>=Bottom2[j1+1])&&(j2<=Top2[j1+1]))
		    {
			offset = this->NodeOffset(3, j1+1, j2, t);
			dataL = (double*) ts.mData + offset;

			if ((j3>=Bottom3[j1+1][j2])&&(j3<=Top3[j1+1][j2]))
				indexStepMax = MAX(indexStepMax,
					           fabs (dataL[j3]-indexValue));	
		    }
		}


		// Along the 2nd-axis, both up and down 
		//
		if(j2 > Bottom2[j1])
		{
			offset = this->NodeOffset(3, j1, j2-1, t);
			dataL = (double*) ts.mData + offset;

			if ((j3>=Bottom3[j1][j2-1])&&(j3<=Top3[j1][j2-1]))
				indexStepMax = MAX(indexStepMax,
					           fabs (dataL[j3]-indexValue));	
		}

		if(j2 < Top2[j1])
		{
			offset = this->NodeOffset(3, j1, j2+1, t);
			dataL = (double*) ts.mData + offset;

			if ((j3>=Bottom3[j1][j2+1])&&(j3<=Top3[j1][j2+1]))
				indexStepMax = MAX(indexStepMax,
					           fabs (dataL[j3]-indexValue));	
		}

		break;
	}  // switch

	return indexStepMax;

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
//
//
 
void
KMrNTree::TSliceNodePut(KTSliceNode &tsNode, ostream &os)
{
	int     tpIdx = TPIdxCurrent(); // Current time point

	int     tsNodeDim;
	int	*tsNodeIdx = (int *)tsNode.mIndex;
 
        if (tsNode.mSlice == NULL ||
            tsNode.mIndex == NULL)
                throw KFailure("Can not print empty slice node.\n");


	// Check NOT end of slice
	//
	if(tsNodeIdx[0])
		os << "End of slice node" << endl;
	else
	{
        	tsNodeDim = (int) tsNode.mSlice->GetSliceDim();
 
        	os << format("TPDate = %s:",
                     	GtoFormatDate(TPDateCurrent())) << endl;
 
        	switch (tsNodeDim) {
        	case 1:
                	os << "j1 = " << tsNodeIdx[1] << endl;
                	break;
 
        	case 2:
                	os << "j1 = " << tsNodeIdx[1] << endl;
                	os << "j2 = " << tsNodeIdx[2] << endl;
                	break;
 
        	case 3:
                	os << "j1 = " << tsNodeIdx[1] << endl;
                	os << "j2 = " << tsNodeIdx[2] << endl;
                	os << "j3 = " << tsNodeIdx[3] << endl;
                	break;
 
        	}
	}
 
        os << endl;
 
}
