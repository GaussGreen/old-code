#include "edginc/config.hpp"
#include "edginc/FD1DMultiSolver.hpp"
#include "edginc/FD1DMulti.hpp"
#include "edginc/FDUtils.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/imsl_fd.h"
#include "edginc/imsl.h"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
FD1DMultiSolver::FD1DMultiSolver( const FD1DMulti & model ) : model( model ) 
{
	nPdes = model.nPdes;
	dim = model.dim;

	bot.resize(nPdes);
	top.resize(nPdes);
	int iPde; 
	for (iPde = 0; iPde < nPdes; ++iPde){
		bot[iPde] = iPde * dim;
		top[iPde] = (iPde + 1) * dim -1;
	}

	xSrc.resize(nPdes*dim);
	xVar.resize(nPdes*dim);
	xValue.resize(nPdes*dim);
    
	int ix;
	for (iPde = 0; iPde < nPdes; ++iPde){
		for (ix = 0; ix < dim; ++ix){
			xVar[iPde*dim + ix] = model.xVar[ix];
			xValue[iPde*dim + ix] = model.xValue[ix];
		}
	}

	int j;
	for(j = 0; j < 2; ++j){
		a[j].resize(nPdes*dim);
        c[j].resize(nPdes*dim);
        f[j].resize(nPdes*dim);
        g[j].resize(nPdes*dim);

		q[j].resize(nPdes*nPdes*dim);
		jump[j].resize(nPdes*nPdes*dim);
		jumpWeight[j].resize(nPdes*nPdes*dim);
		jumpIdx[j].resize(nPdes*nPdes*dim);

		xCoeff[j] = SparseCollection(nPdes*dim, nPdes*dim);
    }

    // Crank-Nicolson scheme by default
    theta = .5;
}

//-----------------------------------------------------------------------------

void FD1DMultiSolver::roll()
{
    static const string method = "FD1DMultiSolver::Roll";

    // starting step
    int direction = model.isFwdInduction ? 1 : -1;
    int step = model.isFwdInduction ? 0 : model.timeLine->NumOfStep;
    
    const FDProductArray & products = model.getProducts();
    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
        products[ prodIndex ]->preCalc( step );

    tempSlices.resize( products.size() );
    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        products[ prodIndex ]->update( step,
            model.isFwdInduction ? FDProduct::FWD_0 : FDProduct::BWD_T );

        // initialize temporary slices
        const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();
        tempSlices[ prodIndex ].resize( slices.size() );
        for( int k = 0; k < (int)slices.size(); ++k )
            tempSlices[ prodIndex ][ k ] = slices[ k ]->clone( false );

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( model.DEBUG_DumpToFile.length() )
        {
            string fileName =
                model.DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "wt" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i ", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }
    }

    // use fully implicit scheme first, so PDE coefficients never used at last step
    theta = 1.;

    // sweep the fd
    for( step += direction; step >= 0 && step <= model.timeLine->NumOfStep; step += direction )
        rollOneStep( step, direction, products );
}

//------------------------------------------------------------------------------

void FD1DMultiSolver::rollOneStep( int step, int direction, const FDProductArray & products )
{
    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
        products[ prodIndex ]->preCalc( step );

    // pv factor between steps
    double pv = 1.;
    if( ! model.isFwdInduction )
    {
        pv = model.discYC->pv(
            model.timeLine->StepDates[ step ], model.timeLine->StepDates[ step + 1 ] );
    }

    // update underlying value based on underlying variable
    model.updateValue( step );

    calcCoeffPde( step );
    calcSource( step );

    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
		const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();
        int sliceCount = slices.size();

        // set boundaries and solve
		for( int k = 0; k < sliceCount; ++k )
        {
			TreeSliceLayer * layer = dynamic_cast< TreeSliceLayer * >( slices[ k ].get() );
            int l = layer ? 0 : k;
            int h = layer ? layer->getSlices().size() - 1 : k;
            const vector< TreeSliceSP > & sliceList = layer ? layer->getSlices() : slices;

            for( int j = l; j <= h; ++j )
            {
				TreeSliceEQ & currSlice = static_cast< TreeSliceEQ & >( *sliceList[ j ] );
                int dimBits = currSlice.dimBits();

                // preserve last step slice in prevSlice
				TreeSliceEQ & prevSlice = static_cast< TreeSliceEQ & >( *sliceList[ j ] );
//                TreeSliceEQ & prevSlice = *sliceCache[ dimBits ];
                prevSlice.swapValues( currSlice );

				setBoundary(
                    step, pv,
                    currSlice, prevSlice );

                euroOneStepWithSource1D(
                    step, false, // products[ prodIndex ]->isNonNegative(),
                    currSlice, prevSlice );

				updateBoundary(
                    step, false, // products[ prodIndex ]->isNonNegative(),
                    currSlice );
            }
        }

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( model.DEBUG_DumpToFile.length() )
        {
            string fileName =
                model.DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "at" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i+", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }

        // update product
        products[ prodIndex ]->update( step, model.isFwdInduction ? FDProduct::FWD : FDProduct::BWD );

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( model.DEBUG_DumpToFile.length() )
        {
            string fileName =
                model.DEBUG_DumpToFile +
                ".p" + Format::toString( prodIndex ) +
                ".csv";

            FILE * file = ::fopen( fileName.c_str(), "at" );

            const TreeSlice & slice = products[ prodIndex ]->getValue( step );

            int bot, top;
            slice.getCalcRange( bot, top );
            double * values = slice.getValues();

            ::fprintf( file, "%.5i-", step );
            for( int i = bot; i <= top; ++i )
                ::fprintf( file, ",\t%16.8f", values[ i ] );
            ::fprintf( file, "\n" );

            ::fclose( file );
        }
    }

    // back to Crank-Nicolson after ferst step
    theta = .5;
}

//-----------------------------------------------------------------------------
void FD1DMultiSolver::setBoundary(
    int step, double pv,
    const TreeSlice & currValue,
    const TreeSlice & prevValue )
{
	for ( int iPde = 0; iPde < nPdes; ++iPde){
//		currValue[ bot[iPde] ] = pv * prevValue[ bot[iPde] ];
//		currValue[ top[iPde] ] = pv * prevValue[ top[iPde] ];
	}
}

//-----------------------------------------------------------------------------
void FD1DMultiSolver::updateBoundary(
    int step, bool nonNegative,
    const TreeSliceEQ & currValue )
{
	for (int iPde = 0; iPde < nPdes; ++iPde){
		interpBoundary( currValue, bot[iPde], +1 );
		interpBoundary( currValue, top[iPde], -1 );

		if( nonNegative )
		{
			currValue[ bot[iPde] ] = Maths::max( currValue[ bot[iPde] ], 0. );
			currValue[ top[iPde] ] = Maths::max( currValue[ top[iPde] ], 0. );
		}
	}
}

//-----------------------------------------------------------------------------
void FD1DMultiSolver::interpBoundary( const TreeSlice & s, int i, int d )
{
    // keep 'd' directional ds/dx constant at 'i'
    double ratio = ( xVar[ i + d ] - xVar[ i ] ) / ( xVar[ i + 2 * d ] - xVar[ i + d ] );
//    s[ i ] = ( 1. + ratio ) * s[ i + d ] - ratio * s[ i + 2 * d ];
}

//-----------------------------------------------------------------------------
void FD1DMultiSolver::euroOneStepWithSource1D(
    int step, bool nonNegative,
    const TreeSliceEQ & currValue,
    const TreeSliceEQ & prevValue ){

    int idx = step % 2;

	double dt = model.timeLine->TradeYrFrac[step+1];

	DoubleArray	source(xSrc);
	FDUtils::euroOneStepWithSource(bot[0], 
									top[nPdes-1], 
									xCoeff[idx],
									xCoeff[1-idx],
									currValue,
									prevValue,
									//&xSrc[0],
									&source[0],
									dt,
									theta);

    if( nonNegative ){
        for( int i = bot[0]; i <= top[nPdes]; ++i )
            currValue[ i ] = Maths::max( currValue[ i ], 0. );
    }
}

//-----------------------------------------------------------------------------    
void FD1DMultiSolver::calcSource( int step )
{
    int idx = step % 2;
	for( int i = bot[0]; i <= top[nPdes-1]; ++i )
		xSrc[i] = g[idx][i];
}

//-----------------------------------------------------------------------------    
void FD1DMultiSolver::calcCoeffPde( int step )
{
    // get PDE coefficients for the current step
    int idx = step % 2;
    model.pdeCoeff(step, a[idx], c[idx], f[idx], g[idx], q[idx], jump[idx]);

	interpJump(step);

    if( model.isVariableGrid )
        calcCoeffPdeVar( step );
    else
        calcCoeffPdeFix( step );
}
//-----------------------------------------------------------------------------    
void FD1DMultiSolver::calcCoeffPdeFix( int step )
{
    int idx = step % 2;

	double value, drift, itoTerm;
	double dx = (xVar[top[0]] - xVar[bot[0]]) / (top[0] - bot[0]);

	xCoeff[idx].empty();

	for (int iRow = 0; iRow < nPdes; ++iRow){
		for (int iCol = 0; iCol < nPdes; ++iCol) {
			int i;
			if (iRow == iCol) {
				// at boundary ...
				i = bot[iRow];
				drift = (a[idx][i] + c[idx][i])/dx;
				xCoeff[idx].push_back(i, i, -drift + f[idx][i]);
				xCoeff[idx].push_back(i, i+1, drift);

				i = top[iRow];
				drift = (a[idx][i] + c[idx][i])/dx;
				xCoeff[idx].push_back(i, i-1, -drift);
				xCoeff[idx].push_back(i, i, drift + f[idx][i]);

				for( i = bot[iRow] + 1; i < top[iRow]; ++i ){
					drift = a[idx][i]/(2.*dx);
					itoTerm = c[idx][i]/(dx*dx);

					value = -drift + itoTerm;
					xCoeff[idx].push_back(i, i-1, value);

					value = -2.*itoTerm + f[idx][i];
					xCoeff[idx].push_back(i, i, value);
						
					value = drift + itoTerm;
					xCoeff[idx].push_back(i, i+1, value);
				}
			} else {
				for ( i = bot[0]; i <= top[0]; ++i ){
					if (Maths::isZero(jump[idx][(iRow * nPdes + iCol) * dim + i])){
						xCoeff[idx].push_back(iRow*dim+i, iCol*dim+i, q[idx][(iRow * nPdes + iCol) * dim + i]);
					} else {
						int iMod = jumpIdx[idx][(iRow * nPdes + iCol) * dim + i];
						xCoeff[idx].push_back(iRow*dim+i, iCol*dim+iMod, 
							q[idx][(iRow * nPdes + iCol) * dim + i]*jumpWeight[idx][(iRow * nPdes + iCol) * dim + i]);
						xCoeff[idx].push_back(iRow*dim+i, iCol*dim+iMod+1, 
							q[idx][(iRow * nPdes + iCol) * dim + i]*(1-jumpWeight[idx][(iRow * nPdes + iCol) * dim + i]));
					}
				}
			}
		}
	}
}

//----------------------------------------------------------------------------- 
void FD1DMultiSolver::calcCoeffPdeVar( int step )
{
 
}

//-----------------------------------------------------------------------------    
void FD1DMultiSolver::interpJump( int step )
{
	int idx = step % 2;

	double dx = (xVar[top[0]] - xVar[bot[0]]) / (top[0] - bot[0]);
	int tmp;
	double jumpSize;
	for (int iRow = 0; iRow < nPdes; ++iRow){
		for (int iCol = 0; iCol < nPdes; ++iCol){
			for (int i = 0; i < dim; ++i){
				if (iRow != iCol){
					jumpSize = jump[idx][(iRow * nPdes + iCol) * dim + i];
					if (Maths::isZero(jumpSize)){
						jumpIdx[idx][(iRow * nPdes + iCol) * dim + i] = i;
						jumpWeight[idx][(iRow * nPdes + iCol) * dim + i] = 1.;
					} else if (Maths::isPositive(jumpSize)){
						tmp = (int)(jump[idx][(iRow * nPdes + iCol) * dim + i]/dx) + i;
						if (tmp > dim - 2) tmp = dim - 2; 
						jumpIdx[idx][(iRow * nPdes + iCol) * dim + i] = tmp;
						// linear interpolation should be done with spot price, however, it seems no significant difference between S and lnS
						//jumpWeight[idx][(iRow * nPdes + iCol) * dim + i] = (xVar[tmp+1]-xVar[i]-jumpSize)/dx;
						jumpWeight[idx][(iRow * nPdes + iCol) * dim + i] = (xValue[tmp+1]-exp(xVar[i]+jumpSize))/(xValue[tmp+1]-xValue[tmp]);
					} else {
						tmp = (int)(jump[idx][(iRow * nPdes + iCol) * dim + i]/dx) + i - 1;
						if (tmp < 0) tmp = 0;
						jumpIdx[idx][(iRow * nPdes + iCol) * dim + i] = tmp;
						// linear interpolation should be done with spot price, however, it seems no significant difference between S and lnS
						//jumpWeight[idx][(iRow * nPdes + iCol) * dim + i] = (xVar[tmp+1]-xVar[i]-jumpSize)/dx;
						jumpWeight[idx][(iRow * nPdes + iCol) * dim + i] = (xValue[tmp+1]-exp(xVar[i]+jumpSize))/(xValue[tmp+1]-xValue[tmp]);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------    
DRLIB_END_NAMESPACE
