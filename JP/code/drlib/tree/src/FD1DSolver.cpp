#include "edginc/config.hpp"
#include "edginc/FD1DSolver.hpp"
#include "edginc/FD1D.hpp"
#include "edginc/FDUtils.hpp"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------

namespace
{
inline
void interpBoundary( int i, int d, double * x, double * y )
{
    // keep 'd' directional ds/dx constant at 'i'
    int i1 = i + d;
    int i2 = i1 + d;;
    double r = ( x[ i1 ] - x[ i ] ) / ( x[ i2 ] - x[ i1 ] );
    y[ i ] = y[ i1 ] - r * ( y[ i2 ] - y[ i1 ] );
}
}

//-----------------------------------------------------------------------------

inline
void FD1DSolver::solveOneStep(
    double pv, bool nonNegative,
    const TreeSliceEQ & prevSlice,
    TreeSliceEQ & currSlice )
{
    int rangeBot1 = model.range->limits.bot[ 0 ];
    int rangeTop1 = model.range->limits.top[ 0 ];

    // get solver boundaries
    int bot1 = currSlice.getBotBound( 0 );
    int top1 = currSlice.getTopBound( 0 );

    // set values at boundaries
    currSlice[ bot1 ] = bot1 > rangeBot1 && currSlice.getBotValue( 0 ) != DBL_MAX ?
        currSlice.getBotValue( 0 ) : pv * prevSlice[ bot1 ];
    currSlice[ top1 ] = top1 < rangeTop1 && currSlice.getTopValue( 0 ) != DBL_MAX ?
        currSlice.getTopValue( 0 ) : pv * prevSlice[ top1 ];

    int idx = step % 2;
    FDUtils::euroOneStepWithSource(
        bot1, top1,
        aX[ idx ], bX[ idx ], cX[ idx ],
        aX[ 1 - idx ], bX[ 1 - idx ], cX[ 1 - idx ],
        currSlice, prevSlice,
        xSrc, true, 0, 0 );

    // update values at boundaries
    if( bot1 == rangeBot1 )
        interpBoundary( bot1, +1, xVar, currSlice );
    if( top1 == rangeTop1 )
        interpBoundary( top1, -1, xVar, currSlice );

    if( nonNegative )
    {
        for( int i = bot1; i <= top1; ++i )
            currSlice[ i ] = Maths::max( currSlice[ i ], 0. );
    }

    // initialize boundaries for this slice
    currSlice.setBotBound( 0, rangeBot1 );
    currSlice.setTopBound( 0, rangeTop1 );
}

//-----------------------------------------------------------------------------    

void FD1DSolver::calcCoeffPdeFix()
{
    int idx = step % 2;

    for( int j = 0; j < 2; ++j )
    {
        aX[ j ][ bot1 ] = aX[ j ][ top1 ] = 0.;
        bX[ j ][ bot1 ] = bX[ j ][ top1 ] = 1.;
        cX[ j ][ bot1 ] = cX[ j ][ top1 ] = 0.;
    }

    double dx = ( xVar[ top1 ] - xVar[ bot1 ] ) / ( top1 - bot1 );
    double opt1 = .5 / dx;
    double opt2 = 1. / ( dx * dx );

    for( int i = bot1 + 1; i < top1; ++i )
    {
        for( int j = 0; j < 2; ++j )
        {
            double drift   = opt1 * a[ idx ][ i ];
            double itoTerm = opt2 * c[ idx ][ i ];

            double weight = theta - ( idx ^ j );

            aX[ j ][ i ] = weight * ( drift - itoTerm );
            cX[ j ][ i ] = - weight * ( drift + itoTerm );
            bX[ j ][ i ] = 1. + weight * ( 2. * itoTerm -  f[ idx ][ i ] );
        }
    }
}

//-----------------------------------------------------------------------------    

void FD1DSolver::calcCoeffPdeVar()
{
    int idx = step % 2;

    for( int j = 0; j < 2; ++j )
    {
        aX[ j ][ bot1 ] = aX[ j ][ top1 ] = 0.;
        bX[ j ][ bot1 ] = bX[ j ][ top1 ] = 1.;
        cX[ j ][ bot1 ] = cX[ j ][ top1 ] = 0.;
    }

    for( int i = bot1 + 1; i < top1; ++i )
    {
        double dxp = xVar[ i ] - xVar[ i - 1 ];
        double dxn = xVar[ i + 1 ] - xVar[ i ];

        double dxp2 = dxp * dxp;
        double dxn2 = dxn * dxn;

        double opt = 1. / ( dxp * dxn2 + dxn * dxp2 );

        for( int j = 0; j < 2; ++j )
        {
            double drift   =      opt * a[ idx ][ i ];
            double itoTerm = 2. * opt * c[ idx ][ i ];

            double weight = theta - ( idx ^ j );

            aX[ j ][ i ] = weight * ( drift * dxn2 - itoTerm * dxn );
            cX[ j ][ i ] = - weight * ( drift * dxp2 + itoTerm * dxp );
            bX[ j ][ i ] = 1. + weight * ( drift * ( dxp2 - dxn2 ) + itoTerm * ( dxp + dxn ) - f[ idx ][ i ] );
        }
    }
}

//-----------------------------------------------------------------------------    

inline
void FD1DSolver::calcCoeffPde()
{
    // get PDE coefficients for the current step
    int idx = step % 2;
    model.pdeCoeff( step, a[ idx ], c[ idx ], f[ idx ], g[ idx ] );

    if( model.isVariableGrid )
        calcCoeffPdeVar();
    else
        calcCoeffPdeFix();
}

//------------------------------------------------------------------------------

inline
void FD1DSolver::calcSource()
{
    int idx = step % 2;
    for( int i = bot1; i <= top1; ++i )
        xSrc[ i ] = g[ idx ][ i ];
}

//-----------------------------------------------------------------------------    

inline
void FD1DSolver::preCalc()
{
    // preCalc might induce node insertion
    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
        products[ prodIndex ]->preCalc( step );

    // set solver range
    bot1 = model.range->limits.bot[ 0 ];
    top1 = model.range->limits.top[ 0 ];

    // insert nodes into solver source
    xSrc->insertNodes( false );

    // insert nodes into solver coefficients
    for( int j = 0; j < 2; ++j )
    {
        a[ j ]->insertNodes( false );
        c[ j ]->insertNodes( false );
        f[ j ]->insertNodes( false );
        g[ j ]->insertNodes( false );

        aX[ j ]->insertNodes( false );
        bX[ j ]->insertNodes( false );
        cX[ j ]->insertNodes( false );
    }

    int maxBits = sliceCache.size();
    for( int b = 0; b < maxBits; ++b )
        sliceCache[ b ]->insertNodes( false );

    // node insertion into product slices done implicitly by slice operations or explicitly in a solver loop
}

//-----------------------------------------------------------------------------

void FD1DSolver::rollOneStep()
{
    preCalc();

    // update underlying value based on underlying variable
    model.updateValue( step );

    calcCoeffPde();
    calcSource();

    // pv factor between steps
    double pv = 1.;
    if( ! model.isFwdInduction )
    {
        pv = model.discYC->pv(
            model.timeLine->StepDates[ step ], model.timeLine->StepDates[ step + 1 ] );
    }

    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        const vector< TreeSliceSP > & slices = products[ prodIndex ]->getSlicesToDEV();
        int sliceCount = slices.size();
        for( int k = 0; k < sliceCount; ++k )
        {
            TreeSliceLayer * layer = dynamic_cast< TreeSliceLayer * >( slices[ k ].get() );
            int l = layer ? 0 : k;
            int h = layer ? layer->getSlices().size() - 1 : k;
            const vector< TreeSliceSP > & sliceList = layer ? layer->getSlices() : slices;

            for( int j = l; j <= h; ++j )
            {
                TreeSliceEQ & currSlice = static_cast< TreeSliceEQ & >( *sliceList[ j ] );
                // insert nodes into the slice
                currSlice.insertNodes();

                // preserve last step slice in prevSlice
                TreeSliceEQ & prevSlice = *sliceCache[ currSlice.dimBits() ];
                prevSlice.swapValues( currSlice );

                // set boundaries and solve
                solveOneStep(
                    pv, false, // products[ prodIndex ]->isNonNegative(),
                    prevSlice, currSlice );
            }
        }

        //!!! FOR DEBUGGING PURPOSES ONLY
        if( model.DEBUG_DumpToFile.length() && ! products[ prodIndex ]->isElementary() )
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
        if( model.DEBUG_DumpToFile.length() && ! products[ prodIndex ]->isElementary() )
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

void FD1DSolver::roll()
{
    static const string method = "FD1DSolver::Roll";

    // starting step
    step = model.isFwdInduction ? 0 : model.timeLine->NumOfStep;

    preCalc();

    for( int prodIndex = 0; prodIndex < (int)products.size(); ++prodIndex )
    {
        products[ prodIndex ]->update( step,
            model.isFwdInduction ? FDProduct::FWD_0 : FDProduct::BWD_T );

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
    int stepInc = model.isFwdInduction ? 1 : -1;
    for( step += stepInc; step >= 0 && step <= model.timeLine->NumOfStep; step += stepInc )
        rollOneStep();
}

//-----------------------------------------------------------------------------

FD1DSolver::FD1DSolver( const FD1D & model ) :
    model( model ),
    products( model.getProducts() ),
    xVar( model.xVar )
{
    xSrc = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );

    for( int j = 0; j < 2; ++j )
    {
        a[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
        c[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
        f[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
        g[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );

        aX[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
        bX[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
        cX[ j ] = STATIC_POINTER_CAST< TreeSliceEQ >( xVar.clone( false ) );
    }

    int maxBits = 1 << model.range->nDim;
    sliceCache.resize( maxBits );
    for( int b = 0; b < maxBits; ++b )
        sliceCache[ b ] = STATIC_POINTER_CAST< TreeSliceEQ >( model.createSlice( b ) );
}

//-----------------------------------------------------------------------------    

DRLIB_END_NAMESPACE
