//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceBasic.cpp
//
//   Description : Basic tree slice implementation.
//
//   Date        : Sep 8, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

TreeSliceBasic::Range::Range(
    bool transposed,
    int bot1, int top1,
    int bot2, int top2,
    int bot3, int top3 )
    :
    transposed( transposed )
{
    bot[0] = bot1;
    top[0] = top1;
    size[0] = top1 - bot1 + 1;
    if( size[0] <= 0 )
        throw ModelException("TreeSliceBasic::Range", "top1 < bot1 during construction");

    bot[1] = bot2;
    top[1] = top2;
    size[1] = top2 - bot2 + 1;
    if( size[1] <= 0 )
        throw ModelException("TreeSliceBasic::Range", "top2 < bot2 during construction");

    bot[2] = bot3;
    top[2] = top3;
    size[2] = top3 - bot3 + 1;
    if( size[2] <= 0 )
        throw ModelException("TreeSliceBasic::Range", "top3 < bot3 during construction");

    int d = SLICE_BASIC_MAX_DIMS;

    while( --d && size[ d ] <= 1 )
        ;
    nDim = d + 1;

    limits.bot1 = bot1;
    limits.top1 = top1;

    if( nDim > 1 )
    {
        limits.bot2 = new int[ size[0] ] - bot[0];
        limits.top2 = new int[ size[0] ] - bot[0];

        for( int i = bot1; i <= top1; ++i )
        {
            limits.bot2[ i ] = bot2;
            limits.top2[ i ] = top2;
        }

        if( nDim > 2 )
        {
            limits.bot3 = new int *[ size[0] ] - bot[0];
            limits.top3 = new int *[ size[0] ] - bot[0];

            for( int i = bot1; i <= top1; ++i )
            {
                limits.bot3[ i ] = new int[ size[1] ] - bot[1];
                limits.top3[ i ] = new int[ size[1] ] - bot[1];

                for( int j = bot2; j <= top2; ++j )
                {
                    limits.bot3[ i ][ j ] = bot3;
                    limits.top3[ i ][ j ] = top3;
                }
            }
        }
    }
}

TreeSliceBasic::Range::~Range()
{
    if( nDim > 1 )
    {
        if( nDim > 2 )
        {
            for( int i = top[0]; i >= bot[0]; --i )
            {
                delete [] ( limits.top3[ i ] + bot[1] );
                delete [] ( limits.bot3[ i ] + bot[1] );
            }

            delete [] ( limits.top3 + bot[0] );
            delete [] ( limits.bot3 + bot[0] );
        }

        delete [] ( limits.top2 + bot[0] );
        delete [] ( limits.bot2 + bot[0] );
    }
}

int TreeSliceBasic::Range::getSize( int dim ) const
{
#ifdef DEBUG
    if( nDim <= dim )
        throw ModelException("TreeSliceBasic::Range::getSize", "dim out of range");
#endif

    return dim < 0 ? size[0] * size[1] * size[2] : size[ dim ];
}

TreeSliceBasic::TreeSliceBasic( const Range & range, int dimBits ) :
    range( range )
{
    subDimBits = 0;
    valueCount = 1;

    for( int d = 0; d < SLICE_BASIC_MAX_DIMS; ++d )
    {
        int sd = range.transposed ? d : SLICE_BASIC_MAX_DIMS - d - 1;

        if( range.size[ sd ] > 1 && ( dimBits & ( 1 << sd ) ) )
        {
            subDimBits |= ( 1 << sd );
            subSize[ sd ] = valueCount;
            valueCount *= range.size[ sd ];
        }
        else
            subSize[ sd ] = 0;
    }

    values = new double[ valueCount ];

    // initialize slice to 0.
    ::memset( values, 0, valueCount * sizeof( *values ) );

    if( subDim( 0 ) + subDim( 1 ) + subDim( 2 ) <= 1 )
        offsetValues = values - subDim( 0 ) * range.bot[0]
                              - subDim( 1 ) * range.bot[1]
                              - subDim( 2 ) * range.bot[2];
    else
        offsetValues = values;
}

TreeSliceBasic::TreeSliceBasic( const TreeSliceBasic & slice, bool copyValues ) :
    range( slice.range )
{
    subDimBits = slice.subDimBits;
    ::memcpy( subSize, slice.subSize, sizeof( subSize ) );

    valueCount = slice.valueCount;
    values = new double[ valueCount ];

    if( copyValues )
        ::memcpy( values, slice.values, valueCount * sizeof( *values ) );
    else
        ::memset( values, 0, valueCount * sizeof( *values ) );

    offsetValues = values - ( slice.values - slice.offsetValues );
}

void TreeSliceBasic::assume( TreeSliceBasic & slice )
{
    subDimBits = slice.subDimBits;
    ::memcpy( subSize, slice.subSize, sizeof( subSize ) );

    valueCount = slice.valueCount;

    // swap values pointers so they get destroyed appropriately
    double * tmpValues = values;
    values = slice.values;
    slice.values = tmpValues;

    offsetValues = slice.offsetValues;
}

void TreeSliceBasic::expand( int dimBits, bool keepValues )
{
    TreeSliceBasic slice( range, subDimBits | dimBits );
    if( keepValues )
        slice.eval( *this );
    assume( slice );
}

TreeSliceSP TreeSliceBasic::calcSmoothStep() const {
    TreeSliceBasicSP r(new TreeSliceBasic(range, subDimBits));
    r->name = "calcSmoothStep("+name+")";
    ((TreeSlice&)*r) = *this; // force the buffer to the right size

    switch( range.nDim )
    {
        case 1:
        {
            int iStart = subDim( 0 ) * range.limits.bot1;
            int iStop  = subDim( 0 ) * range.limits.top1;
            for( int i = iStart; i <= iStop; ++i ) 
            {
                double v, m=0;
                v = operator()(i);
                if (i!=iStart) m = max(m, fabs( operator()(i-1) - v ));
                if (i!=iStop)  m = max(m, fabs( operator()(i+1) - v ));
                (*r)(i) = m;
            }
            break;
        }
        case 2:
        {
            int iStart = subDim( 0 ) * range.limits.bot1;
            int iStop  = subDim( 0 ) * range.limits.top1;
            for( int i = iStart; i <= iStop; ++i ) 
            {
                int jStart = subDim( 1 ) * range.limits.bot2[ i ];
                int jStop  = subDim( 1 ) * range.limits.top2[ i ];
                for( int j = jStart; j <= jStop; ++j ) 
                {
                    double v, m=0;
                    v = operator()(i, j);
                    if (i!=iStart) m = max(m, fabs( operator()(i-1,j) - v ));
                    if (i!=iStop)  m = max(m, fabs( operator()(i+1,j) - v ));
                    if (j!=jStart) m = max(m, fabs( operator()(i,j-1) - v ));
                    if (j!=jStop)  m = max(m, fabs( operator()(i,j+1) - v ));
                    (*r)(i, j) = m;
                }
            }
            break;
        }
        case 3:
        {
            int iStart = subDim( 0 ) * range.limits.bot1;
            int iStop  = subDim( 0 ) * range.limits.top1;
            for( int i = iStart; i <= iStop; ++i ) 
            {
                int jStart = subDim( 1 ) * range.limits.bot2[ i ];
                int jStop  = subDim( 1 ) * range.limits.top2[ i ];
                for( int j = jStart; j <= jStop; ++j ) 
                {
                    int kStart = subDim( 2 ) * range.limits.bot3[ i ][ j ];
                    int kStop  = subDim( 2 ) * range.limits.top3[ i ][ j ];
                    for( int k = kStart; k <= kStop; ++k ) 
                    {
                        double v, m=0;
                        v = operator()(i, j, k);
                        if (i!=iStart) m = max(m, fabs( operator()(i-1,j,k) - v ));
                        if (i!=iStop)  m = max(m, fabs( operator()(i+1,j,k) - v ));
                        if (j!=jStart) m = max(m, fabs( operator()(i,j-1,k) - v ));
                        if (j!=jStop)  m = max(m, fabs( operator()(i,j+1,k) - v ));
                        if (k!=kStart) m = max(m, fabs( operator()(i,j,k-1) - v ));
                        if (k!=kStop)  m = max(m, fabs( operator()(i,j,k+1) - v ));
                        (*r)(i, j, k) = m;
                    }
                }
            }
            break;
        }
        default:
            throw ModelException("calcSmoothStep", "unsupported slice dimensionality");
    }
    return r;
}

DRLIB_END_NAMESPACE
