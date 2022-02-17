//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceEQ.cpp
//
//   Description : Equity specific tree slice implementation.
//
//   Date        : Sep 8, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

TreeSliceEQ::Range::Range(
    bool transposed,
    int bot1, int top1,
    int bot2, int top2,
    int bot3, int top3 )
    :
    transposed( transposed )
{
    bot[ 0 ] = bot1;
    top[ 0 ] = top1;
    size[ 0 ] = top1 - bot1 + 1;
    limits.bot[ 0 ] = bot1;
    limits.top[ 0 ] = top1;
    if( size[ 0 ] <= 0 )
        throw ModelException("TreeSliceEQ::Range", "top1 < bot1 during construction");

    bot[ 1 ] = bot2;
    top[ 1 ] = top2;
    size[ 1 ] = top2 - bot2 + 1;
    limits.bot[ 1 ] = bot2;
    limits.top[ 1 ] = top2;
    if( size[ 1 ] <= 0 )
        throw ModelException("TreeSliceEQ::Range", "top2 < bot2 during construction");

    bot[ 2 ] = bot3;
    top[ 2 ] = top3;
    size[ 2 ] = top3 - bot3 + 1;
    limits.bot[ 2 ] = bot3;
    limits.top[ 2 ] = top3;
    if( size[ 2 ] <= 0 )
        throw ModelException("TreeSliceEQ::Range", "top3 < bot3 during construction");

    int d = SLICE_EQ_MAX_DIMS;

    while( --d && size[ d ] <= 1 )
        ;
    nDim = d + 1;
}

void TreeSliceEQ::Range::insertNode( int dim, const InsNodeSP & insNode )
{
#ifdef DEBUG
    if( insNode->index < limits.bot[ dim ] || limits.top[ dim ] < insNode->index )
        throw ModelException("TreeSliceEQ::Range::insertNode", "insert node index out of range");
#endif
    ++top[ dim ], ++size[ dim ], ++limits.top[ dim ];
    insNodes[ dim ].push_back( insNode );
}

TreeSliceEQ::TreeSliceEQ( const Range & range, int dimBits ) :
    range( range )
{
    subDimBits = 0;
    valueCount = 1;

    for( int d = 0; d < SLICE_EQ_MAX_DIMS; ++d )
    {
        int sd = range.transposed ? d : SLICE_EQ_MAX_DIMS - d - 1;

        if( range.size[ sd ] > 1 && ( ( dimBits >> sd ) & 1 ) )
        {
            subDimBits |= ( 1 << sd );

            subSize[ sd ] = valueCount;
            insNodeCount[ sd ] = range.insNodes[ sd ].size();

            valueCount *= range.size[ sd ];
        }
        else
        {
            subSize[ sd ] = 0;
            insNodeCount[ sd ] = 0;
        }

        botBound[ sd ] = range.limits.bot[ sd ];
        topBound[ sd ] = range.limits.top[ sd ];
        botValue[ sd ] = 0.;
        topValue[ sd ] = 0.;
    }

    allocCount = valueCount + 2; //!!! REVIEW THIS: temporarily allocate 2 more then needed
    values = (double *)::malloc( allocCount * sizeof( double ) );

    // initialize slice to 0.
    ::memset( values, 0, valueCount * sizeof( double ) );

    if( subDim( 0 ) + subDim( 1 ) + subDim( 2 ) <= 1 )
        offsetValues = values - subDim( 0 ) * range.bot[ 0 ]
                              - subDim( 1 ) * range.bot[ 1 ]
                              - subDim( 2 ) * range.bot[ 2 ];
    else
        offsetValues = values;
}

TreeSliceEQ::TreeSliceEQ( const TreeSliceEQ & slice, bool copyValues ) :
    range( slice.range )
{
    subDimBits = slice.subDimBits;

    ::memcpy( subSize, slice.subSize, sizeof( subSize ) );
    ::memcpy( insNodeCount, slice.insNodeCount, sizeof( insNodeCount ) );

    allocCount = slice.allocCount;
    valueCount = slice.valueCount;
    values = (double *)::malloc( allocCount * sizeof( double ) );

    if( copyValues )
        ::memcpy( values, slice.values, valueCount * sizeof( double ) );
    else
        ::memset( values, 0, valueCount * sizeof( double ) );

    offsetValues = values - ( slice.values - slice.offsetValues );

    ::memcpy( botBound, slice.botBound, sizeof( botBound ) );
    ::memcpy( topBound, slice.topBound, sizeof( topBound ) );
    ::memcpy( botValue, slice.botValue, sizeof( botValue ) );
    ::memcpy( topValue, slice.topValue, sizeof( topValue ) );
}

void TreeSliceEQ::insertNodes( int dim, bool keepValues )
{
    // don't insert if dim is not expanded in the slice
    if( ! subDim( dim ) )
        return;

    const vector< Range::InsNodeSP > & insNodes = range.insNodes[ dim ];

    // check the number of nodes to insert
    int insCount = insNodes.size() - insNodeCount[ dim ];
    if( ! insCount )
        return; // no nodes to insert

    ASSERT( insCount > 0 );

    //!!! implemented for the 1st dimension only
    ASSERT( dim == 0 );

    // check if need allocating more space
    if( ( valueCount += insCount ) > allocCount )
    {
        allocCount += valueCount;
        int offset = offsetValues - values;
        values = (double *)::realloc( values, allocCount * sizeof( double ) );
        offsetValues = values + offset;
    }

    for( int i = 0; i < insCount; ++i )
    {
        const Range::InsNode & insNode = *insNodes[ insNodeCount[ dim ] + i ];

        int index = insNode.index;
        if( index <= botBound[ dim ] )
        {
            ++botBound[ dim ], ++topBound[ dim ];

            if( keepValues )
            {
                ::memmove(
                    offsetValues + botBound[ dim ],
                    offsetValues + botBound[ dim ] - 1,
                    ( topBound[ dim ] - botBound[ dim ] + 1 ) * sizeof( double ) );
            }
        }
        else if( index <= topBound[ dim ] )
        {
            ++topBound[ dim ];

            if( keepValues )
            {
                ::memmove(
                    offsetValues + index + 1,
                    offsetValues + index,
                    ( topBound[ dim ] - index ) * sizeof( double ) );

                offsetValues[ index ] =
                    insNode.coeff1 * offsetValues[ index - 1 ] +
                    insNode.coeff2 * offsetValues[ index + 1 ];
            }
        }
    }

    // all nodes have been inserted
    insNodeCount[ dim ] += insCount;
}

TreeSliceSP TreeSliceEQ::calcSmoothStep() const
{
    TreeSliceEQSP r(new TreeSliceEQ(range, subDimBits));
    r->name = "calcSmoothStep("+name+")";
    ((TreeSlice&)*r) = *this; // force the buffer to the right size

    switch( range.nDim )
    {
        case 1:
        {
            int iStart = subDim( 0 ) * range.limits.bot[ 0 ];
            int iStop  = subDim( 0 ) * range.limits.top[ 0 ];
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
            int iStart = subDim( 0 ) * range.limits.bot[ 0 ];
            int iStop  = subDim( 0 ) * range.limits.top[ 0 ];
            for( int i = iStart; i <= iStop; ++i ) 
            {
                int jStart = subDim( 1 ) * range.limits.bot[ 1 ];
                int jStop  = subDim( 1 ) * range.limits.top[ 1 ];
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
            int iStart = subDim( 0 ) * range.limits.bot[ 0 ];
            int iStop  = subDim( 0 ) * range.limits.top[ 0 ];
            for( int i = iStart; i <= iStop; ++i ) 
            {
                int jStart = subDim( 1 ) * range.limits.bot[ 1 ];
                int jStop  = subDim( 1 ) * range.limits.top[ 1 ];
                for( int j = jStart; j <= jStop; ++j ) 
                {
                    int kStart = subDim( 2 ) * range.limits.bot[ 2 ];
                    int kStop  = subDim( 2 ) * range.limits.top[ 2 ];
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
