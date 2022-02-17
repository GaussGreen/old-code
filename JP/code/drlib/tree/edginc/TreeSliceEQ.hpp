//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceEQ.hpp
//
//   Description : Equity specific tree slice implementation.
//
//   Date        : Sep 8, 2005
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_EQ_HPP
#define TREE_SLICE_EQ_HPP

#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

#define SLICE_EQ_MAX_DIMS 3

// holds slice values
// supports up to SLICE_EQ_MAX_DIMS dimensions (somewhat hard-coded)
class TreeSliceEQ;
typedef refCountPtr< TreeSliceEQ > TreeSliceEQSP;

class TREE_DLL TreeSliceEQ : public TreeSlice
{
public:
    // contains initial range and current limits
    class Range;
    typedef refCountPtr< Range > RangeSP;
    class TREE_DLL Range
    {
    public:
        // allocates new 1D slice range
        static RangeSP create(
            int bot1, int top1,
            bool transposed = false )
        {
            return RangeSP(
                new Range( transposed, bot1, top1 ) );
        }

        // allocates new 2D slice range
        static RangeSP create(
            int bot1, int top1, int bot2, int top2,
            bool transposed = false )
        {
            return RangeSP(
                new Range( transposed, bot1, top1, bot2, top2 ) );
        }

        // allocates new 3D slice range
        static RangeSP create(
            int bot1, int top1, int bot2, int top2, int bot3, int top3,
            bool transposed = false )
        {
            return RangeSP(
                new Range( transposed, bot1, top1, bot2, top2, bot3, top3 ) );
        }

        bool transposed;
        int bot[ SLICE_EQ_MAX_DIMS ];  // for each dimension, starting index
        int top[ SLICE_EQ_MAX_DIMS ];  // for each dimension, ending index
        int size[ SLICE_EQ_MAX_DIMS ]; // for each dimension, top - bot + 1;
        int nDim;                      // number of dimensions

        struct Limits
        {
            int bot[ SLICE_EQ_MAX_DIMS ];
            int top[ SLICE_EQ_MAX_DIMS ];
        };
        Limits limits;                 // current range limits

        struct InsNode
        {
            int index;
            double coeff1;
            double coeff2;

            InsNode( int index, double coeff1, double coeff2 ) :
                index( index ), coeff1( coeff1 ), coeff2( coeff2 ) {}
        };
        typedef refCountPtr< InsNode > InsNodeSP;
        vector< InsNodeSP > insNodes[ SLICE_EQ_MAX_DIMS ];

        void insertNode( int dim, const InsNodeSP & insNode );

    private:
        explicit Range(
            bool transposed = false,    // element matrix is transposed
            int bot1 = 0, int top1 = 0, // 1st dimension range
            int bot2 = 0, int top2 = 0, // 2st dimension range
            int bot3 = 0, int top3 = 0  // 3st dimension range
            );
    };

    virtual const string typeName() const { return "TreeSliceEQ"; }

    // creates new slice
    static TreeSliceEQSP create( const Range & range, int dimBits = 0 )
    {
        return TreeSliceEQSP( new TreeSliceEQ( range, dimBits ) );
    }
    
    // clones existing slice
    virtual TreeSliceSP clone( bool copyValues = true ) const
    {
        return TreeSliceSP( new TreeSliceEQ( *this, copyValues ) );
    }

    ~TreeSliceEQ() { ::free( values ); }

    template< typename T >
    TreeSliceEQ & operator =( const SliceMarker< T > & exprMarker )
    {
        eval( static_cast< const T & >( exprMarker ) );
        return *this;
    }
    TreeSliceEQ & operator =( const TreeSliceEQ & slice )
    {
        eval( slice );
        return *this;
    }
    virtual TreeSlice& operator =( double v )
    {
        eval( DoubleOperand( v ) );
        return *this;
    }

    int getDim() const { return subDim( 0 ) + subDim( 1 ) + subDim( 2 ); }

    virtual TreeSliceSP calcSmoothStep() const;

    const Range & getRange() const { return range; }
    int dimBits() const { return subDimBits; }
    int subDim( int dim ) const { return ( subDimBits >> dim ) & 1; }

    // expands slice if needed
    bool expand( int dimBits, bool keepValues = true )
    {
        if( ( subDimBits & dimBits ) == dimBits )
            return false;

        TreeSliceEQ slice( range, subDimBits | dimBits );
        if( keepValues )
            slice.eval( *this );

        // assume slice
        subDimBits = slice.subDimBits;

        ::memcpy( subSize, slice.subSize, sizeof( subSize ) );
        ::memcpy( insNodeCount, slice.insNodeCount, sizeof( insNodeCount ) );

        allocCount = slice.allocCount;
        valueCount = slice.valueCount;

        // swap values pointers so they get destroyed appropriately
        swapT( values, slice.values );

        offsetValues = slice.offsetValues;

        return true;
    }

    // 1-dimensional access
    int offset( int index1 ) const
    {
#ifdef DEBUG
        if( range.nDim != 1 )
            throw ModelException("TreeSliceEQ::operator( 1D )", "slice is not 1D");
        if( subDim( 0 ) && ( index1 < range.limits.bot[ 0 ] || range.limits.top[ 0 ] < index1 ) )
            throw ModelException("TreeSliceEQ::operator( 1D )", "index1 out of range");
#endif

        return subSize[ 0 ] * ( index1 - range.bot[ 0 ] );
    }
    double operator ()( int index1 ) const
    {
        return values[ offset( index1 ) ];
    }
    double & operator ()( int index1 )
    {
        return values[ offset( index1 ) ];
    }

    // 2-dimensional access
    int offset( int index1, int index2 ) const
    {
#ifdef DEBUG
        if( range.nDim != 2 )
            throw ModelException("TreeSliceEQ::operator( 2D )", "slice is not 2D");
        if( subDim( 0 ) && ( index1 < range.limits.bot[ 0 ] || range.limits.top[ 0 ] < index1 ) )
            throw ModelException("TreeSliceEQ::operator( 2D )", "index1 out of range");
        if( subDim( 1 ) && ( index2 < range.limits.bot[ 1 ] || range.limits.top[ 1 ] < index2 ) )
            throw ModelException("TreeSliceEQ::operator( 2D )", "index2 out of range");
#endif

        return subSize[ 0 ] * ( index1 - range.bot[ 0 ] )
             + subSize[ 1 ] * ( index2 - range.bot[ 1 ] );
    }
    double operator ()( int index1, int index2 ) const
    {
        return values[ offset( index1, index2 ) ];
    }
    double & operator ()( int index1, int index2 )
    {
        return values[ offset( index1, index2 ) ];
    }
    
    // 3-dimensional access
    int offset( int index1, int index2, int index3 ) const
    {
#ifdef DEBUG
        if( range.nDim != 3 )
            throw ModelException("TreeSliceEQ::operator( 3D )", "slice is not 3D");
        if( subDim( 0 ) && ( index1 < range.limits.bot[ 0 ] || range.limits.top[ 0 ] < index1 ) )
            throw ModelException("TreeSliceEQ::operator( 3D )", "index1 out of range");
        if( subDim( 1 ) && ( index2 < range.limits.bot[ 1 ] || range.limits.top[ 1 ] < index2 ) )
            throw ModelException("TreeSliceEQ::operator( 3D )", "index2 out of range");
        if( subDim( 2 ) && ( index3 < range.limits.bot[ 2 ] || range.limits.top[ 2 ] < index3 ) )
            throw ModelException("TreeSliceEQ::operator( 3D )", "index3 out of range");
#endif

        return subSize[ 0 ] * ( index1 - range.bot[ 0 ] )
             + subSize[ 1 ] * ( index2 - range.bot[ 1 ] )
             + subSize[ 2 ] * ( index3 - range.bot[ 2 ] );
    }
    double operator ()( int index1, int index2, int index3 ) const
    {
        return values[ offset( index1, index2, index3 ) ];
    }
    double & operator ()( int index1, int index2, int index3 )
    {
        return values[ offset( index1, index2, index3 ) ];
    }
    
    void setIter( int index1 ) const
    {
        iter = values + offset( index1 );
    }
    void setIter( int index1, int index2 ) const
    {
        iter = values + offset( index1, index2 );
    }
    void setIter( int index1, int index2, int index3 ) const
    {
        iter = values + offset( index1, index2, index3 );
    }
    double * nextIter( int dim ) const
    {
        return iter + subSize[ dim ];
    }

    // used by template slice operators
    template< typename ARG >
    void eval( const ARG & arg );

    // Do an operation on slices. Output slices come first
    // T must defile nbSlices, compute() and printDebug()
    template< typename T >
    void loopOnSlices( const T & oper, TreeSliceEQ ** slices, int nbOutput );

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const
    {
        if( getDim() <= 1 )
        {
            bot = subDim( 0 ) * range.limits.bot[ 0 ]
                + subDim( 1 ) * range.bot[ 1 ]  // instead of limist.bot[ 1 ]
                + subDim( 2 ) * range.bot[ 2 ]; // instead of limist.bot[ 2 ]
            top = subDim( 0 ) * range.limits.top[ 0 ]
                + subDim( 1 ) * ( range.top[ 1 ] )  // limist.top[ 1 ]
                + subDim( 2 ) * ( range.top[ 2 ] ); // limist.top[ 2 ]
        }
        else
        {
            bot = 0;
            top = valueCount - 1;
        }
    }
    virtual double * getValues() const
    {
        // insert nodes if any
        const_cast< TreeSliceEQ & >( *this ).insertNodes();

        // if slice is scalar, expand to all dimensions
        if( ! dimBits() )
            const_cast< TreeSliceEQ & >( *this ).expand( ( 1 << range.nDim ) - 1 );

        return offsetValues;
    }
    operator double *() const { return offsetValues; }

    void swapValues( TreeSliceEQ & slice )
    {
#ifdef DEBUG
        ASSERT( subDimBits == slice.subDimBits );
        ASSERT( valueCount == slice.valueCount );
        ASSERT( insNodeCount[ 0 ] == slice.insNodeCount[ 0 ] &&
                insNodeCount[ 1 ] == slice.insNodeCount[ 1 ] &&
                insNodeCount[ 2 ] == slice.insNodeCount[ 2 ] );
#endif

        swapT( allocCount, slice.allocCount );
        swapT( values, slice.values );
        swapT( offsetValues, slice.offsetValues );
    }

    void insertNodes( bool keepValues = true )
    {
        for( int d = 0; d < SLICE_EQ_MAX_DIMS; ++d )
            insertNodes( d, keepValues );
    }

    void setBotBound( int dim, int bound ) { botBound[ dim ] = bound; }
    void setTopBound( int dim, int bound ) { topBound[ dim ] = bound; }
    void setBotValue( int dim, double value ) { botValue[ dim ] = value; }
    void setTopValue( int dim, double value ) { topValue[ dim ] = value; }

    int getBotBound( int dim ) const { return botBound[ dim ]; }
    int getTopBound( int dim ) const { return topBound[ dim ]; }
    double getBotValue( int dim ) const { return botValue[ dim ]; }
    double getTopValue( int dim ) const { return topValue[ dim ]; }

private:
    const Range & range;                   // dimensions & running range limits
    int subDimBits;                        // 1/0 bits to specify used dimensions
    int subSize[ SLICE_EQ_MAX_DIMS ];      // dimension sizes
    int insNodeCount[ SLICE_EQ_MAX_DIMS ]; // number of nodes already inserted into the slice
    int allocCount;                        // number of values allocated
    int valueCount;                        // number of values
    double * values;                       // pointer to slice values
    double * offsetValues;                 // shifted pointer to slice values (for legacy purposes only)
    //!!! REVIEW THIS
    // solve boundaries
    int botBound[ SLICE_EQ_MAX_DIMS ];
    int topBound[ SLICE_EQ_MAX_DIMS ];
    // values at solve boundaries
    double botValue[ SLICE_EQ_MAX_DIMS ];
    double topValue[ SLICE_EQ_MAX_DIMS ];
    
    mutable LoopList< const TreeSliceEQ * > loopSlices;
    mutable LoopList< double * > loopIters;

    explicit TreeSliceEQ( const Range & range, int dimBits = 0 );
    TreeSliceEQ( const TreeSliceEQ & slice, bool copyValues = true );

    void insertNodes( int dim, bool keepValues = true );
};

namespace
{
// helper class to avoid additional dereferencing
template< typename Slice >
class SliceRef : public refCountPtr< Slice >
{
    typedef refCountPtr< Slice > SliceSP;
public:
    SliceRef( Slice * slice = 0 ) : SliceSP( slice ) {}
    SliceRef( const SliceSP & slice ) : SliceSP( slice ) {}
    virtual ~SliceRef() {}

    operator const Slice &() const { return SliceSP::operator *(); }
    operator Slice &() { return SliceSP::operator *(); }

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    void getCalcRange( int & bot, int & top ) const { SliceSP::operator *().getCalcRange( bot, top ); }
    double * getValues() const { return SliceSP::operator *().getValues(); }
    operator double *() const { return SliceSP::operator *(); }
};
}
typedef SliceRef< TreeSliceEQ > TreeSliceEQRef;

#define STATIC_POINTER_CAST boost::static_pointer_cast

namespace
{
template< typename ARG >
class EvalOper
{
    const ARG & arg;
    TreeSliceEQ & slice;

public:
    EvalOper( const ARG & arg, TreeSliceEQ & slice ) : arg( arg ), slice( slice ) {}

    static const int sliceCount = ARG::sliceCount + 1;
    void compute() const { *slice.iter = arg.calc(); }
};
}

template< typename ARG >
inline
void TreeSliceEQ::eval( const ARG & arg )
{
#ifdef DEBUG
    const char* (*print)(void*) = printExpr< ARG >;
    (void)print; // to remove compiler warnings;

    // watch "print((void*)arg)" and "*this->iter" in your debugger
#endif

    const TreeSliceEQ ** slices = loopSlices.reserve( arg.sliceCount + 1 );
    *slices = this;
    const TreeSliceEQ ** end = arg.listSlices( slices + 1 );

    ASSERT( end - slices == arg.sliceCount + 1 );
#ifdef DEBUG
    for( int n = 0; n < arg.sliceCount; ++n )
        ASSERT( slices[ n ] );      
#endif

    loopOnSlices( EvalOper< ARG >( arg, *this ), const_cast< TreeSliceEQ ** >( slices ), 1 );
}

template< typename T >
inline
void TreeSliceEQ::loopOnSlices( const T & oper, TreeSliceEQ ** slices, int nbOutput )
{
    const Range & range = (*slices)->range;

    // insert nodes if any
    for( int n = 0; n < oper.sliceCount; ++n )
        slices[ n ]->insertNodes();

    int dimBits = 0;
    for( int n = 0; n < oper.sliceCount; ++n )
        dimBits |= slices[ n ]->dimBits();

    // expand slices if needed
    for( int n = 0; n < nbOutput; ++n )
        slices[ n ]->expand( dimBits );

    double ** iters = loopIters.reserve( oper.sliceCount );
    switch( range.nDim )
    {
        case 1:
        {
            int subDim1 = ( dimBits >> 0 ) & 1;

            int bot1 = subDim1 * range.limits.bot[ 0 ];
            int top1 = subDim1 * range.limits.top[ 0 ];
            for( int n = 0; n < oper.sliceCount; ++n )
                slices[ n ]->setIter( bot1 );
            for( int i = bot1; i <= top1; ++i )
            {
                oper.compute();
                for( int n = 0; n < oper.sliceCount; ++n )
                    iters[ n ] = slices[ n ]->nextIter( 0 );
                for( int n = 0; n < oper.sliceCount; ++n )
                    slices[ n ]->iter = iters[ n ];
            }
            break;
        }
        case 2:
        {
            int subDim1 = ( dimBits >> 0 ) & 1;
            int subDim2 = ( dimBits >> 1 ) & 1;

            int bot1 = subDim1 * range.limits.bot[ 0 ];
            int top1 = subDim1 * range.limits.top[ 0 ];
            for( int i = bot1; i <= top1; ++i )
            {
                int bot2 = subDim2 * range.limits.bot[ 1 ];
                int top2 = subDim2 * range.limits.top[ 1 ];
                for( int n = 0; n < oper.sliceCount; ++n )
                    slices[ n ]->setIter( i, bot2 );
                for( int j = bot2; j <= top2; ++j )
                {
                    oper.compute();
                    for( int n = 0; n < oper.sliceCount; ++n )
                        iters[ n ] = slices[ n ]->nextIter( 1 );
                    for( int n = 0; n < oper.sliceCount; ++n )
                        slices[ n ]->iter = iters[ n ];
                }
            }
            break;
        }
        case 3:
        {
            int subDim1 = ( dimBits >> 0 ) & 1;
            int subDim2 = ( dimBits >> 1 ) & 1;
            int subDim3 = ( dimBits >> 2 ) & 1;

            int bot1 = subDim1 * range.limits.bot[ 0 ];
            int top1 = subDim1 * range.limits.top[ 0 ];
            for( int i = bot1; i <= top1; ++i )
            {
                int bot2 = subDim2 * range.limits.bot[ 1 ];
                int top2 = subDim2 * range.limits.top[ 1 ];
                for( int j = bot2; j <= top2; ++j )
                {
                    int bot3 = subDim3 * range.limits.bot[ 2 ];
                    int top3 = subDim3 * range.limits.top[ 2 ];
                    for( int n = 0; n < oper.sliceCount; ++n )
                        slices[ n ]->setIter( i, j, bot3 );
                    for( int k = bot3; k <= top3; ++k )
                    {
                        oper.compute();
                        for( int n = 0; n < oper.sliceCount; ++n )
                            iters[ n ] = slices[ n ]->nextIter( 2 );
                        for( int n = 0; n < oper.sliceCount; ++n )
                            slices[ n ]->iter = iters[ n ];
                    }
                }
            }
            break;
        }
        default:
            throw ModelException( "TreeSliceEQ::loopOnSlices", "unsupported slice dimensionality" );
    }
}

DRLIB_END_NAMESPACE

#endif
