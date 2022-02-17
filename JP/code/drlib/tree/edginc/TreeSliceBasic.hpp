//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceBasic.hpp
//
//   Description : Basic tree slice implementation.
//
//   Date        : Sep 8, 2005
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_BASIC_HPP
#define TREE_SLICE_BASIC_HPP

#include "edginc/Format.hpp"

#define SLICE_BASIC_MAX_DIMS 3

DRLIB_BEGIN_NAMESPACE

// holds slice values
// supports up to SLICE_BASIC_MAX_DIMS dimensions (somewhat hard-coded)
class TreeSliceBasic;
typedef refCountPtr< TreeSliceBasic > TreeSliceBasicSP;

class TREE_DLL TreeSliceBasic : public TreeSlice
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

        ~Range();

        int getSize( int dim = -1 ) const;

        struct Limits
        {
            int    bot1,    top1;
            int  * bot2,  * top2;
            int ** bot3, ** top3;
        };

        bool transposed;
        int bot[ SLICE_BASIC_MAX_DIMS ];  // for each dimension, starting index from outside
        int top[ SLICE_BASIC_MAX_DIMS ];  // for each dimension, ending index from outside
        int size[ SLICE_BASIC_MAX_DIMS ]; // for each dimension, top - bot + 1;
        int nDim;                         // number of dimensions
        Limits limits;                    // current range limits
        // example:
        // limits.bot3[index1][index2] is the bottom index for 3rd dimension given indexes in 1st and 2nd dimensions

    private:
        explicit Range(
            bool transposed = false,    // element matrix is transposed
            int bot1 = 0, int top1 = 0, // 1st dimension range
            int bot2 = 0, int top2 = 0, // 2st dimension range
            int bot3 = 0, int top3 = 0  // 3st dimension range
            );
    };

    // creates new slice
    static TreeSliceBasicSP create( const Range & range, int dimBits = 0 )
    {
        return TreeSliceBasicSP( new TreeSliceBasic( range, dimBits ) );
    }
    
    virtual const string typeName() const { return "TreeSliceBasic"; }
    // clones existing slice
    virtual TreeSliceSP clone( bool copyValues = true ) const
    {
        return TreeSliceSP( new TreeSliceBasic( *this, copyValues ) );
    }

    ~TreeSliceBasic() { delete [] values; }

    template< typename T >
    TreeSliceBasic & operator =( const SliceMarker< T > & exprMarker )
    {
        eval( static_cast< const T & >( exprMarker ) );
        return *this;
    }
    TreeSliceBasic & operator =( const TreeSliceBasic & slice )
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

    void expand( int dimBits, bool keepValues = true );

    // 1-dimensional access
    int offset( int index1 ) const
    {
#ifdef DEBUG
        if( range.nDim != 1 )
            throw ModelException("TreeSliceBasic::operator( 1D )", "slice is not 1D");
        if( subDim( 0 ) && ( index1 < range.limits.bot1 || range.limits.top1 < index1 ) )
            throw ModelException("TreeSliceBasic::operator( 1D )", "index1 out of range");
#endif

        return subSize[0] * ( index1 - range.bot[0] );
    }

        // 1-dimensional access
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
            throw ModelException("TreeSliceBasic::operator( 2D )", "slice is not 2D");
        if( subDim( 0 ) && ( index1 < range.limits.bot1 || range.limits.top1 < index1 ) )
            throw ModelException("TreeSliceBasic::operator( 2D )", "index1 out of range");
        if( subDim( 1 ) && ( index2 < range.limits.bot2[ index1 ] || range.limits.top2[ index1 ] < index2 ) )
            throw ModelException("TreeSliceBasic::operator( 2D )", "index2 out of range");
#endif

        return subSize[0] * ( index1 - range.bot[0] )
             + subSize[1] * ( index2 - range.bot[1] );
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
            throw ModelException("TreeSliceBasic::operator( 3D )", "slice is not 3D");
        if( subDim( 0 ) && ( index1 < range.limits.bot1 || range.limits.top1 < index1 ) )
            throw ModelException("TreeSliceBasic::operator( 3D )", "index1 out of range");
        if( subDim( 1 ) && ( index2 < range.limits.bot2[ index1 ] || range.limits.top2[ index1 ] < index2 ) )
            throw ModelException("TreeSliceBasic::operator( 3D )", "index2 out of range");
        if( subDim( 2 ) && ( index3 < range.limits.bot3[ index1 ][ index2 ] || range.limits.top3[ index1 ][ index2 ] < index3 ) )
            throw ModelException("TreeSliceBasic::operator( 3D )", "index3 out of range");
#endif

        return subSize[0] * ( index1 - range.bot[0] )
             + subSize[1] * ( index2 - range.bot[1] )
             + subSize[2] * ( index3 - range.bot[2] );
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
    void eval( const ARG & arg )
    {
#ifdef DEBUG
        const char* (*print)(void*) = printExpr< ARG >;
        (void)print; // to remove compiler warnings;

        // watch "print((void*)arg)" and "*this->iter" in your debugger
#endif

        // get a list of the slices used in the expression
        const TreeSliceBasic ** slices = loopSlices.reserve( arg.sliceCount + 1 );
        const TreeSliceBasic ** end = arg.listSlices( slices );
        slices[ arg.sliceCount ] = this;

        ASSERT( end - slices == arg.sliceCount );
#ifdef DEBUG
        for( int n = 0; n < arg.sliceCount; ++n )
            ASSERT( slices[ n ] );
#endif

        int dimBits = 0;
        for( int n = 0; n < arg.sliceCount; ++n )
            dimBits |= slices[ n ]->dimBits();

        if( ( subDimBits & dimBits ) != dimBits ) 
        {
            // expand slice
            TreeSliceBasic slice( range, subDimBits | dimBits );
            slice.eval( arg );
            assume( slice );
            return;
        }

        double ** iters = loopIters.reserve( arg.sliceCount + 1 );
        switch( range.nDim )
        {
            case 1:
            {
                int bot1 = subDim( 0 ) * range.limits.bot1;
                int top1 = subDim( 0 ) * range.limits.top1;
                for( int n = 0; n <= arg.sliceCount; ++n )
                    slices[ n ]->setIter( bot1 );
                for( int i = bot1; i <= top1; ++i )
                {
                    *iter = arg.calc();
                    for( int n = 0; n <= arg.sliceCount; ++n )
                        iters[ n ] = slices[ n ]->nextIter( 0 );
                    for( int n = 0; n <= arg.sliceCount; ++n )
                        slices[ n ]->iter = iters[ n ];
                }
                break;
            }
            case 2:
            {
                int bot1 = subDim( 0 ) * range.limits.bot1;
                int top1 = subDim( 0 ) * range.limits.top1;
                for( int i = bot1; i <= top1; ++i )
                {
                    int bot2 = subDim( 1 ) * range.limits.bot2[ i ];
                    int top2 = subDim( 1 ) * range.limits.top2[ i ];
                    for( int n = 0; n <= arg.sliceCount; ++n )
                        slices[ n ]->setIter( i, bot2 );
                    for( int j = bot2; j <= top2; ++j )
                    {
                        *iter = arg.calc();
                        for( int n = 0; n <= arg.sliceCount; ++n )
                            iters[ n ] = slices[ n ]->nextIter( 1 );
                        for( int n = 0; n <= arg.sliceCount; ++n )
                            slices[ n ]->iter = iters[ n ];
                    }
                }
                break;
            }
            case 3:
            {
                int bot1 = subDim( 0 ) * range.limits.bot1;
                int top1 = subDim( 0 ) * range.limits.top1;
                for( int i = bot1; i <= top1; ++i )
                {
                    int bot2 = subDim( 1 ) * range.limits.bot2[ i ];
                    int top2 = subDim( 1 ) * range.limits.top2[ i ];
                    for( int j = bot2; j <= top2; ++j )
                    {
                        int bot3 = subDim( 2 ) * range.limits.bot3[ i ][ j ];
                        int top3 = subDim( 2 ) * range.limits.top3[ i ][ j ];
                        for( int n = 0; n <= arg.sliceCount; ++n )
                            slices[ n ]->setIter( i, j, bot3 );
                        for( int k = bot3; k <= top3; ++k )
                        {
                            *iter = arg.calc();
                            for( int n = 0; n <= arg.sliceCount; ++n )
                                iters[ n ] = slices[ n ]->nextIter( 2 );
                            for( int n = 0; n <= arg.sliceCount; ++n )
                                slices[ n ]->iter = iters[ n ];
                        }
                    }
                }
                break;
            }
            default:
                throw ModelException( "TreeSliceBasic::assign", "unsupported slice dimensionality" );
        }
    }

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const
    {
        if( getDim() <= 1 )
        {
            bot = subDim( 0 ) * range.limits.bot1
                + subDim( 1 ) * range.bot[1]  // instead of min(bot2[i])
                + subDim( 2 ) * range.bot[2]; // instead of min(bot3[i][j])
            top = subDim( 0 ) * range.limits.top1
                + subDim( 1 ) * ( range.top[1] )  // instead of max(top2[i])
                + subDim( 2 ) * ( range.top[2] ); // instead of max(top3[i][j])
        }
        else
        {
            bot = 0;
            top = valueCount - 1;
        }
    }
    virtual double * getValues() const { return offsetValues; }

    operator double *() const { return offsetValues; }

private:
    const Range & range;                 // dimensions & running range limits
    int subDimBits;                      // 1/0 bits to specify used dimensions
    int subSize[ SLICE_BASIC_MAX_DIMS ]; // dimension sizes
    int valueCount;                      // number of values allocated
    double * values;                     // pointer to slice values
    double * offsetValues;               // shifted pointer to slice values (for legacy purposes only)
    
    mutable LoopList< const TreeSliceBasic * > loopSlices;
    mutable LoopList< double * > loopIters;

    explicit TreeSliceBasic( const Range & range, int dimBits = 0 );
    TreeSliceBasic( const TreeSliceBasic & slice, bool copyValues = true );

    void assume( TreeSliceBasic & slice );
};

DRLIB_END_NAMESPACE

#endif
