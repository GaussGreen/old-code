//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceGeneral.hpp
//
//   Description : General tree slice implementation.
//
//   Date        : Apr 11, 2006
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_GENERAL_HPP
#define TREE_SLICE_GENERAL_HPP

DRLIB_BEGIN_NAMESPACE

class TreeSliceGeneral;
typedef refCountPtr< TreeSliceGeneral > TreeSliceGeneralSP;
class TreeSliceGeneral : public TreeSlice
{
public:
    class Range;
    typedef refCountPtr< Range > RangeSP;
    class Range
    {
    public:
        // allocates new 1D slice range
        static RangeSP create(
            int count,
            int bot1, int top1,
            bool transposed = false )
        {
            return RangeSP(
                new Range( count, transposed, bot1, top1 ) );
        }

        // allocates new 2D slice range
        static RangeSP create(
            int count,
            int bot1, int top1, int bot2, int top2,
            bool transposed = false )
        {
            return RangeSP(
                new Range( count, transposed, bot1, top1, bot2, top2 ) );
        }

        // allocates new 3D slice range
        static RangeSP create(
            int count,
            int bot1, int top1, int bot2, int top2, int bot3, int top3,
            bool transposed = false )
        {
            return RangeSP(
                new Range( count, transposed, bot1, top1, bot2, top2, bot3, top3 ) );
        }

        const TreeSliceBasic::Range & operator []( int index ) const { return *ranges[ index ]; }
        TreeSliceBasic::Range & operator []( int index ) { return *ranges[ index ]; }

        const TreeSliceBasic::Range & operator *() const { return operator[]( rangeIndex ); }
        TreeSliceBasic::Range & operator *() { return operator[]( rangeIndex ); }

        operator const TreeSliceBasic::Range &() const { return operator *(); }
        operator TreeSliceBasic::Range &() { return operator *(); }

        const TreeSliceBasic::Range * operator ->() const { return &operator *(); }
        TreeSliceBasic::Range * operator ->() { return &operator *(); }

        int count() const { return ranges.size(); }

        void setIndex( int index ) { rangeIndex = index; }
        int getIndex() const { return rangeIndex; }

    private:
        vector< TreeSliceBasic::RangeSP > ranges; // vector of ranges
        int rangeIndex;                           // current range index

        explicit Range(
            int count,                  // number of ranges
            bool transposed = false,    // element matrix is transposed
            int bot1 = 0, int top1 = 0, // 1st dimension range
            int bot2 = 0, int top2 = 0, // 2st dimension range
            int bot3 = 0, int top3 = 0  // 3st dimension range
            ) :
            ranges( count ),
            rangeIndex( 0 )
        {
            for( int i = 0; i < count; ++i )
            {
                ranges[ i ] = TreeSliceBasic::Range::create(
                    bot1, top1, bot2, top2, bot3, top3, transposed );
            }
        }
    };

public:
    // creates new slice
    static TreeSliceGeneralSP create( const Range & range, int dimBits = 0 )
    {
        return TreeSliceGeneralSP( new TreeSliceGeneral( range, dimBits ) );
    }
    
    virtual const string typeName() const { return "TreeSliceGeneral"; }
    // clones existing slice
    virtual TreeSliceSP clone( bool copyValues = true ) const
    {
        return TreeSliceSP( new TreeSliceGeneral( *this, copyValues ) );
    }

    const TreeSliceBasic & operator []( int index ) const { return *slices[ index ]; }
    TreeSliceBasic & operator []( int index ) { return *slices[ index ]; }

    const TreeSliceBasic & operator *() const { return operator []( range.getIndex() ); }
    TreeSliceBasic & operator *() { return operator []( range.getIndex() ); }

    operator const TreeSliceBasic &() const { return operator *(); }
    operator TreeSliceBasic &() { return operator *(); }

    const TreeSliceBasic * operator ->() const { return &operator *(); }
    TreeSliceBasic * operator ->() { return &operator *(); }

    int count() const { return range.count(); }

    const Range & getRange() const { return range; }
    int dimBits() const { return operator *().dimBits(); }

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const { operator *().getCalcRange( bot, top ); }
    virtual double * getValues() const { return operator *().getValues(); }

    void expand( int dimBits, bool keepValues = true )
    {
        int count = slices.size();
        for( int i = 0; i < count; ++i )
            slices[ i ]->expand( dimBits, keepValues );
    }

    virtual int getDim() const { return operator *().getDim(); }


	virtual bool isZero() const { 
		return ( getDim() &&  Maths::isZero(*getValues()) ); 
	}

    // 1-dimensional access
    double operator ()( int index1 ) const
    {
        return operator *().operator ()( index1 );
    }
    double & operator ()( int index1 )
    {
        return operator *().operator ()( index1 );
    }

    // 2-dimensional access
    double operator ()( int index1, int index2 ) const
    {
        return operator *().operator ()( index1, index2 );
    }
    double & operator ()( int index1, int index2 )
    {
        return operator *().operator ()( index1, index2 );
    }

    // 3-dimensional access
    double operator ()( int index1, int index2, int index3 ) const
    {
        return operator *().operator ()( index1, index2, index3 );
    }
    double & operator ()( int index1, int index2, int index3 )
    {
        return operator *().operator ()( index1, index2, index3 );
    }

    // used by template slice operators
    template< typename ARG >
    void eval( const ARG & arg )
    {
        // get a list of the slices used in the expression
        const TreeSliceGeneral ** slices = loopSlices.reserve( arg.sliceCount + 1 );
        const TreeSliceGeneral ** end = arg.listSlices( slices );
        slices[ arg.sliceCount ] = this;

        ASSERT( end - slices == arg.sliceCount );
#ifdef DEBUG
        for( int n = 0; n < arg.sliceCount; ++n )
            ASSERT( slices[ n ] );
#endif

        ExposedSliceOverride * override = loopOverride.reserve( arg.sliceCount + 1 );
        for( int n = 0; n <= arg.sliceCount; ++n )
        {
            new ( override + n ) ExposedSliceOverride;
            override[ n ].set( slices[ n ], slices[ n ]->operator->() );
        }

        operator *().eval( arg );

        for( int n = arg.sliceCount; 0 <= n; --n )
            override[ n ].~ExposedSliceOverride();
    }

    template< typename T >
    TreeSliceGeneral & operator =( const SliceMarker< T > & exprMarker )
    {
        eval( static_cast< const T & >( exprMarker ) );
        return *this;
    }
    TreeSliceGeneral & operator =( const TreeSliceGeneral & slice )
    {
        eval( slice );
        return *this;
    }
    virtual TreeSlice& operator =( double v )
    {
        eval( DoubleOperand( v ) );
        return *this;
    }

    // used by template slice operators
    template< typename ARG >
    TreeSliceGeneral & operator +=( const ARG & arg ) { return *this = *this + arg; }
    template< typename ARG >
    TreeSliceGeneral & operator -=( const ARG & arg ) { return *this = *this - arg; }
    template< typename ARG >
    TreeSliceGeneral & operator *=( const ARG & arg ) { return *this = *this * arg; }
    template< typename ARG >
    TreeSliceGeneral & operator /=( const ARG & arg ) { return *this = *this / arg; }

private:
    const Range & range;               // dimensions & running range limits
    vector< TreeSliceBasicSP > slices; // vector of slices

    mutable LoopList< const TreeSliceGeneral * > loopSlices;
    mutable LoopList< ExposedSliceOverride > loopOverride;

    explicit TreeSliceGeneral( const Range & range, int dimBits = 0 ) :
        range( range )
    {
        int count = this->count();
        slices.resize( count );
        for( int i = 0; i < count; ++i )
            slices[ i ] = TreeSliceBasic::create( range[ i ], dimBits );
    }

    TreeSliceGeneral( const TreeSliceGeneral & slice, bool copyValues = true ) :
        range( slice.range )
    {
        int count = this->count();
        slices.resize( count );
        for( int i = 0; i < count; ++i )
            slices[ i ] = DYNAMIC_POINTER_CAST<TreeSliceBasic>( slice.slices[ i ]->clone( copyValues ) );
    }
};

class TreeSliceGeneralCont;
typedef refCountPtr< TreeSliceGeneralCont > TreeSliceGeneralContSP;
class TreeSliceGeneralCont
{
public:
    // allocates new slice container
    static TreeSliceGeneralContSP create( const TreeSliceGeneral::Range & range, int size, int dimBits = 0 )
    {
        return TreeSliceGeneralContSP( new TreeSliceGeneralCont( range, size, dimBits ) );
    }

    operator double const * const *() const { return &values[ range.getIndex() ][ 0 ]; }
    operator double * const *() { return &values[ range.getIndex() ][ 0 ]; }

    double const * const * operator []( int rangeIndex ) const { return &values[ rangeIndex ][ 0 ]; }
    double * const * operator []( int rangeIndex ) { return &values[ rangeIndex ][ 0 ]; }

    operator const vector< TreeSliceSP > &() const { return slices; }
    const vector< TreeSliceSP > & operator *() const { return slices; }

    int count() const { return range.count(); }

private:
    TreeSliceGeneralCont( const TreeSliceGeneral::Range & range, int size, int dimBits = 0 ) :
        range( range )
    {
        slices.resize( size );
        for( int i = 0; i < size; ++i )
            slices[ i ] = TreeSliceGeneral::create( range, dimBits );

        int count = this->count();
        values.resize( count );
        for( int j = 0; j < count; ++j )
        {
            values[ j ].resize( size );
            for( int i = 0; i < size; ++i )
            {
                values[ j ][ i ] =
                    static_cast< TreeSliceGeneral & >( *slices[ i ] )[ j ].getValues();
            }
        }
    }

    const TreeSliceGeneral::Range & range;
    vector< TreeSliceSP > slices;
    vector< vector< double * > > values;
};

DRLIB_END_NAMESPACE

#endif
