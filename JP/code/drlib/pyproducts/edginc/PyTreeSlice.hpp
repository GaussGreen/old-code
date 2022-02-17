#ifndef QLIB_PY_TREE_SLICE_HPP
#define QLIB_PY_TREE_SLICE_HPP

DRLIB_BEGIN_NAMESPACE

class SliceOperandPy
{
protected:
    //----------------------------------------------------------------------------
    // slice marker
    //----------------------------------------------------------------------------
    class SliceMarkerPy : public SliceMarker< SliceMarkerPy >
    {
    protected:
        SliceMarkerPy( int sliceCount ) : sliceCount( sliceCount ) {}

    public:
        const int sliceCount;

        virtual const TreeSlice** listSlices(const TreeSlice** list) const = 0;
        virtual double calc( const double * = 0 ) const
        {
            throw ModelException( "SliceMarkerPy::calc( double )", "Invalid slice expression" );
        }
        virtual bool calc( const bool * = 0 ) const
        {
            throw ModelException( "SliceMarkerPy::calc( bool )", "Invalid slice expression" );
        }

        void printDebug(char *s) const {}

        virtual ~SliceMarkerPy() {}
    };
    typedef refCountPtr< SliceMarkerPy > SliceMarkerPySP;

    //----------------------------------------------------------------------------
    // simple operand
    //----------------------------------------------------------------------------
    template< typename T >
    class SimpleOperandPy : public SliceMarkerPy
    {
        T v;

    public:
        SimpleOperandPy( T v ) :
            SliceMarkerPy( 0 ),
            v( v ) {}

        const TreeSlice** listSlices(const TreeSlice** list) const { return list; }
        T calc( const T * = 0 ) const { return v; }
    };
    template< typename T >
    static SliceMarkerPySP create( T v )
    {
        return SliceMarkerPySP( new SimpleOperandPy< T >( v ) );
    }

    //----------------------------------------------------------------------------
    // slice operand
    //----------------------------------------------------------------------------
    class TreeSliceOperandPy : public SliceMarkerPy
    {
        const TreeSlice & slice;

    public:
        TreeSliceOperandPy( const TreeSlice & slice ) :
            SliceMarkerPy( 1 ),
            slice( slice ) {}

        const TreeSlice** listSlices(const TreeSlice** list) const { return slice.listSlices( list ); }
        double calc( const double * = 0 ) const { return slice.calc(); }
    };

    //----------------------------------------------------------------------------
    // unary slice expression
    //----------------------------------------------------------------------------
    template< typename T, typename ARG >
    class SliceUnaryExprPy : public SliceMarkerPy
    {
        typedef T (*Oper)( ARG );
        Oper oper;
        SliceMarkerPySP arg;

    public:
        SliceUnaryExprPy( Oper oper, const SliceMarkerPySP & arg ) :
            SliceMarkerPy( arg->sliceCount ),
            oper( oper ), arg( arg ) {}

        const TreeSlice** listSlices(const TreeSlice** list) const { return arg->listSlices( list ); }
        T calc( const T * = 0 ) const { return (*oper)( arg->calc( (const ARG *)0 ) ); }
    };
    template< typename T, typename ARG >
    static SliceMarkerPySP create( T (*oper)( ARG ), const SliceMarkerPySP & arg )
    {
        return SliceMarkerPySP( new SliceUnaryExprPy< T, ARG >( oper, arg ) );
    }

    //----------------------------------------------------------------------------
    // binary slice expression
    //----------------------------------------------------------------------------
    template< typename T, typename L, typename R >
    class SliceBinaryExprPy : public SliceMarkerPy
    {
        typedef T (*Oper)( L, R );
        Oper oper;
        SliceMarkerPySP l;
        SliceMarkerPySP r;

    public:
        SliceBinaryExprPy( Oper oper, const SliceMarkerPySP & l, const SliceMarkerPySP & r ) :
            SliceMarkerPy( l->sliceCount + r->sliceCount ),
            oper( oper ), l( l ), r( r ) {}

        const TreeSlice** listSlices(const TreeSlice** list) const { return r->listSlices( l->listSlices( list ) ); }
        T calc( const T * = 0 ) const { return (*oper)( l->calc( (const L *)0 ), r->calc( (const R *)0 ) ); }
    };
    template< typename T, typename L, typename R >
    static SliceMarkerPySP create( T (*oper)( L, R ), const SliceMarkerPySP & l, const SliceMarkerPySP & r )
    {
        return SliceMarkerPySP( new SliceBinaryExprPy< T, L, R >( oper, l, r ) );
    }

    //----------------------------------------------------------------------------
    // 3 argument slice expression
    //----------------------------------------------------------------------------
    template< typename T, typename S0, typename S1, typename S2 >
    class Slice3ExprPy : public SliceMarkerPy
    {
        typedef T (*Oper)( S0, S1, S2 );
        Oper oper;
        SliceMarkerPySP s0;
        SliceMarkerPySP s1;
        SliceMarkerPySP s2;

    public:
        Slice3ExprPy(
            Oper oper,
            const SliceMarkerPySP & s0,
            const SliceMarkerPySP & s1,
            const SliceMarkerPySP & s2 )
            :
            SliceMarkerPy( s0->sliceCount + s1->sliceCount + s2->sliceCount ),
            oper( oper ), s0( s0 ), s1( s1 ), s2( s2 ) {}

        const TreeSlice** listSlices(const TreeSlice** list) const
        {
            return
                s2->listSlices(
                    s1->listSlices(
                        s0->listSlices( list ) ) );
        }
        T calc( const T * = 0 ) const
        {
            return (*oper)( s0->calc( (const S0 *)0 ), s1->calc( (const S1 *)0 ), s2->calc( (const S2 *)0 ) );
        }
    };
    template< typename T, typename S0, typename S1, typename S2 >
    static SliceMarkerPySP create(
        T (*oper)( S0, S1, S2 ),
        const SliceMarkerPySP & s0,
        const SliceMarkerPySP & s1,
        const SliceMarkerPySP & s2 )
    {
        return SliceMarkerPySP( new Slice3ExprPy< T, S0, S1, S2 >( oper, s0, s1, s2 ) );
    }

public:
    SliceMarkerPySP obj;

    SliceOperandPy( const SliceMarkerPySP & obj ) : obj( obj ) {}
    virtual ~SliceOperandPy() {}

    static SliceOperandPy __smin( const SliceOperandPy & l, const SliceOperandPy & r )
    {
        return create( &oper_min::apply, l.obj, r.obj );
    }
    static SliceOperandPy __sminL( const SliceOperandPy & l, double r )
    {
        return create( &oper_min::apply, l.obj, create( r ) );
    }
    static SliceOperandPy __sminR( double l, const SliceOperandPy & r )
    {
        return create( &oper_min::apply, create( l ), r.obj );
    }

    static SliceOperandPy __smax( const SliceOperandPy & l, const SliceOperandPy & r )
    {
        return create( &oper_max::apply, l.obj, r.obj );
    }
    static SliceOperandPy __smaxL( const SliceOperandPy & l, double r )
    {
        return create( &oper_max::apply, l.obj, create( r ) );
    }
    static SliceOperandPy __smaxR( double l, const SliceOperandPy & r )
    {
        return create( &oper_max::apply, create( l ), r.obj );
    }

    static SliceOperandPy __log( const SliceOperandPy & l )
    {
        return create( &oper_log::apply, l.obj );
    }

    static SliceOperandPy __exp( const SliceOperandPy & l )
    {
        return create( &oper_exp::apply, l.obj );
    }

    static SliceOperandPy __sqrt( const SliceOperandPy & l )
    {
        return create( &oper_sqrt::apply, l.obj );
    }

    static SliceOperandPy __cond( const SliceOperandPy & c, const SliceOperandPy & a, const SliceOperandPy & b )
    {
        return create( &oper_cond::apply, c.obj, a.obj, b.obj );
    }
    static SliceOperandPy __condL( const SliceOperandPy & c, const SliceOperandPy & a, double b )
    {
        return create( &oper_cond::apply, c.obj, a.obj, create( b ) );
    }
    static SliceOperandPy __condR( const SliceOperandPy & c, double a, const SliceOperandPy & b )
    {
        return create( &oper_cond::apply, c.obj, create( a ), b.obj );
    }
    static SliceOperandPy __condB( bool c, const SliceOperandPy & a, const SliceOperandPy & b )
    {
        return create( &oper_cond::apply, create( c ), a.obj, b.obj );
    }
    static SliceOperandPy __condBL( bool c, const SliceOperandPy & a, double b )
    {
        return create( &oper_cond::apply, create( c ), a.obj, create( b ) );
    }
    static SliceOperandPy __condBR( bool c, double a, const SliceOperandPy & b )
    {
        return create( &oper_cond::apply, create( c ), create( a ), b.obj );
    }

    SliceOperandPy __add( const SliceOperandPy & r )
    {
        return create( &oper_add::apply, obj, r.obj );
    }
    SliceOperandPy __addL( double r )
    {
        return create( &oper_add::apply, obj, create( r ) );
    }
    SliceOperandPy __addR( double l )
    {
        return create( &oper_add::apply, create( l ), obj );
    }

    SliceOperandPy __sub( const SliceOperandPy & r )
    {
        return create( &oper_sub::apply, obj, r.obj );
    }
    SliceOperandPy __subL( double r )
    {
        return create( &oper_sub::apply, obj, create( r ) );
    }
    SliceOperandPy __subR( double l )
    {
        return create( &oper_sub::apply, create( l ), obj );
    }

    SliceOperandPy __mul( const SliceOperandPy & r )
    {
        return create( &oper_mul::apply, obj, r.obj );
    }
    SliceOperandPy __mulL( double r )
    {
        return create( oper_mul::apply, obj, create( r ) );
    }
    SliceOperandPy __mulR( double l )
    {
        return create( &oper_mul::apply, create( l ), obj );
    }

    SliceOperandPy __div( const SliceOperandPy & r )
    {
        return create( &oper_div::apply, obj, r.obj );
    }
    SliceOperandPy __divL( double r )
    {
        return create( &oper_div::apply, obj, create( r ) );
    }
    SliceOperandPy __divR( double l )
    {
        return create( &oper_div::apply, create( l ), obj );
    }

    SliceOperandPy __pow( const SliceOperandPy & r )
    {
        return create( &oper_pow::apply, obj, r.obj );
    }
    SliceOperandPy __powL( double r )
    {
        return create( &oper_pow::apply, obj, create( r ) );
    }
    SliceOperandPy __powR( double l )
    {
        return create( &oper_pow::apply, create( l ), obj );
    }

    SliceOperandPy __pos()
    {
        return create( &oper_pos::apply, obj );
    }
    SliceOperandPy __neg()
    {
        return create( &oper_neg::apply, obj );
    }
    SliceOperandPy __abs()
    {
        return create( &oper_abs::apply, obj );
    }

    SliceOperandPy __lt( const SliceOperandPy & r )
    {
        return create( &oper_lt::apply, obj, r.obj );
    }
    SliceOperandPy __ltL( double r )
    {
        return create( &oper_lt::apply, obj, create( r ) );
    }
    SliceOperandPy __ltR( double l )
    {
        return create( &oper_lt::apply, create( l ), obj );
    }

    SliceOperandPy __le( const SliceOperandPy & r )
    {
        return create( &oper_le::apply, obj, r.obj );
    }
    SliceOperandPy __leL( double r )
    {
        return create( &oper_le::apply, obj, create( r ) );
    }
    SliceOperandPy __leR( double l )
    {
        return create( &oper_le::apply, create( l ), obj );
    }

    SliceOperandPy __gt( const SliceOperandPy & r )
    {
        return create( &oper_gt::apply, obj, r.obj );
    }
    SliceOperandPy __gtL( double r )
    {
        return create( &oper_gt::apply, obj, create( r ) );
    }
    SliceOperandPy __gtR( double l )
    {
        return create( &oper_gt::apply, create( l ), obj );
    }

    SliceOperandPy __ge( const SliceOperandPy & r )
    {
        return create( &oper_ge::apply, obj, r.obj );
    }
    SliceOperandPy __geL( double r )
    {
        return create( &oper_ge::apply, obj, create( r ) );
    }
    SliceOperandPy __geR( double l )
    {
        return create( &oper_ge::apply, create( l ), obj );
    }

    SliceOperandPy __or( const SliceOperandPy & r )
    {
        return create( &oper_or::apply, obj, r.obj );
    }
    SliceOperandPy __orL( bool r )
    {
        return create( &oper_or::apply, obj, create( r ) );
    }
    SliceOperandPy __orR( bool l )
    {
        return create( &oper_or::apply, create( l ), obj );
    }

    SliceOperandPy __and( const SliceOperandPy & r )
    {
        return create( &oper_and::apply, obj, r.obj );
    }
    SliceOperandPy __andL( bool r )
    {
        return create( &oper_and::apply, obj, create( r ) );
    }
    SliceOperandPy __andR( bool l )
    {
        return create( &oper_and::apply, create( l ), obj );
    }

    SliceOperandPy __not()
    {
        return create( &oper_not::apply, obj );
    }
};

class TreeSliceConstPy : public SliceOperandPy
{
public:
    const TreeSlice & obj;

    TreeSliceConstPy( const TreeSlice & obj ) :
        SliceOperandPy( SliceMarkerPySP( new TreeSliceOperandPy( obj ) ) ),
        obj( obj ) {}
};

class TreeSlicePy : public SliceOperandPy
{
    class EvalOper
    {
        const SliceMarkerPy & arg;
        TreeSlice & slice;

    public:
        EvalOper( const SliceMarkerPy & arg, TreeSlice & slice ) :
            sliceCount( arg.sliceCount + 1 ), arg( arg ), slice( slice ) {}

        const int sliceCount;

        const TreeSlice** listInputSlices( const TreeSlice** list ) const { return arg.listSlices( list ); }
        TreeSlice** listOutputSlices( TreeSlice** list ) const { *list = &slice; return ++list; }
        void compute() const { *slice.iter = arg.calc( (const double *)0 ); }

        void printDebug( char *s ) const {}
    };

    void eval( const SliceMarkerPy & arg ) { TreeSlice::loopOnSlices( EvalOper( arg, *obj ) ); }

public:
    TreeSliceSP obj;

    TreeSlicePy( const TreeSliceSP & obj ) :
        SliceOperandPy( SliceMarkerPySP( new TreeSliceOperandPy( *obj ) ) ),
        obj( obj ) {}

    void assignSliceMarker( const SliceOperandPy & src ) { eval( *src.obj ); }
    void assignTreeSliceConst( const TreeSliceConstPy & src ) { *obj = src.obj; }
    void assignTreeSlice( const TreeSlicePy & src ) { *obj = *src.obj; }
    void assignDouble( double x ) { *obj = x; }
    bool isZero() { return obj->isZero(); }
};

class TreeSliceSPPy
{
public:
    TreeSliceSP obj;

    TreeSlicePy getObject()
    {
        if( ! obj )
            throw ModelException( "TreeSliceSPPy::getObject()", "Object not initialized" );
        return obj;
    }
};

DRLIB_END_NAMESPACE

#endif
