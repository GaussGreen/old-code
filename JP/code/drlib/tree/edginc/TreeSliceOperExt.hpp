//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : TreeSliceOperExt.hpp
//
//   Description : Special tree slice operation.
//
//   Date        : Jan 24, 2005
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_OPER_EXT_HPP
#define TREE_SLICE_OPER_EXT_HPP

#include "edginc/UtilFuncs.hpp"
#include "edginc/Black.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------
// BS price expression and operator
//----------------------------------------------------------------------------
template< typename FWD >
class SliceBSExpr : public SliceMarker< SliceBSExpr< FWD > >
{
    typedef typename SliceMarker< FWD >::RefType FWDRefType;

    const bool isCall;
    FWDRefType fwd;
    const double strike;
    const double pv;
    const double variance;

    double price( double fwd ) const
    {
        return Black::price( isCall, fwd, strike, pv, variance );
    }

public:
    SliceBSExpr( bool isCall, FWDRefType fwd, double strike, double pv, double variance ) :
        isCall( isCall ),
        fwd( fwd ),
        strike( strike ),
        pv( pv ),
        variance( variance )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = SliceMarker< FWD >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { return fwd.listSlices( list ); }
    inline double calc() const { return price( fwd.calc() ); }
    void printDebug(char *s) const
    {
        char buf[200];
        sprintf(buf, "BS(%s, K=%f, PV=%f, Var=%f, Fwd=", isCall ? "CALL" : "PUT", strike, pv, variance);
        strcat(s, buf);
        this->fwd.printDebug(s);
        strcat(s, ")");
    }
};
template< typename FWD >
inline
SliceBSExpr< FWD > priceBS( bool isCall, const FWD & fwd, double strike, double pv, double variance )
{
    return SliceBSExpr< FWD >( isCall, fwd, strike, pv, variance );
}

//----------------------------------------------------------------------------
// accumulation slice expression
//----------------------------------------------------------------------------
template< typename OPER, typename T, typename ARG >
class SliceAccumExpr : public SliceMarker< SliceAccumExpr< OPER, T, ARG > >
{
    typedef typename SliceMarker< ARG >::RefType ARGRefType;

    double & value;
    ARGRefType arg;

public:
    static double apply( const T & slice, double value, ARGRefType arg )
    {
        TreeSliceSP dummy = slice.clone( false );
        *dummy = SliceAccumExpr< OPER, T, ARG >( value, arg );
        return value;
    }

    SliceAccumExpr( double & value, ARGRefType arg ) : value( value ), arg( arg ) {}

    // TreeSlice "expression template" primitives
    static const int sliceCount =
        SliceMarker< ARG >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { return arg.listSlices( list ); }
    inline double calc() const { return value = OPER::apply( value, arg.calc() ); }
    void printDebug(char *s) const
    {
        strcat(s, "Accum(");
        strcat(s, OPER::symbol());
        strcat(s, "(");
        arg.printDebug(s);
        strcat(s, "))");
    }
};

//----------------------------------------------------------------------------
// accumulation sum - returns sum of values in the slice
//----------------------------------------------------------------------------
template< typename T, typename ARG >
inline
double esum( const T & slice, const SliceMarker< ARG > & arg )
{
    return SliceAccumExpr< oper_add, T, ARG >::apply( slice, 0., static_cast< const ARG & >( arg ) );
}

//----------------------------------------------------------------------------
// accumulation prod - returns product of values in the slice
//----------------------------------------------------------------------------
template< typename T, typename ARG >
inline
double eprod( const T & slice, const SliceMarker< ARG > & arg )
{
    return SliceAccumExpr< oper_mul, T, ARG >::apply( slice, 1., static_cast< const ARG & >( arg ) );
}

//----------------------------------------------------------------------------
// accumulation max - returns maximum value in the slice
//----------------------------------------------------------------------------
template< typename T, typename ARG >
inline
double emax( const T & slice, const SliceMarker< ARG > & arg )
{
    return SliceAccumExpr< oper_max, T, ARG >::apply( slice, -DBL_MAX, static_cast< const ARG & >( arg ) );
}

//----------------------------------------------------------------------------
// accumulation min - returns maximum value in the slice
//----------------------------------------------------------------------------
template< typename T, typename ARG >
inline
double emin( const T & slice, const SliceMarker< ARG > & arg )
{
    return SliceAccumExpr< oper_min, T, ARG >::apply( slice, DBL_MAX, static_cast< const ARG & >( arg ) );
}

//----------------------------------------------------------------------------
// Knock-out-type operators
//----------------------------------------------------------------------------
struct oper_ko_ternary_boundary_in : SliceOper< double >
{
    static const char* symbol(){return "oper_ko_ternary_boundary_in";}
    static Type apply( double value, double lo, double hi, double up, double mid, double down) 
    {
        if (value < lo )
            return down;
        if ( hi < value)
            return up;
        return mid;
    }
};

struct oper_ko_ternary_boundary_out : SliceOper< double >
{
    static const char* symbol(){return "oper_ko_ternary_boundary_out";}
    static Type apply( double value, double lo, double hi, double up, double mid, double down) 
    {
        if (value <= lo )
            return down;
        if ( hi <= value)
            return up;
        return mid;
    }
};

struct oper_ko_smooth_ternary : SliceOper< double >
{
    static const char* symbol(){return "oper_ko_smooth_ternary";}
    static Type apply( double value, double lo, double hi,
                        double up, double mid, double down,
                        double smoothStep) 
    {        
        double x = SmoothValue( mid,  lo + smoothStep, 
                                down, lo - smoothStep,
                                value);

        return SmoothValue( up, hi + smoothStep,
                            x,  hi - smoothStep,
                            value);
    }
};

template <typename S0, typename S1, typename S2>
void koTernary(TreeSlice &result, 
               const TreeSlice &value,
               double lo,
               double hi,
               const S0 &up,
               const S1 &mid,
               const S2 &down,
               bool boundaryIn,
               bool smooth) 
{
    try {
        typedef const TreeSlice T;
        if (smooth) {
            TreeSliceSP smoothStep = value.calcSmoothStep();

            result = Slice7Expr<oper_ko_smooth_ternary, T, double, double, S0, S1, S2, T>
                (value, lo, hi, 
                up, mid, down, *smoothStep);
            return;
        }
        if (boundaryIn) {
            result = Slice6Expr<oper_ko_ternary_boundary_in, T, double, double, S0, S1, S2>(
                value, lo, hi, 
                up, mid, down);
        }
        else {
            result = Slice6Expr<oper_ko_ternary_boundary_out, T, double, double, S0, S1, S2>(
                value, lo, hi, 
                up, mid, down);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

DRLIB_END_NAMESPACE

#endif
