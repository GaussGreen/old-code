//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceExpr.hpp
//
//   Description : Generic tree slice expressions.
//
//   Date        : Apr 12, 2006
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_EXPR_HPP
#define TREE_SLICE_EXPR_HPP

#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

template<class A>
const char* printExpr(void *p) {
    const A *a = (const A*)p;
    static char buf[1000];
    buf[0]=0;
    a->printDebug(buf);
    return buf;
}

//----------------------------------------------------------------------------
// slice marker
//----------------------------------------------------------------------------
template< typename T > 
class SliceMarker
{
public:
    typedef T DefType;
    typedef const T & RefType;
};

//----------------------------------------------------------------------------
// simple operand
//----------------------------------------------------------------------------
template< typename T >
class SimpleOperand : public SliceMarker< SimpleOperand< T > >
{
    T v;

public:
    SimpleOperand( T v ) : v( v ) {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 0;
    template< typename S >
    const S** listSlices(const S** list) const { return list; }
    inline T calc() const { return v; }
    void printDebug(char *s) const;
};

//----------------------------------------------------------------------------
// double operand
//----------------------------------------------------------------------------
typedef SimpleOperand< double > DoubleOperand;
template<>
class SliceMarker< double >
{
public:
    typedef DoubleOperand DefType;
    typedef DoubleOperand RefType;
};
template<>
inline
void DoubleOperand::printDebug(char *s) const
{
    char buf[30];
    sprintf(buf, "%f", v);
    strcat(s, buf);
    strcat(s, "{CST}");
}

//----------------------------------------------------------------------------
// bool operand
//----------------------------------------------------------------------------
typedef SimpleOperand< bool > BoolOperand;
template<>
class SliceMarker< bool >
{
public:
    typedef BoolOperand DefType;
    typedef BoolOperand RefType;
};
template<>
inline
void BoolOperand::printDebug(char *s) const
{
    char buf[10];
    sprintf(buf, "%c", v ? 'Y' : 'N' );
    strcat(s, buf);
    strcat(s, "{CST}");
}

//----------------------------------------------------------------------------
// unary slice expression
//----------------------------------------------------------------------------
template< typename OPER, typename ARG >
class SliceUnaryExpr : public SliceMarker< SliceUnaryExpr< OPER, ARG > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< ARG >::RefType ARGRefType;

    ARGRefType arg;

public:
    SliceUnaryExpr( ARGRefType arg ) : arg( arg ) {}

    // TreeSlice "expression template" primitives
    static const int sliceCount =
        SliceMarker< ARG >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { return arg.listSlices( list ); }
    inline OperType calc() const { return OPER::apply( arg.calc() ); }
    void printDebug(char *s) const
    {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        arg.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// binary slice expression
//----------------------------------------------------------------------------
template< typename OPER, typename L, typename R >
class SliceBinaryExpr : public SliceMarker< SliceBinaryExpr< OPER, L, R > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< L >::RefType LRefType;
    typedef typename SliceMarker< R >::RefType RRefType;

    LRefType l;
    RRefType r;

public:
    SliceBinaryExpr( LRefType l, RRefType r ) : l( l ), r( r ) {}

    // TreeSlice "expression template" primitives
    static const int sliceCount =
        SliceMarker< L >::DefType::sliceCount +
        SliceMarker< R >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { return r.listSlices( l.listSlices( list ) ); }
    inline OperType calc() const { return OPER::apply( l.calc(), r.calc() ); }
    void printDebug(char *s) const
    {
        strcat(s, "(");
        l.printDebug(s);
        strcat(s, OPER::symbol());
        r.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// 3 slice expression
//----------------------------------------------------------------------------
template<
    typename OPER,
    typename S0,
    typename S1,
    typename S2 >
class Slice3Expr : public SliceMarker< Slice3Expr< OPER, S0, S1, S2 > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< S0 >::RefType S0RefType;
    typedef typename SliceMarker< S1 >::RefType S1RefType;
    typedef typename SliceMarker< S2 >::RefType S2RefType;

    S0RefType s0;
    S1RefType s1;
    S2RefType s2;

public:
    Slice3Expr(
        S0RefType s0,
        S1RefType s1,
        S2RefType s2 )
        :
        s0( s0 ), s1( s1 ), s2( s2 )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 
        SliceMarker< S0 >::DefType::sliceCount + 
        SliceMarker< S1 >::DefType::sliceCount +
        SliceMarker< S2 >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { 
        return
            s2.listSlices(
                s1.listSlices(
                    s0.listSlices( list ) ) );
    }
    inline OperType calc() const
    {
        return OPER::apply(
            s0.calc(),
            s1.calc(),
            s2.calc() );
    }
    void printDebug(char *s) const {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        s0.printDebug(s);  strcat(s, ",");
        s1.printDebug(s);  strcat(s, ",");
        s2.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// 4 slice expression
//----------------------------------------------------------------------------
template<
    typename OPER,
    typename S0,
    typename S1,
    typename S2,
    typename S3 >
class Slice4Expr : public SliceMarker< Slice4Expr< OPER, S0, S1, S2, S3 > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< S0 >::RefType S0RefType;
    typedef typename SliceMarker< S1 >::RefType S1RefType;
    typedef typename SliceMarker< S2 >::RefType S2RefType;
    typedef typename SliceMarker< S3 >::RefType S3RefType;

    S0RefType s0;
    S1RefType s1;
    S2RefType s2;
    S3RefType s3;

public:
    Slice4Expr(
        S0RefType s0,
        S1RefType s1,
        S2RefType s2,
        S3RefType s3 )
        :
        s0( s0 ), s1( s1 ), s2( s2 ), s3( s3 )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 
        SliceMarker< S0 >::DefType::sliceCount + 
        SliceMarker< S1 >::DefType::sliceCount +
        SliceMarker< S2 >::DefType::sliceCount +
        SliceMarker< S3 >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { 
        return
            s3.listSlices(
                s2.listSlices(
                    s1.listSlices(
                        s0.listSlices( list ) ) ) );
    }
    inline OperType calc() const
    {
        return OPER::apply(
            s0.calc(),
            s1.calc(),
            s2.calc(),
            s3.calc() );
    }
    void printDebug(char *s) const {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        s0.printDebug(s);  strcat(s, ",");
        s1.printDebug(s);  strcat(s, ",");
        s2.printDebug(s);  strcat(s, ",");
        s3.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// 5 slice expression
//----------------------------------------------------------------------------
template<
    typename OPER,
    typename S0,
    typename S1,
    typename S2,
    typename S3,
    typename S4 >
class Slice5Expr : public SliceMarker< Slice5Expr< OPER, S0, S1, S2, S3, S4 > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< S0 >::RefType S0RefType;
    typedef typename SliceMarker< S1 >::RefType S1RefType;
    typedef typename SliceMarker< S2 >::RefType S2RefType;
    typedef typename SliceMarker< S3 >::RefType S3RefType;
    typedef typename SliceMarker< S4 >::RefType S4RefType;

    S0RefType s0;
    S1RefType s1;
    S2RefType s2;
    S3RefType s3;
    S4RefType s4;

public:
    Slice5Expr(
        S0RefType s0,
        S1RefType s1,
        S2RefType s2,
        S3RefType s3,
        S4RefType s4 )
        :
        s0( s0 ), s1( s1 ), s2( s2 ), s3( s3 ), s4( s4 )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 
        SliceMarker< S0 >::DefType::sliceCount +
        SliceMarker< S1 >::DefType::sliceCount +
        SliceMarker< S2 >::DefType::sliceCount +
        SliceMarker< S3 >::DefType::sliceCount +
        SliceMarker< S4 >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { 
        return
            s4.listSlices(
                s3.listSlices(
                    s2.listSlices(
                        s1.listSlices(
                            s0.listSlices( list ) ) ) ) );
    }
    inline OperType calc() const
    {
        return OPER::apply(
            s0.calc(),
            s1.calc(),
            s2.calc(),
            s3.calc(),
            s4.calc() );
    }
    void printDebug(char *s) const {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        s0.printDebug(s);  strcat(s, ",");
        s1.printDebug(s);  strcat(s, ",");
        s2.printDebug(s);  strcat(s, ",");
        s3.printDebug(s);  strcat(s, ",");
        s4.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// 6 slice expression
//----------------------------------------------------------------------------
template<
    typename OPER,
    typename S0,
    typename S1,
    typename S2,
    typename S3,
    typename S4,
    typename S5 >
class Slice6Expr : public SliceMarker< Slice6Expr< OPER, S0, S1, S2, S3, S4, S5 > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< S0 >::RefType S0RefType;
    typedef typename SliceMarker< S1 >::RefType S1RefType;
    typedef typename SliceMarker< S2 >::RefType S2RefType;
    typedef typename SliceMarker< S3 >::RefType S3RefType;
    typedef typename SliceMarker< S4 >::RefType S4RefType;
    typedef typename SliceMarker< S5 >::RefType S5RefType;

    S0RefType s0;
    S1RefType s1;
    S2RefType s2;
    S3RefType s3;
    S4RefType s4;
    S5RefType s5;

public:
    Slice6Expr(
        S0RefType s0,
        S1RefType s1,
        S2RefType s2,
        S3RefType s3,
        S4RefType s4,
        S5RefType s5 )
        :
        s0( s0 ), s1( s1 ), s2( s2 ), s3( s3 ), s4( s4 ), s5( s5 )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 
        SliceMarker< S0 >::DefType::sliceCount + 
        SliceMarker< S1 >::DefType::sliceCount +
        SliceMarker< S2 >::DefType::sliceCount +
        SliceMarker< S3 >::DefType::sliceCount +
        SliceMarker< S4 >::DefType::sliceCount +
        SliceMarker< S5 >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { 
        return
            s5.listSlices(
                s4.listSlices(
                    s3.listSlices(
                        s2.listSlices(
                            s1.listSlices(
                                s0.listSlices( list ) ) ) ) ) );
    }
    inline OperType calc() const
    {
        return OPER::apply(
            s0.calc(),
            s1.calc(),
            s2.calc(),
            s3.calc(),
            s4.calc(),
            s5.calc() );
    }
    void printDebug(char *s) const {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        s0.printDebug(s);  strcat(s, ",");
        s1.printDebug(s);  strcat(s, ",");
        s2.printDebug(s);  strcat(s, ",");
        s3.printDebug(s);  strcat(s, ",");
        s4.printDebug(s);  strcat(s, ",");
        s5.printDebug(s);
        strcat(s, ")");
    }
};

//----------------------------------------------------------------------------
// 7 slice expression
//----------------------------------------------------------------------------
template<
    typename OPER,
    typename S0,
    typename S1,
    typename S2,
    typename S3,
    typename S4,
    typename S5,
    typename S6 >
class Slice7Expr : public SliceMarker< Slice7Expr< OPER, S0, S1, S2, S3, S4, S5, S6 > >
{
    typedef typename OPER::Type OperType;
    typedef typename SliceMarker< S0 >::RefType S0RefType;
    typedef typename SliceMarker< S1 >::RefType S1RefType;
    typedef typename SliceMarker< S2 >::RefType S2RefType;
    typedef typename SliceMarker< S3 >::RefType S3RefType;
    typedef typename SliceMarker< S4 >::RefType S4RefType;
    typedef typename SliceMarker< S5 >::RefType S5RefType;
    typedef typename SliceMarker< S6 >::RefType S6RefType;

    S0RefType s0;
    S1RefType s1;
    S2RefType s2;
    S3RefType s3;
    S4RefType s4;
    S5RefType s5;
    S6RefType s6;

public:
    Slice7Expr(
        S0RefType s0,
        S1RefType s1,
        S2RefType s2,
        S3RefType s3,
        S4RefType s4,
        S5RefType s5,
        S6RefType s6 )
        :
        s0( s0 ), s1( s1 ), s2( s2 ), s3( s3 ), s4( s4 ), s5( s5 ), s6( s6 )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 
        SliceMarker< S0 >::DefType::sliceCount + 
        SliceMarker< S1 >::DefType::sliceCount +
        SliceMarker< S2 >::DefType::sliceCount +
        SliceMarker< S3 >::DefType::sliceCount +
        SliceMarker< S4 >::DefType::sliceCount +
        SliceMarker< S5 >::DefType::sliceCount +
        SliceMarker< S6 >::DefType::sliceCount;
    template< typename S >
    const S** listSlices(const S** list) const { 
        return
            s6.listSlices(
                s5.listSlices(
                    s4.listSlices(
                        s3.listSlices(
                            s2.listSlices(
                                s1.listSlices(
                                    s0.listSlices( list ) ) ) ) ) ) );
    }
    inline OperType calc() const
    {
        return OPER::apply(
            s0.calc(),
            s1.calc(),
            s2.calc(),
            s3.calc(),
            s4.calc(),
            s5.calc(),
            s6.calc() );
    }
    void printDebug(char *s) const {
        strcat(s, OPER::symbol());
        strcat(s, "(");
        s0.printDebug(s);  strcat(s, ",");
        s1.printDebug(s);  strcat(s, ",");
        s2.printDebug(s);  strcat(s, ",");
        s3.printDebug(s);  strcat(s, ",");
        s4.printDebug(s);  strcat(s, ",");
        s5.printDebug(s);  strcat(s, ",");
        s6.printDebug(s);
        strcat(s, ")");
    }
};

DRLIB_END_NAMESPACE

#endif
