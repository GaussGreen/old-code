//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceOper.hpp
//
//   Description : Generic tree slice operations.
//
//   Date        : Oct 14, 2005
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_OPER_HPP
#define TREE_SLICE_OPER_HPP

#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------
// slice operator base
//----------------------------------------------------------------------------
template< typename T > 
struct SliceOper
{
    typedef T Type;
};

//----------------------------------------------------------------------------
// binary slice operators
//----------------------------------------------------------------------------
// operator +
struct oper_add : SliceOper< double >
{
    static Type apply( double a, double b ) { return a + b; }
    static const char* symbol(){return " + ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_add, L, R > add( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_add, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_add, L, R > operator +( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return add( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_add, L, double > operator +( const SliceMarker< L > & l, double r )
{
    return add( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_add, double, R > operator +( double l, const SliceMarker< R > & r )
{
    return add( l, static_cast< const R & >( r ) );
}

// operator -
struct oper_sub : SliceOper< double >
{
    static Type apply( double a, double b ) { return a - b; }
    static const char* symbol(){return " - ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_sub, L, R > sub( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_sub, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_sub, L, R > operator -( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return sub( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_sub, L, double > operator -( const SliceMarker< L > & l, double r )
{
    return sub( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_sub, double, R > operator -( double l, const SliceMarker< R > & r )
{
    return sub( l, static_cast< const R & >( r ) );
}

// operator *
struct oper_mul : SliceOper< double >
{
    static Type apply( double a, double b ) { return a * b; }
    static const char* symbol(){return " * ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_mul, L, R > mul( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_mul, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_mul, L, R > operator *( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return mul( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_mul, L, double > operator *( const SliceMarker< L > & l, double r )
{
    return mul( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_mul, double, R > operator *( double l, const SliceMarker< R > & r )
{
    return mul( l, static_cast< const R & >( r ) );
}

// operator /
struct oper_div : SliceOper< double >
{
    static Type apply( double a, double b ) { return a / b; }
    static const char* symbol(){return " / ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_div, L, R > div( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_div, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_div, L, R > operator /( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return div( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_div, L, double > operator /( const SliceMarker< L > & l, double r )
{
    return div( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_div, double, R > operator /( double l, const SliceMarker< R > & r )
{
    return div( l, static_cast< const R & >( r ) );
}

// operator min
struct oper_min : SliceOper< double >
{
    static Type apply( double a, double b ) { return a > b ? b : a; }
    static const char* symbol(){return " min ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_min, L, R > smin( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_min, L, R >( l, r );
}

// operator max
struct oper_max : SliceOper< double >
{
    static Type apply( double a, double b ) { return a < b ? b : a; }
    static const char* symbol(){return " max ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_max, L, R > smax( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_max, L, R >( l, r );
}

// operator pow
struct oper_pow : SliceOper< double >
{
    static Type apply( double a, double b ) { return ::pow( a, b ); }
    static const char* symbol(){return " pow ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_pow, L, R > spow( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_pow, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_pow, L, R > pow( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return spow( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_pow, L, double > pow( const SliceMarker< L > & l, double r )
{
    return spow( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_pow, double, R > pow( double l, const SliceMarker< R > & r )
{
    return spow( l, static_cast< const R & >( r ) );
}

//----------------------------------------------------------------------------
// unary slice operators
//----------------------------------------------------------------------------
// operator +
struct oper_pos : SliceOper< double >
{
    static Type apply( double v ) { return +v; }
    static const char* symbol(){return " + ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_pos, ARG > pos( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_pos, ARG >( static_cast< const ARG & >( arg ) );
}
template< typename ARG >
inline
SliceUnaryExpr< oper_pos, ARG > operator +( const SliceMarker< ARG > & arg )
{
    return pos( static_cast< const ARG & >( arg ) );
}

// operator -
struct oper_neg : SliceOper< double >
{
    static Type apply( double v ) { return -v; }
    static const char* symbol(){return " - ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_neg, ARG > neg( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_neg, ARG >( static_cast< const ARG & >( arg ) );
}
template< typename ARG >
inline
SliceUnaryExpr< oper_neg, ARG > operator -( const SliceMarker< ARG > & arg )
{
    return neg( static_cast< const ARG & >( arg ) );
}

// operator abs
struct oper_abs : SliceOper< double >
{
    static Type apply( double v ) { return ::fabs( v ); }
    static const char* symbol(){return " abs ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_abs, ARG > abs( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_abs, ARG >( static_cast< const ARG & >( arg ) );
}

// operator log
struct oper_log : SliceOper< double >
{
    static Type apply( double v ) { return ::log( v ); }
    static const char* symbol(){return " log ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_log, ARG > log( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_log, ARG >( static_cast< const ARG & >( arg ) );
}

// operator exp
struct oper_exp : SliceOper< double >
{
    static Type apply( double v ) { return ::exp( v ); }
    static const char* symbol(){return " exp ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_exp, ARG > exp( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_exp, ARG >( static_cast< const ARG & >( arg ) );
}

// operator sqrt
struct oper_sqrt : SliceOper< double >
{
    static Type apply( double v ) { return ::sqrt( v ); }
    static const char* symbol(){return " sqrt ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_exp, ARG > sqrt( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_sqrt, ARG >( static_cast< const ARG & >( arg ) );
}


//----------------------------------------------------------------------------
// logical slice operators
//----------------------------------------------------------------------------
//!!! REVIEW == AND != BELOW
/*
// operator ==
struct oper_eq : SliceOper< bool >
{
    static Type apply( double a, double b ) { return Maths::equals( a, b ); }
    static const char* symbol(){return " == ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_eq, L, R > eq( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_eq, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_eq, L, R > operator ==( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return eq( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_eq, L, double > operator ==( const SliceMarker< L > & l, double r )
{
    return eq( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_eq, double, R > operator ==( double l, const SliceMarker< R > & r )
{
    return eq( l, static_cast< const R & >( r ) );
}

// operator !=
struct oper_ne : SliceOper< bool >
{
    static Type apply( double a, double b ) { return ! Maths::equals( a, b ); }
    static const char* symbol(){return " != ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_ne, L, R > ne( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_ne, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_ne, L, R > operator !=( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return ne( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_ne, L, double > operator !=( const SliceMarker< L > & l, double r )
{
    return ne( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_ne, double, R > operator !=( double l, const SliceMarker< R > & r )
{
    return ne( l, static_cast< const R & >( r ) );
}
*/

// operator <
struct oper_lt : SliceOper< bool >
{
    static Type apply( double a, double b ) { return a < b; }
    static const char* symbol(){return " < ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_lt, L, R > lt( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_lt, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_lt, L, R > operator <( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return lt( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_lt, L, double > operator <( const SliceMarker< L > & l, double r )
{
    return lt( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_lt, double, R > operator <( double l, const SliceMarker< R > & r )
{
    return lt( l, static_cast< const R & >( r ) );
}

// operator <=
struct oper_le : SliceOper< bool >
{
    static Type apply( double a, double b ) { return a <= b; }
    static const char* symbol(){return " <= ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_le, L, R > le( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_le, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_le, L, R > operator <=( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return le( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_le, L, double > operator <=( const SliceMarker< L > & l, double r )
{
    return le( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_le, double, R > operator <=( double l, const SliceMarker< R > & r )
{
    return le( l, static_cast< const R & >( r ) );
}

// operator >
struct oper_gt : SliceOper< bool >
{
    static Type apply( double a, double b ) { return a > b; }
    static const char* symbol(){return " > ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_gt, L, R > gt( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_gt, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_gt, L, R > operator >( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return gt( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_gt, L, double > operator >( const SliceMarker< L > & l, double r )
{
    return gt( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_gt, double, R > operator >( double l, const SliceMarker< R > & r )
{
    return gt( l, static_cast< const R & >( r ) );
}

// operator >=
struct oper_ge : SliceOper< bool >
{
    static Type apply( double a, double b ) { return a >= b; }
    static const char* symbol(){return " >= ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_ge, L, R > ge( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_ge, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_ge, L, R > operator >=( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return ge( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_ge, L, double > operator >=( const SliceMarker< L > & l, double r )
{
    return ge( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_ge, double, R > operator >=( double l, const SliceMarker< R > & r )
{
    return ge( l, static_cast< const R & >( r ) );
}

// operator ||
struct oper_or : SliceOper< bool >
{
    static Type apply( bool a, bool b ) { return a || b; }
    static const char* symbol(){return " || ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_or, L, R > sor( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_or, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_or, L, R > operator ||( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return sor( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_or, L, bool > operator ||( const SliceMarker< L > & l, bool r )
{
    return sor( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_or, bool, R > operator ||( bool l, const SliceMarker< R > & r )
{
    return sor( l, static_cast< const R & >( r ) );
}

// operator &&
struct oper_and : SliceOper< bool >
{
    static Type apply( bool a, bool b ) { return a && b; }
    static const char* symbol(){return " && ";}
};
template< typename L, typename R >
inline
SliceBinaryExpr< oper_and, L, R > sand( const L & l, const R & r )
{
    return SliceBinaryExpr< oper_and, L, R >( l, r );
}
template< typename L, typename R >
inline
SliceBinaryExpr< oper_and, L, R > operator &&( const SliceMarker< L > & l, const SliceMarker< R > & r )
{
    return sand( static_cast< const L & >( l ), static_cast< const R & >( r ) );
}
template< typename L >
inline
SliceBinaryExpr< oper_and, L, bool > operator &&( const SliceMarker< L > & l, bool r )
{
    return sand( static_cast< const L & >( l ), r );
}
template< typename R >
inline
SliceBinaryExpr< oper_and, bool, R > operator &&( bool l, const SliceMarker< R > & r )
{
    return sand( l, static_cast< const R & >( r ) );
}

// operator !
struct oper_not : SliceOper< bool >
{
    static Type apply( bool v ) { return ! v; }
    static const char* symbol(){return " ! ";}
};
template< typename ARG >
inline
SliceUnaryExpr< oper_not, ARG > snot( const SliceMarker< ARG > & arg )
{
    return SliceUnaryExpr< oper_not, ARG >( static_cast< const ARG & >( arg ) );
}
template< typename ARG >
inline
SliceUnaryExpr< oper_not, ARG > operator !( const SliceMarker< ARG > & arg )
{
    return snot( static_cast< const ARG & >( arg ) );
}

struct oper_cond : SliceOper< double >
{
    static Type apply( bool c, double a, double b )
    {
        return c ? a : b;
    }
    static const char* symbol(){return " cond ";}
};
template< typename C, typename A, typename B >
inline
Slice3Expr< oper_cond, C, A, B >
cond( const C & c, const A & a, const B & b )
{
    return Slice3Expr< oper_cond, C, A, B >( c, a, b );
}

#define IF( x ) cond( x,
#define ELSE ,
#define ENDIF )

DRLIB_END_NAMESPACE

#endif
