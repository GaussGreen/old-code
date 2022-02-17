#ifndef fix123_oper_dot_hpp
#define fix123_oper_dot_hpp

#include <fix123head.h>

typedef double value_type;

// -------------------------------------------------------------------------
//                            Supported Slice Operations                            
// -------------------------------------------------------------------------
struct UnaryOper
{
    virtual void execute(value_type* lhs, value_type const* x, size_t cnt) const=0;
    virtual ~UnaryOper() {}
};

struct UnarySmoothOper
{
    // can turn off smoothing at runtime by passing false
    UnarySmoothOper(bool smooth=true) : m_smooth(smooth) {}
    virtual void execute(value_type* lhs, value_type const* x, value_type const* step, size_t cnt) const=0;
    virtual ~UnarySmoothOper() {}

    bool    m_smooth;
};

struct BinaryOper
{
    virtual void execute(value_type* lhs, value_type const* x, value_type const* y, size_t cnt) const=0;
    virtual ~BinaryOper() {}
};


// To be used in operations that on a set of same offset nodes from different slices
struct XSUnaryOper
{
    virtual void execute(value_type* lhs, value_type const* rhs, size_t cnt, value_type x) const=0;
    virtual ~XSUnaryOper() {}
};

// To be used in operations that on a set of same offset nodes from different slices
struct XSBinaryOper
{
    virtual void execute(value_type* lhs, value_type const* rhs, size_t cnt, value_type x, value_type y) const=0;
    virtual ~XSBinaryOper() {}
};


//------------------------------------------------------------------------------
// Unary slice operation
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, 
               value_type const* arg, UnaryOper const& fun);

//------------------------------------------------------------------------------
// Unary slice smooth operation
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, 
               value_type const* arg, 
               UnarySmoothOper const& fun);

//------------------------------------------------------------------------------
// Binary slice operation - functor version
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, value_type const* arg1, value_type const* arg2, 
               BinaryOper const& fun);


#endif

