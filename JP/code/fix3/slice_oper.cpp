#include <fix123head.h>
#include <fix123oper.hpp>
#include <assert.h>
#include <math.h>
#include <vector>

//------------------------------------------------------------------------------
// Unary slice operation
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, 
               value_type const* arg, UnaryOper const& fun)
{
//    static char const*  routine = "Fix3Tree::sliceOper unary";

    int     i, j;
    int     offset;

    // source and destination slice element pointers
    value_type const*   srcPtr;
    value_type*         trgPtr;

    // Slice boundary
    int     Top1 = tree->Top1   [t];
    int*    Top2 = tree->Top2   [t];
    int**   Top3 = tree->Top3   [t];
    int     Bot1 = tree->Bottom1[t];
    int*    Bot2 = tree->Bottom2[t];
    int**   Bot3 = tree->Bottom3[t];

//std::cout << "UNARY OPER" << std::endl;
//std::cout << "-LHS" << *lhs << std::endl;
//std::cout << "ARG " << *arg << std::endl;
    switch (tree->NbFactor) {
        case 1: 
            offset = Fix3_Node_Offset(1, 0, 0, t, tree) + Bot1;
            srcPtr = arg + offset; 
            trgPtr = lhs + offset;

            fun.execute(trgPtr, srcPtr, Top1 - Bot1 + 1);
//std::cout << "+LHS (" << Top1-Bot1+1 << ")" << *lhs << std::endl;
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Fix3_Node_Offset(2, i, 0, t, tree) + Bot2[i];
                srcPtr = arg + offset; 
                trgPtr = lhs + offset;

                fun.execute(trgPtr, srcPtr, Top2[i] - Bot2[i] + 1);
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Fix3_Node_Offset(3, i, j, t, tree) + Bot3[i][j];
                    srcPtr = arg + offset; 
                    trgPtr = lhs + offset;

                    fun.execute(trgPtr, srcPtr, Top3[i][j] - Bot3[i][j] + 1);
                }
            }
            break;
    } 
}



//------------------------------------------------------------------------------
// Unary slice smooth operation
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, 
               value_type const* arg, 
               UnarySmoothOper const& fun)
{
//    static char const*  routine = "Fix3Tree::sliceOper smooth";

    int     i, j, k;
    int     offset;

    // for smoothing function lhs != arg
    assert(lhs != arg);

    // source and destination slice element pointers
    value_type const*   srcPtr;
    value_type*         trgPtr;
    value_type const*   refPtr;

    // Slice boundary
    int     Top1 = tree->Top1   [t];
    int*    Top2 = tree->Top2   [t];
    int**   Top3 = tree->Top3   [t];
    int     Bot1 = tree->Bottom1[t];
    int*    Bot2 = tree->Bottom2[t];
    int**   Bot3 = tree->Bottom3[t];


    value_type const* ref = arg;

//std::cout << "UNARY SMOOTH OPER" << std::endl;
    switch (tree->NbFactor) {
        case 1: 
            offset = Fix3_Node_Offset(1, 0, 0, t, tree) + Bot1;
            {
                srcPtr = arg + offset; 
                trgPtr = lhs + offset;

                refPtr = ref;

                // build a vector of smoothing steps
                int len = Top1 - Bot1 + 1;
                int idx = 0;
                std::vector<value_type> steps(len);

                if (fun.m_smooth)
                    for (i = Bot1; i <= Top1; i++)
                        steps[idx++] = Fix3_GetIndexStep(refPtr, 1, i, 0, 0, t, tree);

                fun.execute(trgPtr, srcPtr, &steps[0], len);
            }
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Fix3_Node_Offset(2, i, 0, t, tree) + Bot2[i];
                {
                    srcPtr = arg + offset; 
                    trgPtr = lhs + offset;

                    refPtr = ref;

                    // build a vector of smoothing steps
                    int len = Top2[i] - Bot2[i] + 1;
                    int idx = 0;
                    std::vector<value_type> steps(len);

                    if (fun.m_smooth)
                        for (j = Bot2[i]; j <= Top2[i]; j++)
                            steps[idx++] = Fix3_GetIndexStep(refPtr, 2, i, j, 0, t, tree);

                    fun.execute(trgPtr, srcPtr, &steps[0], len);
                }
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Fix3_Node_Offset(3, i, j, t, tree) + Bot3[i][j];
                    {
                        srcPtr = arg + offset; 
                        trgPtr = lhs + offset;

                        refPtr = ref;

                        // build a vector of smoothing steps
                        int len = Top3[i][j] - Bot3[i][j] + 1;
                        int idx = 0;
                        std::vector<value_type> steps(len);

                        if (fun.m_smooth)
                            for (k = Bot3[i][j]; k <= Top3[i][j]; k++)
                                steps[idx++] = Fix3_GetIndexStep(refPtr, 3, i, j, k, t, tree);

                        fun.execute(trgPtr, srcPtr, &steps[0], len);
                    }
                }
            }
            break;
    } 
}


//------------------------------------------------------------------------------
// Binary slice operation - functor version
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type* lhs, value_type const* arg1, value_type const* arg2, 
               BinaryOper const& fun)
{
//    static char const*  routine = "Fix3Tree::sliceOper binary";

    int     i, j;
    int     offset;

    // source and destination slice element pointers
    value_type const*   src1Ptr;
    value_type const*   src2Ptr;
    value_type*         trgPtr;

    // Slice boundary
    int     Top1 = tree->Top1   [t];
    int*    Top2 = tree->Top2   [t];
    int**   Top3 = tree->Top3   [t];
    int     Bot1 = tree->Bottom1[t];
    int*    Bot2 = tree->Bottom2[t];
    int**   Bot3 = tree->Bottom3[t];

/*
std::cout << "BINARY OPER dim=" << lhs->mNbDim << std::endl;
std::cout << "-LHS " << *lhs << std::endl;
std::cout << "ARG1 " << *arg1 << std::endl;
std::cout << "ARG2 " << *arg2 << std::endl;
*/
    switch (tree->NbFactor) {
        case 1: 
            offset = Fix3_Node_Offset(1, 0, 0, t, tree) + Bot1;
            {
                trgPtr =  lhs + offset;
                src1Ptr= arg1 + offset; 
                src2Ptr= arg2 + offset; 

                fun.execute(trgPtr, src1Ptr, src2Ptr, Top1 - Bot1 + 1);
            }
//std::cout << "+LHS (" << Top1-Bot1+1 << ")" << *lhs << std::endl;
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Fix3_Node_Offset(2, i, 0, t, tree) + Bot2[i];
                {
                    trgPtr =  lhs + offset;
                    src1Ptr= arg1 + offset; 
                    src2Ptr= arg2 + offset; 

                    fun.execute(trgPtr, src1Ptr, src2Ptr, Top2[i] - Bot2[i] + 1);
                }
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Fix3_Node_Offset(3, i, j, t, tree) + Bot3[i][j];
                    {
                        trgPtr =  lhs + offset;
                        src1Ptr= arg1 + offset; 
                        src2Ptr= arg2 + offset; 

                        fun.execute(trgPtr, src1Ptr, src2Ptr, Top3[i][j] - Bot3[i][j] + 1);
                    }
                }
            }
            break;
    } 
}




//------------------------------------------------------------------------------
// Unary xsectional slice operation - functor version
// NOTE: Some explanation is in order so here is how this function and its 
// binary version work.
// Functions operate on a set of slices at once. Given the bank of count argument 
// slices (rhs), 'count' of result slices (lhs), one (two for binary) argument 
// slices x (x and y for binary) and a function object (fun) the function does 
// the following:
//
// 1. Create 'count' size input array by selecting same offset nodes from the
//    'rhs' set
// 2. Create 'count' size empty output array
// 3. Select same offset nodes from one or two argument slices
// 4. Invoke user supplied function with two arrays and one or two arguments.
// 5. Copy outpur array into appropriate slots in the 'lhs' set of slices.
// 
// The above operations are performed for all nodes. Result and argument slice 
// sets can be the same.
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type** lhs, value_type* const* rhs, size_t count,
               value_type const* x, XSUnaryOper const& fun)
{
    int     i, j, k, offset;
    size_t  c;

    // source and destination slice element pointers
    value_type const*   srcPtr;
    value_type const*     xPtr;
    value_type*         trgPtr;

    // Slice boundary
    int     Top1 = tree->Top1   [t];
    int*    Top2 = tree->Top2   [t];
    int**   Top3 = tree->Top3   [t];
    int     Bot1 = tree->Bottom1[t];
    int*    Bot2 = tree->Bottom2[t];
    int**   Bot3 = tree->Bottom3[t];

    std::vector<value_type> ixs(count);
    std::vector<value_type> oxs(count);

    switch (tree->NbFactor) {
        case 1: 
            offset = Fix3_Node_Offset(1, 0, 0, t, tree);
            for (i = Bot1; i <= Top1; i++)
            {
                for (c=0; c<count; ++c) 
                {
                    srcPtr = rhs[c] + offset; 
                    ixs[c] = srcPtr[i];
                }

                xPtr = x + offset;
                fun.execute(&oxs[0], &ixs[0], count, *xPtr);

                for (c=0; c<count; ++c) 
                {
                    trgPtr = lhs[c] + offset; 
                    trgPtr[i] = oxs[c];
                }
            }
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Fix3_Node_Offset(2, i, 0, t, tree);
                for (j = Bot2[i]; j <= Top2[i]; j++)
                {
                    for (c=0; c<count; ++c) 
                    {
                        srcPtr = rhs[c] + offset; 
                        ixs[c] = srcPtr[j];
                    }

                    xPtr = x + offset;
                    fun.execute(&oxs[0], &ixs[0], count, *xPtr);

                    for (c=0; c<count; ++c) 
                    {
                        trgPtr = lhs[c] + offset; 
                        trgPtr[j] = oxs[c];
                    }
                }
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Fix3_Node_Offset(3, i, j, t, tree);
                    for (k = Bot3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (c=0; c<count; ++c) 
                        {
                            srcPtr = rhs[c] + offset; 
                            ixs[c] = srcPtr[k];
                        }

                        xPtr = x + offset;
                        fun.execute(&oxs[0], &ixs[0], count, *xPtr);

                        for (c=0; c<count; ++c) 
                        {
                            trgPtr = lhs[c] + offset; 
                            trgPtr[k] = oxs[c];
                        }
                    }
                }
            }
            break;
    } 
}


//------------------------------------------------------------------------------
// Binary xsectional slice operation - functor version
//------------------------------------------------------------------------------
void
Fix3_SliceOper(FIX3_TREE_DATA const* tree, size_t t,
               value_type** lhs, value_type* const* rhs, size_t count,
               value_type const* x, value_type const* y, XSBinaryOper const& fun)
{
    int     i, j, k, offset;
    size_t  c;

    // source and destination slice element pointers
    value_type const*   srcPtr;
    value_type const*     xPtr;
    value_type const*     yPtr;
    value_type*         trgPtr;

    // Slice boundary
    int     Top1 = tree->Top1   [t];
    int*    Top2 = tree->Top2   [t];
    int**   Top3 = tree->Top3   [t];
    int     Bot1 = tree->Bottom1[t];
    int*    Bot2 = tree->Bottom2[t];
    int**   Bot3 = tree->Bottom3[t];

    std::vector<value_type> ixs(count);
    std::vector<value_type> oxs(count);

    switch (tree->NbFactor) {
        case 1: 
            offset = Fix3_Node_Offset(1, 0, 0, t, tree);
            for (i = Bot1; i <= Top1; i++)
            {
                for (c=0; c<count; ++c) 
                {
                    srcPtr = rhs[c] + offset; 
                    ixs[c] = srcPtr[i];
                }

                xPtr = x + offset;
                yPtr = y + offset;
                fun.execute(&oxs[0], &ixs[0], count, *xPtr, *yPtr);

                for (c=0; c<count; ++c) 
                {
                    trgPtr = lhs[c] + offset; 
                    trgPtr[i] = oxs[c];
                }
            }
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Fix3_Node_Offset(2, i, 0, t, tree);
                for (j = Bot2[i]; j <= Top2[i]; j++)
                {
                    for (c=0; c<count; ++c) 
                    {
                        srcPtr = rhs[c] + offset; 
                        ixs[c] = srcPtr[j];
                    }

                    xPtr = x + offset;
                    yPtr = y + offset;
                    fun.execute(&oxs[0], &ixs[0], count, *xPtr, *yPtr);

                    for (c=0; c<count; ++c) 
                    {
                        trgPtr = lhs[c] + offset; 
                        trgPtr[j] = oxs[c];
                    }
                }
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Fix3_Node_Offset(3, i, j, t, tree);
                    for (k = Bot3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (c=0; c<count; ++c) 
                        {
                            srcPtr = rhs[c] + offset; 
                            ixs[c] = srcPtr[k];
                        }

                        xPtr = x + offset;
                        yPtr = y + offset;
                        fun.execute(&oxs[0], &ixs[0], count, *xPtr, *yPtr);

                        for (c=0; c<count; ++c) 
                        {
                            trgPtr = lhs[c] + offset; 
                            trgPtr[k] = oxs[c];
                        }
                    }
                }
            }
            break;
    } 
}








//------------------------------------------------------------------------------
// These should be in a separate file. Since new .cpp modules are combustible,
// in the interest of fire safety they are kept here.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Calculate (smooth) step function depending on the value of the index in 
// relation to the barrier. If 'ba' = 'A'bove, calculate H(x, barrier), 
// if 'ba' = 'B'elow, calculate 1- H(x, barrier). H(x, b) is a Heaviside function 
// defined as: H(x, b) = 0; x<b, 1 otherwise
//------------------------------------------------------------------------------
int Fix3_ExIndicator(
         double*                indicator,    /**<  Result slice           */
         double const*          index,        /**<  Observation index      */
         double                 barrier,      /**<  barrier (strike) value */
         char                   ba,           /**<  'B'elow or 'A'bove     */
         int                    smooth,       /**<  1 for smoothing        */
         int                    t,            /**<  Current time point     */
         FIX3_TREE_DATA const*  tree_data)    /**<  Tree data              */
{
    // -------------------------------------------------------------------------
    // Ex indicator function
    // -------------------------------------------------------------------------
    struct ExIndicator : public UnarySmoothOper
    {
        ExIndicator(double barrier, bool above, bool smooth) 
            : UnarySmoothOper(smooth), m_barrier(barrier), m_above(above) {}

        virtual void execute(double* lhs, double const* x, double const* step, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i) {
                if (m_above) // TODO - should have STRIKE_TOL (like BARRIER_TOL) - check the logic
                    lhs[i] = DrSmoothStep(1.0, 0.0, x[i], m_barrier * (1.0-STRIKE_TOL), step[i]);
                else
                    lhs[i] = DrSmoothStep(0.0, 1.0, x[i], m_barrier * (1.0+STRIKE_TOL), step[i]);
            }
        }
        double m_barrier;
        bool   m_above;
    };

    Fix3_SliceOper(tree_data, t, indicator, index, ExIndicator(barrier, ba=='A', smooth!=0));

    return SUCCESS;
}


//------------------------------------------------------------------------------
// Calculate (smooth) double step function depending on the value of the index 
// in relation to the high and low barriers. If 'io' = 'I'nside, calculate 
// H(x, low) - H(x, high), if 'io' = 'O'utside, calculate 
// 1- [H(x, low) + H(x, high)].  H(x, b) is a Heaviside function defined as:
// H(x, b) = 0; x<b, 1 otherwise
//------------------------------------------------------------------------------
int Fix3_KOIndicator(
         double*                indicator,    /**<  Slice to be modified   */
         double const*          index,        /**<  Observation index      */
         double                 low,          /**<  Low barrier            */
         double                 high,         /**<  High barrier           */
         char                   io,           /**<  'I'nside or 'O'utside  */
         int                    smooth,       /**<  1 for smoothing        */
         int                    t,            /**<  Current time point     */
         FIX3_TREE_DATA const*  tree_data)    /**<  Tree data              */
{
    // -------------------------------------------------------------------------
    // KO indicator function
    // -------------------------------------------------------------------------
    struct KOIndicator : public UnarySmoothOper
    {
        KOIndicator(double low, double high, bool inside, bool smooth) 
            : UnarySmoothOper(smooth), m_low(low), m_high(high), m_inside(inside) {}

        virtual void execute(double* lhs, double const* x, double const* step, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i) {
                if (m_inside) {
                    lhs[i] = DrSmoothStep(1.0, 0.0,    x[i], m_low  * (1.0-BARRIER_TOL), step[i]);
                    lhs[i] = DrSmoothStep(0.0, lhs[i], x[i], m_high * (1.0+BARRIER_TOL), step[i]);
                }
                else {
                    lhs[i] = DrSmoothStep(0.0, 1.0,    x[i], m_low  * (1.0+BARRIER_TOL), step[i]);
                    lhs[i] = DrSmoothStep(1.0, lhs[i], x[i], m_high * (1.0-BARRIER_TOL), step[i]);
                }
            }
        }
        double m_low;
        double m_high;
        bool   m_inside;
    };

    Fix3_SliceOper(tree_data, t, indicator, index, KOIndicator(low, high, io=='I', smooth!=0));

    return SUCCESS;
}





//------------------------------------------------------------------------------
// Calculate expectation given two values and the probability these values take.
// Probability is given as a slice, one value is a slice, the other is a scalar.
// Result is a slice. Of course, to cover all cases we should write 2*2*2=8
// functions where two values and the probability argument are either a scalar
// or a slice. If we used expression templates the compiler will do this for us. 
// Until then, write another overload when you need it.
// Calculates result = x * (1-p) + y * p
//------------------------------------------------------------------------------
int Fix3_Expectation(
         double*                result,   /**<  Result slice                */
         double const*          x,        /**<  Value of the random variable*/
         double                 y,        /**<  Value of the random variable*/
         double const*          p,        /**<  Probability slice           */
         int                    t,        /**<  Current time point          */
         FIX3_TREE_DATA const*  tree_data)/**<  Tree data                   */
{
    // -------------------------------------------------------------------------
    // Expectation function
    // -------------------------------------------------------------------------
    struct Expectation : public BinaryOper
    {
        Expectation(double y) : m_y(y) {}

        virtual void execute(double* lhs, double const* x, double const* p, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i) {
                assert(p[i] >= 0.0 && p[i] <= 1.0);
                lhs[i] = x[i] * (1.0-p[i]) + m_y * p[i];
            }
        }
        double m_y;
    };

    Fix3_SliceOper(tree_data, t, result, x, p, Expectation(y));

    return SUCCESS;
}




//------------------------------------------------------------------------------
// Calculate logical OR of two slices. The arguments are slices with values
// between 0 and 1. 
// Calculates: result = 1 - (1 - x) * (1 - y)
//------------------------------------------------------------------------------
int Fix3_LogicalOr(
         double*                result,   /**<  Result slice                */
         double const*          x,        /**<  Value of the random variable*/
         double const*          y,        /**<  Value of the random variable*/
         int                    t,        /**<  Current time point          */
         FIX3_TREE_DATA const*  tree_data)/**<  Tree data                   */
{
    // -------------------------------------------------------------------------
    // Logical OR function
    // -------------------------------------------------------------------------
    struct LogicalOr : public BinaryOper
    {
        virtual void execute(double* lhs, double const* x, double const* y, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i) {
                // asserts below are known to fail on the slice boundaries
                //assert(x[i] >= 0.0 && x[i] <= 1.0);
                //assert(y[i] >= 0.0 && y[i] <= 1.0);
                lhs[i] = 1.0 - (1.0 - x[i]) * (1.0 - y[i]);
            }
        }
    };

    Fix3_SliceOper(tree_data, t, result, x, y, LogicalOr());

    return SUCCESS;
}




int Fix3_AddScalar2(
         double*                result,       
         double const*          argument,     
         double                 scalar,       
         int                    t,            
         FIX3_TREE_DATA const*  tree_data)    
{
    // -------------------------------------------------------------------------
    // add scalar function: lhs = x + scalar
    // -------------------------------------------------------------------------
    struct AddScalar : public UnaryOper
    {
        AddScalar(double scalar)
            : m_scalar(scalar) {}

        virtual void execute(double* lhs, double const* x, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i)
                lhs[i] = x[i] + m_scalar;
        }
        double m_scalar;
    };

    Fix3_SliceOper(tree_data, t, result, argument, AddScalar(scalar));

    return SUCCESS;
}


int  Fix3_LadderInterpolate_t
             (double**              SVSlices,   /* (I/O) State variable slices */
              double const*         Step,       /* (I) Step for ladder         */
              double const*         Spread,     /* (I) Spread for ladder       */
              double                Sticky,     /* (I) Sticky coefficient      */
              double                FloorSt,    /* (I) Floor spread            */
              double                CapSt,      /* (I) Cap spread              */
              char                  AoM,        /* (I) Additive/multiplicative */
              double                RibFrac,    /* (I) Rib obs in range        */
              int                   NbStates,   /* (I) Nb of states            */
              double const*         CurrStates, /* (I) Curr state levels at i-1*/
              double const*         PrevStates, /* (I) Prev state levels at i  */
              int                   t,          /* (I) Current time point      */
              FIX3_TREE_DATA const* tree_data)  /* (I) Tree data structure     */
{
    // -------------------------------------------------------------------------
    // x-section interpolation function 
    // -------------------------------------------------------------------------
    struct Interpolator : public XSBinaryOper
    {
        Interpolator(double const* curr, double const* prev, double sticky, 
                     double cap, double floor, bool add, double rib)
            : m_curr(curr), m_prev(prev), m_sticky(sticky), 
              m_cap(cap), m_floor(floor), m_add(add), m_rib(rib) {}

        virtual void execute(double* lhs, double const* rhs, size_t cnt, double step, double spread) const
        {
            for (size_t i=0; i<cnt; ++i)
            {
                double X = m_sticky * m_curr[i] + spread;

                if (m_add)
                {
                    X += step;
                    X =  COLLAR(X, m_cap, m_floor);
                }
                else
                {
                    X =  COLLAR(X, m_cap, m_floor);
                    X *= step;
                }

                X *= m_rib;

                Fix3_DoubleQuadraticInterp(m_prev, rhs, cnt, X, &lhs[i]);
            }
        }
        double const*   m_curr;
        double const*   m_prev;
        double          m_sticky;
        double          m_cap;
        double          m_floor;
        bool            m_add;
        double          m_rib;
    };

    Fix3_SliceOper(tree_data, t, SVSlices, SVSlices, NbStates, Step, Spread,
            Interpolator(CurrStates, PrevStates, Sticky, CapSt, FloorSt, AoM=='A', RibFrac));

    return SUCCESS;
}


// This will replace static method in the product
int  XFix3_StickyLadderSwap_t
             (double    **Sticky,         /**< (O) Prices for all states   */
              long        ResetFlagSt,    /**< (I) Reset at timept         */
              long        ResetFlagFund,  /**< (I) Reset at timept         */
              double     *Step,           /**< (I) Step for ladder         */
              double     *Spread,         /**  (I) Spread for ladder       */
              char        AoM,            /**  (I) Additive or mult step   */
              double     *Funding,        /**< (I) Funding                 */
              double     *ZeroToPmtSt,    /**< (I) PmtZero 4 this reset    */
              double      FloorSt,        /**< (I) Floor spread            */
              double      CapSt,          /**< (I) Cap spread              */
              double      StickyCoeff,    /**< (I) Sticky coefficient      */
              double      RibFrac,        /**< (I) Rib obs in range        */
              double      OutsSt,         /**< (I) Outstanding             */
              double      DcfSt,          /**< (I) Day count fract         */
              char        SoZ,            /**< (I) Swap/zero coupon        */
              char        CompSt,         /**< (I) Simple/Compound pmt     */
              int         NbStates,       /**< (I) Nb of states            */
              double     *CurrStates,     /**< (I) Curr state levels at i-1*/
              double     *PrevStates,     /**< (I) Prev state levels at i  */
              int         t,              /**< (I) Current time point      */
              FIX3_TREE_DATA   *tree_data)     /**< (I) Tree data structure     */
{
    // -------------------------------------------------------------------------
    // Add sticky coupon to the bank of slices
    // -------------------------------------------------------------------------
    struct AddLadCoupon : public XSUnaryOper
    {
        AddLadCoupon(double const* prev, double dcf, double notional, bool isSimple, bool isZero)
            : m_prev(prev), m_dcf(dcf), m_notional(notional), m_isSimple(isSimple), m_isZero(isZero) {}

        virtual void execute(double* lhs, double const* rhs, size_t cnt, double zero) const
        {
            for (size_t i=0; i<cnt; ++i)
            {
                double R = ACC_FN(m_prev[i], m_dcf, m_isSimple);

                // Update current value: mulitply and add coupon
                lhs[i] = rhs[i];
                if (m_isZero)
                    lhs[i] *= (1.0 + R);
                lhs[i] += zero * R * m_notional;
            }
        }
        double const*   m_prev;
        double          m_dcf;
        double          m_notional;
        bool            m_isSimple;
        bool            m_isZero;
    };

    // -------------------------------------------------------------------------
    // Subtract funding coupon from the bank of slices
    // -------------------------------------------------------------------------
    struct SubFndCoupon : public XSUnaryOper
    {
        virtual void execute(double* lhs, double const* rhs, size_t cnt, double coupon) const
        {
            for (size_t i=0; i<cnt; ++i)
                lhs[i] = rhs[i] - coupon;
        }
    };

    if (ResetFlagSt)
    {
        // add coupon: stk = prev * dcf * zero * notional
        Fix3_SliceOper(tree_data, t, Sticky, Sticky, NbStates, ZeroToPmtSt, 
                AddLadCoupon(PrevStates, DcfSt, OutsSt, CompSt=='S', SoZ=='Z'));

        // interpolate
        Fix3_LadderInterpolate_t(
                Sticky, 
                Step, 
                Spread, 
                StickyCoeff,
                FloorSt, 
                CapSt,
                AoM,
                RibFrac,
                NbStates, 
                CurrStates,
                PrevStates,
                t, tree_data);
    }

    /* Add floating payoff if it's a reset BUT ONLY if it's ladder date */
    /* This restriction is achieved by setting flags in calc properly   */
    if (ResetFlagFund)
    {
        Fix3_SliceOper(tree_data, t, Sticky, Sticky, NbStates, Funding, SubFndCoupon());

        Fix3_Set_Slice(Funding, 0.0, t, tree_data);
    }

    return SUCCESS;
} 
