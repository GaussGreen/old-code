#include <fix123head.h>
#include <fix123oper.hpp>
#include <assert.h>
#include <vector>

//------------------------------------------------------------------------------
// Unary slice operation
//------------------------------------------------------------------------------
void
SliceOper(TREE_DATA const* tree, size_t t,
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
            offset = Node_Offset(1, 0, 0, t, tree) + Bot1;
            srcPtr = arg + offset; 
            trgPtr = lhs + offset;

            fun.execute(trgPtr, srcPtr, Top1 - Bot1 + 1);
//std::cout << "+LHS (" << Top1-Bot1+1 << ")" << *lhs << std::endl;
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Node_Offset(2, i, 0, t, tree) + Bot2[i];
                srcPtr = arg + offset; 
                trgPtr = lhs + offset;

                fun.execute(trgPtr, srcPtr, Top2[i] - Bot2[i] + 1);
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Node_Offset(3, i, j, t, tree) + Bot3[i][j];
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
SliceOper(TREE_DATA const* tree, size_t t,
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
            offset = Node_Offset(1, 0, 0, t, tree) + Bot1;
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
                        steps[idx++] = GetIndexStep(refPtr, 1, i, 0, 0, t, tree);

                fun.execute(trgPtr, srcPtr, &steps[0], len);
            }
            break;

        case 2:
            for (i = Bot1; i <= Top1; i++) {
                offset = Node_Offset(2, i, 0, t, tree) + Bot2[i];
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
                            steps[idx++] = GetIndexStep(refPtr, 2, i, j, 0, t, tree);

                    fun.execute(trgPtr, srcPtr, &steps[0], len);
                }
            }

            break;
        case 3:
            for (i = Bot1; i <= Top1; i ++) {
                for (j = Bot2[i]; j <= Top2[i]; j++) {
                    offset = Node_Offset(3, i, j, t, tree) + Bot3[i][j];
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
                                steps[idx++] = GetIndexStep(refPtr, 3, i, j, k, t, tree);

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
SliceOper(TREE_DATA const* tree, size_t t,
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
            offset = Node_Offset(1, 0, 0, t, tree) + Bot1;
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
                offset = Node_Offset(2, i, 0, t, tree) + Bot2[i];
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
                    offset = Node_Offset(3, i, j, t, tree) + Bot3[i][j];
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
// These should be in a separate file. Since new .cpp modules are combustible,
// in the interest of fire safety they are kept here.
//------------------------------------------------------------------------------
int KOIndicator(
         double*                indicator,    /**< (O) Slice to be modified   */
         double const*          index,        /**< (I) observation index      */
         double                 low,          /**< (I) Low barrier            */
         double                 high,         /**< (I) High barrier           */
         char                   io,           /**< (I) 'I'nside or 'O'utside  */
         int                    t,            /**< (I) Current time point     */
         TREE_DATA const*  tree_data)    /**< (I) Tree data              */
{
    // -------------------------------------------------------------------------
    // KO indicator function
    // -------------------------------------------------------------------------
    struct KOIndicator : public UnaryOper
    {
        KOIndicator(double low, double high, bool inside) 
            : m_low(low), m_high(high), m_inside(inside) {}

        virtual void execute(double* lhs, double const* x, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i)
            {
                if (m_inside)
                    lhs[i] = (x[i] > m_low  * (1.0-BARRIER_TOL) && x[i] < m_high * (1.0+BARRIER_TOL));
                else
                    lhs[i] = (x[i] < m_low  * (1.0+BARRIER_TOL) || x[i] > m_high * (1.0-BARRIER_TOL));
            }
        }
        double m_low;
        double m_high;
        bool   m_inside;
    };

    SliceOper(tree_data, t, indicator, index, KOIndicator(low, high, io=='I'));

    return SUCCESS;
}


int AddScalar2(
         double*                result,       /**< (O) Slice to be modified   */
         double const*          argument,     /**< (I) observation index      */
         double                 scalar,       /**< (I) Low barrier            */
         int                    t,            /**< (I) Current time point     */
         TREE_DATA const*  tree_data)    /**< (I) Tree data              */
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

    SliceOper(tree_data, t, result, argument, AddScalar(scalar));

    return SUCCESS;
}


