#ifndef TREE_SLICE_RATES_HPP
#define TREE_SLICE_RATES_HPP

// this header file should be included by TreeSlice.hpp only

#include "edginc/Maths.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/RatesSliceRange.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL TreeSliceRates : public TreeSlice { // if the compiler complains about TreeSlice undefined it is because TreeSliceRates.hpp has been included where TreeSlice.hpp should have been.
public:

    /*************************** types ***************************/
    class TREE_DLL Range : public RatesSliceRange {
    public:
        Range(int bot1, int top1, // 1st dimension range
              int bot2, int top2, // 2st dimension range
              int bot3, int top3) // 3st dimension range
        {
            init(bot1, top1, bot2, top2, bot3, top3);
        }

        Range(int bot1, int top1, // 1st dimension range
              int bot2, int top2) // 2st dimension range
        {
            init(bot1, top1, bot2, top2, 0, 0);
        }

        Range(int bot1, int top1) // 1st dimension range
        {
            init(bot1, top1, 0, 0, 0, 0);
        }

    private:

        void init(int bot1, int top1,          // 1st dimension range
                  int bot2, int top2,          // 2st dimension range
                  int bot3, int top3,          //
                  int bot4 = 0, int top4 = 0); // 3st dimension range
        
    };

    DECLARE_REF_COUNT(Range);

    string curveToDEV;
    int curveToDEVIdx;

    virtual int getCurveToDEVIdx() const { return curveToDEVIdx; }

    /*************************** variables ***************************/
protected:
    const Range &range;
    int dim;
    double *values;   // flat array containing the values
    bool deleteWhenDestroyed; // responsible to delete the "values" array

    mutable int incInner; // value to add to "iter" from the inner loop

    /*************************** methods ***************************/
public:

    /******** [begin] must-implement for TreeSlice ***/
    virtual void printDetails(char *s) const;
    virtual const string typeName() const { return "TreeSliceRates"; }
    virtual TreeSliceSP clone( bool copyValues = true ) const;
    virtual TreeSlice& operator=(double constant);
    virtual TreeSlice& operator=(const TreeSlice&);

    virtual double getCentre() const;
    virtual bool isZero() const { return (dim==0) && Maths::isZero(*values); }
    virtual void clear() {allocDim(-1);}

    TreeSliceSP calcSmoothStep() const;

    template<class A> 
    void eval(const A & expr);

    // do an operation on slices. Output slices come first
    template<class T> // T must define sliceCount, compute() and printDebug()
    static void loopOnSlices(const T & oper, TreeSliceRates **s, int nbOutput);

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const;
    virtual double * getValues() const;

    // to be removed
    virtual string getCurveToDEV() const { return curveToDEV; }

    /******** [end] must-implement for TreeSlice ***/


    TreeSliceRates(const RatesSliceRange &range, const string& curveToDEV, int curveToDEVIdx);
    ~TreeSliceRates();

    TreeSliceRates& operator=(const TreeSliceRates& s);

    double* getValuePtr() { return values; }
    void    allocDim(int newDim); // prepare the slice to receive data in dimension newDim
    int     getDim() { return dim; }
    void    expand(int newDim);
    void    testTreeStep() const;

protected:

     // computes the size of the "values" array
    int valuesSize() const;

    // alloc n doubles for "values". if n<0 no alloc is done. See also allocDim
    void resize(int n); 

    // copy the content and takes the ownership of an existing slice
    void takeFrom(TreeSliceRates &s);

    // calc the offet to add to "values" to access an element in a 3D slice
    int offset( int i, int j, int k ) const;
    double val(int i, int j, int k) const {return values[offset(i,j,k)];}
    // calc the offet to add to "values" to access an element in a 2D slice
    int offset( int i, int j ) const;
    double val(int i, int j) const {return values[offset(i,j)];}
    // calc the offet to add to "values" to access an element in a 1D slice
    int offset( int i ) const;
    double val(int i) const {return values[offset(i)];}

    // see usage of iterSeek(...) in calc() below
    void iterSeek(int i) const;             // get ready to loop on dimension #0 (1st)
    void iterSeek(int i, int j) const;        // get ready to loop on dimension #1 (2nd)
    void iterSeek(int i, int j, int k) const; // get ready to loop on dimension #2 (3rd)

    double* iterInc(void) const { return iter+incInner;}

private:
    
    template<class A>
    class ExprOper;

    template<class T> // T must defile sliceCount and compute() and printDebug()
    void loopOnSlicesImpl(const T & oper, TreeSliceRates **s, int dim);

    void prepareSlicesForEval(int nbSlices, TreeSliceRates **s, TreeSliceRates &newSlice);

    template<class T> // T must defile sliceCount
    static void incAllSlices(const T & oper, double** tmpIter, TreeSliceRates **s);

    TreeSliceRates(const TreeSliceRates&);
};
typedef refCountPtr< TreeSliceRates > TreeSliceRatesSP;

/*****************************************************************************/
/************************** template functions definition ********************/
/*****************************************************************************/

// An operation object that computes an expression. For eval() to use on loopOnSlicesImpl().
template<class A>
class TreeSliceRates::ExprOper {
    TreeSliceRates & slice;
    const A & expr;
public:
    static const int sliceCount = A::sliceCount+1; // +1 for newSlice

    inline void compute() const {
        *slice.iter = expr.calc();
    }
    void printDebug(char *s) const { expr.printDebug(s); }
    ExprOper(const A & expr, TreeSliceRates & slice) : slice(slice), expr(expr)
    {}
};

template<class A>
void TreeSliceRates::eval(const A & expr) {
    try {
        const char* (*print)(void*) = printExpr<A>; 
        (void)print; // to remove compiler warnings
        // watch "print((void*)&expr)" and "*this->iter" in your debugger

        TreeSliceRates newSlice(range, curveToDEV, curveToDEVIdx);
        ExprOper<A> exprOper(expr, newSlice);

        const TreeSliceRates** s = (const TreeSliceRates**)::alloca( (expr.sliceCount+1) * sizeof(const TreeSliceRates*) );
        const TreeSliceRates** end = expr.listSlices(s);
        ASSERT((end-s)==expr.sliceCount);

        prepareSlicesForEval(expr.sliceCount, const_cast<TreeSliceRates**>(s), newSlice);

        loopOnSlicesImpl(exprOper, const_cast<TreeSliceRates**>(s), newSlice.dim);

        takeFrom(newSlice); // take the values of newSlice
        treeStep = range.treeStep;
    }
    catch (exception& e) {
        char buf[10000];
        *buf=0;
        expr.printDebug(buf);
        throw ModelException(e, "TreeSliceRates::eval(const A & expr), "+string(buf));
    }
}

template<class T>
void TreeSliceRates::loopOnSlicesImpl(const T & oper, TreeSliceRates **s, int dim) {
    const char* (*print)(void*) = printExpr<T>; (void)print; // for debug print((void*)&oper)
#if !defined(LINUX)
    double** tmpIter = (double**)::alloca( oper.sliceCount * sizeof(double*) );
#else
    double* tmpIter[oper.sliceCount]; // fix gcc-3.3 segfault
#endif

    // loop depending on the output dim
    switch (dim) {
        case 0:
            for (int n=0; n<oper.sliceCount; ++n) {
                s[n]->iterSeek(0);
            }
            oper.compute();
            break;

        case 1: {
            int iStart = range.limits.bot1;
            int nbLoops = range.limits.top1 - iStart + 1;
            for (int n=0; n<oper.sliceCount; ++n) {
                s[n]->iterSeek(iStart);
            }
            while (nbLoops--) { // "while" instead of "for" for efficiency
                oper.compute();
                incAllSlices(oper, tmpIter, s);
            }
            break;
        }
        case 2: {
            int t0 = range.limits.top1;
            for (int i = range.limits.bot1; i <= t0; i++) 
            {
                int jStart = range.limits.bot2[i];
                int nbLoops = range.limits.top2[i] - jStart + 1;
                for (int n=0; n<oper.sliceCount; ++n) {
                    s[n]->iterSeek(i, jStart);
                }
                while (nbLoops--) { // "while" instead of "for" for efficiency
                    oper.compute();
                    incAllSlices(oper, tmpIter, s);
                }
            }
            break;
        }
        case 3: {
            int t0 = range.limits.top1;
            for (int i = range.limits.bot1; i <= t0; i++) 
            {
                int t1 = range.limits.top2[i];
                for (int j = range.limits.bot2[i]; j <= t1; j++) 
                {
                    int kStart = range.limits.bot3[i][j];
                    int nbLoops = range.limits.top3[i][j] - kStart + 1;
                    for (int n=0; n<oper.sliceCount; ++n) {
                        s[n]->iterSeek(i, j, kStart);
                    }
                    while (nbLoops--) { // "while" instead of "for" for efficiency
                        oper.compute();
                        incAllSlices(oper, tmpIter, s);
                    }
                }
            }
            break;
        }
        default: 
            throw ModelException("Invalid dimension");
    }
}

template<class T>
void TreeSliceRates::loopOnSlices(const T & oper, TreeSliceRates **s, int nbOutput) {

    // find max dimension on input slices
    int dim=-1;
    for (int n = nbOutput; n<oper.sliceCount; ++n) {
        if (s[n]->dim<0) {
            throw ModelException(
                "Cannot evaluate an expression using uninitialized argument \""
                + s[n]->name+"\" (index "+Format::toString(n)+")");
        }
        dim = Maths::max(dim, s[n]->dim);
    }
    // set dim on output slices
    for (int n = 0; n<nbOutput; ++n) {
        s[n]->allocDim(dim);
    }
    s[0]->loopOnSlicesImpl(oper, s, dim);
}

template<class T>
void TreeSliceRates::incAllSlices(const T & oper, double** tmpIter, TreeSliceRates **s) {
    // two step increment to avoid incrementing twice if a slice appears twice in the expression
    for (int n=0; n<oper.sliceCount; ++n) {
        tmpIter[n] = s[n]->iterInc();
    }
    for (int n=0; n<oper.sliceCount; ++n) {
        s[n]->iter = tmpIter[n];
    }
}
DRLIB_END_NAMESPACE

#endif
