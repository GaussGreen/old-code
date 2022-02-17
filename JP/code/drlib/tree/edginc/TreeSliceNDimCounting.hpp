//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TreeSliceNDimCounting.hpp
//
//   Author      : Charles Morcom, July 31, 2006
//
//   Description : A general n-dimensional slice with rectangular limits
//                 and with discrete counting states.
//                 
//                 This is designed to be used by CountingTree and successor
//                 classes. While this compiles, it is not yet fully usable,
//                 and will require some debugging. As an early prototype,
//                 you will also notice that things are not yet very carefully
//                 optimized.
//----------------------------------------------------------------------------

#ifndef QR_TREESLICENDIMCOUNTING_HPP
#define QR_TREESLICENDIMCOUNTING_HPP

// this header file should be included by TreeSlice.hpp only
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


class TREE_DLL TreeSliceNDimCounting : public TreeSlice { // if the compiler complains about TreeSlice undefined it is because TreeSliceNDimCounting.hpp has been included where TreeSlice.hpp should have been.
public:

    /* Class defining tree limits. It might seem more natural to have this in CountingTree,
       but it has to be here since TreeSlice causes this to be included before CountingTree,
       and the template eval requires that things from CountingTree needed here have to be
       more than just forward declared. */
    class TREE_DLL Range {
    public:
        /*Counting state and diffusion grid state limits at each time point*/
	    vector< vector<int> > cLInnerLimits;
	    vector< vector<int> > cUInnerLimits;
	    vector< vector<int> > dLInnerLimits;
	    vector< vector<int> > dUInnerLimits;
        /* These are the limits for the (cPrime, dStar) post-c-transition states */
        vector< vector<int> > cLMidLimits;
	    vector< vector<int> > cUMidLimits;
	    vector< vector<int> > dLMidLimits;
	    vector< vector<int> > dUMidLimits;
        /* These are the limits for the (cPrime, dStar) post-c-transition states */
        vector< vector<int> > cLOuterLimits;
	    vector< vector<int> > cUOuterLimits;
	    vector< vector<int> > dLOuterLimits;
	    vector< vector<int> > dUOuterLimits;
        /* These are the max (for U) and min (for L) of the limits across time steps*/
        vector<int> dLGlobalLimits;
        vector<int> dUGlobalLimits;
        vector<int> cLGlobalLimits;
        vector<int> cUGlobalLimits;
        /* These are the global sizes */
        vector<int> dSize; //dUGlobalLimits-dLGlobalLimits+1
        vector<int> cSize; //cUGlobalLimits-cLGlobalLimits+1
        int tIdx; // current time-step
        bool hasCStates;

    };

    /*************************** variables ***************************/
protected:
    const Range& range;       // the limits of the tree that this slice is associated with
	int dim;                  // d-dimension of this slice (1/2 total dim if c-states)
    double *values;           // flat array containing the values
    bool deleteWhenDestroyed; // responsible to delete the "values" array

    string curveToDEV;

    /*************************** methods ***************************/
public:

    /******** [begin] must-implement for TreeSlice ***/
	virtual const string typeName() const { return "TreeSliceNDimCounting"; }
    virtual TreeSliceSP clone( bool copyValues = true ) const;
    virtual TreeSlice& operator=(double constant);
    virtual TreeSlice& operator=(const TreeSlice&);
    template<class A> 
        void eval(const A & expr);
    TreeSliceSP calcSmoothStep() const;
    
	// !!! TO BE REMOVED: ONLY HERE BECAUSE IN INTERFACE
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const
    {
        if( dim == 1 )
        {
            //bot = range.globalMinLimits[0];
            //top = range.globalMaxLimits[0];
        }
        else
        {
            bot = 0;
            top = valuesSize() - 1;
        }
    }
    virtual double * getValues() const
    {
        return values; //dim == 1 ? values - range.globalMinLimits[0] : values;
    }

	/******** [end] must-implement for TreeSlice ***/



    TreeSliceNDimCounting(const Range& range, const string& curveToDEV="");
    ~TreeSliceNDimCounting();

    TreeSliceNDimCounting& operator=(const TreeSliceNDimCounting& s);

	double* getValuePtr() { return values; }

    // prepare the slice to receive data in dimension newDim
    void allocDim(int newDim); 

    // gets the value at the centre of the slice 
    // e.g. coord (0) or (0,0) or (0,0,0), etc. depending on the dimension
	// note that, for counting processes this may not be very interesting
    double getCentre() const;

    int getDim() { return dim; }
    void expand(int newDim);

protected:

     // computes the size of the "values" array
    int valuesSize() const;

    // alloc n doubles for "values". if n<0 no alloc is done. See also allocDim
    void resize(int n); 

    // copy the content and takes the ownership of an existing slice
    void takeFrom(TreeSliceNDimCounting &s);

    /* I have to provide these for TreeSlice, but this is not how
       these should be used. */
    // calc the offet to add to "values" to access an element in a 3D slice
    int offset( int i, int j, int k ) const;
    double val(int i, int j, int k) const {return values[offset(i,j,k)];}
    // calc the offet to add to "values" to access an element in a 2D slice
    int offset( int i, int j ) const;
    double val(int i, int j) const {return values[offset(i,j)];}
    // calc the offet to add to "values" to access an element in a 1D slice
    int offset( int i ) const;
    double val(int i) const {return values[offset(i)];}

    /* These are new offset methods which include the counting-states */
    int offset(const vector<int>& dCoords) const;
    double val(const vector<int>& dCoords) const {return values[offset(dCoords)];}
    int offset(
        const vector<int>& cCoords, // counting-state coordinates
        const vector<int>& dCoords  // diffusion-state coordinates
        ) const;
    double val(
        const vector<int>& cCoords, 
        const vector<int>& dCoords) const {
            return values[offset(cCoords, dCoords)];}


    // see usage of iterSeek(...) in calc() below: sets TreeSlice::iter pointer
    void iterSeek(int i) const;             // get ready to loop on dimension #0 (1st)
    void iterSeek(int i, int j) const;        // get ready to loop on dimension #1 (2nd)
    void iterSeek(int i, int j, int k) const; // get ready to loop on dimension #2 (3rd)

    void iterSeek(const vector<int>& dCoords) const; // n-dimensions
    void iterSeek(
        const vector<int>& cCoords, 
        const vector<int>& dCoords) const; // n-dims with c-states

    double* iterInc(void) const { return iter;}
	void testTreeStep() const;

private:
    TreeSliceNDimCounting(const TreeSliceNDimCounting&);

    friend class CountingTree;
};
typedef refCountPtr< TreeSliceNDimCounting > TreeSliceNDimCountingSP;


/***************** body of the calc function (template) ******************/
template<class A>
void TreeSliceNDimCounting::eval(const A & expr) {
try {
    //const int N = expr.sliceCount();
    const char* (*print)(void*) = printExpr<A>; 
    (void)print; // to remove compiler warnings
    // watch "print((void*)expr)" and "*this->iter" in your debugger

	double** iterators = (double**)::alloca( expr.sliceCount * sizeof(double*) );

    // get a list of the slices used in the expression and convert them to TreeSliceNDimCounting
	const TreeSliceNDimCounting** s = 
        (const TreeSliceNDimCounting**)::alloca( (expr.sliceCount+1) * sizeof(const TreeSliceNDimCounting*) );
    {
        const TreeSliceNDimCounting** end = expr.listSlices(s);
        ASSERT((end-s)==expr.sliceCount);
        for (int i=0; i<expr.sliceCount; ++i) {
			if (!s[i]) {
				throw ModelException("Cannot evaluate an expression using arguments of type "+
					s[i]->typeName());
			}
        }
        s[expr.sliceCount]=this;
    }


    // calc the max dimension of the slices of the expression
    int newDim = 0;
    for (int i=0; i<expr.sliceCount; ++i) {
        newDim = Maths::max(newDim, s[i]->dim);
    }

	// check treeStep of the arguments
	for (int i=0; i<expr.sliceCount; i++) {
		s[i]->testTreeStep();
	}

    // if the output dimension is wrong alloc a new slice for the result
    TreeSliceNDimCounting newSlice(range);
    if (newDim != dim) {
        newSlice.allocDim(newDim);
        newSlice.copyTreeSlice(*this);
    } else {
        newSlice.takeFrom(*this);
    }

    /* NOTE THAT, FOR SLICE EVALUATIONS, YOU ALWAYS WANT THE INNER LIMITS */
    int t = range.tIdx;
    const vector<int>& dLLimits = range.dLInnerLimits[t];
    const vector<int>& dULimits = range.dUInnerLimits[t];
    const vector<int>& cLLimits = range.cLInnerLimits[t];
    const vector<int>& cULimits = range.cUInnerLimits[t];
 
    // loop depending on the output dim
    int vs = newSlice.valuesSize();
    switch (newSlice.dim) {
        case 0:
            // zero dimension is single point
            for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(0);
            newSlice.iterSeek(0);
            ASSERT( (newSlice.iter - newSlice.values) < vs );
            *newSlice.iter = expr.calc();
            break;

        case 1: 
            if (range.hasCStates) {
                /* SINGLE DIMENSION WITH COUNTING STATE */
                static vector<int> c(1); // non-reentrant!
                static vector<int> d(1);
                for (c[0]=cLLimits[0]; c[0]<=cULimits[0]; c[0]++) {
                    int iStart = dLLimits[0];
                    int iStop  = dULimits[0];
                    d[0] = iStart;
                    for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(c,d);
                    newSlice.iterSeek(c,d);

                    // i <=> newSlice.iter - stop + iStop
                    double *stop = newSlice.iter + iStop - iStart;
                    while (newSlice.iter <= stop) { // innermost "for" replaced by while for performance
                        ASSERT( (newSlice.iter - newSlice.values) < vs );
                        *newSlice.iter = expr.calc();
                        // two step increment to avoid incrementing twice if a slice appears twice in the expression
                        for (int n=0; n<expr.sliceCount; ++n) iterators[n] = s[n]->iterInc();
                        for (int n=0; n<expr.sliceCount; ++n) s[n]->iter = iterators[n];
                        ++newSlice.iter;
                    }
                }
            } else {
                /* SINGLE DIFFUSION DIMENSION; NO COUNTING */
                int iStart = dLLimits[0];
                int iStop  = dULimits[0];
                for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(iStart);
                newSlice.iterSeek(iStart);

                // i <=> newSlice.iter - stop + iStop
                double *stop = newSlice.iter + iStop - iStart;
                while (newSlice.iter <= stop) { // innermost "for" replaced by while for performance
                    ASSERT( (newSlice.iter - newSlice.values) < vs );
                    *newSlice.iter = expr.calc();
                    // two step increment to avoid incrementing twice if a slice appears twice in the expression
                    for (int n=0; n<expr.sliceCount; ++n) iterators[n] = s[n]->iterInc();
                    for (int n=0; n<expr.sliceCount; ++n) s[n]->iter = iterators[n];
                    ++newSlice.iter;
                }
            }
            break;
        case 2:
            if (range.hasCStates) {
                /* 2-DIM DIFFUSION; AT LEAST ONE HAS COUNTING STATE */
                static vector<int> d(2);
                static vector<int> c(2);
                for (c[0]=cLLimits[0]; c[0]<=cULimits[0]; c[0]++) {
                    for (c[1]=cLLimits[1]; c[1]<=cULimits[1]; c[1]++) {
                        int t0 = dULimits[0];
                        for (int i = dLLimits[0]; i <= t0; i++) {
                            int jStart = dLLimits[1];
                            int jStop  = dULimits[1];
                            d[0] = i;
                            d[1] = jStart;
                            for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(c,d);
                            newSlice.iterSeek(c,d);

                            // j <=> newSlice.iter - stop + jStop
                            double *stop = newSlice.iter + jStop - jStart;
                            while (newSlice.iter <= stop) // innermost "for" replaced by while for performance
                            {
                                ASSERT( (newSlice.iter - newSlice.values) < vs );
                                *newSlice.iter = expr.calc();
                                // two step increment to avoid incrementing twice if a slice appears twice in the expression
                                for (int n=0; n<expr.sliceCount; ++n) iterators[n] = s[n]->iterInc();
                                for (int n=0; n<expr.sliceCount; ++n) s[n]->iter = iterators[n];
                                ++newSlice.iter;
                            }
                        } 
                    }
                }
                
            } else {
                /* 2-DIM DIFFUSION; NEITHER HAS COUNTING STATE */
                int t0 = dULimits[0];
                for (int i = dLLimits[0]; i <= t0; i++) 
                {
                    int jStart = dLLimits[1];
                    int jStop  = dULimits[1];
                    for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(i, jStart);
                    newSlice.iterSeek(i, jStart);

                    // j <=> newSlice.iter - stop + jStop
                    double *stop = newSlice.iter + jStop - jStart;
                    while (newSlice.iter <= stop) // innermost "for" replaced by while for performance
                    {
                        ASSERT( (newSlice.iter - newSlice.values) < vs );
                        *newSlice.iter = expr.calc();
                        // two step increment to avoid incrementing twice if a slice appears twice in the expression
                        for (int n=0; n<expr.sliceCount; ++n) iterators[n] = s[n]->iterInc();
                        for (int n=0; n<expr.sliceCount; ++n) s[n]->iter = iterators[n];
                        ++newSlice.iter;
                    }
                }
            }
            break;
        case 3:
            if (range.hasCStates) {
                throw ModelException("Invalid dimension: I haven't written eval for dimensions>2 with c-states, yet!");
            } else {
                int t0 = dULimits[0];
                for (int i = dLLimits[0]; i <= t0; i++) 
                {
                    int t1 = dULimits[1];
                    for (int j = dLLimits[1]; j <= t1; j++) 
                    {
                        int kStart = dULimits[2];
                        int kStop  = dLLimits[2];
                        for (int n=0; n<expr.sliceCount; ++n) s[n]->iterSeek(i, j, kStart);
                        newSlice.iterSeek(i, j, kStart);

                        // k <=> newSlice.iter - stop + kStop
                        double *stop = newSlice.iter + kStop - kStart;
                        while (newSlice.iter <= stop) // innermost "for" replaced by while for performance
                        {
                            ASSERT( (newSlice.iter - newSlice.values) < vs );
                            *newSlice.iter = expr.calc();
                            // two step increment to avoid incrementing twice if a slice appears twice in the expression
                            for (int n=0; n<expr.sliceCount; ++n) iterators[n] = s[n]->iterInc();
                            for (int n=0; n<expr.sliceCount; ++n) s[n]->iter = iterators[n];
                            ++newSlice.iter;
                        }
                    }
                }
                break;
            }
        default: 
			throw ModelException("Invalid dimension: I haven't written eval for dimensions>3, yet!");
    }
    // reference the value of newSlice
    takeFrom(newSlice);
} catch (exception& e) {
	char buf[10000];
	*buf=0;
	expr.printDebug(buf);
	throw ModelException(e, "TreeSliceNDimCounting::eval(const A& expr), "+string(buf));
}
}

DRLIB_END_NAMESPACE

#endif
