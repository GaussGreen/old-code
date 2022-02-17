//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSlice.hpp
//
//   Description : Generic tree slice.
//
//   Date        : Apr 11, 2006
//
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_HPP
#define TREE_SLICE_HPP

#include "edginc/ModelException.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/TreeSliceExpr.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include <string.h>
#ifdef _MSC_VER
#include <malloc.h>
#else
#include <alloca.h>
#endif

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(TreeSlice);

class TREE_DLL TreeSlice : public SliceMarker< TreeSlice >
{
    /****************** public variables ******************/
public:
    string name;        // for debug only (meant to stay)
    double * theValues; // for debug purposes (meant to be removed)
    int treeStep;       // to help check slice consistency during operations

    mutable double *iter;

protected:
    template< typename T >
    class LoopList
    {
        int count;
        T * list;
    public:
        LoopList() : count( 0 ), list( 0 ) {}
        ~LoopList() { if( list ) ::free( list ); }
        operator T *() const { return list; }
        T * reserve( int newCount )
        {
            if( count < newCount )
            {
                list = (T *)::realloc( list, newCount * sizeof( T ) );
                count = newCount;
            }
            return list;
        }
    };

    /****************** private variables ******************/
private:
    mutable const TreeSlice *exposedSlice;

    /****************** public methods ******************/
public:
    TreeSlice();
    TreeSlice(const TreeSlice &s) ;
    virtual ~TreeSlice() {}

    // do an operation on slices
    template<class T> // T must define nbSlices, compute(), printDebug(), listInputSlices(), listOutputSlices()
    static void loopOnSlices( const T & oper );

    // resizes a vector filling the gaps with cloned proto slice
    static void resizeVector(
        int newSize,
        const TreeSlice & proto,
        vector< TreeSliceSP > & slices,
        const string &namePrefix,
        bool copyValues = false );

    // expression evaluation
    template< typename ARG >
    void evalMultiType( const SliceMarker< ARG > & arg ); // does the assignment

    template< typename ARG >
    TreeSlice & operator =( const SliceMarker< ARG > & arg )
    {
        evalMultiType( arg );
        return *this;
    }
    template< typename ARG >
    TreeSlice & operator +=( const ARG & arg ) { return *this = *this + arg; }
    template< typename ARG >
    TreeSlice & operator -=( const ARG & arg ) { return *this = *this - arg; }
    template< typename ARG >
    TreeSlice & operator *=( const ARG & arg ) { return *this = *this * arg; }
    template< typename ARG >
    TreeSlice & operator /=( const ARG & arg ) { return *this = *this / arg; }

    /************ expression template primitives ************/
public:
    static const int sliceCount = 1;
    template< typename S >
    const S** listSlices(const S** l) const {*l=dynamic_cast<const S*>(exposedSlice); return ++l;}
    inline double calc() const {return *exposedSlice->iter;}
    void printDebug(char *s) const;

    /****************** virtual methods ******************/
public:
    virtual void printDetails(char *s) const; // to be called by the debugger
    virtual const string typeName() const = 0;
    virtual TreeSliceSP clone( bool copyValues = true ) const = 0;
    virtual TreeSlice& operator=(double x) = 0;
    virtual TreeSlice& operator=(const TreeSlice& s) {
        evalMultiType(s);   
        return *this;
    }

    // gets the value at the centre of the slice 
    // e.g. coord (0) or (0,0) or (0,0,0) depending on the dimension
    virtual double getCentre() const;

    // Returns true if the slice contains all zero values.
    // A return value of false does not imply that one of the values is not zero.
    // The aim of this function is to avoid unnecessary computations when possible.
    virtual bool isZero() const { return false; }

    // makes the slice invalid, as it was at its construction.
    virtual void clear() const {}

    virtual TreeSliceSP calcSmoothStep() const;

    virtual int getCurveToDEVIdx() const { return 0; }

    // Other methods needed:
    // template< typename A > void eval( const A & expr );
    // template< typename T > void loopOnSlices( const T & oper, TreeSliceType ** slices, int nbOutput );

    // !!! TO BE REMOVED
    // direct access to values (for legacy purposes only)
    virtual void getCalcRange( int & bot, int & top ) const = 0;
    virtual double * getValues() const = 0;

    // to be removed
    virtual string getCurveToDEV() const { return ""; }

    /********************* protected  *********************/
protected:
    class TREE_DLL ExposedSliceOverride
    {
        const TreeSlice * slice;
        const TreeSlice * oldValue;
        void clear();
    public:
        void set( const TreeSlice * source, const TreeSlice * target );
        ExposedSliceOverride() : slice( 0 ), oldValue( 0 ) {}
        ~ExposedSliceOverride() { clear(); }
    };

    void copyTreeSlice(const TreeSlice &s); // copies the fields of TreeSlice only
};

namespace
{

// helper function to perform swap
template< typename T >
inline
static void swapT( T & t1, T & t2 )
{
    T t = t1;
    t1 = t2;
    t2 = t;
}

}

DRLIB_END_NAMESPACE

/*****************************************************************************/
/************************** template functions definition ********************/
/*****************************************************************************/

#include "edginc/TreeSliceRates.hpp"
#include "edginc/TreeSliceRatesCompact.hpp"
#include "edginc/TreeSliceNDimCounting.hpp"
#include "edginc/TreeSliceBasic.hpp"
#include "edginc/TreeSliceGeneral.hpp"
#include "edginc/TreeSliceEQ.hpp"
#include "edginc/TreeSliceLayer.hpp"

#include "edginc/TreeSliceOper.hpp"

DRLIB_BEGIN_NAMESPACE

template< typename T >
void TreeSlice::loopOnSlices( const T & oper )
{
    if( ! oper.sliceCount ) {
        oper.compute();
        return;
    }

    // get a list of the slices
    TreeSlice ** slices = (TreeSlice **)::alloca( oper.sliceCount * sizeof( TreeSlice * ) );
    TreeSlice ** end = oper.listOutputSlices( slices );
    const int nbOutput = end - slices;
    end = const_cast< TreeSlice ** >( oper.listInputSlices( const_cast< const TreeSlice ** >( end ) ) );
#ifdef DEBUG
    ASSERT( end - slices == oper.sliceCount );
    for( int n = 0; n < oper.sliceCount; ++n )
        ASSERT( slices[ n ] );
#endif

    TreeSliceRates* treeSliceRates = dynamic_cast<TreeSliceRates*>(*slices);
    if (treeSliceRates) {
        treeSliceRates->loopOnSlices(oper, reinterpret_cast< TreeSliceRates ** >( slices ), nbOutput);
        return;
    }

	TreeSliceRatesCompact* treeSliceRatesCompact = dynamic_cast<TreeSliceRatesCompact*>(*slices);
    if (treeSliceRatesCompact) {
        treeSliceRatesCompact->loopOnSlices(oper, reinterpret_cast< TreeSliceRatesCompact ** >( slices ), nbOutput);
        return;
    }

    TreeSliceEQ* treeSliceEQ = dynamic_cast<TreeSliceEQ*>(*slices);
    if (treeSliceEQ) {
        treeSliceEQ->loopOnSlices(oper, reinterpret_cast< TreeSliceEQ ** >( slices ), nbOutput);
        return;
    }
/*
    TreeSliceLayer* treeSliceLayer = dynamic_cast<TreeSliceLayer*>(*slices);
    if (treeSliceLayer) {
        treeSliceLayer->loopOnSlices(oper, slices, nbOutput);
        return;
    }*/

    string type = slices[0] ? slices[0]->typeName() : "";
    throw ModelException( "TreeSlice::loopOnSlices", "Slice " + type + " not handled" );
}

template< typename ARG >
void TreeSlice::evalMultiType( const SliceMarker< ARG > & arg )
{
    const ARG & expr = static_cast< const ARG & >( arg );

    TreeSliceLayer* treeSliceLayer = dynamic_cast<TreeSliceLayer*>(this);
    if (treeSliceLayer) {
        treeSliceLayer->eval(expr);
        return;
    }

    TreeSliceEQ* treeSliceEQ = dynamic_cast<TreeSliceEQ*>(this);
    if (treeSliceEQ) {
        treeSliceEQ->eval(expr);
        return;
    }

    TreeSliceRates* treeSliceRates = dynamic_cast<TreeSliceRates*>(this);
    if (treeSliceRates) {
        treeSliceRates->eval(expr);
        return;
    }

	TreeSliceRatesCompact* treeSliceRatesCompact = dynamic_cast<TreeSliceRatesCompact*>(this);
    if (treeSliceRatesCompact) {
        treeSliceRatesCompact->eval(expr);
        return;
    }

    TreeSliceNDimCounting* treeSliceNDimCounting = dynamic_cast<TreeSliceNDimCounting*>(this);
    if (treeSliceNDimCounting) {
        treeSliceNDimCounting->eval(expr);
        return;
    }

    TreeSliceBasic* treeSliceBasic = dynamic_cast<TreeSliceBasic*>(this);
    if (treeSliceBasic) {
        treeSliceBasic->eval(expr);
        return;
    }

    TreeSliceGeneral* treeSliceGeneral = dynamic_cast<TreeSliceGeneral*>(this);
    if (treeSliceGeneral) {
        treeSliceGeneral->eval(expr);
        return;
    }

    throw ModelException("template<class A> TreeSlice& operator=", "Slice type not handled");
}

DRLIB_END_NAMESPACE

#endif
