//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Lattice.hpp
//
//   Description : Implementation of a Lattice: "Maturity x Strikes".
//                 Different size across Maturity
//                 Typically used by Local Vols, but Tree & MC could as well
//
//   Author      : JNJ after Regis G.
//
//   Date        : 18 Sep 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_LATTICE_HPP
#define EDR_LATTICE_HPP

#if defined(_MSC_VER)
// disable warning
#pragma warning(disable : 4355)
#endif

#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

// If you use this template for something other than doubles and you get
// undefined symbols when you link then change the code below which stops
// the inclusion of Lattice.cpp

// Template classes representing a Slice and a Lattice.
// A slice is essentially an array.
// A lattice is a collection (array) of slices.
// The added feature is that a slice is also a (1-slice) lattice.
// This means that routines that take a lattice as an argument can also take a slice.

template<class Comp,
         class CompArray,
         class CompArraySP>
class CSlice;


template<class Comp,
         class CompArray = array<Comp>,
         class CompArraySP = smartPtr< array<Comp> > >
class TOOLKIT_DLL CLattice {
private:
    vector< CSlice<Comp, CompArray, CompArraySP>* > m;
public:
    CLattice(const vector<int>& sizes,
             Comp               val = Comp());

    // Copy constructor
    CLattice(const CLattice<Comp, CompArray, CompArraySP>& rhs);

    /** Return the size of each of the slices making up the lattice */
    vector<int> sizes() const;

    /** Returns the number of slices in the lattice */
    int size() const;

    /** Easy access to the i'th slice in the lattice */
    const CSlice<Comp, CompArray, CompArraySP>& operator[](int i) const;

    /** Easy access to the i'th slice in the lattice */
    CSlice<Comp, CompArray, CompArraySP>& operator[](int i);

    /** scales all values in the lattice by supplied factor */
    void scale(double factor);

    virtual ~CLattice();
protected:
    CLattice(int size, CSlice<Comp, CompArray, CompArraySP>* p,
             bool isOwner=true);

    void setToNull();

private:
    void clean();

    bool    isOwner;
};

template<class Comp,
         class CompArray = array<Comp>,
         class CompArraySP = smartPtr< array<Comp> > > 
class TOOLKIT_DLL CSlice: public CLattice<Comp, CompArray, CompArraySP> {
private:
    CompArraySP v_sptr;

    // added to allow input array
    Comp*  v_ptr;
    int    v_size;

    CSlice(const CompArraySP& inarray);
public:
    static CSlice<Comp, CompArray, CompArraySP> fromArray(
        const CompArraySP& inarray);

    const CompArray& toArray() const;

    CSlice<Comp, CompArray, CompArraySP>& operator=(
        const CSlice<Comp, CompArray, CompArraySP>& rhs);

#if defined(__GNUC__) && ( __GNUC__ >= 3)
    // added to allow input array - where iterators aren't typedefd to pointers
    CSlice(typename CompArray::iterator in, int size);
#endif

    // added to allow input array
    CSlice(Comp* in, int size);

    CSlice(int  size,
           Comp val = Comp());

    CSlice(const CSlice<Comp, CompArray, CompArraySP>& rhs);

    int size() const;

    const Comp& operator[](int i) const;

    Comp& operator[](int i);

    /** scales all values in the slice by supplied factor */
    void scale(double factor);

    virtual ~CSlice();
};

typedef CLattice<double> CLatticeDouble;
typedef CSlice  <double> CSliceDouble;

DECLARE_REF_COUNT(CLatticeDouble);

#if !defined(DEBUG) || defined(QLIB_ARRAYMD_CPP)
// currently only used for doubles so avoid cost of including it
#include "edginc/Lattice.inl"
#endif

#ifndef QLIB_ARRAYMD_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL CLattice<double>);
EXTERN_TEMPLATE(class TOOLKIT_DLL CSlice<double>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL CLattice<double>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL CSlice<double>);
#endif

DRLIB_END_NAMESPACE
#endif
