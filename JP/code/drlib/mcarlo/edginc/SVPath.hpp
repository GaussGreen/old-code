//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVPath.hpp
//
//   Description : A Legacy non-type-specific access mode to the Spot Variable 
//
//
//
//----------------------------------------------------------------------------


#ifndef SVPath_HPP
#define SVPath_HPP

#include "edginc/StateVariable.hpp"

DRLIB_BEGIN_NAMESPACE

/** A Legacy non-type-specific access mode to the Spot Variable */
    /** Class for returning equity spot values (for one asset) in a MC
        simulation. This is used by products in the payoff() method. An
        SVPath is obtained by calling the method path() on the
        MCPath::IStateVar object below. Payoffs should obtain a
        MCPath::IStateVar object in the pathGenUpdated() method by calling
        getSpotSV() on this (ie MCPath) object
    */

class MCARLO_DLL SVPath : public IAdvanceableStateVariable {
public:
    virtual ~SVPath();

    /** Constructs a SVPath where the spots are not in a
        contiguous block */
    SVPath(const double*      ppath,
            const vector<int>& offsets,
            int                beginIdx,
            int                endIdx);

    /** Constructs 'trivial' SVPath object where the spots are in a
        contiguous block */
    SVPath(const double* ppath,
            int beginIdx,
            int endIdx);

    /** Default constructor - to be used only together with any initialize() method below */
    SVPath();

    void initialize(const double*      ppath,
                    const vector<int>& offsets,
                    int                beginIdx,
                    int                endIdx);

    void initialize(const double*      ppath,
                    int beginIdx,
                    int endIdx);

    void initialize(int beginIdx,
                    int endIdx);

    /** Return the spot level for the specified index. Inline for
        performance */
    double operator[](int pathIndex) const;

    // for stateless payoffs
    double getValue() const;

    // implementing this interface function
    void reset();

    // implementing this interface function
    void advance();

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    int begin() const;

    //// see  begin()
    int end() const;

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables -
        need to see how hard it is to implement */
    virtual bool doingPast() const{ return false; }

    /** This virtual function allows for dynamic data access 'on-the-fly'
        while preserving the existing interface and its efficient mode */
    virtual double element(int idx) const;

private:
    SVPath(const SVPath& rhs); // do not stucture copy this class
    SVPath& operator=(const SVPath& rhs); // do not stucture copy this class

    ////// fields //////
    const double*  ppath;
    vector<int>    offsets;
    int            beginIdx;
    int            endIdx;
    int            currentIdx;
    bool           isTrivial;
};
typedef smartPtr<SVPath> SVPathSP;
typedef vector<SVPathSP> SVPathArray;

#if !defined(DEBUG) || defined(QLIB_SVPATH_CPP)
//////////////////////////////////////////////////////////////////////////
// inline these functions for optimised build

/** Return the spot level for the specified index. Inline for
    performance */
OPT_INLINE double SVPath::operator[](int pathIndex) const {
    return ppath ? ppath[isTrivial ? pathIndex: offsets[pathIndex]] : element(pathIndex);
}

// for stateless payoffs
OPT_INLINE double SVPath::getValue() const {
    return (*this)[currentIdx];
}

OPT_INLINE void SVPath::reset() {
    currentIdx = beginIdx;
}

OPT_INLINE void SVPath::advance() {
    ++currentIdx;
}

/** Indicates which values in array returned by Path() is being
    simulated. Clients should use a loop of the form:
    "for (iStep = begin(); iStep < end(); iStep++)" */
OPT_INLINE int SVPath::begin() const { return beginIdx; }

//// see  begin()
OPT_INLINE int SVPath::end() const { return endIdx; }

//////////////////////////////////////////////////////////////////////////
#endif

DRLIB_END_NAMESPACE

#endif
