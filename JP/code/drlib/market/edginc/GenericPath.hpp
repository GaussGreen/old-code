//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : GenericPath.hpp
//
// Description : Representing the path of a MonteCarlo run. It is the templatised form of MCPath::Path
//
// Date        : Aug 2006
//
//----------------------------------------------------------------------------



#ifndef GENERICPATH_HPP
#define GENERICPATH_HPP

#include "edginc/StateVariable.hpp"

DRLIB_BEGIN_NAMESPACE

template <class tuple>
class MARKET_DLL GenericPath : public IStateVariable
{
public:
    /** Virtual Destructor */
    virtual ~GenericPath() {};

    /** Constructs a MCPath::Path where the spots are not in a contiguous block */
    GenericPath(
		const tuple*		ppath, 
        const vector<int>&	offsets,
        int					beginIdx, 
        int					endIdx)
		:ppath(ppath), 
		offsets(offsets), 
		beginIdx(beginIdx), 
		endIdx(endIdx), 
		currentIdx(beginIdx), 
		isTrivial(false) {}

    /** Constructs 'trivial' MCPath::GenericPath object where the spots are in a 
        contiguous block */
    GenericPath(
		const tuple* ppath, 
        int beginIdx, 
        int endIdx)
		:ppath(ppath), 
		offsets(0), 
		beginIdx(beginIdx), 
		endIdx(endIdx), 
		currentIdx(beginIdx), 
		isTrivial(true) {}    

    /** Default constructor - to be used only together with any initialize() method below */
    GenericPath():
        ppath(0), 
		offsets(0), 
		beginIdx(0), 
		endIdx(0), 
		currentIdx(0), 
		isTrivial(true) {}           

    /** returns the 'path'  Use
        begin() and end() on MCPath::Path to see which values are being
        simulated. It is valid to access values before begin() but they
        will be zero. (for inherited classes) */
    virtual const GenericPath& path() const 
	{ 
		return *this; 
	}

    /** Initializes the member variables */
    void initialize
		(	const tuple*		ppath, 
			const vector<int>&	offsets,
			int					beginIdx, 
			int					endIdx)
    {
        this->ppath = ppath; 
        this->offsets = offsets;
        this->beginIdx = beginIdx;
        this->endIdx = endIdx;
        this->currentIdx = beginIdx;
        this->isTrivial = false;
    }

    /** Initializes the member variables */
    void initialize(
		const tuple*	ppath, 
        int				beginIdx, 
        int				endIdx)
    {
        this->ppath = ppath; 
        this->offsets.clear();
        this->beginIdx = beginIdx;
        this->endIdx = endIdx;
        this->currentIdx = beginIdx;
        this->isTrivial = true;
    }

    /** Initializes the member variables */
    void initialize(
		int beginIdx, 
        int endIdx)
    {
        this->ppath = 0;
        this->offsets.clear();
        this->beginIdx = beginIdx;
        this->endIdx = endIdx;
        this->currentIdx = beginIdx;
        this->isTrivial = true;
    }
    
    /** Return the spot level for the specified index. Inline for
        performance */
    OPT_INLINE const tuple& operator[](int pathIndex) const
    {
        return this->ppath ? 
							this->ppath[this->isTrivial ? 
								pathIndex
								:this->offsets[pathIndex]] 
							: this->element(pathIndex);                
    }

    /** for stateless payoffs */
    OPT_INLINE const tuple& getValue() const
    {
        return (*this)[this->currentIdx];
    }

    /** implementing this interface function */
    OPT_INLINE void reset()
    {
        this->currentIdx = this->beginIdx;
    }

    /** implementing this interface function */
    OPT_INLINE void advance()
    {
        ++this->currentIdx;
    }

    /** Indicates which values in array returned by Path() is being
        simulated. Clients should use a loop of the form:
        "for (iStep = begin(); iStep < end(); iStep++)" */
    int begin() const
    {
        return this->beginIdx;
    }
    
    /** see comments for begin() */
    int end() const
    {
        return this->endIdx;
    }

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const
	{
		return false; 
	}

    /** This virtual function allows for dynamic data access 'on-the-fly'
        while preserving the existing interface and its efficient mode */
    virtual const tuple& element(int idx) const
    {
        if (this->ppath == 0)
            throw ModelException("GenericPath<tuple>::element", "No elements present");
        if (idx < this->beginIdx) 
            throw ModelException("GenericPath<tuple>::element", "Element idx is less than beginIdx");
        if (idx > this->endIdx)
            throw ModelException("GenericPath<tuple>::element", "Element idx is larger than endIdx");
        return this->ppath[idx];
    }

private:
    GenericPath(const GenericPath& rhs); // do not stucture copy this class
    GenericPath& operator=(const GenericPath& rhs); // do not stucture copy this class
        
// ## Member variables
    const tuple*   ppath;
    vector<int>    offsets;
    int            beginIdx;
    int            endIdx;
    int            currentIdx;
    bool           isTrivial;
};

DRLIB_END_NAMESPACE
#endif
