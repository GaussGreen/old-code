//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ArrayMD.hpp
//
//   Description : Multidimensional Array
//
//   Date        : Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MARRAY_HPP
#define EDR_MARRAY_HPP

#if defined(_MSC_VER)
// disable warning
#pragma warning(disable : 4503)
#pragma warning(disable : 4284)
#endif

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Format.hpp"
#include <iterator>

DRLIB_BEGIN_NAMESPACE

/** Multidimensional array. Implemented as a contiguous array of elementary
    'CompArray' arrays. The elements of the array can be accessed either via:
        - an 'Index' object; or, via
        - an 'Iterator' object.
    An iterator can be 'turned' into an index, but not the other way around */
template<
class Comp,                                            // component type
    class CompArray = array<Comp>,                         // component array type
    class CompArraySP = smartPtr<CompArray>,               // smart ptr to comp array type
#if (defined(__GNUC__) && ( __GNUC__ >= 3)) || (defined(_MSC_VER) && (_MSC_VER >= 1300)) // VC 7 
    class CompArray_Iterator = typename vector<Comp>::iterator>     // comp array iterator type
#else
    class CompArray_Iterator = vector<Comp>::iterator>     // comp array iterator type
#endif

class TOOLKIT_DLL ArrayMD: public CObject{
public:
    static CClassConstSP const TYPE;    

    /** Size of the array. For instance, for a M * N, 2-dimensional MArray, the 'Size'
        of the array would be (M, N)*/
    typedef CIntArray Size;

private:
    typedef vector<CompArraySP>         CompArrayArray;     // array of elementary arrays
    typedef refCountPtr<CompArrayArray> CompArrayArraySP;
    
    struct ContIndex;

public:
    class Index;
    friend class Index;

    /** The 'Index' class supports elementary operators, such as ++, ==, !=, <.
        Given an 'ArrayMD' m and an 'Index' idx, the element pointed to by idx
        can be accessed as follows: m[idx] */
    class Iterator;
    class TOOLKIT_DLL Index{
    public:        
        friend class Iterator;
        friend struct ContIndex;

        Index& operator++(){
            // if we've reached the last index, nothing to do
            if (pos > 0){
                return *this;
            }
            // if before first index, set to first
            if (pos < 0){
                return (*this = IndexFirst(size));
            }
            // otherwise, increment indices starting from the right
            int i = dim - 1;
            for (; i >= 0; --i){
                ++idx[i];
                if (idx[i] < size[i]){
                    break;
                }
                idx[i] = 0;
            }
            if (i >= 0){
                return *this;
            }
            // if here, we've reached the end of the sequence
            return (*this = IndexLast(size));
        }

        bool operator==(const Index& rhs){
            if (dim != rhs.dim){
                throw ModelException(ArrayMD::TYPE->getName() + "::Index::operator==",
                                     "Cannot compare two indices that have different dimensions");
            }
            if (pos != rhs.pos){
                return false;
            }
            if (!pos){
                return (idx == rhs.idx);
            }
            return true;
        }

        bool operator!=(const Index& rhs){
            if (dim != rhs.dim){
                throw ModelException(ArrayMD::TYPE->getName() + "::Index::operator!=",
                                     "Cannot compare two indices that have different dimensions");
            }
            if (pos != rhs.pos){
                return true;
            }
            if (!pos){
                return (idx != rhs.idx);
            }
            return false;
        }

        bool operator<(const Index& rhs){
            if (dim != rhs.dim){
                throw ModelException(ArrayMD::TYPE->getName() + "::Index::operator<",
                                     "Cannot compare two indices that have different dimensions");
            }
            if (pos || rhs.pos){
                return (pos < rhs.pos);
            }
            int i = 0;
            for (; i < dim; ++i){
                if (idx[i] >= rhs.idx[i]){
                    return false;
                }
            }
            return true;
        }

        int operator[](int i) const{
            return idx[i];
        }

        int& operator[](int i){
            return idx[i];
        }

    protected:
        // First index ctor
        Index(const Size& size):
        size(size),
        dim(size.size()),
        idx(dim, 0),
        pos(0){}

        // Last index ctor
        Index(const Size& size,
              bool is_last):
        size(size),
        dim(size.size()),
        idx(calcLastIdx(size)),
        pos(1){}

    private:
        Index& operator=(const Index& rhs){
            if (this == &rhs){
                return *this;
            }
            idx = rhs.idx;
            pos = rhs.pos;
            return *this;
        }

        // Sets the index to the last valid index, + 1
        static vector<int> calcLastIdx(const Size& size){
            int dim = size.size();
            vector<int> idx(dim);
            int i = 0;
            for (; i < dim - 1; ++i){
                idx[i] = size[i] - 1;
            }
            idx[i] = size[i];
            return idx;
        }

#if 0
        // Turns a contiguous index into an 'Index'
        static vector<int> fromCIndex(const Size& size, const ContIndex& cont_idx){
            int dim = size.size();
            vector<int> idx(dim);
            idx[dim-1] = cont_idx.idx_back;
            int i = dim - 2;
            CompArrayArray::size_type j = cont_idx.c_idx;
            for(; i >= 0; --i){
                idx[i] = j % size[i];
                j /= size[i];
            }
            return idx;
        }
#endif

        const Size& size;
        const int   dim;
        vector<int> idx;
        int         pos;   // position: -1 if before sequence, 0 if in sequence, +1 if after sequence
    };

private:
    class TOOLKIT_DLL IndexFirst: public Index{
    public:
        IndexFirst(const Size& size):
        Index(size){}
    };
    typedef refCountPtr<IndexFirst> IndexFirstSP;

    class TOOLKIT_DLL IndexLast: public Index{
    public:
        IndexLast(const Size& size):
        Index(size, true){}
    };
    typedef refCountPtr<IndexLast> IndexLastSP;

private:
    /* Represents a contiguous index. A contiguous index is made of 2 elements:
       i) 'c_idx', which points to an elementary array
       ii) 'idx_back', which points to an element of the elementary array pointed
       to by c_idx */
    struct ContIndex{
        typename CompArrayArray::size_type c_idx;
        int         idx_back;
        int         pos;

        static typename CompArrayArray::size_type calc_c_size(
            const Size& size){
            int dim = size.size();
            int i = dim - 2;
            typename CompArrayArray::size_type p = 1;
            for(; i >= 0; --i){
                p *= size[i];
            }
            return p;
        }
        
        // ctor from Index to ContIndex
        ContIndex(const Index& index):
        c_idx(to_c_index(index)),
        idx_back(index.idx.back()),
        pos(index.pos){}

    private:
        // Convert the dim-1 first elements of index into a contiguous index
        static typename CompArrayArray::size_type to_c_index(
            const Index& index){
            int i = index.dim - 2;
            typename CompArrayArray::size_type c_idx = 0;
            typename CompArrayArray::size_type p = 1;
            for(; i >= 0; --i){
                c_idx += index[i] * p;
                p *= index.size[i];
            }
            return c_idx;
        }
    };

public:
    friend class Iterator;

    /** The 'Iterator' class is a 'forward' iterator. It supports elementary 
        operators, such as ++, ==, !=, *, ->. */
    class TOOLKIT_DLL Iterator{
    public:
        friend class Index;

	    typedef forward_iterator_tag iterator_category;
	    typedef Comp                 value_type;
        typedef Comp&                reference_type;
        typedef Comp*                pointer_type;

        reference_type operator*(){
            return *curr_comp_it;
        }

        pointer_type operator->(){
            return &*curr_comp_it;
        }

        Iterator& operator++(){
            // if finished, nothing to do
            if (curr_sub_seq_it == seq_end){
                return *this;
            }
            // increment index and iterator to current comp
            ++curr_index;
            ++curr_comp_it;
            // if we've moved sub sequence, we must increase iterator to current
            // subsequence
            if (curr_comp_it == sub_seq_end){
                ++curr_sub_seq_it;
                // if we haven't reached the end of the sequence, we must reset 
                // the iterator to the current comp
                if (curr_sub_seq_it != seq_end){
                    curr_comp_it = (*curr_sub_seq_it)->begin();
                }
            }
            return *this;
        }

        bool operator==(const Iterator& rhs){
            if (theArray != rhs.theArray){
                throw ModelException(ArrayMD::TYPE->getName() + "::Iterator::operator==",
                                     "Cannot compare two iterators that do not refer to the same ArrayMD");
            }
            return (curr_index == rhs.curr_index);
        }

        bool operator!=(const Iterator& rhs){
            if (theArray != rhs.theArray){
                throw ModelException(ArrayMD::TYPE->getName() + "::Iterator::operator!=",
                                     "Cannot compare two iterators that do not refer to the same ArrayMD");
            }
            return (curr_index != rhs.curr_index);
        }

        const Index& index() const{
            return curr_index;
        }

        // copy ctor
        Iterator(const Iterator& rhs):
        theArray(rhs.theArray),
        curr_index(rhs.curr_index),
        curr_sub_seq_it(rhs.curr_sub_seq_it),
        seq_end(rhs.seq_end),
        curr_comp_it(rhs.curr_comp_it),
        sub_seq_end(rhs.sub_seq_end){}

        Iterator& operator=(const Iterator& rhs){
            if (this == &rhs){
                return *this;
            }
            theArray  = rhs.theArray;
            curr_index = rhs.curr_index;
            curr_sub_seq_it = rhs.curr_sub_seq_it;
            seq_end = rhs.seq_end;
            curr_comp_it = rhs.curr_comp_it;
            sub_seq_end = rhs.sub_seq_end;
            return *this;
        }

    protected:        
        // 'begin' iterator ctor
        Iterator(ArrayMD<Comp, CompArray, CompArraySP>* theArray):
        theArray(theArray),
        curr_index(theArray->first()),
        curr_sub_seq_it(theArray->v->begin()),
        seq_end(theArray->v->end()),
        curr_comp_it((*curr_sub_seq_it)->begin()),
        sub_seq_end((*curr_sub_seq_it)->end()){}

        // 'end' iterator ctor
        Iterator(ArrayMD<Comp, CompArray, CompArraySP>* theArray,
                 bool is_end):
        theArray(theArray),
        curr_index(theArray->last()),
        curr_sub_seq_it(theArray->v->end()),
        seq_end(theArray->v->end()),
        curr_comp_it(),    // not used
        sub_seq_end((*curr_sub_seq_it)->end()){}

    private:        
        ArrayMD<Comp, CompArray, CompArraySP, CompArray_Iterator>*  theArray;    // the array refered to
        Index                    curr_index;        // current index 
        typename CompArrayArray::iterator curr_sub_seq_it;   // points to current sub sequence
        typename CompArrayArray::iterator seq_end;           // end of sequence
        CompArray_Iterator       curr_comp_it;      // points to current component
        CompArray_Iterator       sub_seq_end;       // end of current sub sequence
    };

private:
    class IteratorBegin: public Iterator{
    public:
        IteratorBegin(ArrayMD<Comp, CompArray, CompArraySP, CompArray_Iterator>* theArray):
        Iterator(theArray){}
    };
    typedef refCountPtr<IteratorBegin> IteratorBeginSP;

    class IteratorEnd: public Iterator{
    public:
        IteratorEnd(ArrayMD<Comp, CompArray, CompArraySP, CompArray_Iterator>* theArray):
        Iterator(theArray, true){}
    };
    typedef refCountPtr<IteratorEnd> IteratorEndSP;

public:
    ArrayMD(const Size& size):
    CObject(TYPE),
    thesize(size),
    dim(size.size()){
        // Get size of contiguous array and allocate mem of that size
        typename CompArrayArray::size_type c_size = 
            ContIndex::calc_c_size(size);
        v = CompArrayArraySP(new CompArrayArray(c_size));
        // Size of back index
        int idx_back_size = size[dim-1];
        // Allocate mem for entire cont array
        typename CompArrayArray::size_type c_idx = 0;
        for (; c_idx < c_size; ++c_idx){
            (*v)[c_idx] = CompArraySP(new CompArray(idx_back_size));
        }
        // Create begin/end iterators + first/last indices
        beginIt = IteratorBeginSP(new IteratorBegin(this));
        endIt = IteratorEndSP(new IteratorEnd(this));
        firstIdx = IndexFirstSP(new IndexFirst(thesize));
        lastIdx = IndexLastSP(new IndexLast(thesize));
    }

    // copy ctor
    ArrayMD(const ArrayMD& rhs):
    CObject(TYPE),
    thesize(rhs.thesize),
    dim(rhs.dim),
    v(rhs.v),               // not a deep copy
    beginIt(rhs.beginIt),   // not a deep copy
    endIt(rhs.endIt),       // not a deep copy
    firstIdx(rhs.firstIdx), // not a deep copy
    lastIdx(rhs.lastIdx){}  // not a deep copy

    ArrayMD& operator=(const ArrayMD& rhs){
        thesize = rhs.thesize;
        dim = rhs.dim;
        v = rhs.v;                // not a deep copy
        beginIt = rhs.beginIt;    // not a deep copy
        endIt = rhs.endIt;        // not a deep copy
        firstIdx = rhs.firstIdx;  // not a deep copy
        lastIdx = rhs.lastIdx;    // not a deep copy
        return *this;
    }

    int dimension() const{
        return dim;
    }

    const Size& size() const{
        return thesize;
    }

    const Comp& operator[](const Index& index) const{
        ContIndex cont_idx(index);
        return (*(*v)[cont_idx.c_idx])[cont_idx.idx_back];
    }

    Comp& operator[](const Index& index){
        ContIndex cont_idx(index);
        return (*(*v)[cont_idx.c_idx])[cont_idx.idx_back];
    }

    /** Returns 'begin' iterator */
    const Iterator& begin(){
        return *beginIt;
    }

    /** Returns 'end' iterator */
    const Iterator& end(){
        return *endIt;
    }

    /** Returns 'first' index */
    const Index& first() const{
        return *firstIdx;
    }

    /** Returns 'last' index (which corresponds to the 'end' iterator).
        For instance, for a 1-dim array [0..n-1] last would return n */
    const Index& last() const{
        return *lastIdx;
    }

    virtual IObject* clone() const{
        static const string method(ArrayMD::TYPE->getName() + "::clone");
        try{
            IObjectSP objcopy(CObject::clone());
            ArrayMD<Comp, CompArray, CompArraySP, CompArray_Iterator>& thecopy
                = dynamic_cast<ArrayMD<Comp, CompArray, CompArraySP, CompArray_Iterator>&>(*objcopy);
            thecopy.v = v;                // not a deep copy
            thecopy.beginIt = beginIt;    // not a deep copy
            thecopy.endIt = endIt;        // not a deep copy
            thecopy.firstIdx = firstIdx;  // not a deep copy
            thecopy.lastIdx = lastIdx;    // not a deep copy
            return &thecopy;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

private:

    // for reflection
    ArrayMD():
    CObject(TYPE),
    dim(0){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ArrayMD, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(thesize, "");
        FIELD_MAKE_TRANSIENT(thesize);
        FIELD(dim, "");
        FIELD_MAKE_TRANSIENT(dim);
        // not public for now
    }

    static IObject* defaultCtor(){
        return new ArrayMD();
    }
    
    Size             thesize;
    int              dim;

    CompArrayArraySP v;         // $unregistered
    IteratorBeginSP  beginIt;   // $unregistered
    IteratorEndSP    endIt;     // $unregistered
    IndexFirstSP     firstIdx;  // $unregistered
    IndexLastSP      lastIdx;   // $unregistered
};

typedef ArrayMD<double> DoubleArrayMD;
typedef smartPtr<DoubleArrayMD> DoubleArrayMDSP;
#ifndef QLIB_ARRAYMD_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL ArrayMD<double>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DoubleArrayMD>);
#endif
DRLIB_END_NAMESPACE
#endif
