//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : Shuffle.cpp
//
//   Description : Generic Knuth shuffle
//
//   Author      : Jon Dee
//
//   Date        : 17 November 2007
//
//----------------------------------------------------------------------------

#ifndef EDR_SHUFFLE_HPP
#define EDR_SHUFFLE_HPP

DRLIB_BEGIN_NAMESPACE

// Knuth shuffle: O(N), works with vectors and QLib arrays.
template<typename T, typename R>
static void shuffle(T& theArray, const R& rnd) {
    int rand;
    int nTop = (int)theArray.size()-1;
    for(int i=nTop;i>=0;i--)
    {
        rand = rnd.getRandomInteger(i);
        std::swap(theArray[i], theArray[rand]);
    }
}

// An example of a type that can be used with shuffle, 
// but beware it'll not be a perfect shuffle: 
// rand() is not great.
class basicRandGen {
public:
    basicRandGen(int seed) { 
        srand(seed);
    }
    int getRandomInteger(int maxInt) const { 
        int rnd = (rand() % (maxInt + 1));
        return rnd;
    }
};


DRLIB_END_NAMESPACE

#endif
