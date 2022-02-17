
//
// C++ Interface: BoostRNGGen
//
// Description: Wrappers around Boost Random Number Generators (found in <boost/random.hpp>)
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//


#ifndef BOOSTRNGGEN_HPP
#define BOOSTRNGGEN_HPP

#include "edginc/coreConfig.hpp"
#include "edginc/IRNGGenerator.h"

#include <boost/random.hpp>

CORE_BEGIN_NAMESPACE


/** Wrapper template to bridge Boost Generators with IUniformRNGGen
    The implementation is largely based on the Boost example file random_demo.cpp
*/

template <class rand_gen>
class BoostUniformGen : public IUniformRNGGen {

public:

    typedef rand_gen base_generator_type;
    typedef boost::uniform_real<> distribution_type; //FIXME: not sure why, but uniform_01 didn't work on the first try; reuse their sample code here.
    typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
    
    DECLARESP(BoostUniformGen);

    static BoostUniformGenSP create(unsigned seed = 0) {
        return BoostUniformGenSP(new BoostUniformGen(seed));
    }
    virtual double fetch() ///< Fetch next random number
    {
        return uni_gen();
    }

    virtual IRNGGeneratorSP  clone() const  ///< Clone the existing state
    {
        return IRNGGeneratorSP(new BoostUniformGen(uni_gen));
    }

protected:
    BoostUniformGen(unsigned seed = 0) :
        generator(seed),
        uni_dist(0.0, 1.0),
        uni_gen(generator, uni_dist)
    {
    }

private:

    base_generator_type generator;
    distribution_type   uni_dist;
    gen_type            uni_gen;

    BoostUniformGen(gen_type _gen) :
        uni_gen(_gen)
    {}
    
};

/**  For the full list of generators and their properties (perf./memory) see documentation on the Boost web site: 
     http://www.boost.org/libs/random/random-generators.html

     The boost::mt19937 should be a good start for many applications. 
 */

typedef BoostUniformGen< boost::minstd_rand > BoostMindStdUniformGen; // initialize with seed != 0
typedef BoostUniformGen< boost::mt19937 > BoostMT19937UniformGen; // mersenne_twister 19937
typedef BoostUniformGen< boost::mt11213b > BoostMT11213bUniformGen; // mersenne_twister 11213b
typedef BoostUniformGen< boost::lagged_fibonacci44497 > BoostLaggedFibonachi44497UniformGen; // initialize with seed != 0




CORE_END_NAMESPACE
#endif
