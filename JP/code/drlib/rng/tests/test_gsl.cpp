#include "edginc/GSLRNG.h"
using namespace core;
int main(int argc, char**argv)
{
    long size = (argc > 1 ) ? atol(argv[1]) : 1000;
    vector<double> randoms(size);
    GSLSuperCubeGenSP gen = GSLSuperCubeGen::create(); // generator/seed can be set via en. variable
    GSLZigguratNormalSuperCubeRNGSP rng = GSLZigguratNormalSuperCubeRNG::create(gen);
    std::generate(randoms.begin(), randoms.end(), *rng);
    cout << randoms[0] << " " << randoms.back() << endl;
}

