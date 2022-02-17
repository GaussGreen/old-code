#include "edginc/BoostRNGGen.h"
#include "edginc/gasdev2.h"

using namespace core;

int main(int argc, char**argv) {
    long size = (argc > 1) ? atol(argv[1]) : 1000;
    vector<double> randoms(size);
    //    BoostUniformGen< boost::lagged_fibonacci44497 > rng(1);

    IUniformRNGGenSP gen = BoostMT19937UniformGen::create(1);
    INormalRNGSP     rng = GasDev2::create(gen);

    for(size_t i =0; i < randoms.size(); ++i)
        cout << (randoms[i]=rng->fetch()) << " ";
    cout << endl;
}
