#include "edginc/ran2.h"
#include "edginc/gasdev2.h"
using namespace core;
int main(int argc, char**argv) {
    long size = (argc > 1) ? atol(argv[1]) : 1000;
    vector<double> randoms(size);
    SC_gasdev2SP rng = SC_gasdev2::create(SC_ran2Gen::create());
    std::generate(randoms.begin(), randoms.end(), *rng);
    cout << randoms[0] << " " << randoms.back() << endl;
}
