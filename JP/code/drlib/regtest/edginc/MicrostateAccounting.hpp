//Microstate acccounting code for accurate and consistent timing
//of execution. Only available on Solaris.
#ifdef sun

#include "edginc/config.hpp"
#include "edginc/ModelException.hpp"
#include <fstream>
#include <exception>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/procfs.h>
#include <unistd.h>
#include <stropts.h>
#include <string>

class MicrostateAccounting  {
private:
    static const string timeExtention;
    static const string difExtention;
    static const string outExtention;
    static const int    maxLoopsNew;
    static const int    maxLoopsExisting;
    static const double thresholdNew;
    static const double thresholdExisting;
    static const char   commentToken;
    static const int    failedTestCode;
	bool previousTimingFound; 
    hrtime_t startCycles;
    hrtime_t oldTiming;
    ifstream inFile;
    ofstream outFile;
    vector<hrtime_t> timings;
    static int toggleMicrostateAccounting(int onoff);
    static double relDiff(double a, double b);

public:
    MicrostateAccounting(string outFilename, int permNum);
    void recordStartTime();
    void recordStopTime(int& loops);
    void reportTime();
    void reportFailedTest();
};

#endif
