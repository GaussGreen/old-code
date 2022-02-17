//Microstate acccounting code for accurate and consistent timing
//of execution. Only available on Solaris.
#ifdef sun
#undef __STDC__  
#define __STDC__ 0
#endif
#include "edginc/config.hpp"
#ifdef sun
#undef __STDC__  
#define __STDC__ 0

#include "edginc/MicrostateAccounting.hpp"
//#include "edginc/TRACE.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
         
// The microstate accounting isn't ANSI C++, so disable it for this file.
// Extention of dif directories. Can't currently calculate this.
const string MicrostateAccounting::timeExtention = ".time";
// Extention of dif directories. Can't currently calculate this.
const string MicrostateAccounting::difExtention = "dif-solaris.opt";
// The equivalent for dif directories.
const string MicrostateAccounting::outExtention = "out";
// The maximum number of loops to try to converge result for new timing.
const int MicrostateAccounting::maxLoopsNew = 30;
// The equivalent if there already is a timing.
const int MicrostateAccounting::maxLoopsExisting = 5;
// The greatest relative difference allowed between the two shortest timings.
const double MicrostateAccounting::thresholdNew = 0.02;
// The equivalent if previous timings exist.
const double MicrostateAccounting::thresholdExisting = 0.10;
// Comment prefix.
const char MicrostateAccounting::commentToken = '#';
// Comment prefix.
const int MicrostateAccounting::failedTestCode = 0;

int MicrostateAccounting::toggleMicrostateAccounting(int onoff) {
    char buffer[80];
    int fd;
    int state = PR_MSACCT;
    sprintf (buffer, "/proc/%d", getpid ());
    if ((fd = open(buffer, O_RDWR)) == -1) {
        perror("open:");
        cerr << "MicrostateAccounting: "
                "Can't open process file, "
                "can't run with -time option." << endl;
        exit(1);
    }
    if (ioctl(fd, (onoff ? PIOCSET : PIOCRESET), &state) == -1) {
        perror("PIOCSET/PIOCRESET");
        close(fd);
        cerr << "MicrostateAccounting: "
                "Microstate accounting not available,"
                "can't run with -time option." << endl;
        exit(1);
    }
    close(fd);
    return(0);
}

double MicrostateAccounting::relDiff(const double a, const double b) {
    return fabs((b-a)/(double)a);
}

 
hrtime_t convertToInteger(const std::string& s)
{
    std::istringstream i(s);
    hrtime_t x;
    if (!(i >> x)) {
        throw ModelException("convertToInteger(\"" + s + "\")");
    }
    return x;
}

MicrostateAccounting::MicrostateAccounting(string outFilename, const int permNum)
    : oldTiming(0), previousTimingFound(false)  {
    toggleMicrostateAccounting(1);
    // Create filenames
    outFilename += timeExtention;
    string inFilename = outFilename;
    //TRACE("Output file name: " + outFilename); 
    inFilename.replace(outFilename.find(difExtention),
                       difExtention.length(), outExtention);
    errno = 0;
    inFile.open(inFilename.c_str());
    if (errno != 0) {
       char const* message = strerror( errno ) ; 
       string errorString= "Cannot open file " +  inFilename + ", error was: " + message;
       //TRACE(errorString);
    }
    int i = 0;
    if (inFile.is_open()) {
        for (i = 0; i < permNum; ++i) {         // Find right permutation
            string buffer; 
            while (inFile.peek() == commentToken) { // Remove comments
                getline(inFile, buffer);
            }
            getline(inFile,buffer);
            oldTiming = convertToInteger(buffer);       
        }
	previousTimingFound = true;
        //TRACE("Previous timing found (" + oldTiming + ") for permutation " + i + " in file " + inFilename) ;
    }
    inFile.close();
    // Open in append mode if this is not the first permutation
    outFile.open(outFilename.c_str(),
                 (permNum > 1) ? (ios::app) : (ios::out));
    if (permNum == 1) { // Write comment header
		outFile << commentToken 
                << " This file contains microstate Solaris timings for the test"
                   " (" << failedTestCode << " meaning a failed test)"
                << endl;
    }
}

void MicrostateAccounting::recordStartTime() {
    startCycles = gethrvtime(); // Get microstate timing
}

void MicrostateAccounting::recordStopTime(int &loops) {
    hrtime_t timeTaken = gethrvtime () - startCycles;
    timings.push_back(timeTaken); // Add timing to our list
    //TRACE("Timed: " + timeTaken);
    sort(timings.begin(), timings.end()); 
    if (previousTimingFound) { // No previous timing found, create new one.
        if ((timings.size() < 2) ||
            (relDiff (timings[1], timings[0]) > thresholdNew)) {
            loops = min(maxLoopsNew, (int)loops + 1);
            if (maxLoopsNew == loops) {
                // If we don't converge the lowest timing is the most likely...
                outFile << "#Comparing to previous: Failed to converge." << endl;
            }
        }
    } else if (oldTiming == failedTestCode) { 
        // Previous timing was a failed test
        // Do nothing
    } else { // Compare to previous timing (is faster than creating new)
        if (relDiff(oldTiming, timings[0]) > thresholdExisting) {
            loops = min (maxLoopsExisting, loops + 1);
            if (maxLoopsExisting == loops) {
                // If we don't converge the lowest timing is the most likely...
                outFile << "#No previous timing: Failed to converge." << endl;
            }
        }
    }
}

void MicrostateAccounting::reportTime(void) {
    toggleMicrostateAccounting(0); 
    outFile << timings[0] << endl;      
    outFile.close();
}

void MicrostateAccounting::reportFailedTest(void) {
    outFile << failedTestCode << endl;
    toggleMicrostateAccounting(0);
    outFile.close();
}
#endif
