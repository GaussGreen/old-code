//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : models.cpp
//
//   Description : regression tester
//
//   Author      : Andrew J Swain
//
//   Date        : 19 January 2001
//
//
//----------------------------------------------------------------------------
#ifdef sun
#undef __STDC__
#define __STDC__ 0
#endif

#include "edginc/config.hpp"
#include <signal.h>
#include <time.h>
#include "edginc/XLAddin.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Library.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/AddinLib.hpp"
#include "edginc/Addin.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/Format.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/XObject.hpp"
#include "edginc/RegressionTestMode.hpp"
#include "edginc/EDGServices.h"
#include "edginc/ToolkitDebug.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/MegaShuffleInterface.hpp"

#ifdef sun
#include "edginc/MicrostateAccounting.hpp"
#endif 

#ifdef USING_VALGRIND_REGTESTER
#include "../../../3rdParty/valgrind-3.0.1/include/valgrind.h"
#endif

#ifdef PURIFY
extern "C"
{
    extern int purify_new_leaks();
    extern int purify_printf(char* fmt, ...);
}
#define PURIFY_NEW_LEAKS purify_new_leaks()
#else
#define PURIFY_NEW_LEAKS
#endif
#ifdef QUANTIFY
extern "C"
{
    int quantify_start_recording_data();
    int quantify_stop_recording_data();
}
#define QUANTIFY_STOP_RECORDING_DATA quantify_stop_recording_data()
#define QUANTIFY_START_RECORDING_DATA quantify_start_recording_data()
#else
#define QUANTIFY_STOP_RECORDING_DATA
#define QUANTIFY_START_RECORDING_DATA
#endif

USING_DRLIB_NAMESPACE
using namespace std;

// define SRM3 to link against srm3wrapper.cpp
#define xSRM3
// too lazy to have separate header files
extern int edgminiMain(void);
extern int edgmaxiMain (int argc, char* argv[]);
extern int driminiMain(void);
extern int drimaxiMain (int argc, char* argv[]);
#ifdef SRM3
extern int srm3Main (int argc, char* argv[]);
#endif
extern void runViaDRI(IObjectSP input, const string& outfile);


static void runTest(int& permNum, IObjectSP input, OutputFile& output);
static void parseAndRunInputFiles();

static const string MULTIPERM = "MULTIPERM";

/* signal handler for overnight on NT (avoid the just in time debugging) */
#if defined(WIN32)
static void signalHandler(int signal){
    exit(signal); /* not allowed to call printf etc, so just give up */
}
#endif

static int quantifyTarget = 0; // drives what to quantify

int main(int argc, char *argv[]){
    static const string routine = "main";
    int exitCode = 1; // failure
    try {
        RegressionTestMode::on(); // flag we are in regression test mode
        // ensure all symbols are linked in
        AddinLib::linkInAllLibs();
        Library::startup();
        // what we need to do for regression tester (no Excel around)
        Addin::registerKit(XLAddin::defaultRegister);

        // easier to deal with edgmaxi and edgmini separately since they
        // have their own set of command line options
        if (argc > 1 && (strcmp(argv[1], "-edgmini") == 0 ||
                         strcmp(argv[1], "-edgmaxi") == 0 ||
                         strcmp(argv[1], "-srm3") == 0 ||
                         strcmp(argv[1], "-drimini") == 0 ||
                         strcmp(argv[1], "-drimaxi") == 0)){
            if (strcmp(argv[1], "-edgmini") == 0) {
                return edgminiMain();
            }
            else if (strcmp(argv[1], "-drimini") == 0) {
                return driminiMain();
            }
            else if (strcmp(argv[1], "-srm3") == 0) {
#ifdef SRM3
                return srm3Main(argc, argv);
#else
                throw ModelException("main", "Recompile this file with SRM3"
                                     "defined to get SRM3 linked in");
#endif
            }
            else {
                bool edr = !strcmp(argv[1], "-edgmaxi");
                // remove first option from inputs
                argc--;
                for (int i = 1; i < argc; i++){
                    argv[i] = argv[i+1];
                }
                return edr ? edgmaxiMain(argc, argv) : drimaxiMain(argc, argv);
            }
        }
        if ( !CommandLineParams::parseOptions(argc, argv) ) {
            return 0;
        }

#       if defined (WIN32)
            if (
#               if defined (DOING_OFF_BUILD)
                    !CommandLineParams::hasParameter(CommandLineParams::JIT)
#               else
                    CommandLineParams::hasParameter(CommandLineParams::NOJIT)
#               endif
               ) {
                if (signal(SIGFPE,  signalHandler)  == SIG_ERR ||
                    signal(SIGILL,  signalHandler)  == SIG_ERR ||
                    signal(SIGSEGV, signalHandler)  == SIG_ERR) {
                    throw ModelException("main", "Failed to set signal handler");
                }
            }
#       endif

        IObjectSP quantifyObj =
            CommandLineParams::getParameterArgs(CommandLineParams::Quantify);
        quantifyTarget = !quantifyObj?
            0: dynamic_cast<CInt&>(*quantifyObj).intValue();
       if (quantifyTarget == 0){
           QUANTIFY_START_RECORDING_DATA;
       }

       parseAndRunInputFiles();
       exitCode = 0; // success
    } catch (ModelException& e){
        // print to the screen - safer as it's not clear whether we've got
        // a valid output file
        e.printStackTrace();
    } catch (exception& e){
        cout << "Caught exception: " << e.what() << endl;
        cout << "main failed" << endl;
    }
    Library::shutdown();
    return exitCode;
}

static string recreateInputFile(const string& infile, bool allowBackReferences) {
    static const string method = "recreateInputFile";

    try {
        // where to write to
        string newfile = infile + ".new";

        XMLWriter xmlOut(newfile);
        xmlOut.allowBackReferences(allowBackReferences);

        // get the root of the document
        XMLReader  xml(infile, true);
        Reader::NodeSP root(xml.root());

        /* if the root node has the 'magic' tag MULTIPERM then
           this is several tests stuck together, so we can't
           read it in & process it directly */
        string rootname = root->name();

        if (rootname == MULTIPERM) {
            xmlOut.objectStart(MULTIPERM, "", 0, true);
            // tests should be children of the root
            Reader::NodeListSP nl(root->children());

            // walk through all root children
            for (unsigned int i = 0; i < nl->size(); i++) {
                Reader::Node* child = (*nl)[i].get();
                IObjectSP obj(xml.read(child));
                obj->write(child->name(), &xmlOut);
            }

            xmlOut.objectEnd(MULTIPERM, 0);
        }
        else {
            // a single test, read it in directly
            IObjectSP obj(xml.read(root.get()));
            obj->write(root->name(), &xmlOut);
        }

        return newfile;
    }
    catch (ModelException& e){
        throw ModelException(e, method);
    }
}

//// Allow EDRInterface access to these functions. A bit hacky but saves
//// a lot of work and (in theory....) the need for these will go soon
// convert a DRObject to an IObject
extern QLIB_DLL_IMPORT IObject* drObject2Object(DRObject drObj);
// convert an IObject to a DRObject
extern QLIB_DLL_IMPORT DRObject object2DRObject(IObject* object,
                                                DRService *svc);
static DRService* svc = 0; // bit hacky

// To be able to use Quantify to tell us the difference between
// the 1st loop iteration and the following ones, we neeed to
// introduce one level of indirection. The result is the same,
// but now we call an intermediate function for iterations after
// the first one. The var args stuff is just a cute way of stopping
// the compiler optimising this function away.
IObjectSP wrapperAroundRunTest (ClientRunnable* runnable, ...) {
    va_list param_pt;
    va_start(param_pt, runnable);

    IRegressionTest* test = va_arg(param_pt, IRegressionTest*);
    va_end(param_pt);

    return  (!test) ? runnable->run(): test->runTest();
}

static void runTest(int& permNum, IObjectSP input, OutputFile& output) {
    static const string method = "runTest";

    try {
        // print the magic dots ...
        permNum++;
        cout << "." << flush;
        if (!(permNum % 20)){
            cout << Format::toString(permNum) << flush;
        }
        IObjectSP results;

        // run the regression test, be it addin, risk manager, or other
        IRegressionTest* test = dynamic_cast<IRegressionTest*>(input.get());
        ClientRunnable* runnable = 0;
        if (!test){
            runnable = dynamic_cast<ClientRunnable*>(input.get());
            if (!runnable){
                throw ModelException(method, "Don't know how to execute "
                                     "regression test for object of type "+
                                     input->getClass()->getName());
            }
        }

        // if megaShuffle has been specified on the commandline, 
        // we'll do something special 
        MegaShuffleInterfaceSP megaShuffleTest;
        if (CommandLineParams::hasParameter(CommandLineParams::MegaShuff)) {

            //Either use a seed specified on the commandline, 
            //or create a seed from today's date (not including the time).
            int seed = MegaShuffleInterface::createShuffleSeed();
            if (CommandLineParams::hasParameter(CommandLineParams::ShuffSeed)) { 
                IObjectSP obj = CommandLineParams::getParameterArgs(CommandLineParams::ShuffSeed);
                seed = dynamic_cast<CInt&>(*obj).intValue();
            }
            megaShuffleTest = MegaShuffleInterfaceSP(new MegaShuffleInterface(input, seed));
            test = megaShuffleTest.get();
        }

        // command line can force multiple pricings (v poor man's quantify)
        IObjectSP loopObj =
            CommandLineParams::getParameterArgs(CommandLineParams::Loop);
        int       loop = 1;
        if (loopObj.get()){
            loop = dynamic_cast<CInt&>(*loopObj).intValue();
        }
        // save exceptions in a string array
        vector<string> errorStack;
        ErrorHandler origErrorHandler(
            ErrorHandler::set(ErrorHandler(errorStack), false));

        // for timings
        clock_t startLoop = clock();
        if (quantifyTarget & 1){
            QUANTIFY_START_RECORDING_DATA;
        }

#ifdef sun
        MicrostateAccounting* myAcc; 
        if (CommandLineParams::hasParameter(CommandLineParams::Time)){
          myAcc = new MicrostateAccounting(output.getFilename(),permNum); 
        }
#endif
        bool viaDRI =
            CommandLineParams::hasParameter(CommandLineParams::DRI_WRAP);
        DRObject drObj = 0;
        try{
            for (int i = 0; i < loop; i++){
#ifdef sun
                if (CommandLineParams::hasParameter(CommandLineParams::Time)){ 
                    myAcc->recordStartTime();
                }
#endif 
                if (viaDRI){
                    // turn IObject into DRObject
                    drObj = object2DRObject(input.get(), svc);
                    // turn into XObject
                    XObjectSP xObj(XObject::toXObject(drObj,
                                                      false /*don't own */));
                    // clone
                    xObj = xObj->recursiveClone();
                    // then execute
                    IDRObjectSP result(xObj->execute());
                    // then recover our object
                    XObjectSP theResult(XObjectSP::dynamicCast(result));
                    theResult = theResult->recursiveClone();// extra testing
                    results =
                        IObjectSP(drObject2Object(theResult->getDRObject()));
                } else {
                    // To be able to use Quantify to tell us the difference between
                    // the 1st loop iteration and the following ones, we neeed to
                    // introduce one level of indirection. The result is the same,
                    // but now we call an intermediate function for iterations after
                    // the first one.
                    if (i == 0) {
                        results = !test ? runnable->run(): test->runTest();
                    }
                    else {
                        results = wrapperAroundRunTest(runnable, test);
                    }
                }
#ifdef sun 
                if (CommandLineParams::hasParameter(CommandLineParams::Time) ){
                   myAcc->recordStopTime(loop);
                }
#endif
            }
#ifdef sun 
            //Solaris microstate accounting
            if (CommandLineParams::hasParameter(CommandLineParams::Time) ){
               myAcc->reportTime();
            }
#endif
        } catch (exception& e){
#ifdef sun 
            if (CommandLineParams::hasParameter(CommandLineParams::Time) ){
                myAcc->reportFailedTest();
            }
#endif
            ModelException modelEx(e);
            modelEx.errorLog();
        }
        if (viaDRI && drObj){
            DRError err = 0;
            if ((err = (*svc->fptr->objectFree)(svc, drObj))) {
                (*svc->fptr->stringFree)(svc, err);
            }
            drObj = 0;
        }

        if (quantifyTarget & 1){
            QUANTIFY_STOP_RECORDING_DATA;
        }
        clock_t endLoop = clock();
        if (loopObj.get()){
            cout << "Average time (" << loop << " runs): " <<
                (double)(endLoop-startLoop)*1000.0/
                (double)(CLOCKS_PER_SEC * loop) << "ms" << endl;
        }
        // restore error handler
        ErrorHandler::set(origErrorHandler, true);

        // write results (including error stack) to output file
        output.write(results.get(), CStringArray(errorStack.begin(), errorStack.end()));
    }
    catch (ModelException& e){
        throw ModelException(e, method);
    }
}
//// can the specified file be opened for reading
static bool fileIsReadable(const string& infile){
    FILE* file = fopen(infile.c_str(), "r");
    if (!file){
        return false;
    }
    fclose(file);
    return true;
}

static void parseAndRunInputFiles() {
    static const string method = "parseAndRunInputFiles";

    const char *err;

    try{
        if (CommandLineParams::hasParameter(CommandLineParams::DRI_WRAP)){
            // fire up DRService
            DRServiceInitArgs args;
            EDGSInitServiceCreateArgs(&args);
            if (!DRCreateService(&args, &svc)) {
                printf("Failed to fire up DRService");
                EDGSGetServiceInvocationError(&err);
                throw ModelException(method, err);
            }
        }
        string       rootname;
        IObjectSP    objInFile  = CommandLineParams::getParameterArgs(
                                  CommandLineParams::InFiles);
        IObjectSP    objOutFile = CommandLineParams::getParameterArgs(
                                  CommandLineParams::OutFiles);
        string       infile;

        CStringArraySP inFilesSP  = CStringArraySP::dynamicCast(objInFile);
        CStringArraySP outFilesSP = CStringArraySP::dynamicCast(objOutFile);

        CStringArray* inFiles  = inFilesSP.get();
        CStringArray* outFiles = outFilesSP.get();

        if ( !inFiles) {
            throw ModelException(method, "Input file name is NULL");
        }

        if ( !outFiles) {
            throw ModelException(method, "Output file name is NULL");
        }

        for ( int fileNo=0 ;  fileNo<inFiles->size() ; ++fileNo )
        {
            // always clear out any existing output file. Don't move this
            // line.
            auto_ptr<OutputFile> output(new OutputFile((*outFiles)[fileNo]));
            try {
                IObjectSP    test;
                int permNum = 0; // just used for magic dots eg ...20...
                infile  = (*inFiles)[fileNo];
#ifdef PURIFY
                purify_printf("*** Processing %s ***\n", infile.c_str());
#endif
#ifdef USING_VALGRIND_REGTESTER
                VALGRIND_PRINTF("<inputfile filename=\"%s\"/>\n", infile.c_str());
#endif
                cout << infile.c_str() << ":" << flush;
                CStringSP outFile(CString::create((*outFiles)[fileNo]));
                CommandLineParams::addParameter(
                    CommandLineParams::CurrentOutFile, outFile);
                if (!fileIsReadable(infile)){
                    throw ModelException(method, "Unable to open/read file "+
                                         infile);
                }
                if (CommandLineParams::hasParameter(CommandLineParams::Canon)) {
                    cout << (" written to " + recreateInputFile(infile, false))
                         << endl;
                    continue;
                }
                if (CommandLineParams::hasParameter(CommandLineParams::Unparp)){
                    recreateInputFile(infile, true);
                }
                if (quantifyTarget & 4){
                    QUANTIFY_START_RECORDING_DATA;
                }
                // get the root of the document
                XMLReader  xml(infile, true);
                Reader::NodeSP root(xml.root());
                /* if the root node has the 'magic' tag MULTIPERM then
                   this is several tests stuck together, so we can't
                   read it in & process it directly */
                if (quantifyTarget & 4){
                    QUANTIFY_STOP_RECORDING_DATA;
                }
                rootname = root->name();
                if (CommandLineParams::hasParameter(CommandLineParams::DRI)){
                    output.reset(); // let go of the output file
                    // a single test, read it in directly & run it via DRI
                    test = IObjectSP(xml.read(root.get()));
                    runViaDRI(test, (*outFiles)[fileNo]);
                } else if (rootname == MULTIPERM) {
                    // tests should be children of the root
                    Reader::NodeListSP nl(root->children());
                    // check if we want to run a specific permutation
                    IObjectSP permObj;
                    permObj = CommandLineParams::getParameterArgs(
                        CommandLineParams::Perm);
                    if ( !(!permObj) ) {
                        // run a specific testcase
                        int permutation =
                            dynamic_cast<CInt&>(*permObj).intValue();
                        unsigned int i            = 0;
                        unsigned int numNodes     = nl->size();
                        unsigned int numElemNodes = 0;
                        bool found                = false;
                        while ( i < numNodes && !found ) {
                            Reader::Node* child = (*nl)[i++].get();
                            ++numElemNodes;
                            if (numElemNodes == (unsigned int) permutation){
                                found = true;
                                // run testcase
                                try {
                                    if (quantifyTarget & 2){
                                        QUANTIFY_START_RECORDING_DATA;
                                    }
                                    test = IObjectSP(xml.read(child));
                                    if (quantifyTarget & 2){
                                        QUANTIFY_STOP_RECORDING_DATA;
                                    }
                                    runTest(permNum, test, *output);
                                }
                                catch (ModelException& e){
                                    string m("continuing with next "
                                             "permutation (if any) ...");
                                    ModelException x(e, method, m);
                                    output->write(x);
                                }
                            } else {
                                CStringArray exceptions(0);
                                CStringSP message(
                                    CString::create("skipped ..."));
                                output->write(message.get(), exceptions);
                            }
                        }
                        if (!found) {
                            string m("Could not find permutation "+
                                     Format::toString(permutation)+
                                     "in file");
                            throw ModelException(method, m);
                        }
                    } else {
                        // walk through all root children
                        for (unsigned int i = 0; i < nl->size(); i++) {
                            Reader::Node* child = (*nl)[i].get();
                            try {
                                if (quantifyTarget & 2){
                                    QUANTIFY_START_RECORDING_DATA;
                                }
                                test = IObjectSP(xml.read(child));
                                if (quantifyTarget & 2){
                                    QUANTIFY_STOP_RECORDING_DATA;
                                }
                                runTest(permNum, test, *output);
                            }
                            catch (ModelException& e){
                                ModelException x(e, method,
                                                 "continuing with next "
                                                 "permutation (if any) ...");
                                output->write(x);
                            }
                        }
                    }
                } else {
                    // a single test, read it in directly & run it
                    test = IObjectSP(xml.read(root.get()));
                    runTest(permNum, test, *output);
                }
                cout << endl << flush;
            } catch (ModelException& e) {
                // write errors to output file
                output->write(e);
            }
            PURIFY_NEW_LEAKS;
//#ifdef USING_VALGRIND_REGTESTER
//          VALGRIND_PRINTF("</inputfile>\n");
//#endif
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }

    if (svc && (err = svc->fptr->serviceFree(svc))){
        svc->fptr->stringFree(svc, err);
        throw ModelException(method, "Failed to shut DRService down");
    }
}
