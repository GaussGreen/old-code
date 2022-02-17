//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CommandLineParams.hpp
//
//   Description : Processed reg test command line options
//
//   Author      : André Segger
//
//   Date        : 04 Jun 2001
//
//----------------------------------------------------------------------------

#ifndef EDG_CL_PARAMS_H
#define EDG_CL_PARAMS_H
#include <set>
#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/Hashtable.hpp"

DRLIB_BEGIN_NAMESPACE

typedef std::set<string>           StringSet;
typedef std::set<string>::iterator StringSetIterator;

/** a hashtable which stores all command line options */
class TOOLKIT_DLL CommandLineParams : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class CommandLineParamsHelper;

    static string EdgMini;
    static string EdgMaxi;
    static string Perm;
    static string PermAlias;
    static string Pong;
    static string InFiles;
    static string OutFiles;
    static string CurrentOutFile;
    static string Quantify; // make quantify only measure certain bits of code
    static string NOJIT;    // disable just in time debugging on NT
                            // (disabled by default in overnight/produc build)
    static string JIT;      // enable just in time debugging on NT
                            // (enabled by default in non-overnight build)
    static string Loop;     // forces pricing to be done n times
#ifdef sun
    static string Time;     // Times execution using microstate accounting and saves to file (Solaris only) 
#endif
    static string AddSens;  // add a greek to the control
    static string Unparp;   // rebuild input file
    static string Canon;    // rebuild input file without back-references & exit
    static string FewIter;  // for Monte Carlo - reduce the iterations
    static string AddReq;   // add a request to the control
    static string Help;     // echo help to stdout
	static string Mega;     // use MegaControl
	static string MegaShuff;// variant of MegaControl: shuffle greeks, compare
    static string ShuffSeed;// seed to use to override seed based on today's date.
	static string Wrapper;  // wrapper --> object file conversion
    static string Doc;      // create HTML doc for a class
    static string DRI;      // rebuild inputs via Global DR Interface
    static string DRI_WRAP;   // treat everything as an XObject
    static string DebugPackets; // report debug outputs in .out summary
    static string DLLPath; /* where to find dll/so which are dynamically opened.
                              Split paths by either colons (unix) or 
                              semi-colons (windows) */
    static string UseStateVars;  // for State Var - overrides model param
     
    /** add a command line option to the hash */
    static void addParameter(const string& option, IObjectSP value);

    /** retrieve parameter for a given command line option 
        returns empty sp if this option has not been set      */
    static IObjectSP getParameterArgs(const string& option);

    /** check whether command line option has been set */
    static bool hasParameter(const string& option);

    /** check whether a given option is valid */
    static bool isValid(const string& option);

    /** parse all input parameters */
    static bool parseOptions(int argc, char *argv[]);

    static IObjectSP parseParameter(const string& option, char *arg);

    /** destructor */
    virtual ~CommandLineParams();

    static bool hasAdditionalParam(const string& option);

    /** gives executable name (less any path and .exe) */
    static string executableName();

    /** echo help to stdout */
    static void help();

    /** Add a command line parameter to the list of recognised options */
    static void addCommandLineParameter(const string& param);
private:
    /** a hashtable of all valid command line options */
    static StringSet   allOptions;

    /** name of the executable */
    static string exeName;

    /** the options passed into the current test */
    static HashtableSP activeOptions;

    /** default constructor */
    CommandLineParams();

    /** remove prefix from command line option if it exists */
    static string stripPrefix(const string& option);

    /** strips out base name and returns executable name (less any .exe) */
    static string parsePathToBaseName(char *fullExeName);
        
    static IObject* defaultCLParams();
    static void load(CClassSP& clazz);

};

DRLIB_END_NAMESPACE

#endif
