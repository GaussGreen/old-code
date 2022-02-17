//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CommandLineParams.cpp
//
//   Description : Processed reg test command line options
//
//   Author      : André Segger
//
//   Date        : 04 Jun 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Version.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Null.hpp"
#include "edginc/HTMLWriter.hpp"

DRLIB_BEGIN_NAMESPACE

// initialise static strings - note there is more to do than just adding a
// string here. It would be nice if you just created a class for each
// option encoding all the relevant information
// Search for EXPLICIT_ACTION to find out where to make changes
string CommandLineParams::EdgMini        = "edgmini";
string CommandLineParams::EdgMaxi        = "edgmaxi";
string CommandLineParams::Perm           = "perm";
string CommandLineParams::PermAlias      = "p";
string CommandLineParams::Pong           = "pong";
string CommandLineParams::InFiles        = "inputFiles";
string CommandLineParams::OutFiles       = "outputFiles";
string CommandLineParams::CurrentOutFile = "CurrentOutFile";
string CommandLineParams::Quantify       = "quantify";
string CommandLineParams::JIT            = "jit";
string CommandLineParams::NOJIT          = "nojit";
string CommandLineParams::Loop           = "loop";
#ifdef sun
string CommandLineParams::Time           = "time";
#endif
string CommandLineParams::AddSens        = "addsens";
string CommandLineParams::Unparp         = "unparp";
string CommandLineParams::Canon          = "canon";
string CommandLineParams::FewIter        = "fewIter";
string CommandLineParams::AddReq         = "addreq";
string CommandLineParams::Help           = "h";
string CommandLineParams::Mega           = "mega";
string CommandLineParams::MegaShuff      = "megaShuffle";
string CommandLineParams::ShuffSeed      = "shuffleSeed";
string CommandLineParams::Wrapper        = "wrapper";
string CommandLineParams::Doc            = "doc";
string CommandLineParams::DRI            = "dri";
string CommandLineParams::DRI_WRAP       = "driwrap";
string CommandLineParams::DebugPackets   = "debugPackets";
string CommandLineParams::DLLPath        = "dllPath";
string CommandLineParams::UseStateVars   = "useStateVars";

StringSet   CommandLineParams::allOptions;
string      CommandLineParams::exeName("");
HashtableSP CommandLineParams::activeOptions(new Hashtable());

/** add a command line option to the hash */
void CommandLineParams::addParameter(const string& option, IObjectSP value)
{
    CommandLineParams::activeOptions->put(option, value);
}

/** retrieve parameter for a given command line option 
    returns false if this option has not been set      */
IObjectSP CommandLineParams::getParameterArgs(const string& option)
{
    return CommandLineParams::activeOptions->get(option);
}

/** check whether command line option is has been set */
bool CommandLineParams::hasParameter(const string& option)
{
    return CommandLineParams::activeOptions->containsKey(option);
}

/** check whether command line option is has been set */
bool CommandLineParams::isValid(const string& option)
{
    StringSetIterator pos;
    pos = CommandLineParams::allOptions.find(option);
    return ( pos != CommandLineParams::allOptions.end());
}

/** remove prefix from command line option if it exists */
string CommandLineParams::stripPrefix(const string& option)
{
    static const string prefix("-");
    string normalisedOption = option;
    if (normalisedOption[0] == '-') {
        normalisedOption.erase(0,1);
    }
    return normalisedOption;
}

bool CommandLineParams::hasAdditionalParam(const string& option)
{
    // currently only allows 1 parameter - could possibly be extended
    // EXPLICIT_ACTION  - possibly
    if ( !option.compare(CommandLineParams::EdgMini) || 
         !option.compare(CommandLineParams::Pong)    ||
         !option.compare(CommandLineParams::Unparp)  ||
         !option.compare(CommandLineParams::Canon)   ||
         !option.compare(CommandLineParams::JIT)     ||
#ifdef sun
         !option.compare(CommandLineParams::Time)    ||
#endif
         !option.compare(CommandLineParams::NOJIT)   ||
         !option.compare(CommandLineParams::Mega)    ||
         !option.compare(CommandLineParams::MegaShuff)    ||
         !option.compare(CommandLineParams::Help)    ||
         !option.compare(CommandLineParams::DRI)     ||
         !option.compare(CommandLineParams::DRI_WRAP)||
         !option.compare(CommandLineParams::FewIter) ||
         !option.compare(CommandLineParams::UseStateVars) )
    {
        return false;
    }
    return true;
}

/** parse all input parameters. returns true if all options are 
    valid and could be parsed successfully */
bool CommandLineParams::parseOptions(int argc, char *argv[])
{
    static const string routine("CommandLineParams::parseOptions");
    string    currentOption;
    string    parameter;
    bool      foundFile  = false; // means have we got an input file to run
    bool      finished = false;

    // store executable name
    exeName = parsePathToBaseName(argv[0]);

    // parse optional parameters
    for (int i = 1; i < argc && !finished; i++) {
        string argument = argv[i];
        bool   isOption = !argument.empty() && argument[0] == '-';
        if (isOption){
            currentOption = CommandLineParams::stripPrefix(argument);
            // EXPLICIT_ACTION  - possibly
            if ( currentOption == CommandLineParams::PermAlias) {
                currentOption = CommandLineParams::Perm;
            }
            if ( CommandLineParams::isValid(currentOption) ) {
                // we know we have a valid option - just have to check for
                // additional parameters 
                IObjectSP value;
                if ( CommandLineParams::hasAdditionalParam(currentOption) ) {
                    if ( i == argc-1 ) {
                        throw ModelException(routine, "Command line switch "+
                                             argument+" requires additional "+
                                             "parameter");
                    }
                    value = CommandLineParams::parseParameter(
                        currentOption,
                        argv[++i]);
                } else {
                    value = CNull::create(); // add something
                }
                // finally add the parameter to the hashtable
                CommandLineParams::addParameter(currentOption, value);
            } else {
                throw ModelException(routine, "Unrecognised command line "
                                     "option: "+argument);
            }
        }
        else if (foundFile){
            throw ModelException(routine, "Parameter "+argument+" unexpected");
        } else {
            CStringArraySP inFiles (new StringArray());
            CStringArraySP outFiles(new StringArray());
            string outFileExtension;
            // the last parameter was not a valid option - we assume it's 
            // the input file name 
            if ( argument[0] == '@' )
            {
                string inFileName;
                string outFileName;
                // multi file - need to parse file
                if ( i  == argc-1 )
                {
                    throw ModelException(routine,
                                         "The input filename and output "
                                         "directory name must be provided.\n");
                }

                string dirName = argv[++i];

                if ( i  == argc-1 )
                {
                    outFileExtension = "out";
                }
                else
                {
                    outFileExtension = argv[++i] ;
                }

                // get rid of magic @
                argument.erase(0,1);

                FILE* inputFile = NULL;
                char line[255];

                if( (inputFile = fopen( argument.c_str(), "r" )) == NULL )
                {
                    throw ModelException(
                        routine,
                        "Could not open file " + argument +
                        "containing list of files to process");
                }

                // parse input file 
                while ( fscanf( inputFile, "%s", line ) != EOF )
                {
                    // read input file name
                    inFileName = line;

                    // create output file name
                    outFileName = dirName;

                    // add a directory slash if it's not there yet
                    if ( outFileName[outFileName.size()-1] != '/' )
                    {
                        outFileName += "/";
                    }
                    // so how do you do a replace with stl?
                    string myInFileName = inFileName;
                    for (size_t i = 0; i < myInFileName.size(); i++){
                        if (myInFileName[i] == '\\'){
                            myInFileName[i] = '/';
                        }
                    }
                    size_t lastSlash = myInFileName.rfind('/');
                    if (lastSlash == string::npos){
                        // and then append the file name
                        outFileName += inFileName;
                    } else {
                        outFileName += myInFileName.substr(lastSlash, 
                                                           string::npos);
                    }
                    size_t dotPos = outFileName.rfind('.');

                    if ( dotPos == string::npos )
                    {
                        outFileName = outFileName + "." + outFileExtension;
                    }
                    else
                    {
                        outFileName = outFileName.substr(0,dotPos) + "." + 
                                      outFileExtension;
                    }

                    inFiles->push_back(inFileName);
                    outFiles->push_back(outFileName);
                }

                if( fclose( inputFile ) )
                {
                    throw ModelException(
                        routine,
                        "Unable to close file " + argument);
                }
            }
            else
            {
                // single file
                string outFile;
                inFiles->push_back(argument);

                if ( i  < argc-1 )
                {
                    // check if there is an output file name
                    ++i;
                    outFile = CommandLineParams::stripPrefix(argv[i]);
                }
                else
                {
                    size_t dotPos = argument.find('.');
                    if ( dotPos == string::npos )
                    {
                        outFile = argument + ".out";
                    }
                    else
                    {
                        outFile = argument.substr(0,dotPos) + ".out";
                    }
                }
                outFiles->push_back(outFile);

                FILE* inputFile = NULL;
                if( (inputFile = fopen( argument.c_str(), "r" )) == NULL )
                {
                    // set output file before throwing an exception
                    addParameter(CurrentOutFile, CStringSP(CString::create(outFile)));
                    throw ModelException(routine,
                                         "Could not open file " + argument);
                }

                fclose(inputFile);
            }

            CommandLineParams::addParameter(CommandLineParams::InFiles,  inFiles);
            CommandLineParams::addParameter(CommandLineParams::OutFiles, outFiles);

            foundFile  = true;
            finished = true;
        }
    }

    // validate arguments
    // CommandLineParams::validateArgs();

    if ( !foundFile ) {
        cout << "Quantitative Research Regression Tester v" << 
                CVersion::DRLibVersion() << endl;
        cout << "Usage: models <infile> <outfile>" << endl;

        if (CommandLineParams::hasParameter(CommandLineParams::Help)){
            CommandLineParams::help();
        }
        if (CommandLineParams::hasParameter(CommandLineParams::Doc)){
            IObjectSP param = CommandLineParams::
                getParameterArgs(CommandLineParams::Doc);
            const string& classname = 
                (dynamic_cast<CString&>(*param)).stringValue();
            string filename = classname + ".html";

            HTMLWriter* html = HTMLWriter::htmlFile(filename);
            CClassConstSP clazz(CClass::forName(classname));
            html->document(clazz);
            cout << "Documented " << classname << " in " << filename << endl;
        }
        return false;
    }

    return true;
}

IObjectSP CommandLineParams::parseParameter(
            const string& option,
            char *arg)
{

    // EXPLICIT_ACTION  - possibly
    if ( !option.compare(CommandLineParams::Perm)      ||
         !option.compare(CommandLineParams::PermAlias) ||
         !option.compare(CommandLineParams::EdgMaxi)   ||
         !option.compare(CommandLineParams::Quantify)  ||
         !option.compare(CommandLineParams::Loop)      ||
         !option.compare(CommandLineParams::ShuffSeed))
    {
        int intParam;
        sscanf( arg, "%d", &intParam );
        return CIntSP(CInt::create(intParam));
    }
    else if (!option.compare(CommandLineParams::AddSens) ||
             !option.compare(CommandLineParams::AddReq)  ||
             !option.compare(CommandLineParams::Doc)     ||
             !option.compare(CommandLineParams::DebugPackets)     ||
             !option.compare(CommandLineParams::DLLPath) ||
             !option.compare(CommandLineParams::Wrapper)) {
        char strParam[128];
        sscanf( arg, "%s", strParam);
        return CStringSP(CString::create(strParam));
    }
    else
    {
        throw ModelException(
                "CommandLineParams::parseParameter", 
                "Option " + option + 
                " does not have an additional parameter");
    }
}

static char UNIX_PATH_SEP_CHAR = '/';
static char PC_PATH_SEP_CHAR   = '\\';
static char PC_DRIVE_SEP_CHAR  = ':';
static char EXTENSION_SEP_CHAR = '.';

/** strips out base name and returns executable name (less any .exe) */
string CommandLineParams::parsePathToBaseName(char *fullExeName) {
    static const string method = "CommandLineParams::parsePathToBaseName";
    try {
        char *baseName;  // final result 
        char *rref;      // readonly ptr into fullExeName 
        char *wref;      // write ptr into baseName 
        int   len = strlen(fullExeName);

        // extravagant but certainly big enough! (not forgetting 
        // string terminator) 
        baseName = new char[len+1];

        /* 2 stages: 
         * 1=read from end-to-start and stop at separator or start of string
         * 2=read from start-to-end and stop at period or end of string
         */
        for (rref=fullExeName+len-1; rref>=fullExeName; rref--) {
            if (*rref == UNIX_PATH_SEP_CHAR || 
                *rref == PC_PATH_SEP_CHAR ||
                *rref == PC_DRIVE_SEP_CHAR) {
                break;
            }
        }
        // skip over the path separator or reset on start of fullExeName
        rref++;
        // read off base name into baseName. rref initialised above ; if no
        // path separator then rref == fullExeName
        // Note may read up to the terminating '\0' 
        for (wref = baseName; rref<=fullExeName+len; rref++, wref++) {
            if (*rref == EXTENSION_SEP_CHAR) {
                // done - if no period then the '\0' terminator is
                // copied by the line below 
                *wref = '\0';
                break;
            }

            // pick out the base name char by char 
            *wref = *rref;
        }

        string baseExeName = baseName;

        delete[] baseName;

        return baseExeName;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** gives executable name (less any path and .exe) */
string CommandLineParams::executableName() {
    return exeName;
}

void CommandLineParams::help() {
    cout << "Command line options:" << endl;
    for (StringSetIterator pos = CommandLineParams::allOptions.begin();
         pos != CommandLineParams::allOptions.end();
         ++pos) {
        cout << *pos << (hasAdditionalParam(*pos) ? " <param>" : "") << endl;
    }
}



CommandLineParams::CommandLineParams(): CObject(TYPE){}

CommandLineParams::~CommandLineParams() {}

/** Add a command line parameter to the list of recognised options */
void CommandLineParams::addCommandLineParameter(const string& param){
    allOptions.insert(param);
}


/** Invoked when Class is 'loaded' */
void CommandLineParams::load(CClassSP& clazz){
    REGISTER(CommandLineParams, clazz);
    SUPERCLASS(CObject);
    /** add all known parameter options to the hash table */
    // EXPLICIT_ACTION
    addCommandLineParameter(EdgMini);
    addCommandLineParameter(EdgMaxi);
    addCommandLineParameter(Perm);
    addCommandLineParameter(Pong);
    addCommandLineParameter(Quantify);
    addCommandLineParameter(JIT);
    addCommandLineParameter(NOJIT);
    addCommandLineParameter(Loop);
#ifdef sun
    addCommandLineParameter(Time);
#endif
    addCommandLineParameter(AddSens);
    addCommandLineParameter(Unparp);
    addCommandLineParameter(Canon);
    addCommandLineParameter(FewIter);
    addCommandLineParameter(AddReq);
    addCommandLineParameter(Help);
    addCommandLineParameter(Mega);
    addCommandLineParameter(MegaShuff);
    addCommandLineParameter(ShuffSeed);
    addCommandLineParameter(Wrapper);
    addCommandLineParameter(Doc);
    addCommandLineParameter(DRI);
    addCommandLineParameter(DRI_WRAP);
    addCommandLineParameter(DebugPackets);
    addCommandLineParameter(DLLPath);
    addCommandLineParameter(UseStateVars);
}

IObject* CommandLineParams::defaultCLParams(){
    return new CommandLineParams();
}

CClassConstSP const CommandLineParams::TYPE = 
CClass::registerClassLoadMethod(
    "CommandLineParams", typeid(CommandLineParams), load);

DRLIB_END_NAMESPACE
