#include <iostream>
#include <esl_error.h>
#include <esl_log.h>


void testLog(int argc, char* argv[])
{
    // trap log mesages that have been printed to the default log stream
    StringLogger acc(eslLog);

    // make the file logger associated with the default stream
    FileLogger fl(eslLog, "./qt.log");

    // make the stream logger so we can see messages on the screen
    StreamLogger sl(eslLog, std::cerr);

    // set log level - anything at or above this level gets logged
    eslLog.setLevel(Log::info);

    // log some messages to the default stream
    ESL_LOG(Log::error) << "This is error on line " << __LINE__ << std::endl;
    ESL_LOG(Log::warn)  << "This is warning" << std::endl;

    // log some messages to the specific stream
    LOG_MSG(Log::warn, eslLog) << "This is warning on line " << __LINE__ << std::endl;

    ESL_LOG(Log::debug) << "This is debug - you should not see it" << std::endl;

    // log some messages using old scheme
    DR_Error("This is DR Error reported from '%s', line %d",__FILE__, __LINE__); 
    DR_Error("This is another DR Error"); 
    
    // test exception throwing
    ESL_LOG(Log::info) << "About to throw" << std::endl;

    try
    {
        throw EslException("main") << " Throwing esl exception from" 
                                   << " file: " << __FILE__ 
                                   << " line: " << __LINE__;

        ESL_LOG(Log::info) << "Exception thrown - you should not see this" << std::endl;
    }
    catch (std::exception& e)
    {
        std::cerr << "got exception: " << e.what() 
                  << "\nLog messages:\n" << acc.str() << std::endl;
    }
}





