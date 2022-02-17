#include <iostream>
#include <esl_error.h>
#include <esl_log.h>


void testLog(int argc, char* argv[]);

void usage(char const* program)
{
    std::cerr << "USAGE: " << program << " "
              << "log"
              << std::endl;
}


int main(int argc, char* argv[])
{

    if (argc < 2)
    {
        usage(argv[0]);
        return -1;
    }

    try
    {
        std::string module(argv[1]);
        if (module == "log")
        {
            testLog(argc - 1, &argv[1]);
        }
        else
        {
            usage(argv[0]);
            return -1;
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    return 0;
}





