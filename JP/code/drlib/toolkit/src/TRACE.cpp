/**
 * @file TRACE.cpp
 */

#include "edginc/config.hpp"
#include <fstream>
#include "edginc/TRACE.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

#ifdef DEBUG

bool TRACE_enabled = false;
static ostream* os = &cerr;

struct TRACE_Init {
    TRACE_Init() {
        TRACE_enabled = false;
        os = &cerr;
        const char* f = getenv("QLIB_TRACE");

        if (f) {
            const char* nonwhite = f;
            while (isspace(*nonwhite)) ++nonwhite;

            if (*nonwhite) {
                static ofstream fs;
                fs.open(f);
                if (!fs) {
                    cerr << "Warning: QLIB_TRACE is set to \"" << f <<
                            "\" but I can't open that file for writing\n";
                }
                else {
                    os = &fs;
                    TRACE_enabled = true;
                }
            }
            else {
                os = &cerr;
                TRACE_enabled = true;
            }
        }
    }
};

static TRACE_Init dummy;

static int indentLevel = 0;

static void indent() {
    for (int i = 0; i < indentLevel; ++i) {
        *os << "   ";
    }
    (*os).flush();
}

ostream& TRACE_stream() {
    indent();
    return *os;
}

void TRACE_setStream(ostream* newos) {
    os = newos ? newos : &cerr;
    TRACE_enabled = newos ? true : false;
}

TRACE_Block::TRACE_Block() {
    ++indentLevel;
}

TRACE_Block::~TRACE_Block() {
    --indentLevel;
}

TRACE_Block::TRACE_Block(const char* functionName) {
    if (TRACE_enabled) {
        if (!strncmp(functionName, "drlib::", 7)) {
            functionName += 7;
        }
        indent();
        *os << "---- " << functionName << "\n";
        (*os).flush();
        ++indentLevel;
    }
}

#endif

DRLIB_END_NAMESPACE
