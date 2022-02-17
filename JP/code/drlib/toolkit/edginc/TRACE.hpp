/**
 * @file TRACE.hpp
 *
 * A very simple-minded mechanism for printing indented traces of program
 * execution to stderr or a file.  See RiskQuantityEvaluator and various
 * other files in riskmgr for examples.
 *
 *    -  TRACE("hello " << 2+2) writes "hello 4"
 *    -  TRACE_SHOW(2+2) writes "2+2 = 4"
 *    -  TRACE_BLOCK("price the options") writes "prices the options"
 *       and increases the indent level from here to the end of the
 *       current C++ block
 *    -  TRACE_METHOD writes "MyClass::mymethod()" and increases the
 *       indent level from here to the end of the current C++ block
 *
 * Tracing is only compiled in if DEBUG is defined, i.e. it's not included in
 * the optimised production code.
 *
 * Even in the debug build it's disabled by default, depending on the
 * value of the QLIB_TRACE environment variable:
 *
 *    -  QLIB_TRACE unset: no tracing
 *    -  QLIB_TRACE empty (or just white space): tracing to stderr
 *    -  QLIB_TRACE set and nonempty: tracing to the file it names
 *
 * Example output:
 * 
 * <PRE>
 *       ...
 *       Recursing into hypothetical world
 *       ---- HypothesisTree::evaluate
 *          We're in a world subjected to 0.002000 tweak of 25301290 VolParallel
 *          Now we price the instrument(s) in this world state
 *          and store the prices in buffers
 *             ---- HypothesisTree::QuantityToCompute::compute
 *                ---- PriceResultsFunction::operator`()'
 *                   results->retrievePrice() = 13.7859
 *                Store that number, for future calculation of VEGA_PARALLEL 25301290
 *       Applying hypothesis "0.001000 tweak of 25301290 Spot"
 *       Recursing into hypothetical world
 *       ---- HypothesisTree::evaluate
 *          We're in a world subjected to 0.001000 tweak of 25301290 Spot
 *          Now we price the instrument(s) in this world state
 *          and store the prices in buffers
 *             ---- HypothesisTree::QuantityToCompute::compute
 *                ---- PriceResultsFunction::operator`()'
 *                   results->retrievePrice() = 13.7288
 *                Store that number, for future calculation of DELTA 25301290
 *                Store that number, for future calculation of GAMMA 25301290
 *       ...
 * </PRE>
 */

#ifndef QLIB_TRACE_H
#define QLIB_TRACE_H

DRLIB_BEGIN_NAMESPACE

#ifdef DEBUG

// 
// -----------
//  Internals
// -----------
// 

// True if the environment variable QLIB_TRACE is set

TOOLKIT_DLL extern bool TRACE_enabled;

// Returns selected output stream having written indent prefix

TOOLKIT_DLL ostream& TRACE_stream();

TOOLKIT_DLL void TRACE_setStream(ostream*);

// Constructor indents, destructor unindents

struct TOOLKIT_DLL TRACE_Block {
    TRACE_Block();
    TRACE_Block(const char* functionName);
    ~TRACE_Block();
};

// Visual C++ turns out to be really bad at token pasting
// See http://msdn2.microsoft.com/en-us/library/b0084kay.aspx

#   define __TRACE_BLOCK(X, L) TRACE(X); TRACE_Block __TRACE_dummy##L
#   define _TRACE_BLOCK(X, L) __TRACE_BLOCK(X, L) // [sic]

#   define __TRACE_METHOD(L) TRACE_Block  __TRACE_dummy_##L(__FUNCTION__)
#   define _TRACE_METHOD(L) __TRACE_METHOD(L)

// 
// ----------------------
//  Programmer interface
// ----------------------
// 

/**
 * Write a message comprising a string or a series of C++ values to the
 * trace
 *
 * E.g. TRACE("Two plus two is " << 2+2 << " so there");
 *
 * You get to write one line only, with no newlines please.  It will be
 * indented by three spaces times the level of TRACE_BLOCK / TRACE_METHOD
 * nesting.
 */

#   define TRACE(X)                                \
        do {                                       \
            if (TRACE_enabled) {                   \
                TRACE_stream() << X << "\n";       \
            }                                      \
        }                                          \
        while (0)

/**
 * Write the value of an expression to the trace,
 * labelled with the code for the expression
 *
 * E.g. TRACE_SHOW(2+2) writes "2+2 = 4"
 */

#   define TRACE_SHOW(X) TRACE(#X << " = " << (X))

/**
 * Write a message to the trace and increase the indent level for
 * subsequent messages, from here to the end of the current C++ block.
 */

#   define TRACE_BLOCK(X) _TRACE_BLOCK(X, __LINE__)

/**
 * Write the name of the current method to the trace, and increase the
 * indent level for subsequent messages, from here to the end of the
 * current C++ block
 */

#   define TRACE_METHOD _TRACE_METHOD(__LINE__)

#else

#   define TRACE_enabled    (false)
#   define TRACE(X)         do {} while (0)
#   define TRACE_SHOW(X)    do {} while (0)
#   define TRACE_BLOCK(X)   do {} while (0)
#   define TRACE_METHOD     do {} while (0)

#endif

DRLIB_END_NAMESPACE

#endif
