//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RegressionTestMode.hpp
//
//   Description : sets up environment used by regression tester
//
//   Author      : Mark A Robson
//
//   Date        : 26 Oct 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_REGRESSIONTESTMODE_HPP
#define QLIB_REGRESSIONTESTMODE_HPP

DRLIB_BEGIN_NAMESPACE

/** sets up environment used by regression tester. Note that although
    the use of global statics should be discouraged in general, if used
    wisely can be beneficial. One such use is to allow extra checks when
    running in regression test mode -- this is not an alternative to using
    assert etc but allows say a replacement fast methodology to be testing
    against the original methodology. */
class TOOLKIT_DLL RegressionTestMode {
public:
    /** turns it on */
    static void on();

    /** turns it off */
    static void off();

    /** is it on or off ? */
    static bool isOn();

private:
    RegressionTestMode(); // not implemented
    static bool inRegressionTester;
};

DRLIB_END_NAMESPACE

#endif




