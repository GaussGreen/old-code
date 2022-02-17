// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 1/12/00 Bruce Broder
//
// $Header$
//

#if ! defined(_CM_GENERAL_GENERAL_)
#include "../General/General.h"
#endif
#if ! defined(_CM_GENERAL_TEST_MANAGER_)
#include "../General/TestManager.h"
#endif

using namespace std;
using namespace CM;

//
//   o b j e c t s
//

static const String theProgramName =
    "CMTest";

static const String theBanner =
    "\n" + String( theProgramName ) + " " + CMVersion() + "\n";

namespace CM {

//
//   t e s t   i n i t i a l i z a t i o n s
//

CM_TEST_DESCRIPTOR_INIT( environment, EnvironmentObjects, 1 )
CM_TEST_DESCRIPTOR_INIT( market, MarketObjects, 2 )
CM_TEST_DESCRIPTOR_INIT( pricer, PricerObjects, 3 )
CM_TEST_DESCRIPTOR_INIT( recovery, RecoveryObjects, 4 )

//
//   f u l l   t e s t s
//

CM_TEST_DESCRIPTOR_FULL(
    basket,
        BasketObjects, 10, BasketTest, true )

CM_TEST_DESCRIPTOR_FULL(
    tranchedBasket,
        TranchedBasketObjects, 10, TranchedBasketTest, true )

CM_TEST_DESCRIPTOR_FULL(
    bootstrapSpreads,
        ParCDSToSpreadObjects, 10, ParCDSToSpreadTest, true )

CM_TEST_DESCRIPTOR_FULL(
    cds,
        CDSObjects, 10, CDSTest, true )

CM_TEST_DESCRIPTOR_FULL(
    dayCount,
        DayCountObjects, 10, DayCountTest, true )

CM_TEST_DESCRIPTOR_FULL(
    badDayAdjustment,
        BadDayAdjustmentObjects, 10, BadDayAdjustmentTest, true )

CM_TEST_DESCRIPTOR_FULL(
    fractionalDates,
        FractionalDatesObjects, 10, FractionalDatesTest, true )

CM_TEST_DESCRIPTOR_FULL(
    riskyFixed,
        RiskyBondFixedObjects, 10, RiskyBondFixedTest, true )

CM_TEST_DESCRIPTOR_FULL(
    riskyFloat,
        RiskyBondFloatObjects, 10, RiskyBondFloatTest, true )

CM_TEST_DESCRIPTOR_FULL(
    riskyZero,
        RiskyZeroObjects, 10, RiskyZeroTest, true )

CM_TEST_DESCRIPTOR_FULL(
    cashflowInstrument,
        CashflowInstrumentObjects, 10, CashflowInstrumentTest, true )

CM_TEST_DESCRIPTOR_FULL(
    parAssetSwap,
        ParAssetSwapObjects, 10, ParAssetSwapTest, true )

CM_TEST_DESCRIPTOR_FULL(
    assetSwapOption,
        AssetSwapOptionObjects, 10, AssetSwapOptionTest, true )

CM_TEST_DESCRIPTOR_FULL(
    zeroCurve,
        ZeroCurveObjects, 10, ZeroCurveTest, true )

CM_TEST_DESCRIPTOR_FULL(
    volCurve,
        VolCurveObjects, 10, VolCurveTest, true )

CM_TEST_DESCRIPTOR_FULL(
    cashflow,
        CashflowObjects, 10, CashflowTest, true )

CM_TEST_DESCRIPTOR_FULL(
    compositeSpread,
        CompositeSpreadObjects, 10, CompositeSpreadTest, true )

CM_TEST_DESCRIPTOR_FULL(
    specialCase,
        SpecialCaseObjects, 10, SpecialCaseTest, true )

CM_TEST_DESCRIPTOR_FULL(
    kapital,
        KapitalObjects, 10, KapitalTest, true )

CM_TEST_DESCRIPTOR_FULL(
    bloomberg,
        BloombergObjects, 10, BloombergTest, true )

CM_TEST_DESCRIPTOR_FULL(
    numerical,
        NumericalAlgorithmsObjects, 10, NumericalAlgorithmsTest, true )

CM_TEST_DESCRIPTOR_FULL(
    timing,
        TimingObjects, 10, TimingTest, false )

//CM_TEST_DESCRIPTOR_FULL(
//    waterfall,
//        WaterfallTestObjects, 10, WaterfallTest, true )

} // CM

//
//   u s a g e
//

static void usage(
    int exitValue = 0 )
{
    Array<String> modules( TheTestManager().getModules() );

    // Print usage message.

    cerr << theProgramName << " "
         << "[-h] [moduleName [iteration number]]\n\n"
         << "\twhere -h produces this message and\n"
         << "\tmoduleName is one of the following:" << endl << endl;

    for( size_t i = 0; i < modules.size(); i++ )
        cerr << "\t\t" << modules[i] << endl;

    exit( exitValue );
}

//
//   c h e c k M o d u l e N a m e
//

static void checkModuleName(
    const char* moduleName )
{
    Array<String> modules( TheTestManager().getModules() );

    for( size_t i = 0; i < modules.size(); i++ ) {
        if( moduleName == modules[i] )
            return;
    }

    cerr << theProgramName
         << ": Invalid module name '" << moduleName << "'"
         << endl << endl;

    usage( -1 );  // this exits
}

//
//   m a i n
//

int main(
    int argc,
    char *argv[] )
{
    SetDebugTrace( false );

    cout << theBanner << endl;

    if( argc > 3 )
        usage( -1 );    // this exits

    if( argc == 2 && CompareNoCase( argv[1], "-h" ) == 0 )
        usage();        // this exits

    // Turn on Alib error logging.

    GtoErrMsgOn();

    try {
        if( argc == 1 ) {
            // No arguments -> run all regression test.

            TheTestManager().testAll();
        }
        else {
            // Regression test for a particular module.

            String module = argv[1];
            checkModuleName( argv[1] );

            if( argc == 2 ) 
                TheTestManager().test( module );
            else 
                TheTestManager().test( module, atoi( argv[2] ) );
        }
    }
    catch( Exception& e )
    {
        cerr << e.what() << endl;
        exit( -1 );
    }
    catch( ... )
    {
        cerr << "Unknown exception encountered." << endl;
        exit( -1 );
    }

    return 0;
}
