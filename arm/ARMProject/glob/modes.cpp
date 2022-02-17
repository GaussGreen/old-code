/*
 * $Log: modes.cpp,v $
 * Revision 1.2  2001/04/27 09:25:12  smysona
 * IsCallOrPut
 *
 * Revision 1.1  2001/01/30 09:52:36  smysona
 * Initial revision
 *
 */


#include "modes.h"

int IsCallOrPut(int Type)
{
    if ((Type & CALIB_PORTFOLIO_MODE::Call) == CALIB_PORTFOLIO_MODE::Call)
    {
        return 1;
    }
    else
    {
        if ((Type & CALIB_PORTFOLIO_MODE::Put) == CALIB_PORTFOLIO_MODE::Put)
        {
            return -1;
        }
    }

    return 0;
}

