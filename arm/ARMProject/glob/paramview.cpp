/*
 * Copyright CDC IXIS CM : Paris July 2003
 *
 * $Log: paramview.cpp,v $
 * Revision 1.6  2004/05/28 11:50:42  jmprie
 * modif du GetMappingName() pour eviter la boucle infinie qd il n'y a pas
 * de match pour un type donne !!
 *
 * Revision 1.5  2003/08/01 14:08:43  mab
 * Just Formatting
 *
 * Revision 1.4  2003/08/01 10:16:34  ekoehler
 * ARM formatting
 *
 * Revision 1.2  2003/07/31 19:40:21  ekoehler
 * ARM style formatting.
 *
 * Revision 1.1  2003/07/30 12:12:41  ekoehler
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*
    paramview.cpp
 
    This file implements the ARM_ParamView class

*----------------------------------------------------------------------------*/


#include "paramview.h"
#include <string.h>
#include <stdio.h>


/*----------------------------------------------------------------------------*/

/*
 * general function to find the name that corresponds to a constant
 */

char* ARM_ParamView::GetMappingName( char* GenericTag, long MethodFlag )
{
    int NbOfRows = sizeof(itsParamMappingTable)/sizeof(itsParamMappingTable[0]);
    char* msg;

    int i = 0;
    bool Found=false;
    while(i < NbOfRows)
    {
        if(!strcmp(itsParamMappingTable[i].itsGenericTag, GenericTag))
        {
            int j = 0;

            while((Found = itsParamMappingTable[i].itsTagArray[j].itsMethodFlag
                           != ARM_PARAMVIEW_ENDOFLINE))
            {
                if(itsParamMappingTable[i].itsTagArray[j].itsMethodFlag == MethodFlag )
                {
                    msg = itsParamMappingTable[i].itsTagArray[j].itsMethodName;
                    break;
                }
                ++j;
            }
            break;
        }
        ++i;
    }

    if(!Found)
        return NULL;

    return(msg);
}



/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

