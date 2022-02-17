/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : armval.cpp                                                   */
/*                                                                            */
/* DESCRIPTION : global utilities                                             */
/*                                                                            */
/* DATE        : Thu Feb 1998                                                 */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include "armval.h"




int ARM_Val::operator == (ARM_Val& srcVal)
{
    if ( valType != srcVal.valType )
    {
       return(0);
    }

    switch(srcVal.valType)
    {
        case ARM_INT :
        {
            return( val.intVal == srcVal.val.intVal );
        };
        break; 
 
        case ARM_DOUBLE :
        {
            return( val.doubleVal == srcVal.val.doubleVal );
        };
        break;
           
        case ARM_STRING:
        {
            return(strcmp(val.stringVal, srcVal.val.stringVal) == 0 );
        };
        break;
 
        default :
                return(0);
                break;
     }
}



int ARM_Val::operator < (ARM_Val& srcVal)
{
    if ( valType != srcVal.valType )
    {
       return(0);
    }

    switch(srcVal.valType)
    {
        case ARM_INT :
        {
            return( val.intVal < srcVal.val.intVal ); 
        };
        break; 
 
        case ARM_DOUBLE :
        {
            return( val.doubleVal < srcVal.val.doubleVal ); 
        };
        break;
          
        case ARM_STRING:
        {
            return(strcmp(val.stringVal, srcVal.val.stringVal) < 0 );
        };
        break;
 
        default :
            return(0);
        break;
     };
}


void ARM_Val::View(char* id, FILE* ficOut)
{
    FILE* fOut;

    char fOutName[200];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n ARM_Val");


    switch(valType)
    {
        case ARM_INT :
        {
			fprintf(fOut, "\n Value : %d\n\n", val.intVal); 
        };
        break;

        case ARM_DOUBLE :
        {
			fprintf(fOut, "\n Value : %10.5lf\n\n", val.doubleVal);
        };
        break;
    
        case ARM_STRING:
        {
			fprintf(fOut, "\n Value : %s\n\n", val.stringVal);
        };
        break;

        default : 
        break;
    };

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
