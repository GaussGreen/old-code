#include "alibfunclog.h"

#include <stdarg.h>
#include <stdio.h>

#include <alib/addinlib.h>
#include <alib/cerror.h>
#include <alib/cfileio.h>
#include <alib/modlcall.h>
#include <alib/strutil.h>

static TVar* makeVarFromInput
(TParamDef *paramDef,
 TAddinLib *addinLib,
 void      *data);

static int writeFunctionCall
(FILE           *fp,
 char           *libName,
 TObjectManager *om,          /* (I) Object Manager */
 TFuncDef       *funcDef,     /* (I) Function Definition*/
 TArgs          *inputs       /* (I) Input Perameters*/
);

int CxAlibFuncLog
(char    *fileName,
 TBoolean append,
 char    *funcName,
 ...)
{
    int status = FAILURE;
    static char routine[] = "CxAlibFuncLog";

    TFuncDef   *funcDef;
    FILE       *fp = NULL;
    size_t      i;
    va_list     ap;
    void      **rawArgs = NULL;
    TArgs      *args = NULL;
    TAddinLib  *addinLib = NULL;

    if (fileName != NULL)
    {
        fp = fopen (fileName, append ? "a" : "w");
        if (!fp)
        {
            GtoErrMsg ("%s: Cannot open file %s\n", routine, fileName);
            goto done; /* failure */
        }
    }

    addinLib = GtoGetAddinLib (GtoGetAddinLibName());
    if (addinLib == NULL)
    {
        GtoErrMsg ("%s: Cannot get addin library details\n", routine);
        goto done; /* failure */
    }

    funcDef = GtoGetFuncDefFromName (funcName,
                                     addinLib->numFuncs,
                                     addinLib->funcNames,
                                     addinLib->funcDefs);
    if (funcDef == NULL)
    {
        GtoErrMsg ("%s: Cannot find function %s\n", routine, funcName);
        goto done; /* failure */
    }

    args = GtoNewArgs (funcDef->numInputs + funcDef->numOutputs);
    
    if (funcDef->numInputs > 0)
    {
        rawArgs = NEW_ARRAY(void*, funcDef->numInputs);
        va_start (ap, funcName);
        for (i = 0; i < funcDef->numInputs; ++i)
        {
            rawArgs[i] = va_arg(ap, void*);
        }
        va_end (ap);

        for (i = 0; i < funcDef->numInputs; ++i)
        {
            TParamDef *paramDef = funcDef->params + i;
            args->vars[i] = makeVarFromInput (paramDef, addinLib, rawArgs[i]);
        }
    }

    writeFunctionCall (fp,
                       addinLib->name,
                       addinLib->om,
                       funcDef,
                       args);
    
    status = SUCCESS;

done:

    GtoFreeArgs (args);
    FREE (rawArgs);
    
    if (fp)
        fclose (fp);

    return status;
}

/* call this function if you are not calling from a regular alib add-in */
void CxAlibFuncLogInit (char *libname, TRegisterObjectTypesFunc objRegFunc)
{
    static TBoolean called = FALSE;

    if (!called)
    {
        extern TFuncName gtoFuncNameList[];
        extern TFuncDef *gtoFuncDefList[];
        extern size_t gtoFuncDefCount;

        TObjectManager *om;

        called = TRUE;
        
        om = GtoNewObjectManager();
        if (objRegFunc != NULL) objRegFunc(om);
        GtoRegisterAddinLib (libname, om, gtoFuncDefCount,
                             gtoFuncDefList, gtoFuncNameList);
    }
}

static TVar* makeVarFromInput
(TParamDef *paramDef,
 TAddinLib *addinLib,
 void      *data)
{
    static char routine[] = "makeVarFromInput";

    void *varData = NULL; /* will be owned by TVar */
    TVar *var = NULL;
    TObjectType *objectType = NULL;

    if (paramDef->isScalar)
    {
        switch (paramDef->varType)
        {
        case GTO_VT_LONG:
        {
            long *p = NEW_ARRAY(long,2);
            p[0] = 1;
            p[1] = *((long*)data);
            varData = (void*)p;
            break;
        }
        case GTO_VT_DATE:
        {
            TDate *p = NEW_ARRAY(TDate,2);
            p[0] = 1;
            p[1] = *((TDate*)data);
            varData = (void*)p;
            break;
        }
        case GTO_VT_DOUBLE:
        {
            double *p = NEW_ARRAY(double,2);
            p[0] = 1;
            p[1] = *((double*)data);
            varData = (void*)p;
            break;
        }
        case GTO_VT_CHAR_BLOCK:
        {
            char *p = NEW_ARRAY(char,129);
            p[0] = 1;
            strcpy (p+1, (char*)data);
            varData = (void*)p;
            break;
        }
        case GTO_VT_STRING:
        {
            char **p = NEW_ARRAY(char*,2);
            p[0] = (char*)1;
            p[1] = GtoStringDuplicate((char*)data);
            varData = (void*)p;
            break;
        }
        case GTO_VT_OBJECT:
        {
            void **p;
            objectType = GtoGetObjectType (addinLib->om,
                                           paramDef->classType);
            if (objectType == NULL)
                goto done; /* failure */
            p = NEW_ARRAY(void*,2);
            p[0] = (void*)1;
            p[1] = data;
            varData = (void*)p;
            break;
        }
        default:
            goto done; /* failure */
        }
    }
    else
    {
        switch (paramDef->varType)
        {
        case GTO_VT_LONG:
        {
            long n = ((long*)data)[0];
            long *p = NEW_ARRAY(long, n+1);
            COPY_ARRAY (p, data, long, n+1);
            varData = (void*)p;
            break;
        }
        case GTO_VT_DATE:
        {
            long n = (long)((TDate*)data)[0];
            TDate *p = NEW_ARRAY(TDate, n+1);
            COPY_ARRAY (p, data, TDate, n+1);
            varData = (void*)p;
            break;
        }
        case GTO_VT_DOUBLE:
        {
            long n = (long)((double*)data)[0];
            double *p = NEW_ARRAY(double, n+1);
            COPY_ARRAY (p, data, double, n+1);
            varData = (void*)p;
            break;
        }
        case GTO_VT_CHAR_BLOCK:
        {
            long n = (long)((char*)data)[0];
            char *p = NEW_ARRAY(char, 128*n+1);
            COPY_ARRAY (p, data, char, 128*n+1);
            varData = (void*)p;
            break;
        }
        case GTO_VT_STRING:
        {
            long n = (long)((char**)data)[0];
            long i;
            char **p = NEW_ARRAY(char*, n+1);
            p[0] = (char*)n;
            for (i = 1; i <= n; ++i)
                p[i] = GtoStringDuplicate (((char**)data)[i]);
            varData = (void*)p;
            break;
        }
        case GTO_VT_OBJECT:
        {
            long n;
            void **p;

            objectType = GtoGetObjectType (addinLib->om,
                                           paramDef->classType);
            if (objectType == NULL)
                goto done; /* failure */
            
            n = (long)((void**)data)[0];
            p = NEW_ARRAY(void*, n+1);
            COPY_ARRAY (p, data, void*, n+1);
            varData = (void*)p;
            break;
        }
        default:
            goto done; /* failure */
        }
    }

    var = GtoNewVar (paramDef->varType,
                     objectType,
                     paramDef->isScalar,
                     varData,
                     NULL);

done:

    if (var == NULL)
        GtoErrMsgFailure (routine);

    return var;
}
    


static int writeFunctionCall
(
    FILE           *fp,
    char           *libName,
    TObjectManager *om,          /* (I) Object Manager */
    TFuncDef       *funcDef,     /* (I) Function Definition*/
    TArgs          *inputs       /* (I) Input Perameters*/
 )
{
    char           routine[] = "writeFunctionCall";
    int            status    = FAILURE; /* Until proven otherwise*/
    TModelCall    *funcCall  = NULL;
    char          *str       = NULL;

    funcCall = GtoModelCallMakeFromArgs(om,libName,funcDef,inputs);
    if (funcCall == NULL)
        goto done;

    str = GtoSerializeObject (funcCall, GtoModelCallObjectType(), NULL, FALSE);
    if (str == NULL)
        goto done;

    if (fp == NULL)
    {
        /* don't use GtoErrMsg because the output is going to be too long */
        GtoErrLogWrite (str);
        GtoErrLogWrite ("\n");
    }
    else
    {
        fprintf (fp, "%s\n", str);
    }

    status = SUCCESS;

done:
    
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Unable to print Function Call.\n", routine);
    }
    FREE (str);
    GtoModelCallDelete (funcCall);
    return status;
}

