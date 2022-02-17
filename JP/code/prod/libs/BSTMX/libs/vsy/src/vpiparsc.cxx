/***************************************************************
 * Module:        BasisTree
 * Submodule:        
 * File:        
 * Function:        
 * Author:        C. Daher and David Liu
 ***************************************************************/
#include "vpbase.h"        /* Definitions and error handling */

#include "vpbundle.h"
#include "vpcashfl.h"
#include "vpoption.h"
#include "vpfleg.h"
#include "vpkio.h"
#include "vpkionew.h"
#include "vpkio2idx.h"
#include "vpfleg2Idx.h"
#include "vpfleg3Idx.h"
#include "vpfleg.h"
#include "vpfleg2Idx.h"
#include "vpfleg3Idx.h"
#include "vpprotleg.h"
#include "vpdefprotect.h"
#include "vpflegnew.h"

#include "kcplxrate.h"      // complex rate index
#include "kindicator.h"     // indicator 

#include "kutilios.h"
#include "kmemutil.h"

#include "vpipars.h"                // prototype consistency

// Global Constant Table

KMap(String,double)                globVpConstTable;



static        KVector(SharedPointer<KVPAtom>)        *vpArray = NULL;
static        SharedPointer<KVPInstr>                vpRoot;
static        String                        *sVpDiscName = NULL;
static        KResetBank                *sVpResetBank = NULL;
static        KMap(String,double)        *sVpConstTable = NULL;        // Constant table



extern "C" {

#include "drlio.h"        /* DrlFPrintf */
#include "drltime.h"        /* TDayCount */
#include "drlmem.h"        /* Mem alloc */

#define        _DINSPAR_SRC
#include "vpiparsh.h"        /* TvpiparsSymbol */
};

/*----------------------------------------------------------------------
 * Type for yacc interface lists
 */

struct        KVpiListItem {
        KVpiListItem(char type)
        { mType = type; }

        char                mType;
        KVector(TDate)        mDtval;
        KVector(double)        mDval;
        KVector(String) mSval;

        size_t        size() const
        {
            switch (mType) {
            case 'D': return mDtval.size(); break;
            case 'd': return mDval.size(); break;
            case 's': return mSval.size(); break;
            default:  throw KFailure("Bad type %c\n", (char) mType);
            }
            return(-1);
        }
};

class        KVpiList : public KVector(KVpiListItem) {
public:
        void        InsertColumn(char type)
        {
                insert(end(), KVpiListItem(type));
        }

        void        CheckType(int idxF, char type, char *routine)
        {
                if (idxF >= size())
                        throw KFailure("%s: idxF %d >= size %d.\n",
                                routine, idxF, size()); 
                if ((*this)[idxF].mType != type)
                        throw KFailure("%s: expected type %c, got %c.\n",
                                routine, type, (*this)[idxF].mType);
        }


        void        AddDouble(int idxF, double value)
        {
                CheckType(idxF, 'd', "KVpiList::AddDouble");
                (*this)[idxF].mDval.insert((*this)[idxF].mDval.end(), value);
        }

        void        AddTDate(int idxF, TDate value)
        {
                CheckType(idxF, 'D', "KVpiList::AddTDate");
                (*this)[idxF].mDtval.insert((*this)[idxF].mDtval.end(), value);
        }

        void        AddString(int idxF, char* value)
        {
                CheckType(idxF, 's', "KVpiList::AddString");
                (*this)[idxF].mSval.insert((*this)[idxF].mSval.end(), 
                            (String) value);
        }


        KVector(double)        VectDouble(int idxF)
        {
                CheckType(idxF, 'd', "KVpiList::VectDouble");
                return (*this)[idxF].mDval;
        }

        KVector(TDate)        VectTDate(int idxF)
        {
                CheckType(idxF, 'D', "KVpiList::VectTDate");
                return (*this)[idxF].mDtval;
        }

        KVector(String) VectString(int idxF)
        {
                CheckType(idxF, 's', "KVpiList::VectString");
                return (*this)[idxF].mSval;
        }

friend        ostream& operator<<(ostream& os, const KVpiList& object)
        {
                int        idxF, numF, numR = 0;

                numF = object.size();
                for (idxF=0; idxF<numF; idxF++) {
                        numR = MAX(numR, object[idxF].size());
                }
                os << format("KVpiList: %3dx%d\n", numR, numF);
                for (idxF=0; idxF<numF; idxF++) {
                        os << "\t`" << object[idxF].mType << "'";
                }
                os << endl;
                return (os);
        }
};

                                // Garbage collection (needed for yacc)
static        KGCollect<KVpiList>        sVpiListGc;









/*----------------------------------------------------------------------
 * Evaluates the input 
 */

SharedPointer<KVPInstr>
VPInsParEvalFile(
        const char* fnam,                // (I) file name
        KResetBank &vpResetBank,        // (O) Resets table
        String &vpDiscName)                // (O) discount curve name
{
        SharedPointer<KVPInstr>        vpRoot;
        String        txtFile;

    try {

        txtFile = DppFileToString(fnam);


        VPInsParEval(
                txtFile.c_str(),
                NULL,
                &vpRoot,
                globVpConstTable,
                vpResetBank,
                vpDiscName);

        return(vpRoot);
    }
    catch (KFailure) {
        throw KFailure("VPInsParEvalFile: failed.\n");
    }
}


/*----------------------------------------------------------------------
 * Evaluates the input 
 */

void
VPInsParEval(
        const char *inp,                           // (I) input string
        KVector(SharedPointer<KVPAtom>) *vpArrayO, // (O) array of instruments 
        SharedPointer<KVPInstr> *vpRootO,                 // (O) root instruments
        KMap(String,double) &vpConstTable,                // (O) Constants table
        KResetBank &vpResetBank,                        // (O) Resets table
        String &vpDiscName)                                // (O) discount curve name

{
static        char        routine[] = "VPInsParEval";
        int        errCode;

        vpiparsSymbTableG = NULL;

    try {

        // Set static variable to communicate with wrappers
        //
        vpArray = NULL;
        vpArray = new KVector(SharedPointer<KVPAtom>);

        sVpConstTable = &vpConstTable;
        sVpResetBank  = &vpResetBank;
        sVpDiscName   = &vpDiscName;


        /* Set static variable */

        /* Allocate global symbol table */
        vpiparsSymbTableG = _vpiparsSymbolTableNew(NUMSYMMAX);
        ASSERT_OR_THROW(vpiparsSymbTableG != NULL)


        /* initialize the lexer */
        IF_FAILED_THROW( vpiparsLexInit(inp));


        /* parse */
        errCode = vpiparsparse();
        if (errCode != 0) throw KFailure();



        /* cleanup */
        _vpiparsSymbolTableClear(vpiparsSymbTableG);
        vpiparsLexClean();

        // 
                dppLog << "------------------------------------------" << endl;
                dppLog << "CONSTANTS TABLE" << endl;
                for(KMap(String,double)::iterator p=sVpConstTable->begin();
                    p != sVpConstTable->end();
                    ++p) {
                        dppLog << format("%10s  %18.10f",
                                (const char*) ((*p).first).c_str(),
                                (*sVpConstTable)[(*p).first]) << endl;
                }

                if (!vpResetBank.IsEmpty()) {
                dppLog << "------------------------------------------" << endl;
                dppLog << vpResetBank;
                }

                int        idx, n = vpArray->size();
                dppLog << "------------------------------------------" << endl;
                dppLog << routine << ": Output Array" << endl;
                for (idx=0; idx<vpArray->size(); idx++) {
                        dppLog << format("----- VP %2d/%2d", idx+1, n);
                        dppLog << format(" %15s", (*vpArray)[idx]->TypeName());
                        dppLog << format(" %15s", (*vpArray)[idx]->GetName());
                        dppLog << endl;
                        dppLog << *((*vpArray)[idx]);
                }
                dppLog << routine << ": Root : ";
                if (vpRoot) dppLog << vpRoot->GetName();
                else dppLog << "NOT SET" << endl;
                dppLog << endl;

                dppLog << "------------------------------------------" << endl;
                dppLog << "Root Dependency Graph: " << endl;
                vpRoot->PutTree(dppLog);
                dppLog << "------------------------------------------" << endl;

        if (vpArrayO) {
                vpArrayO = vpArray;
        } else {
                delete vpArray;
        }
        if (vpRootO) {
                *vpRootO  = vpRoot;
                if (*vpRootO == NULL) {
                        throw KFailure("%s: root instrument undefined.\n",
                        routine);
                }


        }

        vpiListCleanup();

    }
    catch (KFailure) {
        /* cleanup */
        _vpiparsSymbolTableClear(vpiparsSymbTableG);
        vpiparsLexClean();
        delete vpArray;
        
        vpiListCleanup();

        throw KFailure("%s: error parsing.\n", routine);

    }
    catch (...) {
        /* cleanup */
        _vpiparsSymbolTableClear(vpiparsSymbTableG);
        vpiparsLexClean();
        delete vpArray;
        vpiListCleanup();

        throw KFailure("%s: non KFailure exception.\n", routine);

    }
}



//----------------------------------------------------------------------
// Adds a constant to the constant table.
// We need this routine because we can't call the template members
// under C linkage.
//

static        void
addConstantToTable(const char *name, double cValue)
{
        if (sVpConstTable) {
                String        cName = String(name);
                (*sVpConstTable)[cName] = cValue;
        }
}



/*======================================================================
 *
 *        WARNING:      C-LINKAGE SECTION
 *
 *        All the following is defined with C linkage
 *        in order to be called by lex/yacc
 *
 *======================================================================*/

extern "C" {

extern        char*        TSymTypeName(TSymType type);

/*----------------------------------------------------------------------
 * Creates a new symbol table.
 */

TvpiparsSymbolTable*
_vpiparsSymbolTableNew(int numSymMax)
{
        TvpiparsSymbolTable        *symbtab = NULL;
        int                        i;

        symbtab = NEW(TvpiparsSymbolTable);
        if (symbtab == NULL) goto failed;

        symbtab->table = NEW_ARRAY(TvpiparsSymbol, numSymMax);
        if (symbtab->table == NULL) goto failed;

        symbtab->numSymMax = numSymMax;

        /* clear */
        for (i=0; i<=symbtab->numSymMax-1; i++) {
                symbtab->table[i].name[0] = '\0';
                symbtab->table[i].setFlag = FALSE;
                symbtab->table[i].type = SYMNIL;
                symbtab->table[i].ptr = NULL;
        }
        symbtab->numSym = 0;
        return(symbtab);

failed:
        DppErrMsg("_vpiparsSymbolTableNew: failed (numSymMax = %d).\n",
                numSymMax);
        return(NULL);
}

/*----------------------------------------------------------------------
 * Clear symbol table.
 */

int
_vpiparsSymbolTableClear(TvpiparsSymbolTable *symbtab)
{
        int        i;

        for (i=0; i<=symbtab->numSym-1; i++) {
                symbtab->table[i].name[0] = '\0';
                symbtab->table[i].setFlag = FALSE;
        }

        symbtab->numSym = 0;

        return(SUCCESS);
}



/*----------------------------------------------------------------------
 * Checks if symbol "s" already in table and returns it if
 * exists. Creates new symbol otherwise.
 */

int
_vpiparsSymbolTableLookup(
        TvpiparsSymbolTable *symbtab,        /* (I) symbol table */
        char *s,                        /* (I) symbol name */
        TvpiparsSymbol **sym)                /* (O) */
{
static        char        routine[] = "_vpiparsSymbolTableLookup";
        int                i;

        if (strlen(s) >= INSPARMAXS) {
                DppErrMsg("%s: symbol `%s' too long.\n", routine, s);
                return(FAILURE);
        }

        for (i=0; i<=symbtab->numSym-1; i++) {
            /* is it already here? */
            if (symbtab->table[i].name &&
                !strcmp(symbtab->table[i].name, s)) {
                        *sym = &symbtab->table[i];
                        return(SUCCESS);
            }

        }

        /* check enough mem */
        if (symbtab->numSym > symbtab->numSymMax-1) {
                DppErrMsg("%s: too many symbols (max %d).\n",
                        routine, symbtab->numSymMax);
                return(FAILURE);
        }

        /* symb not there, create new one */
        i = symbtab->numSym;
        strcpy(symbtab->table[i].name, s);
        symbtab->table[i].setFlag = FALSE;
        symbtab->table[i].type = SYMNIL;
        symbtab->numSym++;

        *sym = &symbtab->table[i];
        return(SUCCESS);
}



/*----------------------------------------------------------------------
 * Prints the symbol table.
 */

int
_vpiparsSymbolTablePrint(TvpiparsSymbolTable *symbtab)
{
static        char        routine[] = "_vpiparsSymbolTablePrint";
        int                i;
        TvpiparsSymbol        *symb;

        dppLog << "vpiparsSymbolTable: ACTIVE SYMBOL TABLE:\n";

        dppLog << "NO/TOT SYMBOL NAME           "
                "SYMBOL TYPE         STAT  VALUE\n" << endl;

        for (i=0; i<=symbtab->numSym-1; i++) {
            symb = &symbtab->table[i];


            dppLog << format("%2d/%2d  %-20s  %-15s  %6s ",
                i+1, symbtab->numSym,
                symbtab->table[i].name,
                TSymTypeName(symbtab->table[i].type),
                (symbtab->table[i].setFlag ? "SET" : " NOTSET"));

            switch (symb->type) {
            case SYMSTRING:
                dppLog << format("\"%s\"", ((char*)symb->ptr));
                break;
            case SYMDOUBLE:
                dppLog << format("%12.8f", symb->dval);
                break;
            case SYMTDATE:
                dppLog << format("%10s", DrlTDatePrint(NULL, symb->dtval));
                break;
            case SYMVPBASE:
                dppLog << format("%p", symb->ptr);
                break;
            default:
                dppLog << "N/A";
                break;
            }

            dppLog << endl;
        }
        return(SUCCESS);
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*x                                                                    x*/
/*x                                                                    x*/
/*x                                                                    x*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/*----------------------------------------------------------------------
 * 
 */

int
_vpiparsSymbolSetValue(TvpiparsSymbol *symb, TSymType type, void *value)
{

        if ((symb->type != SYMNIL) ||
            (symb->setFlag != FALSE)) {
            DppErrMsg("SetValue: variable `%s' of type %s is already set.\n",
                symb->name, TSymTypeName(symb->type));
            return(FAILURE);
        }

        switch (type) {
        case SYMDOUBLE:
                symb->dval = *((double*) value);

                //
                // Add to table of constants
                //
                addConstantToTable(symb->name, *((double*) value));


                break;

        case SYMTDATE:
                symb->dtval = *((TDate*) value);
                break;

        case SYMSTRING:
                if (strlen(((char*)value)) >= INSPARMAXS)
                {
                        DppErrMsg("%s: string `%s' too long.\n", 
                                  (char*)value);
                        return(FAILURE);
                }

                strcpy(symb->strval, ((char*)value));
                break;
        default:
                DppErrMsg("Cannot set value of type %s in variable `%s'.\n",
                        TSymTypeName(type), symb->name);
                return(FAILURE);
        }

        symb->type = type;
        symb->setFlag = TRUE;

        return(SUCCESS);
}



/*----------------------------------------------------------------------
 * 
 */

int
_vpiparsSymbolGetValue(TvpiparsSymbol *symb, TSymType type, void *value)
{
static        char routine[] = "_vpiparsSymbolGetValue";
        if (symb->setFlag == FALSE) {
            DppErrMsg("%s: undefined variable `%s'.\n", routine, symb->name);
            return(FAILURE);
        }

        if (type != symb->type) {
            DppErrMsg("%s: variable `%s' is of type %s, "
                "expected type %s.\n", routine,
                symb->name, TSymTypeName(symb->type), TSymTypeName(type));
            return(FAILURE);
        }

        switch (type) {
        case SYMDOUBLE:
                *((double*) value) = symb->dval;
                break;
        case SYMTDATE:
                *((TDate*) value) = symb->dtval;
                break;
        case SYMSTRING:
                strncpy(((char*) value), symb->strval, INSPARMAXS);
                break;
        default:
                DppErrMsg("%s: cannot get value of variable `%s' of type %s.\n",
                        symb->name, TSymTypeName(symb->type));
                return(FAILURE);
        }
        symb->setFlag = TRUE;

        return(SUCCESS);
}

/*----------------------------------------------------------------------
 * Sets the root instrument.
 */

int
_vpiparsSymbolSetRoot(
        TvpiparsSymbol *insSymb)                /* (I) instr symb */
{
static char routine[] = "_vpiparsSymbolSetRoot";
    try {


        /* Check all symb represent instruments */
        if (insSymb->type != SYMVPBASE) {
                DppErrMsg("%s: variable `%s' is not an instrument.\n",
                        routine, insSymb->name);
                throw KFailure();
        }

        vpRoot = Raw2SharedPointer(((KVPInstr*)insSymb->ptr));

        return(SUCCESS);
    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Sets a fixed cash flows list.
 */

int
_vpiparsSetDiscName(
        char *discZcName)                /* (I) curve name */
{
static char routine[] = "_vpiparsSymbolSetDiscName";


    try {
        *sVpDiscName = String(discZcName);
        return(SUCCESS);
    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*x                                                                    x*/
/*x     List utility to handle "flows" style inputs                    x*/
/*x                                                                    x*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/



int
vpiListInit(
        void **oListNew,                /* (O) */
        char *types)                        /* (I) Ex: "Dds" for date,double,string */
{
        KVpiList        *listNew = sVpiListGc.New();
        char        *p;

        for (p=types; *p!='\0'; p++)
                (*listNew).InsertColumn(*p);
        *oListNew = (void*) listNew;
        //dppLog << __LINE__ << ":" << *listNew;

        return(SUCCESS);
}

int
vpiListAdd(
        void **oListNew,                /* (O) */
        void **oListOld,                /* (I) old list */
        char *types,                        /* (I) Ex: "Dds" for date,,double,string */
        ...)                                /* (I) all arguments */
{
        va_list        ap;
        char        *p;
        int        idxF;

        if (*oListOld == NULL) {
                DppErrMsg("Received null pointer `%s'.\n", types);
                return (FAILURE);
        }

        KVpiList        &listOld = *((KVpiList*) *oListOld);

        va_start(ap, types);
        for (p=types, idxF=0; *p!='\0'; p++, idxF++) {
            switch (*p) {
            case 'D': 
                listOld.AddTDate(idxF, (TDate) va_arg(ap, TDate));
                break;
            case 'd': 
                listOld.AddDouble(idxF, (double) va_arg(ap, double));
                break;
            case 's': 
                listOld.AddString(idxF, (char*) va_arg(ap, char*));
                break;
            default:
                DppErrMsg("Unknown type %c.\n", *p);
                return(FAILURE);
            }
        }
        va_end(ap);

        *oListOld = (void*) NULL;
        *oListNew = (void*) &listOld;
        //dppLog << __LINE__ << ":" << listOld;

        return(SUCCESS);
}


int
vpiListCleanup()
{
    try {
        sVpiListGc.Free();
        return (SUCCESS);
    }
    catch(KFailure) {
        return (FAILURE);
    }
}


/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*x                                                                    x*/
/*x     KVP classes wrappers                                           x*/
/*x                                                                    x*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
int
_vpiparsSymbolSetKDateInterval(
        TvpiparsSymbol  *symb,                /* (O) */
        TDateInterval        interval,        /* (I) Must be days if isBusDays=TRUE */
        int                dateAdjType,    /* (I) Use GTO_DATE_ADJ_TYPE_... flag */
        char                *holidayFile,        /* (I) Set to "NO_WEEKENDS" for no adj*/
        long                badDayConv)        /* (I) Only applies if isBusDays=FALSE*/
{
static        char routine[] = "vpiparsSymbolSetKDateInterval";

        TDateAdjIntvl        *intvlAdj  = NULL;
        KDateInterval        *kIntvlAdj = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        /* Create object */
        intvlAdj = GtoDateIntervalAdjustedNew(
                                        &interval,
                                        dateAdjType,
                                        holidayFile,
                                        badDayConv);

        if (intvlAdj == NULL) return(FAILURE);

        KDateInterval *kIntvlAdj = new KDateInterval(*intvlAdj);
        if (kIntvlAdj == NULL) return(FAILURE);

        symb->type = SYMDATEINT;
        symb->ptr = (void*) kIntvlAdj;
        symb->setFlag = TRUE;

        GtoDateIntervalAdjustedFree(intvlAdj);

        return(SUCCESS);

    }
    catch (KFailure) {
        GtoDateIntervalAdjustedFree(intvlAdj);
        delete kIntvlAdj;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}





/*----------------------------------------------------------------------
 * Sets a variable to a rate.
 */


int
_vpiparsSymbolSetKRate_Fixed(
        TvpiparsSymbol *symb,                /* (O) */
        double fixedRate)

{
static        char routine[] = "vpiparsSymbolSetKRate";

        KRate                *floatRate= NULL;
        static char discZcName[] = "Default";

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        /* Create object */
        floatRate = new KRate (fixedRate);

        if (floatRate == NULL) return(FAILURE);

        symb->type = SYMRATE;
        symb->ptr = (void*) floatRate;
        symb->setFlag = TRUE;

        /* Store discount curve name 
         */
        if (strlen(discZcName) >= INSPARMAXS) {
            DppErrMsg("%s: curve name `%s' too long.\n", routine, discZcName);
            return(FAILURE);
        }
        strcpy(symb->strval, discZcName);
        dppLog << "Curve name: " << discZcName << endl;


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatRate;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a floating rate with 0 spot offset.
 */

int
_vpiparsSymbolSetKRate_Float1(
        TvpiparsSymbol *symb,                /* (O) */
        TDateInterval mat,                /* (I) */
        TDateInterval freq,                /* (I) */
        long dayCc,                        /* (I) */
        char *curveName)                /* (I) */

{
static        char routine[] = "vpiparsSymbolSetKRate1";

        KRate                *floatRate= NULL;
        KDateInterval        spotOffset(0e0);
        double                spread = 0e0;
        double                weight = 1e0;

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        /* Create object */
        floatRate = new KRate (
                curveName,
                KDateInterval(mat),
                KDateInterval(freq),
                KDayCc(dayCc),
                spotOffset,
                spread,
                weight);

        if (floatRate == NULL) return(FAILURE);

        symb->type = SYMRATE;
        symb->ptr = (void*) floatRate;
        symb->setFlag = TRUE;

        /* Store curve name in symb structure to propagate to the 
         * product
         */
        if (strlen(curveName) >= INSPARMAXS) {
            DppErrMsg("%s: curve name `%s' too long.\n", routine, curveName);
            return(FAILURE);
        }
        strcpy(symb->strval, curveName);
        dppLog << "Curve name: " << curveName << endl;

        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatRate;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Sets a floating rate with given spot offset.
 */


int
_vpiparsSymbolSetKRate_Float2(
        TvpiparsSymbol *symb,                /* (O) */
        TDateInterval mat,                /* (I) */
        TDateInterval freq,                /* (I) */
        long dayCc,                        /* (I) */
        TvpiparsSymbol *spotOffSym,        /* (I) spot offset interval */
        char *curveName)                /* (I) */

{
static        char routine[] = "vpiparsSymbolSetKRate2";

        KRate                *floatRate = NULL;
        KDateInterval        *spotOffset = NULL;
        double                spread = 0e0;
        double                weight = 1e0;

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        if (spotOffSym->type != SYMDATEINT) {
                throw KFailure("%s: variable `%s' is not a KDateInterval.\n",
                        routine, spotOffSym->name);
        }

        spotOffset = (KDateInterval*) spotOffSym->ptr;

        /* Create object */
        floatRate = new KRate (
                curveName,
                KDateInterval(mat),
                KDateInterval(freq),
                KDayCc(dayCc),
                *spotOffset,
                spread,
                weight);

        if (floatRate == NULL) return(FAILURE);

        symb->type = SYMRATE;
        symb->ptr = (void*) floatRate;
        symb->setFlag = TRUE;

        /* Store curve name in symb structure to propagate to the 
         * product
         */
        if (strlen(curveName) >= INSPARMAXS) {
            DppErrMsg("%s: curve name `%s' too long.\n", routine, curveName);
            return(FAILURE);
        }
        strcpy(symb->strval, curveName);
        dppLog << "Curve name: " << curveName << endl;

        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatRate;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Sets a variable to a complex float index
 */
int
_vpiparsSymbolSetKRate_Cplx(
    TvpiparsSymbol *symb,       /* (O)               */
    int            numIns,      /* (I) # of instr    */
    TvpiparsSymbol **insSym,    /* (I) rate indexes  */
    char           *formula)    /* (I) index formula */
{
static        char routine[] = "_vpiparsSymbolSetKRateCplx";
    KCplxRate       *cplxfloat = NULL;
    KVector(KRate*) instrRates;

    try
    {     
        dppLog << "Setting symbol " << symb->name << endl;

        /* creat KRate vector*/
        for (int i = 0; i < numIns; i++)
        {
            if (insSym[i]->type != SYMRATE) 
            {
                    throw KFailure("%s: instrument `%s' is not a KRate.\n",
                            routine, insSym[i]->name);
            }
            instrRates.push_back( (KRate*)insSym[i]->ptr );
        }

        
        /* creat complex rate object*/
        cplxfloat = new KCplxRate( 
                        symb->name,
                        instrRates,
                        formula);

/*        cout << routine << '\t'
             << "Name " << cplxfloat->GetName() 
             << "Curve Name " << cplxfloat->CurveName() << endl;
 */       
        if (cplxfloat == NULL) return (FAILURE);
        
        symb->type = SYMCPLXRATE;
        symb->ptr  = (void*) cplxfloat;
        symb->setFlag = TRUE;

        instrRates.clear();

        return (SUCCESS);
    }
    catch (KFailure)
    {
        delete cplxfloat;
        instrRates.clear();
        DppErrMsg("%s: failed.\n", routine);  
        return (FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Sets a indicator
 */
int 
_vpiparsSymbolSetIndicator(
    TvpiparsSymbol *symb,       /* (O)                  */
    TvpiparsSymbol *rateSym,    /* (I) rate             */
    double          lbarrier,   /* (I) lower barrier    */
    double          hbarrier,   /* (I) higher barrier   */
    char            ioWindow,   /* (I) Inside or outside*/
    char            smooth)     /* (I) Smoothing method */
{
static char routine[] = "_vpiparsSymbolSetIndicator";
    KIndicator     *koIndicator = NULL;
    KCplxRate      *cplxfloatRate;
    KRate          *floatRate;
    KKnockIO        kioWindow;
    KSmooth         ksmooth;

    try
    {     
        dppLog << "Setting symbol " << symb->name << endl;
        cout << symb->name << endl;

        /* creat KRate vector*/


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

            
        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        /* Create floatLeg */
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;
            koIndicator = new KIndicator( 
                        symb->name,
                        cplxfloatRate,
                        lbarrier,
                        hbarrier,
                        kioWindow,
                        ksmooth);
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;
            koIndicator = new KIndicator( 
                        symb->name,
                        floatRate,
                        lbarrier,
                        hbarrier,
                        kioWindow,
                        ksmooth);
        }

        /* creat a Indicator object */
/*        koIndicator = new KIndicator( 
                        symb->name,
                        instrRates,
                        formula,
                        lbarrier,
                        hbarrier,
                        kioWindow,
                        ksmooth);*/
        
        if (koIndicator == NULL) return (FAILURE);
        
        symb->type = SYMINDICATOR;
        symb->ptr  = (void*) koIndicator;
        symb->setFlag = TRUE;


        return (SUCCESS);
    }
    catch (KFailure)
    {
        delete koIndicator;
        DppErrMsg("%s: failed.\n", routine);  
        return (FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a bundle.
 */

int
_vpiparsSymbolSetKVPWBundle(
        TvpiparsSymbol *symb,                /* (O) */
        int numIns,                        /* (I) # of instr */
        TvpiparsSymbol **insSymb,        /* (I) array of instr symb */
        double *insWei)                        /* (I) array of instr weight */
{
static        char routine[] = "_vpiparsSymbolSetBundle";
        KVPWBundle        *wb = new KVPWBundle(symb->name);

    try {
        int                idx;

        dppLog << "Setting symbol " << symb->name << endl;

        /* Check all symb represent instruments */
        for (idx=0; idx<=numIns-1; idx++) {
            if (insSymb[idx]->type != SYMVPBASE) {
                DppErrMsg("%s: variable `%s' is not an instrument.\n",
                        routine, insSymb[idx]->name);
                throw KFailure();
            }
        }

        /* Add instruments in the bundle */
        for (idx=0; idx<=numIns-1; idx++) {
                wb->AddDepWeight(
                        Raw2SharedPointer((KVPInstr*)insSymb[idx]->ptr),
                        insWei[idx]);
        }

        symb->type = SYMVPBASE;
        symb->ptr = (void*) wb;
        symb->setFlag = TRUE;

        // Add instrument to root bundle 
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));

        return(SUCCESS);
    }
    catch (KFailure) {
        delete wb;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Sets a fixed cash flows list.
 */

int
_vpiparsSymbolSetKVPCashFlows(
        TvpiparsSymbol *symb,                /* (O) */
        OVpiList *ocfSchedule,                /* (I) date and amount schedule */
        char *discZcName)                /* (I) curve name */
{
static char routine[] = "_vpiparsSymbolSetKVPCashFlows";


        KVpiList &cfSchedule = *((KVpiList*) *ocfSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        dppLog << "Curve name: " << discZcName << endl;

        //
        // Construct object 
        //
        symb->ptr = new KVPCashFlows(
                symb->name,
                cfSchedule.VectTDate(0),        // cash flow dates
                cfSchedule.VectDouble(1),        // cash flow amounts
                discZcName);
        if (symb->ptr == NULL) throw KFailure();

        symb->setFlag = TRUE;
        symb->type = SYMVPBASE;

        //
        // Add instrument to root bundle
        //
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));

        return(SUCCESS);
    }
    catch (KFailure) {
        delete symb->ptr;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

//**********************************************************************
//
// Option and KO constructors
//
//**********************************************************************



/*----------------------------------------------------------------------
 * Sets a variable to a KVPOption (simple)
 */

int
_vpiparsSymbolSetKVPOptionSimple(
        TvpiparsSymbol *symb,                /* (O) */
        char* type,                        /* (I) call/put */
        TBoolean american,                /* (I) Amer/Euro */
        TDate startDate,                /* (I) This date is not included */
        TDate matDate,                        /* (I) */
        TDateInterval freq,                /* (I) */
        int nDays,                        /* (I) # of notifcation days */
        double strike,                        /* (I) strike */
        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetOptionSimple";


        KDateInterval        notDays = KDateInterval(nDays, FALSE);
        TBoolean        stubAtEnd = FALSE;
        KVPAtom                *under = NULL;
        KVPOption        *option = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create option */
        option = new KVPOption(
                symb->name,
                KVPOptionType(type),
                american,
                startDate,
                matDate,
                freq,
                stubAtEnd,
                notDays,
                strike,
                discZcName);
        if (option == NULL) return(FAILURE);

        /* Add underlying */
        option->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) option;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete option;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Sets a variable to a KVPOption (arbitrary)
 * Strike pay dates are assumed to be the same as
 * settlement dates.
 */

int
_vpiparsSymbolSetKVPOption1(
        TvpiparsSymbol *symb,                /* (O) */
        char* type,                        /* (I) call/put */
        TBoolean american,                /* (I) Amer/Euro */
        OVpiList *ooptSchedule,                /* (I) reset,acc, etc. scehdule  */
        int nDays,                        /* (I) # of notifcation days */
        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetOption1";


        KDateInterval        notDays = KDateInterval(nDays, FALSE);
        KVPAtom                *under = NULL;
        KVPOption        *option = NULL;
        KVpiList &optSchedule = *((KVpiList*) *ooptSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        // Get underlying 
        //
        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create option */
        option = new KVPOption(
                symb->name,
                KVPOptionType(type),
                american,
                optSchedule.VectTDate(0),        // notifDates
                optSchedule.VectTDate(1),        // settleDates
                optSchedule.VectDouble(2),        // strikes
                notDays,
                discZcName);
        if (option == NULL) return(FAILURE);

        /* Add underlying */
        option->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) option;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete option;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a variable to a KVPOption (arbitrary)
 * Include the strike pay dates, which can be different
 * from settlement dates.
 */

int
_vpiparsSymbolSetKVPOption2(
        TvpiparsSymbol *symb,             /* (O) */
        char* type,                       /* (I) call/put */
        TBoolean american,                /* (I) Amer/Euro */
        OVpiList *ooptSchedule,           /* (I) reset,acc, etc. scehdule  */
        int nDays,                        /* (I) # of notifcation days */
        TvpiparsSymbol *underSym,         /* (I) underlying bundle */
        char *discZcName)                 /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetOption2";


        KDateInterval        notDays = KDateInterval(nDays, FALSE);
        KVPAtom                *under = NULL;
        KVPOption        *option = NULL;
        KVpiList &optSchedule = *((KVpiList*) *ooptSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        // Get underlying 
        //
        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create option */
        option = new KVPOption(
                symb->name,
                KVPOptionType(type),
                american,
                optSchedule.VectTDate(0),        // notifDates
                optSchedule.VectTDate(1),        // settleDates
                optSchedule.VectTDate(2),        // strikeDates
                optSchedule.VectDouble(3),        // strikes
                notDays,
                discZcName);
        if (option == NULL) return(FAILURE);

        /* Add underlying */
        option->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) option;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete option;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




//**********************************************************************
//
// Knock In/Out constructors
//
//**********************************************************************


/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 * Effective observation dates are implicitly given in the index rate
 */

int
_vpiparsSymbolSetKVPKnockIO1(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        char ioWindow,                        /* (I) Knock In/Out window */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */

        OVpiList *oknockSchedule,        /* (I) reset,acc, etc. scehdule  */

        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPKnockIO1";

        KKnockIO        kioType, kioWindow;
        KSmooth                ksmooth;

        KRate                *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KVPAtom                *under = NULL;
        KVPKnockIO        *knockIO = NULL;
        KVpiList &knockSchedule = *((KVpiList*) *oknockSchedule);


    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);

        //
        // Create knock IO
        //
        knockIO = new KVPKnockIO(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                ksmooth,
                knockSchedule.VectTDate(0),        // obsDates
                knockSchedule.VectTDate(1),        // settleDates
                knockSchedule.VectDouble(2),        // barrierLos
                knockSchedule.VectDouble(3),        // barrierHis
                knockSchedule.VectDouble(4),        // rebates
                discZcName);

        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 * Effective observation dates are explicitly given in the schedule list
 */
int
_vpiparsSymbolSetKVPKnockIO2(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        char ioWindow,                        /* (I) Knock In/Out window */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */

        OVpiList *oknockSchedule,        /* (I) reset,acc, etc. scehdule  */

        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPKnockIO2";

        KKnockIO        kioType, kioWindow;
        KSmooth                ksmooth;

        KRate                *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KVPAtom                *under = NULL;
        KVPKnockIO        *knockIO = NULL;
        KVpiList &knockSchedule = *((KVpiList*) *oknockSchedule);


    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);

        //
        // Create knock IO
        //
        knockIO = new KVPKnockIO(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                ksmooth,
                knockSchedule.VectTDate(0),        // obsDates
                knockSchedule.VectTDate(1),        // obsEffDates
                knockSchedule.VectTDate(2),        // settleDates
                knockSchedule.VectDouble(3),        // barrierLos
                knockSchedule.VectDouble(4),        // barrierHis
                knockSchedule.VectDouble(5),        // rebates
                discZcName);

        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 */

int
_vpiparsSymbolSetKVPKnockIOSimple(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        char ioWindow,                        /* (I) Knock In/Out window */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */
        TDate startDate,                /* (I) This date is not included */
        TDate matDate,                        /* (I) */
        TDateInterval freq,                /* (I) */
        int nDays,                        /* (I) # of notifcation days */

        double barrierLo,                /* (I) Low barrier */
        double barrierHi,                /* (I) High barrier */
        double rebate,                        /* (I) Rebate */
        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPKnockIOSimple";

        KKnockIO        kioType, kioWindow;
        KSmooth                ksmooth;

        KRate                *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KDateInterval        notDays = KDateInterval(nDays, FALSE);
        TBoolean        stubAtEnd = FALSE;

        TDate        barrierD[2] = {startDate, matDate};
        KVector(TDate)        barrierDates(barrierD, barrierD+2); 

        KVector(double)        barrierLos(2, barrierLo); 
        KVector(double)        barrierHis(2, barrierHi); 
        KVector(double)        rebates(2, rebate); 

        KVPAtom                *under = NULL;
        KVPKnockIO        *knockIO = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);


        /* Create knock IO */
        knockIO = new KVPKnockIO(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                ksmooth,
                startDate,
                matDate,
                freq,
                stubAtEnd,
                notDays,
                barrierDates,
                barrierLos,
                barrierHis,
                rebates,
                discZcName);
        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 * Effective observation dates are implicitly given in the index rate
 */

int
_vpiparsSymbolSetKVPKnockIO1New(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        char ioWindow,                        /* (I) Knock In/Out window */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */

        OVpiList *oknockSchedule,        /* (I) reset,acc, etc. scehdule  */

        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPKnockIO1";

        KKnockIO        kioType, kioWindow;
        KSmooth                ksmooth;

        KRate                *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KVPAtom                *under = NULL;
        KVPKnockIO        *knockIO = NULL;
        KVpiList &knockSchedule = *((KVpiList*) *oknockSchedule);


    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);

        //
        // Create knock IO
        //
        knockIO = new KVPKnockIO(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                ksmooth,
                knockSchedule.VectTDate(0),        // obsDates
                knockSchedule.VectTDate(1),        // settleDates
                knockSchedule.VectDouble(2),        // barrierLos
                knockSchedule.VectDouble(3),        // barrierHis
                knockSchedule.VectDouble(4),        // rebates
                discZcName);

        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 * Effective observation dates are explicitly given in the schedule list
 */
int
_vpiparsSymbolSetKVPKnockIO2New(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        char ioWindow,                        /* (I) Knock In/Out window */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */

        OVpiList *oknockSchedule,        /* (I) reset,acc, etc. scehdule  */

        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPKnockIO2";

        KKnockIO        kioType, kioWindow;
        KSmooth                ksmooth;

        KRate                *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KVPAtom                *under = NULL;
        KVPKnockIO        *knockIO = NULL;
        KVpiList &knockSchedule = *((KVpiList*) *oknockSchedule);


    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);
        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);

        //
        // Create knock IO
        //
        knockIO = new KVPKnockIO(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                ksmooth,
                knockSchedule.VectTDate(0),        // obsDates
                knockSchedule.VectTDate(1),        // obsEffDates
                knockSchedule.VectTDate(2),        // settleDates
                knockSchedule.VectDouble(3),        // barrierLos
                knockSchedule.VectDouble(4),        // barrierHis
                knockSchedule.VectDouble(5),        // rebates
                discZcName);

        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a variable to a KVPKnockIO
 */

int
_vpiparsSymbolSetKVPKnockIOSimpleNew(
        TvpiparsSymbol *symb,                /* (O) */
        char ioType,                        /* (I) Knock In/Out */
        TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index */
        char smooth,                        /* (I) Node smoothing method */
        TDate startDate,                /* (I) This date is not included */
        TDate matDate,                        /* (I) */
        TDateInterval freq,                /* (I) */
        int nDays,                        /* (I) # of notifcation days */

        double rebate,                        /* (I) Rebate */
        TvpiparsSymbol *underSym,        /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "_vpiparsSymbolSetKVPKnockIOSimpleNew";

        KKnockIO        kioType;
        KSmooth         ksmooth;

        KIndicator      *floatRate= NULL;
        char            *curveName = NULL; // rate index curve

        KDateInterval   notDays = KDateInterval(nDays, FALSE);
        TBoolean        stubAtEnd = FALSE;

        TDate           barrierD[2] = {startDate, matDate};
        KVector(TDate)  barrierDates(barrierD, barrierD+2); 

        KVector(double) rebates(2, rebate); 

        KVPAtom         *under = NULL;
        KVPKnockIONew   *knockIO = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine << endl;

        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMINDICATOR) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KIndicator*) rateSym->ptr;

/*        cout << routine << '\t' 
             << "Name: " << floatRate->GetName() << '\t'
             << "rate type: " << rateSym->type << '\t' 
             << "curve name:" << floatRate->CurveName() << endl;
        cout << *floatRate << endl;
*/
        /* Create knock IO */
        knockIO = new KVPKnockIONew(
                symb->name,
                kioType,
                Raw2SharedPointer(floatRate),
                ksmooth,
                startDate,
                matDate,
                freq,
                stubAtEnd,
                notDays,
                barrierDates,
                rebates,
                discZcName);
        if (knockIO == NULL) return(FAILURE);

        /* Add underlying */
        knockIO->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



int	
_vpiparsSymbolSetKVPKnockIO2Idx(
	TvpiparsSymbol *symb,           /* (O) */
	char ioType,                    /* (I) Knock In/Out              */
	char ioWindow,                  /* (I) Knock In/Out window 1     */
	TvpiparsSymbol *rateSym,        /* (I) Knock In/Out rate index 1 */
    char ioWindow2,                 /* (I) Knock In/Out window 2     */
	TvpiparsSymbol *rateSym2,       /* (I) Knock In/Out rate index 2 */
	char smooth,              	    /* (I) Node smoothing method     */
	OVpiList *oknockSchedule,	    /* (I) obs, settle, barrLo, 
					                 *     barrHi, rebate scehdule   */
	TvpiparsSymbol *underSym,       /* (I) underlying bundle         */
	char *discZcName)
{
static        char routine[] = "vpiparsSymbolSetKVPKnockIO2Idx";

        KKnockIO        kioType, kioWindow, kioWindow2;
        KSmooth                ksmooth;

        KRate             *floatRate  = NULL;
        KRate             *floatRate2 = NULL;
        char              *curveName  = NULL; // rate index curve 1
        char              *curveName2 = NULL; // rate index curve 2

        KVPAtom           *under      = NULL;
        KVPKnockIO2Idx    *knockIO2Idx= NULL;
        KVpiList &knockSchedule       = *((KVpiList*) *oknockSchedule);


    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        // Knock in/out type
        //
        if (ioType == 'I')
                kioType = CRX_KNOCK_IN;
        else if (ioType == 'O')
                kioType = CRX_KNOCK_OUT;
        else if (ioType == 'N')
                kioType = CRX_NONE;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioType);

        // Knock in/out window type
        //
        if (ioWindow == 'I')
                kioWindow = CRX_KNOCK_IN;
        else if (ioWindow == 'O')
                kioWindow = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);

        // Same for second index
        if (ioWindow2 == 'I')
                kioWindow2 = CRX_KNOCK_IN;
        else if (ioWindow2 == 'O')
                kioWindow2 = CRX_KNOCK_OUT;
        else
                throw KFailure("%s: invalid knock in/out type (%c).\n",
                                routine,
                                ioWindow);

        // Smoothing type
        //
        if (smooth == 'D')
                ksmooth = DOUBLE_SMOOTH;
        else if (smooth == 'S')
                ksmooth = SINGLE_SMOOTH;
        else if (smooth == 'N')
                ksmooth = NO_SMOOTH;
        else
                throw KFailure("%s: invalid smoothing type (%c).\n",
                                routine,
                                smooth);

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;
        curveName = rateSym->strval;

        floatRate->SetCurveName(curveName);

        if (rateSym2->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym2->name);
        }
        floatRate2 = (KRate*) rateSym2->ptr;
        curveName2 = rateSym2->strval;

        floatRate2->SetCurveName(curveName2);


        //
        // Create knock IO
        //
        knockIO2Idx = new KVPKnockIO2Idx(
                symb->name,
                kioType,
                kioWindow,
                Raw2SharedPointer(floatRate),
                kioWindow2,
                Raw2SharedPointer(floatRate2),
                ksmooth,
                knockSchedule.VectTDate(0),        // obsDates
                knockSchedule.VectTDate(1),        // obsEffDates
                knockSchedule.VectTDate(2),        // settleDates
                knockSchedule.VectDouble(3),        // barrierLos
                knockSchedule.VectDouble(4),        // barrierHis
                knockSchedule.VectDouble(5),        // barrierLos
                knockSchedule.VectDouble(6),        // barrierHis
                knockSchedule.VectDouble(7),        // rebates
                discZcName);

        if (knockIO2Idx == NULL) return(FAILURE);

        /* Add underlying */
        knockIO2Idx->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) knockIO2Idx;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete knockIO2Idx;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



//**********************************************************************
//
// Floating leg constructors
//
//**********************************************************************


/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (simple)
 */

int
_vpiparsSymbolSetKVPFloatLegSimple(
        TvpiparsSymbol *symb,                /* (O) */
        TDate startDate,                /* (I) start date  */
        TDate matDate,                        /* (I) end date */
        TDateInterval freq,                /* (I) frequency */
        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) payment formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegSimple";


        TBoolean        stubAtEnd = FALSE;
        KRate                *floatRate= NULL;
        KVPFloatLeg        *floatLeg = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;

        /* Create floatLeg */
        floatLeg = new KVPFloatLeg(
                symb->name,
                startDate,
                matDate,
                KDateInterval(freq),
                dayCc,
                stubConv,
                stubAtEnd,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets)
 * Reset effective dates are specified in KRate.
 */

int
_vpiparsSymbolSetKVPFloatLeg1(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg1";

        TBoolean        stubAtEnd = FALSE;
        KRate                *floatRate= NULL;
        KVPFloatLeg        *floatLeg = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;

        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset effective
                flpSchedule.VectTDate( 1),        // acc start
                flpSchedule.VectTDate( 2),        // acc end
                flpSchedule.VectTDate( 3),        // pay
                flpSchedule.VectDouble(4),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets)
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg2(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg2";

        TBoolean        stubAtEnd = FALSE;
        KRate                *floatRate= NULL;
        KVPFloatLeg        *floatLeg = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;

        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets) with
 * arbitrary rate spread schedule. 
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg3(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg2";

        TBoolean        stubAtEnd = FALSE;
        KRate                *floatRate= NULL;
        KVPFloatLeg        *floatLeg = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;

        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // spread
                flpSchedule.VectDouble(6),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets) with
 * arbitrary formula schedule
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg4(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. schedule  */
                                        /*     including formula schedule */
        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg4";

        TBoolean        stubAtEnd = FALSE;
        KRate                *floatRate= NULL;
        KVPFloatLeg        *floatLeg = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        floatRate = (KRate*) rateSym->ptr;

        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(6),        // notional
                flpSchedule.VectString(5),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

//**********************************************************************
//
// Floating leg w/ multiple set of index rates
//
//**********************************************************************

/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (simple)
 */

int
_vpiparsSymbolSetKVPFloatLegNewSimple(
        TvpiparsSymbol *symb,                /* (O) */
        TDate startDate,                /* (I) start date  */
        TDate matDate,                        /* (I) end date */
        TDateInterval freq,                /* (I) frequency */
        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) payment formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNewSimple";


        TBoolean        stubAtEnd = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name 
             << " rate type " << rateSym->type << '\t';


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        /* Create floatLeg */
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;
            cout << "Name " << cplxfloatRate->GetName() << '\t'
                 << "Curve Name " << cplxfloatRate->CurveName() << endl;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                startDate,
                matDate,
                KDateInterval(freq),
                dayCc,
                stubConv,
                stubAtEnd,
                Raw2SharedPointer(cplxfloatRate),
                formula,
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                startDate,
                matDate,
                KDateInterval(freq),
                dayCc,
                stubConv,
                stubAtEnd,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve
        }
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets)
 * Reset effective dates are specified in KRate.
 */

int
_vpiparsSymbolSetKVPFloatLegNew1(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNew1";

        TBoolean        stubAtEnd = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        //
        // Create floatLeg 
        //
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset effective
                flpSchedule.VectTDate( 1),        // acc start
                flpSchedule.VectTDate( 2),        // acc end
                flpSchedule.VectTDate( 3),        // pay
                flpSchedule.VectDouble(4),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(cplxfloatRate),
                formula,
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset effective
                flpSchedule.VectTDate( 1),        // acc start
                flpSchedule.VectTDate( 2),        // acc end
                flpSchedule.VectTDate( 3),        // pay
                flpSchedule.VectDouble(4),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve
        }
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets)
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLegNew2(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNew2";

        TBoolean        stubAtEnd = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        //
        // Create floatLeg 
        //
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(cplxfloatRate),
                formula,
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve
        }
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets) with
 * arbitrary rate spread schedule. 
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLegNew3(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *formula,                        /* (I) formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNew3";

        TBoolean        stubAtEnd     = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        //
        // Create floatLeg 
        //
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // spread
                flpSchedule.VectDouble(6),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(cplxfloatRate),
                formula,
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // spread
                flpSchedule.VectDouble(6),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                formula,
                discZcName);   // discount curve

        }
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg (complex arbitrary resets) with
 * arbitrary formula schedule
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLegNew4(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. schedule  */
                                        /*     including formula schedule */
        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNew4";

        TBoolean        stubAtEnd     = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        //
        // Create floatLeg 
        //
        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(6),        // notional
                flpSchedule.VectString(5),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(cplxfloatRate),
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(6),        // notional
                flpSchedule.VectString(5),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                discZcName);   // discount curve

        }
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}
/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLegNew (complex arbitrary resets)
 * Reset effective dates are specified in KRate.
 */

int
_vpiparsSymbolSetKVPFloatLegNew5(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym,        /* (I) rate index */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLegNew5";

        TBoolean        stubAtEnd     = FALSE;
        KCplxRate      *cplxfloatRate = NULL;
        KRate          *floatRate     = NULL;
        KVPFloatLegNew *floatLeg      = NULL;

        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;
        cout << routine <<  " Setting symbol " << symb->name << endl;


        if (rateSym->type != SYMCPLXRATE &&
            rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }

        if (rateSym->type == SYMCPLXRATE)
        {
            cplxfloatRate = (KCplxRate*) rateSym->ptr;

            //
            // Create floatLeg by complex rate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                flpSchedule.VectString(6),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(cplxfloatRate),
                discZcName);   // discount curve
        }
        else
        {
            floatRate = (KRate*) rateSym->ptr;

            //
            // Create floatLeg by KRate
            //
            floatLeg = new KVPFloatLegNew(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                flpSchedule.VectString(6),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate),
                discZcName);   // discount curve

        }

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));

 //       delete floatRate;
        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


//**********************************************************************
//
// Floating leg w/ two set of index rates
//
//**********************************************************************


/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg2Idx (simple)
 */

int
_vpiparsSymbolSetKVPFloatLeg2IdxSimple(
        TvpiparsSymbol *symb,                /* (O) */
        TDate startDate,                /* (I) start date */
        TDate matDate,                        /* (I) end date */
        TDateInterval freq,                /* (I) frequency */
        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym1,        /* (I) rate index 1 */
        TvpiparsSymbol *rateSym2,        /* (I) rate index 2 */
        char *formula,                        /* (I) payment formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg2IdxSimple";


        TBoolean        stubAtEnd   = FALSE;

        KRate                *floatRate1 = NULL,
                        *floatRate2 = NULL;

        KVPFloatLeg2Idx        *floatLeg   = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym1->type != SYMRATE) { 
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym1->name);
        }
        if (rateSym2->type != SYMRATE) { 
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym2->name);
        }

        floatRate1 = (KRate*) rateSym1->ptr;
        floatRate2 = (KRate*) rateSym2->ptr;

        /* Create floatLeg */
        floatLeg = new KVPFloatLeg2Idx(
                symb->name,
                startDate,
                matDate,
                KDateInterval(freq),
                dayCc,
                stubConv,
                stubAtEnd,
                Raw2SharedPointer(floatRate1),
                Raw2SharedPointer(floatRate2),
                formula,
                discZcName);   // discount curve
        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg2Idx (complex arbitrary resets)
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg2Idx(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym1,        /* (I) rate index 1 */
        TvpiparsSymbol *rateSym2,        /* (I) rate index 2 */
        char *formula,                        /* (I) payment formula */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg2Idx";

        TBoolean        stubAtEnd   = FALSE;

        KRate                *floatRate1 = NULL,
                        *floatRate2 = NULL;

        KVPFloatLeg2Idx        *floatLeg   = NULL;


        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym1->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym1->name);
        }
        if (rateSym2->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym2->name);
        }

        floatRate1 = (KRate*) rateSym1->ptr;
        floatRate2 = (KRate*) rateSym2->ptr;

        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg2Idx(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate1),
                Raw2SharedPointer(floatRate2),
                formula,
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg2Idx (complex arbitrary resets)
 * Allows an arbitrary formula schedule
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg2Idx2(
        TvpiparsSymbol *symb,                /* (O) new symbol */

        OVpiList *oflpSchedule,                /* (I) reset,acc, etc. scehdule  */

        long dayCc,                        /* (I) day count convention */
        long stubConv,                        /* (I) stub convention */
        TvpiparsSymbol *rateSym1,        /* (I) rate index 1 */
        TvpiparsSymbol *rateSym2,        /* (I) rate index 2 */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg2Idx2";

        TBoolean        stubAtEnd   = FALSE;

        KRate                *floatRate1 = NULL,
                        *floatRate2 = NULL;

        KVPFloatLeg2Idx        *floatLeg   = NULL;


        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym1->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym1->name);
        }
        if (rateSym2->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym2->name);
        }

        floatRate1 = (KRate*) rateSym1->ptr;
        floatRate2 = (KRate*) rateSym2->ptr;
        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg2Idx(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(6),        // notional
                flpSchedule.VectString(5),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate1),
                Raw2SharedPointer(floatRate2),
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}




/*----------------------------------------------------------------------
 * Sets a variable to a KVPFloatLeg3Idx (complex arbitrary resets)
 * Allows an arbitrary formula schedule
 * Reset effective dates are given explicitly.
 */

int
_vpiparsSymbolSetKVPFloatLeg3Idx(
        TvpiparsSymbol *symb,            /* (O) new symbol */

        OVpiList *oflpSchedule,          /* (I) reset,acc, etc. scehdule  */

        long dayCc,                      /* (I) day count convention */
        long stubConv,                   /* (I) stub convention */
        TvpiparsSymbol *rateSym1,        /* (I) rate index 1 */
        TvpiparsSymbol *rateSym2,        /* (I) rate index 2 */
        TvpiparsSymbol *rateSym3,        /* (I) rate index 3 */
        char *discZcName)                /* (I) discount curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPFloatLeg3Idx";

        TBoolean        stubAtEnd   = FALSE;

        KRate           *floatRate1 = NULL,
                        *floatRate2 = NULL,
                        *floatRate3 = NULL;

        KVPFloatLeg3Idx        *floatLeg   = NULL;


        KVpiList &flpSchedule = *((KVpiList*) *oflpSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (rateSym1->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym1->name);
        }
        if (rateSym2->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym2->name);
        }

        if (rateSym3->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym3->name);
        }

        floatRate1 = (KRate*) rateSym1->ptr;
        floatRate2 = (KRate*) rateSym2->ptr;
        floatRate3 = (KRate*) rateSym3->ptr;
        //
        // Create floatLeg 
        //
        floatLeg = new KVPFloatLeg3Idx(
                symb->name,
                flpSchedule.VectTDate( 0),        // reset
                flpSchedule.VectTDate( 1),        // reset effective
                flpSchedule.VectTDate( 2),        // acc start
                flpSchedule.VectTDate( 3),        // acc end
                flpSchedule.VectTDate( 4),        // pay
                flpSchedule.VectDouble(5),        // notional
                flpSchedule.VectString(6),        // formula
                dayCc,
                stubConv,
                Raw2SharedPointer(floatRate1),
                Raw2SharedPointer(floatRate2),
                Raw2SharedPointer(floatRate3),
                discZcName);   // discount curve

        if (floatLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) floatLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete floatLeg;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Sets a variable to a KVPProtLeg (simple)
 */

int
_vpiparsSymbolSetKVPProtLegSimple(
        TvpiparsSymbol *symb,        /* (O) */
        TDate startDate,             /* (I) start date  */
        TDate endDate,               /* (I) end date */
        double notional,             /* (I) notional */
        char   *recovery,            /* (I) recovery */
        long   defConv,              /* (I) stub convention */
        char *discZcName)            /* (I) CDS curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPProtLegSimple";


        KProtPayConv    payType;
        KVPProtLeg      *protLeg = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        if (defConv == 0L)
            payType = PAY_DEF;
        else
            payType = PAY_MAT;

        /* Create floatLeg */
        protLeg = new KVPProtLeg(
                symb->name,
                startDate,
                endDate,
                notional,
                recovery,
                payType,
                discZcName);   // discount curve

        if (protLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) protLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Sets a variable to a KVPProtLeg with a schedule of notionals and dates
 */

int
_vpiparsSymbolSetKVPProtLegList1(
        TvpiparsSymbol *symb,        /* (O) */
        OVpiList *ontlSchedule,       /* (I) notional and dates, scehdule  */
        char   *recovery,            /* (I) recovery */
        long   defConv,              /* (I) pay at default or end */
        char *discZcName)            /* (I) CDS curve */
{
static        char routine[] = "_vpiparsSymbolSetKVPProtLegList1";


        KProtPayConv    payType;
        KVPProtLeg      *protLeg = NULL;

    try {
        KVpiList &ntlSchedule = *((KVpiList*) *ontlSchedule);
        dppLog << "Setting symbol " << symb->name << endl;

        if (defConv == 0L)
            payType = PAY_DEF;
        else
            payType = PAY_MAT;

        /* Create floatLeg */
        protLeg = new KVPProtLeg(
                symb->name,
                ntlSchedule.VectTDate( 0),         // start dates
                ntlSchedule.VectTDate( 1),         // end dates
                ntlSchedule.VectDouble(2),         // notional
                recovery,
                payType,
                discZcName);   // discount curve

        if (protLeg == NULL) return(FAILURE);

        symb->type = SYMVPBASE;
        symb->ptr = (void*) protLeg;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Sets a variable to a KVPDefProtect (simple)
 */

int
_vpiparsSymbolSetKVPDefExposureSimple(
        TvpiparsSymbol *symb,           /* (O) */
        TDate startDate,                /* (I) start date          */
        TDate endDate,                  /* (I) end date            */
        TDateInterval freq,             /* (I) freq                */
        double   recovery,              /* (I) recovery rate       */
        TvpiparsSymbol *underSym,       /* (I) underlying bundle   */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetDefExposureSimple";

        KVPAtom         *under    = NULL;
        KVPDefProtect   *defProt  = NULL;

    try {
        dppLog << "Setting symbol " << symb->name << endl;


        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create default protection */
        defProt = new KVPDefProtect(
                symb->name,
                startDate,
                endDate,
                freq,
                DEF_EXPOSURE,
                recovery,
                discZcName);

        if (defProt == NULL) return(FAILURE);

        /* Add underlying */
        defProt->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) defProt;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete defProt;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}


/*----------------------------------------------------------------------
 * Default exposure with with arbitrary schedule
 */

int
_vpiparsSymbolSetKVPDefExposureGeneral(
        TvpiparsSymbol *symb,           /* (O) */
        OVpiList *oDefProtSchedule,     /* (I) start,end,settle scehdule */
        double   recovery,              /* (I) recovery rate             */
        TvpiparsSymbol *underSym,       /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPDefExposureGeneral";

        KVPAtom              *under = NULL;
        KVPDefProtect        *defProt = NULL;
        KVpiList &defProtSchedule = *((KVpiList*) *oDefProtSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        // Get underlying 
        //
        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create option */
        defProt = new KVPDefProtect(
                symb->name,
                defProtSchedule.VectTDate(0),    // startDates
                defProtSchedule.VectTDate(1),    // endDates
                defProtSchedule.VectTDate(2),    // settleDates
                DEF_EXPOSURE,
                recovery,
                discZcName);
        if (defProt == NULL) return(FAILURE);

        /* Add underlying */
        defProt->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) defProt;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete defProt;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



/*----------------------------------------------------------------------
 * Default knock-in with arbitrary schedule and rebates)
 */

int
_vpiparsSymbolSetKVPDefKnockIn(
        TvpiparsSymbol *symb,           /* (O) */
        OVpiList *oDefProtSchedule,     /* (I) start,end,settle,rebate sch */
        TvpiparsSymbol *underSym,       /* (I) underlying bundle */
        char *discZcName)               /* (I) discount curve name */

{
static        char routine[] = "vpiparsSymbolSetKVPDefKnockIn";


        KVPAtom              *under = NULL;
        KVPDefProtect        *defProt = NULL;
        KVpiList &defProtSchedule = *((KVpiList*) *oDefProtSchedule);

    try {
        dppLog << "Setting symbol " << symb->name << endl;

        // Get underlying 
        //
        if (underSym->type != SYMVPBASE) {
                throw KFailure("%s: variable `%s' is not a stream.\n",
                        routine, underSym->name);
        }
        under  = (KVPAtom*) underSym->ptr;

        /* Create option */
        defProt = new KVPDefProtect(
                symb->name,
                defProtSchedule.VectTDate(0),    // startDates
                defProtSchedule.VectTDate(1),    // endDates
                defProtSchedule.VectTDate(2),    // settleDates
                defProtSchedule.VectDouble(3),   // rebates
                DEF_KNOCKIN,
                0e0,           // Not used for def knock-in
                discZcName);
        if (defProt == NULL) return(FAILURE);

        /* Add underlying */
        defProt->AddDep(Raw2SharedPointer(under));

        symb->type = SYMVPBASE;
        symb->ptr = (void*) defProt;
        symb->setFlag = TRUE;

        // Add instrument to root bundle
        vpArray->insert(vpArray->end(), 
                        Raw2SharedPointer((KVPAtom*)symb->ptr));


        return(SUCCESS);

    }
    catch (KFailure) {
        delete defProt;
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}



//**********************************************************************
//
// Resets
//
//**********************************************************************

/*----------------------------------------------------------------------
 * Adds a rate reset
 */

int
_vpiparsAddRateResetRank(
        TvpiparsSymbol *rateSym,        /* (I) underlying bundle */
        TDate resetDate,                /* (I) */
        double value)                        /* (I) */
{
static        char routine[] = "_vpiparsAddRateResetRank";

    try {

        if (rateSym->type != SYMRATE) {
                throw KFailure("%s: variable `%s' is not a rate.\n",
                        routine, rateSym->name);
        }
        KRate        *floatRate = (KRate*) rateSym->ptr;
        KRate        floatRateNS;
        
        // Create a copy
        floatRateNS = *floatRate;
        floatRateNS.SetSpread(0e0);

        // Add reset to bank
        sVpResetBank->Insert(floatRateNS, resetDate, value);

        return(SUCCESS);

    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}








/*----------------------------------------------------------------------
 * Frees.
 */

int
_vpiparsSymbolFree(
        TvpiparsSymbol *symb)                /* (I) */
{
    static char routine[] = "vpiparsSymbolFree";

    try {
        switch (symb->type) {
        case DRL_STRING_T:
                delete [] ((char*)symb->ptr);
                break;
        case DRL_DOUBLE_T:
                /* Nothing to do! */
                break;
        default:
                DppErrMsg("_vpiparsSymbolDelete: undefined.\n");
                return(FAILURE);
        }
        
        return (TRUE);

    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
}

/*----------------------------------------------------------------------
 * Prints a symbol.
 */

int
_vpiparsSymbolPrint(TvpiparsSymbol *symb)
{
static        char routine[] = "_vpiparsSymbolPrint";
        KVPAtom                *vp;

    try {

        if (symb->setFlag == FALSE) {
            DppErrMsg("%s: undefined variable `%s'.\n", routine, symb->name);
            throw KFailure();
        }

        dppLog << format("Symbol `%s' [%s] = ",
                        symb->name,
                        TSymTypeName(symb->type));

        switch (symb->type) {
        case SYMSTRING:
                dppLog << format("\"%s\"\n", ((char*)symb->ptr));
                break;
        case SYMDOUBLE:
                dppLog << format("%12.8f\n", symb->dval);
                break;
        case SYMTDATE:
                dppLog << format("%10s\n", DrlTDatePrint(NULL, symb->dtval));
                break;

        case SYMDATEINT:
                {
                        KDateInterval *vpDtIntvl = ((KDateInterval*) symb->ptr);
                        dppLog << *vpDtIntvl << endl;
                }
                break;
        case SYMRATE:
                {
                        KRate *vpRate = ((KRate*) symb->ptr);
                        dppLog << *vpRate << endl;
                        dppLog << format("\"%s\"\n", symb->strval); 
                }
                break;
        case SYMCPLXRATE:
                {
                        KCplxRate *vpCplxRate = ((KCplxRate*) symb->ptr);
                        dppLog << *vpCplxRate << endl;
                }
                break;
        case SYMINDICATOR:
                {
                        KIndicator *vpIndicator = ((KIndicator*) symb->ptr);
                        dppLog << *vpIndicator << endl;
                }
                break;
        case SYMVPBASE:
                vp = ((KVPAtom*) symb->ptr);
                dppLog << vp->GetName() << "(" << vp->TypeName() << ')' << endl;
                dppLog << *vp << endl;
                break;

        default:
                DppErrMsg("%s: cannot print variable `%s' is of type %s.\n",
                        routine, symb->name,
                        TSymTypeName(symb->type));
                throw KFailure();
        }
    }
    catch (KFailure) {
        DppErrMsg("%s: failed.\n", routine);
        return(FAILURE);
    }
    return(SUCCESS);
}

/*----------------------------------------------------------------------
 * 
 */

char*
TSymTypeName(TSymType type)
{
        switch (type) {
        case SYMINT:        return("SYMINT");
        case SYMDOUBLE:     return("SYMDOUBLE");
        case SYMTDATE:      return("SYMTDATE");
        case SYMLONG:       return("SYMLONG");
        case SYMSTRING:     return("SYMSTRING");
        case SYMDATEINT:    return("SYMDATEINT");
        case SYMRATE:       return("SYMRATE");
        case SYMCPLXRATE:   return("SYMCPLXRATE");
        case SYMINDICATOR:  return("SYMINDICATOR");
        case SYMVPBASE:     return("SYMVPBAS");
        case SYMNIL:        return("SYMNIL");
        }
        return("unknown");
}


/*----------------------------------------------------------------------*/
};        /* extern "C" */




