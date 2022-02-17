//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLTest.hpp
//
//   Description : Creating/running Tests for Addins and XL Kit
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLAddin.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/XLTest.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/Format.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ErrorHandler.hpp"

DRLIB_BEGIN_NAMESPACE

/* spoof type for handles - note larger than a short */
#define xltypeEdgHandle     0x10000

/** Class for wrapping the opers supplied from excel to our addin
    function together with the other data needed in order to be able
    to create a regression test */
class ExcelParams: public CObject,
                   public IRegressionTest{
public:
    /** wrapper class for XL OPER */
    class EDROper: public CObject{
    public:
        static CClassConstSP const TYPE;
    private: //fields
        XL_OPER oper; // $unregistered
    private: // methods
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            REGISTER(EDROper, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultEDROper);
        }

        static IObject* defaultEDROper(){
            // would like to do: OPER tmp = {}; but doesn't compile
            XL_OPER tmp;
            memset(&tmp, 0, sizeof(XL_OPER)); // clear memory - avoid UMR
            return new EDROper(tmp);
        }

        static void operArrayFree(XL_OPER* opers, int length){
            if (opers){
                for (int i = 0; i < length; i++){
                    XLConvert::operFreeSimple(opers[i]);
                }
                FREE(opers);
            }
        }

        /** write array of OPERs out  */
        static void write(const XL_OPER*   opers,
                          int              numOpers,
                          const string&    tag,
                          Writer*          writer) {
            try{
                EDROperArray array(numOpers);
                for (int i = 0; i < numOpers; i++){
                    array[i] = EDROperSP(new EDROper(opers[i]));
                }
                array.write(tag, writer);
            } catch (exception& e){
                throw ModelException(e, "ExcelParams::write");
            }
        }


        /** read an array of OPERs from reader */
        static XL_OPER* read(Reader::Node*   elem,
                             int&            length,
                             Reader*         reader){
            static const string method("ExcelParams::read");
            XL_OPER*  opers = 0;
            length = 0;
            try{
                IObjectSP input(reader->read(elem));
                EDROperArray& array = dynamic_cast<EDROperArray&>(*input);
                length = array.size();
                opers = NEW_ARRAY(XL_OPER, length);
                for (int i = 0; i < length; i++){
                    operCopy(array[i]->getOper(), &opers[i]);
                }
            } catch (exception &e) {
                operArrayFree(opers, length);
                throw ModelException(&e, method);
            }
            return opers;
        }

    public:
        /** 'simple' copy of an OPER. Assumes newOper is 'empty' */
        static void operCopy(const XL_OPER *oper, XL_OPER *newOper){
            static const string  routine("ExcelParams::operCopy");
            try{
                *newOper = *oper; /* structure copy - need to be careful */
                if (oper->type & xltypeMulti){
                    int numOper = oper->val.xlarray.rows *
                        oper->val.xlarray.columns;
                    newOper->val.xlarray.lparray = 0; // avoid memory problems
                    newOper->val.xlarray.lparray = NEW_ARRAY(XL_OPER, numOper);
                    for (int idx = 0; idx < numOper; idx++) {
                        XL_OPER *currSubOper = &oper->val.xlarray.lparray[idx];
                        XL_OPER *newSubOper  = &newOper->val.xlarray.lparray[idx];
                        if (currSubOper->type & xltypeMulti) {
                            throw ModelException(routine, "Arrays of arrays"
                                                 " not supported");
                        }
                        *newSubOper = *currSubOper;  /* structure copy */
                        if (currSubOper->type & xltypeStr) {
                            /* copy counted string */
                            newSubOper->val.str =
                                NEW_ARRAY(char, currSubOper->val.str[0] + 1);
                            memcpy(newSubOper->val.str, currSubOper->val.str,
                                   currSubOper->val.str[0] + 1);
                        }
                    }
                } else if (oper->type & xltypeStr) {
                    newOper->val.str = 0; // avoid memory problems
                    /* copy counted string */
                    newOper->val.str = NEW_ARRAY(char, oper->val.str[0] + 1);
                    memcpy(newOper->val.str, oper->val.str,oper->val.str[0]+1);
                } else if (oper->type & xltypeRef){
                    newOper->val.mref.lpmref = 0; // avoid memory problems
                    newOper->val.mref.lpmref = NEW(XLMREF);
                    // then simple structure copy
                    *newOper->val.mref.lpmref = *oper->val.mref.lpmref;
                }
            } catch (exception& e){
                XLConvert::operFreeSimple(*newOper);
                throw ModelException(e, routine);
            }
        }

        EDROper(const XL_OPER& oper): CObject(TYPE){
            operCopy(&oper, &this->oper);
        }

        virtual IObject* clone() const{
            return new EDROper(oper);
        }

        const XL_OPER* getOper() const{
            return &oper;
        }

        ~EDROper(){
            XLConvert::operFreeSimple(oper);
        }

        /** write object out to writer */
        virtual void write(const string& tag, Writer* writer) const{
            static const string routine("ExcelParams::write");
            try{
                unsigned long  xlType = oper.type; // default
                IObjectSP object;
                writer->objectStart(tag, "", this, true); /* can igore return -
                                                             no SPs in opers */
                if (XLConvert::operContainsHandle(oper)){
                    // is it really a handle
                    if (Handle::exists(
                        XLConvert::operStringToCPPString(oper)))
                    {
                        object = XLConvert::getObjectFromHandle(oper, NULL);
                        /* xltypeEdgHandle is a spoof type */
                        xlType = xltypeEdgHandle | (oper.type & xlbitDLLFree);
                    }
                }
                CIntSP eqXlType(CInt::create(xlType));
                eqXlType->write("xlType", writer);
                xlType = xlType & (~xlbitDLLFree); // remove xlbitDLLFree bit
                IObjectSP subObj;
                switch (xlType) {
                case xltypeNum:
                    subObj = IObjectSP(CDouble::create(oper.val.num));
                    break;
                case xltypeBool:
                    subObj = IObjectSP(CInt::create(oper.val.boolVal));
                    break;
                case xltypeErr:
                    subObj = IObjectSP(CInt::create(oper.val.err));
                    break;
                case xltypeMissing:
                case xltypeNil:
                    // do nothing
                    break;
                case xltypeEdgHandle:
                {
                    object->write("handle", writer);
                    /* save user handle part of name */
                    string name = XLConvert::operStringToCPPString(oper);
                    subObj = IObjectSP(CString::create(
                        Handle::getUserHandleName(name)));
                    break;
                }
                case xltypeStr:
                {
                    string value = XLConvert::operStringToCPPString(oper);
                    subObj = IObjectSP(CString::create(value));
                    break;
                }
                case xltypeMulti:
                {
                    IObjectSP rows(CInt::create(oper.val.xlarray.rows));
                    IObjectSP columns(CInt::create(oper.val.xlarray.columns));
                    rows->write("rows", writer);
                    columns->write("columns", writer);
                    write(oper.val.xlarray.lparray,
                          oper.val.xlarray.rows*oper.val.xlarray.columns,
                          "lparray", writer);
                    break;
                }
                case xltypeRef:
                case xltypeSRef:
                {
                    if (xlType == xltypeRef){
                        IObjectSP idSheet(CInt::create(oper.val.mref.idSheet));
                        idSheet->write("idSheet", writer);
                    }
                    const XLREF& xlRef = xlType == xltypeSRef?
                        oper.val.sref.ref: *oper.val.mref.lpmref->reftbl;
                    IObjectSP rwFirst(CInt::create(xlRef.rwFirst));
                    IObjectSP rwLast(CInt::create(xlRef.rwLast));
                    IObjectSP colFirst(CInt::create(xlRef.colFirst));
                    IObjectSP colLast(CInt::create(xlRef.colLast));
                    rwFirst->write("rwFirst", writer);
                    rwLast->write("rwLast", writer);
                    colFirst->write("colFirst", writer);
                    colLast->write("colLast", writer);
                    break;
                }
                default:
                    throw ModelException(routine, "Unrecognised xl type "+
                                         Format::toString(oper.type));
                }

                if (subObj.get()){
                    subObj->write("val", writer);
                }
                writer->objectEnd(tag, this);
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }

        /** populate an empty object from reader */
        virtual void import(Reader::Node* elem, Reader* reader){
            static const string method("ExcelParams::import");
            string              fieldname;
            try {
                XL_OPER*  operArray = 0;
                int    numOpers = 0;
                HashtableSP dd(new Hashtable());
                try {
                    Reader::NodeListSP nl(elem->children());
                    for (unsigned int i = 0; i < nl->size(); i++) {
                        Reader::Node*  child = (*nl)[i].get();
                        fieldname = child->name();
                        if (fieldname == "lparray"){
                            operArray = read(child, numOpers, reader);
                        }
                        IObjectSP obj(reader->read(child));
                        dd->put(fieldname, obj);
                    }
                }
                catch (exception& e) {
                    throw ModelException(e, method, "Failed whilst "
                                         "processing: "+fieldname);
                }
                // now build up oper from components
                CIntSP typeObject(CIntSP::dynamicCast(dd->get("xlType")));
                int xlType = typeObject->intValue();
                oper.type = (unsigned short) xlType;
                xlType &= ~xlbitDLLFree; /* switch off xlbitDLLFree */

                IObjectSP val;
                // pull out the val field if we ought to have one
                if (xlType == xltypeNum || xlType == xltypeStr ||
                    xlType == xltypeBool || xlType == xltypeErr ||
                    xlType == xltypeEdgHandle){
                    val = dd->get("val");
                    if (!val){
                        throw ModelException(method, "Component 'val' is "
                                             "missing!");
                    }
                }

                switch (xlType){
                case xltypeNum:
                {
                    CDoubleSP data(CDoubleSP::dynamicCast(val));
                    oper.val.num = data->doubleValue();
                    break;
                }
                case xltypeEdgHandle:
                {
                    int   origType   = oper.type;
                    // retrieve the actual handle
                    IObjectSP object = dd->get("handle");
                    /* actual name doesn't matter too much as long we
                       avoid clashes - this is our attempt to do that */
                    static int handleCount = 0;
                    string handleSuffix(Format::toString(handleCount));
                    handleCount = (handleCount + 1) %  10000;
                    CStringSP data(CStringSP::dynamicCast(val));
                    XLConvert::createHandle(
                        data->stringValue(),
                        handleSuffix, // append count on end
                        object, 0, -1, // don't append a number
                        oper);
                    oper.type |= (origType & xlbitDLLFree); // recover dllfree
                    break;
                }
                case xltypeStr:
                {
                    CStringSP data(CStringSP::dynamicCast(val));
                    XLConvert::populateOperFromString(data->stringValue(),
                                                      oper);
                    break;
                }
                case xltypeBool:
                {
                    CIntSP data(CIntSP::dynamicCast(val));
                    oper.val.boolVal = (unsigned short) data->intValue();
                    break;
                }
                case xltypeErr:
                {
                    CIntSP data(CIntSP::dynamicCast(val));
                    oper.val.err = (unsigned short) data->intValue();
                    break;
                }
                case xltypeMulti:
                {
                    CIntSP rows(CIntSP::dynamicCast(dd->get("rows")));
                    CIntSP columns(CIntSP::dynamicCast(dd->get("columns")));
                    oper.val.xlarray.rows = (unsigned short) rows->intValue();
                    oper.val.xlarray.columns = (unsigned short)
                        columns->intValue();
                    oper.val.xlarray.lparray = operArray;
                    if (numOpers !=
                        oper.val.xlarray.columns * oper.val.xlarray.rows) {
                        throw ModelException("Mismatch between actual number"
                                             "in array and given number");
                    }
                    break;
                }
                case xltypeRef:
                case xltypeSRef:
                {

                    if (xlType == xltypeRef){
                        CIntSP idSheet(CIntSP::dynamicCast(dd->
                                                           get("idSheet")));
                        oper.val.mref.idSheet = idSheet->intValue();
                        oper.val.mref.lpmref = NEW(XLMREF);
                        oper.val.mref.lpmref->count = 1;
                    } else {
                        oper.val.sref.count = 1;
                    }
                    XLREF& xlRef = xlType == xltypeSRef?
                        oper.val.sref.ref: *oper.val.mref.lpmref->reftbl;
                    CIntSP rwFirst(CIntSP::dynamicCast(dd->get("rwFirst")));
                    CIntSP rwLast(CIntSP::dynamicCast(dd->get("rwLast")));
                    CIntSP colFirst(CIntSP::dynamicCast(dd->get("colFirst")));
                    CIntSP colLast(CIntSP::dynamicCast(dd->get("colLast")));
                    xlRef.rwFirst = rwFirst->intValue();
                    xlRef.rwLast = rwLast->intValue();
                    xlRef.colFirst = colFirst->intValue();
                    xlRef.colLast = colLast->intValue();
                    break;
                }
                case xltypeMissing:
                case xltypeNil:
                    /* do nothing */
                    break;
                default:
                    throw ModelException("Unrecognised xl type "+
                                         Format::toString(oper.type));
                }
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

        /** write object out in 'output' format - ie suitable for comparing
            regression files with */
        virtual void outputWrite(const string& linePrefix,
                                 const string& prefix,
                                 ostream&      stream) const{
            static const string routine("ExcelParams::outputWrite");
            try{
                IObjectSP subObj;
                string prefixToUse = prefix.empty()? "OPER": prefix;
                if (XLConvert::operContainsHandle(oper)){
                    subObj = XLConvert::getObjectFromHandle(oper, NULL);
                    prefixToUse = prefix.empty()?
                        subObj->getClass()->getName():
                        prefix+"_"+subObj->getClass()->getName();
                } else {
                    switch (oper.type & (~xlbitDLLFree)) {
                    case xltypeNum:
                        subObj = IObjectSP(CDouble::create(oper.val.num));
                        break;
                    case xltypeBool:
                        subObj = IObjectSP(CInt::create(oper.val.boolVal));
                        break;
                    case xltypeErr:
                        subObj = IObjectSP(CString::create(
                            "Error code "+Format::toString(oper.val.err)));
                        break;
                    case xltypeMissing:
                        subObj = IObjectSP(CString::create("xlType Missing"));
                        break;
                    case xltypeNil:
                        subObj = IObjectSP(CString::create("xlType Nil"));
                        break;
                    case xltypeStr:
                    {
                        string value = XLConvert::operStringToCPPString(oper);
                        subObj = IObjectSP(CString::create(value));
                        break;
                    }

                    case xltypeMulti:
                    {
                        for (int i = 0; i < oper.val.xlarray.columns; i++)
                        {
                            for (int j = 0; j < oper.val.xlarray.rows; j++)
                            {
                                EDROper tmp(*(
                                    oper.val.xlarray.lparray +
                                    XLConvert::offset(i, j,
                                                      oper.val.xlarray.columns,
                                                      oper.val.xlarray.rows)));
                                string label("["+Format::toString(i)+"]"+
                                             "["+Format::toString(j)+"]");
                                tmp.outputWrite(linePrefix, prefixToUse + "_" +
                                                label, stream);
                            }
                        }
                        break;
                    }
                    case xltypeRef:
                    case xltypeSRef:
                    {
                        if (oper.type == xltypeRef){
                            IObjectSP idSheet(CInt::
                                              create(oper.val.mref.idSheet));
                            idSheet->outputWrite(linePrefix, prefixToUse+
                                                 "_idSheet", stream);
                        }
                        const XLREF& xlRef = oper.type == xltypeSRef?
                            oper.val.sref.ref: *oper.val.mref.lpmref->reftbl;
                        IObjectSP rwFirst(CInt::create(xlRef.rwFirst));
                        IObjectSP rwLast(CInt::create(xlRef.rwLast));
                        IObjectSP colFirst(CInt::create(xlRef.colFirst));
                        IObjectSP colLast(CInt::create(xlRef.colLast));
                        rwFirst->outputWrite(linePrefix,
                                             prefixToUse+"_rwFirst", stream);
                        rwLast->outputWrite(linePrefix,
                                            prefixToUse+"_rwLast", stream);
                        colFirst->outputWrite(linePrefix,
                                              prefixToUse+"_colFirst", stream);
                        colLast->outputWrite(linePrefix,
                                             prefixToUse+"_colLast", stream);
                        break;
                    }
                    default:
                        throw ModelException(routine, "Unrecognised xl type "+
                                             Format::toString(oper.type));
                    }
                }
                if (subObj.get()){
                    subObj->outputWrite(linePrefix, prefixToUse, stream);
                }

            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
    };

    // smart pointer support
    typedef smartPtr<EDROper> EDROperSP;
    // support for arrays of EDROper (note array of smart pointers)
    typedef array<EDROperSP, EDROper> EDROperArray;



private:  ////// fields ///////

    string           addinName;          /* what addin to invoke */
    EDROperArray     params;             /* the parameters */
    bool             calledFromVB;       /* have we been called from VB */
    EDROperArray     resolvedRefParams;  /* for xltypeSRef/xltypeRef opers.
                                            Length same as params */

    // allow ability to emulate excel when running through regression tester
    static bool          staticCalledFromVB;
    static EDROperArray  staticRefParams; /* opers which are
                                             xltypeSRef/xltypeRef */
    static EDROperArray  staticResolvedRefParams; // "xlcoerce(refParams)"
    static vector<BoolArray> staticRecalculated;  /* whether each cell in
                                                     staticRefParams has been
                                                     'calculated' */
private:     /// methods ///

    static bool callerIsVB(){
        return staticCalledFromVB;
    }

    /** for reflection */
    ExcelParams(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ExcelParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultExcelParams);
        FIELD(addinName, "Identifies addin function");
        FIELD(params, "Excel input opers");
        FIELD(calledFromVB, "Whether called from VB or not");
        FIELD(resolvedRefParams, "Resolved references");
        FIELD_MAKE_OPTIONAL(resolvedRefParams); // backwards compatibility
        // set default method for determining (this is subsequently overwritten
        // in xlAutoOpen)
        XLAddin::setCallerIsVBMethod(callerIsVB);
        // same thing for coerce methods
        XLConvert::setXLCoerceMethods(coerceToMulti, xlFreeCoerce);
    }

    static IObject* defaultExcelParams(){
        return new ExcelParams();
    }

    /* wrapper around XLAddin::execute() - turns parameter list into a
       va_list */
    static XL_OPER * xlAddinWrap(const XLAddin* xlAddin,
                                 const XL_OPER *a0, ...) {
        XL_OPER       *output;
        va_list     args;
        va_start(args, a0);
        output = xlAddin->execute(a0, args);
        va_end(args);
        return output;
    }

public:
    static CClassConstSP const TYPE;

    /** Creates new ExcelParams, takes deep copy of supplied OPERs. The
        return value of cachedRefs must be checked - if false then there are
        uncalculated cells and control must be returned to Excel */
    ExcelParams(const string&  addinName,
                va_list        paramsIn,    /* the parameters */
                int            numParams,
                bool           calledFromVB,
                bool&          cachedRefs):
        CObject(TYPE), addinName(addinName), params(numParams),
        calledFromVB(calledFromVB), resolvedRefParams(numParams){
        try{
            cachedRefs = true;
            for (int i = 0; i < numParams; i++){
                XL_OPER *currentOper = va_arg(paramsIn, XL_OPER *);
                EDROperSP oper(new EDROper(*currentOper));
                params[i] = oper;
                EDROperSP resolvedOper;
                // if oper is a ref/sref then look up value and store
                if ((currentOper->type & xltypeRef) ||
                    (currentOper->type & xltypeSRef)){
                    XL_OPER coerceOutput;
                    if (XLConvert::coerceToMulti(*currentOper,
                                                 coerceOutput)==xlretUncalced){
                        cachedRefs = false;
                        return; // exit immediately
                    }
                    resolvedOper = EDROperSP(new EDROper(coerceOutput));
                    XLConvert::freeCoerce(coerceOutput);
                }
                resolvedRefParams[i] = resolvedOper;
            }
        } catch (exception& e){
            throw ModelException(e, "ExcelParams");
        }
    }

    /** Run regression test for inputs stored within this class */
    IObjectSP runTest() const{
        static const string routine("ExcelParams::runTest");
        XL_OPER     *output    = 0;
        IObjectSP result;
        try{
            /* find addin type */
            const XLAddin* xlAddin = XLAddin::lookUp(addinName);
            // check we've got the right number of parameters
            int numParamsPassed = params.size();
            int numParamsInAddin = xlAddin->numParams();
            //int numParams = xlAddin->numParams();
            int numOptionalParams = xlAddin->numOptionalParams();
            if ((numParamsPassed < (numParamsInAddin - numOptionalParams)) ||
                (numParamsPassed > numParamsInAddin))
            {
                throw ModelException(
                    routine, "Function " +addinName+ " takes "+
                    Format::toString(numParamsInAddin - numOptionalParams)+
                    " to " +Format::toString(numParamsInAddin)+
                    " parameters but "+Format::toString(params.size())+
                    " were supplied");
            }

            /* set static 'VB' flag */
            staticCalledFromVB = calledFromVB;
            staticRefParams.clear();
            staticResolvedRefParams.clear();
            staticRecalculated.clear();
            // now need to arrange OPERs in same was as excel. This
            // supports 33 - more than Excel apparently.
            const XL_OPER* inputs[33] = {0};
            int maxNumCalls = 1; // how many cells are in references + 1
            for (int i = 0; i < numParamsPassed; i++){
                inputs[i] = params[i]->getOper();
                if (inputs[i]->type & (xltypeSRef | xltypeRef)){
                    // set up static cache (for test harness)
                    staticRefParams.push_back(params[i]);
                    staticResolvedRefParams.push_back(resolvedRefParams[i]);
                    // update maxNumCalls
                    const XLREF& xlRef = (inputs[i]->type & xltypeSRef)?
                        inputs[i]->val.sref.ref:
                        *inputs[i]->val.mref.lpmref->reftbl;
                    int numCells = (xlRef.rwLast-xlRef.rwFirst+1)*
                        (xlRef.colLast-xlRef.colFirst+1);
                    staticRecalculated.push_back(BoolArray(numCells));
                    maxNumCalls +=  numCells;
                }
            }
            XL_OPER missingVal;
            missingVal.type = xltypeMissing;
            for (int j = numParamsPassed; j < numParamsInAddin; j++){
                inputs[j] = &missingVal;
            }

            /* now a little subtlety. We now pretend that cells are
               uncalculated and wait until the kit calls us back to
               explicitly calculate a cell. So we loop until we get
               a result (or too many loops) */
            int numCalls = 0;
            do{
                /* Note we are deliberately ignoring the numParams field
                   and we are almost certainly passing illegal addresses
                   to ExcelMainWrap - but EdgExcelMain knows how many
                   parameters it is allowed to read so it's okay - it
                   won't read the bogus addresses */
                output = xlAddinWrap(xlAddin,
                                     inputs[0],  inputs[1], inputs[2],
                                     inputs[3], inputs[4], inputs[5],
                                     inputs[6], inputs[7], inputs[8],
                                     inputs[9], inputs[10], inputs[11],
                                     inputs[12], inputs[13], inputs[14],
                                     inputs[15], inputs[16], inputs[17],
                                     inputs[18], inputs[19], inputs[20],
                                     inputs[21], inputs[22], inputs[23],
                                     inputs[24], inputs[25], inputs[26],
                                     inputs[27], inputs[28], inputs[29],
                                     inputs[30], inputs[31], inputs[32]);
                numCalls++;
            } while (!output && numCalls <= maxNumCalls);
            /* incidentally ,replacing each of the 'inputs[<n>]' with
               '*(inputs++)' doesn't work */
            if (!output){
                // shouldn't happen unless no memory
                throw ModelException(routine, "XLAddin::execute returned null"
                                     " oper!");
            }
            result = IObjectSP(new EDROper(*output));
        } catch (exception& e){
            if (output && (output->type & xlbitDLLFree)) {
                XLConvert::operFreeSimple(*output);
                FREE(output);
                // avoid calling xlAutoFree(output); - EAS link issues
            }
            throw ModelException(e, routine);
        }
        if (output && (output->type & xlbitDLLFree)) {
            XLConvert::operFreeSimple(*output);
            FREE(output);
            // avoid calling xlAutoFree(output); - EAS link issues
        }
        return result;
    }

    /** Requests conversion of XLOPER to type xlmulti. If successful, output
        must be passed to XLFree */
    static int coerceToMulti(const XLOPER& input, XLOPER& output){
        // identify which oper in staticRefParams corresponds to request
        for (int i = 0; i < staticRefParams.size(); i++){
            const XL_OPER* oper = staticRefParams[i]->getOper();
            if ((input.type & oper->type & xltypeSRef) ||
                ((input.type & oper->type & xltypeRef) &&
                 input.val.mref.idSheet == oper->val.mref.idSheet)){
                // now see if requested range is within oper
                const XLREF& xlRef = (oper->type & xltypeSRef)?
                    oper->val.sref.ref: *oper->val.mref.lpmref->reftbl;
                const XLREF& xlRefInput = (oper->type & xltypeSRef)?
                    input.val.sref.ref: *input.val.mref.lpmref->reftbl;
                if (xlRefInput.rwFirst >= xlRef.rwFirst &&
                    xlRefInput.colFirst >= xlRef.colFirst &&
                    xlRefInput.rwLast <= xlRef.rwLast &&
                    xlRefInput.colLast <= xlRef.colLast){
                    // have found a match - look up resolved oper
                    const XL_OPER* multiOper =
                        staticResolvedRefParams[i]->getOper();
                    // then create xltypeMulti oper of right dimension
                    output.type = xltypeMulti;
                    output.val.xlarray.rows = xlRefInput.rwLast -
                        xlRefInput.rwFirst + 1;
                    output.val.xlarray.columns = xlRefInput.colLast -
                        xlRefInput.colFirst + 1;
                    int size = output.val.xlarray.rows *
                        output.val.xlarray.columns;
                    output.val.xlarray.lparray = size == 0? 0:
                        NEW_ARRAY(XL_OPER, size);
                    // work out offset into multi oper
                    int colOffset = xlRefInput.colFirst-xlRef.colFirst;
                    int rowOffset = xlRefInput.rwFirst - xlRef.rwFirst;
                    // then shallow copy over oper's
                    // and see if the opers have been 'calculated' ?
                    bool hitNoCalc = false;
                    BoolArray& recalc = staticRecalculated[i];
                    int width = multiOper->val.xlarray.columns;
                    for (int j = 0; j < output.val.xlarray.columns; j++){
                        for (int k = 0; k < output.val.xlarray.rows; k++){
                            XLConvert::offset(j, k, output) =
                                XLConvert::offset(colOffset+j, rowOffset+k,
                                                  *multiOper);
                            if (!recalc[(k+rowOffset)*width+(j+colOffset)]){
                                recalc[(k+rowOffset)*width+(j+colOffset)]=true;
                                hitNoCalc = true;
                            }
                        }
                    }
                    if (hitNoCalc){
                        xlFreeCoerce(output);
                        return xlretUncalced; // and escape
                    }
                    return xlretSuccess;
                }
            }
        }
        // something's either wrong with the tester or the addin kit
        throw ModelException("XLTest::coerceToMulti", "Couldn't find cached "
                             "value for reference oper");
    }

    /** Free the contents of the OPER created in XLCoerceToMulti method */
    static void xlFreeCoerce (XLOPER& toFreeContentsOf){
        FREE(toFreeContentsOf.val.xlarray.lparray);
    }
};

bool ExcelParams::staticCalledFromVB = false;
ExcelParams::EDROperArray ExcelParams::staticRefParams;
ExcelParams::EDROperArray ExcelParams::staticResolvedRefParams;
vector<BoolArray> ExcelParams::staticRecalculated;

CClassConstSP const ExcelParams::EDROper::TYPE =
CClass::registerClassLoadMethod(
    "ExcelParams::EDROper", typeid(ExcelParams::EDROper), load);

// work around for msvc 7
typedef ExcelParams::EDROperArray ExcelParamsEDROperArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("ExcelParams::EDROperArray", ExcelParamsEDROperArray);

CClassConstSP const ExcelParams::TYPE = CClass::registerClassLoadMethod(
    "ExcelParams", typeid(ExcelParams), load);

/* Run any addin functin where the name is specified by the first OPER -
   addin function name needs to be without the EDG_ prefix */
XL_OPER* XLTest::generic(
    const XL_OPER*  addinNameOper, /* (I) NB type 'P' */
    const XL_OPER*  fileNameOper,  /* (I) NB type 'P', may be NULL */
    const XL_OPER*  a0,            /* (I) first real parameter */
    va_list      args)          /* (I) remaining args */
{
    static const string routine  = "XLTest::generic";
    XL_OPER*               xlOutput = 0;

    try{
        bool        createFile;
        /* convert addinName into string */
        if (!(addinNameOper->type & xltypeStr)) {
            // this is probably more helpful that, say, converting a number
            // to a filename
            throw ModelException(routine,
                                 "First cell must contain addin name");
        }
        string addinName = XLConvert::operToString(*addinNameOper);
        /* if fileNameOper a string, then use it - otherwise
           don't do dump to file */
        string inpFileName;
        string outFileName;
        if (fileNameOper && fileNameOper->type & xltypeStr) {
            inpFileName = XLConvert::operToString(*fileNameOper);
        }
        createFile = !inpFileName.empty();
        /* find addin type */
        const XLAddin* xlAddin = XLAddin::lookUp(addinName);

        /* construct output file name */
        if (createFile) {
            outFileName = OutputFile::createOutputFileName(inpFileName);
        }
        /* and call main function */
        if (!(xlOutput = xlAddin->execute(createFile? &inpFileName: 0,
                                          createFile? &outFileName: 0,
                                          a0, args))){
            // shouldn't happen unless no memory
            throw ModelException(routine, "XLAddin::execute returned null"
                                 " oper!");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        e.errorLog();
    } catch (exception& e){
        ModelException e2(e, routine);
        e2.addMsg(routine);
        e2.errorLog();
    }

    if (!xlOutput && (xlOutput = NEW(XL_OPER))) {
        xlOutput->type = xltypeErr;
        xlOutput->val.err = xlerrValue;
        xlOutput->type |= xlbitDLLFree;
    }
    return xlOutput;
}

/* Run any addin functin where the name is specified by the first OPER -
   addin function name needs to be without the EDG_ prefix */
XL_OPER* XLTest::genericXL(
    bool         calledFromVB,  // (I) true if called from VB
    const XL_OPER*  addinNameOper, /* (I) NB type 'P' */
    const XL_OPER*  fileNameOper,  /* (I) NB type 'P' */
    va_list      args)          /* (I) remaining args */
{
    static const string routine  = "XLTest::genericXL";
    XL_OPER*               xlOutput = 0;

    try{
        bool        createFile;
        /* convert addinName into string */
        if (!(addinNameOper->type & xltypeStr)) {
            // this is probably more helpful that, say, converting a number
            // to a filename
            throw ModelException(routine,
                                 "First cell must contain addin name");
        }
        string addinName = XLConvert::operToString(*addinNameOper);
        /* if fileNameOper a string, then use it - otherwise
           don't do dump to file */
        string inpFileName;
        string outFileName;
        ErrorHandler origErrorHandler;
        vector<string> exceptionStack; // holds error messages

        if (fileNameOper && fileNameOper->type & xltypeStr) {
            inpFileName = XLConvert::operToString(*fileNameOper);
        }
        createFile = !inpFileName.empty();
        /* find addin type */
        const XLAddin* xlAddin = XLAddin::lookUp(addinName);
        // find number of parameters this function takes
        int numParams = xlAddin->numParams();

        // now start process of catching inputs and writing to file
        bool cachedRefsOK;
        ExcelParams xlParams(addinName, args, numParams,
                             calledFromVB, cachedRefsOK);
        if (!cachedRefsOK){
            return 0; // what MS tells us to do
        }
        /* serialise out */
        if (createFile){
            XMLWriter    xml(inpFileName);
            xlParams.write("XL-TEST", &xml);
            outFileName = OutputFile::createOutputFileName(inpFileName);
            // create and use new error handler & save old one.
            // (send errors to regression test output file)
            origErrorHandler = ErrorHandler::set(ErrorHandler(exceptionStack),
                                                 false);
        }

        /* and call main function */
        IObjectSP result;
        try{
            result = xlParams.runTest();
        } catch (exception&){
            if (createFile){
                 // restore error handler
                ErrorHandler::set(origErrorHandler, true);
            }
            throw;
        }
        // serialise out results
        if (createFile){
            // restore error handler
            ErrorHandler::set(origErrorHandler, true);
            OutputFile fileOut(outFileName);
            fileOut.write(result.get(), CStringArray(exceptionStack.begin(), exceptionStack.end()));
        }
        /* extract return outputs */
        ExcelParams::EDROperSP
            oper(ExcelParams::EDROperSP::dynamicCast(result));
        xlOutput = NEW(XL_OPER);
        ExcelParams::EDROper::operCopy(oper->getOper(), xlOutput);

    } catch (ModelException& e){
        e.addMsg(routine);
        e.errorLog();
    } catch (exception& e){
        ModelException e2(e, routine);
        e2.addMsg(routine);
        e2.errorLog();
    }

    if (!xlOutput && (xlOutput = NEW(XL_OPER))) {
        xlOutput->type = xltypeErr;
        xlOutput->val.err = xlerrValue;
        xlOutput->type |= xlbitDLLFree;
    }
    return xlOutput;
}

DRLIB_END_NAMESPACE
