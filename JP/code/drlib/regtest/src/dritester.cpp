//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : dritester.cpp
//
//   Description : regression tester for dr interface
//
//   Author      : Andrew J Swain (from models.cpp)
//
//   Date        : Autumn 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4786)
#endif
#include "edginc/EDGServices.h"
#include "edginc/DataDictionary.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Null.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Addin.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/Library.hpp"
#include "edginc/DRUtil.hpp"

USING_DRLIB_NAMESPACE

// Invoke only if (x) is a function that returns DRError or is a DRError type.
#define CHECK(x, userMsg)                                       \
    DRI_CHECK(x, method.c_str(), userMsg, handleError, 0, ;, ;, ;)

static void handleError(const char* method, const char *err, void *cbParam) {
    ModelException e(method, err);
    throw e;
}

static void toDRValue(DRService* svc, IObjectSP object, DRValue* val);

// turn EDR objects into DRObjects
static DRObject toDRObject(DRService* svc, IObjectSP object) {
    static const string method = "toDRObject";
    DRMap    map = 0;

    try {
        DRObject dro = 0;
        
        // 2 routes to test here - for regular classes and "exported" ones
        // where the interface is fixed for all time and you can assume
        // the order of arguments
        if (Modifier::isExport(object->getClass()->getModifiers())) {
            // if class is exported, then interface is fixed, so stupid 
            // to use proxies here
            
            // turn object into dictionary version
            CDataDictionarySP dd(CDataDictionary::pop2DataDict(object));
            CClassConstSP clazz = dd->getType();
            const char* typeName = clazz->getName().c_str();

            CFieldArray fields(Addin::getDataClassFields(clazz));
            int numFields = fields.size();
            int i;
            DRValue* params = new DRValue[numFields];
            for (i = 0; i < numFields; i++) {
                IObjectSP elem(dd->get(fields[i]->getName()));
                toDRValue(svc, elem, &params[i]);
            }

            CHECK(svc->fptr->mapNewFromValues(svc, typeName, params, numFields,
                                              &map), 
                  "map failed");
            for (i = 0; i < numFields; i++) {
                CHECK(svc->fptr->valueClear(svc, &params[i]), 0);
            }            
            delete[] params;
        }
        else {
            // here we might want to use a proxy (in interface terms) for the
            // class we want to build
            // typeName is what we want to build

            // bit broken right now - if we have a public & private interface
            // (NOT proxy) then getting the type from the object gives the 
            // private type, but from the equivalent dictionary gives the 
            // public type
            IObjectSP objPublic(CObject::convertToPublicRep(object));

            // make a dictionary of (any) proxy to use its interface
            IObjectSP o(CObject::toProxy(object));

            // if this object is a "true map" then as far as DR interface is
            // concerned, it's a map
            IMap* omap = dynamic_cast<IMap*>(o.get());    
            CDataDictionarySP dd;
            IMap::IIteratorSP iter;
            if (omap && omap->isTrueMap()) {
                iter = IMap::IIteratorSP(omap->createIterator());
            }
            else {
                dd = CDataDictionarySP(CDataDictionary::pop2DataDict(o));
                iter = IMap::IIteratorSP(dd->createIterator());
            }
                
            const char* typeName = objPublic->getClass()->getName().c_str();

            CHECK(svc->fptr->mapNew(svc, typeName, &map), "map failed");

            while (iter->hasMoreElements()) {
                DRValue obj;
                const char* field = iter->getKey().c_str(); 
                IObjectSP   elem(iter->getElement());

                toDRValue(svc, elem, &obj);
               
                CHECK(svc->fptr->mapAddItem(svc, map, field, &obj),
                      "map add failed");

                CHECK(svc->fptr->valueClear(svc, &obj), "valueClear failed");

                iter->increment();
            }
            
        }
        CHECK(svc->fptr->mapToObject(svc, map, &dro), "mapToObject failed");

        if (map) {
            CHECK(svc->fptr->objectFree(svc, map), "objectFree failed");
        }

        map = 0;
        return dro;
    } 
    catch (ModelException& e){
        if (map) {
            CHECK(svc->fptr->objectFree(svc, map), "objectFree failed");
        }
        throw ModelException(
            e, method, "Failed to rebuild " + object->getClass()->getName());
    }
}

// convert EDR matrix to DRI version
static DRMatrix toDRMatrix(DRService* svc, IObjectSP object) {
    static const string method = "toDRMatrix";

    try {
        DRMatrix matrix = 0;
        CDoubleMatrixSP mtx(CDoubleMatrixSP::dynamicCast(object));

        int rows = mtx->numRows();
        int cols = mtx->numCols();

        CHECK(svc->fptr->matrixNew(svc, cols, rows, &matrix),
              "matrixNew failed");

        for (int i = 0; i < cols; i++){
            for (int j = 0; j < rows; j++){
                CHECK(svc->fptr->matrixSet(svc, matrix, i, j, (*mtx)[i][j]),
                      "arraySet failed");
            }
        }   
        return matrix;
    } 
    catch (ModelException& e){
        throw ModelException(e, method);
    }  
}

// convert EDR arrays to DRI ones
static DRArray toDRArray(DRService* svc, IObjectSP object) {
    static const string method = "toDRArray";

    try {
        DRArray dra = 0;
        IArraySP array(IArraySP::dynamicCast(object));

        int    length = array->getLength();
        string className = array->getClass()->getComponentType()->getName();

        CHECK(svc->fptr->arrayNew(svc, length, className.c_str(), &dra),
              "arrayNew failed");
 
        for (int i = 0; i < length; i++){
            DRValue obj;
            CHECK(EDGSDRValue(&obj, DR_UNDEFINED), 0);

            toDRValue(svc, array->get(i), &obj);

            CHECK(svc->fptr->arraySet(svc, dra, i, &obj),
                  "arraySet failed");

            CHECK(svc->fptr->valueClear(svc, &obj), 0);
        }

        return dra;
    } 
    catch (ModelException& e){
        throw ModelException(e, method);
    }  
}

// turn primitive types into DRValues
static void toDRValue(DRService* svc, IObjectSP object, DRValue* val) {
    static const string method = "toDRValue";
    try {
        if (object->getClass()->isPrimitive()) {
            if (object->getClass()->isAssignableFrom(CBool::TYPE)) {
                CBool* bp = dynamic_cast<CBool*>(object.get());
                CHECK(EDGSDRValue(val, DR_BOOL, bp->boolValue()), 0);
            }
            else if (object->getClass()->isAssignableFrom(CInt::TYPE)) {
                CInt* ip = dynamic_cast<CInt*>(object.get());
                CHECK(EDGSDRValue(val, DR_INT, ip->intValue()), 0);
            }
            else if (object->getClass()->isAssignableFrom(CDouble::TYPE)) {
                CDouble* dp = dynamic_cast<CDouble*>(object.get());
                CHECK(EDGSDRValue(val, DR_DOUBLE, dp->doubleValue()), 0);
            }
            else if (object->getClass()->isAssignableFrom(CString::TYPE)) {
                CString* sp = dynamic_cast<CString*>(object.get());
                CHECK(EDGSDRValue(val, DR_STRING, sp->stringValue().c_str()), 
                      0);
            }
            else {
                CHECK(EDGSDRValue(val, DR_UNDEFINED), 0);
            }
        }
        else {
            DRObject dro = 0;

            if (object->getClass()->isArray()){
                dro = toDRArray(svc, object); //IObjectSP::attachToRef(object));
            }
            else if (CDoubleMatrix::TYPE->isInstance(object)){
                dro = toDRMatrix(svc, object);
            }
            else if (!CNull::TYPE->isInstance(object)){
                dro = toDRObject(svc, object); //IObjectSP::attachToRef(object));
            }

            CHECK(EDGSDRValue(val, DR_OBJECT, dro), 0);
        }
    } 
    catch (ModelException& e){
        throw ModelException(e, method);
    }  
}

// pull inputs to bits and rebuild via Global DR Interface
// then run it and see what happens
REGTEST_DLL void runViaDRI(IObjectSP input, const string& outfile) {
    static const string method = "runViaDRI";
    
    // fire up DRService
    DRServiceInitArgs args;
    EDGSInitServiceCreateArgs(&args);

    DRService *svc = 0;

    if (!DRCreateService(&args, &svc)) {
        printf("Failed to create EDR Service!\n");
        const char *err = 0;
        EDGSGetServiceInvocationError(&err);
        throw ModelException(method, err);
    }
 
    DRObject   dro = 0;
    DRValue    output;

    try {
        // gotta have that dot
        cout << "." << flush;

        CHECK(EDGSDRValue(&output, DR_UNDEFINED), 0);

        dro = toDRObject(svc, input);
        
        CHECK(svc->fptr->execute(svc, dro, &output), "failed to execute");
            
        if (output.type != DR_OBJECT) {
            throw ModelException(method, "result is not an object");
        }

        CHECK(EDGSWriteResultsToFile(svc, outfile.c_str(), 
                                     output.value.object),
              "failed to write results");
    } 
    catch (ModelException& e){
        ModelException error(method, e.stackTrace());
        OutputFile out(outfile);
        out.write(error);
    }
    if (svc){
        DRI_FREE(svc, dro);
        CHECK(svc->fptr->valueClear(svc, &output), 0);
        CHECK(svc->fptr->serviceFree(svc), 0);
    }
    // bring the library back to life ...
    Library::startup();
}
