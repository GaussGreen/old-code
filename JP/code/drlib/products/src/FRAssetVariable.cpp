//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRAssetVariable.cpp
//
//   Description : Variables Relating to Market Data for Flex Rules
//
//   Author      : Mark A Robson
//
//   Date        : 9 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRController.hpp"
#include "edginc/FRFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Used to 'read' the spot price of the asset in the simulation */
class FRAssetVariable: public FR::LValueBase,
                       public FRIfaces::IVarBarrierLevelAssist {
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns the input to define the variable in IMS 
        asset variables are define at the level of the basket definition, not in the list of variables */
    // assetVarNames needed for FlexInstrument interface
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Asset";
        result[1] = name;
        result[2] = "";

        return result;
    }

    /** Class for caching random number path. */
    class Path: public virtual FRController::INotifyPayoffCall{
        int                    assetIdx;
        int                    pathIdx;
        bool                   calcRefLevel;
        mutable const double*  pathArray;
        mutable double         refLevel;
    public:
        ~Path(){}
        // simple constructor
        Path(int                 assetIdx,
             int                 pathIdx,
             bool                calcRefLevel): 
            assetIdx(assetIdx), pathIdx(pathIdx), calcRefLevel(calcRefLevel){}
        /** invoked every time we price a new path */
        virtual void update(const FRController* ctrl){
            const MCPathGenerator* pathGen =
                ctrl->getPathGenerator();
            pathArray = pathGen->Path(assetIdx, pathIdx);
            if (calcRefLevel){
                refLevel = pathGen->refLevel(assetIdx, pathIdx);
            }
        }

        // returns the path as an array of doubles
        inline const double* path() const{
            return pathArray;
        }
        inline double perf(int simIdx) const{
            return (pathArray[simIdx]/refLevel);
        }

    };

    class Spot: public FR::RValueDouble, 
                public virtual FRIfaces::ILValueDouble{
    private:
    struct MyRT{
        TGetValue*         func;
        TSetValue*         setFunc;
        const Path*        path;  // reference to 'Path' class - defined above
        int                simIdx;   // which date we're on

        explicit MyRT(const Path*        path,
                      int                simIdx):
            func(&getValue), setFunc(&setValue), path(path), simIdx(simIdx){}

        static void setValue(void* structure, double value){
            throw ModelException("FRAssetVariable", "Cannot assign a"
                                 " value to a const variable");
        }
        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->path->path()[rt->simIdx]);
        }
        
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT* rt;
    public:
        // simple constructor
        Spot(const Path* path, int  simIdx):
            rt(new MyRT(path, simIdx)){}

        ~Spot(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

    class Perf: public FR::RValueDouble,
                public virtual FRIfaces::ILValueDouble{
    private:
    struct MyRT{
        TGetValue*         func;
        TSetValue*         setFunc;
        const Path*        path;  // reference to 'Path' class - defined above
        int                simIdx;   // which date we're on

        explicit MyRT(const Path*        path,
                      int                simIdx):
            func(&getValue), setFunc(&setValue), path(path), simIdx(simIdx){}

        static void setValue(void* structure, double value){
            throw ModelException("FRAssetVariable", "Cannot assign a"
                                 " value to a const variable");
        }
        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->path->perf(rt->simIdx));
        }
        
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT* rt;
    public:
        // simple constructor
        Perf(const Path* path, int  simIdx):
            rt(new MyRT(path, simIdx)){}

        ~Perf(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

    /** Class for caching random number path when using state variables */
    class PathSV: public virtual FRController::INotifyPayoffCall,
                  public virtual FRController::IStateVarClient{
        const SVPath*       spotPath; /* holds pathSV->path(assetIdx) for 
                                               current path */
        SVGenSpot::IStateVarSP       spotSV;
        IRefLevel::IStateVarSP    refLevelSV;
        int                       assetIdx;
        mutable double            refLevel;
        SVGenSpotSP                  spotGen;      //!< Generator for spot
        IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    public:
        ~PathSV(){}
        // simple constructor
        PathSV(SVGenSpotSP                  spotGen,      //!< Generator for spot
               IRefLevel::IStateVarGenSP refLevelGen,
               int                       assetIdx): 
            spotPath(0), assetIdx(assetIdx), 
            spotGen(spotGen), refLevelGen(refLevelGen) {}

        /** invoked every time we price a new path */
        virtual void update(const FRController* ctrl){
            spotPath = &spotSV->path(assetIdx);
            if (refLevelGen.get()){
                refLevel = refLevelSV->refLevel(assetIdx);
            }
        }

        /** Populates the collector with the state variables required by the
            various assets */
        virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
            svCollector->append(spotGen.get()); // ask for spot SV
            if (refLevelGen.get()){
                svCollector->append(refLevelGen.get()); // ask for ref Level SV
            }
        }

        /** To be called when the path generator changes (currently before doing
            the past, and before doing the future) */
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
            spotSV = spotGen->getSpotSV(newPathGen);
            if (refLevelGen.get()){
                refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            }
        }

        //// returns the path for specified simulation index
        inline double path(int simIdx) const{
            return ((*spotPath)[simIdx]);
        }
        inline double perf(int simIdx) const{
            return ((*spotPath)[simIdx]/refLevel);
        }
    };

    class SpotSV: public FR::RValueDouble, 
                  public virtual FRIfaces::ILValueDouble{
    private:
    struct MyRT{
        TGetValue*         func;
        TSetValue*         setFunc;
        const PathSV*      path;  // reference to 'Path' class - defined above
        int                simIdx;   // which date we're on

        explicit MyRT(const PathSV*      path,
                      int                simIdx):
            func(&getValue), setFunc(&setValue), path(path), simIdx(simIdx){}

        static void setValue(void* structure, double value){
            throw ModelException("FRAssetVariable", "Cannot assign a"
                                 " value to a const variable");
        }
        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->path->path(rt->simIdx));
        }
        
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT* rt;
    public:
        // simple constructor
        SpotSV(const PathSV* path, int  simIdx):
            rt(new MyRT(path, simIdx)){}

        ~SpotSV(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

    class PerfSV: public FR::RValueDouble,
                  public virtual FRIfaces::ILValueDouble{
    private:
    struct MyRT{
        TGetValue*         func;
        TSetValue*         setFunc;
        const PathSV*      path;  // reference to 'Path' class - defined above
        int                simIdx;   // which date we're on

        explicit MyRT(const PathSV*      path,
                      int                simIdx):
            func(&getValue), setFunc(&setValue), path(path), simIdx(simIdx){}

        static void setValue(void* structure, double value){
            throw ModelException("FRAssetVariable", "Cannot assign a"
                                 " value to a const variable");
        }
        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->path->perf(rt->simIdx));
        }
        
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT* rt;
    public:
        // simple constructor
        PerfSV(const PathSV* path, int  simIdx):
            rt(new MyRT(path, simIdx)){}

        ~PerfSV(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return name;
    }

    /** for FRIfaces::IVarBarrierLevelAssist */
    string getAssetName() const {
        return assetName;
    }
    double getRefLevel(FRController* frCtrl) const {
        double refLevel;
        if (frCtrl->stateVarUsed()){
            throw ModelException("FRAssetVariable::getRefLevel", "Not yet supported with SV");
        } else {
            const MCPathGenerator* pathGen =
                frCtrl->getPathGenerator();
            const FRIfaces::IProductView* productView = frCtrl->productView();
            int assetIdx = productView->findAssetIndex(assetName, ccyTreatment);
            // This is messy : findPathIndex() actually updates the instrument's
            // internal list of volRequests, which we do not wish to do here.
            // It is safe to simply use 0 since a check is in place to 
            // only use this after the last ref level date, and so there
            // is no longer any dependency on the vol.
#if 0
            // const problems: stems from interface to LN path generator
            CVolRequestSP req(CVolRequestSP::constCast(volRequest));
            // NB Need to call this method currently (even for SV) as this
            // passes the vol request to the MC product
            int pathIdx = productView->findPathIndex(assetIdx, req);
#else
            int pathIdx = 0;
#endif
            refLevel = pathGen->refLevel(assetIdx, pathIdx);
        }
        if (isPerf) {
            return refLevel;
        } else {
            return 1.0;
        }
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        // indentify our asset idx and path idx
        const FRIfaces::IProductView* productView = frCtrl->productView();
        int assetIdx = productView->findAssetIndex(assetName, ccyTreatment);
        // const problems: stems from interface to LN path generator
        CVolRequestSP req(CVolRequestSP::constCast(volRequest));
        // NB Need to call this method currently (even for SV) as this
        // passes the vol request to the MC product
        int pathIdx = productView->findPathIndex(assetIdx, req);
        int numDates = productView->getSimDates()->size();
        // create instance of Spot for all sim dates - there's no real need
        // to do this - just thought (incorrectly?) that the cost of looking
        // up the assetIdx and pathIdx would be costly if we kept on doing it
        // perhaps we should always do the 1st sim date and then look that up
        // to get the assetIdx and pathIdx.
        if (frCtrl->stateVarUsed()){
            if (pathIdx != 0){
                throw ModelException("FRAssetVariable::createRValue", "Only "
                                     "one vol request per asset supported");
            }
            IRefLevel::IStateVarGenSP refLevelGen;
            if (isPerf){
                refLevelGen = productView->getRefLevelGen();
            }
            SVGenSpotSP spotGen(productView->getSpotGen());
            // create common 'Path' class for sim time performance
            PathSV* path = new PathSV(spotGen, refLevelGen, assetIdx);
            frCtrl->addToNotifyPayoffCall(path, true); // avoid leak
            frCtrl->addToStateVars(path, false); // register with controller
            if (isPerf) {
                for (int i = 0; i < numDates; i++){
                    FRIfaces::IRValueSP rValue(new PerfSV(path, i));
                    frCtrl->setRValue(i, this, rValue);
                }
            } else {
                for (int i = 0; i < numDates; i++){
                    FRIfaces::IRValueSP rValue(new SpotSV(path, i));
                    frCtrl->setRValue(i, this, rValue);
                }
            }
        } else {
            // create common 'Path' class for sim time performance
            Path* path = new Path(assetIdx, pathIdx, isPerf);
            frCtrl->addToNotifyPayoffCall(path, true); // avoid leak
            if (isPerf) {
                for (int i = 0; i < numDates; i++){
                    FRIfaces::IRValueSP rValue(new Perf(path, i));
                    frCtrl->setRValue(i, this, rValue);
                }
            } else {
                for (int i = 0; i < numDates; i++){
                    FRIfaces::IRValueSP rValue(new Spot(path, i));
                    frCtrl->setRValue(i, this, rValue);
                }
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string               name;
    const string               assetName;
    const string               ccyTreatment; /* optional - may
                                                specialise assetName
                                                default is to use that
                                                given in MultiAsset */
    const CVolRequestConstSP   volRequest;   // optional
    const bool                 isPerf;       // optional, default false

    FRAssetVariable(): LValueBase(TYPE), isPerf(false) {}

    FRAssetVariable(const string& name,
                    const string& assetName,
                    bool          isPerf): 
        LValueBase(TYPE), name(name),
        assetName(assetName), isPerf(isPerf) {}

    static IObject* defaultFRAssetVariable(){
        return new FRAssetVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRAssetVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not appropriate here");
        }
        // Here "value" means assetName
        // Default to isPerf==true, but should allow to be false.
        // Possibly parse value to recognise "a1,n", for example?
        size_t n = strcspn(varName.c_str(), ",");
        size_t len = varName.length();
        bool   isPerf = true;
        string myVarName;
        if (n == len) {
            // no comma, so just treat directly
            myVarName = varName;
        } else {
            // interpret trailing "n"/"N"/"y"/"Y" as isPerf
            // and extract the first part as the actual var name
            if (n+2 != len) {
                throw ModelException(routine,
                                     "Allow format <varName,p> where p is n/N/y/Y but given " + varName);
            }
            if (varName[n+1]=='n' || varName[n+1]=='N') {
                isPerf = false;
            } else if (varName[n+1]=='y' || varName[n+1]=='Y') {
                isPerf = true;
            } else {
                throw ModelException(routine,
                                     "Allow format <varName,p> where p is n/N/y/Y but given " + varName);
            }
            myVarName = string(varName.c_str(), n);
        }
        return new FRAssetVariable(myVarName, value[0], isPerf);
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRAssetVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(assetName, "Asset Name");
        FIELD(ccyTreatment, "Currency Treatment");
        FIELD_MAKE_OPTIONAL(ccyTreatment);
        FIELD(volRequest,            "Vol Request");
        FIELD_MAKE_OPTIONAL(volRequest);
        FIELD(isPerf,         "true=measures perf not absolute level");
        FIELD_MAKE_OPTIONAL(isPerf);
        EMPTY_SHELL_METHOD(defaultFRAssetVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_ASSET",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI asset spot Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRAssetVariable::create, "Asset");
    }
    
};

CClassConstSP const FRAssetVariable::TYPE =
CClass::registerClassLoadMethod("FRAssetVariable", typeid(FRAssetVariable),
                                load);

bool loadFRMarketVariables() {
    return (FRAssetVariable::TYPE != 0);
}

DRLIB_END_NAMESPACE

