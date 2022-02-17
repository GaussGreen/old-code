//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRBarrierVariable.cpp
//
//   Description : Barrier Variables Flex Rules. 
//                 Allows semantic extension of Bool to report events.
//
//   Date        : July 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FR.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Format.hpp"
#include "edginc/FRController.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/Events.hpp"

DRLIB_BEGIN_NAMESPACE

class FRBarrierVariable: public FR::LValueBase {

    class IBarrierEvent {
    public:
        virtual bool getValue() = 0;
        virtual FlexBarrierBreach* createEvent(const string& monVarName,
                                               const DateTime& barDate) = 0; 
        virtual ~IBarrierEvent(){};
    };

    // Works for both known and other cases via IBarrierEvent
    class BarrierEventer : public virtual FRController::IFREvent {
    private:
        IBarrierEvent*      var; // const ...?
        const string&       monVarName;
        const DateTime&     barDate;
        
    public:
        BarrierEventer(IBarrierEvent*      var,
                       const string&       monVarName,
                       const DateTime&     barDate):
            var(var), monVarName(monVarName), barDate(barDate) {}

        virtual bool hasEvent() {
            return var->getValue();
        };
        virtual FlexBarrierBreach* createEvent() { 
            return var->createEvent(monVarName, barDate);  
        };
        
    };

    // FRIfaces::ILValueBool and not IRValueBool so that 
    // instances of this type can be placed into BoolArrays
    class Barrier : public FR::RValueBool,
                    public virtual FRIfaces::ILValueBool,
                    public virtual IBarrierEvent {
    private:
        struct MyRT{
            TGetValue*                    func;
            TSetValue*                    setFunc;
            FRIfaces::IRValueDouble::RT*  monVar;
            bool                          isUp;
            FRIfaces::IRValueDouble::RT*  levelVar;

            explicit MyRT(FRIfaces::ILValueDouble* monVar,
                          bool                     isUp,
                          FRIfaces::ILValueDouble* levelVar): 
                func(&getValue),
                setFunc(&setValue),
                monVar(monVar->getRT()), isUp(isUp),
                levelVar(levelVar->getRT()) {}

            static void setValue(void* structure, bool value){
                throw ModelException("FRBarrierVariable", "Cannot assign a"
                                     " value to a read-only variable");
            }
            static bool getValue(void* structure) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVar->func(rt->monVar) >
                            rt->levelVar->func(rt->levelVar));
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        FRIfaces::IRValueDouble*  monVar;
        FRIfaces::IRValueDouble*  levelVar;
        MyRT* rt;
    public:
        Barrier(FRIfaces::ILValueDouble*  monVar,
                bool                      isUp,
                FRIfaces::ILValueDouble*  levelVar):
            monVar(monVar), levelVar(levelVar),
            rt(new MyRT(monVar, isUp, levelVar)) {}

        ~Barrier() {
            delete rt;
        }

        virtual bool isKnown() const {
            return (monVar->isKnown() && levelVar->isKnown());
        }

        virtual bool getValue() {
            return MyRT::getValue(rt);
        }

        virtual void setValue(bool value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} 
        
        virtual FRIfaces::IRValueBool::RT* getRT() {
            return (FRIfaces::IRValueBool::RT*)rt;
        }
        virtual FlexBarrierBreach* createEvent(const string& monVarName,
                                               const DateTime& barDate) { 
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->levelVar->func(rt->levelVar), 
                                         rt->monVar->func(rt->monVar), 
                                         rt->isUp);  
        };
    };

    class BarrierKnownLevel : public FR::RValueBool,
                              public virtual FRIfaces::ILValueBool,
                              public virtual IBarrierEvent  {
    private:
        struct MyRT{
            TGetValue*                    func;
            TSetValue*                    setFunc;
            FRIfaces::IRValueDouble::RT*  monVar;
            bool                          isUp;
            double                        level;

            explicit MyRT(FRIfaces::ILValueDouble* monVar,
                          bool                     isUp,
                          double                   level): 
                func(&getValue),
                setFunc(&setValue),
                monVar(monVar->getRT()), isUp(isUp), level(level) {}
            

            static void setValue(void* structure, bool value){
                throw ModelException("FRBarrierVariable", "Cannot assign a"
                                     " value to a read-only variable");
            }
            static bool getValue(void* structure) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVar->func(rt->monVar) > rt->level);
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        FRIfaces::IRValueDouble*  monVar;
        MyRT* rt;
    public:
        BarrierKnownLevel(FRIfaces::ILValueDouble*  monVar,
                          bool                      isUp,
                          double                    level) :
            monVar(monVar), rt(new MyRT(monVar, isUp, level)) {}

        ~BarrierKnownLevel() {
            delete rt;
        }

        virtual bool isKnown() const {
            return monVar->isKnown();
        }

        virtual bool getValue() {
            return MyRT::getValue(rt);
        }

        virtual void setValue(bool value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} 

        virtual FRIfaces::IRValueBool::RT* getRT() {
            return (FRIfaces::IRValueBool::RT*)rt;
        }
        virtual FlexBarrierBreach* createEvent(const string& monVarName,
                                               const DateTime& barDate) {
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->level, 
                                         rt->monVar->func(rt->monVar),
                                         rt->isUp);  
        };
    };

public:
    static CClassConstSP const TYPE;

    static const string UP;
    static const string DOWN;

    ~FRBarrierVariable(){}

    virtual void validatePop2Object() {
        isUp = CString::equalsIgnoreCase(barrierType, UP) ||
                (barrierType.length()==1 &&
                 CString::equalsIgnoreCase(barrierType, UP, 1)); // allow 'U'

        bool isDown = CString::equalsIgnoreCase(barrierType, DOWN) ||
                (barrierType.length()==1 &&
                 CString::equalsIgnoreCase(barrierType, DOWN, 1)); // allow 'D'
        
        if (!isUp && !isDown) {
            throw ModelException("FRBarrierVariable::validatePop2Object", 
                                 "Barrier type must indicate " + UP + " or " + DOWN +
                                 " but given " + barrierType);
        }
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Barrier";
        result[1] = name;
        result[2] = monitorVarName+","+(isUp?"Up":"Down")+","+ level;

        return result;
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRBarrierVariable::getLValue", "Variable "+
                             name +" is read-only");
    }

    FRBarrierVariable(): LValueBase(TYPE),
                         isUp(false)  {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint.
    The monitoring is implicitly performed only at time points where this
    variable is referenced. */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRBarrierVariable::createRValue";
        FRIfaces::IRValueBool* bar = 0;
        double barLevel = 0.;
        bool   haveBarLevel = false;

        if (!frCtrl->getRValue(this, index)){
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            FRIfaces::ILValueExpression* lValueExp = 
                frCtrl->getLValueExpression(monitorVarName);
            if (!lValueExp){
                string m("No variable with name '"+monitorVarName+"'\n"
                         "When defining '"+name+"' you must supply the "
                         "name of the variable whose level is monitored. ");
                throw ModelException(method, m);
            }
            FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                             frCtrl);
            // cast to double
            FRIfaces::ILValueDouble* monVar = 
                dynamic_cast<FRIfaces::ILValueDouble*>(lValue);
            if (!monVar){
                throw ModelException(method, "Monitoring variable is '"+
                                     monitorVarName+"' which must be a "
                                     "double variable");
            }
            // level 
            FRIfaces::ILValueExpression* lValueExpLevel = 
                frCtrl->getLValueExpression(level);
            if (lValueExpLevel){
                // it's a variable name
                FRIfaces::IRValue* rValueLevel = 
                    lValueExpLevel->getRValue(index, frCtrl);
                // cast to double
                FRIfaces::IRValueDouble* levelVar = 
                    dynamic_cast<FRIfaces::IRValueDouble*>(rValueLevel);
                if (!levelVar){
                    throw ModelException(method, "Barrier level variable is '"+
                                         level+"' which must be a "
                                         "double variable");
                }
                if (levelVar->isKnown()) {
                    barLevel = levelVar->getValue();
                    haveBarLevel = true;
                    bar = new BarrierKnownLevel(monVar, isUp, barLevel);
                } else {
                    FRIfaces::ILValueDouble* lValueLevel = 
                        dynamic_cast<FRIfaces::ILValueDouble*>(lValueExpLevel->getLValue(index, frCtrl));

                    bar = new Barrier(monVar, isUp, lValueLevel);
                }
            } else {
                // not a variable name - see if it can be parsed as a number
                char*   endPos;
                double  dVal = strtod(level.c_str(), &endPos);
                if (*endPos == '\0') {
                    bar = new BarrierKnownLevel(monVar, isUp, dVal);
                } else {
                    string m("No variable with name '"+level+"'\n"
                             "When defining '"+name+"' the level must be given "
                             "as the name of a variable or directly as a value. ");
                    throw ModelException(method, m);
                }
                barLevel = dVal;
                haveBarLevel = true;
            }
            // keep track of the event if we need to 
            if (index <= frCtrl->getEventIndex()) {
                IBarrierEvent* barEv = dynamic_cast<IBarrierEvent*>(bar);
                frCtrl->addEvent(new BarrierEventer(barEv, monitorVarName,
                                                    frCtrl->getDate(index)));
            }
            // Consider reporting a BARRIER_LEVEL
            const DateTime& today = frCtrl->productView()->getValueDate();
            DateTime upperDate = BarrierLevel::barrierWindow(today);
            const DateTime& thisDate = frCtrl->getDate(index);
            if (thisDate>=today &&
                thisDate<=upperDate &&
                haveBarLevel) {
                // check mon var is a simple asset
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression(monitorVarName);
                if (lValueExp) {
                    FRIfaces::IVarBarrierLevelAssist* barLevelAssist = 
                        dynamic_cast<FRIfaces::IVarBarrierLevelAssist*>(lValueExp);
                    if (barLevelAssist) {
                        string assetName = barLevelAssist->getAssetName();
                        frCtrl->addBarrierLevel(assetName,
                                                isUp,
                                                barLevel,
                                                thisDate,
                                                barLevelAssist);
                    }
                }
            }
        }
        return bar;
    }

private:
    // field
    const string name;
    const string monitorVarName; // which var is checked against level
    const string barrierType;    // "Up" or "Down". String rather than bool allows easier extension
                                 // and greater clarity in the interface
    const string level;          // var name or double value

    // transient
    bool isUp;                   // translated from barrierType

    FRBarrierVariable(const string varName,
                      const string monitorVarName,
                      const string barrierType,
                      const string level): 
        LValueBase(TYPE),
        name(varName), monitorVarName(monitorVarName),
        barrierType(barrierType), level(level) {
            validatePop2Object();
    }

    static IObject* defaultFRBarrierVariable(){
        return new FRBarrierVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRBarrierVariable::create";
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        if (value[0].empty() || value[0] == "-" ) {
            throw ModelException(routine, "Value should contain monitorVarName, "
                                 "U/D and level");
        }
        
        // Format of 'value' is these params in this order, comma separated
        string monitorVarName; // var name we check against 'level'
        string barrierType; // up/down
        string levelStr; // var name or double value that defines 'level' 

        const char* pos = value[0].c_str();
        int paramId = 0;
        while (*pos) {
            char   c;
            // skip whitespace
            while ((c = *pos) == ' ' || c == '\t'){
                pos++;
            }
            if (c == ','){
                // skip over our chosen separator and continue
                pos++;
                // count params
                paramId++;
                if (paramId>2) {
                    throw ModelException(routine, 
                                         "Too many commas (must have 2) with "
                                         "barrier variable " + varName + 
                                         " and value of " + value[0]);
                }
            } else {
                // next contiguous block is the text we want
                const char* varBegin = pos;
                do{
                    pos++; /* Get another character. */
                    c = *pos;
                } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                const char* varEnd = pos;
                // done with this one
                switch (paramId) {
                case 0:
                    monitorVarName = string(varBegin, varEnd - varBegin);
                    break;
                case 1:
                    barrierType = string(varBegin, varEnd - varBegin);
                    break;
                case 2:
                    // Perhaps offer that this could be interpreted as
                    // a double directly?
                    levelStr = string(varBegin, varEnd - varBegin);
                    break;
                }
            } 
        }

        if (paramId != 2) {
            throw ModelException(routine, "For barrier variable " + varName + 
                                 " require name of monitored var, U/D and "
                                 "level in that order and comma "
                                 "separated but given " + value[0]);
        }

        // A simple "return new ..." gives gcc compilation warning 
        FR::LValueBase* bar = new FRBarrierVariable(varName,
                                                    monitorVarName,
                                                    barrierType,
                                                    levelStr);
        return bar;
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBarrierVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRBarrierVariable);
        FIELD(name, "name of variable");
        FIELD(monitorVarName, "Name of variable to monitor versus "
                     "level");
        FIELD(barrierType, "Up/U or Down/D");
        FIELD(level, "Name of variable, or double value, that defines "
                     "monitoring level");
        FIELD_NO_DESC(isUp);
        FIELD_MAKE_TRANSIENT(isUp);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_BARRIER",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Barrier Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRBarrierVariable::create, 
                                           "Barrier");
    }
    
};

CClassConstSP const FRBarrierVariable::TYPE = CClass::registerClassLoadMethod(
    "FRBarrierVariable", typeid(FRBarrierVariable), FRBarrierVariable::load);

const string FRBarrierVariable::UP = "Up";
const string FRBarrierVariable::DOWN = "Down";

bool  FRBarrierVariableLoad() {
    return ( FRBarrierVariable::TYPE != 0);
}
bool  loadFRBarrierVariable() {
    return ( FRBarrierVariable::TYPE != 0);
}

/************************************************************/
/* A convenience which has a single up/down and level, but 
   multiple monitoring vars. This saves creation of many 
   BarrierVars which are then put into an array. */
/************************************************************/
class FRBarrierArrayVariable: public FR::LValueBase {

    class IBarrierArrayEvent {
    public:
        virtual bool getValue(unsigned int index) = 0;
        virtual FlexBarrierBreach* createEvent(const string& monVarName,
                                               const DateTime& barDate,
                                               unsigned int  index) = 0; 
        virtual ~IBarrierArrayEvent(){};
    };

    // Works for both known and other cases via IBarrierArrayEvent
    class BarrierArrayEventer : public virtual FRController::IFREvent {
    private:
        IBarrierArrayEvent* var; // const ...?
        const string&       monVarName;
        const DateTime&     barDate;
        unsigned int        index;
        
    public:
        BarrierArrayEventer(IBarrierArrayEvent* var,
                            const string&       monVarName,
                            const DateTime&     barDate,
                            unsigned int        index):
            var(var), monVarName(monVarName), barDate(barDate), index(index) {}

        virtual bool hasEvent() {
            return var->getValue(index);
        };
        virtual FlexBarrierBreach* createEvent() { 
            return var->createEvent(monVarName, barDate, index);  
        };
        
    };

    class BarrierArray : public FR::RValueBoolArray,
                         public virtual FRIfaces::IRValueBoolArray,
                         public virtual IBarrierArrayEvent {
    private:
        struct MyRT : public FRIfaces::IRValueBoolArray::RT{
            vector<FRIfaces::IRValueDouble::RT*>  monVars;
            bool                                  isUp;
            FRIfaces::IRValueDouble::RT*          levelVar;

            explicit MyRT(vector<FRIfaces::ILValueDouble*> monVars,
                          bool                             isUp,
                          FRIfaces::ILValueDouble*         levelVar): 
                monVars(monVars.size()), isUp(isUp),
                levelVar(levelVar->getRT()) {
                size = &getSize;
                func = &getValue;
                for(unsigned int i=0; i<monVars.size(); i++) {
                    this->monVars[i] = monVars[i]->getRT();
                }
            }
            
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->monVars.size();
            }

            static bool getValue(void* structure, int index) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVars[index]->func(rt->monVars[index]) >
                            rt->levelVar->func(rt->levelVar));
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };

        vector<FRIfaces::ILValueDouble*>  monVars;
        FRIfaces::IRValueDouble*          levelVar;
        MyRT*                             rt;
    public:
        BarrierArray(vector<FRIfaces::ILValueDouble*>  monVars,
                     bool                              isUp,
                     FRIfaces::ILValueDouble*          levelVar):
            monVars(monVars), levelVar(levelVar),
            rt(new MyRT(monVars, isUp, levelVar)) {}

        ~BarrierArray() {
            delete rt;
        }

        virtual bool isKnown() const {
            bool known = levelVar->isKnown();
            for(unsigned int i=0; !known && i<monVars.size(); i++) {
                known = monVars[i]->isKnown();
            }
            return known;
        }

        virtual IObjectConstSP get(){
            CBoolArraySP bools(new BoolArray(rt->monVars.size()));
            for (unsigned int i = 0; i < rt->monVars.size(); i++){
                (*bools)[i] = rt->getValue(rt, i);
            }
            return IObjectSP(bools);
        }

        virtual void setReset(char* reset){} 
        
        virtual FRIfaces::IRValueBoolArray::RT* getRT() {
            return (FRIfaces::IRValueBoolArray::RT*)rt;
        }

        // IBarrierArrayEvent
        bool getValue(unsigned int index) {
            return rt->getValue(rt, index);
        }

        FlexBarrierBreach* createEvent(const string& monVarName,
                                       const DateTime& barDate,
                                       unsigned int  index) { 
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->levelVar->func(rt->levelVar),
                                         rt->monVars[index]->func(rt->monVars[index]),
                                         rt->isUp);  
        };
        
    };

    class BarrierArrayNative : public FR::RValueBoolArray,
                               public virtual FRIfaces::IRValueBoolArray,
                               public virtual IBarrierArrayEvent {
    private:
        struct MyRT : public FRIfaces::IRValueBoolArray::RT{
            FRIfaces::IRValueDoubleArray::RT*     monVars;
            bool                                  isUp;
            FRIfaces::IRValueDouble::RT*          levelVar;

            explicit MyRT(FRIfaces::ILValueDoubleArray* monVars,
                          bool                          isUp,
                          FRIfaces::ILValueDouble*      levelVar): 
                monVars(monVars->getRT()), isUp(isUp),
                levelVar(levelVar->getRT()) {
                size = &getSize;
                func = &getValue;
            }
            
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->monVars->size(rt->monVars);
            }

            static bool getValue(void* structure, int index) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVars->func(rt->monVars,index) >
                            rt->levelVar->func(rt->levelVar));
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };

        FRIfaces::ILValueDoubleArray*  monVars;
        FRIfaces::IRValueDouble*       levelVar;
        MyRT*                          rt;
    public:
        BarrierArrayNative(FRIfaces::ILValueDoubleArray*  monVars,
                           bool                           isUp,
                           FRIfaces::ILValueDouble*       levelVar):
            monVars(monVars), levelVar(levelVar),
            rt(new MyRT(monVars, isUp, levelVar)) {}

        ~BarrierArrayNative() {
            delete rt;
        }

        virtual bool isKnown() const {
            bool known = levelVar->isKnown() &&
                monVars->isKnown();
            return known;
        }

        virtual IObjectConstSP get(){
            int sz = rt->monVars->size(rt->monVars);
            CBoolArraySP bools(new BoolArray(sz));
            for (int i = 0; i < sz; i++){
                (*bools)[i] = rt->getValue(rt, i);
            }
            return IObjectSP(bools);
        }

        virtual void setReset(char* reset){} 
        
        virtual FRIfaces::IRValueBoolArray::RT* getRT() {
            return (FRIfaces::IRValueBoolArray::RT*)rt;
        }

        // IBarrierArrayEvent
        bool getValue(unsigned int index) {
            return rt->getValue(rt, index);
        }

        FlexBarrierBreach* createEvent(const string& monVarName,
                                       const DateTime& barDate,
                                       unsigned int  index) { 
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->levelVar->func(rt->levelVar),
                                         rt->monVars->func(rt->monVars,index),
                                         rt->isUp);  
        };
        
    };

    class BarrierArrayKnownLevel : public FR::RValueBoolArray,
                                   public virtual FRIfaces::IRValueBoolArray,
                                   public virtual IBarrierArrayEvent  {
    private:
        struct MyRT : public FRIfaces::IRValueBoolArray::RT{
            vector<FRIfaces::IRValueDouble::RT*>  monVars;
            bool                                  isUp;
            double                                level;

            explicit MyRT(vector<FRIfaces::ILValueDouble*> monVars,
                          bool                             isUp,
                          double                           level): 
                monVars(monVars.size()), isUp(isUp), level(level) {
                size = &getSize;
                func = &getValue;
                for(unsigned int i=0; i<monVars.size(); i++) {
                    this->monVars[i] = monVars[i]->getRT();
                }
            }
            
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->monVars.size();
            }

            static bool getValue(void* structure, int index) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVars[index]->func(rt->monVars[index]) > rt->level);
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        vector<FRIfaces::ILValueDouble*>  monVars;
        MyRT*                             rt;
    public:
        BarrierArrayKnownLevel(vector<FRIfaces::ILValueDouble*>  monVars,
                               bool                              isUp,
                               double                            level):
            monVars(monVars), rt(new MyRT(monVars, isUp, level)) {}

        ~BarrierArrayKnownLevel() {
            delete rt;
        }

        virtual bool isKnown() const {
            bool known = false;
            for(unsigned int i=0; !known && i<monVars.size(); i++) {
                known = monVars[i]->isKnown();
            }
            return known;
        }

        virtual IObjectConstSP get(){
            CBoolArraySP bools(new BoolArray(rt->monVars.size()));
            for (unsigned int i = 0; i < rt->monVars.size(); i++){
                (*bools)[i] = rt->getValue(rt, i);
            }
            return IObjectSP(bools);
        }

        virtual void setReset(char* reset){} 

        virtual FRIfaces::IRValueBoolArray::RT* getRT() {
            return (FRIfaces::IRValueBoolArray::RT*)rt;
        }

        // IBarrierArrayEvent
        bool getValue(unsigned int index) {
            return rt->getValue(rt, index);
        }

        FlexBarrierBreach* createEvent(const string& monVarName,
                                       const DateTime& barDate,
                                       unsigned int  index) { 
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->level,
                                         rt->monVars[index]->func(rt->monVars[index]),
                                         rt->isUp);  
        };

    };

    class BarrierArrayNativeKnownLevel : public FR::RValueBoolArray,
                                         public virtual FRIfaces::IRValueBoolArray,
                                         public virtual IBarrierArrayEvent  {
    private:
        struct MyRT : public FRIfaces::IRValueBoolArray::RT{
            FRIfaces::IRValueDoubleArray::RT*  monVars;
            bool                               isUp;
            double                             level;

            explicit MyRT(FRIfaces::ILValueDoubleArray* monVars,
                          bool                          isUp,
                          double                        level): 
                monVars(monVars->getRT()), isUp(isUp), level(level) {
                size = &getSize;
                func = &getValue;
            }
            
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->monVars->size(rt->monVars);
            }

            static bool getValue(void* structure, int index) {
                MyRT* rt = (MyRT*)structure;
                bool ret = (rt->monVars->func(rt->monVars,index) > rt->level);
                return (rt->isUp? ret : !ret);
            }

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        FRIfaces::ILValueDoubleArray*  monVars;
        MyRT*                          rt;
    public:
        BarrierArrayNativeKnownLevel(FRIfaces::ILValueDoubleArray*  monVars,
                                     bool                           isUp,
                                     double                         level):
            monVars(monVars), rt(new MyRT(monVars, isUp, level)) {}

        ~BarrierArrayNativeKnownLevel() {
            delete rt;
        }

        virtual bool isKnown() const {
            return monVars->isKnown();
        }

        virtual IObjectConstSP get(){
            int sz = rt->monVars->size(rt->monVars);
            CBoolArraySP bools(new BoolArray(sz));
            for (int i = 0; i < sz; i++){
                (*bools)[i] = rt->getValue(rt, i);
            }
            return IObjectSP(bools);
        }

        virtual void setReset(char* reset){} 

        virtual FRIfaces::IRValueBoolArray::RT* getRT() {
            return (FRIfaces::IRValueBoolArray::RT*)rt;
        }

        // IBarrierArrayEvent
        bool getValue(unsigned int index) {
            return rt->getValue(rt, index);
        }

        FlexBarrierBreach* createEvent(const string& monVarName,
                                       const DateTime& barDate,
                                       unsigned int  index) { 
            return new FlexBarrierBreach(barDate, monVarName, 
                                         rt->level,
                                         rt->monVars->func(rt->monVars,index),
                                         rt->isUp);  
        };

    };

public:
    static CClassConstSP const TYPE;

    static const string UP;
    static const string DOWN;

    ~FRBarrierArrayVariable(){}

    virtual void validatePop2Object() {
        static const string routine("FRBarrierArrayVariable::validatePop2Object");
        isUp = CString::equalsIgnoreCase(barrierType, UP) ||
                (barrierType.length()==1 &&
                 CString::equalsIgnoreCase(barrierType, UP, 1)); // allow 'U'

        bool isDown = CString::equalsIgnoreCase(barrierType, DOWN) ||
                (barrierType.length()==1 &&
                 CString::equalsIgnoreCase(barrierType, DOWN, 1)); // allow 'D'
        
        if (!isUp && !isDown) {
            throw ModelException(routine, 
                                 "BarrierArray type must indicate " + UP + " or " + DOWN +
                                 " but given " + barrierType);
        }
        if (level.empty()) {
            throw ModelException(routine, 
                                 "No level supplied!");
        }
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolArrayType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);

        if (monitorVarNames.size()==0) {
            throw ModelException("FRBarrierArrayVariable::getIMSInput",
                                 "No monitorVarNames!");
        }
        result[0] = "BarrierArray";
        result[1] = name;
        // if there's only one monitor var name dispense with '{' & '}'
        string vars;
        if (monitorVarNames.size()>1) {
            vars += "{";
        }
        vars += monitorVarNames[0];
        for(int i=1; i<monitorVarNames.size(); i++) {
            vars += ","+monitorVarNames[i];
        }
        if (monitorVarNames.size()>1) {
            vars += "}";
        }
        result[2] = vars+","+(isUp?"Up":"Down")+","+ level;

        return result;
    }

    FRBarrierArrayVariable(): LValueBase(TYPE),
        isUp(false)  {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint.
        The monitoring is implicitly performed only at time points where this
        variable is referenced. */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRBarrierArrayVariable::createRValue";
        FRIfaces::IRValueBoolArray* bar = 0;
        int v = 0; // indexes through variables, on several occasions

        if (!frCtrl->getRValue(this, index)){
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')

            // Process monitoring variables...
            bool isNativeArray = false;
            FRIfaces::ILValueDoubleArray* monVarArray = 0;
            if (monitorVarNames.size() == 1) {
                // allow a DoubleArray (identified here) or a single Double (picked up below)
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression(monitorVarNames[0]);
                if (!lValueExp){
                    string m("No variable with name '"+monitorVarNames[0]+"'\n"
                             "When defining '"+name+"' you have indicated that "
                             "this is the name of a variable whose level is monitored. ");
                    throw ModelException(method, m);
                }
                FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                                 frCtrl);
                // try as double array
                monVarArray = dynamic_cast<FRIfaces::ILValueDoubleArray*>(lValue);
                if (monVarArray) {
                    isNativeArray = true;
                }
                // any other case is caught below (i.e. single Double is ok, other types fail)
            } 
            vector<FRIfaces::ILValueDouble*> monVars(monitorVarNames.size());
            if (!isNativeArray) {
                // allow an array of Doubles
                for(v=0; v<monitorVarNames.size(); v++) {
                    FRIfaces::ILValueExpression* lValueExp = 
                        frCtrl->getLValueExpression(monitorVarNames[v]);
                    if (!lValueExp){
                        string m("No variable with name '"+monitorVarNames[v]+"'\n"
                                 "When defining '"+name+"' you have indicated that "
                                 "this is the name of a variable whose level is monitored. ");
                        throw ModelException(method, m);
                    }
                    FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                                     frCtrl);
                    // cast to double
                    monVars[v] = dynamic_cast<FRIfaces::ILValueDouble*>(lValue);
                    if (!monVars[v]){
                        throw ModelException(method, "Monitoring variable '"+
                                             monitorVarNames[v]+"' must be a "
                                             "double variable");
                    }
                }
            }

            // level 
            double barLevel = 0.;
            bool   haveBarLevel = false;
            FRIfaces::ILValueExpression* lValueExpLevel = 
                frCtrl->getLValueExpression(level);
            if (lValueExpLevel){
                // it's a variable name
                FRIfaces::IRValue* rValueLevel = 
                    lValueExpLevel->getRValue(index, frCtrl);
                // cast to double
                FRIfaces::IRValueDouble* levelVar = 
                    dynamic_cast<FRIfaces::IRValueDouble*>(rValueLevel);
                if (!levelVar){
                    throw ModelException(method, "BarrierArray level variable is '"+
                                         level+"' which must be a "
                                         "double variable");
                }
                if (levelVar->isKnown()) {
                    barLevel = levelVar->getValue();
                    haveBarLevel = true;
                    if (isNativeArray) {
                        bar = new BarrierArrayNativeKnownLevel(monVarArray, isUp, barLevel);
                    } else {
                        bar = new BarrierArrayKnownLevel(monVars, isUp, barLevel);
                    }
                } else {
                    FRIfaces::ILValueDouble* lValueLevel = 
                        dynamic_cast<FRIfaces::ILValueDouble*>(lValueExpLevel->getLValue(index, frCtrl));
                    if (isNativeArray) {
                        bar = new BarrierArrayNative(monVarArray, isUp, lValueLevel);
                    } else {
                        bar = new BarrierArray(monVars, isUp, lValueLevel);
                    }
                }
            } else {
                // not a variable name - see if it can be parsed as a number
                char*   endPos;
                double  dVal = strtod(level.c_str(), &endPos);
                if (*endPos == '\0') {
                    barLevel = dVal;
                    haveBarLevel = true;
                    if (isNativeArray) {
                        bar = new BarrierArrayNativeKnownLevel(monVarArray, isUp, dVal);
                    } else {
                        bar = new BarrierArrayKnownLevel(monVars, isUp, dVal);
                    }
                } else {
                    string m("No variable with name '"+level+"'\n"
                             "When defining '"+name+"' the level must be given "
                             "as the name of a variable or directly as a value. ");
                    throw ModelException(method, m);
                }
            }
            // keep track of the event if we need to 
            if (index <= frCtrl->getEventIndex()) {
                IBarrierArrayEvent* barEv = dynamic_cast<IBarrierArrayEvent*>(bar);
                if (isNativeArray) {
                    frCtrl->addEvent(new BarrierArrayEventer(barEv, monitorVarNames[0], 
                                                            frCtrl->getDate(index), v));
                } else {
                    for(v=0; v<monitorVarNames.size(); v++) {
                        frCtrl->addEvent(new BarrierArrayEventer(barEv, monitorVarNames[v], 
                                                            frCtrl->getDate(index), v));
                    }
                }
            }
            // Consider reporting BARRIER_LEVELs
            const DateTime& today = frCtrl->productView()->getValueDate();
            DateTime upperDate = BarrierLevel::barrierWindow(today);
            const DateTime& thisDate = frCtrl->getDate(index);
            if (thisDate>=today &&
                thisDate<=upperDate &&
                haveBarLevel) {
                if (!isNativeArray) {
                    // allow an array of Doubles
                    for(v=0; v<monitorVarNames.size(); v++) {
                        FRIfaces::ILValueExpression* lValueExp = 
                            frCtrl->getLValueExpression(monitorVarNames[v]);
                        // we know lValueExp is ok from above check
                        FRIfaces::IVarBarrierLevelAssist* barLevelAssist = 
                            dynamic_cast<FRIfaces::IVarBarrierLevelAssist*>(lValueExp);
                        if (barLevelAssist) {
                            string assetName = barLevelAssist->getAssetName();
                            frCtrl->addBarrierLevel(assetName,
                                                    isUp,
                                                    barLevel,
                                                    thisDate,
                                                    barLevelAssist);
                        }
                    }
                } else {
                    FRIfaces::ILValueExpression* lValueExp = 
                        frCtrl->getLValueExpression(monitorVarNames[0]);
                    FRIfaces::IVarArrayBarrierLevelAssist* barrayLevelAssist = 
                        dynamic_cast<FRIfaces::IVarArrayBarrierLevelAssist*>(lValueExp);
                    if (barrayLevelAssist) {
                        vector<FRIfaces::IVarBarrierLevelAssist*> compAssists = 
                            barrayLevelAssist->getComponentAssists(frCtrl);
                        for (unsigned int i=0; i<compAssists.size(); i++) {
                            if (compAssists[i]) {
                                string assetName = compAssists[i]->getAssetName();
                                frCtrl->addBarrierLevel(assetName,
                                                        isUp,
                                                        barLevel,
                                                        thisDate,
                                                        compAssists[i]);
                            }
                        }
                    }
                }
            }
        }
        return bar;
    }

private:
    // field
    const string       name;
    const CStringArray monitorVarNames;// which var is checked against level
    const string       barrierType;    // "Up" or "Down". String rather than bool allows easier extension
                                       // and greater clarity in the interface
    const string       level;                // var name or double value

    // transient
    bool isUp;                   // translated from barrierType

    FRBarrierArrayVariable(const string       varName,
                           const CStringArray monitorVarNames,
                           const string       barrierType,
                           const string       level): 
        LValueBase(TYPE),
        name(varName), monitorVarNames(monitorVarNames),
        barrierType(barrierType), level(level) {
            validatePop2Object();
    }

    static IObject* defaultFRBarrierArrayVariable(){
        return new FRBarrierArrayVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRBarrierArrayVariable::create";
        try {
            if (value.size()!=1) {
                throw ModelException(routine, 
                                     "Per sim date initialisation not yet supported here");
            }
            if (value[0].empty() || value[0] == "-") {
                throw ModelException(routine, "No definition provided!");
            }

            // We allow the case where there's a single var name to be without braces
            // var names we check against 'level'
            StringArraySP monitorVarNames;
            size_t openBrace = value[0].find("{");
            size_t closeBrace = value[0].find("}");
            unsigned int remIdx = 0;
            if (openBrace == string::npos && closeBrace == string::npos) {
                // no braces - must have a single var name
                monitorVarNames = StringArraySP(new StringArray(0));
                string nm(value[0],0,value[0].find(","));
                monitorVarNames->push_back(nm);
                remIdx = value[0].find(",");
            } else if (openBrace != string::npos && closeBrace != string::npos) {
                // both open and close brace supplied
                // Note the length is 2 short in order to strip off the braces
                string monNamesValue(value[0],openBrace+1,closeBrace-openBrace-1);
                if (monNamesValue.empty()) {
                    throw ModelException(routine, "No monitoring variables specified!");;
                }
                IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                               StringArray::TYPE, monNamesValue));
                monitorVarNames = StringArraySP(StringArraySP::dynamicCast(vars));
                remIdx = value[0].find("}")+1; // +1 skips the brace
            } else {
                string msg("Unpaired brace - missing ");
                msg += (openBrace == string::npos)?"opening":"closing";
                msg += " brace.";
                throw ModelException(routine, msg);
            }

            string remValue(value[0],remIdx,string::npos);
            string barrierType; // up/down
            string levelStr; // var name or double value that defines 'level' 

            const char* pos = remValue.c_str();
            int paramId = 0;
            while (*pos) {
                char   c;
                // skip whitespace
                while ((c = *pos) == ' ' || c == '\t'){
                    pos++;
                }
                if (c == ','){
                    // skip over our chosen separator and continue
                    pos++;
                    // count params
                    paramId++;
                    if (paramId>2) {
                        throw ModelException(routine, 
                                             "Too many parameters for "
                                             "barrier array variable " + varName + 
                                             ". Count the commas!");
                    }
                } else {
                    // next contiguous block is the text we want
                    const char* varBegin = pos;
                    do{
                        pos++; /* Get another character. */
                        c = *pos;
                    } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                    const char* varEnd = pos;
                    // done with this one
                    switch (paramId) {
                    case 0:
                        // done above
                        break;
                    case 1:
                        barrierType = string(varBegin, varEnd - varBegin);
                        break;
                    case 2:
                        levelStr = string(varBegin, varEnd - varBegin);
                        break;
                    }
                } 
            }

            if (paramId != 2) {
                throw ModelException(routine, "Wrong number of parameters for barrier variable " + varName);
            }

            // A simple "return new ..." gives gcc compilation warning 
            FR::LValueBase* bar = new FRBarrierArrayVariable(varName,
                                                             *monitorVarNames.get(),
                                                             barrierType,
                                                             levelStr);
            return bar;
        } catch (exception& e) {
            throw ModelException(e, routine,
                                 "Failed to parse value string: " + value[0] + 
                                 "\nFormat required : "
                                 "{monitorVarNames},U/D,level\n"
                                 "monitorVarNames being a comma-separated list."
                                 "Braces can be omitted with a single monitorVarName.");
        }
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBarrierArrayVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRBarrierArrayVariable);
        FIELD(name, "name of variable");
        FIELD(monitorVarNames, "Array of names of variables to monitor versus "
                     "level");
        FIELD(barrierType, "Up/U or Down/D");
        FIELD(level, "Name of variable, or double value, that defines "
                     "monitoring level");
        FIELD_NO_DESC(isUp);
        FIELD_MAKE_TRANSIENT(isUp);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_BARRIER_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX BarrierArray Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRBarrierArrayVariable::create, 
                                           "BarrierArray");
    }
    
};

CClassConstSP const FRBarrierArrayVariable::TYPE = CClass::registerClassLoadMethod(
    "FRBarrierArrayVariable", typeid(FRBarrierArrayVariable), FRBarrierArrayVariable::load);

const string FRBarrierArrayVariable::UP = "Up";
const string FRBarrierArrayVariable::DOWN = "Down";



DRLIB_END_NAMESPACE
