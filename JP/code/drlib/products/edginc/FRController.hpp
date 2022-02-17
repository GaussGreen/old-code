//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRController.hpp
//
//   Description : Defines common concrete classes used internally
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FRCONTROLLER_HPP
#define EDR_FRCONTROLLER_HPP

#include <memory>
#include "edginc/MonteCarlo.hpp"
#include "edginc/FRIfaces.hpp"
#include "edginc/BarrierLevel.hpp"

DRLIB_BEGIN_NAMESPACE
class Hashtable;
class FlexBarrierBreach;
class EventResults;

class PRODUCTS_DLL FRController {
public:
    /** interface for object that want to be notified every time payoff is
        called */
    class PRODUCTS_DLL INotifyPayoffCall{
    public:
        virtual void update(const FRController* ctrl) = 0;
        virtual ~INotifyPayoffCall();
    };

    /** add an object to the list of objects that get called when the
        payoff is invoked. If memoryManage is true the supplied object
        will be freed when this object is deleted. */
    void addToNotifyPayoffCall(INotifyPayoffCall* obj, bool memoryManage);

    // interface for flex objects which can handle flex breach events
    class PRODUCTS_DLL IFREvent {
    public:
        virtual FlexBarrierBreach* createEvent() = 0;
        virtual bool hasEvent() = 0;
        virtual ~IFREvent() {};
    };
    typedef refCountPtr<IFREvent> IFREventSP;

    /** Returns true if use of the state var framework is requested */
    bool stateVarUsed() const;

    /** Flex variables that use state variables need to create an instance
        of this class and pass it to addToStateVars */
    class PRODUCTS_DLL IStateVarClient{ // possibly should derive from IMCStateVarClient
    public:
        virtual ~IStateVarClient();
        /** Populates the collector with the state variables required by the
            various assets */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const = 0;

        /** To be called when the path generator changes (currently before doing
            the past, and before doing the future) */
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) = 0;
    };
        
    /** add an object to the list of flex variables that use state variables.
        If memoryManage is true the supplied object
        will be freed when this object is deleted. */
    void addToStateVars(IStateVarClient* obj, bool memoryManage);

    /** Builds an FRIfaces::IAssignment for the given l and r values eg maps
        double r value to double l value. If the r value is known then the
        lvalue is assined to and null is returned */
    static FRIfaces::IAssignment* createAssignment(FRIfaces::ILValue* lvalue,
                                                   FRIfaces::IRValue* rvalue);

    /** returns the index of the current point */
    int getIndex() const;

    /** Returns how many simulation dates there are (ie max value of
        getIndex()) */
    int numDates() const;

    /** Returns the path generator for the current path. Value
        undefined until payoff() invoked */
    const MCPathGenerator*  getPathGenerator() const;

    /** returns the date corresponding to index */
    const DateTime& getDate(int index) const;

    /** returns view of product */
    const FRIfaces::IProductView* productView() const;

    /** Returns the IRValue object for supplied variable at 
        specified index. Should only be used by implementations of
        IRValueExpression::getRValue as it returns null if the corresponding
        IRValueExpression has not set it yet */
    FRIfaces::IRValue* getRValue(
        const FRIfaces::IRValueExpression* expression,
        int                                index) const;
    
    /** Returns the IRValue object for supplied variable at 
        the current index. Should only be used by implementations of
        IRValueExpression::getRValue as it returns null if the corresponding
        IRValueExpression has not set it yet  */
    FRIfaces::IRValue* getRValue(
        const FRIfaces::IRValueExpression* expression) const;

    /** Sets/stores the RValueSP for current index for the specified 
        expression */
    void setRValue(const FRIfaces::IRValueExpression* expression,
                   const FRIfaces::IRValueSP&         rValue);

     /** Sets/stores the RValueSP for the specified index for the specified 
        expression */
    void setRValue(int                                index,
                   const FRIfaces::IRValueExpression* expression,
                   const FRIfaces::IRValueSP&         rValue);

    /** Given variable name returns LValueExpression representing that
     * variable. Note returns NULL if no such variable */
    FRIfaces::ILValueExpression* getLValueExpression(
        const string& varName) const;

    /** Returns a vector of IRValues which reflect how the given
        variable is calculated at each simulation date */
    const vector<FRIfaces::RValUnion>& getRValuesForLValueExpression(
        const FRIfaces::ILValueExpression* lValueExp);

   /** Creates FRController object ready for simulation run. The 3rd
        parameter controls whether const variables can be resolved
        immediately or must wait until the simulation (this is for
        testing purposes)  */
    FRController(const FRIfaces::IProductView* productView,
                 bool                          useStateVars, // new style MC?
                 bool                          allowResolutionOfConstVars,
                 bool                          triggerEvents);

    virtual ~FRController();

    /** evaluates one path */
    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices);

    /** evaluates one path gathering debug info */
    void payoffWithDebug(const MCPathGenerator*  pathGen,
                         IMCPrices&              prices);

    /** returns a hashtable holding the value of each variable in the latest
        simulated path - for debugging purposes as s/sheet */
    Hashtable* getCurrentValues() const;

    /** More debug ...
     */
    void updateDebugValues(bool doMin,
                           bool doMax,
                           bool doAvg);
    Hashtable* getDebugValues(const string& id) const;

    /** easiest to handle recording and reporting internally to 
        controller where the info all resides. We need to be able
        to trigger it from the instrument though. */
    void recordPaymentDates(CControl*     control,
                            Results*      results);
    void recordKnownCashFlows(CControl*     control,
                              Results*      results,
                              const string& ccyName) const;
    
    /** tester. Ignores value of productView->getPayVariable() and
        uses supplied value. Can cope with any types of variable for
        'pay' variable. Returns 2D array of values of 'pay' variable
        for each date for each pay variable. Runs simulation across
        all dates. */
    ObjectArraySP tester(
        const FRIfaces::ILValueExpressionArray* payVars);

    /** memory management utility - stores pointer and frees it when 
        FRController goes out of scope. Repeated calls with the same
        pointer are ok */
    void store(FRIfaces::IRValue*  object);

    /** memory management utility - stores pointer and frees it when 
        FRController goes out of scope */
    void store(const FRIfaces::IRValueSP&  rValueSP);

    /** state management utility for objects which store state which needs to
        be reset between each simulation */
    void registerStateObject(FRIfaces::IHoldsState* obj);

    /** Populates the collector with the state variables required by the
        various assets. To be called by the containing IMCProduct */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** To be called when the path generator changes (currently before doing
        the past, and before doing the future).
        To be called by the containing IMCProduct */
    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    bool getTriggerEvents() const;

    int getEventIndex() const;

    void addEvent(IFREvent* event);

    void retrieveEvents(EventResults* events);

    void addBarrierLevel(const string&   assetName,
                         bool            isUp,
                         double          barLevel,
                         const DateTime& barrierDate,
                         const FRIfaces::IVarBarrierLevelAssist* barLevelAssist);
    
    void getBarrierLevelReports(vector<string>&              assetNames,
                                vector<BarrierLevelArraySP>& levels);

private:
    void assignValues(int start, int end);
    static void populateRValueUnion(
        const FRIfaces::ILValueExpression* lValue,
        FRIfaces::IRValue*                 rValue,
        FRIfaces::RValUnion&               rValUnion);
    
   
    class AssignmentObject;
    class AssignmentDouble;
    template <class L, class R> class AssignmentArray;
    //class AssignmentDoubleArray;
    class AssignmentDoubleSpecial1;
    class AssignmentDoubleSpecial2;
    class AssignmentDoubleFromInt;
    class AssignmentInt;
    class AssignmentIntSpecial1;
    class AssignmentBool;
    class AssignmentBoolSpecial1;
    class AssignmentDate;
    class AssignmentDateSpecial1;
    class Imp;
    friend class Imp;
    Imp*    my; // hides implementation
};

DRLIB_END_NAMESPACE
#endif
