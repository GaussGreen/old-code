//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceLayer.hpp
//
//   Description : tree slice layer, one state variable per layer
//
//   Date        : May 19, 2006
//----------------------------------------------------------------------------

#ifndef TREE_SLICE_LAYER_HPP
#define TREE_SLICE_LAYER_HPP

DRLIB_BEGIN_NAMESPACE

class FDModel;

class TreeSliceLayer;
typedef refCountPtr< TreeSliceLayer > TreeSliceLayerSP;

// ********* TreeSliceLayer class ************
class TREE_DLL TreeSliceLayer : public TreeSlice
{
public:
    ///////// StateSupport [begin]: state variable support using slice layers ////////
    class TREE_DLL StateSupport
    {
        //************ types defined **************/
    public:
        typedef enum { INTERP_LINEAR, INTERP_QUADRATIC } InterpMethod;

        //************ variables defined ************/
    protected:
        double todayValue;
        vector< double > prevGrid;
        vector< double > currGrid;

    private:
        InterpMethod interpMethod;
        bool isFwdTransition;
        // isFwdTransition == true: transition grid from t- to t+
        // isFwdTransition == false: transition grid from t+ to t-

        TreeSliceLayerSP gridLayer;

        const FDModel * model;
        int level; /* When more than one state variable are stacked altogether
                      the one with the smallest "level" is closer to the concrete
                      slice implementation. If level==0, then slices[0] is never
                      a TreeSliceLayer. Set temporarily to -1 by constructor. */

        int currInterpStep; // time step at which the latest slice interpolation has been performed
        int prevInterpStep; // previous value of currInterpStep. Used only between interpPrepare() and interpPerform()

        string name; // for debug purposes

        //************* functions defined ***********/
    public:

        // Transitions the 'grid' using 'srcSlices'
        // To be implemented by inherited classes
        virtual void transition(
            const vector< const TreeSlice * > & srcSlices,
            const vector< double > & gridIn,
            vector< double > & gridOut ) const = 0;

        // Called from FDModel::registerStateSupport
        // to set model and layer level 
        void init( const FDModel * model, int level );

        // Builds a slice representing the currGrid
        void populateGrid( int step );

        // Returns current 'gridLayer' populated during last call to 'populateGrid'
        const TreeSlice & getGridSlice() const { return *gridLayer; }

        // Interpolates all the 'layers'
        // 'srcSlices' will be iterated along and passed into 'transition'
        void interpolate( int step, const vector< const TreeSlice * > & srcSlices );

        // Interpolates at 'todayValue' and returns corresponding 'layer' value 
        double getPrice0( const TreeSliceLayer & layer ) const;

        // Return the last time step at which an interpolation has been performed.
        // If a slice has a lastInterpStep greater than getCurrInterpStep() it 
        // means that something went wrong with it.
        int getCurrInterpStep() { return currInterpStep; }

        // Checks that lastInterpStep has an allowed value.
        // INPUT: currInterpStep, layer, interpStep
        // OUTPUT: exception thrown if invalid lastInterpStep
        void checkInterpStep( const TreeSliceLayer & layer, bool usePrev = false ) const;

        // constructor / destructor
        StateSupport(
            const string & name,
            double todayValue,
            InterpMethod interpMethod,
            bool isFwdTransition = true );
        virtual ~StateSupport() {}
    };
    ///////// StateSupport [end] ////////

    /****** TreeSlice virtual function implementation ******/
    virtual void printDetails(char *s) const;
    virtual const string typeName() const { return "TreeSliceLayer"; }
    virtual TreeSliceSP clone( bool copyValues = true ) const;
    virtual TreeSlice& operator =( double v );
    virtual double getCentre() const;
    virtual bool isZero() const;
    double getPrice0() const;

    template< typename A >
    inline void eval( const A & expr );

//    template<class T> // T must define sliceCount, compute() and printDebug()
//    void loopOnSlices(const T & oper, TreeSliceRates **s, int nbOutput);

    // !!! TO BE REMOVED
    virtual void getCalcRange( int & bot, int & top ) const;
    virtual double * getValues() const;

    /****** TreeSliceLayer-specific ******/

    TreeSliceLayer(
        StateSupport * owner,
        const TreeSlice & slice,
        bool copyValues,
        int lastInterpSatet = -1 );

    /** collects all the 'layers' of a given 'owner'
        if 'expandedOnly' then skip layers that are not expanded in state variable dimension  */
    static void collectLayers(
        const StateSupport * owner,
        const vector< TreeSliceSP > & slices,
        vector< TreeSliceLayer * > & layers,
        bool expandedOnly = true );

    // accessor function for StateSupport and Model (to DEV what is inside)
    const vector< TreeSliceSP > & getSlices() const { return slices; }
    vector< TreeSliceSP >       & getSlices()       { return slices; }

    vector< TreeSliceSP >       & getTempSlices()   { return tempSlices; }

private:

    vector< TreeSliceSP > tempSlices; // temporary slices for interpolation
    vector< TreeSliceSP > slices; // vector of tree slices, one for each state variable value
    int lastInterpStep;           // step at which the last interpolation (update) occurred
    StateSupport * owner;

    mutable LoopList< const TreeSlice * > loopSlices;
    mutable LoopList< const TreeSliceLayer * > loopLayers;
    mutable LoopList< ExposedSliceOverride > loopOverride;
};

#if 0
template<class T>
void loopOnSlices(const T & oper, TreeSlice **s, int nbOutput) 
{
    const TreeSliceLayer** layerArgs =  loopLayers.reserve( oper.sliceCount );
    int outputLayerLevel = dynamic_cast<const TreeSliceLayer&>(*s[0]).level;
    int newNbSlices = 1;
    for (int n = nbOutput; n<oper.sliceCount; ++n) {
        
        const TreeSliceLayer *layer = dynamic_cast<const TreeSliceLayer*>(s[i]);
        if (layer) {
            if (layer->level > outputLayerLevel ) {
                throw ModelException("Output slice "+Format::toString(i)
                +" cannot support the state variable provided by input #"+laver->owner->name);
            }
            if (layer->level < outputLayerLevel ) 
                layer=0;
        }
        layerArgs[i] = layer;
        if (layer)
            newNbSlices = max(newNbSlices, (int)layer->slices.size());
    }
    for (int n = 0; n<nbOutput; ++n) {
        if (s[i]->level != outputLayerLevel ) {
            throw ModelException("Output slices "+Format::toString(i)
            +" are not all the same type");
        }
        TreeSlice::resizeVector( newNbSlices, *s[i]->slices[0], s[i]->slices, s[i]->name, true);
    }

    for (int i=0; i<newNbSlices; ++i) // for each value of the state variable
    {
        vector< TreeSlice::ExposedSliceOverride > override( expr.sliceCount );
        for (int j=0; j<expr.sliceCount; ++j) // for each slice appearing in the expression
        {
            const TreeSliceLayer *layer = layerArgs[j];
            // If the slice appearing at the j-th position in the expression is of type
            // TreeSliceLayer then replace it with the underlying slice corresponding
            // to the same value of the state variable
            if( layer )
                override[j].set( layer, layer->slices[ layer->slices.size() > 1 ? i : 0 ].get() );
            // If the slice is not of type layer (or refers to another state variable)
            // do not override (e.g. use the current slice). This provides the 
            // auto-expansion feature
        }
        // do the operation on the underlying slices
        loopOnSlices(oper, s, nbOutput);
        slices[i]->evalMultiType(expr);
    }    
    
}
#endif

template< typename A >
inline void TreeSliceLayer::eval( const A & expr )
{
    try {
        const char* (*print)(void*) = printExpr<A>; 
        (void)print; // to remove compiler warnings

        // watch "print((void*)expr)" and "*this->iter" in your debugger

        // get a list of the slices used in the expression
        const TreeSlice** treeSlices = loopSlices.reserve( expr.sliceCount );
        const TreeSlice** end = expr.listSlices(treeSlices);
        ASSERT((end-treeSlices)==expr.sliceCount);

        const TreeSliceLayer** layerArgs = loopLayers.reserve( expr.sliceCount );
        int newNbSlices = 1;
        int newInterpStep = -1;
        int newTreeStep = -1;

        for (int i=0; i<expr.sliceCount; ++i) { // for each slice appearing in the expression
            const TreeSliceLayer *layer = dynamic_cast<const TreeSliceLayer*>(treeSlices[i]);
            if( layer && ! ( layer->owner == owner ) )
            {
                // If layerId are different it means that the layers are
                // dealing with different state variables, so ignore.
                layer=0;
            }
            layerArgs[i] = layer;
            if (layer) {
                newNbSlices = max(newNbSlices, (int)layer->slices.size());
                layer->owner->checkInterpStep(*layer);
                
                if (newInterpStep<0)
                    newInterpStep = layer->lastInterpStep;
                    
                if (layer->lastInterpStep>=0)
                    newInterpStep = min (newInterpStep, layer->lastInterpStep);
            }
            if (treeSlices[i]->treeStep >= 0) {
                if (newTreeStep <0 ) {
                    newTreeStep = treeSlices[i]->treeStep;
                }
                // need more work to pass this test
//                else if (treeSlices[i]->treeStep != newTreeStep) {
//                    throw ModelException("Inconsistent treeStep");
//                }
            }
        }
        resizeVector( newNbSlices, *slices[ 0 ], slices, name, true);
        lastInterpStep = newInterpStep;
        treeStep = newTreeStep;

        ExposedSliceOverride * override = loopOverride.reserve( expr.sliceCount );
        for (int i=0; i<(int)slices.size(); ++i) // for each value of the state variable
        {
            int k = 0;
            for (int j=0; j<expr.sliceCount; ++j) // for each slice appearing in the expression
            {
                const TreeSliceLayer *layer = layerArgs[j];
                // If the slice appearing at the j-th position in the expression is of type
                // TreeSliceLayer then replace it with the underlying slice corresponding
                // to the same value of the state variable
                if( layer ) {
                    ASSERT((size_t)i<layer->slices.size() || layer->slices.size()==1);
                    new (override + k) ExposedSliceOverride;
                    override[k].set( layer, layer->slices[ layer->slices.size() > 1 ? i : 0 ].get() );
                    ++k;
                }
                // If the slice is not of type layer (or refers to another state variable)
                // do not override (e.g. use the current slice). This provides the 
                // auto-expansion feature
            }
            // do the operation on the underlying slices
            slices[i]->evalMultiType(expr);
            while( --k >= 0 )
                override[k].~ExposedSliceOverride();
        }
    }
    catch (exception& e) {
        char buf[10000];
        *buf=0;
        expr.printDebug(buf);
        throw ModelException(e, "TreeSliceLayer::eval(const A & expr), "+string(buf));
    }
}

DRLIB_END_NAMESPACE

#endif
