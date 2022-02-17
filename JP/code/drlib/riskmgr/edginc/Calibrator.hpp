//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Calibrator.hpp
//
//   Description : Calibrator
//
//   Author      : regis Guichard
//
//   Date        : 21 May 02
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CALIBRATOR_H
#define EDG_CALIBRATOR_H
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Optimizer.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

// Forward declaration of IInstanceIDBootstrapper
class IInstanceIDBootstrapper;
typedef smartPtr<IInstanceIDBootstrapper> IInstanceIDBootstrapperSP;
typedef smartConstPtr<IInstanceIDBootstrapper> IInstanceIDBootstrapperConstSP;
typedef array<IInstanceIDBootstrapperSP, IInstanceIDBootstrapper> IInstanceIDBootstrapperArray;
typedef smartPtr<IInstanceIDBootstrapperArray> IInstanceIDBootstrapperArraySP;

// Forward declaration of IBootstrapper
class IBootstrapper;
typedef smartPtr<IBootstrapper> IBootstrapperSP;
typedef smartConstPtr<IBootstrapper> IBootstrapperConstSP;
typedef array<IBootstrapperSP, IBootstrapper> IBootstrapperArray;
typedef smartPtr<IBootstrapperArray> IBootstrapperArraySP;

/** Calibrator Sensitivity. */
/** Calibrator Sensitivity. */
class RISKMGR_DLL Calibrator: public CObject{
public:
    friend class CalibratorHelper;
    static CClassConstSP const TYPE;
    /** Classes implement this interface to support their fields being
        'calibrated'. The 'reflection' mechanism is used to get/set the
        value of the fields. Developers need to be aware of this 'behind
        the scenes' alterations to the object data members. In particular
        they want to override the default implementation of
        IObject::fieldsUpdated() */
    class RISKMGR_DLL IAdjustable: public virtual IObject{
    public:
        static CClassConstSP const TYPE;
        /** Returns the name of this object. This name will be used for
            returning tweaks (and for calibration will be input data) so
            the name should reflect the data's market data name. */
        virtual string getName() const = 0;

        virtual ~IAdjustable();

        /** Registers the specified field for the given class with the
            Calibrator. Only registered fields will considered by the
            calibrator. The method takes ownership of the Range supplied */
        static void registerField(const CClassConstSP& clazz,
                                  const string&        fieldName,
                                  Range*               range);

        /** Registers the specified field for the given class with the
            Calibrator. Only registered fields will considered by the
            calibrator. The registered field is associated with an array
            of expirires via the method getExpiries.
            The method takes ownership of the Range supplied */
        typedef ExpiryArraySP (*TGetExpiriesMethod)(const IObject* obj);
        static void registerBootstrappableField(
            const CClassConstSP& clazz,
            const string&        fieldName,
            Range*               range,
            TGetExpiriesMethod   getExpiriesMethod);

        /** Returns all fields that are registered with the calibrator
            for the specified class */
        static CFieldArray getFields(const CClassConstSP& clazz);

        /** Checks fields (which can be calibrated) in object are
            within range */
        static void checkRange(const IAdjustable* adjustable);

        /** Returns the recorded range for the specified field */
        static const Range& getRange(const CFieldConstSP& field);

        static TGetExpiriesMethod getGetExpiriesMethod(const CFieldConstSP& field);

        static bool hasGetExpiriesMethod(const CFieldConstSP& field);

       

        //private: // MS can't cope with this being private
        class FieldInfo;
    };

    DECLARE(IAdjustable)

    class InstanceIDDb; // derived class for doubles
    class InstanceIDDbArray; // derived class for DoubleArrays
	typedef smartPtr<InstanceIDDb> InstanceIDDbSP;
	typedef smartPtr<InstanceIDDbArray> InstanceIDDbArraySP;

    class ObjFunc;
	typedef smartPtr<ObjFunc> ObjFuncSP;

	// InstanceID is a polymorphic type and we want to handle arrays of them
    class InstanceID;
    typedef smartPtr<InstanceID> InstanceIDSP;
    typedef array<InstanceIDSP, InstanceID> InstanceIDArray;
    typedef smartPtr<InstanceIDArray> InstanceIDArraySP;
    /** Identifies field within object that is being calibrated */
    class RISKMGR_DLL InstanceID: public CObject{
    protected:
        string             typeToCalibrate;
        string             name;             // name of instance
        string             fieldName;        // field to calibrate
        // looked up upon construction
        CClassConstSP      classToCalibrate; // what type we're looking for $unregistered
        CFieldConstSP      field;            // field to calibrate $unregistered

        InstanceID(const CClassConstSP&  clazz,
                   const string&  typeToCalibrate,
                   const string&  name,
                   const string&  fieldName);

        // for default constructor
        InstanceID(const CClassConstSP&  clazz);

    private:
        static void load(CClassSP& clazz);

    public:

        /** 
        * Returns the object corresponding to typeToCalibrate for this InstanceID
        */
        IObjectConstSP getObject(const IObjectSP& adjGroup) const;
        
        /** access to field name*/
        string  getFieldName() const;   

        static CClassConstSP const TYPE;

        virtual void validatePop2Object();

        virtual IObject* clone() const;

        // for initial guess ////////
        static InstanceIDSP dbArray(const string&         typeToCalibrate,
                                    const string&         name,
                                    const string&         field);

		static void initialise(InstanceIDSP ids, ObjFuncSP objFunc);

		static ExpiryArraySP getExpiries(InstanceIDSP ids, ObjFuncSP objFunc);

		//					////////

        virtual ~InstanceID(){}

        /** Must be called before any other methods - allows object to see the
            market data (eg length of arrays) */
        virtual void initialise(const IObjectSP&  adjGroup) = 0;

        /** Returns the number of variables */
        virtual int numVariables() const = 0;

        /** Populates the names array from 'index' to
            'index' + numVariables() */
        virtual void getNames(StringArray& names, int index) const = 0;

        /** Populates the ranges array from 'index' to
            'index' + numVariables() */
        virtual void getRanges(RangeArray& ranges, int index) const = 0;

        /** Populates the values array from 'index' to
            'index' + numVariables() */
        virtual void getInitialValues(
            DoubleArray&         values,
            int                  index,
            const IObjectSP&     adjGroup) const = 0;

        /** Populates the values array from 'index' to
            'index' + numVariables() */
        virtual void getValues(
            DoubleArray&         values,
            int                  index,
            const IObjectSP&     adjGroup) const = 0;

        /** Uses the values in x from 'index' to
            'index' + numVariables() to adjust the field being calibrated */
        virtual void applyAdjustment(const IObjectSP&     adjGroup,
                                     const CDoubleArray&  x,
                                     int                  index) const = 0;

        /** Put the calibrated results into some sort of context -
            not clear what is the best way of doing this at the moment */
        virtual IObjectSP writeResults(const DoubleArray& calibratedVals,
                                       int                index) const = 0;

        /** Read off the 'values' object and populate calibratedVars array */
        virtual void readResults(const IObject& values,
                                 DoubleArray&   calibratedVals,
                                 int            index) const = 0;

        /** Returns true if scalar type; false if vector type */
        virtual bool isScalar() const = 0;

        /** If vector type, returns indices of elts that got calibrated.
            Returns 0, otherwise. */
        virtual IObjectSP getArrayIndices() const = 0;

        /** Applies adjustments to array of InstanceID's */
        static void applyAdjustment(const InstanceIDArray& ids,
                                    const IObjectSP&       adjGroup,
                                    const CDoubleArray&    x);

        /** Gets names from array of InstanceID's */
        static void getNames(const InstanceIDArray& ids,
                             StringArray&           names);

        /** Gets ranges from array of InstanceID's */
        static void getRanges(const InstanceIDArray& ids,
                              RangeArray&            ranges);

        /** Gets initial values from array of InstanceID's */
        static void getInitialValues(const InstanceIDArray& ids,
                                     const IObjectSP&       adjGroup,
                                     DoubleArray&           x);

        /** Gets current values from array of InstanceID's */
        static void getValues(const InstanceIDArray& ids,
                              const IObjectSP&       adjGroup,
                              DoubleArray&           x);

        enum TResultIndices{
            OBJ_TYPE = 0,
            NAME,
            FIELD_NAME,
            VALUE,
            DB_TYPE,
            ARRAY_INDEX,
            NB_RES
        };

        /** Put the calibrated results into some sort of context -
            not clear what is the best way of doing this at the moment */
        static IObjectSP writeResults(const InstanceIDArray& ids,
                                      const DoubleArray&     calibratedVals);

        /** Read off the res object and return a pair of populated ids and
            calibrated values */
        static pair<InstanceIDArray, DoubleArray>
                         readResults(const IObjectConstSP& res);

	    /** InstanceIDDbArray class implements this interface to support the use of the
            bootstrapping routine for the calibration. */
		class IBootstrappable;
		typedef smartPtr<IBootstrappable> IBootstrappableSP;
		typedef array<IBootstrappableSP, IBootstrappable> IBootstrappableArray;

	    class RISKMGR_DLL IBootstrappable: public virtual IObject{
	    public:
		    static CClassConstSP const TYPE;

            static void load(CClassSP& clazz);

		    /** Extracts the terms associated with that IBoostrappable object */
            virtual ExpiryArraySP getExpiries(const IObjectConstSP& adjGroup) const = 0;

		    /** Extracts the InstanceID corresponding to the index idxMat */
		    virtual InstanceIDSP getInstanceID(int idxMat) const = 0;
	    };
	};

    /** Set of methods that need to be implemented by the objective function
        and will be invoked by the calibrator */
    class RISKMGR_DLL ObjFunc: public CObject{
    public:
        static CClassConstSP const TYPE;

        /** Give the derived class a chance to do s'thing with the market */
        virtual void getMarket(MarketData* market) = 0;

        /** Helper class */
        class RISKMGR_DLL Helper{
        public:
            /** Helper function to get the market */
            static void getMarket(CInstrument* inst,
                                  IModel*      model,
                                  CControl*    control,
                                  MarketData*  market);
        };

        /** Gives the derived class a chance to do additional validation
            irrespective of getMarket being called or not */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup() = 0;

        /** Calculates the value of the objective function */
        virtual double calcValue() const = 0;

	    /** ObjFunc class implements this interface to support the use of the
            bootstrapping routine for the calibration. */
        class RISKMGR_DLL IBootstrappable: public virtual IObject{
        public:
            static CClassConstSP const TYPE;

            static void load(CClassSP& clazz);

            /** Makes additional validations for calibration with bootstrapping */
            virtual void validate(const InstanceID::IBootstrappableArray& ids) const = 0;

            /** Updates the instrument for calibration with bootstrapping
                before each calibrator run */
		    virtual void update(int idxMat) = 0;

            /** Reset the instrument for calibration with bootstrapping
                after each calibrator run */
            virtual void reset() = 0;
        };
		typedef smartPtr<IBootstrappable> IBootstrappableSP;

		/**
		 * Interface an objective function should implement to support
		 * "generic bootstrapping"
		 * ["generic" means independent of eg asset class; for the moment
		 * ObjFunc::IBootstrappable is still used for equity vol and SRM calibration]
		 *
		 
		 * */
		 // TODO opportunity to derive from IBootstrapper ?
        class RISKMGR_DLL IGenericBootstrappable: public virtual IObject{
        public:

        	/** Used for reflection */
            static CClassConstSP const TYPE;

        	/** Used for reflection */
            static void load(CClassSP& clazz);

            /** main method used to retrieve a bootstrapper capable of bootstrapping */
            virtual IBootstrapperSP getBootstrapper(
                const InstanceIDArray & ids) = 0;
			
		};
		typedef smartPtr<IGenericBootstrappable> IGenericBootstrappableSP;

    protected:
        /** for reflection */
        ObjFunc(const CClassConstSP& clazz);

    private:
        static void load(CClassSP& clazz);
    };
    //    typedef smartPtr<ObjFunc> ObjFuncSP;
    typedef array<ObjFuncSP, ObjFunc> ObjFuncArray;

    /** Additional set of methods that need to be implemented by the objective
        function if the calibrator is to operate with a least square optimizer */
    class RISKMGR_DLL ObjFuncLeastSquare: public ObjFunc{
    public:
        static CClassConstSP const TYPE;

        /** Returns the number of functions */
        virtual int getNbFuncs() const = 0;

        /** Calculates the values of the objective functions */
        virtual void calcValue(CDoubleArray&   funcvals) const = 0;

        /** Calculates the value of the objective function */
        virtual double calcValue() const;

        static const string DO_NOT_NORMALIZE;

        void isObjFuncNormalized(string doNormalize);

        string getNormalizedWeights();


    protected:
        /** for reflection */
        ObjFuncLeastSquare(const CClassConstSP& clazz);

        string objFuncNormalized;

    private:
        static void load(CClassSP& clazz);

        /* Transient field */
        mutable CDoubleArraySP funcs;
    };
    typedef smartPtr<ObjFuncLeastSquare> ObjFuncLeastSquareSP;
    typedef array<ObjFuncLeastSquareSP, ObjFuncLeastSquare> ObjFuncLeastSquareArray;

    /** Linear combination of a set of ObjFunc's */
    class RISKMGR_DLL ObjFuncCombo: public ObjFunc,
                        public virtual ObjFunc::IBootstrappable{
    public:
        static CClassConstSP const TYPE;

        void validatePop2Object();

        /** Give the derived class a chance to do s'thing with the market */
        virtual void getMarket(MarketData* market);

        /** Give the derived class a chance to do additional validation
            irrespective of getMarket being called or not */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup();

        /** Calculates the value of the objective function */
        virtual double calcValue() const;

        /** Makes additional validations for calibration with bootstrapping */
        virtual void validate(const InstanceID::IBootstrappableArray& ids) const;

        /** Updates the instrument for calibration with bootstrapping
            before each calibrator run */
		virtual void update(int idxMat);

        /** Reset the instrument for calibration with bootstrapping
            after each calibrator run */
		virtual void reset();

    protected:
        /** for reflection */
        ObjFuncCombo();

    private:
        static void load(CClassSP& clazz);
        static IObject* defaultCtor();

        // registered
        ObjFuncArray objFuncs;
        DoubleArray  weights;

        // transient
        int           nbFuncs;
        ObjectArray   adjGroupArray;
    };

    /** Linear combination of a set of ObjFuncLeastSquare's */
    class RISKMGR_DLL ObjFuncLeastSquareCombo: public ObjFuncLeastSquare,
                                   public virtual ObjFunc::IBootstrappable{
    public:
        static CClassConstSP const TYPE;

        void validatePop2Object();

        /** Give the derived class a chance to do s'thing with the market */
        virtual void getMarket(MarketData* market);

        /** Give the derived class a chance to do additional validation
            irrespective of getMarket being called or not */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup();

        /** Returns the number of functions */
        virtual int getNbFuncs() const;

		/** Returns the array of objective functions */
        virtual ObjFuncLeastSquareArray getObjFuncArray();

        /** Calculates the values of the objective functions */
        virtual void calcValue(CDoubleArray& funcvals) const;

        /** Makes additional validations for calibration with bootstrapping */
        virtual void validate(const InstanceID::IBootstrappableArray& ids) const;

        /** Updates the instrument for calibration with bootstrapping
            before each calibrator run */
		virtual void update(int idxMat);

        /** Reset the instrument for calibration with bootstrapping
        after each calibrator run */
		virtual void reset();

    protected:
        /** for reflection */
        ObjFuncLeastSquareCombo();

    private:
        static void load(CClassSP& clazz);
        static IObject* defaultCtor();

        // registered
        ObjFuncLeastSquareArray  objFuncs;
        DoubleArray              weights;

        // transient
        int                      nbFuncs;
        int                      totalNbFuncs;
        mutable DoubleArrayArray vals;
        ObjectArray              adjGroupArray;
    };
	typedef smartPtr<ObjFuncLeastSquareCombo> ObjFuncLeastSquareComboSP;

    /**  PenaltyFunc penalizes a set of variables for moving away
         from a set of initial values too much. */
    class RISKMGR_DLL PenaltyFunc: public ObjFuncLeastSquare{
    public:
        static CClassConstSP const TYPE;
        friend class Calibrator_PenaltyFuncHelper;

        PenaltyFunc();

        virtual void validatePop2Object();

        /** Gives this class the chance to use the market */
        virtual void getMarket(MarketData* market);

        /** Gives this class a chance to do additional validation
            after getMarket has been called */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup();

        /** Returns the number of functions */
        virtual int getNbFuncs() const;

        /** Calculates the values of the objective functions */
        virtual void calcValue(CDoubleArray& funcvals) const;

    private:

        /* Registered fields */
        InstanceIDArray     ids;
        DoubleArray         weights;
        IObjectSP           object;
        string              objectType;     // needed if object is a market object

        /* Transient fields */
        int                 nbVars;
        DoubleArray         initvals;
        mutable DoubleArray currvals;
        DoubleArray         usedweights;
    };

    // InstanceID for doubles
    class RISKMGR_DLL InstanceIDDb: public InstanceID {
    private:
        bool               useOverride;
        double             override; // overrides initial value in object
        int                arrayIndex;// if field is an array, the elt to calibrate
		RangeSP			   rangeOverride;	// optional override for Range in object
        bool               skipTransient;   // optional 
        class ReadWriteVal;
        friend class ReadWriteVal;

        InstanceIDDb();

        static IObject* defaultCtor();

        static void load(CClassSP& clazz);

    public:
        static CClassConstSP const TYPE;

        virtual void validatePop2Object();

        InstanceIDDb(const string&  typeToCalibrate,
                     const string&  name,
                     const string&  field,
                     bool           useOverride,
                     double         override,
                     int            arrayIndex);

		/** override constructor to allow for rangeOverride */
		InstanceIDDb(const string&  typeToCalibrate,
			const string&  name,
			const string&  field,
			bool           useOverride,
			double         override,
			int            arrayIndex,
			RangeSP		   rangeOverride);

        /** override constructor to allow for rangeOverride */
		InstanceIDDb(const string&  typeToCalibrate,
			const string&  name,
			const string&  field,
			bool           useOverride,
			double         override,
			int            arrayIndex,
			bool		   skipTransient);

        /** Must be called before any other methods - allows object to see the
            market data (eg length of arrays) */
        virtual void initialise(const IObjectSP&  adjGroup);

        /** Returns the number of variables */
        virtual int numVariables() const;

        /** Populates the names array from 'index' to
            'index' + numVariables() */
        virtual void getNames(StringArray& names, int index) const;

        /** Populates the ranges array from 'index' to
            'index' + numVariables() */
        virtual void getRanges(RangeArray& ranges, int index) const;

        /** Returns true if scalar type; false if vector type */
        virtual bool isScalar() const;

        /** If vector type, returns indices of elts that got calibrated.
            Returns 0, otherwise. */
        virtual IObjectSP getArrayIndices() const;

        /** Populates the values array from 'index' to
            'index' + numVariables() */
        virtual void getInitialValues(
            DoubleArray&         values,
            int                  index,
            const IObjectSP&  adjGroup) const;

        /** Populates the values array from 'index' to
            'index' + numVariables() */
        virtual void getValues(
                DoubleArray&         values,
                int                  index,
                const IObjectSP&  adjGroup) const;

        /** Uses the values in x from 'index' to
            'index' + numVariables() to adjust the field being calibrated */
        void applyAdjustment(const IObjectSP&  adjGroup,
                             const CDoubleArray&  x,
                             int                  index) const;

	    /** Return the calibratedVals[index] to
            calibratedVals[index+numVariables] as an object */
        IObjectSP writeResults(const DoubleArray& calibratedVals,
                               int                index) const;

        /** Read off the 'values' object and populate calibratedVars array */
        virtual void readResults(const IObject& values,
                                 DoubleArray&   calibratedVals,
                                 int            index) const;
    };

    /** Run calibrator */
    CResultsSP run(ObjFunc&               objFunc,
                   const InstanceIDArray& instanceIDs) const;

    CResultsSP run(
        ObjFunc&               objFunc,
        const InstanceIDArray& instanceIDs,
        DoubleArray xguess,
        double objFuncVal) const;


    /** method to extract values from results used in bootstrapping */
    static void getResultValues(
        CResultsSP currentResults,  // (I)
        DoubleArray & aggregCalibValues // (O)
        );

    //// Constructor
    Calibrator(const OptimizerNDSP& optimizer);

    //// Utility: Returns all distinct names of objects of specified type
    static StringArray getNames(CClassConstSP  clazz,
                                IObjectConstSP adjGroup);

private:
    Calibrator();
    Calibrator(const Calibrator &rhs);
    Calibrator& operator=(const Calibrator& rhs);
    class NameCollector;

    /* Registered vars */
    OptimizerNDSP       optimizer;
};

typedef smartConstPtr<Calibrator> CalibratorConstSP;
typedef smartPtr<Calibrator> CalibratorSP;

DRLIB_END_NAMESPACE
#endif

