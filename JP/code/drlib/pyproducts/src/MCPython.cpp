//----------------------------------------------------------------------------
//
//   Group       : Global QR
//
//   Filename    : MCPython.cpp
//   Date        : 1 Nov 2006
//   Description : python flexible payoff for Monte Carlo engine
//----------------------------------------------------------------------------

#include "edginc/PyInterface.hpp"

DRLIB_BEGIN_NAMESPACE

void MCPythonProd::payoff(const IPathGenerator* pathGenerator, IMCPrices& prices) {
    if (useSrcCode){ // add ! when we actually store code text, only for testing performance
        int    iAsset;
        double payoff;
        // Any past samples
        sumOut = sumSoFar;
        
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            const SVPath& path = spotSV->path(iAsset);
            for (int iStep=path.begin(); iStep<path.end(); iStep++) {
                sumOut[iAsset] += path[iStep];
            }
        }
        if (spotSV->doingPast()){ // preserve values from past calc
            sumSoFar = sumOut;
        }

        // Continue using sumOut as a working area - now for perfs
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            // calc average
            sumOut[iAsset] /= refLevelSV->refLevel(iAsset) * inst->averageOutDates.size();
        }
        
        // Finding best (or worst) does not need a sort - just a
        // single pass thru
        int    pickedAsset = 0;

        if (inst->isBest) {
            for (iAsset=1; iAsset<nbAssets; iAsset++) {
                if (sumOut[iAsset] > 
                    sumOut[pickedAsset]) {
                    pickedAsset = iAsset;
                } 
            }
        } else {
            for (iAsset=1; iAsset<nbAssets; iAsset++) {
                if (sumOut[iAsset] < 
                    sumOut[pickedAsset]) {
                    pickedAsset = iAsset;
                } 
            }
        }

        payoff = CoP*(sumOut[pickedAsset] - inst->strike);

        // cout << "C++ payoff = " << payoff << endl;

        payoff = inst->notional * prices.maxWithZero(payoff);

        prices.add(payoff); 
    }
    else{
        try {
            pyInterfaceInit();
            boost::python::object module( boost::python::handle<>( PyImport_ImportModule( "MCPython" ) ) );
            boost::python::object payoffFunction = module.attr( "payoff" );
            MCPythonProdPy productPy( *this );
            IPathGeneratorPy pathGeneratorPy( *pathGenerator );
            IMCPricesPy pricesPy( prices );
            payoffFunction( productPy, pathGeneratorPy, pricesPy );
        } catch ( boost::python::error_already_set ) {
            PyErr_Print();
        }
    }
}

void MCPythonProd::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    static const string routine = "MCPython::pathGenUpdated";
    try {
        spotSV = spotGen->getSpotSV(newPathGen);
        refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
        dfSV = dfGen->getSVDiscFactor(newPathGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* MCPython::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new MCPythonProd(this, simSeries);
}

CClassConstSP const MCPython::TYPE = CClass::registerClassLoadMethod(
    "MCPython", typeid(MCPython), MCPython::load);

bool PYPRODUCTS_DLL MCPythonLoad() {
    return (MCPython::TYPE != 0);
}

DRLIB_END_NAMESPACE