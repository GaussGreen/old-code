#include "../credit/src/pch-credit.cxx"

#include "../convolution/edginc/ABSCDOConvolution.hpp"
#include "../convolution/edginc/BaseCorrelationLossCalculator.hpp"
#include "../convolution/edginc/BCConvolutionModelConfig.hpp"
#include "../convolution/edginc/BCPortfolioNameConvolutionLossGen.hpp"
#include "../convolution/edginc/BCTrancheConvolutionLossGen.hpp"
#include "../convolution/edginc/BetaConvolutor.hpp"
#include "../convolution/edginc/CCMCalibration.hpp"
#include "../convolution/edginc/CCMConvolution.hpp"
#include "../convolution/edginc/CCMFastLossCalculator.hpp"
#include "../convolution/edginc/CCMLossCalculator.hpp"
#include "../convolution/edginc/CCMLossCalculatorBase.hpp"
#include "../convolution/edginc/CCMLossUnit.hpp"
#include "../convolution/edginc/CCMRFLFastLossCalculator.hpp"
#include "../convolution/edginc/CCMRFLLossCalculator.hpp"
#include "../convolution/edginc/CCMSkew.hpp"
#include "../convolution/edginc/CCMTrancheCalculatorLegacy.hpp"
#include "../convolution/edginc/CCMTrancheFastLossCalculatorLegacy.hpp"
#include "../convolution/edginc/CCMTrancheUtils.hpp"
#include "../convolution/edginc/ConvolutionEngine.hpp"
#include "../convolution/edginc/ConvolutionLib.hpp"
#include "../convolution/edginc/ConvolutionModel.hpp"
#include "../convolution/edginc/ConvolutionModelConfig.hpp"
#include "../convolution/edginc/ConvolutionProduct.hpp"
//#include "../convolution/edginc/CounterPartyCredit.hpp"
#include "../convolution/edginc/CreditMetricsABSCDOLossCalculator.hpp"
#include "../convolution/edginc/CreditMetricsDefaultsModel.hpp"
#include "../convolution/edginc/CreditMetricsFastLossCalculator.hpp"
#include "../convolution/edginc/CreditMetricsLossCalculator.hpp"
#include "../convolution/edginc/CreditMetricsLossCalculatorBase.hpp"
#include "../convolution/edginc/CreditMetricsRFLFastLossCalculator.hpp"
#include "../convolution/edginc/CreditMetricsRFLLossCalculator.hpp"
#include "../convolution/edginc/FixedTrancheLossCalculator.hpp"
#include "../convolution/edginc/GaussianMarketFactorModel.hpp"
#include "../convolution/edginc/IBCCondLossDistributionsGen.hpp"
#include "../convolution/edginc/IConditionalDefaultsModel.hpp"
#include "../convolution/edginc/ICondLossDistributionsGen.hpp"
#include "../convolution/edginc/IConvolutor.hpp"
#include "../convolution/edginc/ILossDistributionsGen.hpp"
#include "../convolution/edginc/IMarketFactorModel.hpp"
#include "../convolution/edginc/IMarketFactorValue.hpp"
#include "../convolution/edginc/LondonFloorLossCalculator.hpp"
#include "../convolution/edginc/LossDistribution.hpp"
//#include "../convolution/edginc/MarketFactorValue.hpp"
#include "../convolution/edginc/PortfolioNameConvolutionLossGen.hpp"
#include "../convolution/edginc/RecursiveConvolutor.hpp"
#include "../convolution/edginc/SaddlePoint.hpp"
#include "../convolution/edginc/TrancheConvolutionLossGen.hpp"
#include "../convolution/edginc/TrancheLossCalculator.hpp"
#include "../convolution/edginc/TrancheLossCalculatorLegacy.hpp"
