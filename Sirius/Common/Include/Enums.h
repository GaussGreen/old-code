///////////////////////////////////////////////////////////////////////////////
//	Enums.h : Contains all enumerator values that are common across the whole
//		      of Sirius
//

#ifndef _ENUMS_H
#define _ENUMS_H

#if !defined(__midl)
#pragma once
#endif

typedef enum DataSourceEnum
// NOTES: 1) don't make the enums any longer than 16 characters or stored proedures
//			 (e.g. sp_user_get_zerocurve) will break.
//		  2) Pick the data source names carefully so they will not clash with location
//			 codes. (Or, for example, the Location / DataSource specification in the 
//			 volatility structure put_Value will become ambiguous).
{
	NoDataSource = 0,
	Last = 1,
	PL = 2	
} DataSourceEnum;

typedef enum PerformanceTypeEnum	// for Variance Swap
{ 
	NoPerformanceType = 0,
	LogarithmicType = 1,
	StraddleType = 2
} PerformanceTypeEnum;

typedef enum SpotOrForwardEnum
{
	NoSpotOrForward = 0,
	Spot = 1,
	Forward = 2
} SpotOrForwardEnum;

typedef enum DividendTypeEnum
{
	NoDividendType = 0,
	Discrete = 1,
	Continuous = 2
} DividendTypeEnum;

typedef enum FileSystemEnum 
{
	fsUnknown = 0,
	fsSQLServer = 1,
	fsNTFS = 2
} FileSystemEnum;
	
typedef enum DatabaseModeEnum
{
	NoDatabaseConnection = 0,
	DirectDatabaseConnection = 1,
	UseSiriusServer = 2,
	RetrieveSpotSchedulesOnly = 3,
	DirectExceptPositions = 4,
	UseTaurus = 5,							// Don't change this value since BwlData relies on it
	TaurusExceptPositions = 6
} DatabaseModeEnum;

typedef enum PrincipalExchangedTypeEnum
{
	NoPrincipalExchange = 0,
	Both = 1,
	Final = 2,
	Initial = 3
} PrincipalExchangedTypeEnum;

typedef enum NotionalTypeEnum
{
	NoNotionalType = 0,
	ConstantNotional = 1,
	ConstantShare = 2,
	RAMConstantNotional = 3
} NotionalTypeEnum;

typedef enum SwapLegEnum
{
	NoSwapLeg = 0,
	A_Leg = 1,
	B_Leg = 2
} SwapLegEnum;

typedef enum SwapLegTypeEnum
{
	NoSwapLegType = 0,
	FixedLeg = 1,
	FloatingLeg = 2
} SwapLegTypeEnum;

typedef enum YesNoEnum
{
	No = 0,
	Yes = 1	
} YesNoEnum;

typedef enum MinMaxEnum
{
	Minimum = 0,
	Maximum = 1
} MinMaxEnum;

typedef enum PayoffTypeEnum
{
	NoPayoffType = 0,
	Call = 1,
	Put = 2,
	Straddle = 3
} PayoffTypeEnum;

typedef enum RollDateConventionEnum
{
	NoRollDateConvention = 0,
	EOM = 1,
	FRN = 2,
	IMM = 3,
	Standard = 4
} RollDateConventionEnum;

typedef enum StubTypeEnum
{
	NoStubType = 0,
	Auto = 1,
	Long = 2,
	Short = 2
} StubTypeEnum;

typedef enum InterpolateTypeEnum
{
	Volatility = 0,
	Variance = 1,
	AsVolData = 2
} InterpolateTypeEnum;

typedef enum ZeroCurveTypeEnum
{
	Natural			= 0,
	Pay				= 1,
	Composite		= 2
} ZeroCurveTypeEnum;

typedef enum AssetTypeEnum
{
	InvalidAsset = 0,	
	CurrencyAsset = 1,
	PrototypeNaturalAsset = 2,	// This is a natural asset without currencies defined
	NaturalAsset = 3,
	SingleUnderlying = 4,
	QuantoBasket = 5,
	CompositeBasket = 6
} AssetTypeEnum;

typedef enum DayCountConventionEnum
{	
	NoDayCountConvention	= 0,
	ActualActual			= 1,
	Act365					= 2,
	Act360					= 3,	
	Thirty360				= 4,
	Act365F					= 5,
	ActActISDA				= 6,
	Actual365				= 7,
	Actual360				= 8,
	Actual365F				= 9,
	ActualActualISDA		= 10
} DayCountConventionEnum;

typedef enum BusinessDayConventionEnum
{	
	NoBusinessDayConvention = 0,
	ModifiedFollowing	= 1,
	Following			= 2,
	ModifiedPreceding	= 3,
	Preceding			= 4,
	NoChange			= 5
} BusinessDayConventionEnum;

typedef enum WeekdayEnum
{ 
	Sunday    = 1,
	Monday    = 2,
	Tuesday   = 3,
	Wednesday = 4,
	Thursday  = 5,
	Friday    = 6,
	Saturday  = 7
} WeekdayEnum;

typedef enum MonthEnum
{ 
	January   = 1,
	February  = 2,
	March     = 3,
	April     = 4,
	May       = 5,
	June      = 6,
	July      = 7,
	August    = 8,
	September = 9,
	October   = 10,
	November  = 11,
	December  = 12 
} MonthEnum;

typedef enum TimeUnitEnum
{ 
	Days   = 0,
	Weeks  = 1,
	Months = 2,
	Years  = 3 
} TimeUnitEnum;

typedef enum InterpolatorTypeEnum
{
	NoInterpolatorType		= 0,
	Linear					= 1,
	CubicSpline				= 2,
	MonotonicSpline			= 3,
	Constant				= 4,
	TwoDimensional			= 5,
	FitPolynomial			= 6,
	FitPolynomialExlicit	= 7
} InterpolatorTypeEnum;



typedef enum StrikesTypeEnum
{
	NoStrikesType	= 0,
	Fixed			= 1,		
	Normalised		= 2,
	ForwardBased	= 3,
	SpotBased		= 4,
	LogBased		= 5	
} StrikesTypeEnum;

typedef enum MonteCarloTypeEnum
{
	NoMonteCarloType			= 0,
	ForwardSkew					= 1,
	General						= 2,
	HybridForwardSkew			= 3,
	Quasi						= 6,
	Hermite						= 7,
	LocalVolMonteCarlo			= 8
} MonteCarloTypeEnum;

typedef enum PdeTypeEnum
{
	NoPdeType					= 0,
	BlackScholes				= 1,
	LocalVol					= 2
} PdeTypeEnum;

typedef enum BidAskMidEnum
{
	Bid				= 1,
	Ask				= 2,
	Middle			= 3			// Don't change this value since BwlData relies on it, cannot have Mid since it's a VB keywork
} BidAskMidEnum;

typedef enum VolatilityDataTypeEnum
{
	TermVolatilities = 1,
	ForwardVolatilities = 2
} VolatilityDataTypeEnum;


typedef enum CurrencyEnum 
{
    NoCurrency = 0,
	ARS = 1,    //!< Argentinian Peso
    ATS = 2,    //!< Austrian Schillings
    AUD = 3,    //!< Australian Dollar
    BDT = 4,    //!< Bangladesh Taka
    BEF = 5,    //!< Belgian Franc
    BGL = 6,    //!< Bulgarian Lev
    BRL = 7,    //!< Brazilian Real
    BYB = 8,    //!< Belarusian Ruble
    CAD = 9,    //!< Canadian Dollar
    CHF = 10,    //!< Swiss Franc
    CLP = 11,    //!< Chilean Peso
    CNY = 12,    //!< Chinese Yuan
    COP = 13,    //!< Colombian Peso
    CYP = 14,    //!< Cyprus Pound
    CZK = 15,    //!< Czech Koruna
    DEM = 16,    //!< German Mark
    DKK = 17,    //!< Danish Krone
    EEK = 18,    //!< Estonian Kroon
    EUR = 19,    //!< Euro
    GBP = 20,    //!< British Pound
    GRD = 21,    //!< Greek Drachma
    HKD = 22,    //!< Hong Kong Dollar
    HUF = 23,    //!< Hungarian Forint
    ILS = 24,    //!< Israeli Shekel
    INR = 25,    //!< Indian Rupee
    IQD = 26,    //!< Iraqi Dinar
    IRR = 27,    //!< Iranian Rial
    ISK = 28,    //!< Iceland Krona
    ITL = 29,    //!< Italian Lira
    JPY = 30,    //!< Japanese Yen
    KRW = 31,    //!< South-Korean Won
    KWD = 32,    //!< Kuwaiti dinar
    LTL = 33,    //!< Lithuanian Litas
    LVL = 34,    //!< Latvian Lats
    MTL = 35,    //!< Maltese Lira
    MXP = 36,    //!< Mexican Peso
    NOK = 37,    //!< Norwegian Kroner
    NPR = 38,    //!< Nepal Rupee
    NZD = 39,    //!< New Zealand Dollar
    PKR = 40,    //!< Pakistani Rupee
    PLN = 41,    //!< New Polish Zloty
    ROL = 42,    //!< Romanian Leu
    SAR = 43,    //!< Saudi Riyal
    SEK = 44,    //!< Swedish Krona
    SGD = 45,    //!< Singapore Dollar
    SIT = 46,    //!< Slovenian Tolar
    SKK = 47,    //!< Slovak Koruna
    THB = 48,    //!< Thai Baht
    TRL = 49,    //!< Turkish Lira
    TTD = 50,    //!< Trinidad & Tobago dollar
    TWD = 51,    //!< Taiwan Dollar
    USD = 52,    //!< US Dollar
    VEB = 53,    //!< Venezuelan Bolivar
    ZAR = 54     //!< South African Rand
} CurrencyEnum;

#endif
