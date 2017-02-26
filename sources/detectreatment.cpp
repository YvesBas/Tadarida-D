#include "detec.h"
#include "detec.h"
class Detec;

// ParamToSave: the settings of ".ta" files are described by objects of this class
ParamToSave::ParamToSave(int arrayNumber,int paramaterNumber,QString columnTitle)
{
    ArrayNumber=arrayNumber;
    ParameterNumber=paramaterNumber;
    ColumnTitle=columnTitle;
    FromVersion = 0;
	ToVersion = -1;
}

ParamToSave::ParamToSave(int arrayNumber,int parameterNumber,QString columnTitle,int needVer)
{
    ArrayNumber=arrayNumber;
    ParameterNumber=parameterNumber;
    ColumnTitle=columnTitle;
    FromVersion = needVer;
	ToVersion = -1;
}

ParamToSave::ParamToSave(int arrayNumber,int parameterNumber,QString columnTitle,int needVer,int limVer)
{
    ArrayNumber=arrayNumber;
    ParameterNumber=parameterNumber;
    ColumnTitle=columnTitle;
    FromVersion = needVer;
	ToVersion = limVer;
}

ParamToSave::ParamToSave()
{
}

ParamToSave::~ParamToSave()
{
}

// each thread of detec class creates an object of DetecTreatment class
// which supports the processing of sound files
DetecTreatment::DetecTreatment(Detec *pDet)
{
    _detec = pDet;
    _detec->LogStream << "detec->IThread=" << _detec->IThread << endl;
    // global variables allowing the computefft method to be thread safe
    _complexInput = _detec->PMainWindow->ComplexInput[_detec->IThread];
    _fftRes = _detec->PMainWindow->FftRes[_detec->IThread];
    _resultSuffix = QString("ta");
    if(LINWIN==1) _resultCompressedSuffix = QString("tac");
    else _resultCompressedSuffix = QString("ta.gz");
    // version of the list of settings of sound events in the ".ta" files
    _paramVersion = 1;
    // initializes the definitions of settings
    initVectorParams();
}

DetecTreatment::DetecTreatment()
{
    _paramVersion = 1;
    initVectorParams();
}

DetecTreatment::~DetecTreatment()
{
}

void DetecTreatment::InitializeDetecTreatment()
{
    // this method creates the variables and arrays used by treatments of sound files
    _data				= ( float* ) fftwf_malloc( sizeof( float ) * FFT_HEIGHT_MAX );
    _coeff = new double[FFT_HEIGHT_HALF_MAX];
    SonogramArray = new qint16*[MAXHEIGHT];
    PointFlagsArray = new char *[MAXHEIGHT];
    _charSonogramArray = new char[MAXHEIGHT*SONOGRAM_WIDTH_MAX*sizeof(qint16)];
    _charPointFlagsArray = new char[MAXHEIGHT*LD8];
    for ( int i=0 ; i <MAXHEIGHT ; i++)
    {
        SonogramArray[i] = (qint16 *)((char *)(_charSonogramArray+(i*SONOGRAM_WIDTH_MAX*sizeof(qint16))));
        PointFlagsArray[i] = (char *)(_charPointFlagsArray+(i*LD8));
        memset(SonogramArray[i],0,SONOGRAM_WIDTH_MAX*2);
    }
    EnergyMin = 0.0f;
    _charParamsArray = new char[MAXCRI*NBTAB*NBPAR*sizeof(float)];
    _paramsArray = new float**[MAXCRI];
    _numberCallParameters = VectPar.size();
    for (int i=0;i < MAXCRI; i++)
    {
        _paramsArray[i] = new float *[NBTAB];
        for(int j=0;j < NBTAB; j++)
            _paramsArray[i][j] = (float *)((char *)(_charParamsArray+(i*NBTAB*NBPAR+j*NBPAR)*sizeof(float)));
    }
    LowSlope = new int[NCRETES*MAXCRI*2];
    Inflexion1 = new int[NCRETES*MAXCRI*2];
    Inflexion3 = new int[NCRETES*MAXCRI*2];
    _harmonic=new int *[2];
    _dpm=new float *[2];
    _dypm=new int *[2];
    _ypm=new int *[2];
    for(int k=0;k<2;k++)
    {
        _harmonic[k]=new int [MAXCRI];
        _dpm[k]=new float[MAXCRI];
        _dypm[k]=new int[MAXCRI];
        _ypm[k]=new int[MAXCRI];
    }
    _tabY = new float[MAXHEIGHT];
    _numberPixelsPerY = new int[MAXHEIGHT];
    _numberPixelsPerX = new int[MAXLARCRI];
    _averagePerX = new float[MAXLARCRI];
    _xMinPerY = new quint16[MAXHEIGHT];
    _xSecondWestRidgePerY = new quint16[MAXHEIGHT];
    _xMaxPerY = new int[MAXHEIGHT];
    _yMinPerX = new quint16[MAXLARCRI];
    _yMaxPerX = new quint16[MAXLARCRI];
    _yEbarPerX = new float[MAXLARCRI];
    _totYEPerX = new float[MAXLARCRI];
    _eMaxPerX = new float[MAXLARCRI];
    int maxOfTwo = qMax(MAXLARCRI,MAXHEIGHT);
    _slope = new float[maxOfTwo];
    _charTabX = new char[MAXCRI*MAXLARCRI*sizeof(float)];
    _tabX = new float*[MAXCRI];
    for(int i=0 ; i < MAXCRI ; i++) _tabX[i] = (float *)((char *)(_charTabX+(i*MAXLARCRI*sizeof(float))));
    _charYEmaxPerX = new char[MAXCRI*MAXLARCRI*sizeof(qint16)];
    _yEmaxPerX = new quint16*[MAXCRI];
    for(int i=0;i<MAXCRI;i++)  _yEmaxPerX[i] = (quint16 *)((char *)(_charYEmaxPerX+(i*MAXLARCRI*sizeof(quint16))));
    _charTabYX = new char[MAXHEIGHT*MAXLARCRI*sizeof(float)];
    _tabYX = new float*[MAXHEIGHT];
    for(int i=0;i<MAXHEIGHT;i++) _tabYX[i] = (float *)((char *)(_charTabYX+(i*MAXLARCRI*sizeof(float))));
    EnergyColumAverage = new float[SONOGRAM_WIDTH_MAX];
    FlagGoodCol = new bool[SONOGRAM_WIDTH_MAX];
    FlagGoodColInitial = new bool[SONOGRAM_WIDTH_MAX];
    _sortMp = new int[MAXCRI];
    _invMp = new int[MAXCRI];
    _xMp = new int[MAXCRI];
    _detec->LogStream << "dtt modefreq=" << _modeFreq << endl;
}

void DetecTreatment::SetDirParameters(QString wavPath,QString txtPath,bool imageData,
                                                                       QString imagePath,QString datPath)
{
    // this method affects the paths of the directories of input and output files
    _wavPath = wavPath;
    _txtPath = txtPath;
    _imageData = imageData;
    _imagePath = imagePath;
    _datPath = datPath;
}

// assignment of settings for processing files
void DetecTreatment::SetGlobalParameters(int modeFreq,int timeExpansionLeft,int timeExpansionRight,int detectionThreshold,int stopThreshold,
                                 int freqMin,int nbo,bool useValflag,
                                int jumpThreshold,int widthBigControl,int widthLittleControl,
                                int highThresholdJB,int lowThresholdJB,int lowThresholdC,int highThresholdC,int qR,int qN,int parVer,bool desactiveCorrectNoise)
{
    // multifreq :
    _modeFreq = modeFreq;
    _timeExpansionLeft = timeExpansionLeft;
    _timeExpansionRight = timeExpansionRight;
    _detectionThreshold = detectionThreshold;
    _stopThreshold = stopThreshold;
    _freqMin = freqMin;
    _nbo = nbo;
    _useValflag = useValflag;
    _jumpThreshold = jumpThreshold;
    _widthBigControl = widthBigControl;
    _widthLittleControl = widthLittleControl;
    _highThresholdJB = highThresholdJB;
    _lowThresholdJB = lowThresholdJB;
    _lowThresholdC = lowThresholdC;
    _highThresholdC = highThresholdC;
    _qR = qR;
    _qN = qN;
    _paramVersion = parVer;
    //
    if(_modeFreq==1) _freqCallMin = 8.0f;
   else _freqCallMin=0.8f;
    _desactiveCorrectNoise = desactiveCorrectNoise;
}

// this method defines the settings which will be computed by the method detectsparameter2
// they are recorded in a vector of objects of the class: ParamToSave
// according to the values of parameters needVer and limVer, the settings will be recorded
// in the ".ta" file
void DetecTreatment::initVectorParams()
{

    QString prefix[] = {"","CM_","CS_","CN_","CO_","CO2_"};
    VectPar.push_back(ParamToSave(SH,StTime,"StTime"));
    VectPar.push_back(ParamToSave(SH,Dur,"Dur"));
    VectPar.push_back(ParamToSave(SH,PrevSt,"PrevSt"));
    VectPar.push_back(ParamToSave(SH,Fmax,"Fmax",0,1));
    VectPar.push_back(ParamToSave(SH,Fmin,"Fmin"));
    VectPar.push_back(ParamToSave(SH,BW,"BW"));
    VectPar.push_back(ParamToSave(SH,FreqMP,"FreqMP",0,1));
    VectPar.push_back(ParamToSave(SH,PosMP,"PosMP"));
    VectPar.push_back(ParamToSave(SH,FreqPkS,"FreqPkS",0,1));
    VectPar.push_back(ParamToSave(SH,FreqPkM,"FreqPkM",0,1));
    VectPar.push_back(ParamToSave(SH,PosPkS,"PosPkS",0,1));
    VectPar.push_back(ParamToSave(SH,PosPkM,"PosPkM",0,1));
    VectPar.push_back(ParamToSave(SH,FreqPkS2,"FreqPkS2",0,1));
    VectPar.push_back(ParamToSave(SH,FreqPkM2,"FreqPkM2",0,1));
    VectPar.push_back(ParamToSave(SH,PrevMP1,"PrevMP1"));
    VectPar.push_back(ParamToSave(SH,PrevMP2,"PrevMP2"));
    VectPar.push_back(ParamToSave(SH,NextMP1,"NextMP1"));
    VectPar.push_back(ParamToSave(SH,NextMP2,"NextMP2"));
    VectPar.push_back(ParamToSave(SH,Amp1,"Amp1"));
    VectPar.push_back(ParamToSave(SH,Amp2,"Amp2"));
    VectPar.push_back(ParamToSave(SH,Amp3,"Amp3"));
    VectPar.push_back(ParamToSave(SH,Amp4,"Amp4"));
    VectPar.push_back(ParamToSave(SH,NoisePrev,"NoisePrev"));
    VectPar.push_back(ParamToSave(SH,NoiseNext,"NoiseNext"));
    VectPar.push_back(ParamToSave(SH,NoiseDown,"NoiseDown"));
    VectPar.push_back(ParamToSave(SH,NoiseUp,"NoiseUp"));
    VectPar.push_back(ParamToSave(SH,CVAmp,"CVAmp"));
    VectPar.push_back(ParamToSave(CO,Dur,"CO_Dur",0,1));
    VectPar.push_back(ParamToSave(CO2,Dur,"CO2_Dur",0,1));
    VectPar.push_back(ParamToSave(CM,Fmax,"CM_Fmax",0,1));
    VectPar.push_back(ParamToSave(CS,Fmax,"CS_Fmax",0,1));
    VectPar.push_back(ParamToSave(CM,Fmin,"CM_Fmin",0,1));
    VectPar.push_back(ParamToSave(CN,Fmin,"CN_Fmin",0,1));
    VectPar.push_back(ParamToSave(CM,BW,"CM_BW",0,1));
    VectPar.push_back(ParamToSave(CS,BW,"CS_BW",0,1));
    VectPar.push_back(ParamToSave(CN,BW,"CN_BW",0,1));
    VectPar.push_back(ParamToSave(CO2,FPk,"CO2_FPk",0,1));
    VectPar.push_back(ParamToSave(CO2,FPkD,"CO2_FPkD"));
    VectPar.push_back(ParamToSave(CO2,TPk,"CO2_TPk"));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,Slope,prefix[i]+"Slope"));
    for(int i=CO;i<=CO2;i++)
        if(i==CO) VectPar.push_back(ParamToSave(i,ISlope,prefix[i]+"ISlope",1,1));
       else VectPar.push_back(ParamToSave(i,ISlope,prefix[i]+"ISlope",1));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,HCF,prefix[i]+"HCF",0,1));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,THCF,prefix[i]+"THCF"));
    for(int i=CM;i<=CO2;i++)
        if(i>=CO) VectPar.push_back(ParamToSave(i,FIF,prefix[i]+"FIF",0,1));
        else VectPar.push_back(ParamToSave(i,FIF,prefix[i]+"FIF"));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,LCF,prefix[i]+"LCF",0,1));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,UpSl,prefix[i]+"UpSl"));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,LoSl,prefix[i]+"LoSl"));
    VectPar.push_back(ParamToSave(CM,StF,"CM_StF",0,1));
    VectPar.push_back(ParamToSave(CM,EnF,"CM_EnF",0,1));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,StSl,prefix[i]+"StSl"));
    for(int i=CM;i<=CO2;i++) VectPar.push_back(ParamToSave(i,EnSl,prefix[i]+"EnSl"));
    for(int i=CM;i<=CO2;i++)
        if(i==CM || i==CO2) VectPar.push_back(ParamToSave(i,FISl,prefix[i]+"FPSl",0,1));
        else VectPar.push_back(ParamToSave(i,FISl,prefix[i]+"FPSl"));
    VectPar.push_back(ParamToSave(CM,FISl,"CM_FISl"));
    VectPar.push_back(ParamToSave(CO2,FISl,"CO2_FISl"));
    for(int i=CM;i<=CN;i++) VectPar.push_back(ParamToSave(i,CeF,prefix[i]+"CeF",0,1));
    VectPar.push_back(ParamToSave(CM,B5dBBF,"CM_5dBBF",0,1));
    VectPar.push_back(ParamToSave(CM,B5dBAF,"CM_5dBAF",0,1));
    VectPar.push_back(ParamToSave(CM,B5dBBW,"CM_5dBBW"));
    VectPar.push_back(ParamToSave(CM,B5dBDur,"CM_5dBDur"));
    VectPar.push_back(ParamToSave(CO2,B5dBBF,"CO2_5dBBF",0,1));
    VectPar.push_back(ParamToSave(CO2,B5dBAF,"CO2_5dBAF",0,1));
    VectPar.push_back(ParamToSave(CO2,B5dBBW,"CO2_5dBBW"));
    VectPar.push_back(ParamToSave(CO2,B5dBDur,"CO2_5dBDur"));
    VectPar.push_back(ParamToSave(SH,Hup_RFMP,"Hup_RFMP"));
    VectPar.push_back(ParamToSave(SH,Hup_PosMP,"Hup_PosMP",0,1));
    VectPar.push_back(ParamToSave(SH,Hup_PosSt,"Hup_PosSt",0,1));
    VectPar.push_back(ParamToSave(SH,Hup_PosEn,"Hup_PosEn",0,1));
    VectPar.push_back(ParamToSave(SH,Hup_AmpDif,"Hup_AmpDif"));
    VectPar.push_back(ParamToSave(SH,Hup_RSlope,"Hup_RSlope",0,1));
    VectPar.push_back(ParamToSave(SH,Hlo_RFMP,"Hlo_RFMP",0,1));
    VectPar.push_back(ParamToSave(SH,Hlo_PosMP,"Hlo_PosMP",0,1));
    VectPar.push_back(ParamToSave(SH,Hlo_PosSt,"Hlo_PosSt",0,1));
    VectPar.push_back(ParamToSave(SH,Hlo_PosEn,"Hlo_PosEn"));
    VectPar.push_back(ParamToSave(SH,Hlo_AmpDif,"Hlo_AmpDif"));
    VectPar.push_back(ParamToSave(SH,Hlo_RSlope,"Hlo_RSlope",0,1));
    VectPar.push_back(ParamToSave(SH,Ramp_2_1,"Ramp_2_1"));
    VectPar.push_back(ParamToSave(SH,Ramp_3_1,"Ramp_3_1"));
    VectPar.push_back(ParamToSave(SH,Ramp_3_2,"Ramp_3_2"));
    VectPar.push_back(ParamToSave(SH,Ramp_1_2,"Ramp_1_2"));
    VectPar.push_back(ParamToSave(SH,Ramp_4_3,"Ramp_4_3",0,1));
    VectPar.push_back(ParamToSave(SH,Ramp_2_3,"Ramp_2_3"));
    VectPar.push_back(ParamToSave(SH,RAN_2_1,"RAN_2_1"));
    VectPar.push_back(ParamToSave(SH,RAN_3_1,"RAN_3_1"));
    VectPar.push_back(ParamToSave(SH,RAN_3_2,"RAN_3_2"));
    VectPar.push_back(ParamToSave(SH,RAN_1_2,"RAN_1_2"));
    VectPar.push_back(ParamToSave(SH,RAN_4_3,"RAN_4_3"));
    VectPar.push_back(ParamToSave(SH,RAN_2_3,"RAN_2_3"));
    VectPar.push_back(ParamToSave(SH,HetX,"HetX"));
    VectPar.push_back(ParamToSave(SH,HetY,"HetY",0,1));
    VectPar.push_back(ParamToSave(SH,Dbl8,"Dbl8"));
    VectPar.push_back(ParamToSave(SH,Stab,"Stab"));
    VectPar.push_back(ParamToSave(SH,HeiET,"HeiET"));
    VectPar.push_back(ParamToSave(SH,HeiEM,"HeiEM"));
    VectPar.push_back(ParamToSave(SH,HeiRT,"HeiRT"));
    VectPar.push_back(ParamToSave(SH,HeiRM,"HeiRM"));
    VectPar.push_back(ParamToSave(SH,HeiETT,"HeiETT",0,1));
    VectPar.push_back(ParamToSave(SH,HeiEMT,"HeiEMT"));
    VectPar.push_back(ParamToSave(SH,HeiRTT,"HeiRTT"));
    VectPar.push_back(ParamToSave(SH,HeiRMT,"HeiRMT"));
    VectPar.push_back(ParamToSave(SH,MedInt,"MedInt",0,1));
    VectPar.push_back(ParamToSave(SH,Int25,"Int25"));
    VectPar.push_back(ParamToSave(SH,Int75,"Int75"));
    VectPar.push_back(ParamToSave(SH,RInt1,"RInt1"));
    VectPar.push_back(ParamToSave(SH,IntDev,"IntDev",0,1));
    VectPar.push_back(ParamToSave(SH,SmIntDev,"SmIntDev"));
    VectPar.push_back(ParamToSave(SH,LgIntDev,"LgIntDev"));
    VectPar.push_back(ParamToSave(SH,VarInt,"VarInt"));
    VectPar.push_back(ParamToSave(SH,VarSmInt,"VarSmInt"));
    VectPar.push_back(ParamToSave(SH,VarLgInt,"VarLgInt"));
    VectPar.push_back(ParamToSave(SH,RIntDev1,"RIntDev1"));
    VectPar.push_back(ParamToSave(SH,EnStabSm,"EnStabSm"));
    VectPar.push_back(ParamToSave(SH,EnStabLg,"EnStabLg"));
    VectPar.push_back(ParamToSave(SH,HetXr,"HetXr",0,1));
    VectPar.push_back(ParamToSave(SH,HetYr,"HetYr"));
    VectPar.push_back(ParamToSave(SH,HetYr2,"HetYr2",0,1));
    VectPar.push_back(ParamToSave(SH,HetCMC,"HetCMC"));
    VectPar.push_back(ParamToSave(SH,HetCMD,"HetCMD"));
    VectPar.push_back(ParamToSave(SH,HetCTC,"HetCTC"));
    VectPar.push_back(ParamToSave(SH,HetCTD,"HetCTD"));
    VectPar.push_back(ParamToSave(SH,HetCMnP,"HetCMnP",0,1));
    VectPar.push_back(ParamToSave(SH,HetCMfP,"HetCMfP"));
    VectPar.push_back(ParamToSave(SH,HetCTnP,"HetCTnP",0,1));
    VectPar.push_back(ParamToSave(SH,HetCTfP,"HetCTfP"));
    VectPar.push_back(ParamToSave(SH,HetPicsMAD,"HetPicsMAD",0,1));
    VectPar.push_back(ParamToSave(SH,HetPicsMALD,"HetPicsMALD"));
    VectPar.push_back(ParamToSave(SH,HetPicsMABD,"HetPicsMABD"));
    VectPar.push_back(ParamToSave(SH,HetPicsMRBLD,"HetPicsMRLBD"));
    VectPar.push_back(ParamToSave(SH,HetPicsTAD,"HetPicsTAD",0,1));
    VectPar.push_back(ParamToSave(SH,HetPicsTALD,"HetPicsTALD",0,1));
    VectPar.push_back(ParamToSave(SH,HetPicsTABD,"HetPicsTABD"));
    VectPar.push_back(ParamToSave(SH,HetPicsTRBLD,"HetPicsTRLBD"));
    VectPar.push_back(ParamToSave(SH,VDPicsM,"VDPicsM",0,1));
    VectPar.push_back(ParamToSave(SH,VLDPicsM,"VLDPicsM",0,1));
    VectPar.push_back(ParamToSave(SH,VBDPicsM,"VBDPicsM",0,1));
    VectPar.push_back(ParamToSave(SH,VDPPicsM,"VDPPicsM",0,1));
    VectPar.push_back(ParamToSave(SH,VLDPPicsM,"VLDPPicsM"));
    VectPar.push_back(ParamToSave(SH,VBDPPicsM,"VBDPPicsM"));
    VectPar.push_back(ParamToSave(SH,VDPicsT,"VDPicsT",0,1));
    VectPar.push_back(ParamToSave(SH,VLDPicsT,"VLDPicsT",0,1));
    VectPar.push_back(ParamToSave(SH,VBDPicsT,"VBDPicsT",0,1));
    VectPar.push_back(ParamToSave(SH,VDPPicsT,"VDPPicsT",0,1));
    VectPar.push_back(ParamToSave(SH,VLDPPicsT,"VLDPPicsT"));
    VectPar.push_back(ParamToSave(SH,VBDPPicsT,"VBDPPicsT"));
    {
        for(int i=CM;i<=CO2;i++)
        {
            VectPar.push_back(ParamToSave(i,SDC,prefix[i]+"SDC",1,1));
            VectPar.push_back(ParamToSave(i,SDCR,prefix[i]+"SDCR",1));
        }
        VectPar.push_back(ParamToSave(CM,SDCRY,"CM_SDCRY",1,1));
        VectPar.push_back(ParamToSave(CS,SDCRY,"CS_SDCRY",1,1));
        VectPar.push_back(ParamToSave(CM,SDCRXY,"CM_SDCRXY",1));
        VectPar.push_back(ParamToSave(CS,SDCRXY,"CS_SDCRXY",1));
        for(int i=CM;i<=CS;i++)
        {
            if(i==CS) VectPar.push_back(ParamToSave(i,SDCL,prefix[i]+"SDCL",1,1));
            else VectPar.push_back(ParamToSave(i,SDCL,prefix[i]+"SDCL",1));
            VectPar.push_back(ParamToSave(i,SDCLR,prefix[i]+"SDCLR",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRY,prefix[i]+"SDCLRY",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRXY,prefix[i]+"SDCLRXY",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRXY2,prefix[i]+"SDCLRXY2",1,1));
            VectPar.push_back(ParamToSave(i,SDCLOP,prefix[i]+"SDCLOP",1));
            VectPar.push_back(ParamToSave(i,SDCLROP,prefix[i]+"SDCLROP",1));
            if(i==CM) VectPar.push_back(ParamToSave(i,SDCLRYOP,prefix[i]+"SDCLRYOP",1,1));
            else VectPar.push_back(ParamToSave(i,SDCLRYOP,prefix[i]+"SDCLRYOP",1));
            VectPar.push_back(ParamToSave(i,SDCLRXYOP,prefix[i]+"SDCLRXYOP",1,1));
            if(i==CM) VectPar.push_back(ParamToSave(i,SDCLWB,prefix[i]+"SDCLWB",1,1));
            else VectPar.push_back(ParamToSave(i,SDCLWB,prefix[i]+"SDCLWB",1));
            if(i==CS) VectPar.push_back(ParamToSave(i,SDCLRWB,prefix[i]+"SDCLRWB",1,1));
            else VectPar.push_back(ParamToSave(i,SDCLRWB,prefix[i]+"SDCLRWB",1));
            VectPar.push_back(ParamToSave(i,SDCLRYWB,prefix[i]+"SDCLRYWB",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRXYWB,prefix[i]+"SDCLRXYWB",1,1));
            VectPar.push_back(ParamToSave(i,SDCLOPWB,prefix[i]+"SDCLOPWB",1,1));
            VectPar.push_back(ParamToSave(i,SDCLROPWB,prefix[i]+"SDCLROPWB",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRYOPWB,prefix[i]+"SDCLRYOPWB",1,1));
            if(i==CS) VectPar.push_back(ParamToSave(i,SDCLRXYOPWB,prefix[i]+"SDCLRXYOPWB",1,1));
            else VectPar.push_back(ParamToSave(i,SDCLRXYOPWB,prefix[i]+"SDCLRXYOPWB",1));
            VectPar.push_back(ParamToSave(i,SDCL_DNP,prefix[i]+"SDCL_DNP",1,1));
            VectPar.push_back(ParamToSave(i,SDCLR_DNP,prefix[i]+"SDCLR_DNP",1));
            if(i==CM) VectPar.push_back(ParamToSave(i,SDCLRY_DNP,prefix[i]+"SDCLRY_DNP",1,1));
            else VectPar.push_back(ParamToSave(i,SDCLRY_DNP,prefix[i]+"SDCLRY_DNP",1));
            VectPar.push_back(ParamToSave(i,SDCLRXY_DNP,prefix[i]+"SDCLRXY_DNP",1,1));
            VectPar.push_back(ParamToSave(i,SDCLRXY2_DNP,prefix[i]+"SDCLRXY2_DNP",1,1));
        }
        VectPar.push_back(ParamToSave(CM,ELBPOS,"CM_ELBPOS",1));
        VectPar.push_back(ParamToSave(CS,ELBPOS,"CS_ELBPOS",1));
        VectPar.push_back(ParamToSave(CM,ELBSB,"CM_ELBSB",1));
        VectPar.push_back(ParamToSave(CS,ELBSB,"CS_ELBSB",1));
        VectPar.push_back(ParamToSave(CM,ELB2POS,"CM_ELB2POS",1));
        VectPar.push_back(ParamToSave(CS,ELB2POS,"CS_ELB2POS",1));
        VectPar.push_back(ParamToSave(CM,ELB2SB,"CM_ELB2SB",1));
        VectPar.push_back(ParamToSave(CS,ELB2SB,"CS_ELB2SB",1));
        VectPar.push_back(ParamToSave(CM,RAF,"CM_RAF",1,1));
        VectPar.push_back(ParamToSave(CM,RAE,"CM_RAE",1,1));
        VectPar.push_back(ParamToSave(CM,RAFE,"CM_RAFE",1));
        VectPar.push_back(ParamToSave(CM,RAFP,"CM_RAFP",1,1));
        VectPar.push_back(ParamToSave(CM,RAFP2,"CM_RAFP2",1,1));
        VectPar.push_back(ParamToSave(CM,RAFP3,"CM_RAFP3",1));
        VectPar.push_back(ParamToSave(CM,SBMP,"CM_SBMP",1));
        VectPar.push_back(ParamToSave(CM,SAMP,"CM_SAMP",1));
        VectPar.push_back(ParamToSave(CM,SBAR,"CM_SBAR",1));
        VectPar.push_back(ParamToSave(CM,RAHP2,"RAHP2",1,1));
        VectPar.push_back(ParamToSave(CM,RAHP4,"RAHP4",1,1));
        VectPar.push_back(ParamToSave(CM,RAHP8,"RAHP8",1,1));
        VectPar.push_back(ParamToSave(CM,RAHP16,"RAHP16",1,1));
        VectPar.push_back(ParamToSave(CM,RAHE2,"RAHE2",1,1));
        VectPar.push_back(ParamToSave(CM,RAHE4,"RAHE4",1));
        VectPar.push_back(ParamToSave(CM,RAHE8,"RAHE8",1,1));
        VectPar.push_back(ParamToSave(CM,RAHE16,"RAHE16",1,1));
    } 
}

// this method deletes dynamically created objects and frees memory allocated
// in the method:InitializeDetecTreatment()
void DetecTreatment::EndDetecTreatment()
{
    fftwf_free(_data);
    delete _coeff;
    delete[] _charSonogramArray;
    delete[] _charPointFlagsArray;
    delete[] SonogramArray;
    delete[] PointFlagsArray;
    for (int i=0;i < MAXCRI; i++) delete[] _paramsArray[i];
    delete[] _paramsArray;
    delete[] _charParamsArray;
    delete[] LowSlope;
    delete[] Inflexion1;
    delete[] Inflexion3;
    for(int k=0;k<2;k++)
    {
        delete[] _harmonic[k];
        delete[] _dpm[k];
        delete[] _dypm[k];
        delete[] _ypm[k];
    }
    delete[] _harmonic;
    delete[] _dpm;
    delete[] _dypm;
    delete[] _ypm;
    delete[] _tabY;
    delete[] _numberPixelsPerY;
    delete[] _numberPixelsPerX;
    delete[] _averagePerX;
    delete[] _xMinPerY;
    delete[] _xMaxPerY;
    delete[] _yMinPerX;
    delete[] _yMaxPerX;
    delete[] _yEbarPerX;
    delete[] _totYEPerX;
    delete[] _eMaxPerX;
    delete[] _xSecondWestRidgePerY;
    delete[] _slope;
    delete[] _tabX;
    delete[] _charTabX;
    delete[] _yEmaxPerX;
    delete[] _charYEmaxPerX;
    delete[] _tabYX;
    delete[] _charTabYX;
    delete[] EnergyColumAverage;
    delete[] FlagGoodCol;
    delete[] FlagGoodColInitial;
    delete[] _sortMp;
    delete[] _invMp;
    delete[] _xMp;
}

// CallTreatmentsForOneFile: a very important method which is called by the object
// of detec class to treat one file
// this method calls the main methods of treatment for each file
//     openWavFile
//     correctNoise
//     shapesDetects
//     detectsParameter2
//     saveParameters
bool DetecTreatment::CallTreatmentsForOneFile(QString& wavFile,QString &pathFile)
{
    _detec->LogStream << "c:" << endl;
    qint64 d0 = QDateTime::currentDateTime().toMSecsSinceEpoch();
   _wavFile = wavFile;
    int d[5];
    clearVars();
    NError = -1;
    if (openWavFile(pathFile))
    {
        if(computeFFT(pathFile))
        {
            d[0]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->LogStream <<   "cF"<< endl;
            correctNoise();
            d[1]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->LogStream <<   "cN"<< endl;
            shapesDetects();
            d[2]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->LogStream <<   "sD: " << endl;
            _callsNumber = (int)CallsArray.size();
            detectsParameter2();
            d[3]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->LogStream <<   "dP" << endl;
            saveParameters(wavFile);
            if(_detec->MustCompress) saveCompressedParameters(wavFile);
            d[4]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->LogStream <<   "sP" << endl;
            if(_detec->TimeFileOpen)
            {
                _detec->TimeStream  << wavFile << '\t' << d[0] << '\t';
                for(int j=1;j<=4;j++) _detec->TimeStream  << d[j] - d[j-1] << '\t';
                _detec->TimeStream  << d[4] << endl;
            }
        }
        else
        {
            return(false);
        }
    }
    else
    {
        return(false);
    }
    return(true);
}

// this method cleans variables which must be cleaned for each treatment of a file
void DetecTreatment::clearVars()
{
    CallsArray.clear();
    _vectorXMin.clear();
    _vectorXMax.clear();
    MasterPoints.clear();
    if(_imageData)
    {
        CallMasterRidgeArray.clear();
        CallSouthArray.clear();
        CallNorthRidgeArray.clear();
        CallWestRidgeArray.clear();
        CallSecondWestRidgeArray.clear();
    }
}

// debug informations recorded in debug mode
void DetecTreatment::aff(QString name,qint64 address,int size)
{
    /*
    _detec->LogStream << "P:"<< name << " : " << address ;
    if(size==0) _detec->LogStream << " fr" << endl;
    else _detec->LogStream << " s=" << size << endl;
    */
}

// this method determines if a .wav file corresponds
// to a left / mono or right channel of a recorder
// from conventions on the file name in the vigiechiro protocol
// to be modified for other uses
bool DetecTreatment::determineLeftOrRight(QString& wavFile)
{
    int ct = wavFile.count("_");
    QString toLook;
    if(ct>2)
    {
        for(int i=ct-1;i>=0;i--)
        {
            toLook = wavFile.section("_",i,i);
            if(toLook.length()==1)  break;
        }
    }
    else toLook = "-9";
    if(toLook.toInt()==1) return(false);
    return(true);
}

// this method uses the sndfile library to verify the validity of the .wav file
// very important variables for the treatment are computed : _fftHeight and _iH
bool DetecTreatment::openWavFile(QString& pathFile)
{

    if (! (_soundFile = sf_open(pathFile.toStdString().c_str(), SFM_READ, &_soundFileInfo)))
    {
        if(_detec->ErrorFileOpen) _detec->ErrorStream << pathFile << ": " << sf_strerror (NULL) << endl;
        NError=FNREC;
        return  false;
    }
    if (_soundFileInfo.channels > 1)
    {
        sf_close (_soundFile);
        if(_detec->ErrorFileOpen) _detec->ErrorStream << pathFile << ": " << "multi-channel is not managed" << endl;
        NError=MCNT;
        return  false;
    }
    if(_detec->ReprocessingMode)
    {
        if(_detec->NumVer > 19)
        {
            TimeExpansion = _detec->TE;
        }
        else
        {
            float ate = (_detec->NumtE) / ((float)_soundFileInfo.samplerate);
            if(ate < 2) TimeExpansion = 1; else TimeExpansion = 10;
        }
    }
    else
    {

        if(determineLeftOrRight(_wavFile))
        {
            TimeExpansion = _timeExpansionLeft;
            _detec->LogStream << "g=" << TimeExpansion << endl;

        }
        else
        {
            TimeExpansion = _timeExpansionRight;
            _detec->LogStream << "d=" << TimeExpansion << endl;
        }
    }
    if(TimeExpansion<=0)
    {
        sf_close (_soundFile);
        if(_detec->ErrorFileOpen) _detec->ErrorStream << _wavFile << ": undefined time factor" << endl;
        _detec->LogStream << _wavFile << ": ft ndef" << endl;

        NError=TNT;
        return(false);
    }
    _detec->LogStream << "sfte=" << _soundFileInfo.samplerate*TimeExpansion << endl;
    if(_modeFreq==2)
    {
        if (_soundFileInfo.samplerate*TimeExpansion >= 480000 ) {_fftHeight = 4096; _iH =5;}
        else {
            if (_soundFileInfo.samplerate*TimeExpansion >= 240000 ) {_fftHeight = 2048;  _iH =4;}

            else {

                if (_soundFileInfo.samplerate*TimeExpansion >= 120000 ) {_fftHeight = 1024;  _iH =3;}

                else {
                    if (_soundFileInfo.samplerate*TimeExpansion >= 60000 ) {_fftHeight = 512;  _iH =2;}
                    else {
                        if (_soundFileInfo.samplerate*TimeExpansion >= 30000 ) {_fftHeight = 256;  _iH =1;}
                        else
                            { _fftHeight = 128;  _iH =0;}
                    }
                }
            }
        }

    }
    else
    {
        if (_soundFileInfo.samplerate*TimeExpansion >= 2400000 ) {_fftHeight = 4096; _iH =5;}
        else {
            if (_soundFileInfo.samplerate*TimeExpansion >= 1200000 ) {_fftHeight = 2048;  _iH =4;}

            else {

                if (_soundFileInfo.samplerate*TimeExpansion >= 600000 ) {_fftHeight = 1024;  _iH =3;}

                else {
                    if (_soundFileInfo.samplerate*TimeExpansion >= 300000 ) {_fftHeight = 512;  _iH =2;}
                    else {
                        if (_soundFileInfo.samplerate*TimeExpansion >= 150000 ) {_fftHeight = 256;  _iH =1;}
                        else
                            { _fftHeight = 128;  _iH =0;}
                    }
                }
            }
        }
    }
    return true;
}

// the computeFFT method uses fftw3 library (Fast Fourier transforms algorithms)
// at the end of this treatment, the array: SonogramArray contains a matrix: frequency/time/energy
// used for the detection of soundevents and the calculation of their parameters
bool DetecTreatment::computeFFT(QString &wavFile)
{
    int iCount;
    int readcount;
    double a = 0.0d;
    double epsilon = 0.0001d;
    EnergyMax        = 0.0d;
    EnergyMin        = 0.0d;
    FftHeightHalf		= (int)ceil((float)_fftHeight/2.0f);
    _iOverlapMoving	= (int)ceil((float)_fftHeight/(float)(_nbo*2));
    SonogramWidth		= (int)ceil(_nbo*2*(float)_soundFileInfo.frames/(float)_fftHeight)+1;
    MsPerX =(double)(FftHeightHalf*1000)/(_nbo*_soundFileInfo.samplerate*TimeExpansion); //Time: msec
    KhzPerY =(double)(_soundFileInfo.samplerate*TimeExpansion)/(double)(_fftHeight*1000); //Freq:khz
    if((double)SonogramWidth*MsPerX < 10.0d)
    {
        if(_detec->ErrorFileOpen) _detec->ErrorStream << wavFile << ": too short duration" << endl;
        _detec->LogStream << "too short duration:" << SonogramWidth*MsPerX << " ms" << endl;
        NError=DTP;
        return  false;
    }
    if(SonogramWidth > SONOGRAM_WIDTH_MAX)
    {
        if(_detec->ErrorFileOpen) _detec->ErrorStream << wavFile << ": too long duration" << endl;
        _detec->LogStream <<  "too long" << endl;
        NError=DTG;
        return  false;
    }
    sf_seek(_soundFile, 0, SEEK_END);
    _pPlan = &(_detec->PMainWindow->Plan[_detec->IThread][_iH]);
    double fact1=2.0d*PI;
    double fact2=4.0f*PI;
    double quot1=FftHeightHalf-1;
    LimY = qMin(FftHeightHalf,MAXHEIGHT);
    for (int i = 0 ; i < FftHeightHalf ; i++)
    {
        _coeff[i] = 0.435d - 0.5d*cos(fact1*i/quot1)+ 0.065d*cos(fact2*i/quot1);
    }
    for (int iLoop = 0 ; iLoop < _nbo; iLoop++)
    {
        iCount = 0;
        sf_seek(_soundFile, iLoop * _iOverlapMoving, SEEK_SET);
        while ((readcount = (int)sf_read_float(_soundFile, _data, FftHeightHalf)))
        {
            for (int i = 0 ; i < FftHeightHalf ; i++)
            {
                _complexInput[i][0] =_data[i] * _coeff[i];
                _complexInput[i][1] = 0.0f;
            }
            for (int i = FftHeightHalf ; i < _fftHeight ; i++)
            {
                _complexInput[i][1] = 0.0f;
                _complexInput[i][0] = 0.0f;
            }
            fftwf_execute( *_pPlan );
            qint16 *sml;
            int jc=iCount*_nbo+iLoop;
            double b;
            for(int i =0; i < LimY;  i++)
            {
                sml=SonogramArray[i];
                a = pow(_fftHeight*_fftRes[i][0], 2) + pow(_fftHeight*_fftRes[i][1], 2);
                if (i > _freqMin && a !=0 )
                {
                    b=10*log10(a)+epsilon;
                    sml[jc] = (qint16)(b*100.0d);
                    if(b>EnergyMax) EnergyMax=b;
                    else
                    {
                        if(b<EnergyMin) EnergyMin =b;
                    }
                }
                else 
				{
				sml[jc] = -5000;
				}
            }
            iCount++;
        }
    }
    for(int i =0; i < LimY; i++) SonogramArray[i][SonogramWidth-1]=0;
    EnergyMax = qMax(EnergyMax, (double)0);
    sf_close (_soundFile);
    return true;
}

// this method contains algorithms to eliminate background noise
// to improve the detection of cries...
void DetecTreatment::correctNoise()
{
    _minY=(int)(((double)_freqMin)/KhzPerY);
    if(_modeFreq == 2)
        _maxY=(int)(((double)FREQ_MAX)/(10.0d*KhzPerY));
    else
        _maxY=(int)(((double)FREQ_MAX)/KhzPerY);
    if(_maxY>LimY-1) _maxY = LimY-1;
    if(_detec->IDebug)
	{
        _detec->LogStream << "_minY = " << _minY << endl;
        _detec->LogStream << "_maxY = " << _maxY << endl;
        _detec->LogStream << "_freqMin = " << _freqMin << endl;
        _detec->LogStream << "FREQ_MAX = " << FREQ_MAX << endl;
	}
    if(_desactiveCorrectNoise)
    {
        qint16 *fc;
        int retrait = _stopThreshold * 100;
        for(int y = _minY; y <= _maxY ; y++)
        {
            fc=SonogramArray[y];
            for (int x = 0 ; x < SonogramWidth ; x++) fc[x] -= retrait;
        }
        return;
    }
    int son_min = EMIN,son_max = EMIN+199;
    int minEmc = qMax((int)(20.0f/KhzPerY),_minY);
    int maxEmc = qMin((int)(80.0f/KhzPerY),_maxY);
    // neutralizes silent columns
    for (int x = 0 ; x < SonogramWidth ; x++) EnergyColumAverage[x] = 0.0f;
    //for(int y = _minY; y < _maxY ; y++)
    int decalThreshold = 0;
    if(minEmc < maxEmc)
    {
        double totEmc = 0.0d;
        for(int y = minEmc; y < maxEmc ; y++)
        {
            qint16 *fc = SonogramArray[y];
            for (int x = 0 ; x < SonogramWidth ; x++) EnergyColumAverage[x]  += (float)(*fc++);
        }
        double diviseur = 100.0d * ( (double) (maxEmc - minEmc) );
        float *pemc = EnergyColumAverage;
        for (int x = 0 ; x < SonogramWidth ; x++) {*pemc /= diviseur; totEmc += (*pemc++);}
        decalThreshold = (int)( ((totEmc/SonogramWidth) - ((float)(_lowThresholdJB+_highThresholdJB)/2.0f))/4      );
    }
    bool unSaut = false;
    WithSilence = false;
    int highThresholdJB = _highThresholdJB+decalThreshold;
    int lowThresholdJB = _lowThresholdJB+decalThreshold;
    int lowThresholdC = _lowThresholdC+decalThreshold;
    int highThresholdC = _highThresholdC+decalThreshold;
    if(SonogramWidth>10 && _useValflag)
    {
        bool valFlag = true;
        int patience = 0;
        int jumpThreshold = _jumpThreshold;
        int widthBigControl = _widthBigControl;
        int widthLittleControl = _widthLittleControl;
        for (int x = 0 ; x <= widthBigControl ; x++) FlagGoodColInitial[x]=valFlag;
        double dif2Col,averageLittleBefore,averageLittleNext,averageBigBefore,averageBigNext;
        bool notYet = true;
        if(widthBigControl > SonogramWidth/3)
        {
            widthBigControl = SonogramWidth/3;
            if(widthLittleControl>widthBigControl/2) widthLittleControl=widthBigControl/2;
        }
        for (int x = widthBigControl+1 ; x < SonogramWidth-widthBigControl ; x++)
        {
            averageLittleBefore = 0.0f;
            averageLittleNext   = 0.0f;
            for(int j=0;j<widthLittleControl;j++)
            {
                averageLittleNext   += EnergyColumAverage[x+j];
                averageLittleBefore += EnergyColumAverage[x-1-j];
            }
            averageLittleBefore /= widthLittleControl;
            averageLittleNext   /= widthLittleControl;
            dif2Col = averageLittleNext - averageLittleBefore;
            // -- -- --
            if(dif2Col>jumpThreshold && patience>widthBigControl/3
                    && averageLittleNext >highThresholdJB)
            {
                averageBigBefore = 0.0f;
                averageBigNext   = 0.0f;
                for(int j=0;j<widthBigControl;j++)
                {
                    averageBigNext   += EnergyColumAverage[x+j];
                    averageBigBefore += EnergyColumAverage[x-1-j];
                }
                averageBigBefore /= widthBigControl;
                averageBigNext   /= widthBigControl;
                if((averageBigNext-averageBigBefore)>jumpThreshold
                        && averageBigNext >highThresholdJB)
                {

                    patience = 0;
                    if(valFlag==false)
                    {
                        notYet =false;
                        unSaut = true;
                    }

                    if(valFlag==true && notYet==true
                            && averageLittleBefore < lowThresholdJB
                            && averageBigBefore < lowThresholdJB)
                    {
                        for(int j=x-1;j>=0;j--)
                        {
                            if(FlagGoodColInitial[j]==true) FlagGoodColInitial[j]=false;
                            else break;
                        }
                        unSaut = true;
                        notYet =false;
                    }
                    valFlag=true;
                }
            }
            // -- -- --
            if(dif2Col<-jumpThreshold && patience>widthBigControl/3
                    && averageLittleNext<lowThresholdJB)
            {
                averageBigBefore = 0.0f;
                averageBigNext   = 0.0f;
                for(int j=0;j<widthBigControl;j++)
                {
                    averageBigNext   += EnergyColumAverage[x+j];
                    averageBigBefore += EnergyColumAverage[x-1-j];
                }
                averageBigBefore /= widthBigControl;
                averageBigNext   /= widthBigControl;
                if((averageBigNext-averageBigBefore)<-jumpThreshold
                        && averageBigNext<lowThresholdJB)
                {
                    unSaut = true;
                    patience = 0;
                    valFlag=false;
                    notYet =false;
                }
            }
            FlagGoodColInitial[x] = valFlag;
            patience ++;
        }
        for (int x = SonogramWidth-widthBigControl;x < SonogramWidth;x++) FlagGoodColInitial[x]=valFlag;
    }
    int nff = 0;
    if(unSaut)
    {
        for (int x = 0;x < SonogramWidth;x++)
        {
            if(FlagGoodColInitial[x]==false)
            {
                if(EnergyColumAverage[x]<highThresholdC) {nff++; FlagGoodCol[x]=false;}
                else FlagGoodCol[x] = true;
            }
            else
            {
                if(EnergyColumAverage[x]<lowThresholdC) {nff++; FlagGoodCol[x] = false;}
                else FlagGoodCol[x]=true;
            }
        }
        WithSilence = true;
        if(nff*2 >SonogramWidth)
        {
            for(int jpha=0;jpha<2;jpha++)
            {
                int totTrue = 0;
                int maxTrue = 0;
                int actualWidth = 0;
                for (int x = 0;x < SonogramWidth;x++)
                {
                    if(jpha==1 && FlagGoodCol[x]==false)
                    {
                        if(EnergyColumAverage[x]>lowThresholdC) {nff--; FlagGoodCol[x]=true;}
                    }
                    if(FlagGoodCol[x]==true)
                    {
                        actualWidth++;
                        totTrue++;
                        if(actualWidth>maxTrue) maxTrue = actualWidth;
                    }
                    else actualWidth = 0;
                }
                if(maxTrue*10>SonogramWidth)
                {
                    if(jpha==1)
                    {
                        if(nff==0)
                        {
                            WithSilence = false;
                        }
                    }
                    break;
                }
                else
                {
                    if(jpha==1)
                    {
                        for (int x = 0;x < SonogramWidth;x++) FlagGoodCol[x]=true;
                        nff = 0;
                        WithSilence = false;
                    }
                }
            }
        }
    }
    else
    {
        for (int x = 0;x < SonogramWidth;x++) FlagGoodColInitial[x]=true;
        for (int x = 0;x < SonogramWidth;x++) FlagGoodCol[x]=true;
    }

    if(nff>0)
    {
        _detec->LogStream << "nff = " << nff << " sur " << SonogramWidth << endl;
    }

    int tval[200];
    qint16 *fc;
    int largeurRectifiee;
    for(int y = _minY; y <= _maxY ; y++)
    {
        fc=SonogramArray[y];
        for(int k=0;k<200;k++) tval[k]=0;

        largeurRectifiee = 0;
        for (int x = 0 ; x < SonogramWidth ; x++)
        {
            int son = (fc[x]+50)/100;
            if(son>=son_min && son<=son_max && FlagGoodCol[x])
            {
                tval[son-son_min]++;
                largeurRectifiee++;
            }
            if(son < son_min) fc[x] = son_min*100;
            if(son > son_max) fc[x] = son_max*100;
        }

        int cumul = 0, q5 = largeurRectifiee*_qR/100;
        if(q5 < _qN) q5 = _qN;
        for(int j=0;j<=son_max-son_min;j++)
       {
            cumul += tval[j];
            if(cumul>=q5)
            {
                int retrait = (j + son_min + _stopThreshold) * 100;
                for (int x = 0 ; x < SonogramWidth ; x++) fc[x] -= retrait;
                break;
            }
        }
    }
}

// this method determines sound events
// its algorithm detects points of the time/frequency matrix
// whose energy is greater than a threshold (_detectionThreshold)
// and associates all the adjoining points whose energy is greater than a second
// threshold (_stopThreshold)
void DetecTreatment::shapesDetects()
{
    double epsilon = 0.0001d;
    for(int j=0;j<MAXHEIGHT;j++) memset(PointFlagsArray[j],0,LD8);
    _maxCallWidth = 0;
    _maxCallHeight = 0;
    QPoint Point;
    int ix,iy;
    int curseur;
    int ncontour=0;
    EnergyShapeThreshold = (double)_detectionThreshold-_stopThreshold;
    EnergyStopThreshold = 0.0f;
    int nbcont=0;
    qint16 *smy;
    bool onsarrete = false;
    for(int y = (_maxY-1); y >= _minY ; y--)
    {
        smy = SonogramArray[y];
        int digitPos = 0;
        char *pBoolChar = PointFlagsArray[y];
        for (int x = 0 ; x < SonogramWidth ; x++)
        {
            if(((*pBoolChar) & (1 << digitPos))==0)
                if (smy[x] > EnergyShapeThreshold*100)
                {
                    nbcont++;
                    ncontour++;
                    _vectorCallPoints.clear();
                    Point.setX(x);
                    Point.setY(y);
                    _callEnergyMax = (double)smy[x] / 100.0d;
                    _callEnergyMaxIndex = 0;
                    _vectorCallPoints.push_back(Point);
                    _xMin=x; _xMax=x;
                    _yMin=y; _yMax=y;
                    curseur = 0;
                    *pBoolChar |= (1 << digitPos);
                    while(curseur < _vectorCallPoints.size())
                    {
                        Point=_vectorCallPoints.at(curseur);
                        ix=Point.x();
                        iy=Point.y();
                        for(int jy=iy-1;jy<=iy+1;jy++)
                            // edit yves - broadening spectrum
                        {
                            if(jy>=_minY && jy <= _maxY)
                            {
                                int retx = 5;
                                if(ix<retx) retx = ix;
                                int posBC = (ix-retx)/8;
                                int dP = ix-retx-posBC*8;
                                char *pBC = PointFlagsArray[jy]+posBC;
                                char bC = *pBC;
                                for(int jx=ix-retx;jx<=ix+5;jx++)
                                {
                                    if(jx!=ix || jy!=iy)
                                    {
                                        if(jx>= 0 && jx < SonogramWidth)
                                        {
                                            if((bC & (1 << dP))==0)
                                            {
                                                double val = (double)SonogramArray[jy][jx]/100.0d;
                                                if (val > EnergyStopThreshold)
                                                {
                                                    Point.setX(jx);
                                                    Point.setY(jy);
                                                    *pBC |= (1 << dP);
                                                    bC = *pBC;

                                                    _vectorCallPoints.push_back(Point);
                                                    if(jx<_xMin) _xMin=jx; else {if(jx>_xMax) _xMax=jx;}
                                                    if(jy<_yMin) _yMin=jy; else {if(jy>_yMax) _yMax=jy;}
                                                    if(val > _callEnergyMax)
                                                    {
                                                        _callEnergyMax = val;
                                                        _callEnergyMaxIndex = _vectorCallPoints.size()-1;
                                                    }
                                                }
                                            }

                                        }

                                    }
                                    if(jx<SonogramWidth-1)
                                    {
                                        dP++;
                                        if(dP==8) {pBC++; bC = *pBC; dP=0;}
                                    }
                                }
                            }
                        }
                        curseur++;
                        if(curseur > 300000)
                        {
                            _detec->LogStream << "limited to 300000 pixels" << endl;
                            break;
                        }
                    }
                    double freqpm=((double)_vectorCallPoints.at(_callEnergyMaxIndex).y())*KhzPerY;

                    if(freqpm>(_freqCallMin+epsilon) && (_xMax-_xMin+1)<MAXLARCRI
                            && (_yMax-_yMin+1)<MAXHEIGHT
                            && !onsarrete)
                    {
                        MasterPoints.push_back(_vectorCallPoints.at(_callEnergyMaxIndex));
                        CallsArray.push_back(_vectorCallPoints);
                        _vectorXMin.push_back(_xMin);
                        _vectorXMax.push_back(_xMax);
                        if((_xMax-_xMin+1)>_maxCallWidth) _maxCallWidth=_xMax-_xMin+1;
                        if((_yMax-_yMin+1)>_maxCallHeight) _maxCallHeight=_yMax-_yMin+1;
                        if(CallsArray.size()==MAXCRI-1) onsarrete = true;
                        if(_xMax-_xMin+1 > 5000)
                        {
                            _detec->LogStream << "cryWidth=" << _xMax-_xMin+1 << endl;
                        }
                    }
                }

            if(onsarrete) break;
            if(x<SonogramWidth-1)
            {
                digitPos++;
                if(digitPos==8) {pBoolChar++; digitPos=0;}
            }
        }
        if(onsarrete) break;
    }
    if(CallsArray.size() > MAXCRI)
        _detec->LogStream << "nCries = " << CallsArray.size() << "   " << _wavFile << endl;

    sortWaves();
}

// this method sorts the sound events according to their start time
// it is necessary for some settings computing depending on those of other
// preceding or following sound events
void DetecTreatment::sortWaves()
{
    bool ontrie = true;
    int nbcris = CallsArray.size();
    while(ontrie)
    {
        ontrie = false;
        for(int i=0;i<nbcris-1;i++)
        {
            int x1 = _vectorXMin.at(i);
            int x2 = _vectorXMin.at(i+1);
            int x3 = _vectorXMax.at(i);
            int x4 = _vectorXMax.at(i+1);
            bool permuter = false;
            if(x1>x2) permuter = true;
            if(permuter)
            {
                QVector<QPoint> mvco1=CallsArray.at(i);
                QVector<QPoint> mvco2=CallsArray.at(i+1);
                CallsArray.replace(i,mvco2);
                CallsArray.replace(i+1,mvco1);
                QPoint pm1=MasterPoints.at(i);
                QPoint pm2=MasterPoints.at(i+1);
                MasterPoints.replace(i,pm2);
                MasterPoints.replace(i+1,pm1);
                _vectorXMin.replace(i,x2);
                _vectorXMin.replace(i+1,x1);
                _vectorXMax.replace(i,x4);
                _vectorXMax.replace(i+1,x3);
                ontrie = true;
            }
        }
    }
}

// this very long method computes all the settings to be recorded in the ".ta" file
// read the description of settings in tadarida manual
void DetecTreatment::detectsParameter2()
{
    float *oParam;
    float *oParamCrete[NCRETES];
    double epsilon = 0.0001d;
    //
    //if(_detec->IDebug) _detec->_logText << "dp2-1" << endl ;
    int nbcris = CallsArray.size();
    if(nbcris< 1 || _maxCallWidth < 1 || _maxCallHeight < 1) return;
    int maxlarhau = _maxCallWidth;
    if(_maxCallHeight>maxlarhau) maxlarhau=_maxCallHeight;
    for(int i=0;i<nbcris;i++) {_sortMp[i]=i; _xMp[i]=MasterPoints.at(i).x();}
     sortIntArrays(_sortMp,nbcris,_xMp);
     for(int i=0;i<nbcris;i++)  _invMp[_sortMp[i]]=i;
     //if(_detec->IDebug) _detec->_logText << "-3" << "nbCries="<< nbcris << endl;
    for (int icri = 0 ; icri < nbcris ; icri++) //Execute for each call
    {
        //if(_detec->IDebug) _detec->_logText << endl << "-3-icri=" << icri << endl ; //+++
        oParam = _paramsArray[icri][SH];
        for(int j=0;j<NCRETES;j++) oParamCrete[j] = _paramsArray[icri][j+1];
        //if(_detec->IDebug) _detec->_logText << "-4" << endl;  //+++
        QVector<QPoint> unemat = CallsArray.at(icri);
        int tailleforme = unemat.size(); // The number of vectors of each call
        int xmin = SonogramWidth-1;
        int xmax = 0;
        int ymin=_maxY;
        int ymax=_minY;
        for(int j=0;j<tailleforme;j++)
        {
            int x=unemat.at(j).x();
            int y=unemat.at(j).y();
            if(x>xmax) xmax = x;
            if(x<xmin) xmin = x;
            if(y>ymax) ymax = y;
            if(y<ymin) ymin = y;
        }
        for(int k=0;k<=ymax-ymin;k++)
        {
            _numberPixelsPerY[k]=0;
            _tabY[k]=0.0f;
            _xMinPerY[k]=xmax+1;
            _xMaxPerY[k]=0;
            for(int l=0;l<=xmax-xmin;l++) _tabYX[k][l]=0.0f;
        }
        for(int k=0;k<=xmax-xmin;k++)
        {
            _tabX[icri][k]=0.0f;
            _numberPixelsPerX[k]=0;
            _yMinPerX[k]=0;
            _yMaxPerX[k]=0;
            _yEbarPerX[k]=0.0f;
            _totYEPerX[k]=0.0f;
            _yEmaxPerX[icri][k]=0;
            _eMaxPerX[k]=0.0f;
        }
        double eTot = 0.0d;
        double eMoy = 0.0d;

        for(int j=0;j<tailleforme;j++)
        {
            int x=unemat.at(j).x();
            int y=unemat.at(j).y();
            double e = (double)SonogramArray[y][x]/100.0f;
            _tabY[y-ymin]+=e;
            eTot += e;
            _numberPixelsPerY[y-ymin]++;
            _tabX[icri][x-xmin]+=e;

            _numberPixelsPerX[x-xmin]++;
            _totYEPerX[x-xmin]+=e*(double)y;
            if(y>_yMaxPerX[x-xmin]) _yMaxPerX[x-xmin]=y;
            if(y<_yMinPerX[x-xmin] || _yMinPerX[x-xmin]==0) _yMinPerX[x-xmin]=y;
            if(x>_xMaxPerY[y-ymin]) _xMaxPerY[y-ymin]=x;
            if(x<_xMinPerY[y-ymin]) _xMinPerY[y-ymin]=x;

            if(e>(_eMaxPerX[x-xmin]+epsilon))
            {
                _eMaxPerX[x-xmin]=e;
                _yEmaxPerX[icri][x-xmin]=y;
            }
            _tabYX[y-ymin][x-xmin]=e;

        }
        eMoy = eTot / tailleforme;
        double erec; int xc;
        for(int k=0;k<=ymax-ymin;k++)
        {
            xc=_xMinPerY[k];
            //if(xc>=xmin)
            //{
                erec=_tabYX[k][xc-xmin];
                int patience=0;
                if(xmax>xc)
                {
                    for(int l=xc-xmin+1;l<=xmax-xmin;l++)
                    {
                        if(_tabYX[k][l]>=erec)
                        {
                            erec=_tabYX[k][l];
                            xc=xmin+l;
                            patience=0;
                        }
                        else
                        {
                            patience++;
                            if(patience>3) break;
                        }
                    }
                }
            //}
            _xSecondWestRidgePerY[k]=xc;
        }
        int yfds = 0;
        int yfdm = 0;
        double recfds= 0.0d;
        double recfdm= 0.0d;
        for(int k=0;k<=ymax-ymin;k++)
        {
            if(_tabY[k]>(recfds+epsilon))
            {
                yfds=ymin+k;
                recfds=_tabY[k];
            }
            if(_numberPixelsPerY[k]>0)
                if(_tabY[k]/_numberPixelsPerY[k]>(recfdm+epsilon))
                {
                    yfdm=ymin+k;
                    recfdm=_tabY[k]/_numberPixelsPerY[k];
                }
        }
        int xpps = 0;
        int xppm = 0;
        double recpps= 0.0d;
        double recppm= 0.0d;
        for(int k=0;k<=xmax-xmin;k++)
        {
            if(_tabX[icri][k]>(recpps+epsilon))
            {
                xpps=xmin+k;
                recpps=_tabX[icri][k];
            }
            if(_numberPixelsPerX[k]>0)
            {
                if(_tabX[icri][k]/_numberPixelsPerX[k]>(recppm+epsilon))
                {
                    xppm=xmin+k;
                    recppm=_tabX[icri][k]/_numberPixelsPerX[k];
                }
                _yEbarPerX[k]=_totYEPerX[k]/_tabX[icri][k];
            }

        }
        double dStTime = (double)xmin*(double)MsPerX;
        double dDur = (double)(xmax-xmin)*(double)MsPerX;
        oParam[StTime] = dStTime;
        oParam[Dur] = dDur;
        int criprec=icri;
        if(icri>0)
            for(int jcri=icri-1;jcri>=0;jcri--)
            {
                double ddcomp = dStTime + epsilon - (double)(_paramsArray[jcri][SH][StTime]) - (double)(_paramsArray[jcri][SH][Dur]);
                if(ddcomp > 0)
                {
                    criprec=jcri;
                    break;
                }
            }

        if(criprec==icri) oParam[PrevSt] = 9999.0f;
        else oParam[PrevSt] = oParam[StTime] - _paramsArray[criprec][SH][StTime];

        oParam[Fmax]=(double)ymax*KhzPerY;
        oParam[Fmin]=(double)ymin*KhzPerY;
        oParam[BW]=oParam[Fmax]-oParam[Fmin];
        oParam[FreqMP] = (double)MasterPoints.at(icri).y()*KhzPerY;
        oParam[PosMP] = (double)(MasterPoints.at(icri).x()-xmin+0.5f)/(double)(xmax-xmin+1.0f);
        oParam[FreqPkS]  = yfds * KhzPerY;
        oParam[FreqPkM] = yfdm * KhzPerY;
        oParam[PosPkS]   = (double)(xpps-xmin+0.5f)/(double)(xmax-xmin+1.0f);
        oParam[PosPkM]  = (double)(xppm-xmin+0.5f)/(double)(xmax-xmin+1.0f);
        oParam[FreqPkS2] = (double) _yEmaxPerX[icri][xpps-xmin] * KhzPerY;
        oParam[FreqPkM2] = (double) _yEmaxPerX[icri][xppm-xmin] * KhzPerY;
        oParam[PrevMP1] = 9999.0f;
        oParam[PrevMP2] = 9999.0f;
        oParam[NextMP1] = 9999.0f;
        oParam[NextMP2] = 9999.0f;
        double famp[4];
        double namp[4];
        for(int j=0;j<4;j++)
        {
            famp[j] = 0.0d;
            namp[j] = 0.0d;
        }
        double prorata,bdeb,bfin,larco,lartot;
        int ideb,ifin;
        //if(_detec->IDebug) _detec->_logText << "-6" << endl ;  //+++
        lartot=(double)(xmax-xmin+1);

        for(int k=0;k<xmax-xmin+1;k++)
        {
            bdeb=(((double)k)*4.0d)/lartot;
            bfin=(((double)(k+1))*4.0d)/lartot;
            ideb=(int)bdeb;
            ifin=(int)bfin;
            larco=bfin-bdeb;
            prorata = larco/lartot;
            for(int j=ideb;j<=ifin;j++)
            {
                if(j==4) break;
                if(j==ideb)
                {
                    if(ifin==ideb) prorata=1.0d/larco;
                    else prorata = ((double)(j+1)-bdeb)/larco;
                }
                else
                {
                    if(bfin < (j+1)) prorata=(bfin-(double)j)/larco;
                    else prorata = 1.0d/larco;
                }
                famp[j]+=(double)_tabX[icri][k]*prorata;
                namp[j]+=prorata;
            }
        }
        if(namp[0]>0) oParam[Amp1] = famp[0]*KhzPerY/namp[0]+epsilon;
        else oParam[Amp1] = 0.0f;
        if(namp[1]>0) oParam[Amp2] = famp[1]*KhzPerY/namp[1];
        else oParam[Amp2] = 0.0f;
        if(namp[2]>0) oParam[Amp3] = famp[2]*KhzPerY/namp[2];
        else oParam[Amp3] = 0.0f;
        if(namp[3]>0) oParam[Amp4] = famp[3]*KhzPerY/namp[3];
        else oParam[Amp4] = 0.0f;
        for(int inoise=0;inoise<4;inoise++)
        {
            int jdeb,jfin,x,y;
            int nps=0;
            double bruit=0.0d,bruit_moyen=0.0d;
            if(inoise<2) {jdeb=ymin;jfin=ymax;}
            else {jdeb=xmin;jfin=xmax;}
            for(int j=jdeb;j<=jfin;j++)
            {
                x=j;y=j;
                for(int k=1;k<=3;k++)
                {
                    if((inoise<2 && _xMaxPerY[j-ymin]>0) || (inoise>1 &&  _yMaxPerX[j-xmin]>0))
                    {
                        if(inoise==0) x=_xMinPerY[j-ymin]-k;
                        if(inoise==1) x=_xMaxPerY[j-ymin]+k;
                        if(inoise==2) y=_yMinPerX[j-xmin]-k;
                        if(inoise==3) y=_yMaxPerX[j-xmin]+k;
                        if(x>=0 && x<SonogramWidth && y>0 && y<FftHeightHalf)
                        {
                            bruit += (double)SonogramArray[y][x];
                            nps++;
                        }
                    }
                }
            }
            if(nps>0) bruit_moyen=bruit/ ((double)nps*100.0d);
            if(inoise==0)oParam[NoisePrev]  = bruit_moyen;
            if(inoise==1)oParam[NoiseNext] = bruit_moyen;
            if(inoise==2)oParam[NoiseDown]  = bruit_moyen;
            if(inoise==3)oParam[NoiseUp]    = bruit_moyen;
        }
        int nms=0;
        double sdm=0.0d;
        for(int k=0;k<xmax-xmin+1;k++)
        {
            if(_numberPixelsPerX[k]>0)
            {
                _averagePerX[k]=_tabX[icri][k]/_numberPixelsPerX[k];
                sdm+=_averagePerX[k];
                nms++;
            }
        }
        double mdm = sdm/nms;
        double sdc = 0.0d;
        for(int k=0;k<xmax-xmin+1;k++)
        {
            if(_numberPixelsPerX[k]>0) sdc += (_averagePerX[k]-mdm)*(_averagePerX[k]-mdm);
        }
        oParam[CVAmp]=pow(sdc/(double)nms,0.5);
        int tempsOuest[2];
        tempsOuest[0] = qMax(_xMinPerY[0],_xMinPerY[ymax-ymin])-xmin+1;
        oParamCrete[3][Dur] = (double)tempsOuest[0]*(double)MsPerX;
        tempsOuest[1] = qMax(_xSecondWestRidgePerY[0],_xSecondWestRidgePerY[ymax-ymin])-xmin+1;
        oParamCrete[4][Dur] = (double)tempsOuest[1]*(double)MsPerX;
        int fmin=_yEmaxPerX[icri][0],fmax=fmin;
        for(int k=1;k<=xmax-xmin;k++)
        {
            if(_yEmaxPerX[icri][k]>0)
            {
                if(_yEmaxPerX[icri][k]>fmax) fmax=_yEmaxPerX[icri][k];
                if(_yEmaxPerX[icri][k]<fmin) fmin=_yEmaxPerX[icri][k];
            }
        }
        oParamCrete[0][Fmax]=(double)fmax*(double)KhzPerY;
        oParamCrete[0][Fmin]=(double)fmin*(double)KhzPerY;
        oParamCrete[0][BW]=oParamCrete[0][Fmax]-oParamCrete[0][Fmin];
        quint16 *pc; 
        int alasuite,meilleuresuite,meilleur;
        double meilleurepente;
        int xmcmax[NCRETES],ymcmax[NCRETES];
        double amp[10],ampmoy[10],ran[10],ranmoy[10];
        int num[10],denom[10];
        num[0]=1; denom[0]=1;
        num[1]=2; denom[1]=1;
        num[2]=3; denom[2]=1;
        num[3]=3; denom[3]=2;
        num[4]=1; denom[4]=2;
        num[5]=4; denom[5]=3;
        num[6]=2; denom[6]=3;
        int np,y,xn,nn,yn;
        amp[0]=0.0d;

        for(int x=xmin;x<=xmax;x++) amp[0]+=_tabX[icri][x-xmin];
        ampmoy[0]=amp[0]/((double)tailleforme);
        for(int ih=1;ih<=6;ih++)
        {
            amp[ih]=0.0d;
            ampmoy[ih]=0.0d;
            np=0;
            ran[ih]=0.0d;
            ranmoy[ih]=0.0d;
            nn=0;
            for(int j=ymin;j<=ymax;j++)
            {
                if(ih>0) y=(j*num[ih])/denom[ih];
                if(y>_maxY || y<1) continue;
                for(int x=_xMinPerY[j-ymin];x<=_xMaxPerY[j-ymin];x++)
                {
                    //if(_tabYX[j-ymin][x-xmin]>0.0d)
                    if(_tabYX[j-ymin][x-xmin]>epsilon)
                    {
                        amp[ih] += (double)SonogramArray[y][x];
                        np++;
                    }
                }
                for(int gd=0;gd<2;gd++) for(int k=1;k<=3;k++)
                {
                    if(gd==0)xn=_xMinPerY[j-ymin]-k;
                    else xn=_xMaxPerY[j-ymin]+k;
                    if(xn>=0 && xn<SonogramWidth)
                    {
                        ran[ih] += (double)SonogramArray[y][xn];
                        nn++;
                    }
                }
            }
            for(int x=xmin;x<=xmax;x++)
            {
                for(int gd=0;gd<2;gd++) for(int k=1;k<=3;k++)
                {
                    if(_yMaxPerX[x-xmin]>0)
                    {
                        if(gd==0)yn=((_yMinPerX[x-xmin]*num[ih])/denom[ih])-k;
                        else yn=((_yMaxPerX[x-xmin]*num[ih])/denom[ih])+k;
                        if(yn<=_maxY && yn>0)
                        {
                            ran[ih] += (double)SonogramArray[yn][x];
                            nn++;
                        }
                    }
                }
            }
            if(np>0) ampmoy[ih]=amp[ih]/((double)np*100.0d);
            if(nn>0) ranmoy[ih]=ran[ih]/((double)nn*100.0d);
        } 
        if(ampmoy[0]>0.0d)
        {
            oParam[Ramp_2_1] = ampmoy[1]/ampmoy[0];
            oParam[Ramp_3_1] = ampmoy[2]/ampmoy[0];
            oParam[Ramp_3_2] = ampmoy[3]/ampmoy[0];
            oParam[Ramp_1_2] = ampmoy[4]/ampmoy[0];
            oParam[Ramp_4_3] = ampmoy[5]/ampmoy[0];
            oParam[Ramp_2_3] = ampmoy[6]/ampmoy[0];
        }
        else
        {
            oParam[Ramp_2_1] = 0.0f;
            oParam[Ramp_3_1] = 0.0f;
            oParam[Ramp_3_2] = 0.0f;
            oParam[Ramp_1_2] = 0.0f;
            oParam[Ramp_4_3] = 0.0f;
            oParam[Ramp_2_3] = 0.0f;
        }
        if(ranmoy[1]>0) oParam[RAN_2_1] =  ampmoy[1]/ranmoy[1];
        else oParam[RAN_2_1] =  0.0f;
        if(ranmoy[2]>0) oParam[RAN_3_1] =  ampmoy[2]/ranmoy[2];
        else oParam[RAN_3_1] =  0.0f;
        if(ranmoy[3]>0) oParam[RAN_3_2] =  ampmoy[3]/ranmoy[3];
        else oParam[RAN_3_2] =  0.0f;
        if(ranmoy[4]>0) oParam[RAN_1_2] =  ampmoy[4]/ranmoy[4];
        else oParam[RAN_1_2] =  0.0f;
        if(ranmoy[5]>0) oParam[RAN_4_3] =  ampmoy[5]/ranmoy[5];
        else oParam[RAN_4_3] =  0.0f;
        if(ranmoy[6]>0) oParam[RAN_2_3] =  ampmoy[6]/ranmoy[6];
        else oParam[RAN_2_3] =  0.0f;
        int xmaitr = MasterPoints.at(icri).x();
        int ymaitr = MasterPoints.at(icri).y();
        int npar,ntr,nen,xdeb,xfin,ydeb,yfin;
        double e1,e2,e3;
        for(int r=0;r<2;r++)
        {
            if(r==0) npar = HetX; else npar = HetXr;
            oParam[npar] = 0.0f;
            ntr = 0;
            nen = 0;
            double e1,e2,e3;
            if(r==0)
            {
                ydeb=ymin;
                yfin=ymax;
            }
            else
            {
                ydeb = ymaitr-3;
                if(ydeb < ymin) ydeb = ymin;
                yfin = ymaitr+3;
                if(yfin > ymax) yfin = ymax;
            }
            for(int y=ydeb;y<=yfin;y++)
            {
                for(int x=_xMinPerY[y-ymin];x<=_xMaxPerY[y-ymin];x++)
                {
                    e2=_tabYX[y-ymin][x-xmin];
                    // if(e2>0.0f && x>0 && x<SonogramWidth-1)
                    if(e2>epsilon && x>0 && x<SonogramWidth-1)
                    {
                        ntr++;
                        e1 = (double)SonogramArray[y][x-1]/100.0d;
                        e3 = (double)SonogramArray[y][x+1]/100.0d;
                        if((e2>(e1+epsilon) && e2>(e3+epsilon)) || ((e2+epsilon)<e1 && (e2+epsilon)<e3)) nen++;
                    }
                }
            }
            if(ntr>0) oParam[npar] = ((double) nen)/ ((double) ntr);
        }
        for(int r=0;r<3;r++)
        {
            double dift = 0.0d;
            if(r==0) npar = HetY;
            else
            {
                if(r==1) npar = HetYr;
                else npar = HetYr2;
            }
            oParam[npar] = 0.0d;
            ntr=0;
            nen=0;
            if(r==0 || r==2)
            {
                xdeb=xmin;
                xfin=xmax;
            }
            else
            {
                xdeb = xmaitr-3;
                if(xdeb < xmin) xdeb = xmin;
                xfin = xmaitr+3;
                if(xfin > xmax) xfin = xmax;
            }
            for(int x=xdeb;x<=xfin;x++)
            {
                if(_yMaxPerX[x-xmin]>0)
                {
                    if(r==0 || r==1)
                    {
                        ydeb=_yMinPerX[x-xmin];
                        yfin=_yMaxPerX[x-xmin];
                    }
                    else
                    {
                        ydeb = ymaitr-3;
                        if(ydeb < ymin) ydeb = ymin;
                        yfin = ymaitr+3;
                        if(yfin > ymax) yfin = ymax;
                    }
                    for(int y=ydeb;y<=yfin;y++)
                    {
                        e2=_tabYX[y-ymin][x-xmin];
                        //if(e2>0.0f && y>1 && y<_maxY)
                        if(e2>epsilon && y>1 && y<_maxY)
                        {
                            ntr++;
                            e1=(double)SonogramArray[y-1][x]/100.0d;
                            e3=(double)SonogramArray[y+1][x]/100.0d;
                            dift+=qAbs(e3-e2);
                            if((e2>(e1+epsilon) && e2>(e3+epsilon)) || ((e2+epsilon)<e1 && (e2+epsilon)<e3)) nen++;
                        }
                    }
                }
            }
            if(ntr>0)
            {
                if(npar<2) oParam[npar] = ((double) nen)/ ((double) ntr);
                else oParam[npar] = (double)dift/(double)ntr;
            }
        }
        oParam[HetCMC] = 0.0f;
        oParam[HetCMD] = 0.0f;
        oParam[HetCTC] = 0.0f;
        oParam[HetCTD] = 0.0f;
        //
        oParam[HetCMnP] = 0.0f;
        oParam[HetCMfP] = 0.0f;
        oParam[HetCTnP] = 0.0f;
        oParam[HetCTfP] = 0.0f;
        //
        int pics[2][(MAXLARCRI+1)/2];
        int nbpics[2];
        int modepc;
        int memopr;
        double seuilpc[2];
        seuilpc[0] = 1.0d; seuilpc[1] = 10.0d;
        double valpr,difppr;
        for(int mt=0;mt<2;mt++)
        {
            int sensd=1,dersensd = 0;
            int nch=0,ncc=0;
            double moytot= 0.0d,dermoytot = 0.0d;
            double difC=0.0d,totDifC=0.0d,derDifC=0.0d;
            nbpics[mt]=0;
            modepc=0;
            memopr=0;
            valpr=0.0d;
            for(int x=xmin;x<=xmax;x++)
            {
                if(_numberPixelsPerX[x-xmin]>0)
                {                    if(mt==0) moytot = _averagePerX[x-xmin]; else moytot= _tabX[icri][x-xmin];
                    if(x==xmin) valpr = moytot;
                }
                else moytot = 0.0d;
                difppr=moytot-valpr;
                if(difppr>epsilon)
                {
                    if(difppr>=(seuilpc[mt]-epsilon))
                    {
                        if(modepc==0 || modepc==2)
                        {
                            modepc=1;
                        }
                        memopr=x-xmin;
                        valpr = moytot;
                    }
                    else
                    {
                        if(modepc==1) {memopr=x-xmin; valpr=moytot;}
                    }
                }
                else
                {
                    if(difppr<-epsilon)
                    {
                        if(qAbs(difppr)>=(seuilpc[mt]-epsilon))
                        {
                            if(modepc==0 || modepc==1)
                            {
                                pics[mt][nbpics[mt]++]=memopr;

                                modepc=2;
                            }
                            memopr=x-xmin; valpr=moytot;
                        }
                        else
                        {
                            if(modepc==2) {memopr=x-xmin; valpr=moytot;}
                        }
                    }
                }
                if(x==xmax)
                {
                    if(modepc==1)
                    {
                        pics[mt][nbpics[mt]++] = memopr;
                    }
                }
                difC = moytot - dermoytot;
                if(difC>epsilon) sensd = 1;
                else
                {
                    if(difC < -epsilon) sensd = 0;
                }
                //if(difC!=0 || derDifC!=0) ncc++;
                if(qAbs(difC)>epsilon || qAbs(derDifC)>epsilon) ncc++;
                if(sensd != dersensd) nch++;
                totDifC+=qAbs(difC);
                dermoytot = moytot;
                dersensd = sensd;
                derDifC = difC;
            }
            if(ncc>0)
            {
                if(mt==0)
                {
                    oParam[HetCMC] = ((double) nch)/ ((double) ncc);

                    oParam[HetCMD] = (totDifC)/ ((double) ncc);
                }
                else
                {
                    oParam[HetCTC] = ((double) nch)/ ((double) ncc);

                    oParam[HetCTD] = (totDifC)/ ((double) ncc);
                }

            }
        }
        oParam[HetCMnP] = nbpics[0];
        oParam[HetCMfP] = (double)nbpics[0]/((xmax-xmin+1)*MsPerX);
        oParam[HetCTnP] = nbpics[1];
        oParam[HetCTfP] = (double)nbpics[1]/((xmax-xmin+1)*MsPerX);
        double interc[MAXLARCRI/2];
        double variationPicsInter[MAXLARCRI/2];

        for(int mt=0;mt<2;mt++)
        {
            double medianPicsDistance = ((double)(xmax-xmin+2)/2)*MsPerX;
            double medianPicsLittleDistance =((double)(xmax-xmin+2)/2)*MsPerX;
            double medianPicsBigDistance = ((double)(xmax-xmin+2)/2)*MsPerX;
            int nbi = nbpics[mt]-1;
            int halfnbi = nbi/2;
            int quartnbi =halfnbi/2;
            if(nbi>1)
            {
                for(int j=0;j<nbi;j++)
                {
                    interc[j] = (qAbs((pics[mt][j+1]-pics[mt][j])))*MsPerX;
                }
                sortDoubleArray(interc,nbi);
                if(halfnbi*2<nbi) medianPicsDistance = interc[halfnbi];
                else medianPicsDistance= (interc[halfnbi-1] + interc[halfnbi])/2;
                if(nbi<4)
                {
                    medianPicsLittleDistance = interc[0];
                    medianPicsBigDistance = interc[nbi-1];
                }
                else
                {
                    if(quartnbi*2<halfnbi)
                    {
                        medianPicsLittleDistance = interc[quartnbi];
                        medianPicsBigDistance = interc[nbi-1-quartnbi];
                    }
                    else
                    {
                        medianPicsLittleDistance = (interc[quartnbi-1] + interc[quartnbi])/2;
                        medianPicsBigDistance = (interc[nbi-1-quartnbi] + interc[nbi-quartnbi])/2;
                    }
                }
            }
            if(mt==0)
            {
                _paramsArray[icri][SH][HetPicsMAD]   = medianPicsDistance;
                _paramsArray[icri][SH][HetPicsMALD]  = medianPicsLittleDistance;
                _paramsArray[icri][SH][HetPicsMABD]  = medianPicsBigDistance;
                _paramsArray[icri][SH][HetPicsMRBLD] = medianPicsBigDistance / medianPicsLittleDistance;
            }
            else
            {
                _paramsArray[icri][SH][HetPicsTAD]   = medianPicsDistance;
                _paramsArray[icri][SH][HetPicsTALD]  = medianPicsLittleDistance;
                _paramsArray[icri][SH][HetPicsTABD]  = medianPicsBigDistance;
                _paramsArray[icri][SH][HetPicsTRBLD] = medianPicsBigDistance / medianPicsLittleDistance;
            }
            double medianPicsDistanceVariation = 0;
            if(nbi==1)
            {


            }
            for(int j=0;j<nbi;j++) variationPicsInter[j] = qAbs(interc[j]-medianPicsDistance);

            double medianPicsLittleDistanceVariation = 0.0d;
            double medianPicsBigDistanceVariation = 0.0d;
            if(nbi>1)
            {
                medianPicsLittleDistanceVariation = variationPicsInter[0];
                medianPicsBigDistanceVariation = variationPicsInter[nbi-1];
            }
            if(nbi>1)
            {
                sortDoubleArray(variationPicsInter,nbi);
                if(halfnbi*2<nbi) medianPicsDistanceVariation = variationPicsInter[halfnbi];
                else medianPicsDistanceVariation = (variationPicsInter[halfnbi-1] + variationPicsInter[halfnbi])/2;
                for(int j=0;j<nbi;j++)
                {
                    if(j<halfnbi)
                    {
                        variationPicsInter[j] = qAbs(interc[j]-medianPicsLittleDistance);
                    }
                    else
                    {
                        variationPicsInter[j] = qAbs(interc[j]-medianPicsBigDistance);
                    }
                }
                if(nbi<4)
                {
                    medianPicsLittleDistanceVariation = variationPicsInter[0];
                    medianPicsBigDistanceVariation = variationPicsInter[nbi-1];
                }
                else
                {
                    sortDoubleArray(variationPicsInter,halfnbi);
                    sortDoubleArray(variationPicsInter+nbi-halfnbi,halfnbi);
                    if(quartnbi*2<halfnbi)
                    {
                        medianPicsLittleDistanceVariation = variationPicsInter[quartnbi];
                        medianPicsBigDistanceVariation = variationPicsInter[nbi-1-quartnbi];
                    }
                    else
                    {
                        medianPicsLittleDistanceVariation = (variationPicsInter[quartnbi-1] + variationPicsInter[quartnbi])/2;
                        medianPicsBigDistanceVariation = (variationPicsInter[nbi-1-quartnbi] + variationPicsInter[nbi-quartnbi])/2;
                    }
                }
            }
            if(mt==0)
            {
                _paramsArray[icri][SH][VDPicsM]   = medianPicsDistanceVariation;
                _paramsArray[icri][SH][VLDPicsM]  = medianPicsLittleDistanceVariation;
                _paramsArray[icri][SH][VBDPicsM]  = medianPicsBigDistanceVariation;
                _paramsArray[icri][SH][VDPPicsM]  = medianPicsDistanceVariation/medianPicsDistance;
                _paramsArray[icri][SH][VLDPPicsM]  = medianPicsLittleDistanceVariation/medianPicsLittleDistance;
                _paramsArray[icri][SH][VBDPPicsM]  = medianPicsBigDistanceVariation/medianPicsBigDistance;
            }
            else
            {
                _paramsArray[icri][SH][VDPicsT]   = medianPicsDistanceVariation;
                _paramsArray[icri][SH][VLDPicsT]  = medianPicsLittleDistanceVariation;
                _paramsArray[icri][SH][VBDPicsT]  = medianPicsBigDistanceVariation;
                _paramsArray[icri][SH][VDPPicsT]  = medianPicsDistanceVariation/medianPicsDistance;
                _paramsArray[icri][SH][VLDPPicsT]  = medianPicsLittleDistanceVariation/medianPicsLittleDistance;
                _paramsArray[icri][SH][VBDPPicsT]  = medianPicsBigDistanceVariation/medianPicsBigDistance;
            }
        }
        // -----------------------------------------------------------------------------------
        //if(_detec->IDebug) _detec->_logText << "-8" << endl ;  //+++
        oParam[Dbl8] = 0.0f;
        double e8 = 0.0d;
        int n8 = 0;
        int Y8 = (int)(_freqCallMin / KhzPerY);
        for(int y=0;y<Y8;y++)
        {
            qint16 *fc=SonogramArray[y];
            for(int x=xmin;x<=xmax;x++)
            {
                e8 =+ (double)(*fc++);
                n8++;
            }
        }
        if(n8>0) oParam[Dbl8] = eMoy - (e8/((double)n8*100.0d));
        oParam[Stab] = 0.0f;
        double dif,p;
        double ponderTot = 0.0d;
        double difTot = 0.0d;
        double dxmax = ((double)qMax(xmaitr-xmin,xmax-xmaitr))*MsPerX;
        double dymax = ((double)qMax(ymaitr-ymin,ymax-ymaitr))*KhzPerY;
        double distLim = pow(dxmax,2)+pow(dymax,2);
        //
        if(distLim>epsilon)
        {
            for(int y=ymin;y<=ymax;y++)
            {
                for(int x=_xMinPerY[y-ymin]-1;x<=_xMaxPerY[y-ymin];x++)
                {
                    if(x>=0 && x<SonogramWidth-1)
                    {
                        dif = (double)(qAbs(SonogramArray[y][x+1]-SonogramArray[y][x]));
                        p = 1.0d - (((    pow(   ((double)(x-xmaitr))*MsPerX , 2)
                                          +pow(   ((double)(y-ymaitr))*KhzPerY,2))/distLim) * 0.75f);
                        difTot += dif * p;
                        ponderTot += p;
                    }
                }
            }
            for(int x=xmin;x<=xmax;x++)
            {
                if(_yMaxPerX[x-xmin]>0)
                {
                    for(int y=_yMinPerX[x-xmin]-1;y<=_yMaxPerX[x-xmin];y++)
                    {
                        if(y>=_minY && y<_maxY-1)
                        {
                            dif = (double)(qAbs(SonogramArray[y+1][x]-SonogramArray[y][x]));
                            p = 1.0d - (((    pow(   ((double)(x-xmaitr))*MsPerX , 2)
                                              +pow(   ((double)(y-ymaitr))*KhzPerY,2))/distLim) * 0.75f);
                            difTot += dif * p;
                            ponderTot += p;
                        }
                    }
                }
            }
            if(ponderTot!=0) oParam[Stab] = difTot / (ponderTot*100.0d);
        }
        oParam[EnStabSm] = 0.0f;
        oParam[EnStabLg] = 0.0f;
        int radius,nbp;
        //
        for(int lbr=0;lbr<2;lbr++)
        {
            difTot = 0.0d;
            nbp=0;
            if(lbr==0) radius = 3; else radius = 10;
            ydeb = ymaitr-radius; yfin = ymaitr + radius;
            xdeb = xmaitr - radius ; xfin = xmaitr + radius;
            if(ydeb < _minY) ydeb = _minY;
            if(yfin > _maxY-1) yfin = _maxY-1;
            if(xdeb < 0) xdeb = 0;
            if(xfin > SonogramWidth-1) xfin = SonogramWidth-1;
            for(int y=ydeb;y<=yfin;y++)
            {
                for(int x=xdeb;x<xfin;x++)
                {
                    difTot += (double)(qAbs(SonogramArray[y][x+1]-SonogramArray[y][x]));
                    nbp++;
                }
            }
            for(int x=xdeb;x<=xfin;x++)
            {
                for(int y=ydeb;y<yfin;y++)
                {
                    difTot += (double)(qAbs(SonogramArray[y+1][x]-SonogramArray[y][x]));
                    nbp++;
                }
            }
            if(lbr==0) npar = EnStabSm; else npar = EnStabLg;
            if(nbp>0) oParam[npar] = difTot/((double)nbp*100.0d);
        }
        double enerMaster[2];
        enerMaster[0] = (double) _tabY[ymaitr-ymin];
        enerMaster[1] = (double) _tabY[ymaitr-ymin]/_numberPixelsPerY[ymaitr-ymin];
        int enerHeight[2][2];
        double enerY[2];
        int incr;
        for(int tm=0;tm<2;tm++)
        {
            for(int sens = 0;sens <=1;sens++)
            {
                enerHeight[tm][sens] = 0;
                if(sens) incr = 1; else incr = -1;
                for(int y=ymaitr+incr;y>=ymin && y<=ymax;y+=incr)
                {
                    if(tm==0) enerY[tm] = (double)_tabY[y-ymin];
                    else enerY[tm] = (double) _tabY[y-ymin] / _numberPixelsPerY[y-ymin];
                    if(enerY[tm] > (enerMaster[tm] * 0.8f + epsilon)) enerHeight[tm][sens]++;
                    else break;
                }
            }
        }
        oParam[HeiET] = (double) (enerHeight[0][0]+enerHeight[0][1]);
        oParam[HeiEM] = (double) (enerHeight[1][0]+enerHeight[1][1]);
        oParam[HeiRT] = oParam[HeiET]/(ymax - ymin+1);
        oParam[HeiRM] = oParam[HeiEM]/(ymax - ymin+1);
        oParam[HeiETT] = (double) enerHeight[0][1];
        oParam[HeiEMT] = (double) enerHeight[1][1];
        oParam[HeiRTT] = oParam[HeiETT]/(ymax - ymin+1);
        oParam[HeiRMT] = oParam[HeiEMT]/(ymax - ymin+1);
        if(_imageData)
        {
            _callMasterRidge.clear();
            _callSouthRidge.clear();
            _callNorthRidge.clear();
            _callWestRidge.clear();
            _callSecondWestRidge.clear();
        }
        int lowfc,hifc,fc3;
        double yCmax=0.0d;
        double yCmin=9999.0d;
        for(int jcrete=0;jcrete<NCRETES;jcrete++)
        {
            if(jcrete==0)  pc=_yEmaxPerX[icri];
            if(jcrete==1)  pc=_yMinPerX;
            if(jcrete==2)  pc=_yMaxPerX;
            if(jcrete==3)  pc=_xMinPerY;
            if(jcrete==4)  pc=_xSecondWestRidgePerY;
            oParamCrete[jcrete][Slope] = 9999.0f;
            double emcmax = -50.0d;
            xmcmax[jcrete] = 0;
            ymcmax[jcrete]=0;
            if(jcrete<3)
            {
                alasuite = 0; meilleuresuite = 0; meilleur = xmin; meilleurepente = 100000.0d;
                yCmax=0.0d;
                yCmin=9999.0d;
                for(int x=xmin;x<=xmax;x++)
                {
                    if(pc[x-xmin]>yCmax)  yCmax = pc[x-xmin];
                    if(pc[x-xmin]<yCmin && pc[x-xmin]>0)  yCmin = pc[x-xmin];
                    if(_imageData)
                    {
                        if(pc[x-xmin]>0)
                        {
                            QPoint p;
                            p.setX(x);
                            p.setY(pc[x-xmin]);
                            if(jcrete==0) _callMasterRidge.push_back(p);
                            if(jcrete==1) _callSouthRidge.push_back(p);
                            if(jcrete==2) _callNorthRidge.push_back(p);
                        }
                    }
                    _slope[x-xmin]=100000.0;
                    if(x>xmin+1 && x < xmax-1)
                        if(pc[x-xmin+2]>0 && pc[x-xmin-2]>0)
                            _slope[x-xmin] = qAbs(((double)pc[x-xmin+2] - (double)pc[x-xmin-2])/4.0d);
                    if(_slope[x-xmin]==0.0d)
                    {
                        meilleurepente = 0;
                        alasuite++;
                        if(alasuite > 1 && (meilleur == xmin || alasuite > meilleuresuite))
                        {meilleur = x; meilleuresuite = alasuite;}
                    }
                    else
                    {
                        alasuite = 0;
                        if(_slope[x-xmin] < (meilleurepente-epsilon))
                        {
                            meilleur = x;
                            meilleurepente = _slope[x-xmin];
                        }
                    }
                    double e=(double)SonogramArray[pc[x-xmin]][x]/100.0d;
                    if(e>(emcmax+epsilon))
                    {
                        emcmax = e;
                        xmcmax[jcrete]=x;
                        ymcmax[jcrete] = pc[x-xmin];
                    }
                }
                oParamCrete[jcrete][Fmax] = yCmax * KhzPerY;
                oParamCrete[jcrete][Fmin]  = yCmin  * KhzPerY;
                oParamCrete[jcrete][BW]  = oParamCrete[jcrete][Fmax]-oParamCrete[jcrete][Fmin];
                if(xmax>xmin) oParamCrete[jcrete][Slope] = (double)(pc[xmax-xmin]-pc[0])*KhzPerY/((double)(xmax-xmin)*MsPerX);
                if(meilleurepente ==0 && meilleuresuite>1)
                {
                    int bonmeilleur = meilleur;
                    for(int k=1;k<meilleuresuite/2;k++)
                        if(pc[meilleur-k-xmin]>0) bonmeilleur = meilleur-k;
                    meilleur = bonmeilleur;
                }
                lowfc = meilleur-xmin;
                LowSlope[(icri*NCRETES+jcrete)*2] = meilleur; // (xmin+lowfc)
                LowSlope[(icri*NCRETES+jcrete)*2+1] = pc[lowfc];
                oParamCrete[jcrete][FIF]= pc[lowfc]*KhzPerY;
                double pe1=0.0d,pe2=0.0d,dif;
                double meilleuredist;
                int meilleurk[2];
                int debb,finb;
                meilleurk[0]=0;
                meilleurk[1]=lowfc;
                for(int jp=0;jp<2;jp++)
                {
                    meilleurk[jp]=0;
                    meilleuredist = 1000000.0d;
                    if(jp==0) {debb=0;finb=lowfc;}
                    else {debb=lowfc;finb=xmax-xmin;}
                    if((jp==0 && lowfc > 1) || (jp==1 && lowfc < xmax - xmin-1))
                    {
                        for (int k = debb+1 ; k<= finb-1 ; k++)
                        {
                            dif=0.0d;
                            if(pc[k]>0)
                            {
                                pe1=(double)(pc[k]-pc[debb])/((double)k);
                                pe2=((double)(pc[finb]-pc[k]))/((double)(lowfc-k));
                            }
                            for(int s = debb+1 ; s<k ; s++)
                                if(pc[s]>0)
                                    dif+=(double)pow(pe1*((double)s)+(double)pc[debb]-(double)pc[s],2);

                            for(int s = k+1 ; s<=finb-1 ; s++)
                                if(pc[s]>0)
                                    dif+=(double)pow(pe2*((double)(finb-s))+(double)pc[s]-(double)pc[finb],2);
                            if(dif<(meilleuredist-epsilon))
                            {
                                meilleurk[jp]=k;
                                meilleuredist=dif;
                            }
                        }

                    }
                }
                hifc=meilleurk[0];
                Inflexion1[(icri*NCRETES+jcrete)*2]=hifc+xmin;
                Inflexion1[(icri*NCRETES+jcrete)*2+1]=pc[hifc];
                oParamCrete[jcrete][HCF] = pc[hifc]*KhzPerY;
                oParamCrete[jcrete][THCF]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][THCF] = ((double)meilleurk[0]+0.5f)/((double)(xmax-xmin+1));
                fc3=meilleurk[1];
                Inflexion3[(icri*NCRETES+jcrete)*2]=fc3+xmin;
                Inflexion3[(icri*NCRETES+jcrete)*2+1]=pc[fc3];
                oParamCrete[jcrete][LCF] = pc[fc3]*KhzPerY;

                oParamCrete[jcrete][UpSl] = 0.0f;
                if(meilleurk[0]>0)
                    oParamCrete[jcrete][UpSl] = ((double)(pc[meilleurk[0]]-pc[0])*KhzPerY)/((double)meilleurk[0]*MsPerX);
                oParamCrete[jcrete][LoSl] = 0.0f;
                if(lowfc-meilleurk[0]>0)
                    oParamCrete[jcrete][LoSl] = ((double)(pc[lowfc]-pc[meilleurk[0]])*KhzPerY)/((double)(lowfc-meilleurk[0])*MsPerX);
                int debp,finp,ecart;
                double pj[4];
                for(int j2=0;j2<4;j2++)
                {
                    pj[j2]=9999.0d;
                    if(j2<2)
                    {
                        ecart=(xmax-xmin+10)/20;
                        if(ecart==0 && xmax>xmin) ecart = 1;
                        if(j2==0) {debp=0; finp=ecart;}
                        else  {debp=xmax-xmin-ecart; finp=xmax-xmin;}
                    }
                    else
                    {
                        if(j2==2)
                        {
                        debp=xmcmax[jcrete]-xmin-2;
                        finp=xmcmax[jcrete]-xmin+2;
                        }
                        else
                        {
                            debp = lowfc-2;
                            finp = lowfc+2;
                        }
                        if(debp<0) debp =0;
                        if(finp>xmax-xmin) finp = xmax-xmin;

                        ecart=finp-debp;
                    }
                    if(ecart>0)
                    {
                        pj[j2]=((double)(pc[finp]-pc[debp])*KhzPerY)/((double)ecart*MsPerX);
                    }
                }
                oParamCrete[jcrete][StSl] = pj[0];
                oParamCrete[jcrete][EnSl]   = pj[1];
                oParamCrete[jcrete][FISl]  = pj[3];
                oParamCrete[jcrete][FPSl]  = pj[2];
                //
                if(_paramVersion>=1)
                {
                    oParamCrete[jcrete][SDC]  = 0;
                    oParamCrete[jcrete][SDCR]  = 0;
                    oParamCrete[jcrete][SDCL]  = 0;
                    oParamCrete[jcrete][SDCLR]  = 0;
                    oParamCrete[jcrete][SDCLRY]  = 0;
                    oParamCrete[jcrete][SDCLRXY]  = 0;
                    oParamCrete[jcrete][SDCLRXY2]  = 0;
                    oParamCrete[jcrete][SDCLOP]  = 0;
                    oParamCrete[jcrete][SDCLROP]  = 0;
                    oParamCrete[jcrete][SDCLRYOP]  = 0;
                    oParamCrete[jcrete][SDCLRXYOP]  = 0;
                    oParamCrete[jcrete][SDCLWB]  = 0;
                    oParamCrete[jcrete][SDCLRWB]  = 0;
                    oParamCrete[jcrete][SDCLRYWB]  = 0;
                    oParamCrete[jcrete][SDCLRXYWB]  = 0;
                    oParamCrete[jcrete][SDCLOPWB]  = 0;
                    oParamCrete[jcrete][SDCLROPWB]  = 0;
                    oParamCrete[jcrete][SDCLRYOPWB]  = 0;
                    oParamCrete[jcrete][SDCLRXYOPWB]  = 0;
                    oParamCrete[jcrete][SDCL_DNP]  = 0;
                    oParamCrete[jcrete][SDCLR_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRY_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRXY_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRXY2_DNP]  = 0;
                    double derdif1=0.0d;
                    double dif1,dif2;
                    double totdif2 = 0.0d;
                    double totdif2Op = 0.0d;
                    double totdif2wb = 0.0d;
                    double totdif2opwb = 0.0d;
                    // -----------------------------------------------------------------------------------------
                    bool findSdcl = false;
                    bool findSdclOp = false;
                    bool seenMasterPoint = false;
                    int prem = -1;
                    int xdOp = xmin;
                    int ncs = 1;
                    bool firstdif2 = true;
                    for(int x=xmin+1;x<=xmax;x++)
                    {
                        if(x==xmaitr) seenMasterPoint = true;
                        if(pc[x-xmin]>0 && pc[x-xmin-1]>0) dif1 = pc[x-xmin]-pc[x-xmin-1];
                        else continue;
                        if(dif1 * derdif1 < 0) ncs++;
                        if(dif1>1 && !seenMasterPoint && !findSdclOp)
                        {
                            totdif2Op = 0.0d;
                            totdif2opwb = 0.0d;
                            xdOp=x;
                        }
                        if((jcrete == (CS-1) || jcrete == (CM-1)) && (dif1>1 || x==xmax) && findSdcl==false)
                        {
                            prem = x;
                            oParamCrete[jcrete][SDCL]  = totdif2*KhzPerY;
                            oParamCrete[jcrete][SDCLR]  = (   ((double)totdif2) * KhzPerY) / ( ((double)(x-xmin)+0.5f)*MsPerX);
                            oParamCrete[jcrete][SDCLWB]  = totdif2wb * KhzPerY;
                            oParamCrete[jcrete][SDCLRWB]  = (((double)totdif2wb)*KhzPerY) / (((double)(x-xmin)+0.5f)*MsPerX);
                            if(qAbs(pc[x-xmin]-pc[0])>0)
                            {
                                double ratioY = (totdif2*10.0d)/(double)(qAbs(pc[x-xmin]-pc[0]));
                                oParamCrete[jcrete][SDCLRY]  = ratioY/MsPerX;
                                oParamCrete[jcrete][SDCLRXY]  = ratioY / ((double)(x-xmin)+0.5f);
                                oParamCrete[jcrete][SDCLRXY2]  = ratioY / ((double)(x-xmin)+0.5f);
                                //
                                double ratioYwb = (totdif2wb*10.0d)/(double)(qAbs(pc[x-xmin]-pc[0]));
                                oParamCrete[jcrete][SDCLRYWB]  = ratioYwb / MsPerX;
                                oParamCrete[jcrete][SDCLRXYWB]  = ratioYwb / ((double)(x-xmin)+0.5f);
                            }
                            findSdcl = true;
                        }
                        if((jcrete == (CS-1) || jcrete == (CM-1)) && (dif1>1 || x==xmax) && findSdclOp==false && seenMasterPoint)
                        {
                            oParamCrete[jcrete][SDCLOP]  = totdif2Op*KhzPerY;
                            oParamCrete[jcrete][SDCLROP]  = (((double)totdif2Op)*KhzPerY) / (((double)(x-xdOp)+0.5f)*MsPerX);
                            oParamCrete[jcrete][SDCLOPWB]  = totdif2opwb * KhzPerY;
                            oParamCrete[jcrete][SDCLROPWB]  = (((double)totdif2opwb)*KhzPerY) / (((double)(x-xdOp)+0.5f)*MsPerX);
                            if(qAbs(pc[x-xmin]-pc[xdOp-xmin])>0)
                            {
                                double ratioYOp = (totdif2Op*100.0d)/(double)(qAbs(pc[x-xmin]-pc[xdOp-xmin]));
                                oParamCrete[jcrete][SDCLRYOP]  = ratioYOp/MsPerX;
                                oParamCrete[jcrete][SDCLRXYOP]  = ratioYOp / ((double)(x-xdOp)+0.5f);
                                double ratioYopwb = (totdif2opwb*10.0d)/(double)(qAbs(pc[x-xmin]-pc[xdOp-xmin]));
                                oParamCrete[jcrete][SDCLRYOPWB]  = ratioYopwb/MsPerX;
                                oParamCrete[jcrete][SDCLRXYOPWB]  = ratioYopwb / ((double)(x-xdOp)+0.5f);
                            }
                            findSdclOp = true;
                        }
                        dif2 = qAbs(dif1-derdif1);
                        if(firstdif2) {dif2 = 0; firstdif2=false;}

                        totdif2 += dif2;

                        totdif2Op += dif2;
                        if(x>=(xmin+xmaitr)/2) totdif2wb += dif2;
                        if(x>=(xdOp+xmaitr)/2) totdif2opwb += dif2;
                        derdif1 = dif1;
                    }
                    oParamCrete[jcrete][SDC]  = totdif2 * KhzPerY;
                    oParamCrete[jcrete][SDCR]  = ((double)totdif2*KhzPerY) / (((double)(xmax-xmin)+0.5f)*MsPerX);
                    double labs = qAbs(pc[xmax-xmin]-pc[0])+1;

                    double ratioY = (double)totdif2*10.0d/labs;
                    oParamCrete[jcrete][SDCRY]  = ratioY/MsPerX;


                    oParamCrete[jcrete][SDCRXY]  = ratioY / ((double)(xmax-xmin)+0.5f);
                    double fncs = (double) ncs/100.0d;
                    oParamCrete[jcrete][SDCL_DNP]  = oParamCrete[jcrete][SDCL]/fncs;
                    oParamCrete[jcrete][SDCLR_DNP]  = oParamCrete[jcrete][SDCLR]/fncs;
                    oParamCrete[jcrete][SDCLRY_DNP]  = oParamCrete[jcrete][SDCLRY]/fncs;
                    oParamCrete[jcrete][SDCLRXY_DNP]  = oParamCrete[jcrete][SDCLRXY]/fncs;
                    oParamCrete[jcrete][SDCLRXY2_DNP]  = oParamCrete[jcrete][SDCLRXY2]/fncs;
                    // -----------------------------------------------------------------------------------------
                    int em=3;
                    oParamCrete[jcrete][ELBPOS]  = 9999;
                    oParamCrete[jcrete][ELBSB]  = 0;
                    for(int x=xmin+em;x<=xmax-em;x++)
                    {
                        if(x>=prem || (xmax-xmin) < (em*2)) break;
                        int x0 = x-em,x2=x+em;
                        int y0=pc[x0-xmin],y1=pc[x-xmin],y2=pc[x2-xmin];
                        if(y1>0 && y0>0 && y2>0)
                        {
                            bool unecond = false;
                            if(y2>y1) unecond = true;
                            else
                            {
                                if(y1!=y0)
                                {
                                    if ((((double)(y2-y1)) / ((double)(y1-y0))) < 0.6) unecond = true;
                                }
                            }
                            if(unecond)
                            {
                                if(prem>xmin)
                                {
                                    oParamCrete[jcrete][ELBPOS]  = ((double)(x-xmin))/((double)(prem-xmin));
                                }
                                if(x>xmin)
                                {
                                    oParamCrete[jcrete][ELBSB]  = (   (  (double)(y1-pc[0])  ) *KhzPerY)     / ( ((double)(x-xmin))*MsPerX );
                                }
                                break;
                            }
                        }
                    }
                    //
                    oParamCrete[jcrete][ELB2POS]  = 9999;
                    oParamCrete[jcrete][ELB2SB]  = 0;
                    int yd=pc[0];
                    for(int x=xmin+em;x<=xmax-em;x++)
                    {
                        if(x>=prem || (xmax-xmin) < (em*2)) break;
                        int x2=x+em;
                        int y1=pc[x-xmin],y2=pc[x2-xmin];
                        if(y1>0 && y2>0)
                        {
                            bool unecond = false;
                            if(y2>y1)
                            {
                                unecond = true;
                            }
                            else
                            {
                                if(y1>yd) if  ( (((double)(y2-y1)) / ((double)(y1-yd))) < 0.6)
                                {
                                    unecond = true;
                                }
                            }
                            if(unecond)
                            {
                                if(prem>xmin)
                                {
                                    oParamCrete[jcrete][ELB2POS]  = ((double)(x-xmin))/((double)(prem-xmin));
                                }
                                if(x>xmin)
                                {
                                    oParamCrete[jcrete][ELB2SB]  = (  ( (double)(y1-pc[0]))*KhzPerY)      / ( ((double)(x-xmin))*MsPerX );
                                }
                                break;
                            }
                        }
                    }
                    // -----------------------------------------------------------------------------------------
                    if(jcrete==0)
                    {
                        oParamCrete[jcrete][RAF]  = 0;
                        oParamCrete[jcrete][RAE]  = 0;
                        oParamCrete[jcrete][RAFE]  = 0;
                        oParamCrete[jcrete][RAFP]  = 0;
                        oParamCrete[jcrete][RAFP2]  = 0;
                        oParamCrete[jcrete][RAFP3]  = 0;
                        oParamCrete[jcrete][SBMP]  = 0;
                        oParamCrete[jcrete][SAMP]  = 0;
                        oParamCrete[jcrete][SBAR]  = 0;
                        oParamCrete[jcrete][RAHP2]  = 0;
                        oParamCrete[jcrete][RAHP4]  = 0;
                        oParamCrete[jcrete][RAHP8]  = 0;
                        oParamCrete[jcrete][RAHP16]  = 0;
                        oParamCrete[jcrete][RAHE2]  = 0;
                        oParamCrete[jcrete][RAHE4]  = 0;
                        oParamCrete[jcrete][RAHE8]  = 0;
                        oParamCrete[jcrete][RAHE16]  = 0;
                        double ratioF = 0.0d;
                        double ratioE = 0.0d;
                        double ratioFP = 0.0d;
                        double ratioFP2 = 0.0d;
                        double ratioFP3 = 0.0d;
                        double totF[2],totE[2],np[2],totFP[2],totFP2[2],np2[2],totFP3[2],np3[2];
                        double totE2[2],totE4[2],totE8[2],totE16[2];
                        int hx2[2],hx4[2],hx8[2],hx16[2];
                        int nra2[2],nra4[2],nra8[2],nra16[2];
                        double fc2,fc4,fc8,fc16;
                        int avap2,avap4,avap8,avap16;
                        for(int j=0;j<2;j++)
                        {
                            totF[j]=0.0d; totE[j]=0.0d; np[j]=0.0d; totFP[j]=0.0d;
                            np2[j]=0.0d; totFP2[j]=0.0d;  np3[j]=0.0d; totFP3[j]=0.0d;
                            totE2[j]=0.0d; totE4[j]=0.0d; totE8[j]=0.0d; totE16[j]=0.0d;
                            hx2[j]=0; hx4[j]=0; hx8[j]=0; hx16[j]=0;
                            nra2[j]=0; nra4[j]=0; nra8[j]=0; nra16[j]=0;
                        }
                        int avap;
                        int xdep,xfin;
                        xdep=xmin; xfin = xmax;
                        bool ccsra = true;
                        bool alreadyRestarted = false;
                        for(int x=xmin;x<=xmax;x++)
                        {
                            double e = SonogramArray[pc[x-xmin]][x];
                            if(pc[x-xmin]>0)
                            {
                                bool condRestart = false;
                                if(alreadyRestarted==false)
                                {
                                    if((x-xmin)<(xmaitr-xmin)/4)
                                    {
                                        if(e<(_detectionThreshold - _stopThreshold)) condRestart = true;
                                        if(x>xmin+2 && pc[x-xmin-1]==0 && pc[x-xmin-2]==0) condRestart = true;
                                    }
                                }
                                if(condRestart)
                                {
                                    alreadyRestarted = true;
                                    for(int j=0;j<2;j++)
                                    {
                                        totF[j]=0.0d; totE[j]=0.0d; np[j]=0.0d; totFP[j]=0.0d;
                                        np2[j]=0.0d; totFP2[j]=0.0d;  np3[j]=0.0d; totFP3[j]=0.0d;
                                        totF[j]=0.0d; totE[j]=0.0d; np[j]=0.0d; totFP[j]=0.0d;
                                        np2[j]=0.0d; totFP2[j]=0.0d;  np3[j]=0.0d; totFP3[j]=0.0d;
                                        totE2[j]=0.0d; totE4[j]=0.0d; totE8[j]=0.0d; totE16[j]=0.0d;
                                        hx2[j]=0; hx4[j]=0; hx8[j]=0; hx16[j]=0;
                                        nra2[j]=0; nra4[j]=0; nra8[j]=0; nra16[j]=0;
                                    }
                                    xdep = x;
                                }
                                //
                                if(x<=xmaitr) avap=0; else avap=1;
                                totF[avap]+=(double)pc[x-xmin];
                                totE[avap]+=e;
                                totFP[avap]+= _yEbarPerX[x-xmin];
                                totFP2[avap]+= _yEbarPerX[x-xmin]*_tabX[icri][x-xmin];
                                np[avap]++;
                                np2[avap]+=_tabX[icri][x-xmin];
                                totFP3[avap]+= _yEbarPerX[x-xmin]*e;
                                np3[avap]+=e;
                                //
                                if(x>xmaitr) if( ((x-xmaitr)*4) >(xmaitr-xdep)) ccsra = false;
                                //if(x>=xmaitr)  ccsra = false;
                                if(ccsra)
                                {
                                    if(((x-xdep)*16)<=(xmaitr-xdep)) avap16=0; else avap16=1;
                                    if(((x-xdep)*8)<=(xmaitr-xdep)) avap8=0; else avap8=1;
                                    if(((x-xdep)*4)<=(xmaitr-xdep)) avap4=0; else avap4=1;
                                    if(((x-xdep)*2)<=(xmaitr-xdep)) avap2=0; else avap2=1;
                                    fc2=1.0d; fc4=1.0d; fc8=1.0d; fc16=1.0d;
                                    if(avap16==0)
                                    {
                                        if(nra16[0]<3) fc16=4.0d-nra16[0];
                                    }
                                    if(avap8==0)
                                    {
                                        if(nra8[0]<4) fc8=5.0d-nra8[0];
                                    }
                                    if(avap4==0)
                                    {
                                        if(nra4[0]<6) fc4=7.0d-nra4[0];
                                    }
                                    if(avap2==0)
                                    {
                                        if(nra2[0]<9) fc2=10.0d-nra2[0];
                                    }

                                    nra2[avap2]+=fc2; nra4[avap4]+=fc4; nra8[avap8]+=fc8; nra16[avap16]+=fc16;
                                    double ecol = _tabX[icri][x-xmin];
                                    totE2[avap2]+=ecol*fc2; totE4[avap4]+=ecol*fc4; totE8[avap8]+=ecol*fc8; totE16[avap16]+=ecol*fc16;
                                    double hpix = _numberPixelsPerX[x-xmin];
                                    hx2[avap2]+=hpix*fc2; hx4[avap4]+=hpix*fc4; hx8[avap8]+=hpix*fc8; hx16[avap16]+=hpix*fc16;
                                }
                                //
                                if(x<xmax-2 && pc[x-xmin+1]==0 && pc[x-xmin+2]==0)
                                {
                                    if(np[1]>np[0]*3)
                                    {
                                        xfin = x;
                                        break;
                                    }
                                }
                            }
                        }
                        if(np[0]>0 && np[1]>0)
                        {
                            ratioF = (totF[0]*100.0d/np[0]) / (totF[1]/np[1]);
                            ratioE = (totE[0]*100.0d/np[0]) / (totE[1]/np[1]);
                            ratioFP = (totFP[0]*100.0d/np[0]) / (totFP[1]/np[1]);
                            ratioFP2 = (totFP2[0]*10.0d/np2[0]) / (totFP2[1]/np2[1]);
                            ratioFP3 = (totFP3[0]*100.0d/np3[0]) / (totFP3[1]/np3[1]);
                            oParamCrete[jcrete][RAF]  = ratioF;
                            oParamCrete[jcrete][RAE]  = ratioE;
                            oParamCrete[jcrete][RAFE]  = ratioF+ratioE;
                            oParamCrete[jcrete][RAFP]  = ratioFP*KhzPerY;
                            oParamCrete[jcrete][RAFP2]  = ratioFP2*KhzPerY;
                            oParamCrete[jcrete][RAFP3]  = ratioFP3*KhzPerY;
                            if(xmaitr!=xdep)
                            oParamCrete[jcrete][SBMP]  = ( ((double)  (pc[xmaitr-xmin] - pc[xdep-xmin] ))*KhzPerY)  / (((double)(xmaitr-xdep))*MsPerX);
                            if(xmaitr!=xfin)
                            oParamCrete[jcrete][SAMP]  = ( ((double)  (pc[xfin-xmin] - pc[xmaitr-xmin]))*KhzPerY)   / (((double)(xfin-xmaitr))*MsPerX);
                            oParamCrete[jcrete][SBAR]  = oParamCrete[jcrete][SAMP]  - oParamCrete[jcrete][SBMP];
                            if(nra2[0]>0 && nra2[1]>0 )
                            {
                                if(hx2[1]>0)
                                    oParamCrete[jcrete][RAHP2]  = ((((double)hx2[0]*10.0d)/((double)nra2[0]))  / (((double)hx2[1])/((double)nra2[1]))) * KhzPerY;
                                if(totE2[1] > epsilon)
                                    oParamCrete[jcrete][RAHE2]  = (((totE2[0]*10.0d)/((double)nra2[0])) / ((totE2[1])/((double)nra2[1]))) * KhzPerY;
                            }
                            if(nra4[0]>0 && nra4[1]>0)
                            {
                                if(hx4[1]>0)
                                    oParamCrete[jcrete][RAHP4]     = ((((double)hx4[0]*10.0d)/((double)nra4[0])) / (((double)hx4[1])/((double)nra4[1])))*KhzPerY;
                                if(totE4[1] > epsilon)
                                    oParamCrete[jcrete][RAHE4]  = (((totE4[0]*10.0d)/((double)nra4[0])) / ((totE4[1])/((double)nra4[1])))*KhzPerY;

                            }
                            if(nra8[0]>0 && nra8[1]>0 && hx8[1]>0)
                            {
                                if(hx8[1]>0)
                                    oParamCrete[jcrete][RAHP8]  = ((((double)hx8[0]*10.0d)/((double)nra8[0])) / (((double)hx8[1])/((double)nra8[1])))*KhzPerY;
                                if(totE8[1] > epsilon)
                                    oParamCrete[jcrete][RAHE8]  = (((totE8[0]*10.0d)/((double)nra8[0])) / (totE8[1]/((double)nra8[1])))*KhzPerY;
                            }
                            if(nra16[0]>0 && nra16[1]>0 && hx16[1]>0)
                            {
                                if(hx16[1]>0)
                                    oParamCrete[jcrete][RAHP16]  = ((((double)hx16[0]*10.0d)/((double)nra16[0])) / (((double)hx16[1])/((double)nra16[1])))*KhzPerY;
                                if(totE16[1] >0)
                                    oParamCrete[jcrete][RAHE16]  = (((totE16[0]*10.0d)/((double)nra16[0])) / (totE16[1]/((double)nra16[1])))*KhzPerY;
                            }
                        }
                    }
                }
                // -----------------------------------------------------------------------------------------
                int milieu = (xmin+xmax)/2;
                for(int k=milieu;k<=xmax;k++)
                    if(pc[k-xmin]>0)
                    {
                        milieu=k;
                        break;
                    }
                oParamCrete[jcrete][CeF] = (double)pc[milieu-xmin]*KhzPerY;
                if(jcrete==0)
                {
                    int pos1=xmcmax[jcrete];
                    int pos2=xmax;
                    int borne = ((int)emcmax - 5)*100;
                    if(xmin<xmcmax[jcrete])
                    {
                        for(int k=xmin;k<xmcmax[jcrete];k++)
                        {
                            if(pc[k-xmin]>0)
                            {
                                if(SonogramArray[pc[k-xmin]][k]>borne+epsilon)
                                {
                                    pos1=k;
                                    break;
                                }
                            }
                        }
                    }
                    if(xmcmax[jcrete]<xmax)
                    {
                        for(int k=xmcmax[jcrete];k<xmax;k++)
                        {
                            if(pc[k-xmin]>0)
                            {
                                if(SonogramArray[pc[k-xmin]][k]<(borne-epsilon))
                                {
                                    pos2=k;
                                    break;
                                }
                            }
                        }
                    }
                    oParamCrete[jcrete][B5dBBF]  = pc[pos1-xmin]*KhzPerY;
                    oParamCrete[jcrete][B5dBAF]  = pc[pos2-xmin]*KhzPerY;
                    oParamCrete[jcrete][B5dBBW]  = oParamCrete[jcrete][B5dBAF]-oParamCrete[jcrete][B5dBBF];
                    oParamCrete[jcrete][B5dBDur] = (pos2-pos1)*MsPerX;
                }
            }
            else
            {
                meilleur = ymin; meilleurepente = 100000.0d;
                int dx;
                for(int y=ymin;y<=ymax;y++)
                {
                    if(_imageData)
                    {
                        QPoint p;
                        p.setX(pc[y-ymin]);
                        p.setY(y);
                        if(jcrete==3) _callWestRidge.push_back(p);
                        if(jcrete==4) _callSecondWestRidge.push_back(p);
                    }
                    _slope[y-ymin]=9999.0;
                    if(y>ymin+1 && y < ymax-1)
                    {
                        if(pc[y-ymin+2]>0 && pc[y-ymin-2]>0)
                        {

                            dx = qAbs(pc[y-ymin+2]>0 - pc[y-ymin-2]);
                            if(dx>0)
                                _slope[y-ymin] = 4.0d/(double)dx;
                        }
                        if(_slope[y-ymin] < (meilleurepente-epsilon))
                        {
                            meilleur = y;
                            meilleurepente = _slope[y-ymin];
                        }
                    }
                    //
                    double e = (double)SonogramArray[y][pc[y-ymin]]/100.0d;
                    if(e>(emcmax+epsilon))
                    {
                        emcmax = e;
                        xmcmax[jcrete]=pc[y-ymin];
                        ymcmax[jcrete] = y;
                    }

                }
                if(pc[ymax-ymin]!=pc[0])
                {
                    oParamCrete[jcrete][Slope] = ((double)(ymax-ymin)*KhzPerY)/((double)(pc[ymax-ymin]-pc[0])*MsPerX);
                }
                else oParamCrete[jcrete][Slope] = 9999.0f;
                if(ymin != ymax)
                {
                    oParamCrete[jcrete][ISlope] = ((double)(pc[ymax-ymin]-pc[0])*MsPerX) /  ((double)(ymax-ymin)*KhzPerY);
                }
                else oParamCrete[jcrete][ISlope] = 9999.0f;
                lowfc = meilleur;
                LowSlope[(icri*NCRETES+jcrete)*2]=pc[lowfc-ymin];
                LowSlope[(icri*NCRETES+jcrete)*2+1]=lowfc;
                oParamCrete[jcrete][FIF] = lowfc*KhzPerY;
                double pe1,pe2,dif;
                double meilleuredist;
                int meilleurk[2];
                int debb,finb;
                meilleurk[0]=ymin;
                bool adeb = _wavFile.contains("300711");
                meilleurk[1]=lowfc;
                int sens=1;
                if(pc[ymax-ymin]>=pc[0]) sens=1;
                else sens=-1;
                for(int jp=0;jp<2;jp++)
                {
                    meilleuredist = 1000000.0d;
                    if(jp==0)
                    {
                        if(sens==1) {debb=ymin;finb=lowfc;}
                        else {debb=lowfc;finb=ymax;}
                    }
                    else
                    {
                        if(sens==1) {debb=lowfc;finb=ymax;}
                        else {debb=ymin;finb=lowfc;}
                    }
                    if(debb < finb-1)
                    {
                        for (int k = debb+1 ; k<= finb-1 ; k++)
                        {
                            dif=0.0d;
                            if(pc[k-ymin]>0)
                            {
                                if(pc[k-ymin]-pc[debb-ymin] != 0)
                                    pe1=(double)(k-debb)/((double)(pc[k-ymin]-pc[debb-ymin]));
                                else pe1=100000;
                                if(pc[finb-ymin]-pc[k-ymin] != 0)
                                    pe2=((double)(finb-k))/((double)(pc[finb-ymin]-pc[k-ymin]));
                                else pe2=100000;
                                for(int s = debb+1 ; s<k ; s++)
                                    if(pc[s-ymin]>0)
                                    {
                                        if(pe1>99999)
                                        {
                                            dif +=(double)pow((double)(pc[s-ymin]-pc[debb-ymin]),2);

                                        }
                                        else
                                        {
                                            dif += (double)pow((double)(pc[s-ymin]-pc[debb-ymin])-(double)(s-debb)/pe1,2);
                                        }
                                    }

                                for(int s = k+1 ; s<=finb-1 ; s++)
                                    if(pc[s-ymin]>0)
                                    {
                                        if(pe2>99999)
                                        {
                                            dif +=(double)pow((double)(pc[s-ymin]-pc[finb-ymin]),2);
                                        }
                                        else
                                        {
                                            dif += (double)pow((double)(pc[s-ymin]-pc[k-ymin])-(double)(s-k)/pe2,2);
                                        }
                                    }

                                if((double)epsilon<(meilleuredist-dif))
                                {
                                    meilleurk[jp]=k;
                                    meilleuredist=dif;
                                }
                            }

                        }

                    }
                }
                if(meilleurk[0]==0) hifc=ymin; else hifc=meilleurk[0];
                Inflexion1[(icri*NCRETES+jcrete)*2]=pc[hifc-ymin];
                Inflexion1[(icri*NCRETES+jcrete)*2+1]=hifc;
                oParamCrete[jcrete][HCF] = hifc*KhzPerY;

                oParamCrete[jcrete][THCF]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][THCF] = ((double)(pc[meilleurk[0]-ymin]-xmin))/((double)(xmax-xmin));
                fc3=meilleurk[1];
                Inflexion3[(icri*NCRETES+jcrete)*2]=pc[fc3-ymin];
                Inflexion3[(icri*NCRETES+jcrete)*2+1]=fc3;
                oParamCrete[jcrete][LCF] = fc3*KhzPerY;


                oParamCrete[jcrete][UpSl] = 0.0f;
                if((pc[meilleurk[0]-ymin]-pc[0]!=0 && sens==1)
                        || (pc[meilleurk[0]-ymin]-pc[ymax-ymin]!=0 && sens==-1))
                {
                    if(sens==1)
                        oParamCrete[jcrete][UpSl] = ((double)(meilleurk[0]-ymin)*KhzPerY)
                                /( ( (double)(pc[meilleurk[0]-ymin]-pc[0]) ) * MsPerX);
                    else
                        oParamCrete[jcrete][UpSl] = ((double)(meilleurk[0]-ymax)*KhzPerY)
                                /( ( (double)(pc[meilleurk[0]-ymin]-pc[ymax-ymin]) ) * MsPerX);

                }


                oParamCrete[jcrete][LoSl] = 9999.0f;
                if(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]!=0)
                    oParamCrete[jcrete][LoSl] = ((double)(lowfc-meilleurk[0])*KhzPerY)
                            /( ( (double)(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]) ) * MsPerX);

                int debp,finp,ecart;
                double pj[4];
                for(int j2=0;j2<4;j2++)
                {
                    pj[j2]=9999.0d;
                    if(j2<2)
                    {
                        ecart=(ymax-ymin+1)/20;
                        if(j2==0) {debp=0; finp=ecart;}
                        else  {debp=ymax-ymin-ecart; finp=ymax-ymin;}
                    }
                    else
                    {
                        if(j2==2)
                        {
                            debp=ymcmax[jcrete]-ymin-2;
                            finp=ymcmax[jcrete]-ymin+2;
                        }
                        else
                        {
                            debp = lowfc-ymin-2;
                            finp = lowfc-ymin+2;
                        }
                        if(debp<0) debp =0;
                        if(finp>ymax-ymin) finp = ymax-ymin;
                        ecart=finp-debp;
                    }
                    double dx=(double)(pc[finp]-pc[debp]);
                    if(dx!=0)
                    {
                        pj[j2]=((double)ecart*KhzPerY)/(dx*MsPerX);
                    }
                }
                oParamCrete[jcrete][StSl] = pj[0];
                oParamCrete[jcrete][EnSl]   = pj[1];
                oParamCrete[jcrete][FISl]  = pj[3];
                oParamCrete[jcrete][FPSl]  = pj[2];

                if(_paramVersion>=1)
                {
                    oParamCrete[jcrete][SDC]  = 0;
                    oParamCrete[jcrete][SDCR]  = 0;
                    int derdif1=0;
                    int dif1,dif2;
                    int totdif2 = 0;
                    for(int y=ymin+1;y<=ymax;y++)
                    {
                        dif1 = pc[y-ymin]-pc[y-ymin-1];
                        dif2 = qAbs(dif1-derdif1);
                        totdif2 += dif2;
                        derdif1 = dif1;
                    }
                    oParamCrete[jcrete][SDC]  = totdif2*MsPerX;
                    double labs = qAbs(pc[ymax-ymin]-pc[0])+1;
                    oParamCrete[jcrete][SDCR]  = ((double)totdif2/labs)*(MsPerX/KhzPerY);
                }
                if(jcrete==4)
                {
                    int pos1=ymcmax[jcrete];
                    int pos2=ymax;
                    int borne = ((int)emcmax - 5)*100;
                    if(ymin<ymcmax[jcrete])
                    {
                        for(int k=ymin;k<ymcmax[jcrete];k++)
                        {
                            if(pc[k-ymin]>0)
                            {
                                if(SonogramArray[k][pc[k-ymin]]>borne+epsilon)
                                {
                                    pos1=k;
                                    break;
                                }
                            }
                        }
                    }
                    if(ymcmax[jcrete]<ymax)
                    {
                        for(int k=ymcmax[jcrete];k<ymax;k++)
                        {
                            if(pc[k-ymin]>0)
                            {
                                if(SonogramArray[k][pc[k-ymin]]<(borne - epsilon))
                                {
                                    pos2=k;
                                    break;
                                }
                            }
                        }
                    }
                    if(pc[pos2-ymin]<pc[pos1-ymin])
                    {
                        int cpos2 = pos2;
                        pos2 = pos1;
                        pos1 = cpos2;
                    }
                    oParamCrete[jcrete][B5dBBF]  = pos1*KhzPerY;
                    oParamCrete[jcrete][B5dBAF]  = pos2*KhzPerY;
                    oParamCrete[jcrete][B5dBBW]  = oParamCrete[jcrete][B5dBAF]-oParamCrete[jcrete][B5dBBF];
                    oParamCrete[jcrete][B5dBDur] = (pc[pos2-ymin]-pc[pos1-ymin])*MsPerX;
                }
            }
            oParamCrete[jcrete][FPk]=(double)ymcmax[jcrete]*(double)KhzPerY;
            oParamCrete[jcrete][FPkD] = 0.0f;
            if(criprec<icri) oParamCrete[jcrete][FPkD] = oParamCrete[jcrete][FPk] - _paramsArray[criprec][jcrete+1][FPk];
            int intervTemps;
            if(jcrete<3) intervTemps = xmax-xmin+1;
            else intervTemps = tempsOuest[jcrete-3];
            oParamCrete[jcrete][TPk] = (double)(xmcmax[jcrete]-xmin+0.5f)/(double)intervTemps;
        } // end of crests loop 
        oParamCrete[0][StF] = (double)_yEmaxPerX[icri][0]*KhzPerY;
        oParamCrete[0][EnF] = (double)_yEmaxPerX[icri][xmax-xmin]*KhzPerY;
        if(_imageData)
        {
            CallMasterRidgeArray.push_back(_callMasterRidge);
            CallSouthArray.push_back(_callSouthRidge);
            CallNorthRidgeArray.push_back(_callNorthRidge);
            CallWestRidgeArray.push_back(_callWestRidge);
            CallSecondWestRidgeArray.push_back(_callSecondWestRidge);
        }
        oParam[Hup_RFMP] = 0.0f;
        oParam[Hup_PosMP] = 9999.0f;
        oParam[Hup_PosSt] = 9999.0f;
        oParam[Hup_PosEn] = 9999.0f;
        oParam[Hup_AmpDif] = 0.0f;
        oParam[Hup_RSlope] = -9999.0f;
        oParam[Hlo_RFMP] = 0.0f;
        oParam[Hlo_PosMP] = 9999.0f;
        oParam[Hlo_PosSt] = 9999.0f;
        oParam[Hlo_PosEn] = 9999.0f;
        oParam[Hlo_AmpDif] = 0.0f;
        oParam[Hlo_RSlope] = -9999.0f;
    }
    for(int icri=0;icri<nbcris-1;icri++)
    {
        double freqmp = _paramsArray[icri][SH][FreqMP];
        int pmp = _invMp[icri];
        if(pmp<nbcris-1)
        {
            _paramsArray[icri][SH][NextMP1] = (MasterPoints.at(_sortMp[pmp+1]).x()-MasterPoints.at(icri).x())*MsPerX;
            for(int k=pmp+1;k<nbcris;k++)
            {
                if(_paramsArray[_sortMp[k]][SH][Fmin] < (freqmp-epsilon)
                        && _paramsArray[_sortMp[k]][SH][Fmax] > (freqmp+epsilon))
                {
                    _paramsArray[icri][SH][NextMP2] = (MasterPoints.at(_sortMp[k]).x() -MasterPoints.at(icri).x())*MsPerX;
                    break;
                }
            }
        }
        if(pmp>0)
        {

            _paramsArray[icri][SH][PrevMP1]=(MasterPoints.at(icri).x()-MasterPoints.at(_sortMp[pmp-1]).x())*MsPerX;
            for(int k=pmp-1;k>=0;k--)
            {
                if(_paramsArray[_sortMp[k]][SH][Fmin] < (freqmp-epsilon)
                        && _paramsArray[_sortMp[k]][SH][Fmax] > (freqmp+epsilon))
                {
                    _paramsArray[icri][SH][PrevMP2] = (MasterPoints.at(icri).x()-MasterPoints.at(_sortMp[k]).x())*MsPerX;
                    break;
                }
            }
        }
    }
    // --------------------------------------------------------------------
    for(int k=0;k<2;k++)
        for(int i=0;i<nbcris;i++)
        {_harmonic[k][i] = -1; _dpm[k][i]=0.0d; _dypm[k][i]=0; _ypm[k][i]=0;}

    for(int i=0;i<nbcris;i++)
    {
        int xmin1=_vectorXMin.at(i);
        int xmax1=_vectorXMax.at(i);
        int xmaitr1=MasterPoints.at(i).x();
        int ymaitr1=MasterPoints.at(i).y();
        for(int k=0;k<2;k++)
        {
            for(int j=0;j<nbcris;j++)
            {
                if(j==i) continue;
                int xmax2=_vectorXMax.at(j);
                if(xmax2<xmin1) continue;
                int ymaitr2=MasterPoints.at(j).y();
                if((k==0 && ymaitr2<ymaitr1) || (k==1 && ymaitr2>ymaitr1)) continue;
                int xmin2=_vectorXMin.at(j);
                if(xmaitr1>=xmin2 && xmaitr1<=xmax2)
                {
                    int xmaitr2=MasterPoints.at(j).x();
                    if(xmaitr2>=xmin1 && xmaitr2<=xmax1)
                    {
                        double dm =pow(xmaitr2-xmaitr1,2)+pow(ymaitr2-ymaitr1,2);
                        if(_harmonic[k][i]<0 || dm<_dpm[k][i]
                                || (dm==_dpm[k][i] && abs(ymaitr2-ymaitr1)<<_dypm[k][i])
                                || (dm==_dpm[k][i] && abs(ymaitr2-ymaitr1)<<_dypm[k][i]
                                    && ymaitr2<_ypm[k][i])
                                )
                        {
                            _harmonic[k][i]=j;
                            _dpm[k][i]=dm;
                            _dypm[k][i]=abs(ymaitr2-ymaitr1);
                            _ypm[k][i]=ymaitr2;
                        }
                    }
                }
                if(xmin2>xmax1) break;
            }
        }
    }
    //if(_detec->IDebug) _detec->_logText << "-10" << endl;
    int nhsup,nhinf;
    for(int icri=0;icri<nbcris;icri++)
    {
        nhsup=_harmonic[0][icri];
        if(nhsup>=0)
        {
            if(_paramsArray[icri][SH][FreqMP] > epsilon)
                _paramsArray[icri][SH][Hup_RFMP] = _paramsArray[nhsup][SH][FreqMP]
                        /_paramsArray[icri][SH][FreqMP];
            int xmin1=_vectorXMin.at(icri);
            int xmin2=_vectorXMin.at(nhsup);
            int xmax1=_vectorXMax.at(icri);
            int xmax2=_vectorXMax.at(nhsup);
            double larcri=(double)(xmax1-xmin1);
            if(larcri>epsilon)
            {
                _paramsArray[icri][SH][Hup_PosMP]=(double)(MasterPoints.at(nhsup).x()-xmin1)/larcri;
                _paramsArray[icri][SH][Hup_PosSt]=(double)(xmin2-xmin1)/larcri;
                _paramsArray[icri][SH][Hup_PosEn]=(double)(_vectorXMax.at(nhsup)-xmin1)/larcri;
            }
            int xdeb = qMax(xmin1,xmin2);
            int xfin = qMin(xmax1,xmax2);
            double famp1=0.0d;
            double famp2=0.0d;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=(double)_tabX[icri][i-xmin1];
                famp2+=(double)_tabX[nhsup][i-xmin2];
                double alire = (double)_tabX[icri][i-xmin1];
            }
            _paramsArray[icri][SH][Hup_AmpDif] = (famp2 - famp1)*KhzPerY;


            int nouveauxdeb=xdeb;
            int nouveauxfin=xfin;
            for(int i=xdeb;i<=xfin;i++)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhsup][i-xmin2]<=0) nouveauxdeb++;
                else break;

            }
            for(int i=xfin;i>=nouveauxdeb;i--)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhsup][i-xmin2]<=0) nouveauxfin--;
                else break;
            }
            xdeb=nouveauxdeb; xfin=nouveauxfin;
            if(xfin>xdeb)
            {
                double CM_Slope1=((double)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((double)(xfin-xdeb));
                double CM_Slope2=((double)(_yEmaxPerX[nhsup][xfin-xmin2]
                                 -_yEmaxPerX[nhsup][xdeb-xmin2]))
                        /((double)(xfin-xdeb));
                if(CM_Slope1!=0.0d)
                    _paramsArray[icri][SH][Hup_RSlope] = CM_Slope2 / CM_Slope1;
            }
        }
        nhinf=_harmonic[1][icri];
        if(nhinf>=0)
        {
            if(_paramsArray[icri][SH][FreqMP] > epsilon)
                _paramsArray[icri][SH][Hlo_RFMP] = _paramsArray[nhinf][SH][FreqMP]
                        /_paramsArray[icri][SH][FreqMP];
            int xmin1=_vectorXMin.at(icri);
            int xmax1=_vectorXMax.at(icri);
            int xmin2=_vectorXMin.at(nhinf);
            int xmax2=_vectorXMax.at(nhinf);
            double larcri=(double)(xmax1-xmin1);
            if(larcri>epsilon)
            {
                _paramsArray[icri][SH][Hlo_PosMP]=(double)(MasterPoints.at(nhinf).x()-xmin1)/larcri;
                _paramsArray[icri][SH][Hlo_PosSt]=(double)(xmin2-xmin1)/larcri;
                _paramsArray[icri][SH][Hlo_PosEn]=(double)(xmax2-xmin1)/larcri;
            }
            int xdeb = qMax(xmin1,xmin2);
            int xfin = qMin(xmax1,xmax2);
            double famp1=0.0d;
            double famp2=0.0d;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=_tabX[icri][i-xmin1];
                famp2+=_tabX[nhinf][i-xmin2];
            }
            _paramsArray[icri][SH][Hlo_AmpDif] = (famp2-famp1)*KhzPerY;
            int nouveauxdeb=xdeb;
            int nouveauxfin=xfin;
            for(int i=xdeb;i<=xfin;i++)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhinf][i-xmin2]<=0) nouveauxdeb++;
                else break;

            }
            for(int i=xfin;i>=nouveauxdeb;i--)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhinf][i-xmin2]<=0) nouveauxfin--;
                else break;
            }
            xdeb=nouveauxdeb; xfin=nouveauxfin;
            if(xfin>xdeb)
            {
                double CM_Slope1=((double)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((double)(xfin-xdeb));
                double CM_Slope2=((double)(_yEmaxPerX[nhinf][xfin-xmin2]
                                 -_yEmaxPerX[nhinf][xdeb-xmin2]))
                        /((double)(xfin-xdeb));
                if(CM_Slope1!=0.0d)
                    _paramsArray[icri][SH][Hlo_RSlope] = CM_Slope2 / CM_Slope1;
            }
        }
    }
    //
    int nbb;
    double tsono = SonogramWidth *  MsPerX;
    double interv[MAXCRI];
    double variation[MAXCRI];
    double proxiFreq = 2.0d;
    int rband[MAXCRI]; //+
    bool alreadyTreated[MAXCRI];
    for(int i=0;i<nbcris;i++) alreadyTreated[i]=false;

    for(int icri=0;icri<nbcris;icri++)
    {
        if(alreadyTreated[icri]==false)
        {
            _paramsArray[icri][SH][MedInt]   = tsono/2;
            _paramsArray[icri][SH][Int25]  = tsono/2;
            _paramsArray[icri][SH][Int75]  = tsono/2;
            _paramsArray[icri][SH][RInt1] = 1.0f;
            _paramsArray[icri][SH][IntDev]   = 0.0f;
            _paramsArray[icri][SH][SmIntDev]  = 0.0f;
            _paramsArray[icri][SH][LgIntDev]  = 0.0f;
            _paramsArray[icri][SH][VarInt]  = 0.0f;
            _paramsArray[icri][SH][VarSmInt] = 0.0f;
            _paramsArray[icri][SH][VarLgInt] = 0.0f;
            nbb=0;
            double freqmp = _paramsArray[icri][SH][FreqMP];
            for(int k=0;k<nbcris;k++)
            {
                double freqmp2 = _paramsArray[k][SH][FreqMP];
                if(qAbs(freqmp2-freqmp) < (proxiFreq-epsilon)) rband[nbb++]=k;
            }
            if(nbb>1)
            {
                int rbandx[MAXCRI];
                for(int j=0;j<nbb;j++)  rbandx[j] = MasterPoints.at(rband[j]).x();
                sortIntArrays(rband,nbb,rbandx);
            }
            double medianDistance = tsono/2;
            double medianLittleDistance = tsono/2;
            double medianBigDistance = tsono/2;
            int nbi = nbb+1;
            int halfnbi = nbi/2;
            int quartnbi =halfnbi/2;
            if(nbi>1)
            {
                interv[0] = MasterPoints.at(rband[0]).x() * MsPerX;
                interv[nbi-1] = (SonogramWidth - MasterPoints.at(rband[nbb-1]).x()) * MsPerX;
                for(int j=1;j<nbi-1;j++)
                {
                    interv[j] = (MasterPoints.at(rband[j]).x()-MasterPoints.at(rband[j-1]).x()) * MsPerX;
                }
                sortDoubleArray(interv,nbi);
                if(halfnbi*2<nbi) medianDistance = interv[halfnbi];
                else medianDistance= (interv[halfnbi-1] + interv[halfnbi])/2;
                if(nbi<4)
                {
                    medianLittleDistance = interv[0];
                    medianBigDistance = interv[nbi-1];
                }
                else
                {
                    if(quartnbi*2<halfnbi)
                    {
                        medianLittleDistance = interv[quartnbi];
                        medianBigDistance = interv[nbi-1-quartnbi];
                    }
                    else
                    {
                        medianLittleDistance = (interv[quartnbi-1] + interv[quartnbi])/2;
                        medianBigDistance = (interv[nbi-1-quartnbi] + interv[nbi-quartnbi])/2;
                    }
                }
                //
            } // fin if nbi>1
            _paramsArray[icri][SH][MedInt]   = medianDistance;
            _paramsArray[icri][SH][Int25]  = medianLittleDistance;
            _paramsArray[icri][SH][Int75]  = medianBigDistance;
            if(medianLittleDistance>epsilon) _paramsArray[icri][SH][RInt1] = medianBigDistance / medianLittleDistance;
            else _paramsArray[icri][SH][RInt1] = 9999;
            //
            double medianDistanceVariation = 0;
            // a) �cart m�dian
            for(int j=0;j<nbi;j++) variation[j] = qAbs(interv[j]-medianDistance);
            double medianLittleDistanceVariation = variation[0];
            double medianBigDistanceVariation = variation[nbi-1];
            if(nbi>1)
            {
                sortDoubleArray(variation,nbi);
                if(halfnbi*2<nbi) medianDistanceVariation = variation[halfnbi];
                else medianDistanceVariation = (variation[halfnbi-1] + variation[halfnbi])/2;
                // b) �carts sur petits intervalles et sur grands intervalles
                for(int j=0;j<nbi;j++)
                {
                    if(j<halfnbi) variation[j] = qAbs(interv[j]-medianLittleDistance);
                    else variation[j] = qAbs(interv[j]-medianBigDistance);
                }
                if(nbi<4)
                {
                    medianLittleDistanceVariation = variation[0];
                    medianBigDistanceVariation = variation[nbi-1];
                }
                else
                {
                    sortDoubleArray(variation,halfnbi);
                    sortDoubleArray(variation+nbi-halfnbi,halfnbi);
                    if(quartnbi*2<halfnbi)
                    {
                        medianLittleDistanceVariation = variation[quartnbi];
                        medianBigDistanceVariation = variation[nbi-1-quartnbi];
                    }
                    else
                    {
                        medianLittleDistanceVariation = (variation[quartnbi-1] + variation[quartnbi])/2;
                        medianBigDistanceVariation = (variation[nbi-1-quartnbi] + variation[nbi-quartnbi])/2;
                    }
                }
            }
            _paramsArray[icri][SH][IntDev]   = medianDistanceVariation;
            _paramsArray[icri][SH][SmIntDev]  = medianLittleDistanceVariation;
            _paramsArray[icri][SH][LgIntDev]  = medianBigDistanceVariation;
            _paramsArray[icri][SH][VarInt]  = medianDistanceVariation/medianDistance;
            if(medianLittleDistance>epsilon) _paramsArray[icri][SH][VarSmInt]  = medianLittleDistanceVariation/medianLittleDistance;
            _paramsArray[icri][SH][VarLgInt]  = medianBigDistanceVariation/medianBigDistance;
            if(medianLittleDistanceVariation==0) _paramsArray[icri][SH][RIntDev1] = 0.0f;
            else _paramsArray[icri][SH][RIntDev1] = medianBigDistanceVariation / medianLittleDistanceVariation;
            int ymaitr = MasterPoints.at(icri).y();
            for(int jb =0;jb<nbb;jb++)
            {
                int jcri = rband[jb];
                if(MasterPoints.at(jcri).y()==ymaitr && jcri != icri)
                {
                    _paramsArray[jcri][SH][MedInt]   = _paramsArray[icri][SH][MedInt] ;
                    _paramsArray[jcri][SH][Int25]  = _paramsArray[icri][SH][Int25];
                    _paramsArray[jcri][SH][Int75]  = _paramsArray[icri][SH][Int75];
                    _paramsArray[jcri][SH][RInt1] = _paramsArray[icri][SH][RInt1];
                    _paramsArray[jcri][SH][IntDev]   = _paramsArray[icri][SH][IntDev];
                    _paramsArray[jcri][SH][SmIntDev]  = _paramsArray[icri][SH][SmIntDev];
                    _paramsArray[jcri][SH][LgIntDev]  = _paramsArray[icri][SH][LgIntDev];
                    _paramsArray[jcri][SH][VarInt]  = _paramsArray[icri][SH][VarInt];
                    _paramsArray[jcri][SH][VarSmInt] = _paramsArray[icri][SH][VarSmInt];
                    _paramsArray[jcri][SH][VarLgInt] = _paramsArray[icri][SH][VarLgInt];
                    _paramsArray[jcri][SH][RIntDev1]  = _paramsArray[icri][SH][RIntDev1];
                    alreadyTreated[jcri] = true;
                }
            }
        }
    }
}

// methods used to sort arrays of numbers for computing of settings
void DetecTreatment::sortDoubleArray(double *pf,int nbf)
{
    double conserv;
    double epsilon = 0.00001d;
    bool ontrie = true;
    while(ontrie)
    {
        ontrie = false;
        for(int j=0;j<nbf-1;j++)
        {
            bool permuter = false;
            if(pf[j]>(pf[j+1]+epsilon)) permuter = true;
            if(permuter)
            {
                conserv = pf[j];
                pf[j] = pf[j+1];
                pf[j+1] = conserv;
                ontrie = true;
            }
        }
    }
}

void DetecTreatment::SortFloatIndArray(float *pf,int nbf,int *pos)
{
    float epsilon = 0.00001f;
    float conserv;
    int cpos;
    bool ontrie = true;
    for(int j=0;j<nbf;j++) pos[j] = j;
    while(ontrie)
    {
        ontrie = false;
        for(int j=0;j<nbf-1;j++)
        {
            bool permuter = false;
            if(pf[j]>(pf[j+1]+epsilon)) permuter = true;
            if(permuter)
            {
                conserv = pf[j];
                pf[j] = pf[j+1];
                pf[j+1] = conserv;
                cpos = pos[j];
                pos[j] = pos[j+1];
                pos[j+1] = cpos;
                ontrie = true;
            }
        }
    }
}

void DetecTreatment::sortIntArrays(int *pf,int nbf,int *pf2)
{
    int conserv;
    bool ontrie = true;
    while(ontrie)
    {
        ontrie = false;
        for(int j=0;j<nbf-1;j++)
        {
            bool permuter = false;
            if(pf2[j]>pf2[j+1]) permuter = true;
            if(permuter)
            {
                conserv = pf[j];
                pf[j] = pf[j+1];
                pf[j+1] = conserv;
                conserv = pf2[j];
                pf2[j] = pf2[j+1];
                pf2[j+1] = conserv;
                ontrie = true;
            }
        }
    }
}

// this method saves the computed settings in the ".ta" file
void DetecTreatment::saveParameters(const QString& wavFile)
{
    double epsilon = 0.0001d;
    QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ _resultSuffix;
    QFile txtFile;
    txtFile.setFileName(txtFilePath);
    if(txtFile.open(QIODevice::WriteOnly | QIODevice::Text)==false)
    {
        _detec->LogStream  << " impossible to open the file " << txtFilePath << " !" << endl;
        return;
    }
    QTextStream fileStream;
    fileStream.setDevice(&txtFile);
    fileStream.setRealNumberNotation(QTextStream::FixedNotation);
    fileStream.setRealNumberPrecision(2);
    fileStream << "Filename"<< '\t' << "CallNum"
               << '\t' << "Version"<< '\t' << "FileDur"<< '\t' << "SampleRate";
    for(int j=0;j<_numberCallParameters;j++)
        if(VectPar[j].FromVersion<=_paramVersion && (VectPar[j].ToVersion< 0 || VectPar[j].ToVersion>=_paramVersion))  
		fileStream << '\t' << VectPar[j].ColumnTitle;
    fileStream << endl;
    float dur = (double)SonogramWidth*MsPerX/1000.0d;
    float sr = (double)_soundFileInfo.samplerate*TimeExpansion;
    for(int i=0;i<_callsNumber;i++)
    {
        float **callArray = _paramsArray[i];
        fileStream << wavFile << '\t' << QString::number(i)
                    << '\t' << QString::number(_paramVersion) << '\t'  <<  dur << '\t' << sr;
        for(int j=0;j<_numberCallParameters;j++)
        {
			if(VectPar[j].FromVersion<=_paramVersion && (VectPar[j].ToVersion< 0 || VectPar[j].ToVersion>=_paramVersion))  
            fileStream << '\t' <<  (callArray[VectPar[j].ArrayNumber][VectPar[j].ParameterNumber]+epsilon);
        }
        fileStream << endl;
    }
    txtFile.close();
}

// this method compresses the ".ta" file
void DetecTreatment::saveCompressedParameters(const QString& wavFile)
{
    QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ _resultSuffix;
    QString compressedParametersPath = _txtPath+"/"+wavFile.left(wavFile.length()-3) + _resultCompressedSuffix;
    if(QFile::exists(compressedParametersPath)) QFile::remove(compressedParametersPath);
    int resu;
    if(LINWIN==1)
    {
        QString program = "7za";
        QStringList  arguments;
        arguments << "a" << "-tgzip" << compressedParametersPath <<  txtFilePath;
        resu = QProcess::execute(program,arguments);
        if(resu==0) QFile::remove(txtFilePath);
    }
    else
    {
        QString program = "gzip";
        QStringList  arguments;
        arguments << txtFilePath;
        resu = QProcess::execute(program,arguments);
    }
}



