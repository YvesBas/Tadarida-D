#include "detec.h"
class Detec;

// 21/4/2015
// retiré les paramètres ajoutés récemment


ParamToSave::ParamToSave(int numTableau,int numPar,QString columnTitle)
{
    NumTableau=numTableau;
    NumPar=numPar;
    ColumnTitle=columnTitle;
}

ParamToSave::ParamToSave()
{
}

ParamToSave::~ParamToSave()
{
}

DetecTreatment::DetecTreatment(Detec *pDet)
{
    _detec = pDet;
    _freqCallMin=8.0f;
    _paramVersion = 1;
    _compressedVersion = 1;
    ResultSuffix = QString("ta");
    ResultCompressedSuffix = QString("tac");
    initVectorParams();
    //InitializeDetecTreatment();
}

DetecTreatment::DetecTreatment(Recherche *r)
{
    initVectorParams();
}

DetecTreatment::~DetecTreatment()
{
}

void DetecTreatment::InitializeDetecTreatment()
{
    //£ _sonogramArray = new float*[FFT_HEIGHT_HALF_MAX];
    //µ _sonogramArray = new qint16*[FFT_HEIGHT_HALF_MAX];
    //µ _pointFlagsArray = new char *[FFT_HEIGHT_HALF_MAX];
    _sonogramArray = new qint16*[MAXHEIGHT];
    _pointFlagsArray = new char *[MAXHEIGHT];
    //£ float *smy;
	//§ :
    //int ld8 = (SONOGRAM_WIDTH_MAX+7)/8;
    //µ for ( int i=0 ; i < FFT_HEIGHT_HALF_MAX ; i++)
    for ( int i=0 ; i <MAXHEIGHT ; i++)
    {
        //£ _sonogramArray[i] = new float[SONOGRAM_WIDTH_MAX];
        _sonogramArray[i] = new qint16[SONOGRAM_WIDTH_MAX];
        //§ _pointFlagsArray[i] = new char[SONOGRAM_WIDTH_MAX];
        _pointFlagsArray[i] = new char[LD8];
        //£ for(int j = 0; j < SONOGRAM_WIDTH_MAX; j++) *smy++ = 0.0f;
        //for(int j = 0; j < SONOGRAM_WIDTH_MAX; j++) *smy++ = 0;
        memset(_sonogramArray[i],0,SONOGRAM_WIDTH_MAX*2);
    }
    _energyMin = 0.0f;
    _charParamsArray = new char[MAXCRI*NBTAB*NBPAR*sizeof(float)];
    _paramsArray = new float**[MAXCRI];
    //_simpleParamsArray = new float*[MAXCRI];
    //_valuesToCompressArray = new float*[MAXCRI];
    _numberCallParameters = _vectPar.size();

    for (int i=0;i < MAXCRI; i++)
    {
        _paramsArray[i] = new float *[NBTAB];
       // _simpleParamsArray[i] = new float[_numberCallParameters];
       // _valuesToCompressArray[i] = new float[_numberCallParameters];
        for(int j=0;j < NBTAB; j++)
            _paramsArray[i][j] = (float *)((char *)(_charParamsArray+(i*NBTAB*NBPAR+j*NBPAR)*sizeof(float)));
    }
    _lowSlope = new int[NCRETES*MAXCRI*2];
    _inflexion1 = new int[NCRETES*MAXCRI*2];
    _inflexion3 = new int[NCRETES*MAXCRI*2];
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
    //£ _xMinPerY = new int[MAXHEIGHT];
    _xMinPerY = new quint16[MAXHEIGHT];
    //£ _xSecondWestRidgePerY = new int[MAXHEIGHT];
    _xSecondWestRidgePerY = new quint16[MAXHEIGHT];
    _xMaxPerY = new int[MAXHEIGHT];
    //£ _yMinPerX = new int[MAXLARCRI];
    _yMinPerX = new quint16[MAXLARCRI];
    //£ _yMaxPerX = new int[MAXLARCRI];
    _yMaxPerX = new quint16[MAXLARCRI];
    _yEbarPerX = new float[MAXLARCRI];
    _totYEPerX = new float[MAXLARCRI];
    _eMaxPerX = new float[MAXLARCRI];
    _slope    = new float[qMax(MAXLARCRI,MAXHEIGHT)];
    _charTabX = new char[MAXCRI*MAXLARCRI*sizeof(float)];
    _tabX = new float*[MAXCRI];
    for(int i=0 ; i < MAXCRI ; i++) _tabX[i] = (float *)((char *)(_charTabX+(i*MAXLARCRI*sizeof(float))));
    //
    _charYEmaxPerX = new char[MAXCRI*MAXLARCRI*sizeof(int)];
    //£ _yEmaxPerX = new int*[MAXCRI];
    _yEmaxPerX = new quint16*[MAXCRI];
    //£ for(int i=0;i<MAXCRI;i++)  _yEmaxPerX[i] = (int *)((char *)(_charYEmaxPerX+(i*MAXLARCRI*sizeof(int))));
    for(int i=0;i<MAXCRI;i++)  _yEmaxPerX[i] = (quint16 *)((char *)(_charYEmaxPerX+(i*MAXLARCRI*sizeof(int))));
    _charTabYX = new char[MAXHEIGHT*MAXLARCRI*sizeof(float)];
    _tabYX = new float*[MAXHEIGHT];
    for(int i=0;i<MAXHEIGHT;i++) _tabYX[i] = (float *)((char *)(_charTabYX+(i*MAXLARCRI*sizeof(float))));
    //
    _energyMoyCol = new float[SONOGRAM_WIDTH_MAX];
    _flagGoodCol = new bool[SONOGRAM_WIDTH_MAX];
    _tabr1 = new int[MAXCRI];
}

void DetecTreatment::SetDirParameters(QString wavPath,QString txtPath,bool imageData,
                                                                       QString imagePath,QString datPath)
{
    _wavPath = wavPath;
    _txtPath = txtPath;
    _imageData = imageData;
    _imagePath = imagePath;
    _datPath = datPath;
}

void DetecTreatment::SetGlobalParameters(int timeExpansion,int seuilDetect,int seuilStop,
                                 int freqMin,int nbo,bool useValflag,
                                int jumpThreshold,int widthBigControl,int widthLittleControl,
                                int highThreshold,int lowThreshold,int qR,int qN,bool withNewParams)
{
    _timeExpansion = timeExpansion;
    _detectionThreshold = seuilDetect;
    _stopThreshold = seuilStop;
    _freqMin = freqMin;
    _nbo = nbo;
    _useValflag = useValflag;
    _jumpThreshold = jumpThreshold;
    _widthBigControl = widthBigControl;
    _widthLittleControl = widthLittleControl;
    _highThreshold = highThreshold;
    _lowThreshold = lowThreshold;
    _qR = qR;
    _qN = qN;
    _withNewParams = withNewParams;
}

void DetecTreatment::initVectorParams()
{

    QString prefix[] = {"","CM_","CS_","CN_","CO_","CO2_"};
    _vectPar.push_back(ParamToSave(SH,StTime,"StTime"));
    _vectPar.push_back(ParamToSave(SH,Dur,"Dur"));
    _vectPar.push_back(ParamToSave(SH,PrevSt,"PrevSt"));
    _vectPar.push_back(ParamToSave(SH,Fmax,"Fmax"));
    _vectPar.push_back(ParamToSave(SH,Fmin,"Fmin"));
    _vectPar.push_back(ParamToSave(SH,BW,"BW"));
    _vectPar.push_back(ParamToSave(SH,FreqMP,"FreqMP"));
    _vectPar.push_back(ParamToSave(SH,PosMP,"PosMP"));
    _vectPar.push_back(ParamToSave(SH,FreqPkS,"FreqPkS"));
    _vectPar.push_back(ParamToSave(SH,FreqPkM,"FreqPkM"));
    _vectPar.push_back(ParamToSave(SH,PosPkS,"PosPkS"));
    _vectPar.push_back(ParamToSave(SH,PosPkM,"PosPkM"));
    _vectPar.push_back(ParamToSave(SH,FreqPkS2,"FreqPkM1"));
    _vectPar.push_back(ParamToSave(SH,FreqPkM2,"FreqPkM2"));
    _vectPar.push_back(ParamToSave(SH,PrevMP1,"PrevMP1"));
    _vectPar.push_back(ParamToSave(SH,PrevMP2,"PrevMP2"));
    _vectPar.push_back(ParamToSave(SH,NextMP1,"NextMP1"));
    _vectPar.push_back(ParamToSave(SH,NextMP2,"NextMP2"));
    _vectPar.push_back(ParamToSave(SH,Amp1,"Amp1"));
    _vectPar.push_back(ParamToSave(SH,Amp2,"Amp2"));
    _vectPar.push_back(ParamToSave(SH,Amp3,"Amp3"));
    _vectPar.push_back(ParamToSave(SH,Amp4,"Amp4"));
    _vectPar.push_back(ParamToSave(SH,NoisePrev,"NoisePrev"));
    _vectPar.push_back(ParamToSave(SH,NoiseNext,"NoiseNext"));
    _vectPar.push_back(ParamToSave(SH,NoiseDown,"NoiseDown"));
    _vectPar.push_back(ParamToSave(SH,NoiseUp,"NoiseUp"));
    _vectPar.push_back(ParamToSave(SH,CVAmp,"CVAmp"));
    _vectPar.push_back(ParamToSave(CO,Dur,"CO_Dur"));
    _vectPar.push_back(ParamToSave(CO2,Dur,"CO2_Dur"));
    _vectPar.push_back(ParamToSave(CM,Fmax,"CM_Fmax"));
    _vectPar.push_back(ParamToSave(CS,Fmax,"CS_Fmax"));
    //_vectPar.push_back(ParamToSave(CN,Fmax,"CN_Fmax"));
    _vectPar.push_back(ParamToSave(CM,Fmin,"CM_Fmin"));
    //_vectPar.push_back(ParamToSave(CS,Fmin,"CS_Fmin"));
    _vectPar.push_back(ParamToSave(CN,Fmin,"CN_Fmin"));
    _vectPar.push_back(ParamToSave(CM,BW,"CM_BW"));
    _vectPar.push_back(ParamToSave(CS,BW,"CS_BW"));
    _vectPar.push_back(ParamToSave(CN,BW,"CN_BW"));
    // _vectPar.push_back(ParamToSave(CM,FPk,"CM_FPk"));
    _vectPar.push_back(ParamToSave(CO2,FPk,"CO2_FPk"));
    // _vectPar.push_back(ParamToSave(CM,FPkD,"CM_FPkD"));
    _vectPar.push_back(ParamToSave(CO2,FPkD,"CO2_FPkD"));
    //_vectPar.push_back(ParamToSave(CM,TPk,"CM_Ldom"));
    _vectPar.push_back(ParamToSave(CO2,TPk,"CO2_TPk"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Slope,prefix[i]+"Slope"));
    // ajouté le 30/03
    // RAR : for(int i=CO;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,ISlope,prefix[i]+"ISlope"));
    if(_withNewParams) for(int i=CO;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,ISlope,prefix[i]+"ISlope"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,HCF,prefix[i]+"HCF"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,THCF,prefix[i]+"THCF"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,FIF,prefix[i]+"FIF"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,LCF,prefix[i]+"LCF"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,UpSl,prefix[i]+"UpSl"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,LoSl,prefix[i]+"LoSl"));
    _vectPar.push_back(ParamToSave(CM,StF,"CM_StF"));
    _vectPar.push_back(ParamToSave(CM,EnF,"CM_EnF"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,StSl,prefix[i]+"StSl"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,EnSl,prefix[i]+"EnSl"));
    for(int i=CM;i<=CO2;i++)  _vectPar.push_back(ParamToSave(i,FISl,prefix[i]+"FPSl"));
    _vectPar.push_back(ParamToSave(CM,FISl,"CM_FISl"));
    _vectPar.push_back(ParamToSave(CO2,FISl,"CO2_FISl"));
    for(int i=CM;i<=CN;i++) _vectPar.push_back(ParamToSave(i,CeF,prefix[i]+"CeF"));
    _vectPar.push_back(ParamToSave(CM,B5dBBF,"CM_5dBBF"));
    _vectPar.push_back(ParamToSave(CM,B5dBAF,"CM_5dBAF"));
    _vectPar.push_back(ParamToSave(CM,B5dBBW,"CM_5dBBW"));
    _vectPar.push_back(ParamToSave(CM,B5dBDur,"CM_5dBDur"));

    _vectPar.push_back(ParamToSave(CO2,B5dBBF,"CO2_5dBBF"));
    _vectPar.push_back(ParamToSave(CO2,B5dBAF,"CO2_5dBAF"));
    _vectPar.push_back(ParamToSave(CO2,B5dBBW,"CO2_5dBBW"));
    _vectPar.push_back(ParamToSave(CO2,B5dBDur,"CO2_5dBDur"));

    _vectPar.push_back(ParamToSave(SH,Hup_RFMP,"Hup_RFMP"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosMP,"Hup_PosMP"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosSt,"Hup_PosSt"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosEn,"Hup_PosEn"));
    _vectPar.push_back(ParamToSave(SH,Hup_AmpDif,"Hup_AmpDif"));
    _vectPar.push_back(ParamToSave(SH,Hup_RSlope,"Hup_RSlope"));
    _vectPar.push_back(ParamToSave(SH,Hlo_RFMP,"Hlo_RFMP"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosMP,"Hlo_PosMP"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosSt,"Hlo_PosSt"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosEn,"Hlo_PosEn"));
    _vectPar.push_back(ParamToSave(SH,Hlo_AmpDif,"Hlo_AmpDif"));
    _vectPar.push_back(ParamToSave(SH,Hlo_RSlope,"Hlo_RSlope"));

    _vectPar.push_back(ParamToSave(SH,Ramp_2_1,"Ramp_2_1"));
    _vectPar.push_back(ParamToSave(SH,Ramp_3_1,"Ramp_3_1"));
    _vectPar.push_back(ParamToSave(SH,Ramp_3_2,"Ramp_3_2"));
    _vectPar.push_back(ParamToSave(SH,Ramp_1_2,"Ramp_1_2"));
    _vectPar.push_back(ParamToSave(SH,Ramp_4_3,"Ramp_4_3"));
    _vectPar.push_back(ParamToSave(SH,Ramp_2_3,"Ramp_2_3"));
    _vectPar.push_back(ParamToSave(SH,RAN_2_1,"RAN_2_1"));
    _vectPar.push_back(ParamToSave(SH,RAN_3_1,"RAN_3_1"));
    _vectPar.push_back(ParamToSave(SH,RAN_3_2,"RAN_3_2"));
    _vectPar.push_back(ParamToSave(SH,RAN_1_2,"RAN_1_2"));
    _vectPar.push_back(ParamToSave(SH,RAN_4_3,"RAN_4_3"));
    _vectPar.push_back(ParamToSave(SH,RAN_2_3,"RAN_2_3"));

    _vectPar.push_back(ParamToSave(SH,HetX,"HetX"));
    _vectPar.push_back(ParamToSave(SH,HetY,"HetY"));
    _vectPar.push_back(ParamToSave(SH,Dbl8,"Dbl8"));
    _vectPar.push_back(ParamToSave(SH,Stab,"Stab"));

    _vectPar.push_back(ParamToSave(SH,HeiET,"HeiET"));
    _vectPar.push_back(ParamToSave(SH,HeiEM,"HeiEM"));
    _vectPar.push_back(ParamToSave(SH,HeiRT,"HeiRT"));
    _vectPar.push_back(ParamToSave(SH,HeiRM,"HeiRM"));
    _vectPar.push_back(ParamToSave(SH,HeiETT,"HeiETT"));
    _vectPar.push_back(ParamToSave(SH,HeiEMT,"HeiEMT"));
    _vectPar.push_back(ParamToSave(SH,HeiRTT,"HeiRTT"));
    _vectPar.push_back(ParamToSave(SH,HeiRMT,"HeiRMT"));
    _vectPar.push_back(ParamToSave(SH,MedInt,"MedInt"));
    _vectPar.push_back(ParamToSave(SH,Int25,"Int25"));
    _vectPar.push_back(ParamToSave(SH,Int75,"Int75"));
    _vectPar.push_back(ParamToSave(SH,RInt1,"RInt1"));
    _vectPar.push_back(ParamToSave(SH,IntDev,"IntDev"));
    _vectPar.push_back(ParamToSave(SH,SmIntDev,"SmIntDev"));
    _vectPar.push_back(ParamToSave(SH,LgIntDev,"LgIntDev"));
    _vectPar.push_back(ParamToSave(SH,VarInt,"VarInt"));
    _vectPar.push_back(ParamToSave(SH,VarSmInt,"VarSmInt"));
    _vectPar.push_back(ParamToSave(SH,VarLgInt,"VarLgInt"));
    _vectPar.push_back(ParamToSave(SH,RIntDev1,"RIntDev1"));
    _vectPar.push_back(ParamToSave(SH,EnStabSm,"EnStabSm"));
    _vectPar.push_back(ParamToSave(SH,EnStabLg,"EnStabLg"));
    _vectPar.push_back(ParamToSave(SH,HetXr,"HetXr"));
    _vectPar.push_back(ParamToSave(SH,HetYr,"HetYr"));
    _vectPar.push_back(ParamToSave(SH,HetYr2,"HetYr2"));
    _vectPar.push_back(ParamToSave(SH,HetCMC,"HetCMC"));
    _vectPar.push_back(ParamToSave(SH,HetCMD,"HetCMD"));
    _vectPar.push_back(ParamToSave(SH,HetCTC,"HetCTC"));
    _vectPar.push_back(ParamToSave(SH,HetCTD,"HetCTD"));
    _vectPar.push_back(ParamToSave(SH,HetCMnP,"HetCMnP"));
    _vectPar.push_back(ParamToSave(SH,HetCMfP,"HetCMfP"));
    _vectPar.push_back(ParamToSave(SH,HetCTnP,"HetCTnP"));
    _vectPar.push_back(ParamToSave(SH,HetCTfP,"HetCTfP"));

    _vectPar.push_back(ParamToSave(SH,HetPicsMAD,"HetPicsMAD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsMALD,"HetPicsMALD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsMABD,"HetPicsMABD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsMRBLD,"HetPicsMRLBD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsTAD,"HetPicsTAD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsTALD,"HetPicsTALD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsTABD,"HetPicsTABD"));
    _vectPar.push_back(ParamToSave(SH,HetPicsTRBLD,"HetPicsTRLBD"));
    _vectPar.push_back(ParamToSave(SH,VDPicsM,"VDPicsM"));
    _vectPar.push_back(ParamToSave(SH,VLDPicsM,"VLDPicsM"));
    _vectPar.push_back(ParamToSave(SH,VBDPicsM,"VBDPicsM"));
    _vectPar.push_back(ParamToSave(SH,VDPPicsM,"VDPPicsM"));
    _vectPar.push_back(ParamToSave(SH,VLDPPicsM,"VLDPPicsM"));
    _vectPar.push_back(ParamToSave(SH,VBDPPicsM,"VBDPPicsM"));
    _vectPar.push_back(ParamToSave(SH,VDPicsT,"VDPicsT"));
    _vectPar.push_back(ParamToSave(SH,VLDPicsT,"VLDPicsT"));
    _vectPar.push_back(ParamToSave(SH,VBDPicsT,"VBDPicsT"));
    _vectPar.push_back(ParamToSave(SH,VDPPicsT,"VDPPicsT"));
    _vectPar.push_back(ParamToSave(SH,VLDPPicsT,"VLDPPicsT"));
    _vectPar.push_back(ParamToSave(SH,VBDPPicsT,"VBDPPicsT"));
    if(_withNewParams)
    {
        for(int i=CM;i<=CO2;i++)
        {
            _vectPar.push_back(ParamToSave(i,SDC,prefix[i]+"SDC"));
            _vectPar.push_back(ParamToSave(i,SDCR,prefix[i]+"SDCR"));
        }
        _vectPar.push_back(ParamToSave(CM,SDCRY,"CM_SDCRY"));
        _vectPar.push_back(ParamToSave(CS,SDCRY,"CS_SDCRY"));
        _vectPar.push_back(ParamToSave(CM,SDCRXY,"CM_SDCRXY"));
        _vectPar.push_back(ParamToSave(CS,SDCRXY,"CS_SDCRXY"));
        //

        for(int i=CM;i<=CS;i++)
        {
            _vectPar.push_back(ParamToSave(i,SDCL,prefix[i]+"SDCL"));
            _vectPar.push_back(ParamToSave(i,SDCLR,prefix[i]+"SDCLR"));
            _vectPar.push_back(ParamToSave(i,SDCLRY,prefix[i]+"SDCLRY"));
            _vectPar.push_back(ParamToSave(i,SDCLRXY,prefix[i]+"SDCLRXY"));
            _vectPar.push_back(ParamToSave(i,SDCLRXY2,prefix[i]+"SDCLRXY2"));
            //
            _vectPar.push_back(ParamToSave(i,SDCLOP,prefix[i]+"SDCLOP"));
            _vectPar.push_back(ParamToSave(i,SDCLROP,prefix[i]+"SDCLROP"));
            _vectPar.push_back(ParamToSave(i,SDCLRYOP,prefix[i]+"SDCLRYOP"));
            _vectPar.push_back(ParamToSave(i,SDCLRXYOP,prefix[i]+"SDCLRXYOP"));
            //
            _vectPar.push_back(ParamToSave(i,SDCLWB,prefix[i]+"SDCLWB"));
            _vectPar.push_back(ParamToSave(i,SDCLRWB,prefix[i]+"SDCLRWB"));
            _vectPar.push_back(ParamToSave(i,SDCLRYWB,prefix[i]+"SDCLRYWB"));
            _vectPar.push_back(ParamToSave(i,SDCLRXYWB,prefix[i]+"SDCLRXYWB"));
            //
            _vectPar.push_back(ParamToSave(i,SDCLOPWB,prefix[i]+"SDCLOPWB"));
            _vectPar.push_back(ParamToSave(i,SDCLROPWB,prefix[i]+"SDCLROPWB"));
            _vectPar.push_back(ParamToSave(i,SDCLRYOPWB,prefix[i]+"SDCLRYOPWB"));
            _vectPar.push_back(ParamToSave(i,SDCLRXYOPWB,prefix[i]+"SDCLRXYOPWB"));
            //
            _vectPar.push_back(ParamToSave(i,SDCL_DNP,prefix[i]+"SDCL_DNP"));
            _vectPar.push_back(ParamToSave(i,SDCLR_DNP,prefix[i]+"SDCLR_DNP"));
            _vectPar.push_back(ParamToSave(i,SDCLRY_DNP,prefix[i]+"SDCLRY_DNP"));
            _vectPar.push_back(ParamToSave(i,SDCLRXY_DNP,prefix[i]+"SDCLRXY_DNP"));
            _vectPar.push_back(ParamToSave(i,SDCLRXY2_DNP,prefix[i]+"SDCLRXY2_DNP"));
        }
        _vectPar.push_back(ParamToSave(CM,ELBPOS,"CM_ELBPOS"));
        _vectPar.push_back(ParamToSave(CS,ELBPOS,"CS_ELBPOS"));
        _vectPar.push_back(ParamToSave(CM,ELBSB,"CM_ELBSB"));
        _vectPar.push_back(ParamToSave(CS,ELBSB,"CS_ELBSB"));
        //
        _vectPar.push_back(ParamToSave(CM,ELB2POS,"CM_ELB2POS"));
        _vectPar.push_back(ParamToSave(CS,ELB2POS,"CS_ELB2POS"));
        _vectPar.push_back(ParamToSave(CM,ELB2SB,"CM_ELB2SB"));
        _vectPar.push_back(ParamToSave(CS,ELB2SB,"CS_ELB2SB"));
        //
        _vectPar.push_back(ParamToSave(CM,RAF,"CM_RAF"));
        _vectPar.push_back(ParamToSave(CM,RAE,"CM_RAE"));
        _vectPar.push_back(ParamToSave(CM,RAFE,"CM_RAFE"));
        _vectPar.push_back(ParamToSave(CM,RAFP,"CM_RAFP"));
        _vectPar.push_back(ParamToSave(CM,RAFP2,"CM_RAFP2"));
        _vectPar.push_back(ParamToSave(CM,RAFP3,"CM_RAFP3"));
        //
        _vectPar.push_back(ParamToSave(CM,SBMP,"CM_SBMP"));
        _vectPar.push_back(ParamToSave(CM,SAMP,"CM_SAMP"));
        _vectPar.push_back(ParamToSave(CM,SBAR,"CM_SBAR"));
        //
        _vectPar.push_back(ParamToSave(CM,RAHP2,"RAHP2"));
        _vectPar.push_back(ParamToSave(CM,RAHP4,"RAHP4"));
        _vectPar.push_back(ParamToSave(CM,RAHP8,"RAHP8"));
        _vectPar.push_back(ParamToSave(CM,RAHP16,"RAHP16"));
        _vectPar.push_back(ParamToSave(CM,RAHE2,"RAHE2"));
        _vectPar.push_back(ParamToSave(CM,RAHE4,"RAHE4"));
        _vectPar.push_back(ParamToSave(CM,RAHE8,"RAHE8"));
        _vectPar.push_back(ParamToSave(CM,RAHE16,"RAHE16"));
    }

}

void DetecTreatment::EndDetecTreatment()
{
    //µ for ( int i=0 ; i < FFT_HEIGHT_HALF_MAX ; i++)
    for ( int i=0 ; i < MAXHEIGHT ; i++)
    {
        delete[] _sonogramArray[i];
        delete[] _pointFlagsArray[i];
    }
    delete[] _sonogramArray;
    delete[] _pointFlagsArray;
    for (int i=0;i < MAXCRI; i++) delete[] _paramsArray[i];
    delete[] _paramsArray;
    //delete[] _simpleParamsArray;
    //delete[] _valuesToCompressArray;

    delete[] _charParamsArray;
    delete[] _lowSlope;
    delete[] _inflexion1;
    delete[] _inflexion3;
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
    //
    delete[] _energyMoyCol;
    delete[] _flagGoodCol;
    delete[] _tabr1;
}


bool DetecTreatment::CallTreatmentsForOneFile(QString& wavFile,QString &pathFile)
{
    //_detec->_logText << " DetecTreatment  wavFile=" << wavFile << " pathFile=" << pathFile << endl;
    qint64 d0 = QDateTime::currentDateTime().toMSecsSinceEpoch();
   _wavFile = wavFile;
    int d[5];
    clearVars();
    if (openWavFile(pathFile))
    {
        if(computeFFT(pathFile))
        {

            d[0]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->_logText <<   "Après computeFFT  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            correctNoise();
            //_detec->_logText <<   "entre les 2 méthodes noise  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            // if(_detec->_imageData) calculateMedianNoise();
            d[1]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
             _detec->_logText <<   "Après correctNoise...  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            shapesDetects();
            d[2]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->_logText <<   "Après shapesDetects...  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            _callsNumber = (int)_callsArray.size();
            detectsParameter2();
            d[3]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->_logText <<   "Après detectsParameter2...  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            saveParameters(wavFile);
            d[4]=(int)(QDateTime::currentDateTime().toMSecsSinceEpoch()-d0);
            _detec->_logText <<   "Après saveParam... fin.  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            if(_detec->_timeFileOpen)
            {
                _detec->_timeStream  << wavFile << '\t' << d[0] << '\t';
                for(int j=1;j<=4;j++) _detec->_timeStream  << d[j] - d[j-1] << '\t';
                _detec->_timeStream  << d[4] << endl;
            }

            if(_detec->MustCompress) saveCompressedParameters(wavFile);
            // expandParameters(wavFile);
        }
        else
        {
            // erreur dans computefft
            return(false);
        }
    }
    else
    {
        // erreur dans openwavfile
        return(false);
    }
    return(true);
}

void DetecTreatment::clearVars()
{
    _callsArray.clear();
    _vectorXMin.clear();
    _vectorXMax.clear();
    _masterPoints.clear();
    if(_imageData)
    {
        _callMasterRidgeArray.clear();
        _callSouthArray.clear();
        _callNorthRidgeArray.clear();
        _callWestRidgeArray.clear();
        _callSecondWestRidgeArray.clear();
    }
}


bool DetecTreatment::openWavFile(QString& wavFile)
{
    if (! (_soundFile = sf_open(wavFile.toStdString().c_str(), SFM_READ, &_soundFileInfo)))
    {
        if(_detec->_errorFileOpen) _detec->_errorStream << wavFile << ": " << sf_strerror (NULL) << endl;
        return  false;
    }
    if (_soundFileInfo.channels > 1)
    {
        sf_close (_soundFile);
        if(_detec->_errorFileOpen) _detec->_errorStream << wavFile << ": " << "multi-channel non traité" << endl;
        return  false;
    }
    // edit yves - prise en compte tx ech Vigie Chiro
    if (_soundFileInfo.samplerate*_timeExpansion >= 2400000 ) _fftHeight = 4096;

    else {
        if (_soundFileInfo.samplerate*_timeExpansion >= 1200000 ) _fftHeight = 2048;

        else {

            if (_soundFileInfo.samplerate*_timeExpansion >= 600000 ) _fftHeight = 1024;

            else {
                if (_soundFileInfo.samplerate*_timeExpansion >= 300000 ) _fftHeight = 512;
                else {
                    if (_soundFileInfo.samplerate*_timeExpansion >= 150000 ) _fftHeight = 256;
                    else
                        _fftHeight = 128;
                }
            }
        }
    }
    return true;
}

bool DetecTreatment::computeFFT(QString &wavFile)
{
    int iCount;
    int readcount;
    float a = 0.0f;
    _energyMax        = 0.0f;
    _energyMin        = 0.0f;
    _fftHeightHalf		= (int)ceil((float)_fftHeight/(float)2);

    _iOverlapMoving	= (int)ceil((float)_fftHeight/(float)(_nbo*2));
    _sonogramWidth		= (int)ceil(_nbo*2*(float)_soundFileInfo.frames/(float)_fftHeight)+1;
    _data				= ( float* ) fftwf_malloc( sizeof( float ) * _fftHeight );
    _fftRes 		= ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * _fftHeight );
    _complexInput        = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * _fftHeight );
    _msPerX =(float)(_fftHeightHalf*1000)/(_nbo*_soundFileInfo.samplerate*_timeExpansion); //Time: msec
    //_fIStream << _sonogramWidth*_msPerX/1000 << '\t';
    _khzPerY =(float)(_soundFileInfo.samplerate*_timeExpansion)/(float)(_fftHeight*1000); //Freq:khz
    //_detec->_logText << "_sonogramWidth = " << _sonogramWidth << endl;
    if(_sonogramWidth*_msPerX < 10.0f)
    {
        if(_detec->_errorFileOpen) _detec->_errorStream << wavFile << ": fichier trop petit" << endl;
        _detec->_logText << "durée trop petite : " << _sonogramWidth*_msPerX << " ms" << endl;
        return  false;
    }
    if(_sonogramWidth > SONOGRAM_WIDTH_MAX)
    {
        if(_detec->_errorFileOpen) _detec->_errorStream << wavFile << ": fichier trop grand" << endl;
        _detec->_logText << wavFile << ": fichier trop grand" << endl;
        return  false;
    }
    sf_seek(_soundFile, 0, SEEK_END);
    _plan = fftwf_plan_dft_1d( _fftHeight, _complexInput, _fftRes, FFTW_FORWARD, FFTW_ESTIMATE );
    float fact1=2.0f*PI;
    float fact2=4.0f*PI;
    float quot1=_fftHeightHalf-1;
    //µ :
    _limY = qMin(_fftHeightHalf,MAXHEIGHT);
    _coeff = new float[_fftHeightHalf];
    for (int i = 0 ; i < _fftHeightHalf ; i++)
    {
        _coeff[i] = 0.435f - 0.5f*cos(fact1*i/quot1)+ 0.065f*cos(fact2*i/quot1);
    }
    for (int iLoop = 0 ; iLoop < _nbo; iLoop++)
    {
        iCount = 0;
        sf_seek(_soundFile, iLoop * _iOverlapMoving, SEEK_SET);
        while ((readcount = (int)sf_read_float(_soundFile, _data, _fftHeightHalf)))
        {
            for (int i = 0 ; i < _fftHeightHalf ; i++)
            {
                _complexInput[i][0] =_data[i] * _coeff[i];
                _complexInput[i][1] = 0.0f;
            }

            for (int i = _fftHeightHalf ; i < _fftHeight ; i++)
            {
                _complexInput[i][1] = 0.0f;
                _complexInput[i][0] = 0.0f;
            }
            fftwf_execute( _plan );
            //£ float *sml;
            qint16 *sml;
            int jc=iCount*_nbo+iLoop;
            float b;
            //µ for(int i =0; i < _fftHeightHalf; i++)
            for(int i =0; i < _limY;  i++)
            {
                sml=_sonogramArray[i];
                a = pow(_fftHeight*_fftRes[i][0], 2) + pow(_fftHeight*_fftRes[i][1], 2);
                if (i > _freqMin && a !=0 )
                {
                    b=10*log10(a);
                    //£ sml[jc] = b;
                    sml[jc] = (qint16)(b*100.0f);
                    if(b>_energyMax) _energyMax=b;
                    else
                    {
                        if(b<_energyMin) _energyMin =b;
                    }
                }
                else 
				{
				//£ sml[jc] = -50;
				sml[jc] = -5000;
				}
            }
            iCount++;
        }
    }
    //µ for(int i =0; i < _fftHeightHalf; i++) _sonogramArray[i][_sonogramWidth-1]=0;
    for(int i =0; i < _limY; i++) _sonogramArray[i][_sonogramWidth-1]=0;
    _energyMax = qMax(_energyMax, (double)0);
    sf_close (_soundFile);
    fftwf_free(_data);
    fftwf_free(_fftRes);
    fftwf_free(_complexInput);
    delete _coeff;
    return true;
}

/*
void DetecTreatment::calculateMedianNoise()
{
    // _detec->_logText << "cmn debut _fftHeightHalf = " << _fftHeightHalf << endl;
    // _detec->_logText << "_maxY = " << _maxY << endl;
    _medianNoise = 0;
    int son_min = EMIN,son_max = EMIN+199;
    int tval[200];
    for(int i=0;i<200;i++) tval[i]=0;
    //£ float *fc;
    qint16 *fc;
    _energyMin = son_max;
    _energyMax = son_min;
    //µ for(int y = 0; y < _fftHeightHalf ; y++)
    for(int y = 0; y < _maxY ; y++)
    {
        //£ fc=_sonogramArray[y];
        fc=_sonogramArray[y];
        for (int x = 0 ; x < _sonogramWidth ; x++)
        {
            //£ int son = (int)fc[x];
            int son = (*fc++)/100;
            if(son>_energyMax) _energyMax = son;
            if(son<_energyMin) _energyMin = son;
            if(son>=son_min && son<=son_max) tval[son-son_min]++;
        }
    }
    int cumul = 0;
    //µ int moitie = (_sonogramWidth * _fftHeightHalf)/2;
    int moitie = (_sonogramWidth * _maxY)/2;
    bool oncherche = true;
    for(int i=0;i<=son_max-son_min;i++)
    {
        cumul += tval[i];
        if(cumul>=moitie && oncherche)
        {
            _medianNoise = son_min+i;
            oncherche=false; // provisoire
            break;
        }

    }
}
*/

void DetecTreatment::correctNoise()
{
    _minY=(int)(((float)_freqMin)/((float)_khzPerY));
    _maxY=(int)(((float)FREQ_MAX)/((float)_khzPerY));
    //µ if(_maxY>_fftHeightHalf-1) _maxY = _fftHeightHalf-1;
    if(_maxY>_limY-1) _maxY = _limY-1;

    int son_min = EMIN,son_max = EMIN+199;
    // _detec->_logText <<   "Correct noise 1ère partie :  "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    // 1) neutralisation des colonnes de silence
    /*
    float totCol;
    for (int x = 0 ; x < _sonogramWidth ; x++)
    {
        totCol = 0.0f;
        //for(int y = _minY; y < _maxY ; y++) totCol += _sonogramArray[y][x];
        for(int y = _minY; y < _maxY ; y++) totCol += ((float)_sonogramArray[y][x])/100.0f;
        _energyMoyCol[x] = totCol / (_maxY - _minY);
    }
    */
    for (int x = 0 ; x < _sonogramWidth ; x++) _energyMoyCol[x] = 0.0f;
    for(int y = _minY; y < _maxY ; y++)
    {
        qint16 *fc = _sonogramArray[y];
        for (int x = 0 ; x < _sonogramWidth ; x++) _energyMoyCol[x]  += (float)(*fc++);
    }
    float diviseur = 100.0f * ( (float) (_maxY - _minY) );
    float *pemc = _energyMoyCol;
    for (int x = 0 ; x < _sonogramWidth ; x++) *pemc++ /= diviseur;



    // //

    if(_sonogramWidth>10 && _useValflag)
    {
        bool valFlag = true;
        int patience = 0;
        int jumpThreshold = _jumpThreshold;
        int widthBigControl = _widthBigControl;
        int widthLittleControl = _widthLittleControl;
        int highThreshold = _highThreshold;
        int lowThreshold = _lowThreshold;
        for (int x = 0 ; x <= widthBigControl ; x++) _flagGoodCol[x]=valFlag;
        float dif2Col,averageLittleBefore,averageLittleNext,averageBigBefore,averageBigNext;
        bool notYet = true;
        if(widthBigControl > _sonogramWidth/3)
        {
            //_detec->_logText << "cas particulier : _sonogramWidth = " << _sonogramWidth << endl;
            widthBigControl = _sonogramWidth/3;
            if(widthLittleControl>widthBigControl/2) widthLittleControl=widthBigControl/2;
        }
        for (int x = widthBigControl+1 ; x < _sonogramWidth-widthBigControl ; x++)
        {
            averageLittleBefore = 0.0f;
            averageLittleNext   = 0.0f;
            for(int j=0;j<widthLittleControl;j++)
            {
                averageLittleNext   += _energyMoyCol[x+j];
                averageLittleBefore += _energyMoyCol[x-1-j];
            }
            averageLittleBefore /= widthLittleControl;
            averageLittleNext   /= widthLittleControl;
            dif2Col = averageLittleNext - averageLittleBefore;
            // -- -- --
            if(dif2Col>jumpThreshold && patience>widthBigControl/3
                    && averageLittleNext >highThreshold)
            {
                averageBigBefore = 0.0f;
                averageBigNext   = 0.0f;
                for(int j=0;j<widthBigControl;j++)
                {
                    averageBigNext   += _energyMoyCol[x+j];
                    averageBigBefore += _energyMoyCol[x-1-j];
                }
                averageBigBefore /= widthBigControl;
                averageBigNext   /= widthBigControl;
                if((averageBigNext-averageBigBefore)>jumpThreshold
                        && averageBigNext >highThreshold)
                {

                    patience = 0;
                    if(valFlag==false)
                    {
                    }

                    if(valFlag==true && notYet==true
                            && averageLittleBefore < lowThreshold
                            && averageBigBefore < lowThreshold)
                    {
                        for(int j=x-1;j>=0;j--)
                        {
                            if(_flagGoodCol[j]==true) _flagGoodCol[j]=false;
                            else break;
                        }
                    }
                    notYet =false;
                    valFlag=true;
                }
            }
            // -- -- --
            if(dif2Col<-jumpThreshold && patience>widthBigControl/3
                    && averageLittleNext<lowThreshold)
            {
                averageBigBefore = 0.0f;
                averageBigNext   = 0.0f;
                for(int j=0;j<widthBigControl;j++)
                {
                    averageBigNext   += _energyMoyCol[x+j];
                    averageBigBefore += _energyMoyCol[x-1-j];
                }
                averageBigBefore /= widthBigControl;
                averageBigNext   /= widthBigControl;
                if((averageBigNext-averageBigBefore)<-jumpThreshold
                        && averageBigNext<lowThreshold)
                {
                    patience = 0;
                    valFlag=false;
                    notYet =false;
                }
            }
            _flagGoodCol[x] = valFlag;
            patience ++;
        }
        for (int x = _sonogramWidth-widthBigControl;x < _sonogramWidth;x++) _flagGoodCol[x]=valFlag;
    }
    else
    {
        for (int x = 0;x < _sonogramWidth;x++) _flagGoodCol[x]=true;
    }
    // -----------------------------------------------------------------------------------
    // 2)
    // _detec->_logText <<   "Correct noise 2ème partie :  "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    int tval[200];
    //£ float *fc;
    qint16 *fc;
    int largeurRectifiee;
    for(int y = _minY; y <= _maxY ; y++)
    {
        fc=_sonogramArray[y];
        for(int k=0;k<200;k++) tval[k]=0;

        largeurRectifiee = 0;
        for (int x = 0 ; x < _sonogramWidth ; x++)
        {
            //£ int son = (int)fc[x];
            // int son = qRound((float)fc[x]/100.0f);
            int son = (fc[x]+50)/100;
            if(son>=son_min && son<=son_max && _flagGoodCol[x])
            {
                tval[son-son_min]++;
                largeurRectifiee++;
            }
            //£ if(son < son_min) fc[x] = son_min;
            //£ if(son > son_max) fc[x] = son_max;
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
                //£ ajouté :
                int retrait = (j + son_min + _stopThreshold) * 100;
                //£ for (int x = 0 ; x < _sonogramWidth ; x++) fc[x] = fc[x] - j - son_min - _stopThreshold;
                for (int x = 0 ; x < _sonogramWidth ; x++) fc[x] -= retrait;
                break;
            }
        }
    }
    //_detec->_logText <<   "Correct noise fin :  "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

}

void DetecTreatment::shapesDetects()
{
    //_detec->_logText <<   "shapesDetects 1:  "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    //! :
    //int ld8 = (SONOGRAM_WIDTH_MAX+15)/8;
    //£ for(int j=0;j<FFT_HEIGHT_HALF_MAX;j++) memset(_pointFlagsArray[j],0,SONOGRAM_WIDTH_MAX);
    //µ for(int j=0;j<FFT_HEIGHT_HALF_MAX;j++) memset(_pointFlagsArray[j],0,ld8);
    for(int j=0;j<MAXHEIGHT;j++) memset(_pointFlagsArray[j],0,LD8);
    _maxCallWidth = 0;
    _maxCallHeight = 0;
    //
    QPoint Point;
    int ix,iy;
    int curseur;
    //***_logText << "début DetectsContours" << endl ;
    int ncontour=0;
    // _energyShapeThreshold = (double)_detectionThreshold;
    // _energyStopThreshold = (double)_stopThreshold;
    _energyShapeThreshold = (double)_detectionThreshold-_stopThreshold;
    _energyStopThreshold = 0.0f;
    int nbcont=0;
    //£ float *smy;
    qint16 *smy;
    //§ char *zcy;
    bool onsarrete = false;
    //_detec->_logText <<   "shapesDetects 2:  "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    //_detec->_logText <<   "_maxY =  "<< _maxY  << "_minY =  "<< _minY  << endl;
    //_detec->_logText <<   "_sonogramWidth =  "<< _sonogramWidth << endl;
    for(int y = (_maxY-1); y >= _minY ; y--)
    {
        //_detec->_logText << "y="<<y<<endl;
        smy = _sonogramArray[y];
        //§ zcy=_pointFlagsArray[y];
        //§ :
        int digitPos = 0;
        char *pBoolChar = _pointFlagsArray[y];
        //char boolChar = *pBoolChar;
        for (int x = 0 ; x < _sonogramWidth ; x++)
        {
            if(((*pBoolChar) & (1 << digitPos))==0)
                // £if (smy[x] > _energyShapeThreshold)
                if (smy[x] > _energyShapeThreshold*100)
                {
                    nbcont++;
                    ncontour++;
                    _vectorCallPoints.clear();
                    Point.setX(x);
                    Point.setY(y);
                    //£ _callEnergyMax = smy[x];
                    _callEnergyMax = (float)smy[x] / 100.0f;
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
                            // edit yves - elargir spectre
                        {
                            if(jy>=_minY && jy <= _maxY)
                            {
                                int retx = 5;
                                if(ix<retx) retx = ix;
                                int posBC = (ix-retx)/8;
                                int dP = ix-retx-posBC*8;
                                char *pBC = _pointFlagsArray[jy]+posBC;
                                char bC = *pBC;
                                for(int jx=ix-retx;jx<=ix+5;jx++)
                                {
                                    if(jx!=ix || jy!=iy)
                                    {
                                        if(jx>= 0 && jx < _sonogramWidth)
                                        {
                                            if((bC & (1 << dP))==0)
                                            {
                                                float val = (float)_sonogramArray[jy][jx]/100.0f;
                                                if (val > _energyStopThreshold)
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
                                                } // fin if val
                                            } // fin if bc db...

                                        } // fin 2ème condition sur jx

                                    } // fin 1ère condition sur jx
                                    //§ :
                                    if(jx<_sonogramWidth-1)
                                    {
                                        dP++;
                                        if(dP==8) {pBC++; bC = *pBC; dP=0;}
                                    }
                                } // next jx
                            }

                        } // next jy
                        curseur++;
                        if(curseur > 300000)
                        {
                            _detec->_logText << "atteint la limite des 300000" << endl;
                            break;
                        }
                    } // next jy

                    float freqpm=((float)_vectorCallPoints.at(_callEnergyMaxIndex).y())*_khzPerY;

                    if(freqpm>_freqCallMin && (_xMax-_xMin+1)<MAXLARCRI
                            && (_yMax-_yMin+1)<MAXHEIGHT
                            && !onsarrete)
                    {
                        _masterPoints.push_back(_vectorCallPoints.at(_callEnergyMaxIndex));
                        _callsArray.push_back(_vectorCallPoints);
                        _vectorXMin.push_back(_xMin);
                        _vectorXMax.push_back(_xMax);
                        if((_xMax-_xMin+1)>_maxCallWidth) _maxCallWidth=_xMax-_xMin+1;
                        if((_yMax-_yMin+1)>_maxCallHeight) _maxCallHeight=_yMax-_yMin+1;
                        if(_callsArray.size()==MAXCRI-1) onsarrete = true;
                        if(_xMax-_xMin+1 > 5000)
                        {
                            _detec->_logText << "!!! LARCRI = " << _xMax-_xMin+1 << "   " << _wavFile << endl;
                        }
                    }
                }

            if(onsarrete) break;
            if(x<_sonogramWidth-1)
            {
                digitPos++;
                if(digitPos==8) {pBoolChar++; digitPos=0;}
            }
        } // next x
        if(onsarrete) break;
    }
    if(_callsArray.size() > 1500)
        _detec->_logText << "!!! NB CRIs = " << _callsArray.size() << "   " << _wavFile << endl;

    sortWaves();
}

void DetecTreatment::sortWaves()
{
    bool ontrie = true;
    int nbcris = _callsArray.size();
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
                QVector<QPoint> mvco1=_callsArray.at(i);
                QVector<QPoint> mvco2=_callsArray.at(i+1);
                _callsArray.replace(i,mvco2);
                _callsArray.replace(i+1,mvco1);
                QPoint pm1=_masterPoints.at(i);
                QPoint pm2=_masterPoints.at(i+1);
                _masterPoints.replace(i,pm2);
                _masterPoints.replace(i+1,pm1);
                _vectorXMin.replace(i,x2);
                _vectorXMin.replace(i+1,x1);
                _vectorXMax.replace(i,x4);
                _vectorXMax.replace(i+1,x3);
                ontrie = true;
            }
        }
    }
}

void DetecTreatment::detectsParameter2()
{
    float *oParam;
    float *oParamCrete[NCRETES];
    //_detec->_logText << "dp2-debut";
    int nbcris = _callsArray.size();
    if(nbcris< 1 || _maxCallWidth < 1 || _maxCallHeight < 1) return;
    int maxlarhau = _maxCallWidth;
    if(_maxCallHeight>maxlarhau) maxlarhau=_maxCallHeight;
    //_logText << "-2" ;
    //
    // 23/3/2015 :
    int *sortMp = new int[nbcris];
    int *invMp = new int[nbcris];
    int *xMp = new int[nbcris];
    for(int i=0;i<nbcris;i++) {sortMp[i]=i; xMp[i]=_masterPoints.at(i).x();}
     sortIntArrays(sortMp,nbcris,xMp);
     for(int i=0;i<nbcris;i++)  invMp[sortMp[i]]=i;
     //
    for (int icri = 0 ; icri < nbcris ; icri++) //Execute for each call
    {
        //_logText << endl << "-3-icri=" << icri ;
        oParam = _paramsArray[icri][SH];
        for(int j=0;j<NCRETES;j++) oParamCrete[j] = _paramsArray[icri][j+1];
        //_logText << "-4";
        QVector<QPoint> unemat = _callsArray.at(icri);
        int tailleforme = unemat.size(); // The number of vectors of each call
        int xmin = _sonogramWidth-1;
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
        float eTot = 0.0f;
        float eMoy = 0.0f;
        for(int j=0;j<tailleforme;j++)
        {
            int x=unemat.at(j).x();
            int y=unemat.at(j).y();
            //£ float e = _sonogramArray[y][x];
            float e = (float)_sonogramArray[y][x]/100.0f;
            _tabY[y-ymin]+=e;
            eTot += e;
            _numberPixelsPerY[y-ymin]++;
            _tabX[icri][x-xmin]+=e;
            _numberPixelsPerX[x-xmin]++;
            _totYEPerX[x-xmin]+=e*(float)y;
            if(y>_yMaxPerX[x-xmin]) _yMaxPerX[x-xmin]=y;
            if(y<_yMinPerX[x-xmin] || _yMinPerX[x-xmin]==0) _yMinPerX[x-xmin]=y;
            if(x>_xMaxPerY[y-ymin]) _xMaxPerY[y-ymin]=x;
            if(x<_xMinPerY[y-ymin]) _xMinPerY[y-ymin]=x;
            if(e>_eMaxPerX[x-xmin])
            {
                _eMaxPerX[x-xmin]=e;
                _yEmaxPerX[icri][x-xmin]=y;
            }
            _tabYX[y-ymin][x-xmin]=e;

        }
        eMoy = eTot / tailleforme;
        float erec; int xc;
        for(int k=0;k<=ymax-ymin;k++)
        {
            xc=_xMinPerY[k];
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
            _xSecondWestRidgePerY[k]=xc;
        }
        int yfds = 0;
        int yfdm = 0;
        float recfds= 0.0f;
        float recfdm= 0.0f;
        for(int k=0;k<=ymax-ymin;k++)
        {
            if(_tabY[k]>recfds)
            {
                yfds=ymin+k;
                recfds=_tabY[k];
            }
            if(_numberPixelsPerY[k]>0)
                if(_tabY[k]/_numberPixelsPerY[k]>recfdm)
                {
                    yfdm=ymin+k;
                    recfdm=_tabY[k]/_numberPixelsPerY[k];
                }
        }
        int xpps = 0;
        int xppm = 0;
        float recpps= 0.0f;
        float recppm= 0.0f;
        for(int k=0;k<=xmax-xmin;k++)
        {
            if(_tabX[icri][k]>recpps)
            {
                xpps=xmin+k;
                recpps=_tabX[icri][k];
            }
            if(_numberPixelsPerX[k]>0)
            {
                if(_tabX[icri][k]/_numberPixelsPerX[k]>recppm)
                {
                    xppm=xmin+k;
                    recppm=_tabX[icri][k]/_numberPixelsPerX[k];
                }
                _yEbarPerX[k]=_totYEPerX[k]/_tabX[icri][k];
            }

        }
        oParam[StTime] = (float)xmin*(float)_msPerX;
        oParam[Dur] = (xmax-xmin)*(float)_msPerX;
        int criprec=icri;
        if(icri>0)
            for(int jcri=icri-1;jcri>=0;jcri--)
            {
                if(_paramsArray[jcri][SH][StTime] + _paramsArray[jcri][SH][Dur]<oParam[StTime])
                {
                    criprec=jcri;
                    break;
                }
            }
        if(criprec==icri) oParam[PrevSt] = 9999.0f;
        else oParam[PrevSt] = oParam[StTime] - _paramsArray[criprec][SH][StTime];
        oParam[Fmax]=ymax*(float)_khzPerY;
        oParam[Fmin]=ymin*(float)_khzPerY;
        oParam[BW]=oParam[Fmax]-oParam[Fmin];
        oParam[FreqMP] = (float)_masterPoints.at(icri).y()*_khzPerY;
        oParam[PosMP] = (float)(_masterPoints.at(icri).x()-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[FreqPkS]  = yfds * _khzPerY;
        oParam[FreqPkM] = yfdm * _khzPerY;
        oParam[PosPkS]   = (float)(xpps-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[PosPkM]  = (float)(xppm-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[FreqPkS2] = (float) _yEmaxPerX[icri][xpps-xmin] * _khzPerY;
        oParam[FreqPkM2] = (float) _yEmaxPerX[icri][xppm-xmin] * _khzPerY;
        float prevmp = 9999.0f;
        float prevsmart1 =9999.0f;
        //
        //_logText << "-5";

        /*
        if(icri>0)
        {
            float freqmp = oParam[FreqMP];
            if(criprec<icri)
            {
                prevmp=(_masterPoints.at(icri).x()-_masterPoints.at(criprec).x())*_msPerX;
                for(int k=criprec;k>=0;k--)
                {
                    if(_paramsArray[k][SH][Fmin] < freqmp
                            && _paramsArray[k][SH][Fmax] > freqmp)
                    {
                        prevsmart1 = (_masterPoints.at(icri).x()-_masterPoints.at(k).x())*_msPerX;
                        break;
                    }
                }
            }

        }
        */
        int pmp = invMp[icri];
        if(pmp>0)
        {
            float freqmp = oParam[FreqMP];
            prevmp=(_masterPoints.at(icri).x()-_masterPoints.at(sortMp[pmp-1]).x())*_msPerX;
            for(int k=pmp-1;k>=0;k--)
            {
                if(_paramsArray[sortMp[k]][SH][Fmin] < freqmp
                        && _paramsArray[sortMp[k]][SH][Fmax] > freqmp)
                {
                    prevsmart1 = (_masterPoints.at(icri).x()-_masterPoints.at(sortMp[k]).x())*_msPerX;
                    break;
                }
            }
        }
        oParam[PrevMP1] = prevmp;
        oParam[PrevMP2] = prevsmart1;
        oParam[NextMP1] = 9999.0f;
        oParam[NextMP2] = 9999.0f;
        float famp[4];
        float namp[4];
        for(int j=0;j<4;j++)
        {
            famp[j] = 0.0f;
            namp[j] = 0.0f;
        }
        float prorata,bdeb,bfin,larco,lartot;
        int ideb,ifin;
        //_logText << "-6" ;
        lartot=(float)(xmax-xmin+1);
        for(int k=0;k<xmax-xmin+1;k++)
        {
            bdeb=(((float)k)*4.0f)/lartot;
            bfin=(((float)(k+1))*4.0f)/lartot;
            ideb=(int)bdeb;
            ifin=(int)bfin;
            larco=bfin-bdeb;
            prorata = larco/lartot;
            for(int j=ideb;j<=ifin;j++)
            {
                if(j==4) break;
                if(j==ideb)
                {
                    if(ifin==ideb) prorata=1.0f/larco;
                    else prorata = ((float)(j+1)-bdeb)/larco;
                }
                else
                {
                    if(bfin < (j+1)) prorata=(bfin-(float)j)/larco;
                    else prorata = 1.0f/larco;
                }
                famp[j]+=_tabX[icri][k]*prorata;
                namp[j]+=prorata;
            }
        }
        // modifié le 23/3/2015
        if(namp[0]>0) oParam[Amp1] = famp[0]*_khzPerY/namp[0];
        else oParam[Amp1] = 0.0f;
        if(namp[1]>0) oParam[Amp2] = famp[1]*_khzPerY/namp[1];
        else oParam[Amp2] = 0.0f;
        if(namp[2]>0) oParam[Amp3] = famp[2]*_khzPerY/namp[2];
        else oParam[Amp3] = 0.0f;
        if(namp[3]>0) oParam[Amp4] = famp[3]*_khzPerY/namp[3];
        else oParam[Amp4] = 0.0f;
        for(int inoise=0;inoise<4;inoise++)
        {
            int jdeb,jfin,x,y;
            int nps=0;
            float bruit=0.0f,bruit_moyen=0.0f;
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
                        if(x>=0 && x<_sonogramWidth && y>0 && y<_fftHeightHalf)
                        {
                            //£ bruit += _sonogramArray[y][x];
                            //bruit += (float)_sonogramArray[y][x]/100.0f;
                            bruit += (float)_sonogramArray[y][x];
                            nps++;
                        }
                    }
                }
            }
            if(nps>0) bruit_moyen=bruit/ ((float)nps*100.0f);
            if(inoise==0)oParam[NoisePrev]  = bruit_moyen;
            if(inoise==1)oParam[NoiseNext] = bruit_moyen;
            if(inoise==2)oParam[NoiseDown]  = bruit_moyen;
            if(inoise==3)oParam[NoiseUp]    = bruit_moyen;
        }
        //_logText << "-7";
        int nms=0;
        float sdm=0.0f;
        for(int k=0;k<xmax-xmin+1;k++)
        {
            if(_numberPixelsPerX[k]>0)
            {
                _averagePerX[k]=_tabX[icri][k]/_numberPixelsPerX[k];
                sdm+=_averagePerX[k];
                nms++;
            }
        }
        float mdm = sdm/nms;
        float sdc = 0.0f;
        for(int k=0;k<xmax-xmin+1;k++)
        {
            if(_numberPixelsPerX[k]>0) sdc += (_averagePerX[k]-mdm)*(_averagePerX[k]-mdm);
        }
        oParam[CVAmp]=pow(sdc/(float)nms,0.5);
        int tempsOuest[2];
        tempsOuest[0] = qMax(_xMinPerY[0],_xMinPerY[ymax-ymin])-xmin+1;
        oParamCrete[3][Dur] = (float)tempsOuest[0]*(float)_msPerX;

        tempsOuest[1] = qMax(_xSecondWestRidgePerY[0],_xSecondWestRidgePerY[ymax-ymin])-xmin+1;
        oParamCrete[4][Dur] = (float)tempsOuest[1]*(float)_msPerX;
        int fmin=_yEmaxPerX[icri][0],fmax=fmin;
        for(int k=1;k<=xmax-xmin;k++)
        {
            if(_yEmaxPerX[icri][k]>0)
            {
                if(_yEmaxPerX[icri][k]>fmax) fmax=_yEmaxPerX[icri][k];
                if(_yEmaxPerX[icri][k]<fmin) fmin=_yEmaxPerX[icri][k];
            }
        }
        oParamCrete[0][Fmax]=(float)fmax*(float)_khzPerY;
        oParamCrete[0][Fmin]=(float)fmin*(float)_khzPerY;
        oParamCrete[0][BW]=oParamCrete[0][Fmax]-oParamCrete[0][Fmin];
        //
        //£int *pc; // pointeur sur la crete
        quint16 *pc; // pointeur sur la crete
        int alasuite,meilleuresuite,meilleur;
        float meilleurepente;
        int xmcmax[NCRETES],ymcmax[NCRETES];
        float amp[10],ampmoy[10],ran[10],ranmoy[10];
        int num[10],denom[10];
        num[0]=1; denom[0]=1;
        num[1]=2; denom[1]=1;
        num[2]=3; denom[2]=1;
        num[3]=3; denom[3]=2;
        num[4]=1; denom[4]=2;
        num[5]=4; denom[5]=3;
        num[6]=2; denom[6]=3;
        int np,y,xn,nn,yn;
        amp[0]=0.0f;

        for(int x=xmin;x<=xmax;x++) amp[0]+=_tabX[icri][x-xmin];
        ampmoy[0]=amp[0]/((float)tailleforme);
        for(int ih=1;ih<=6;ih++)
        {
            amp[ih]=0.0f;
            ampmoy[ih]=0.0f;
            np=0;
            ran[ih]=0.0f;
            ranmoy[ih]=0.0f;
            nn=0;
            for(int j=ymin;j<=ymax;j++)
            {
                if(ih>0) y=(j*num[ih])/denom[ih];
                if(y>_maxY || y<1) continue;
                for(int x=_xMinPerY[j-ymin];x<=_xMaxPerY[j-ymin];x++)
                {
                    if(_tabYX[j-ymin][x-xmin]>0.0f)
                    {
                        //£ amp[ih] += _sonogramArray[y][x];
                        //amp[ih] += (float)_sonogramArray[y][x]/100.0f;
                        amp[ih] += (float)_sonogramArray[y][x];
                        np++;
                    }
                }
                for(int gd=0;gd<2;gd++) for(int k=1;k<=3;k++)
                {
                    if(gd==0)xn=_xMinPerY[j-ymin]-k;
                    else xn=_xMaxPerY[j-ymin]+k;
                    if(xn>=0 && xn<_sonogramWidth)
                    {
                        //£ ran[ih]+=_sonogramArray[y][xn];
                        //ran[ih] += (float)_sonogramArray[y][xn]/100.0f;
                        ran[ih] += (float)_sonogramArray[y][xn];
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
                            //£ ran[ih]+=_sonogramArray[yn][x];
                            //ran[ih] += (float)_sonogramArray[yn][x]/100.0f;
                            ran[ih] += (float)_sonogramArray[yn][x];
                            nn++;
                        }
                    }
                }
            }
            //if(np>0) ampmoy[ih]=amp[ih]/((float)np);
            //if(nn>0) ranmoy[ih]=ran[ih]/((float)nn);
            if(np>0) ampmoy[ih]=amp[ih]/((float)np*100.0f);
            if(nn>0) ranmoy[ih]=ran[ih]/((float)nn*100.0f);
        } // next ih
        if(ampmoy[0]>0.0f)
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
        //
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
        //
        int xmaitr = _masterPoints.at(icri).x();
        int ymaitr = _masterPoints.at(icri).y();

        // ---------------------------------------------------
        //_logText << "Avant calcul des hetx, hety"  << endl;
        int npar,ntr,nen,xdeb,xfin,ydeb,yfin;
        float e1,e2,e3;
        for(int r=0;r<2;r++)
        {
            if(r==0) npar = HetX; else npar = HetXr;
            oParam[npar] = 0.0f;
            ntr = 0;
            nen = 0;
            float e1,e2,e3;
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
                    if(e2>0.0f && x>0 && x<_sonogramWidth-1)
                    {
                        ntr++;
                        //£ e1=_sonogramArray[y][x-1];
                        e1 = (float)_sonogramArray[y][x-1]/100.0f;
                        //£ e3=_sonogramArray[y][x+1];
                        e3 = (float)_sonogramArray[y][x+1]/100.0f;
                        if((e2>e1 && e2>e3) || (e2<e1 && e2<e3)) nen++;
                    }
                }
            }
            if(ntr>0) oParam[npar] = ((float) nen)/ ((float) ntr);
        }
        //
        for(int r=0;r<3;r++)
        {
            float dift = 0.0f;
            if(r==0) npar = HetY;
            else
            {
                if(r==1) npar = HetYr;
                else npar = HetYr2;
            }
            oParam[npar] = 0.0f;
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
                        if(e2>0.0f && y>1 && y<_maxY)
                        {
                            ntr++;
                            //£ e1=_sonogramArray[y-1][x];
                            e1=(float)_sonogramArray[y-1][x]/100.0f;
                            //£ e3=_sonogramArray[y+1][x];
                            e3=(float)_sonogramArray[y+1][x]/100.0f;
                            dift+=qAbs(e3-e2);
                            if((e2>e1 && e2>e3) || (e2<e1 && e2<e3)) nen++;
                        }
                    }
                }
            }
            if(ntr>0)
            {
                if(npar<2) oParam[npar] = ((float) nen)/ ((float) ntr);
                else oParam[npar] = (float)dift/(float)ntr;
            }
        }
        //_logText << "Après calcul des hetx, hety"  << endl;
        // ---------------------------------------------------
        // HetCMC et HeCMD et HetCTC et HetCTD
        // ajouté recherche des pics et creux
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
        //int creux[2][(MAXLARCRI+1)/2];
        int nbpics[2];
        //int nbcreux[2];
        int modepc;
        int memopr;
        float seuilpc[2];
        seuilpc[0] = 1.0f; seuilpc[1] = 10.0f;
        float valpr,difppr;
        for(int mt=0;mt<2;mt++)
        {
            int sensd=1,dersensd = 0;
            int nch=0,ncc=0;
            float moytot= 0.0f,dermoytot = 0.0f;
            float difC=0.0f,totDifC=0.0f,derDifC=0.0f;
            nbpics[mt]=0;
            //nbcreux[mt]=0;
            modepc=0;
            memopr=0;
            valpr=0.0f;
            for(int x=xmin;x<=xmax;x++)
            {
                if(_numberPixelsPerX[x-xmin]>0)
                {
                    if(mt==0) moytot = _averagePerX[x-xmin]; else moytot= _tabX[icri][x-xmin];
                    if(x==xmin) valpr = moytot;
                }
                else moytot = 0.0f;
                //
                difppr=moytot-valpr;
                if(difppr>0)
                {
                    if(difppr>=seuilpc[mt])
                    {
                        if(modepc==0 || modepc==2)
                        {
                            //creux[mt][nbcreux[mt]++]=memopr;
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
                    if(difppr<0)
                    {
                        if(qAbs(difppr)>=seuilpc[mt])
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
                    if(modepc==1) pics[mt][nbpics[mt]++] = memopr;
                    //if(modepc==2) creux[mt][nbcreux[mt]++] = memopr;

                }
                // ---
                difC = moytot - dermoytot;
                if(difC>0.0f) sensd = 1;
                else
                {
                    if(difC < 0.0f) sensd = 0;
                }
                if(difC!=0 || derDifC!=0) ncc++;
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
                    oParam[HetCMC] = ((float) nch)/ ((float) ncc);

                    oParam[HetCMD] = (totDifC)/ ((float) ncc);
                }
                else
                {
                    oParam[HetCTC] = ((float) nch)/ ((float) ncc);

                    oParam[HetCTD] = (totDifC)/ ((float) ncc);
                }

            }
        }
        oParam[HetCMnP] = nbpics[0];
        oParam[HetCMfP] = (float)nbpics[0]/((xmax-xmin+1)*_msPerX);
        oParam[HetCTnP] = nbpics[1];
        oParam[HetCTfP] = (float)nbpics[1]/((xmax-xmin+1)*_msPerX);
        // -----------------------------------------------------------------------------------
        float interc[MAXLARCRI/2];
        float variationPicsInter[MAXLARCRI/2];
        for(int mt=0;mt<2;mt++)
        {
            float medianPicsDistance = ((float)(xmax-xmin+2)/2)*_msPerX;
            float medianPicsLittleDistance =((float)(xmax-xmin+2)/2)*_msPerX;
            float medianPicsBigDistance = ((float)(xmax-xmin+2)/2)*_msPerX;
            int nbi = nbpics[mt]-1;
            int halfnbi = nbi/2;
            int quartnbi =halfnbi/2;
            if(nbi>1)
            {
                for(int j=0;j<nbi;j++) interc[j] = (qAbs((pics[mt][j+1]-pics[mt][j])))*_msPerX;
                sortFloatArray(interc,nbi);
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
                //
            } // fin if nbi>1
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
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            float medianPicsDistanceVariation = 0;
            // a) écart médian
            for(int j=0;j<nbi;j++) variationPicsInter[j] = qAbs(interc[j]-medianPicsDistance);
            float medianPicsLittleDistanceVariation = variationPicsInter[0];
            float medianPicsBigDistanceVariation = variationPicsInter[nbi-1];
            if(nbi>1)
            {
                sortFloatArray(variationPicsInter,nbi);
                if(halfnbi*2<nbi) medianPicsDistanceVariation = variationPicsInter[halfnbi];
                else medianPicsDistanceVariation = (variationPicsInter[halfnbi-1] + variationPicsInter[halfnbi])/2;
                // b) écarts sur petits intervalles et sur grands intervalles
                for(int j=0;j<nbi;j++)
                {
                    if(j<halfnbi) variationPicsInter[j] = qAbs(interc[j]-medianPicsLittleDistance);
                    else variationPicsInter[j] = qAbs(interc[j]-medianPicsBigDistance);
                }
                if(nbi<4)
                {
                    medianPicsLittleDistanceVariation = variationPicsInter[0];
                    medianPicsBigDistanceVariation = variationPicsInter[nbi-1];
                }
                else
                {
                    sortFloatArray(variationPicsInter,halfnbi);
                    sortFloatArray(variationPicsInter+nbi-halfnbi,halfnbi);
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
                } // fin nbi>=4
            } // fin if nbi>1
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
        } // next mt
        // -----------------------------------------------------------------------------------
        //_logText << "-8" ;
        oParam[Dbl8] = 0.0f;
        float e8 = 0.0f;
        int n8 = 0;
        int Y8 = (int)(_freqCallMin / _khzPerY);
        for(int y=0;y<Y8;y++)
        {
            qint16 *fc=_sonogramArray[y];
            for(int x=xmin;x<=xmax;x++)
            {
                //£ e8 =+ _sonogramArray[y][x];
                e8 =+ (float)(*fc++);
                n8++;
            }
        }
        oParam[Dbl8] = eMoy - (e8/((float)n8*100.0f));
        // ùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùù
        oParam[Stab] = 0.0f;
        float dif,p;
        float ponderTot = 0.0f;
        float difTot = 0.0f;
        float dxmax = ((float)qMax(xmaitr-xmin,xmax-xmaitr))*_msPerX;
        float dymax = ((float)qMax(ymaitr-ymin,ymax-ymaitr))*_khzPerY;
        float distLim = pow(dxmax,2)+pow(dymax,2);
        //
        if(distLim>0.0f)
        {
            for(int y=ymin;y<=ymax;y++)
            {
                for(int x=_xMinPerY[y-ymin]-1;x<=_xMaxPerY[y-ymin];x++)
                {
                    if(x>=0 && x<_sonogramWidth-1)
                    {
                        //£ dif = qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]);
                        //dif = (float)(qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]))/100.0f;
                        dif = (float)(qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]));
                        p = 1.0f - (((    pow(   ((float)(x-xmaitr))*_msPerX , 2)
                                          +pow(   ((float)(y-ymaitr))*_khzPerY,2))/distLim) * 0.75f);
                        difTot += dif * p;
                        ponderTot += p;
                    }
                }
            }
            // _logText << "1) diftot=" << difTot << "  -  pondertot=" << ponderTot << endl;
            for(int x=xmin;x<=xmax;x++)
            {
                if(_yMaxPerX[x-xmin]>0)
                {
                    for(int y=_yMinPerX[x-xmin]-1;y<=_yMaxPerX[x-xmin];y++)
                    {
                        if(y>=_minY && y<_maxY-1)
                        {
                            //£ dif = qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]);
                            //dif = (float)(qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]))/100.0f;
                            dif = (float)(qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]));
                            p = 1.0f - (((    pow(   ((float)(x-xmaitr))*_msPerX , 2)
                                              +pow(   ((float)(y-ymaitr))*_khzPerY,2))/distLim) * 0.75f);
                            difTot += dif * p;
                            ponderTot += p;
                        }
                    }
                }
            }
            //if(ponderTot!=0) oParam[Stab] = difTot / ponderTot;
            if(ponderTot!=0) oParam[Stab] = difTot / (ponderTot*100.0f);
        }
        //_logText << "Avant calcul stablr et stabbr"  << endl;
        oParam[EnStabSm] = 0.0f;
        oParam[EnStabLg] = 0.0f;
        int radius,nbp;
        //
        for(int lbr=0;lbr<2;lbr++)
        {
            difTot = 0.0f;
            nbp=0;
            if(lbr==0) radius = 3; else radius = 10;
            ydeb = ymaitr-radius; yfin = ymaitr + radius;
            xdeb = xmaitr - radius ; xfin = xmaitr + radius;
            if(ydeb < _minY) ydeb = _minY;
            if(yfin > _maxY-1) yfin = _maxY-1;
            if(xdeb < 0) xdeb = 0;
            if(xfin > _sonogramWidth-1) xfin = _sonogramWidth-1;
            for(int y=ydeb;y<=yfin;y++)
            {
                for(int x=xdeb;x<xfin;x++)
                {
                    //£ difTot += qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]);
                    //difTot += (float)(qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]))/100.0f;
                    difTot += (float)(qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]));
                    nbp++;
                }
            }
            for(int x=xdeb;x<=xfin;x++)
            {
                for(int y=ydeb;y<yfin;y++)
                {
                    //£ difTot += qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]);
                    //difTot += (float)(qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]))/100.0f;
                    difTot += (float)(qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]));
                    nbp++;
                }
            }
            if(lbr==0) npar = EnStabSm; else npar = EnStabLg;
            //if(nbp>0) oParam[npar] = difTot/nbp;
            if(nbp>0) oParam[npar] = difTot/((float)nbp*100.0f);
        }
        //_logText << "Après calcul stablr et stabbr"  << endl;
        float enerMaster[2];
        enerMaster[0] = (float) _tabY[ymaitr-ymin];
        enerMaster[1] = (float) _tabY[ymaitr-ymin]/_numberPixelsPerY[ymaitr-ymin];
        int enerHeight[2][2];
        float enerY[2];
        int incr;
        for(int tm=0;tm<2;tm++)
        {
            for(int sens = 0;sens <=1;sens++)
            {
                enerHeight[tm][sens] = 0;
                if(sens) incr = 1; else incr = -1;
                for(int y=ymaitr+incr;y>=ymin && y<=ymax;y+=incr)
                {
                    if(tm==0) enerY[tm] = (float)_tabY[y-ymin];
                    else enerY[tm] = (float) _tabY[y-ymin] / _numberPixelsPerY[y-ymin];
                    if(enerY[tm] > enerMaster[tm] * 0.8f) enerHeight[tm][sens]++;
                    else break;
                }
            }
        }
        oParam[HeiET] = (float) (enerHeight[0][0]+enerHeight[0][1]);
        oParam[HeiEM] = (float) (enerHeight[1][0]+enerHeight[1][1]);
        oParam[HeiRT] = oParam[HeiET]/(ymax - ymin+1);
        oParam[HeiRM] = oParam[HeiEM]/(ymax - ymin+1);
        oParam[HeiETT] = (float) enerHeight[0][1];
        oParam[HeiEMT] = (float) enerHeight[1][1];
        oParam[HeiRTT] = oParam[HeiETT]/(ymax - ymin+1);
        oParam[HeiRMT] = oParam[HeiEMT]/(ymax - ymin+1);
        //_logText << "Après calcul des Hei..."  << endl;
        // ùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùùù
        if(_imageData)
        {
            _callMasterRidge.clear();
            _callSouthRidge.clear();
            _callNorthRidge.clear();
            _callWestRidge.clear();
            _callSecondWestRidge.clear();
        }
        int lowfc,hifc,fc3;
        float yCmax=0.0f;
        float yCmin=9999.0f;
        //
        for(int jcrete=0;jcrete<NCRETES;jcrete++)
        {
            if(jcrete==0)  pc=_yEmaxPerX[icri];
            if(jcrete==1)  pc=_yMinPerX;
            if(jcrete==2)  pc=_yMaxPerX;
            if(jcrete==3)  pc=_xMinPerY;
            if(jcrete==4)  pc=_xSecondWestRidgePerY;
            oParamCrete[jcrete][Slope] = 9999.0f;
            float emcmax = -50.0f;
            xmcmax[jcrete] = 0;
            ymcmax[jcrete]=0;
            if(jcrete<3)
            {
                alasuite = 0; meilleuresuite = 0; meilleur = xmin; meilleurepente = 100000.0f;
                yCmax=0.0f;
                yCmin=9999.0f;
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
                            _slope[x-xmin] = qAbs(((float)pc[x-xmin+2] - (float)pc[x-xmin-2])/4.0f);
                    if(_slope[x-xmin]==0.0f)
                    {
                        meilleurepente = 0;
                        alasuite++;
                        if(alasuite > 1 && (meilleur == xmin || alasuite > meilleuresuite))
                        {meilleur = x; meilleuresuite = alasuite;}
                    }
                    else
                    {
                        alasuite = 0;
                        if(_slope[x-xmin] < meilleurepente)
                        {
                            meilleur = x;
                            meilleurepente = _slope[x-xmin];
                        }
                    }
                    //£ float e=_sonogramArray[pc[x-xmin]][x];
                    float e=(float)_sonogramArray[pc[x-xmin]][x]/100.0f;
                    if(e>emcmax)
                    {
                        emcmax = e;
                        xmcmax[jcrete]=x;
                        ymcmax[jcrete] = pc[x-xmin];
                    }
                } // fin boucle sur x
                oParamCrete[jcrete][Fmax] = yCmax * _khzPerY;
                oParamCrete[jcrete][Fmin]  = yCmin  * _khzPerY;
                oParamCrete[jcrete][BW]  = oParamCrete[jcrete][Fmax]-oParamCrete[jcrete][Fmin];
                //
                if(xmax>xmin) oParamCrete[jcrete][Slope] = (float)(pc[xmax-xmin]-pc[0])*_khzPerY/((float)(xmax-xmin)*_msPerX);
                if(meilleurepente ==0 && meilleuresuite>1)
                {
                    int bonmeilleur = meilleur;
                    for(int k=1;k<meilleuresuite/2;k++)
                        if(pc[meilleur-k-xmin]>0) bonmeilleur = meilleur-k;
                    meilleur = bonmeilleur;
                }

                lowfc = meilleur-xmin;
                _lowSlope[(icri*NCRETES+jcrete)*2] = meilleur; // (xmin+lowfc)
                _lowSlope[(icri*NCRETES+jcrete)*2+1] = pc[lowfc];
                oParamCrete[jcrete][FIF]= pc[lowfc]*_khzPerY;
                float pe1=0.0f,pe2=0.0f,dif;
                float meilleuredist;
                int meilleurk[2];
                int debb,finb;
                meilleurk[0]=0;
                meilleurk[1]=lowfc;
                for(int jp=0;jp<2;jp++)
                {
                    meilleurk[jp]=0;
                    meilleuredist = 1000000.0f;
                    if(jp==0) {debb=0;finb=lowfc;}
                    else {debb=lowfc;finb=xmax-xmin;}
                    if((jp==0 && lowfc > 1) || (jp==1 && lowfc < xmax - xmin-1))
                    {
                        for (int k = debb+1 ; k<= finb-1 ; k++)
                        {
                            dif=0.0f;
                            if(pc[k]>0)
                            {
                                pe1=(float)(pc[k]-pc[debb])/((float)k);
                                pe2=((float)(pc[finb]-pc[k]))/((float)(lowfc-k));
                            }
                            for(int s = debb+1 ; s<k ; s++)
                                if(pc[s]>0)
                                    dif+=(float)pow(pe1*((float)s)+(float)pc[debb]-(float)pc[s],2);

                            for(int s = k+1 ; s<=finb-1 ; s++)
                                if(pc[s]>0)
                                    dif+=(float)pow(pe2*((float)(finb-s))+(float)pc[s]-(float)pc[finb],2);
                            if(dif<meilleuredist)
                            {
                                meilleurk[jp]=k;
                                meilleuredist=dif;
                            }
                        }

                    }
                }
                hifc=meilleurk[0];
                _inflexion1[(icri*NCRETES+jcrete)*2]=hifc+xmin;
                _inflexion1[(icri*NCRETES+jcrete)*2+1]=pc[hifc];
                oParamCrete[jcrete][HCF] = pc[hifc]*_khzPerY;
                oParamCrete[jcrete][THCF]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][THCF] = ((float)meilleurk[0]+0.5f)/((float)(xmax-xmin+1));
                fc3=meilleurk[1];
                _inflexion3[(icri*NCRETES+jcrete)*2]=fc3+xmin;
                _inflexion3[(icri*NCRETES+jcrete)*2+1]=pc[fc3];
                oParamCrete[jcrete][LCF] = pc[fc3]*_khzPerY;
                oParamCrete[jcrete][UpSl] = 0.0f;
                if(meilleurk[0]>0)
                    oParamCrete[jcrete][UpSl] = ((float)(pc[meilleurk[0]]-pc[0])*_khzPerY)/((float)meilleurk[0]*_msPerX);
                oParamCrete[jcrete][LoSl] = 0.0f;
                if(lowfc-meilleurk[0]>0)
                    oParamCrete[jcrete][LoSl] = ((float)(pc[lowfc]-pc[meilleurk[0]])*_khzPerY)/((float)(lowfc-meilleurk[0])*_msPerX);
                int debp,finp,ecart;
                float pj[4];
                for(int j2=0;j2<4;j2++)
                {
                    pj[j2]=9999.0f;
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
                        pj[j2]=((float)(pc[finp]-pc[debp])*_khzPerY)/((float)ecart*_msPerX);
                    }
                }
                oParamCrete[jcrete][StSl] = pj[0];
                oParamCrete[jcrete][EnSl]   = pj[1];
                oParamCrete[jcrete][FISl]  = pj[3];
                oParamCrete[jcrete][FPSl]  = pj[2];
                //
                // RAR : début
                if(_withNewParams)
                {
                    oParamCrete[jcrete][SDC]  = 0;
                    oParamCrete[jcrete][SDCR]  = 0;
                    //
                    oParamCrete[jcrete][SDCL]  = 0;
                    oParamCrete[jcrete][SDCLR]  = 0;
                    oParamCrete[jcrete][SDCLRY]  = 0;
                    oParamCrete[jcrete][SDCLRXY]  = 0;
                    oParamCrete[jcrete][SDCLRXY2]  = 0;
                    //
                    oParamCrete[jcrete][SDCLOP]  = 0;
                    oParamCrete[jcrete][SDCLROP]  = 0;
                    oParamCrete[jcrete][SDCLRYOP]  = 0;
                    oParamCrete[jcrete][SDCLRXYOP]  = 0;
                    //
                    oParamCrete[jcrete][SDCLWB]  = 0;
                    oParamCrete[jcrete][SDCLRWB]  = 0;
                    oParamCrete[jcrete][SDCLRYWB]  = 0;
                    oParamCrete[jcrete][SDCLRXYWB]  = 0;
                    //
                    oParamCrete[jcrete][SDCLOPWB]  = 0;
                    oParamCrete[jcrete][SDCLROPWB]  = 0;
                    oParamCrete[jcrete][SDCLRYOPWB]  = 0;
                    oParamCrete[jcrete][SDCLRXYOPWB]  = 0;
                    //
                    //
                    oParamCrete[jcrete][SDCL_DNP]  = 0;
                    oParamCrete[jcrete][SDCLR_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRY_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRXY_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRXY2_DNP]  = 0;
                    //
                    oParamCrete[jcrete][SDCLOP_DNP]  = 0;
                    oParamCrete[jcrete][SDCLROP_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRYOP_DNP]  = 0;
                    oParamCrete[jcrete][SDCLRXYOP_DNP]  = 0;

                    float derdif1=0.0f;
                    float dif1,dif2;
                    float totdif2 = 0.0f;
                    float totdif2Op = 0.0f;
                    float totdif2wb = 0.0f;
                    float totdif2opwb = 0.0f;
                    // -----------------------------------------------------------------------------------------
                    // Paramètres de la série SDC
                    // float pente = 0;
                    //_detec->_logText << endl << "filename=" <<  wavFile << "   cri " << icri << endl;
                    //_detec->_logText << "jcrete=" << jcrete << " serie SDC" << endl;
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
                            totdif2Op = 0.0f;
                            totdif2opwb = 0.0f;
                            xdOp=x;
                        }
                        if((jcrete == (CS-1) || jcrete == (CM-1)) && (dif1>1 || x==xmax) && findSdcl==false)
                        {
                            prem = x;
                            oParamCrete[jcrete][SDCL]  = totdif2;
                            oParamCrete[jcrete][SDCLR]  = ((float)totdif2) /  ((float)(x-xmin)+0.5f);
                            oParamCrete[jcrete][SDCLWB]  = totdif2wb;
                            oParamCrete[jcrete][SDCLRWB]  = (float)totdif2wb/((float)(x-xmin)+0.5f);
                            if(qAbs(pc[x-xmin]-pc[0])>0)
                            {
                                float ratioY = (totdif2*10.0f)/(float)(qAbs(pc[x-xmin]-pc[0]));
                                oParamCrete[jcrete][SDCLRY]  = ratioY;
                                oParamCrete[jcrete][SDCLRXY]  = ratioY / ((float)(x-xmin)+0.5f);
                                oParamCrete[jcrete][SDCLRXY2]  = ratioY / ((float)(x-xmin)+0.5f);
                                //
                                float ratioYwb = (totdif2wb*10.0f)/(float)(qAbs(pc[x-xmin]-pc[0]));
                                oParamCrete[jcrete][SDCLRYWB]  = ratioYwb;
                                oParamCrete[jcrete][SDCLRXYWB]  = ratioYwb / ((float)(x-xmin)+0.5f);
                            }
                            findSdcl = true;
                        }
                        if((jcrete == (CS-1) || jcrete == (CM-1)) && (dif1>1 || x==xmax) && findSdclOp==false && seenMasterPoint)
                        {
                            oParamCrete[jcrete][SDCLOP]  = totdif2Op;
                            oParamCrete[jcrete][SDCLROP]  = (float)totdif2Op/((float)(x-xdOp)+0.5f);
                            oParamCrete[jcrete][SDCLOPWB]  = totdif2opwb;
                            oParamCrete[jcrete][SDCLROPWB]  = (float)totdif2opwb/((float)(x-xdOp)+0.5f);
                            if(qAbs(pc[x-xmin]-pc[xdOp-xmin])>0)
                            {
                                float ratioYOp = (totdif2Op*100.0f)/(float)(qAbs(pc[x-xmin]-pc[xdOp-xmin]));
                                oParamCrete[jcrete][SDCLRYOP]  = ratioYOp;
                                oParamCrete[jcrete][SDCLRXYOP]  = ratioYOp / ((float)(x-xdOp)+0.5f);
                                float ratioYopwb = (totdif2opwb*10.0f)/(float)(qAbs(pc[x-xmin]-pc[xdOp-xmin]));
                                oParamCrete[jcrete][SDCLRYOPWB]  = ratioYopwb;
                                oParamCrete[jcrete][SDCLRXYOPWB]  = ratioYopwb / ((float)(x-xdOp)+0.5f);
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
                    oParamCrete[jcrete][SDC]  = totdif2;
                    oParamCrete[jcrete][SDCR]  = (float)totdif2/((float)(xmax-xmin)+0.5f);
                    float labs = qAbs(pc[xmax-xmin]-pc[0])+1;
                    float ratioY = (float)totdif2*10.0f/labs;
                    oParamCrete[jcrete][SDCRY]  = ratioY;
                    oParamCrete[jcrete][SDCRXY]  = ratioY / ((float)(xmax-xmin)+0.5f);

                    float fncs = (float) ncs/100.0f;
                    oParamCrete[jcrete][SDCL_DNP]  = oParamCrete[jcrete][SDCL]/fncs;
                    oParamCrete[jcrete][SDCLR_DNP]  = oParamCrete[jcrete][SDCLR]/fncs;
                    oParamCrete[jcrete][SDCLRY_DNP]  = oParamCrete[jcrete][SDCLRY]/fncs;
                    oParamCrete[jcrete][SDCLRXY_DNP]  = oParamCrete[jcrete][SDCLRXY]/fncs;
                    oParamCrete[jcrete][SDCLRXY2_DNP]  = oParamCrete[jcrete][SDCLRXY2]/fncs;
                    //
                    oParamCrete[jcrete][SDCLOP_DNP]  = oParamCrete[jcrete][SDCLOP]/fncs;
                    oParamCrete[jcrete][SDCLROP_DNP]  = oParamCrete[jcrete][SDCLROP]/fncs;
                    oParamCrete[jcrete][SDCLRYOP_DNP]  = oParamCrete[jcrete][SDCLRYOP]/fncs;
                    oParamCrete[jcrete][SDCLRXYOP_DNP]  = oParamCrete[jcrete][SDCLRXYOP]/fncs;
                    // -----------------------------------------------------------------------------------------
                    // Paramètres de la série Coudes
                    // _detec->_logText << "jcrete=" << jcrete << " série coudes" << endl;
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
                            //_detec->_logText << "arret en prem = " << x-xmin << "  soit " << x*_msPerX<<"ms et "<<pc[x-xmin]*_khzPerY<< "kHz"<< endl;
                            bool unecond = false;
                            if(y2>y1) unecond = true;
                            else
                            {
                                if(y1!=y0)
                                {
                                    if ((((float)(y2-y1)) / ((float)(y1-y0))) < 0.6) unecond = true;
                                }
                            }
                            if(unecond)
                            {
                                if(prem>xmin)
                                {
                                    oParamCrete[jcrete][ELBPOS]  = ((float)(x-xmin))/((float)(prem-xmin));
                                }
                                if(x>xmin)
                                {
                                    oParamCrete[jcrete][ELBSB]  = ((float)(y1-pc[0]))/((float)(x-xmin));
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
                                if(y1>yd) if  ( (((float)(y2-y1)) / ((float)(y1-yd))) < 0.6)
                                {
                                    unecond = true;
                                }
                            }
                            if(unecond)
                            {
                                if(prem>xmin)
                                {
                                    oParamCrete[jcrete][ELB2POS]  = ((float)(x-xmin))/((float)(prem-xmin));
                                }
                                if(x>xmin)
                                {
                                    oParamCrete[jcrete][ELB2SB]  = ((float)(y1-pc[0]))/((float)(x-xmin));
                                    //_detec->_logText << "elb2sb = " << ((float)(y1-pc[0]))/((float)(x-xmin)) << endl;

                                }
                                break;
                            }
                        }
                    }
                   // _detec->_logText << "jcrete=" << jcrete << "fin série coudes" << endl << endl;
                    // -----------------------------------------------------------------------------------------
                    // nouveau paramètres (à réintégrer dans boucles du début de la méthode)
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
                        float ratioF = 0.0f;
                        float ratioE = 0.0f;
                        float ratioFP = 0.0f;
                        float ratioFP2 = 0.0f;
                        float ratioFP3 = 0.0f;
                        float totF[2],totE[2],np[2],totFP[2],totFP2[2],np2[2],totFP3[2],np3[2];
                        float totE2[2],totE4[2],totE8[2],totE16[2];
                        int hx2[2],hx4[2],hx8[2],hx16[2];
                        int nra2[2],nra4[2],nra8[2],nra16[2];
                        float fc2,fc4,fc8,fc16;
                        int avap2,avap4,avap8,avap16;
                        for(int j=0;j<2;j++)
                        {
                            totF[j]=0.0f; totE[j]=0.0f; np[j]=0.0f; totFP[j]=0.0f;
                            np2[j]=0.0f; totFP2[j]=0.0f;  np3[j]=0.0f; totFP3[j]=0.0f;
                            totE2[j]=0.0f; totE4[j]=0.0f; totE8[j]=0.0f; totE16[j]=0.0f;
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
                            float e = _sonogramArray[pc[x-xmin]][x];
                            if(pc[x-xmin]>0)
                            {
                                bool condRestart = false;
                                if(alreadyRestarted==false)
                                {
                                    if((x-xmin)<(xmaitr-xmin)/4)
                                    {
                                        if(e<(_highThreshold - _lowThreshold)) condRestart = true;
                                        if(x>xmin+2 && pc[x-xmin-1]==0 && pc[x-xmin-2]==0) condRestart = true;
                                    }
                                }
                                if(condRestart)
                                {
                                    alreadyRestarted = true;
                                    for(int j=0;j<2;j++)
                                    {
                                        totF[j]=0.0f; totE[j]=0.0f; np[j]=0.0f; totFP[j]=0.0f;
                                        np2[j]=0.0f; totFP2[j]=0.0f;  np3[j]=0.0f; totFP3[j]=0.0f;
                                        totF[j]=0.0f; totE[j]=0.0f; np[j]=0.0f; totFP[j]=0.0f;
                                        np2[j]=0.0f; totFP2[j]=0.0f;  np3[j]=0.0f; totFP3[j]=0.0f;
                                        totE2[j]=0.0f; totE4[j]=0.0f; totE8[j]=0.0f; totE16[j]=0.0f;
                                        hx2[j]=0; hx4[j]=0; hx8[j]=0; hx16[j]=0;
                                        nra2[j]=0; nra4[j]=0; nra8[j]=0; nra16[j]=0;
                                    }
                                    xdep = x;
                                }
                                //
                                if(x<=xmaitr) avap=0; else avap=1;
                                totF[avap]+=(float)pc[x-xmin];
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
                                    fc2=1.0f; fc4=1.0f; fc8=1.0f; fc16=1.0f;
                                    if(avap16==0)
                                    {
                                        if(nra16[0]<3) fc16=4.0f-nra16[0];
                                    }
                                    if(avap8==0)
                                    {
                                        if(nra8[0]<4) fc8=5.0f-nra8[0];
                                    }
                                    if(avap4==0)
                                    {
                                        if(nra4[0]<6) fc4=7.0f-nra4[0];
                                    }
                                    if(avap2==0)
                                    {
                                        if(nra2[0]<9) fc2=10.0f-nra2[0];
                                    }

                                    nra2[avap2]+=fc2; nra4[avap4]+=fc4; nra8[avap8]+=fc8; nra16[avap16]+=fc16;
                                    float ecol = _tabX[icri][x-xmin];
                                    totE2[avap2]+=ecol*fc2; totE4[avap4]+=ecol*fc4; totE8[avap8]+=ecol*fc8; totE16[avap16]+=ecol*fc16;
                                    float hpix = _numberPixelsPerX[x-xmin];
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
                        } // next x
                        if(np[0]>0 && np[1]>0)
                        {
                            ratioF = (totF[0]*100.0f/np[0]) / (totF[1]/np[1]);
                            ratioE = (totE[0]*100.0f/np[0]) / (totE[1]/np[1]);
                            ratioFP = (totFP[0]*100.0f/np[0]) / (totFP[1]/np[1]);
                            ratioFP2 = (totFP2[0]*10.0f/np2[0]) / (totFP2[1]/np2[1]);
                            ratioFP3 = (totFP3[0]*100.0f/np3[0]) / (totFP3[1]/np3[1]);
                            oParamCrete[jcrete][RAF]  = ratioF;
                            oParamCrete[jcrete][RAE]  = ratioE;
                            oParamCrete[jcrete][RAFE]  = ratioF+ratioE;
                            oParamCrete[jcrete][RAFP]  = ratioFP;
                            oParamCrete[jcrete][RAFP2]  = ratioFP2;
                            oParamCrete[jcrete][RAFP3]  = ratioFP3;
                            oParamCrete[jcrete][SBMP]  = ( (float)  (pc[xmaitr-xmin] - pc[xdep-xmin] ))  /  ((float)(xmaitr-xdep));
                            oParamCrete[jcrete][SAMP]  = ( (float)  (pc[xfin-xmin] - pc[xmaitr-xmin]     ))   /  ((float)(xfin-xmaitr));
                            oParamCrete[jcrete][SBAR]  = oParamCrete[jcrete][SAMP]  - oParamCrete[jcrete][SBMP];
                            //
                            if(nra2[0]>0 && nra2[1]>0 )
                            {
                                if(hx2[1]>0)
                                    oParamCrete[jcrete][RAHP2]  = (((float)hx2[0]*10.0f)/((float)nra2[0]))  / (((float)hx2[1])/((float)nra2[1]));
                                //
                                if(totE2[1] >0)
                                    oParamCrete[jcrete][RAHE2]  = ((totE2[0]*10.0f)/((float)nra2[0])) / ((totE2[1])/((float)nra2[1]));
                            }
                            if(nra4[0]>0 && nra4[1]>0)
                            {
                                if(hx4[1]>0)
                                    oParamCrete[jcrete][RAHP4]     = (((float)hx4[0]*10.0f)/((float)nra4[0])) / (((float)hx4[1])/((float)nra4[1]));
                                if(totE4[1] >0)
                                    oParamCrete[jcrete][RAHE4]  = ((totE4[0]*10.0f)/((float)nra4[0])) / ((totE4[1])/((float)nra4[1]));

                            }
                            if(nra8[0]>0 && nra8[1]>0 && hx8[1]>0)
                            {
                                if(hx8[1]>0)
                                    oParamCrete[jcrete][RAHP8]  = (((float)hx8[0]*10.0f)/((float)nra8[0])) / (((float)hx8[1])/((float)nra8[1]));
                                if(totE8[1] >0)
                                    oParamCrete[jcrete][RAHE8]  = ((totE8[0]*10.0f)/((float)nra8[0])) / (totE8[1]/((float)nra8[1]));
                            }
                            if(nra16[0]>0 && nra16[1]>0 && hx16[1]>0)
                            {
                                if(hx16[1]>0)
                                    oParamCrete[jcrete][RAHP16]  = (((float)hx16[0]*10.0f)/((float)nra16[0])) / (((float)hx16[1])/((float)nra16[1]));
                                if(totE16[1] >0)
                                    oParamCrete[jcrete][RAHE16]  = ((totE16[0]*10.0f)/((float)nra16[0])) / (totE16[1]/((float)nra16[1]));
                            }
                            //
                        }
                    }
                }
                //
                // -----------------------------------------------------------------------------------------
                int milieu = (xmin+xmax)/2;
                for(int k=milieu;k<=xmax;k++)
                    if(pc[k-xmin]>0)
                    {
                        milieu=k;
                        break;
                    }

                oParamCrete[jcrete][CeF] = (float)pc[milieu-xmin]*_khzPerY;
                if(jcrete==0)
                {
                    int pos1=xmcmax[jcrete];
                    int pos2=xmax;
                    //£ float borne = emcmax - 5;
                    int borne = ((int)emcmax - 5)*100;
                    if(xmin<xmcmax[jcrete])
                    {
                        for(int k=xmin;k<xmcmax[jcrete];k++)
                        {
                            if(pc[k-xmin]>0)
                            {
                                if(_sonogramArray[pc[k-xmin]][k]>borne)
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
                                if(_sonogramArray[pc[k-xmin]][k]<borne)
                                {
                                    pos2=k;
                                    break;
                                }
                            }
                        }
                    }
                    oParamCrete[jcrete][B5dBBF]  = pc[pos1-xmin]*_khzPerY;
                    oParamCrete[jcrete][B5dBAF]  = pc[pos2-xmin]*_khzPerY;
                    oParamCrete[jcrete][B5dBBW]  = oParamCrete[jcrete][B5dBAF]-oParamCrete[jcrete][B5dBBF];
                    oParamCrete[jcrete][B5dBDur] = (pos2-pos1)*_msPerX;
                }
            }
            else
            {
                meilleur = ymin; meilleurepente = 100000.0f;
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
                                _slope[y-ymin] = 4.0f/(float)dx;
                        }
                        if(_slope[y-ymin] < meilleurepente)
                        {
                            meilleur = y;
                            meilleurepente = _slope[y-ymin];
                        }
                    }
                    //
                    //£ float e=_sonogramArray[y][pc[y-ymin]];
                    float e = (float)_sonogramArray[y][pc[y-ymin]]/100.0f;
                    if(e>emcmax)
                    {
                        emcmax = e;
                        xmcmax[jcrete]=pc[y-ymin];
                        ymcmax[jcrete] = y;
                    }

                }
                //
                //
                if(pc[ymax-ymin]!=pc[0])
                {
                    oParamCrete[jcrete][Slope] = ((float)(ymax-ymin)*_khzPerY)/((float)(pc[ymax-ymin]-pc[0])*_msPerX);
                }
                else oParamCrete[jcrete][Slope] = 9999.0f;
                if(ymin != ymax)
                {
                    oParamCrete[jcrete][ISlope] = ((float)(pc[ymax-ymin]-pc[0])*_msPerX) /  ((float)(ymax-ymin)*_khzPerY);
                }
                else oParamCrete[jcrete][ISlope] = 9999.0f;
                lowfc = meilleur;
                _lowSlope[(icri*NCRETES+jcrete)*2]=pc[lowfc-ymin];
                _lowSlope[(icri*NCRETES+jcrete)*2+1]=lowfc;
                oParamCrete[jcrete][FIF] = lowfc*_khzPerY;
                float pe1,pe2,dif;
                float meilleuredist;
                int meilleurk[2];
                int debb,finb;
                meilleurk[0]=ymin;
                meilleurk[1]=lowfc;
                int sens=1;
                if(pc[ymax-ymin]>=pc[0]) sens=1;
                else sens=-1;
                for(int jp=0;jp<2;jp++)
                {
                    meilleuredist = 1000000.0f;
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
                            dif=0.0f;
                            if(pc[k-ymin]>0)
                            {
                                if(pc[k-ymin]-pc[debb-ymin] != 0)
                                    pe1=(float)(k-debb)/((float)(pc[k-ymin]-pc[debb-ymin]));
                                else pe1=100000;
                                if(pc[finb-ymin]-pc[k-ymin] != 0)
                                    pe2=((float)(finb-k))/((float)(pc[finb-ymin]-pc[k-ymin]));
                                else pe2=100000;
                            }
                            for(int s = debb+1 ; s<k ; s++)
                                if(pc[s-ymin]>0)
                                {
                                    if(pe1>99999)
                                        dif +=(float)pow((float)(pc[s-ymin]-pc[debb-ymin]),2);
                                    else
                                        dif += (float)pow((float)(pc[s-ymin]-pc[debb-ymin])-(float)(s-debb)/pe1,2);
                                }

                            for(int s = k+1 ; s<=finb-1 ; s++)
                                if(pc[s-ymin]>0)
                                {
                                    if(pe2>99999)
                                        dif +=(float)pow((float)(pc[s-ymin]-pc[finb-ymin]),2);
                                    else
                                        dif += (float)pow((float)(pc[s-ymin]-pc[k-ymin])-(float)(s-k)/pe1,2);
                                }



                            if(dif<meilleuredist)
                            {
                                meilleurk[jp]=k;
                                meilleuredist=dif;
                            }
                        }

                    }
                }
                if(meilleurk[0]==0) hifc=ymin; else hifc=meilleurk[0];
                _inflexion1[(icri*NCRETES+jcrete)*2]=pc[hifc-ymin];
                _inflexion1[(icri*NCRETES+jcrete)*2+1]=hifc;
                oParamCrete[jcrete][HCF] = hifc*_khzPerY;
                oParamCrete[jcrete][THCF]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][THCF] = ((float)(pc[meilleurk[0]-ymin]-xmin))/((float)(xmax-xmin));
                fc3=meilleurk[1];
                _inflexion3[(icri*NCRETES+jcrete)*2]=pc[fc3-ymin];
                _inflexion3[(icri*NCRETES+jcrete)*2+1]=fc3;
                oParamCrete[jcrete][LCF] = fc3*_khzPerY;
                oParamCrete[jcrete][UpSl] = 0.0f;
                if((pc[meilleurk[0]-ymin]-pc[0]!=0 && sens==1)
                        || (pc[meilleurk[0]-ymin]-pc[ymax-ymin]!=0 && sens==-1))
                {
                    if(sens==1)
                        oParamCrete[jcrete][UpSl] = ((float)(meilleurk[0]-ymin)*_khzPerY)
                                /( ( (float)(pc[meilleurk[0]-ymin]-pc[0]) ) * _msPerX);
                    else
                        oParamCrete[jcrete][UpSl] = ((float)(meilleurk[0]-ymax)*_khzPerY)
                                /( ( (float)(pc[meilleurk[0]-ymin]-pc[ymax-ymin]) ) * _msPerX);

                }
                oParamCrete[jcrete][LoSl] = 9999.0f;
                if(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]!=0)
                    oParamCrete[jcrete][LoSl] = ((float)(lowfc-meilleurk[0])*_khzPerY)
                            /( ( (float)(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]) ) * _msPerX);

                //
                int debp,finp,ecart;
                float pj[4];
                for(int j2=0;j2<4;j2++)
                {
                    pj[j2]=9999.0f;
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
                    float dx=(float)(pc[finp]-pc[debp]);
                    if(dx!=0)
                    {
                        pj[j2]=((float)ecart*_khzPerY)/(dx*_msPerX);
                    }
                }
                oParamCrete[jcrete][StSl] = pj[0];
                oParamCrete[jcrete][EnSl]   = pj[1];
                oParamCrete[jcrete][FISl]  = pj[3];
                oParamCrete[jcrete][FPSl]  = pj[2];

                if(_withNewParams)
                {
                    // calcul du nouveau paramètre SDC
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
                    oParamCrete[jcrete][SDC]  = totdif2;
                    float labs = qAbs(pc[ymax-ymin]-pc[0])+1;
                    oParamCrete[jcrete][SDCR]  = (float)totdif2/labs;
                }
                //
                // calcul des paramètres ...5db pour crête ouest 2
                if(jcrete==4)
                {
                    int pos1=ymcmax[jcrete];
                    int pos2=ymax;
                    //£ float borne = emcmax - 5;
                    int borne = ((int)emcmax - 5)*100;
                    if(ymin<ymcmax[jcrete])
                    {
                        for(int k=ymin;k<ymcmax[jcrete];k++)
                        {
                            if(pc[k-ymin]>0)
                            {
                                if(_sonogramArray[k][pc[k-ymin]]>borne)
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
                                if(_sonogramArray[k][pc[k-ymin]]<borne)
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
                    oParamCrete[jcrete][B5dBBF]  = pos1*_khzPerY;
                    oParamCrete[jcrete][B5dBAF]  = pos2*_khzPerY;
                    oParamCrete[jcrete][B5dBBW]  = oParamCrete[jcrete][B5dBAF]-oParamCrete[jcrete][B5dBBF];
                    oParamCrete[jcrete][B5dBDur] = (pc[pos2-ymin]-pc[pos1-ymin])*_msPerX;
                } // fin if(jcrete==4)
                //

            } // fin else jcrete....
            oParamCrete[jcrete][FPk]=(float)ymcmax[jcrete]*(float)_khzPerY;
            oParamCrete[jcrete][FPkD] = 0.0f;
            if(criprec<icri) oParamCrete[jcrete][FPkD] = oParamCrete[jcrete][FPk] - _paramsArray[criprec][jcrete+1][FPk];
            int intervTemps;
            if(jcrete<3) intervTemps = xmax-xmin+1;
            else intervTemps = tempsOuest[jcrete-3];
            oParamCrete[jcrete][TPk] = (float)(xmcmax[jcrete]-xmin+0.5f)/(float)intervTemps;
        } // fin boucle jcrete
        oParamCrete[0][StF] = (float)_yEmaxPerX[icri][0]*_khzPerY;
        oParamCrete[0][EnF] = (float)_yEmaxPerX[icri][xmax-xmin]*_khzPerY;
        if(_imageData)
        {
            _callMasterRidgeArray.push_back(_callMasterRidge);
            _callSouthArray.push_back(_callSouthRidge);
            _callNorthRidgeArray.push_back(_callNorthRidge);
            _callWestRidgeArray.push_back(_callWestRidge);
            _callSecondWestRidgeArray.push_back(_callSecondWestRidge);
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
    //_detec->_logText << endl << "-9" ;
    for(int icri=0;icri<nbcris-1;icri++)
    {
        float freqmp = _paramsArray[icri][SH][FreqMP];
        int pmp = invMp[icri];
        if(pmp<nbcris-1)
        {
            _paramsArray[icri][SH][NextMP1] = (_masterPoints.at(sortMp[pmp+1]).x()-_masterPoints.at(icri).x())*_msPerX;
            for(int k=pmp+1;k<nbcris;k++)
            {
                if(_paramsArray[sortMp[k]][SH][Fmin] < freqmp
                        && _paramsArray[sortMp[k]][SH][Fmax] > freqmp)
                {
                    _paramsArray[icri][SH][NextMP2] = (_masterPoints.at(sortMp[k]).x() -_masterPoints.at(icri).x())*_msPerX;
                    break;
                }
            }
        }

    }
    for(int k=0;k<2;k++)
        for(int i=0;i<nbcris;i++)
        {_harmonic[k][i] = -1; _dpm[k][i]=0.0f; _dypm[k][i]=0; _ypm[k][i]=0;}

    for(int i=0;i<nbcris;i++)
    {
        int xmin1=_vectorXMin.at(i);
        int xmax1=_vectorXMax.at(i);
        int xmaitr1=_masterPoints.at(i).x();
        int ymaitr1=_masterPoints.at(i).y();
        for(int k=0;k<2;k++)
        {
            for(int j=0;j<nbcris;j++)
            {
                if(j==i) continue;
                int xmax2=_vectorXMax.at(j);
                if(xmax2<xmin1) continue;
                int ymaitr2=_masterPoints.at(j).y();
                if((k==0 && ymaitr2<ymaitr1) || (k==1 && ymaitr2>ymaitr1)) continue;
                int xmin2=_vectorXMin.at(j);
                if(xmaitr1>=xmin2 && xmaitr1<=xmax2)
                {
                    int xmaitr2=_masterPoints.at(j).x();
                    if(xmaitr2>=xmin1 && xmaitr2<=xmax1)
                    {
                        float dm =pow(xmaitr2-xmaitr1,2)+pow(ymaitr2-ymaitr1,2);
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
                if(xmin2>xmax1) break; // puisque les cris sont classés par xmin
            }
        }
    }

    //_detec->_logText << endl << "début traitement nhsup-nhinf" << endl;
    int nhsup,nhinf;
    for(int icri=0;icri<nbcris;icri++)
    {
         //_detec->_logText <<  "icri=" << icri << endl;
        nhsup=_harmonic[0][icri];
        if(nhsup>=0)
        {
             //_detec->_logText <<  "nhsup=" << nhsup << endl;
            if(_paramsArray[icri][SH][FreqMP] > 0.0f)
                _paramsArray[icri][SH][Hup_RFMP] = _paramsArray[nhsup][SH][FreqMP]
                        /_paramsArray[icri][SH][FreqMP];
            int xmin1=_vectorXMin.at(icri);
            int xmin2=_vectorXMin.at(nhsup);
            int xmax1=_vectorXMax.at(icri);
            int xmax2=_vectorXMax.at(nhsup);
            float larcri=(float)(xmax1-xmin1);
            if(larcri>0.0f)
            {
                _paramsArray[icri][SH][Hup_PosMP]=(float)(_masterPoints.at(nhsup).x()-xmin1)/larcri;
                _paramsArray[icri][SH][Hup_PosSt]=(float)(xmin2-xmin1)/larcri;
                _paramsArray[icri][SH][Hup_PosEn]=(float)(_vectorXMax.at(nhsup)-xmin1)/larcri;
            }
            int xdeb = qMax(xmin1,xmin2);
            int xfin = qMin(xmax1,xmax2);
            float famp1=0.0f;
            float famp2=0.0f;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=_tabX[icri][i-xmin1];
                famp2+=_tabX[nhsup][i-xmin2];
            }
            // if(famp1>0.0f) _paramsArray[icri][SH][Hup_AmpDif] = famp2/famp1;
            _paramsArray[icri][SH][Hup_AmpDif] = (famp2 - famp1)*_khzPerY;

           // _detec->_logText <<  "  zz1" << endl;
            int nouveauxdeb=xdeb;
            int nouveauxfin=xfin;
            for(int i=xdeb;i<=xfin;i++)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhsup][i-xmin2]<=0) nouveauxdeb++;
                else break;

            }
           // _detec->_logText <<  "  zz1,5" << endl;
            for(int i=xfin;i>=nouveauxdeb;i--)
            {
                if(_yEmaxPerX[icri][i-xmin1]<=0 || _yEmaxPerX[nhsup][i-xmin2]<=0) nouveauxfin--;
                else break;
            }
            xdeb=nouveauxdeb; xfin=nouveauxfin;
           //  _detec->_logText <<  "  zz2" << endl;

            if(xfin>xdeb)
            {

             //    _detec->_logText <<  "  zz3" << endl;
                float CM_Slope1=((float)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((float)(xfin-xdeb));
                float CM_Slope2=((float)(_yEmaxPerX[nhsup][xfin-xmin2]
                                 -_yEmaxPerX[nhsup][xdeb-xmin2]))
                        /((float)(xfin-xdeb));
                if(CM_Slope1!=0.0f)
                    _paramsArray[icri][SH][Hup_RSlope] = CM_Slope2 / CM_Slope1;
               //  _detec->_logText <<  "  zz4" << endl;
            }
        }
        nhinf=_harmonic[1][icri];
        if(nhinf>=0)
        {
//             _detec->_logText <<  "nhinf=" << nhinf << endl;
            if(_paramsArray[icri][SH][FreqMP] > 0.0f)
                _paramsArray[icri][SH][Hlo_RFMP] = _paramsArray[nhinf][SH][FreqMP]
                        /_paramsArray[icri][SH][FreqMP];
            int xmin1=_vectorXMin.at(icri);
            int xmax1=_vectorXMax.at(icri);
            int xmin2=_vectorXMin.at(nhinf);
            int xmax2=_vectorXMax.at(nhinf);
            float larcri=(float)(xmax1-xmin1);
            if(larcri>0.0f)
            {
                _paramsArray[icri][SH][Hlo_PosMP]=(float)(_masterPoints.at(nhinf).x()-xmin1)/larcri;
                _paramsArray[icri][SH][Hlo_PosSt]=(float)(xmin2-xmin1)/larcri;
                _paramsArray[icri][SH][Hlo_PosEn]=(float)(xmax2-xmin1)/larcri;
            }
            int xdeb = qMax(xmin1,xmin2);
            int xfin = qMin(xmax1,xmax2);
            float famp1=0.0f;
            float famp2=0.0f;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=_tabX[icri][i-xmin1];
                famp2+=_tabX[nhinf][i-xmin2];
            }
            // if(famp1>0.0f) _paramsArray[icri][SH][Hlo_AmpDif] = famp2/famp1;
            _paramsArray[icri][SH][Hlo_AmpDif] = (famp2-famp1)*_khzPerY;
  //           _detec->_logText <<  "  zz11" << endl;
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
    //         _detec->_logText <<  "  zz12" << endl;

            if(xfin>xdeb)
            {
      //           _detec->_logText <<  "  zz13" << endl;
                float CM_Slope1=((float)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((float)(xfin-xdeb));
                float CM_Slope2=((float)(_yEmaxPerX[nhinf][xfin-xmin2]
                                 -_yEmaxPerX[nhinf][xdeb-xmin2]))
                        /((float)(xfin-xdeb));
                if(CM_Slope1!=0.0f)
                    _paramsArray[icri][SH][Hlo_RSlope] = CM_Slope2 / CM_Slope1;
        //         _detec->_logText <<  "  zz14" << endl;
            }
        }
    }
    //
    // _detec->_logText << "fin traitement nhsup nhinf."  << endl;
    int nbb;
    float tsono = _sonogramWidth *  _msPerX;
    //_detec->_logText << "tsono=" << tsono << endl;
    float interv[MAXCRI];
    float variation[MAXCRI];
    float proxiFreq = 2.0f;
    int rband[MAXCRI]; //+
    // ZZZZZZZZ
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
            //
            _paramsArray[icri][SH][VarInt]  = 0.0f;
            _paramsArray[icri][SH][VarSmInt] = 0.0f;
            _paramsArray[icri][SH][VarLgInt] = 0.0f;
            // 1) alimentation de rband
            nbb=0;
            float freqmp = _paramsArray[icri][SH][FreqMP];
            for(int k=0;k<nbcris;k++)
            {
                float freqmp2 = _paramsArray[k][SH][FreqMP];
                if(qAbs(freqmp2-freqmp) < proxiFreq) rband[nbb++]=k;
            }
            // 2) On retrie la suite de cris de la bande dans l'ordre des PM
            if(nbb>1)
            {
                int rbandx[MAXCRI];
                for(int j=0;j<nbb;j++)  rbandx[j] = _masterPoints.at(rband[j]).x();
                sortIntArrays(rband,nbb,rbandx);
            }
            // 3) tri des intervalles de distance
            float medianDistance = tsono/2;
            float medianLittleDistance = tsono/2;
            float medianBigDistance = tsono/2;
            //int nbi = nbb-1;
            // 18-02-2015 :
            int nbi = nbb+1;
            int halfnbi = nbi/2;
            int quartnbi =halfnbi/2;
            if(nbi>1)
            {
                interv[0] = _masterPoints.at(rband[0]).x() * _msPerX;
                interv[nbi-1] = (_sonogramWidth - _masterPoints.at(rband[nbb-1]).x()) * _msPerX;
                for(int j=1;j<nbi-1;j++)
                {
                    interv[j] = (_masterPoints.at(rband[j]).x()-_masterPoints.at(rband[j-1]).x()) * _msPerX;
                }
                sortFloatArray(interv,nbi);
                // 4) calcul de la moyenne des petits et des grands intervalles
                // et des paramètres liés
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
            //_detec->_logText << "medianLittleDistance = " << medianLittleDistance << endl;
            //_detec->_logText << "medianBigDistance = " << medianBigDistance << endl;
            _paramsArray[icri][SH][MedInt]   = medianDistance;
            _paramsArray[icri][SH][Int25]  = medianLittleDistance;
            _paramsArray[icri][SH][Int75]  = medianBigDistance;
            if(medianLittleDistance>0.0f) _paramsArray[icri][SH][RInt1] = medianBigDistance / medianLittleDistance;
            else _paramsArray[icri][SH][RInt1] = 9999;
            //
            // 5) Variations des intervalles
            float medianDistanceVariation = 0;
            // a) écart médian
            for(int j=0;j<nbi;j++) variation[j] = qAbs(interv[j]-medianDistance);
            float medianLittleDistanceVariation = variation[0];
            float medianBigDistanceVariation = variation[nbi-1];
            if(nbi>1)
            {
                sortFloatArray(variation,nbi);
                if(halfnbi*2<nbi) medianDistanceVariation = variation[halfnbi];
                else medianDistanceVariation = (variation[halfnbi-1] + variation[halfnbi])/2;
                // b) écarts sur petits intervalles et sur grands intervalles
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
                    sortFloatArray(variation,halfnbi);
                    sortFloatArray(variation+nbi-halfnbi,halfnbi);
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
                } // fin nbi>=4
            } // fin if nbi>1
            _paramsArray[icri][SH][IntDev]   = medianDistanceVariation;
            _paramsArray[icri][SH][SmIntDev]  = medianLittleDistanceVariation;
            _paramsArray[icri][SH][LgIntDev]  = medianBigDistanceVariation;
            //
            _paramsArray[icri][SH][VarInt]  = medianDistanceVariation/medianDistance;
            if(medianLittleDistance>0.0f) _paramsArray[icri][SH][VarSmInt]  = medianLittleDistanceVariation/medianLittleDistance;
            _paramsArray[icri][SH][VarLgInt]  = medianBigDistanceVariation/medianBigDistance;
            //
            if(medianLittleDistanceVariation==0) _paramsArray[icri][SH][RIntDev1] = 0.0f;
            else _paramsArray[icri][SH][RIntDev1] = medianBigDistanceVariation / medianLittleDistanceVariation;
            // affecter les paramètres liés à la bande rband
            // sur le même axe y
            int ymaitr = _masterPoints.at(icri).y();
            for(int jb =0;jb<nbb;jb++)
            {
                int jcri = rband[jb];
                if(_masterPoints.at(jcri).y()==ymaitr && jcri != icri)
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
        } // fin if alreadytreated = false
    } // next icri
    //
    //_logText << "-dp2 fin" << endl;
    delete[] sortMp;
    delete[] invMp;
    delete[] xMp;

}

void DetecTreatment::sortFloatArray(float *pf,int nbf)
{
    float conserv;
    bool ontrie = true;
    while(ontrie)
    {
        ontrie = false;
        for(int j=0;j<nbf-1;j++)
        {
            bool permuter = false;
            if(pf[j]>pf[j+1]) permuter = true;
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

void DetecTreatment::sortFloatIndArray(float *pf,int nbf,int *pos)
{
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
            if(pf[j]>pf[j+1]) permuter = true;
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

void DetecTreatment::saveParameters(const QString& wavFile)
{
    QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ ResultSuffix;
    QFile txtFile;
    txtFile.setFileName(txtFilePath);
    if(txtFile.open(QIODevice::WriteOnly | QIODevice::Text)==false)
    {
        _detec->_logText  << "ouverture en écriture du fichier " << txtFilePath << " impossible !" << endl;
        return;
    }
    QTextStream fileStream;
    fileStream.setDevice(&txtFile);
    fileStream.setRealNumberNotation(QTextStream::FixedNotation);
    fileStream.setRealNumberPrecision(2);
    fileStream << "Filename"<< '\t' << "CallNum"
               << '\t' << "Version"<< '\t' << "FileDur"<< '\t' << "SampleRate";
    for(int j=0;j<_numberCallParameters;j++) fileStream << '\t' << _vectPar[j].ColumnTitle;
    fileStream << endl;
    //float **parArray;
    //float u_f;
    float dur = (float)_sonogramWidth*_msPerX/1000;
    float sr = (float)_soundFileInfo.samplerate;
    for(int i=0;i<_callsNumber;i++)
    {
        float **callArray = _paramsArray[i];
        fileStream << wavFile << '\t' << QString::number(i)
                    << '\t' << QString::number(_paramVersion) << '\t'  <<  dur << '\t' << sr;
        for(int j=0;j<_numberCallParameters;j++)
        {
            //u_f = callArray[_vectPar[j].NumTableau][_vectPar[j].NumPar];
            //u_f = ((float)qRound(u_f*100.0f))/100.0f;
            //_simpleParamsArray[i][j]=u_f;
            fileStream << '\t' <<  callArray[_vectPar[j].NumTableau][_vectPar[j].NumPar];
        }
        fileStream << endl;
    }
    txtFile.close();
}

void DetecTreatment::saveCompressedParameters(const QString& wavFile)
{
    QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ ResultSuffix;
    QString compressedParametersPath = _txtPath+"/"+wavFile.left(wavFile.length()-3) + ResultCompressedSuffix;
    QString program = "7z";
    QStringList  arguments;
    arguments << "a" << "-tgzip" << compressedParametersPath <<  txtFilePath;
    QProcess::execute(program,arguments);
    //_detec->_logText << "résultat d'appel 7z=" << a << endl;

    /*
    _compressedParametersFile.setFileName(compressedParametersPath);
    if(_compressedParametersFile.open(QIODevice::WriteOnly)==false)
    {
        _detec->_logText  << "Ouverture en écriture du fichier " << compressedParametersPath << " impossible !" << endl;
        return;
    }
    // entête
    _compressedParametersStream.setDevice(&_compressedParametersFile);
    _compressedParametersStream << (qint8)_compressedVersion;
    //_compressedParametersStream << wavFile;
    _compressedParametersStream << (qint16)_callsNumber;
    if(_numberCallParameters>=255)
    {
        _compressedParametersStream << (qint8)255;
        _compressedParametersStream << (qint8)(_numberCallParameters-255);
    }
    else _compressedParametersStream << (qint8)_numberCallParameters;
    // compression des paramètres

    for(int i=0;i<_callsNumber;i++)
        for(int j=0;j<_numberCallParameters;j++)
        if(i==0) _valuesToCompressArray[i][j]=_simpleParamsArray[i][j];
        else _valuesToCompressArray[i][j]=_simpleParamsArray[i][j]-_simpleParamsArray[i-1][j];
    //
    qint8 octet1,octet2,decimal;
    int casrepet,nsepar,nsuite;
    int signe;
    float treated;
    int entier;
    bool avoir = false;
    int c=0;
    for(int i=0;i<_callsNumber;i++) for(int j=0;j<_numberCallParameters;j++)
    {
       // if(i==83 && j>160) avoir = true; else avoir = false;
        //if(i==92 && j>156) avoir = true; else avoir = false;
        // recherche si paramètre répété
        casrepet=0;
        treated=_simpleParamsArray[i][j];
        if(avoir) _detec->_logText << "i,j="<<i<<","<<j<<" valeur traitée = "<<treated << endl;
        if(j>=0)
        {
            if(treated==_simpleParamsArray[i][j-1])
            {
                casrepet=1;
                nsuite=0;
                for(int j2=j+1;j<qMin(_numberCallParameters-1,j+32);j2++)
                    if(treated==_simpleParamsArray[i][j2]) nsuite++;
                    else break;
            }
        }
        else
        {
            if(casrepet==0 && i>0)
            {
                for(int i2=i-1;i2>=qMax(0,i-64);i2--)
                    if(treated==_simpleParamsArray[i2][j])
                    {
                        casrepet=3; nsepar=i-1-i2; break;
                    }
            }
            else
            {
                if(casrepet==0 && j>0)
                {
                    for(int j2=j-1;j2>=qMax(0,j-32);j2--)
                        if(treated==_simpleParamsArray[i][j2])
                        {
                            casrepet=2; nsepar=j-1-j2; break;
                        }
                }
            }
        } // fin détermination de cas de répétition
        // codage de la répétition
        if(casrepet>0)
        {
            if(casrepet==1) octet1=224+nsuite;
            if(casrepet==2) octet1=192+nsepar;
            if(casrepet==3) octet1=128+nsepar;
 if(avoir) _detec->_logText << "i,j="<<i<<","<<j<<" casrepet="<<casrepet
            << " octet1=" << octet1 << endl;
            _compressedParametersStream << octet1;
        }
        // codage du cas non répétition
        else
        {
            treated = _valuesToCompressArray[i][j];
            if(treated>=0) octet1=64;
            else {octet1=0; treated = - treated;}
            entier=(int)treated;
            decimal=qRound((treated-entier)*100.0f);
            if(decimal==100) {decimal=0;entier++;}
            if(entier<32)
            {
                octet1+=32+entier;
                if(avoir) _detec->_logText << "C__ i,j="<<i<<","<<j<<" cas simple < 32"
                         << "octet1=" << octet1 << endl;
                _compressedParametersStream << octet1;
            }
            else
            {
                int acoder=entier-32;
                int reste1=0,reste2=0;
                if(acoder>=8191) {reste1=acoder-8191;acoder=8191;}
                octet1+=acoder>>8;
                octet2 = acoder & 255;
                _compressedParametersStream << octet1;
                _compressedParametersStream << octet2;
                if(avoir)
                {
                    _detec->_logText << "C__ i,j="<<i<<","<<j<<" cas simple > 32"
                             << "octet1=" << octet1
                             << "octet2" << octet2 << endl;
                }
                if(acoder==8191)
                {
                    if(reste1>=65535) {reste2=reste1-65535;reste1=65535;}
                    _compressedParametersStream << (qint16)reste1;
                    if(avoir)
                    {
                        _detec->_logText << "C__ reste1="<<reste1<<endl;
                    }
                    if(reste1==65535)
                    {
                        _compressedParametersStream << (qint32)reste2;
                        if(avoir)
                        {
                            _detec->_logText << "C__ reste2="<<reste2<<endl;
                        }
                    }
                }
            } // fin du cas entier>=32
            // codage de la partie décimale
            if(avoir) _detec->_logText << "C__  decimal = " << decimal << endl;
            _compressedParametersStream << (qint8)decimal;
        } // fin du codage cas non répétition
        if(casrepet==1) j+=nsuite;
        // if(i>0 || j>13) avoir=false;
        if(avoir) _detec->_logText << endl;
    } // fin boucles i,j
    //
    _compressedParametersFile.close();
    */

}

/*
void DetecTreatment::expandParameters(const QString& wavFile)
{
    //QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ ResultSuffix;
    QString txtFilePath = _txtPath+"/"+wavFile.left(wavFile.length()-3)+ "ta2";
    QFile txtFile;
    txtFile.setFileName(txtFilePath);
    if(txtFile.open(QIODevice::WriteOnly | QIODevice::Text)==false)
    {
        _detec->_logText  << "ouverture en écriture du fichier " << _txtFilePath << " impossible !" << endl;
        return;
    }
    QTextStream fileStream;
    fileStream.setDevice(&txtFile);
    fileStream.setRealNumberNotation(QTextStream::FixedNotation);
    fileStream.setRealNumberPrecision(2);
    // décomprresssion
    QString compressedParametersPath = _txtPath+"/"+wavFile.left(wavFile.length()-3) + ResultCompressedSuffix;
    _compressedParametersFile.setFileName(compressedParametersPath);
    if(_compressedParametersFile.open(QIODevice::ReadOnly)==false)
    {
        _detec->_logText  << "Ouverture en lecture du fichier " << compressedParametersPath << " impossible !" << endl;
        return;
    }
    _compressedParametersStream.setDevice(&_compressedParametersFile);
    // lecture entete
    qint8 compressVersion;
    _compressedParametersStream >> compressVersion;
    if(compressVersion !=1)
    {
        _detec->_logText  << "Version de compression non gérée !" << endl;
        return;
    }
    //_compressedParametersStream << wavFile;
    qint16 nbcris;
    _compressedParametersStream >> nbcris;
    if(nbcris <0 || nbcris>MAXCRI)
    {
        _detec->_logText  << "Nombre de cris incorrect !" << endl;
        return;
    }
    int nbparam;
    qint8 nbp1,nbp2;
    _compressedParametersStream >> nbp1;
    nbparam=(unsigned char)nbp1;
    if(nbparam==255)
    {
        _compressedParametersStream >> nbp2;
        nbparam+=(unsigned char)nbp2;
    }
    if(nbparam!=_numberCallParameters)
    {
        _detec->_logText  << "Nombre de parametres incorrect : " << nbparam << endl;
        return;
    }
    // décompression des paramètres
    for(int i=0;i<nbcris;i++)
        for(int j=0;j<nbparam;j++)
            _simpleParamsArray[i][j]=0.0f;
    qint8 octet1,octet2,decimal,valRepet;
    qint16 sisuppl;
    qint32 lisuppl;
    int casrepet;
    int entier,signe;
    float lastValue;
    bool avoir=false;
    for(int i=0;i<_callsNumber;i++)
    {
        for(int j=0;j<_numberCallParameters;j++)
        {
            //if(i==83 && j>160) avoir = true; else avoir = false;
            //if(i==0) avoir = true; else avoir = false;
            //if(i==92 && j>156) avoir = true; else avoir = false;
            _compressedParametersStream >> octet1;
            if(avoir) _detec->_logText << "D__ i,j="<<i<<","<<j
                     << " octet1=" << octet1 << endl;
            // cas repetition
            if(((octet1>>7)&1)==1)
            {
                if(avoir) _detec->_logText << "D__ cas repet"<<endl;
                if(((octet1>>6)&1)==1)
                {
                    valRepet = octet1 & 31;
                    if(((octet1>>5)&1)==1)
                    {
                        lastValue = _simpleParamsArray[i][j-1];
                        for(int j2=j;j2<=j+valRepet;j2++)
                            _simpleParamsArray[i][j2]=lastValue;
if(avoir) _detec->_logText << "D__  série de répétitions valRepet =" << valRepet << endl;
                        j+=valRepet;
                    }
                    else
                    {
                        if(j-1-valRepet>=0)
              _simpleParamsArray[i][j]=_simpleParamsArray[i][j-1-valRepet];
if(avoir) _detec->_logText << "D__  récupère le paramètre précédant de "
                           << valRepet << endl;
                    }
                }
                else
                {
                    valRepet = octet1 & 63;
                    if(i-1-valRepet>=0)
             _simpleParamsArray[i][j]=_simpleParamsArray[i-valRepet][j];
if(avoir) _detec->_logText << "D__  récupère le param. du cri précédant de "
                       << valRepet << endl;
                }
            } // fin des cas répétition
            // cas valeur simple (non repetition)
            else
            {

                signe = (octet1>>6)&1;
                if(((octet1>>5)&1)==1)
                {
                    entier=octet1 & 31;
                    if(avoir) _detec->_logText << "D__ cas < 32  entier="<<entier <<endl;
                }
                else
                {
                    entier=32;
                    _compressedParametersStream >> octet2;
                    entier+=(int)((octet1 & 31)*256)+(unsigned char)octet2;
                    if(avoir) _detec->_logText << "D__ cas>32  "
                             << " octet2=" << octet2
                             << " entier=" << entier
                             << endl;
                    if(entier==8223)
                    {
                        _compressedParametersStream >> sisuppl;
                        entier+=(unsigned qint16)sisuppl;
                        if(avoir) _detec->_logText << "D__ cas>8223  "
                                 << " sisuppl=" << sisuppl
                                 << " entier=" << entier
                                 << endl;
                        if(entier==73758)
                        {
                            _compressedParametersStream >> lisuppl;
                            entier+=(unsigned long int)lisuppl;
                        }
                    }
                }
                //
                _compressedParametersStream >> decimal;
                if(avoir) _detec->_logText << "D__  decimal = " << decimal << endl;
                float value = (float) entier +(float)decimal/100.0f;
                if(signe==0) value = -value;
                if(avoir) _detec->_logText << "D__  value = " << value << endl;

                if(i==0) _simpleParamsArray[i][j]=value;
                else _simpleParamsArray[i][j]=_simpleParamsArray[i-1][j]+value;
            } // fin cas valeur normale (non repetition)
            //if(i>0 || j>13) avoir=false;
            if(avoir) _detec->_logText << endl;
        } // fin boucle j
    } // fin boucle i
    // enregisstrement dans le fichier
    for(int i=0;i<nbcris;i++)
    {
        fileStream << wavFile << '\t' << QString::number(i);
        //_simplesParamsArray
        for(int j=0;j<_numberCallParameters;j++)
        {
            fileStream << '\t' <<  _simpleParamsArray[i][j];
        }
        fileStream << endl;
    }
    _compressedParametersFile.close();
    txtFile.close();
}
*/

void DetecTreatment::createImage(QString wavFile)
{
    if(_sonogramWidth > 32767) {_xmoitie = true; _imaWidth=(_sonogramWidth+1)/2;}
    else  {_xmoitie = false; _imaWidth=_sonogramWidth;}
    QImage ima = QImage(_imaWidth, _fftHeightHalf,QImage::Format_RGB32);
    // initBvrvb(_energyMin,(double)_medianNoise,_energyMax);
    initBvrvb(_energyMin,_energyMax);
    // refuser le traitement si imax > 2000
    int imax=((int)(_energyMax-_energyMin)*5)+1;
    uint crgb;
    for(int i=0;i<imax;i++)
    {
        _tvaleur[i]=calculateRGB((double)i/5.0f);
    }
    //£ float *ydl;
    qint16 *ydl;
    //& char *mzc;
    int exceptions = 0;
    //µ :
    int lyi = qMin(_fftHeightHalf,_limY);
    //µ for(int y = 0; y < _fftHeightHalf ; y++)
    for(int y = 0; y < lyi ; y++)
    {
        ydl=_sonogramArray[y];
        //& mzc=_pointFlagsArray[y];
        int digitPos = 0;
        char *pBoolChar = _pointFlagsArray[y];;
        char boolChar = *pBoolChar;

        for (int x = 0 ; x < _imaWidth ; x++)
        {
            //£ int valeur=(int)(           ( (*ydl)-_energyMin ) * 5.0f   );
            int valeur=(int)(((float)(*ydl)/20.0f)  - _energyMin*5.0f);
            if(valeur>=0 && valeur<imax)crgb=_tvaleur[valeur];
            else {crgb=0; exceptions++;}
            //& if(*mzc) crgb |= 224 << 16;
            if((boolChar & (1 << digitPos))!=0) crgb |= 224 << 16;
            ima.setPixel(x, _fftHeightHalf-y-1,crgb);
            //& ydl++;
            //& mzc++;
            //& if(_xmoitie) {ydl++; mzc++;}
            for(int k=0;k<1+_xmoitie;k++)
            {
                ydl++;
                digitPos++;
                if(digitPos==8) {pBoolChar++; boolChar = *pBoolChar; digitPos=0;}
            }
        }
    }
    //
    QString imageName = _imagePath + '/' + wavFile.replace(QString(".wav"), QString(".jpg"), Qt::CaseInsensitive);
    ima.save(imageName,0,100); // save image
}

// void DetecTreatment::initBvrvb(double bornemin,double median,double bornemax)
void DetecTreatment::initBvrvb(double bornemin,double bornemax)
{
    /*
    _bRGB[0][0]=bornemin;
    _bRGB[1][0]=median;
    _bRGB[2][0]=median+(double)_stopThreshold;
    _bRGB[3][0]=median+(double)_detectionThreshold;
    _bRGB[4][0]=bornemax+1;
    */
    _bRGB[0][0]=bornemin;
    // _bRGB[1][0]=0;
    // _bRGB[2][0]=(double)_stopThreshold;
    // _bRGB[3][0]=(double)_detectionThreshold;
    _bRGB[1][0]=bornemin/2;
    _bRGB[2][0]=(double)_energyStopThreshold;
    _bRGB[3][0]=(double)_energyShapeThreshold;
    _bRGB[4][0]=bornemax+1;
    for(int i=0;i<5;i++) _bRGB[i][0]-=bornemin;
    _bRGB[0][1]=0;   _bRGB[0][2]=0;   _bRGB[0][3]=32;
    _bRGB[1][1]=0;   _bRGB[1][2]=32;  _bRGB[1][3]=64;
    _bRGB[2][1]=0;   _bRGB[2][2]=64; _bRGB[2][3]=96;
    _bRGB[3][1]=0;   _bRGB[3][2]=96; _bRGB[3][3]=128;
    _bRGB[4][1]=0;   _bRGB[4][2]=128; _bRGB[4][3]=160;
}

uint DetecTreatment::calculateRGB(double value)
{
    double rgb[3];
    for(int j=0;j<3;j++) rgb[j] =0.0f;
    for(int i=0;i<4;i++)
    {
        if(value>=_bRGB[i][0] && value<_bRGB[i+1][0])
        {
            for(int j=0;j<3;j++)
                rgb[j]=_bRGB[i][j+1]
                        +((_bRGB[i+1][j+1]-_bRGB[i][j+1])*(value-_bRGB[i][0]))/(_bRGB[i+1][0]-_bRGB[i][0]);
            break;
        }
    }
    uint urgb = qRgb(qRound(rgb[0]),qRound(rgb[1]),qRound(rgb[2]));
    return(urgb);
}

void DetecTreatment::saveDatFile(QString wavFile)
{
    //_detec->_logText << "savedatfile début" << endl;
    QString da2file = wavFile.replace(QString(".wav"),QString(".da2"), Qt::CaseInsensitive);
    _callMatrixName = _datPath + '/' + da2file;
    _callMatrixFile.setFileName(_callMatrixName);
    _callMatrixFile.open(QIODevice::WriteOnly);
    _callMatrixStream.setDevice(&_callMatrixFile);
    _callMatrixStream << (int)_detec->_logVersion;
    _callMatrixStream << _detec->_userVersion;
    int nbcris = (int)_callsArray.size();
    _callMatrixStream << nbcris;
    _detec->_logText << "   nbcris = " << nbcris << endl;

    //µ :
    int lyi = qMin(_fftHeightHalf,_limY);

    //µ _callMatrixStream << _fftHeightHalf;
    _callMatrixStream << lyi;

    _callMatrixStream << (int)_xmoitie;
    _callMatrixStream << (float)_msPerX;
    _callMatrixStream << (float)_khzPerY;
   // _detec->_logText << "savedatfile 2" << endl;
    for(int i=0;i<nbcris;i++)
    {
        //_logText << "  cri num. = " << i << endl;
        _callMatrixStream << _masterPoints[i].x();
        _callMatrixStream << _masterPoints[i].y();
        int nbPixel=_callsArray[i].size();
        _callMatrixStream << nbPixel;
		int x,y; 
		float e;
        for(int j=0;j<nbPixel;j++)
        {
            x = _callsArray[i].at(j).x();
            y = _callsArray[i].at(j).y();
            //£ e = _sonogramArray[y][x];
            e = (float)_sonogramArray[y][x]/100.0f;
            _callMatrixStream << x;
            _callMatrixStream << y;
            _callMatrixStream << e;
        }
        int longCri=_callMasterRidgeArray[i].size();
        _callMatrixStream << longCri;
        for(int j=0;j<longCri;j++)
        {
            _callMatrixStream << _callMasterRidgeArray[i][j].x();
            _callMatrixStream << _callMasterRidgeArray[i][j].y();
        }
        longCri=_callSouthArray[i].size();
        _callMatrixStream << longCri;
        for(int j=0;j<longCri;j++)
        {
            _callMatrixStream << _callSouthArray[i][j].x();
            _callMatrixStream << _callSouthArray[i][j].y();
        }
        longCri=_callNorthRidgeArray[i].size();
        _callMatrixStream << longCri;
        for(int j=0;j<longCri;j++)
        {
            _callMatrixStream << _callNorthRidgeArray[i][j].x();
            _callMatrixStream << _callNorthRidgeArray[i][j].y();
        }
        longCri=_callWestRidgeArray[i].size();
        _callMatrixStream << longCri;
        for(int j=0;j<longCri;j++)
        {
            _callMatrixStream << _callWestRidgeArray[i][j].x();
            _callMatrixStream << _callWestRidgeArray[i][j].y();
        }
        longCri=_callSecondWestRidgeArray[i].size();
        _callMatrixStream << longCri;
        for(int j=0;j<longCri;j++)
        {
            _callMatrixStream << _callSecondWestRidgeArray[i][j].x();
            _callMatrixStream << _callSecondWestRidgeArray[i][j].y();
        }
        for(int jcrete=0;jcrete<NCRETES;jcrete++)
        {
           _callMatrixStream << _lowSlope[(i*NCRETES+jcrete)*2] << _lowSlope[(i*NCRETES+jcrete)*2+1];
           _callMatrixStream << _inflexion1[(i*NCRETES+jcrete)*2]  << _inflexion1[(i*NCRETES+jcrete)*2+1];
           _callMatrixStream << _inflexion3[(i*NCRETES+jcrete)*2]   << _inflexion3[(i*NCRETES+jcrete)*2+1];
        }
    } // next i
    //_logText << "fin de saveDatFile " << endl;
    // -------------------
    _callMatrixFile.close();
}
