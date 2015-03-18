#include "detec.h"

using namespace std;

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

Detec::Detec(QString dirPath): QThread()
{
    _wavPath = dirPath;
    _timeExpansion=10;
    _detectionThreshold = 26;
    _stopThreshold = 18;
    _freqMin = 0;
    _nbo = 4;
    _freqCallMin=8.0f;
    initVectorParams();
    _initPassed = InitializeDetec();
}

Detec::~Detec()
{
    _logFile.close();
}

bool Detec::InitializeDetec()
{
    IsRunning = true;
    QString logFilePath("detec.log");
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    _logText.setDevice(&_logFile);
    //
    QDir sdir(_wavPath);
    if(!sdir.exists())
    {
        _logText << "répertoire "<< _wavPath << " non trouvé !" << endl;
        return(false);
    }
    _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
    if(_wavFileList.isEmpty()) return(false);
    //
    _sonogramArray = new float*[FFT_HEIGHT_MAX];
    _pointFlagsArray = new char *[FFT_HEIGHT_MAX];
    float *smy;
    for ( int i=0 ; i < FFT_HEIGHT_MAX ; i++)
    {
        _sonogramArray[i] = new float[SONOGRAM_WIDTH_MAX];
        _pointFlagsArray[i] = new char[SONOGRAM_WIDTH_MAX];
        smy=_sonogramArray[i];
        for(int j = 0; j < SONOGRAM_WIDTH_MAX; j++) *smy++ = 0.0f;
    }
    //_txtPath = _wavPath+"/txt";
    _txtPath = sdir.absolutePath()+"/txt";
    QDir reptxt(_txtPath);
    if(!reptxt.exists())
    {
        if(!reptxt.mkdir(_txtPath))
        {
            _logText << "création du sous-répertoire "<< _txtPath << " impossible !" << endl;
            return(false);
        }
    }
    _energyMin = 0.0f;
    QString fileInfo = _txtPath+"/fileInfo.txt";
    _fileInfo.setFileName(fileInfo);
    _fileInfo.open(QIODevice::WriteOnly | QIODevice::Text);
    _fIStream.setDevice(&_fileInfo);
    _fIStream.setRealNumberNotation(QTextStream::FixedNotation);
    _fIStream.setRealNumberPrecision(6);
    _fIStream << "Filename\tSR\tDur\tMAX1\tMED1\tMAX2\tMED2\tMAX3\tMED3\tMAX4\tMED4\tMAX5\tMED5\tMAX6\tMED6\tMAX7\tMED7\tMAX8\tMED8\tMAX9\tMED9\tMAX10\tMED10\n";
    QString errorFilePath(_txtPath+"/error.log");
    _errorFile.setFileName(errorFilePath);
    _errorFile.open(QIODevice::WriteOnly | QIODevice::Text);
    _errorStream.setDevice(&_errorFile);
    _errorStream << "Error List" << endl;
    _charParamsArray = new char[MAXCRI*NBTAB*NBPAR*sizeof(float)];
    _paramsArray = new float**[MAXCRI];
    for (int i=0;i < MAXCRI; i++)
    {
        _paramsArray[i] = new float *[NBTAB];
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
        _ypm[k]=new int[MAXCRI];
    }
    _tabY = new float[MAXHAUCRI];
    _numberPixelsPerY = new int[MAXHAUCRI];
    _numberPixelsPerX = new int[MAXLARCRI];
    _averagePerX = new float[MAXLARCRI];
    _xMinPerY = new int[MAXHAUCRI];
    _xSecondWestRidgePerY = new int[MAXHAUCRI];
    _xMaxPerY = new int[MAXHAUCRI];
    _yMinPerX = new int[MAXLARCRI];
    _yMaxPerX = new int[MAXLARCRI];
    _eMaxPerX = new float[MAXLARCRI];
    _slope    = new float[qMax(MAXLARCRI,MAXHAUCRI)];
    _charTabX = new char[MAXCRI*MAXLARCRI*sizeof(float)];
    _tabX = new float*[MAXCRI];
    for(int i=0 ; i < MAXCRI ; i++) _tabX[i] = (float *)((char *)(_charTabX+(i*MAXLARCRI*sizeof(float))));
    _charYEmaxPerX = new char[MAXCRI*MAXLARCRI*sizeof(int)];
    _yEmaxPerX = new int*[MAXCRI];
    for(int i=0;i<MAXCRI;i++) _yEmaxPerX[i] = (int *)((char *)(_charYEmaxPerX+(i*MAXLARCRI*sizeof(int))));
    _charTabYX = new char[MAXHAUCRI*MAXLARCRI*sizeof(float)];
    _tabYX = new float*[MAXHAUCRI];
    for(int i=0;i<MAXHAUCRI;i++) _tabYX[i] = (float *)((char *)(_charTabYX+(i*MAXLARCRI*sizeof(float))));
    return true;
}

void Detec::initVectorParams()
{
    QString prefix[] = {"","CM_","CN_","CS_","CO_","CO2_"};
    _vectPar.push_back(ParamToSave(SH,StartTime,"St(ms)"));
    _vectPar.push_back(ParamToSave(SH,Dur,"Dur(ms)"));
    _vectPar.push_back(ParamToSave(SH,Prev,"Prev(ms)"));
    _vectPar.push_back(ParamToSave(SH,Fmax,"Fmax(kHz)"));
    _vectPar.push_back(ParamToSave(SH,Fmin,"Fmin(kHz)"));
    _vectPar.push_back(ParamToSave(SH,BW,"BW(kHz)"));
    _vectPar.push_back(ParamToSave(SH,FreqMP,"FreqMP(kHz)"));
    _vectPar.push_back(ParamToSave(SH,PosMP,"PosMP(%)"));
    _vectPar.push_back(ParamToSave(SH,FreqDomSum,"FreqDomSum(kHz)"));
    _vectPar.push_back(ParamToSave(SH,FreqDomMean,"FreqDomMean(kHz)"));
    _vectPar.push_back(ParamToSave(SH,PosPeakSum,"PosPeakSum(kHz)"));
    _vectPar.push_back(ParamToSave(SH,PosPeakMean,"PosPeakMean(kHz)"));
    _vectPar.push_back(ParamToSave(SH,FreqPeakSum,"FreqPeakSum(kHz)"));
    _vectPar.push_back(ParamToSave(SH,FreqPeakMean,"FreqPeakMean(kHz)"));
    _vectPar.push_back(ParamToSave(SH,PrevMP,"PrevMP (ms)"));
    _vectPar.push_back(ParamToSave(SH,PrevSmart1,"PrevSmart1 (ms)"));
    _vectPar.push_back(ParamToSave(SH,NextMP,"NextMP (ms)"));
    _vectPar.push_back(ParamToSave(SH,NextMP,"NextSmart1 (ms)"));
    _vectPar.push_back(ParamToSave(SH,Amp1,"Amp1"));
    _vectPar.push_back(ParamToSave(SH,Amp2,"Amp2"));
    _vectPar.push_back(ParamToSave(SH,Amp3,"Amp3"));
    _vectPar.push_back(ParamToSave(SH,Amp4,"Amp4"));
    _vectPar.push_back(ParamToSave(SH,NoiseLeft,"NoiseLeft"));
    _vectPar.push_back(ParamToSave(SH,NoiseRight,"NoiseRight"));
    _vectPar.push_back(ParamToSave(SH,NoiseDown,"NoiseDown"));
    _vectPar.push_back(ParamToSave(SH,NoiseUp,"NoiseUp"));
    _vectPar.push_back(ParamToSave(SH,CVAmp,"CVAmp"));
    _vectPar.push_back(ParamToSave(CO,Dur,"CO_Dur(ms)"));
    _vectPar.push_back(ParamToSave(CO2,Dur,"CO2_Dur(ms)"));
    _vectPar.push_back(ParamToSave(CM,Fmax,"CM_Fmax(kHz)"));
    _vectPar.push_back(ParamToSave(CM,Fmin,"CM_Fmin(kHz)"));
    _vectPar.push_back(ParamToSave(CM,BW,"CM_BW(kHz)"));
    _vectPar.push_back(ParamToSave(CM,Fdom,"CM_Fdom(kHz)"));
    _vectPar.push_back(ParamToSave(CM,FdomDiff,"CM_FdomDiff(kHz)"));
    _vectPar.push_back(ParamToSave(CM,Ldom,"CM_Ldom(Percent)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Slope,prefix[i]+"Slope(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,HiFc,prefix[i]+"HiFc(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Lhi,prefix[i]+"Lhi(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,LowFc,prefix[i]+"LowFc(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Fc3,prefix[i]+"Fc3(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Su,prefix[i]+"Su(kHz/s)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,Sl,prefix[i]+"Sl(kHz/s)"));
    _vectPar.push_back(ParamToSave(CM,StartF,"CM_StartF(kHz)"));
    _vectPar.push_back(ParamToSave(CM,EndF,"CM_EndF(kHz)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,StartSlope,prefix[i]+"StartSlope(kHz/ms)"));
    for(int i=CM;i<=CO2;i++) _vectPar.push_back(ParamToSave(i,EndSlope,prefix[i]+"EndSlope(kHz/ms)"));
    _vectPar.push_back(ParamToSave(CM,SlopeAtFc,"CM_SlopeAtFc (kHz)"));
    for(int i=CM;i<=CN;i++) _vectPar.push_back(ParamToSave(i,FreqCtr,prefix[i]+"FreqCtr (kHz)"));
    _vectPar.push_back(ParamToSave(CM,FBak5dB,"CM_FBak5dB (kHz)"));
    _vectPar.push_back(ParamToSave(CM,FFwd5dB,"CM_FFwd5dB (kHz)"));
    _vectPar.push_back(ParamToSave(CM,Bndw5dB,"CM_Bndw5dB (kHz)"));
    _vectPar.push_back(ParamToSave(CM,DurOf5dB,"CM_DurOf5dB (kHz)"));
    _vectPar.push_back(ParamToSave(SH,Hup_RFMP,"Hup_RFMP"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosMP,"Hup_PosMP"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosSt,"Hup_PosSt"));
    _vectPar.push_back(ParamToSave(SH,Hup_PosEn,"Hup_PosEn"));
    _vectPar.push_back(ParamToSave(SH,Hup_RAmp,"Hup_RAmp"));
    _vectPar.push_back(ParamToSave(SH,Hup_RSlope,"Hup_RSlope"));
    _vectPar.push_back(ParamToSave(SH,Hlo_RFMP,"Hlo_RFMP"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosMP,"Hlo_PosMP"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosSt,"Hlo_PosSt"));
    _vectPar.push_back(ParamToSave(SH,Hlo_PosEn,"Hlo_PosEn"));
    _vectPar.push_back(ParamToSave(SH,Hlo_RAmp,"Hlo_RAmp"));
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
    _vectPar.push_back(ParamToSave(SH,RPM8,"RPM8"));
    _vectPar.push_back(ParamToSave(SH,Stab,"Stab"));

}

void Detec::endDetec()
{
    if(_initPassed)
    {
        _errorFile.close();
        _fileInfo.close();
        for ( int i=0 ; i < FFT_HEIGHT_MAX ; i++)
        {
            delete[] _sonogramArray[i];
            delete[] _pointFlagsArray[i];
        }
        delete[] _sonogramArray;
        delete[] _pointFlagsArray;
        for (int i=0;i < MAXCRI; i++) delete[] _paramsArray[i];
        delete[] _paramsArray;
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
        delete[] _eMaxPerX;
        delete[] _xSecondWestRidgePerY;
        delete[] _slope;
        delete[] _tabX;
        delete[] _charTabX;
        delete[] _yEmaxPerX;
        delete[] _charYEmaxPerX;
        delete[] _tabYX;
        delete[] _charTabYX;
    }
    IsRunning = false;
}


void Detec::run()
{
    if(_initPassed)
    {
        foreach(QString wavFile,_wavFileList) treatOneFile(wavFile);
    }
    endDetec();
}

void Detec::treatOneFile(QString wavFile)
{
    _logText << "Début de traitement : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    _logText << "wavfile=" << wavFile << endl;
    QString nomfic = wavFile;
    for(int j=0;j<FFT_HEIGHT_MAX;j++) memset(_pointFlagsArray[j],0,SONOGRAM_WIDTH_MAX);
    QString pathFile = _wavPath + '/' + wavFile;
    if (openWavFile(pathFile))
    {
        _logText << "Avant computefft : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        if(computeFFT(pathFile))
        {
            _logText << "Avant CorrectNoise : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            correctNoise();
            /*
            _logText << "Avant CalculeBruitMedian : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            calculateMedianNoise();
            */
            _logText << "Avant detectsContours : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            shapesDetects();
            _callsNumber = (int)_callsArray.size();
            _logText << "Avant detectsparameter2 : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            detectsParameter2();
            _logText << "Avant enregistrement param2.txt : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            saveParameters(wavFile);
            _logText << "Fin de traitement : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        }
        else
        {
            // TODO : gérer erreur dans computefft
            return;
        }
    }
    else
    {
        // TODO : gérer erreur dans computefft
        return;
    }
    // -------------------------------------------------------
    _logText << "Fin de traitement d'un fichier : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
}

bool Detec::openWavFile(QString& wavFile)
{

    if (! (_soundFile = sf_open(wavFile.toStdString().c_str(), SFM_READ, &_soundFileInfo)))
    {
        _errorStream << wavFile << ": " << sf_strerror (NULL) << endl;
        return  false;
    }
    if (_soundFileInfo.channels > 1)
    {
        sf_close (_soundFile);
        _errorStream << wavFile << ": " << "multi-channel non traité" << endl;
        return  false;
    }
    _fIStream << wavFile << '\t' << _soundFileInfo.samplerate << '\t';
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

bool Detec::computeFFT(QString &wavFile)
{
    int iCount;
    int readcount;
    float a = 0.0f;
    _energyMax        = 0.0f;
    _energyMin        = 0.0f;
    _fftHeightHalf		= (int)ceil((float)_fftHeight/(float)2);

    _coeff = new float[_fftHeightHalf];
    _iOverlapMoving	= (int)ceil((float)_fftHeight/(float)(_nbo*2));
    _sonogramWidth		= (int)ceil(_nbo*2*(float)_soundFileInfo.frames/(float)_fftHeight)+1;
    _data				= ( float* ) fftwf_malloc( sizeof( float ) * _fftHeight );
    _fftRes 		= ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * _fftHeight );
    _complexInput        = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * _fftHeight );
    _msPerX =(float)(_fftHeightHalf*1000)/(_nbo*_soundFileInfo.samplerate*_timeExpansion); //Time: msec
    _fIStream << _sonogramWidth*_msPerX/1000 << '\t';
    _khzPerY =(float)(_soundFileInfo.samplerate*_timeExpansion)/(float)(_fftHeight*1000); //Freq:khz
    _logText << "_sonogramWidth = " << _sonogramWidth << endl;
    if(_sonogramWidth*_msPerX < 10.0f)
    {
        _errorStream << wavFile << ": fichier trop petit" << endl;
        _logText << "durée trop petite : " << _sonogramWidth*_msPerX << " ms" << endl;
        return  false;
    }
    if(_sonogramWidth > SONOGRAM_WIDTH_MAX)
    {
        _errorStream << wavFile << ": fichier trop grand" << endl;
        return  false;
    }
    _plan = fftwf_plan_dft_1d( _fftHeight, _complexInput, _fftRes, FFTW_FORWARD, FFTW_ESTIMATE );
    float fact1=2.0f*PI;
    float fact2=4.0f*PI;
    float quot1=_fftHeightHalf-1;
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
            float *sml;
            int jc=iCount*_nbo+iLoop;
            float b;
            for(int i =0; i < _fftHeightHalf; i++)
            {
                sml=_sonogramArray[i];
                a = pow(_fftHeight*_fftRes[i][0], 2) + pow(_fftHeight*_fftRes[i][1], 2);
                if (i > _freqMin && a !=0 )
                {
                    b=10*log10(a);
                    sml[jc]=b;
                    if(b>_energyMax) _energyMax=b;
                    else
                    {
                        if(b<_energyMin) _energyMin =b;
                    }
                }
                else sml[jc] = -50;
            }
            iCount++;
        }
    }
    for(int i =0; i < _fftHeightHalf; i++) _sonogramArray[i][_sonogramWidth-1]=0;
    _energyMax = qMax(_energyMax, (double)0);
    sf_close (_soundFile);
    fftwf_free(_data);
    fftwf_free(_fftRes);
    fftwf_free(_complexInput);
    delete _coeff;
    return true;
}

void Detec::shapesDetects()
{
    QPoint Point;
    int ix,iy;
    int curseur;
    _logText << "début DetectsContours" << endl ;
    bool on_en_a_un = false;int ncontour=0;

    _logText << "début DetectsContours" << endl ;

    _maxCallWidth = 0;
    _maxCallHeight = 0;

    _minY=(int)(((float)_freqMin)/((float)_khzPerY));
    _maxY=(int)(((float)FREQ_MAX)/((float)_khzPerY));
    if(_maxY>_fftHeightHalf-1) _maxY = _fftHeightHalf-1;
    /*
    _energyShapeThreshold = (double)_medianNoise+(double)_detectionThreshold;
    _energyStopThreshold = (double)_medianNoise+(double)_stopThreshold;
    */
    _energyShapeThreshold = (double)_detectionThreshold;
    _energyStopThreshold = (double)_stopThreshold;


    int nbcont=0;
    float *smy;
    char *zcy;
    for(int y = _maxY; y >= _minY ; y--)
    {
        smy=_sonogramArray[y];
        zcy=_pointFlagsArray[y];
        for (int x = 0 ; x < _sonogramWidth ; x++)
        {
            if(zcy[x]==0)
                if (smy[x] > _energyShapeThreshold)
                {
                    on_en_a_un = true;
                    nbcont++;
                    ncontour++;
                    _vectorCallPoints.clear();
                    Point.setX(x);
                    Point.setY(y);
                    _callEnergyMax = smy[x];
                    _callEnergyMaxIndex = 0;
                    _vectorCallPoints.push_back(Point);
                    _xMin=x; _xMax=x;
                    _yMin=y; _yMax=y;
                    curseur = 0;
                    zcy[x] = 1;

                    while(curseur < _vectorCallPoints.size())
                    {
                        Point=_vectorCallPoints.at(curseur);
                        ix=Point.x();
                        iy=Point.y();
                        for(int jy=iy-1;jy<=iy+1;jy++)
                            // edit yves - elargir spectre
                            for(int jx=ix-5;jx<=ix+5;jx++)
                                if(jx!=ix || jy!=iy)
                                {
                                    if(jx>= 0 && jx < _sonogramWidth && jy >= _minY && jy <= _maxY)
                                    {
                                        if(_pointFlagsArray[jy][jx]==0)
                                        {
                                            float val = _sonogramArray[jy][jx];
                                            if (val > _energyStopThreshold)
                                            {
                                                Point.setX(jx);
                                                Point.setY(jy);
                                                _pointFlagsArray[jy][jx] = 1;
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
                        curseur++;
                        if(curseur > 100000) break;
                    }
                    float freqpm=((float)_vectorCallPoints.at(_callEnergyMaxIndex).y())*_khzPerY;

                    if(freqpm>_freqCallMin && (_xMax-_xMin+1)<MAXLARCRI
                            && (_yMax-_yMin+1)<MAXHAUCRI)
                    {
                        _masterPoints.push_back(_vectorCallPoints.at(_callEnergyMaxIndex));
                        _callsArray.push_back(_vectorCallPoints);
                        _vectorXMin.push_back(_xMin);
                        _vectorXMax.push_back(_xMax);
                        if((_xMax-_xMin+1)>_maxCallWidth) _maxCallWidth=_xMax-_xMin+1;
                        if((_yMax-_yMin+1)>_maxCallHeight) _maxCallHeight=_yMax-_yMin+1;
                    }
                }
        }
    }
    _logText << "Avant retrieCris : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    sortWaves();
}

void Detec::sortWaves()
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


// edit yves - calcul bruit par bandes
// edit yves - corrige l'image / heterogeneité du bruit / fréquence
// calcule pour chaque tranche de frequence un quantile indicateur du bruit de fond
// on soustrait ce bruit de fond tout en multipliant par ce niveau de bruit de fond pour tenir compte partiellement de l'echelle log des dB
void Detec::correctNoise()
{
    int son_min = -100,son_max = 99;
    int tval[200];
    float *fc;
    for(int y = 0; y < _fftHeightHalf ; y++)
    {
        fc=_sonogramArray[y];
        for(int k=0;k<200;k++) tval[k]=0;

        for (int x = 0 ; x < _sonogramWidth ; x++)
        {
            int son = qRound(fc[x]);
            if(son>=son_min && son<=son_max)
                tval[son-son_min]++;
        }

        int cumul = 0, q5 = _sonogramWidth*5/100;
        bool oncherche = true;
        for(int j=0;j<=son_max-son_min;j++)
        {
            cumul += tval[j];
            if(cumul>=q5 && oncherche)
            {
                for (int x = 0 ; x < _sonogramWidth ; x++) fc[x] = fc[x]-j-son_min;
                oncherche=false;
                break;
            }
        }
    }
}

void Detec::detectsParameter2()
{
    float *oParam;
    float *oParamCrete[NCRETES];
    //_logText << "dp2-debut";
    int nbcris = _callsArray.size();
    if(nbcris< 1 || _maxCallWidth < 1 || _maxCallHeight < 1) return;
    int maxlarhau = _maxCallWidth;
    if(_maxCallHeight>maxlarhau) maxlarhau=_maxCallHeight;
    //_logText << "-2" ;
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
            _yEmaxPerX[icri][k]=0;
            _eMaxPerX[k]=0.0f;
        }
        float eTot = 0.0f;
        float eMoy = 0.0f;
        for(int j=0;j<tailleforme;j++)
        {
            int x=unemat.at(j).x();
            int y=unemat.at(j).y();
            float e = _sonogramArray[y][x];
            _tabY[y-ymin]+=e;
            eTot += e;
            _numberPixelsPerY[y-ymin]++;
            _tabX[icri][x-xmin]+=e;
            _numberPixelsPerX[x-xmin]++;
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
                if(_tabX[icri][k]/_numberPixelsPerX[k]>recppm)
                {
                    xppm=xmin+k;
                    recppm=_tabX[icri][k]/_numberPixelsPerX[k];
                }
        }
        oParam[StartTime] = (float)xmin*(float)_msPerX;
        oParam[Dur] = (xmax-xmin)*(float)_msPerX;
        int criprec=icri;
        if(icri>0)
            for(int jcri=icri-1;jcri>=0;jcri--)
            {
                if(_paramsArray[jcri][SH][StartTime] + _paramsArray[jcri][SH][Dur]<oParam[StartTime])
                {
                    criprec=jcri;
                    break;
                }
            }
        if(icri == 0 || criprec==icri) oParam[Prev] = 9999.0f;
        else oParam[Prev] = oParam[StartTime] - _paramsArray[criprec][SH][StartTime];
        oParam[Fmax]=ymax*(float)_khzPerY;
        oParam[Fmin]=ymin*(float)_khzPerY;
        oParam[BW]=oParam[Fmax]-oParam[Fmin];
        oParam[FreqMP] = (float)_masterPoints.at(icri).y()*_khzPerY;
        oParam[PosMP] = (float)(_masterPoints.at(icri).x()-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[FreqDomSum]  = yfds * _khzPerY;
        oParam[FreqDomMean] = yfdm * _khzPerY;
        oParam[PosPeakSum]   = (float)(xpps-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[PosPeakMean]  = (float)(xppm-xmin+0.5f)/(float)(xmax-xmin+1.0f);
        oParam[FreqPeakSum] = (float) _yEmaxPerX[icri][xpps-xmin] * _khzPerY;
        oParam[FreqPeakMean] = (float) _yEmaxPerX[icri][xppm-xmin] * _khzPerY;
        float prevmp = 9999.0f;
        float prevsmart1 =9999.0f;
        //
        //_logText << "-5";
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
        oParam[PrevMP] = prevmp;
        oParam[PrevSmart1] = prevsmart1;
        oParam[NextMP] = 9999.0f;
        oParam[NextSmart1] = 9999.0f;
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
                    if(ifin==ideb) prorata=1;
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
        if(namp[0]>0) oParam[Amp1] = famp[0]/namp[0];
        else oParam[Amp1] = 0.0f;
        if(namp[1]>0) oParam[Amp2] = famp[1]/namp[1];
        else oParam[Amp2] = 0.0f;
        if(namp[2]>0) oParam[Amp3] = famp[2]/namp[2];
        else oParam[Amp3] = 0.0f;
        if(namp[3]>0) oParam[Amp4] = famp[3]/namp[3];
        else oParam[Amp4] = 0.0f;
        for(int inoise=0;inoise<4;inoise++)
        {
            int jdeb,jfin,x,y;
            int nps=0;
            float bruit=0.0f,bruit_moyen=0.0f;

            bool ontraite = true;
            if(ontraite)
            {
                if(inoise<2) {jdeb=ymin;jfin=ymax;}
                else {jdeb=xmin;jfin=xmax;}
                for(int j=jdeb;j<=jfin;j++)
                {
                    for(int k=1;k<=3;k++)
                    {
                        if((inoise<2 && _xMaxPerY[j-ymin]>0) || (inoise>1 &&  _yMaxPerX[j-xmin]>0))
                        {
                            if(inoise==0){x=_xMinPerY[j-ymin]-k;y=j;}
                            if(inoise==1){x=_xMaxPerY[j-ymin]+k;y=j;}
                            if(inoise==2){y=_yMinPerX[j-xmin]-k;x=j;}
                            if(inoise==3){y=_yMaxPerX[j-xmin]+k;x=j;}
                            if(x>=0 && x<_sonogramWidth && y>0 && y<_fftHeightHalf)
                            {
                                bruit+=_sonogramArray[y][x];
                                nps++;
                            }
                        }
                    }
                }

                if(nps>0) bruit_moyen=bruit/((float)nps);
                //
            }
            if(inoise==0)oParam[NoiseLeft]  = bruit_moyen;
            if(inoise==1)oParam[NoiseRight] = bruit_moyen;
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
        int *pc; // pointeur sur la crete
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
                        amp[ih] += _sonogramArray[y][x];
                        np++;
                    }
                }
                for(int gd=0;gd<2;gd++) for(int k=1;k<=3;k++)
                {
                    if(gd==0)xn=_xMinPerY[j-ymin]-k;
                    else xn=_xMaxPerY[j-ymin]+k;
                    if(xn>=0 && xn<_sonogramWidth)
                    {
                        ran[ih]+=_sonogramArray[y][xn];
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
                            ran[ih]+=_sonogramArray[yn][x];
                            nn++;
                        }
                    }
                }
            }
            if(np>0) ampmoy[ih]=amp[ih]/((float)np);
            if(nn>0) ranmoy[ih]=ran[ih]/((float)nn);
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
        oParam[HetX] = 0.0f;
        int ntr = 0;
        int nen = 0;
        float e1,e2,e3;
        for(int y=ymin;y<=ymax;y++)
        {
            for(int x=_xMinPerY[y-ymin];x<=_xMaxPerY[y-ymin];x++)
            {
                e2=_tabYX[y-ymin][x-xmin];
                if(e2>0.0f && x>0 && x<_sonogramWidth-1)
                {
                    ntr++;
                    e1=_sonogramArray[y][x-1];
                    e3=_sonogramArray[y][x+1];
                    if((e2>e1 && e2>e3) || (e2<e1 && e2<e3)) nen++;
                }
            }
        }
        if(ntr>0) oParam[HetX] = ((float) nen)/ ((float) ntr);
        oParam[HetY] = 0.0f;
        ntr=0;nen=0;
        for(int x=xmin;x<=xmax;x++)
        {
            if(_yMaxPerX[x-xmin]>0)
            {
                for(int y=_yMinPerX[x-xmin];y<=_yMaxPerX[x-xmin];y++)
                {
                    e2=_tabYX[y-ymin][x-xmin];
                    if(e2>0.0f && y>1 && y<_maxY)
                    {
                        ntr++;
                        e1=_sonogramArray[y-1][x];
                        e3=_sonogramArray[y+1][x];
                        if((e2>e1 && e2>e3) || (e2<e1 && e2<e3)) nen++;
                    }
                }
            }
        }
        //_logText << "-8" ;
        if(ntr>0) oParam[HetY] = ((float) nen)/ ((float) ntr);
        //
        oParam[RPM8] = 0.0f;
        float e8 = 0.0f;
        int n8 = 0;
        int Y8 = (int)(_freqCallMin / _khzPerY);
        for(int x=xmin;x<=xmax;x++)
            for(int y=0;y<Y8;y++)
            {
                e8 =+ _sonogramArray[y][x];
                n8++;
            }
        oParam[RPM8] = eMoy - (e8/n8);
        //
        oParam[Stab] = 0.0f;
        int xmaitr = _masterPoints.at(icri).x();
        int ymaitr = _masterPoints.at(icri).y();
        float dif,p;
        float ponderTot = 0.0f;
        float difTot = 0.0f;
        int dxmax = qMax(xmaitr-xmin,xmax-xmaitr);
        int dymax = qMax(ymaitr-ymin,ymax-ymaitr);
        int distLim = pow(dxmax,2)+pow(dymax,2);
        //
        for(int y=ymin;y<=ymax;y++)
        {
            for(int x=_xMinPerY[y-ymin]-1;x<=_xMaxPerY[y-ymin];x++)
            {
                if(x>=0 && x<_sonogramWidth-1)
                {
                    dif = qAbs(_sonogramArray[y][x+1]-_sonogramArray[y][x]);
                    p = 1.0f - (((pow(x-xmaitr,2)+pow(y-ymaitr,2))/distLim) * 0.75f);
                    difTot += dif * p;
                    ponderTot += p;
                }
            }
        }
        for(int x=xmin;x<=xmax;x++)
        {
            for(int y=_yMinPerX[x-xmin]-1;y<=_yMaxPerX[x-xmin];y++)
            {
                if(y>=_minY && y<_maxY-1)
                {
                    dif = qAbs(_sonogramArray[y+1][x]-_sonogramArray[y][x]);
                    p = 1.0f - (((pow(x-xmaitr,2)+pow(y-ymaitr,2))/distLim) * 0.75f);
                    difTot += dif * p;
                    ponderTot += p;
                }
            }
        }
        oParam[Stab] = difTot / ponderTot;
        //
        if(_imageData)
        {
            _callMasterRidge.clear();
            _callSouthRidge.clear();
            _callNorthRidge.clear();
            _callWestRidge.clear();
            _callSecondWestRidge.clear();
        }
        int lowfc,hifc,fc3;
        for(int jcrete=0;jcrete<NCRETES;jcrete++)
        {
            if(jcrete==0) {pc=_yEmaxPerX[icri];}
            if(jcrete==1) {pc=_yMinPerX;}
            if(jcrete==2) {pc=_yMaxPerX;}
            if(jcrete==3) {pc=_xMinPerY;}
            if(jcrete==4) {pc=_xSecondWestRidgePerY;}
            oParamCrete[jcrete][Slope] = 9999.0f;
            float emcmax = -50.0f;
            xmcmax[jcrete] = 0;
            ymcmax[jcrete]=0;
            if(jcrete<3)
            {
                alasuite = 0; meilleuresuite = 0; meilleur = xmin; meilleurepente = 100000.0f;
                int dmeil = 0;
                for(int x=xmin;x<=xmax;x++)
                {
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
                        {meilleur = x; meilleuresuite = alasuite; dmeil =1;}
                    }
                    else
                    {
                        alasuite = 0;
                        if(_slope[x-xmin] < meilleurepente)
                        {
                            meilleur = x;
                            meilleurepente = _slope[x-xmin];
                            dmeil = 2;
                        }
                    }
                    float e=_sonogramArray[pc[x-xmin]][x];
                    if(e>emcmax)
                    {
                        emcmax = e;
                        xmcmax[jcrete]=x;
                        ymcmax[jcrete] = pc[x-xmin];
                    }
                }
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
                oParamCrete[jcrete][LowFc]= pc[lowfc]*_khzPerY;
                float pe1,pe2,dif;
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
                oParamCrete[jcrete][HiFc] = pc[hifc]*_khzPerY;
                oParamCrete[jcrete][Lhi]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][Lhi] = ((float)meilleurk[0]+0.5f)/((float)(xmax-xmin+1));
                fc3=meilleurk[1];
                _inflexion3[(icri*NCRETES+jcrete)*2]=fc3+xmin;
                _inflexion3[(icri*NCRETES+jcrete)*2+1]=pc[fc3];
                oParamCrete[jcrete][Fc3] = pc[fc3]*_khzPerY;
                oParamCrete[jcrete][Su] = 0.0f;
                if(meilleurk[0]>0)
                    oParamCrete[jcrete][Su] = ((float)(pc[meilleurk[0]]-pc[0])*_khzPerY)/((float)meilleurk[0]*_msPerX);
                oParamCrete[jcrete][Sl] = 0.0f;
                if(lowfc-meilleurk[0]>0)
                    oParamCrete[jcrete][Sl] = ((float)(pc[lowfc]-pc[meilleurk[0]])*_khzPerY)/((float)(lowfc-meilleurk[0])*_msPerX);
                int debp,finp,ecart;
                float pj[3];
                for(int j2=0;j2<3;j2++)
                {
                    pj[j2]=0.0f;
                    if(j2<2)
                    {
                        ecart=(xmax-xmin+10)/20;
                        if(ecart==0 && xmax>xmin) ecart = 1;
                        if(j2==0) {debp=0; finp=ecart;}
                        else  {debp=xmax-xmin-ecart; finp=xmax-xmin;}
                    }
                    else
                    {
                        debp=xmcmax[jcrete]-xmin-2;
                        finp=xmcmax[jcrete]-xmin+2;
                        if(debp<0) debp =0;
                        if(finp>xmax-xmin) finp = xmax-xmin;
                        ecart=finp-debp;
                    }
                    if(ecart>0)
                    {
                        pj[j2]=((float)(pc[finp]-pc[debp])*_khzPerY)/((float)ecart*_msPerX);
                    }
                }
                oParamCrete[jcrete][StartSlope] = pj[0];
                oParamCrete[jcrete][EndSlope]   = pj[1];
                oParamCrete[jcrete][SlopeAtFc]  = pj[2];
                int milieu = (xmin+xmax)/2;
                for(int k=milieu;k<=xmax;k++)
                    if(pc[k-xmin]>0)
                    {
                        milieu=k;
                        break;
                    }

                oParamCrete[jcrete][FreqCtr] = (float)pc[milieu-xmin]*_khzPerY;
                if(jcrete==0)
                {
                    int pos1=xmcmax[jcrete];
                    int pos2=xmax;
                    float borne = emcmax - 5;
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
                    oParamCrete[jcrete][FBak5dB]  = pc[pos1-xmin]*_khzPerY;
                    oParamCrete[jcrete][FFwd5dB]  = pc[pos2-xmin]*_khzPerY;
                    oParamCrete[jcrete][Bndw5dB]  = oParamCrete[jcrete][FFwd5dB]-oParamCrete[jcrete][FBak5dB];
                    oParamCrete[jcrete][DurOf5dB] = (pos2-pos1)*_msPerX;
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
                    float e=_sonogramArray[y][pc[y-ymin]];
                    if(e>emcmax)
                    {
                        emcmax = e;
                        xmcmax[jcrete]=pc[y-ymin];
                        ymcmax[jcrete] = y;
                    }

                }
                if(pc[ymax-ymin]!=pc[0])
                {
                    oParamCrete[jcrete][Slope] = ((float)(ymax-ymin)*_khzPerY)/((float)(pc[ymax-ymin]-pc[0])*_msPerX);
                }
                else oParamCrete[jcrete][Slope] = 9999.0f;
                lowfc = meilleur;
                _lowSlope[(icri*NCRETES+jcrete)*2]=pc[lowfc-ymin];
                _lowSlope[(icri*NCRETES+jcrete)*2+1]=lowfc;
                oParamCrete[jcrete][LowFc] = lowfc*_khzPerY;
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
                oParamCrete[jcrete][HiFc] = hifc*_khzPerY;
                oParamCrete[jcrete][Lhi]  = 0.5f;
                if(xmax-xmin>0) oParamCrete[jcrete][Lhi] = ((float)(pc[meilleurk[0]-ymin]-xmin))/((float)(xmax-xmin));
                fc3=meilleurk[1];
                _inflexion3[(icri*NCRETES+jcrete)*2]=pc[fc3-ymin];
                _inflexion3[(icri*NCRETES+jcrete)*2+1]=fc3;
                oParamCrete[jcrete][Fc3] = fc3*_khzPerY;
                oParamCrete[jcrete][Su] = 0.0f;
                if((pc[meilleurk[0]-ymin]-pc[0]!=0 && sens==1)
                        || (pc[meilleurk[0]-ymin]-pc[ymax-ymin]!=0 && sens==-1))
                {
                    if(sens==1)
                        oParamCrete[jcrete][Su] = ((float)(meilleurk[0]-ymin)*_khzPerY)
                                /( ( (float)(pc[meilleurk[0]-ymin]-pc[0]) ) * _msPerX);
                    else
                        oParamCrete[jcrete][Su] = ((float)(meilleurk[0]-ymax)*_khzPerY)
                                /( ( (float)(pc[meilleurk[0]-ymin]-pc[ymax-ymin]) ) * _msPerX);

                }
                oParamCrete[jcrete][Sl] = 9999.0f;
                if(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]!=0)
                    oParamCrete[jcrete][Sl] = ((float)(lowfc-meilleurk[0])*_khzPerY)
                            /( ( (float)(pc[lowfc-ymin]-pc[meilleurk[0]-ymin]) ) * _msPerX);

                //
                int debp,finp,ecart;
                float pj[3];
                for(int j2=0;j2<3;j2++)
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
                        debp=ymcmax[jcrete]-ymin-2;
                        finp=ymcmax[jcrete]-ymin+2;
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
                oParamCrete[jcrete][StartSlope] = pj[0];
                oParamCrete[jcrete][EndSlope]   = pj[1];
                oParamCrete[jcrete][SlopeAtFc]  = pj[2];
            }
            oParamCrete[jcrete][Fdom]=(float)ymcmax[jcrete]*(float)_khzPerY;
            int intervTemps;
            if(jcrete<3) intervTemps = xmax-xmin+1;
            else intervTemps = tempsOuest[jcrete-3];
            oParamCrete[jcrete][Ldom] = (float)(xmcmax[jcrete]-xmin+0.5f)/(float)intervTemps;
            oParamCrete[jcrete][FdomDiff] = 9999.0f;
            if(criprec<icri) oParamCrete[jcrete][FdomDiff] = oParamCrete[jcrete][Fdom]
                    - _paramsArray[criprec][CM][Fdom];

        }
        oParamCrete[0][StartF] = (float)_yEmaxPerX[icri][0]*_khzPerY;
        oParamCrete[0][EndF] = (float)_yEmaxPerX[icri][xmax-xmin]*_khzPerY;
        if(_imageData)
        {
            _callMasterRidgeArray.push_back(_callMasterRidge);
            _callSouthArray.push_back(_callSouthRidge);
            _callNorthRidgeArray.push_back(_callNorthRidge);
            _callWestRidgeArray.push_back(_callWestRidge);
            _callSecondWestRidgeArray.push_back(_callSecondWestRidge);
        }
        oParam[Hup_RFMP] = 0.0f;
        oParam[Hup_PosMP] = 0.0f;
        oParam[Hup_PosSt] = 0.0f;
        oParam[Hup_PosEn] = 0.0f;
        oParam[Hup_RAmp] = 0.0f;
        oParam[Hup_RSlope] = 0.0f;
        oParam[Hlo_RFMP] = 0.0f;
        oParam[Hlo_PosMP] = 0.0f;
        oParam[Hlo_PosSt] = 0.0f;
        oParam[Hlo_PosEn] = 0.0f;
        oParam[Hlo_RAmp] = 0.0f;
        oParam[Hlo_RSlope] = 0.0f;
    }
    //_logText << endl << "-9" ;
    for(int icri=0;icri<nbcris-1;icri++)
    {
        float freqmp = _paramsArray[icri][SH][FreqMP];
        float fincri = _paramsArray[icri][SH][StartTime]+_paramsArray[icri][SH][Dur];
        bool faitpournextmp = false;
        bool aremplacer = false;
        for(int k=icri+1;k<nbcris;k++)
        {
            if(_paramsArray[k][SH][StartTime] < fincri) continue;
            if(!faitpournextmp)
            {
                _paramsArray[icri][SH][NextMP] =
                        (_masterPoints.at(k).x()-_masterPoints.at(icri).x())*_msPerX;
                faitpournextmp = true;
                aremplacer = true;
            }
            if(_paramsArray[k][SH][Fmin] < freqmp && _paramsArray[k][SH][Fmax] > freqmp)
            {
                _paramsArray[icri][SH][NextSmart1] =
                        (_masterPoints.at(k).x()-_masterPoints.at(icri).x())*_msPerX;
                break;
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
    int nhsup,nhinf;
    for(int icri=0;icri<nbcris;icri++)
    {
        nhsup=_harmonic[0][icri];
        if(nhsup>=0)
        {
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
            int xdeb = max(xmin1,xmin2);
            int xfin = min(xmax1,xmax2);
            float famp1=0.0f;
            float famp2=0.0f;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=_tabX[icri][i-xmin1];
                famp2+=_tabX[nhsup][i-xmin2];
            }
            if(famp1>0.0f) _paramsArray[icri][SH][Hup_RAmp] = famp2/famp1;
            int nouveauxdeb=xdeb;
            int nouveauxfin=xfin;
            for(int i=xdeb;i<=xfin;i++)
            {
                if(_yEmaxPerX[icri]<=0 || _yEmaxPerX[nhsup]<=0) nouveauxdeb++;
                else break;

            }
            for(int i=xfin;i>=nouveauxdeb;i++)
            {
                if(_yEmaxPerX[icri]<=0 || _yEmaxPerX[nhsup]<=0) nouveauxfin--;
                else break;
            }
            xdeb=nouveauxdeb; xfin=nouveauxfin;

            if(xfin>xdeb)
            {

                float CM_Slope1=((float)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((float)(xfin-xdeb));
                float CM_Slope2=((float)(_yEmaxPerX[nhsup][xfin-xmin2]
                                 -_yEmaxPerX[nhsup][xdeb-xmin2]))
                        /((float)(xfin-xdeb));
                if(CM_Slope1!=0.0f)
                    _paramsArray[icri][SH][Hup_RSlope] = CM_Slope2 / CM_Slope1;
            }
        }
        nhinf=_harmonic[1][icri];
        if(nhinf>=0)
        {
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
            int xdeb = max(xmin1,xmin2);
            int xfin = min(xmax1,xmax2);
            float famp1=0.0f;
            float famp2=0.0f;
            for(int i=xdeb;i<=xfin;i++)
            {
                famp1+=_tabX[icri][i-xmin1];
                famp2+=_tabX[nhinf][i-xmin2];
            }
            if(famp1>0.0f) _paramsArray[icri][SH][Hlo_RAmp] = famp2/famp1;
            int nouveauxdeb=xdeb;
            int nouveauxfin=xfin;
            for(int i=xdeb;i<=xfin;i++)
            {
                if(_yEmaxPerX[icri]<=0 || _yEmaxPerX[nhinf]<=0) nouveauxdeb++;
                else break;

            }
            for(int i=xfin;i>=nouveauxdeb;i++)
            {
                if(_yEmaxPerX[icri]<=0 || _yEmaxPerX[nhinf]<=0) nouveauxfin--;
                else break;
            }
            xdeb=nouveauxdeb; xfin=nouveauxfin;

            if(xfin>xdeb)
            {
                float CM_Slope1=((float)(_yEmaxPerX[icri][xfin-xmin1]
                                 -_yEmaxPerX[icri][xdeb-xmin1]))
                        /((float)(xfin-xdeb));
                float CM_Slope2=((float)(_yEmaxPerX[nhinf][xfin-xmin2]
                                 -_yEmaxPerX[nhinf][xdeb-xmin2]))
                        /((float)(xfin-xdeb));
                if(CM_Slope1!=0.0f)
                    _paramsArray[icri][SH][Hlo_RSlope] = CM_Slope2 / CM_Slope1;

            }
        }
    }
    _logText << "-dp2 fin" << endl;
}


void Detec::saveParameters(const QString& wavFile)
{
    _logText << "saveParameters  début  wavFile=" << wavFile << endl;
    int nbps = _vectPar.size();
    _txtFilePath2 = _txtPath+"/"+wavFile.left(wavFile.length()-3)+"csv";
    _txtFile2.setFileName(_txtFilePath2);
    if(_txtFile2.open(QIODevice::WriteOnly | QIODevice::Text)==false)
    {
        _logText << "ouverture en écriture du fichier " << _txtFilePath2 << " impossible !" << endl;
        return;
    }
    _logText << "saveParameters écriture dans fichier " << _txtFilePath2 << endl;
    _fileStream2.setDevice(&_txtFile2);
    _fileStream2.setRealNumberNotation(QTextStream::FixedNotation);
    _fileStream2.setRealNumberPrecision(6);
    _fileStream2 << "Filename"<< '\t' << "CallNum";
    for(int j=0;j<nbps;j++) _fileStream2 << '\t' << _vectPar[j].ColumnTitle;
    _fileStream2 << endl;
    for(int i=0;i<_callsNumber;i++)
    {
        float **callArray = _paramsArray[i];
        _fileStream2 << wavFile << '\t' << QString::number(i);
        for(int j=0;j<nbps;j++)
            _fileStream2 << '\t' <<  callArray[_vectPar[j].NumTableau][_vectPar[j].NumPar] ;
        _fileStream2 << endl;
    }
    _txtFile2.close();
    _callsArray.clear();
    _vectorXMin.clear();
    _vectorXMax.clear();
    _masterPoints.clear();
    _logText << "saveParameters  fin" << endl;
}

