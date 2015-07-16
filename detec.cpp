#include "detec.h"
#include "detectreatment.h"

using namespace std;


Detec::Detec(DetecLaunch *pdl,QString processSuffixe,int iThread,QString threadSuffixe,int modeDirFile,QString wavPath,QStringList wavFileList,QStringList wavRepList,int timeExpansion,bool withTimeCsv,int parVer,bool iDebug): QThread()
{
    PDL = pdl;
    IThread = iThread;
    _processSuffixe = processSuffixe;
    _threadSuffixe = threadSuffixe;
    _modeDirFile = modeDirFile;
    _wavPath = wavPath;
    _wavFileList = wavFileList;
    _wavRepList = wavRepList;
    _timeExpansion = timeExpansion;
    _withTimeCsv = withTimeCsv;
    _paramVersion = parVer;
    IDebug = iDebug;

    ReprocessingMode = false;
    //
    _txtPath = _wavPath + "/txt";
    ResultSuffix = QString("ta");
    MustCompress = false;
    _imageData = false;
    //
    _detectionThreshold = 26;
    _stopThreshold = 20;
    _freqMin = 0;
    _nbo = 4;
    _useValflag = true;
    _jumpThreshold = 30;
    _widthBigControl = 60;
    _widthLittleControl = 5;
    _highThreshold = 10;
    _lowThreshold = -4;
    _qR = 5;
    _qN = 5;
    _freqCallMin=8.0f;
    InitializeDetec();
    _logText << "idebug = " << IDebug << endl;
    _detecTreatment = new DetecTreatment(this);
    _logText << "_timeExpansion = " << _timeExpansion << endl;
    _detecTreatment->SetGlobalParameters(_timeExpansion,_detectionThreshold,_stopThreshold,
                                         _freqMin,_nbo,
                                         _useValflag,_jumpThreshold,_widthBigControl,_widthLittleControl,
                                         _highThreshold,_lowThreshold,_qR,_qN,_paramVersion);
    if(_modeDirFile==DIRECTORYMODE) _detecTreatment->SetDirParameters(_wavPath,_txtPath,false,"","");
    _logText << "cd _wathPath="   << _wavPath << "     _txtPath="   << _txtPath << endl;
    _detecTreatment->InitializeDetecTreatment();
}

Detec::~Detec()
{
    delete _detecTreatment;
}

bool Detec::InitializeDetec()
{
    QString logDirPath = QDir::currentPath()+"/log";
    QString logFilePath(logDirPath + QString("/detec")+_processSuffixe+_threadSuffixe+".log");
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    _logText.setDevice(&_logFile);
    _logText << "Lancement Detec" <<_processSuffixe<<_threadSuffixe << " : " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    //
    //
    QString errorFilePath(logDirPath + QString("/error")+_processSuffixe+_threadSuffixe+".log");
    _errorFile.setFileName(errorFilePath);
    if(_errorFile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        _errorStream.setDevice(&_errorFile);
        _errorFileOpen = true;
    }
    else _errorFileOpen = false;
    //
    _timeFileOpen = false;
    if(_withTimeCsv)
    {
        QString timePath(logDirPath + QString("/time")+_processSuffixe+_threadSuffixe+".csv");
        _timeFile.setFileName(timePath);
        if(_timeFile.open(QIODevice::WriteOnly | QIODevice::Text)==true)
        {
            _timeFileOpen = true;
            _timeStream.setDevice(&_timeFile);
            _timeStream.setRealNumberNotation(QTextStream::FixedNotation);
            _timeStream.setRealNumberPrecision(2);
            _timeStream << "filename" << '\t' << "computefft" << '\t' << "noisetreat" << '\t' << "shapesdetects" << '\t' << "parameters" << '\t'
                        << "save - end" << '\t' << "total time(ms)" << endl;
        }
    }
    //
    return true;
}

void Detec::endDetec()
{
    if(_errorFileOpen) _errorFile.close();
    if(_timeFileOpen) _timeFile.close();
    _detecTreatment->EndDetecTreatment();
    _logText << "Fin de traitement Detec" << _processSuffixe  << _threadSuffixe << " : " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    _logFile.close();
    // IsRunning = false;
}

void Detec::run()
{
    QString dirPath,wavFile;
    for(int i=0;i<_wavFileList.size();i++)
    {
        if(_modeDirFile == DIRECTORYMODE) dirPath=_wavPath;
        else dirPath = _wavRepList.at(i);
        treatOneFile(_wavFileList.at(i),dirPath);
    }
    endDetec();
}

void Detec::treatOneFile(QString wavFile,QString dirPath)
{
    _logText << "Deb:" << wavFile << " : " << " rep:" << dirPath
                   << "   -   " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QString pathFile = dirPath + '/' + wavFile;
    if(_modeDirFile == FILESMODE)
    {
        if(!createTxtFile(dirPath)) return;
        _detecTreatment->SetDirParameters(dirPath,_txtPath,false,"","");
    }
    _detecTreatment->CallTreatmentsForOneFile(wavFile,pathFile);
    _logText << "Fin:"<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
}

bool Detec::createTxtFile(QString dirPath)
{
    // à améliorer : ne pas chercher à recréer à chaque fois
    _txtPath = dirPath+"/txt";
    QDir reptxt(_txtPath);
    if(!reptxt.exists())
    {
        if(!reptxt.mkdir(_txtPath))
        {
            _logText << "création du sous-répertoire "<< _txtPath << " impossible !" << endl;
            return(false);
        }
    }
    return(true);
}



