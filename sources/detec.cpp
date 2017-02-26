#include "detec.h"
#include "detectreatment.h"

using namespace std;


Detec::Detec(DetecLaunch *pdl,int iThread,QString threadSuffixe,int modeDirFile,QString wavPath,QStringList wavFileList,QStringList wavRepList,int timeExpansion,bool withTimeCsv,int parVer,bool iDebug,bool withSox,bool mustCompress,int modeFreq): QThread((QObject *)pdl)
{
    // detec class : deteclaunch treatment launches detec objects as threads which share sound files to treat
    // initialization of variables :
    PMainWindow = pdl;
    IThread = iThread;
    _threadSuffixe = threadSuffixe;
    _modeDirFile = modeDirFile;
    _wavPath = wavPath;
    _wavFileList = wavFileList;
    _wavRepList = wavRepList;
    _timeExpansion = timeExpansion;
    _withTimeCsv = withTimeCsv;
    _paramVersion = parVer;
    IDebug = iDebug;
    // this variable means case of live recordings with sox
    _withSox = withSox;
    MustCompress = mustCompress;
    _modeFreq =  modeFreq;
    ReprocessingMode = false;
    _txtPath = _wavPath + "/txt";
    _resultSuffix = QString("ta");
    _detectionThreshold = 26;
    _stopThreshold = 20;
    _freqMin = 0;
    _nbo = 4;
    _useValflag = true;
    _jumpThreshold = 30;
    _widthBigControl = 60;
    _widthLittleControl = 5;
    _highThreshold = 10;
    _lowThreshold = 8;
    _highThreshold2 = 0;
    _lowThreshold2 = 10;
    _qR = 5;
    _qN = 5;
    _freqCallMin=8.0f;
    //  call of initializeDetec : initialization of log files
    initializeDetec();
    // creation of an object of the detectreatment class dedicated to the processing of sound files
    _detecTreatment = new DetecTreatment(this);
    LogStream << "LINWIN =  " << LINWIN << endl;
    LogStream << "VERQT =  " << VERQT << endl;
    LogStream << "_timeExpansion = " << _timeExpansion << endl;
    bool desactiveCorrectNoise = false;
    // assignment of settings for processing files
    _detecTreatment->SetGlobalParameters(modeFreq,_timeExpansion,_timeExpansion,_detectionThreshold,_stopThreshold,
                                         _freqMin,_nbo,
                                         _useValflag,_jumpThreshold,_widthBigControl,_widthLittleControl,
                                         _highThreshold,_lowThreshold,_highThreshold2,_lowThreshold2,_qR,_qN,_paramVersion,desactiveCorrectNoise);
    if(_modeDirFile==DIRECTORYMODE) _detecTreatment->SetDirParameters(_wavPath,_txtPath,false,"","");
    LogStream << "_wathPath="   << _wavPath << "     _txtPath="   << _txtPath << endl;
    _detecTreatment->InitializeDetecTreatment();
}

Detec::~Detec()
{
    delete _detecTreatment;
}

bool Detec::initializeDetec()
{
    // initialization of log files
    QString logDirPath = QDir::currentPath()+"/log";
    QString logFilePath(logDirPath + QString("/detec")+_threadSuffixe+".log");
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    LogStream.setDevice(&_logFile);
    LogStream.setRealNumberNotation(QTextStream::FixedNotation);
    LogStream.setRealNumberPrecision(22);

    LogStream << "Lancement Detec" <<_threadSuffixe << " : " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QString errorFilePath(logDirPath + QString("/error")+_threadSuffixe+".log");
    _errorFile.setFileName(errorFilePath);
    if(_errorFile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        ErrorStream.setDevice(&_errorFile);
        ErrorFileOpen = true;
    }
    else ErrorFileOpen = false;
    //
    TimeFileOpen = false;
    if(_withTimeCsv)
    {
        QString timePath(logDirPath + QString("/time")+_threadSuffixe+".csv");
        _timeFile.setFileName(timePath);
        if(_timeFile.open(QIODevice::WriteOnly | QIODevice::Text)==true)
        {
            TimeFileOpen = true;
            TimeStream.setDevice(&_timeFile);
            TimeStream.setRealNumberNotation(QTextStream::FixedNotation);
            TimeStream.setRealNumberPrecision(2);
            TimeStream << "filename" << '\t' << "computefft" << '\t' << "noisetreat" << '\t' << "shapesdetects" << '\t' << "parameters" << '\t'
                        << "save - end" << '\t' << "total time(ms)" << endl;
        }
    }
    return true;
}

void Detec::endDetec()
{
    if(ErrorFileOpen) _errorFile.close();
    if(TimeFileOpen) _timeFile.close();
    _detecTreatment->EndDetecTreatment();
    LogStream << "End of treatments by Detec object" << _threadSuffixe << " : " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    _logFile.close();
}

void Detec::run()
{
    // the run method manages the whole treatment by the thread
    QString dirPath,wavFile;
    if(_withSox) wavCut();
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
    LogStream << "Begin:" << wavFile << " : " << " r:" << dirPath
                   << "   -   " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QString pathFile = dirPath + '/' + wavFile;
    emit info1(wavFile,1);
    if(_modeDirFile == FILESMODE)
    {
        if(!createTxtDir(dirPath)) return;
        _detecTreatment->SetDirParameters(dirPath,_txtPath,false,"","");
    }
    bool resu = _detecTreatment->CallTreatmentsForOneFile(wavFile,pathFile);
    emit info1(wavFile,(int)resu);
    SLEEP(50);
    LogStream << "End:"<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
}

bool Detec::createTxtDir(QString dirPath)
{
    // Creation of txt directory
    _txtPath = dirPath+"/txt";
    QDir reptxt(_txtPath);
    if(!reptxt.exists())
    {
        if(!reptxt.mkdir(_txtPath))
        {
            LogStream << "cannot create directory:"<< _txtPath << " !" << endl;
            return(false);
        }
    }
    return(true);
}

void Detec::wavCut()
{
    // method used in the case of live recordings
    // this method cuts the file of sounds recorded by sox
    // into files of duration of 5 seconds maximum
    QString adecouper = _wavFileList.at(0);
    QString grosFichierWav = _wavPath + "/" + adecouper;
    QString resuWav = _wavPath + "/r" + adecouper;
    QString program = "sox";
    QStringList  arguments;
    arguments << grosFichierWav << resuWav << "trim" << "0" << "5"  << ":" << "newfile" <<  "restart";
    LogStream << "decoupe par sox de fichierwav = " << grosFichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QProcess p;
    p.execute(program,arguments);
    QDir sdir(_wavPath);
    if(!sdir.exists()) return;
    sdir.remove(adecouper);
    _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
}

