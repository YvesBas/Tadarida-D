#include "detec.h"
#include "detectreatment.h"
#define PARAMNO 0
#define PARAMEXPANSION 1
#define PARAMCOMPRESS 2
#define PARAMHELP 3

using namespace std;


Detec::Detec(int argc,char **argv): QThread()
{
    // IsRunning = true;
    _timeExpansion=10;
    ResultSuffix = QString("ta");
    MustCompress = false;
    // lecture des paramètres
    bool waitValue=false;
    int paramParam=PARAMNO;
    bool on_a_determine = false;
    QString determine;
    for(int i=1;i<argc;i++)
    {
        QString alire = QString(argv[i]);
        if(alire.left(1)=="-")
        {
            waitValue = false;
            if(alire.length()==2)
            {
                paramParam = PARAMNO;
                if(alire.right(1)=="x") {paramParam = PARAMEXPANSION; waitValue=true;}
                if(alire.right(1)=="c") paramParam = PARAMCOMPRESS; // TODO...
                if(alire.right(1)=="h") paramParam = PARAMHELP; // TODO...
            }
        }
        else
        {
            if(waitValue)
            {
                if(paramParam==PARAMEXPANSION)
                {
                    bool ok;
                    int valueExpansion = alire.toInt(&ok,10) ;
                    if(ok==true && (valueExpansion == 1 || valueExpansion == 10)) _timeExpansion = valueExpansion;
                }
            }
            else
            {
                // QString firstArgument = QString(argv[1]);
                if(!on_a_determine)
                {
                    if(alire.length()>4) determine = alire.right(4).toLower();
                    if(determine == ".wav" ) _modeDirFile = FILESMODE;
                    else
                    {
                        _modeDirFile = DIRECTORYMODE;
                        _wavPath = alire;
                    }
                }
                if(_modeDirFile == FILESMODE) _firstList.append(alire);
            } // fin du else waitvalue = true
            waitValue=false;
        } // fin du cas non paramètre "-"
    }
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
    _initPassed = InitializeDetec();
    if(_initPassed)
    {
        _detecTreatment = new DetecTreatment(this);
        _logText << "_timeExpansion = " << _timeExpansion << endl;
        _detecTreatment->SetGlobalParameters(_timeExpansion,_detectionThreshold,_stopThreshold,
                                             _freqMin,_nbo,
                                             _useValflag,_jumpThreshold,_widthBigControl,_widthLittleControl,
                                             _highThreshold,_lowThreshold,_qR,_qN);
        if(_modeDirFile==DIRECTORYMODE) _detecTreatment->SetDirParameters(_wavPath,_txtPath,false,"","");
        _detecTreatment->InitializeDetecTreatment();
    }
}

/*
Detec::Detec(QString dirPath): QThread()
{
    _wavPath = dirPath;
    _timeExpansion=10;
    _detectionThreshold = 26;
    _stopThreshold = 18;
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
    _initPassed = InitializeDetec();
    if(_initPassed)
    {
        _detecTreatment = new DetecTreatment(this);
        _detecTreatment->SetGlobalParameters(_timeExpansion,_detectionThreshold,_stopThreshold,
                                             _freqMin,_nbo,
                                             _useValflag,_jumpThreshold,_widthBigControl,_widthLittleControl,
                                             _highThreshold,_lowThreshold,_qR,_qN);
        _detecTreatment->SetDirParameters(_wavPath,_wavFileList,_txtPath,false,"","");
        _detecTreatment->InitializeDetecTreatment();
    }
}
*/

Detec::~Detec()
{
}

bool Detec::InitializeDetec()
{
    //IsRunning = true;
    QString logFilePath("detec.log");
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    _logText.setDevice(&_logFile);
    _logText << "Lancement TadaridaD - " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

    QString errorFilePath("error.log");
    _errorFile.setFileName(errorFilePath);
    if(_errorFile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
    _errorStream.setDevice(&_errorFile);
    _errorFileOpen = true;
    }
    else _errorFileOpen = false;
    //
    if(_modeDirFile == DIRECTORYMODE)
    {
        QDir sdir(_wavPath);
        if(!sdir.exists())
        {
            _logText << "répertoire "<< _wavPath << " non trouvé !" << endl;
            return(false);
        }
        _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
        if(_wavFileList.isEmpty())
        {
            _logText << "aucun fichier wav trouvé dans le répertoire "<< _wavPath << " !" << endl;
            return(false);
        }
        QString dirPath = sdir.absolutePath();
        if(!createTxtFile(dirPath)) return(false);
    }
    else
    {
        // _modeDirFile == FILESMODE
        QFile f;
        QDir d;
        QString determine;
        int nwf=0;
        foreach(QString wf,_firstList)
        {
            if(wf.length()>4) determine = wf.right(4).toLower();
            else determine = "";
            if(determine == ".wav" )
            {
                f.setFileName(wf);
                if(f.exists())
                {
                    QFileInfo finf(wf);
                    d = finf.dir();
                    if(d.exists())
                    {
                        _wavFileList.append(finf.fileName());
                        _wavRepList.append(d.absolutePath());
                        nwf++;
                    }
                }
            }
        }
        if(nwf==0)
        {
            _logText << "aucun fichier wav trouvé !" << endl;
            return(false);
        }
        _txtPath = ""; // à déterminerpour chaque fichier
    }
    return true;
}

void Detec::endDetec()
{
    _errorFile.close();
    if(_initPassed)
    {
        _detecTreatment->EndDetecTreatment();
    }
    _logText << "Fin d'exécution TadaridaD - " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    _logFile.close();
    // IsRunning = false;
}

void Detec::run()
{
    if(_initPassed)
    {
        //foreach(QString wavFile,_wavFileList) treatOneFile(wavFile);
        QString dirPath,wavFile;
        for(int i=0;i<_wavFileList.size();i++)
        {
            if(_modeDirFile == DIRECTORYMODE) dirPath=_wavPath;
            else dirPath = _wavRepList.at(i);
            treatOneFile(_wavFileList.at(i),dirPath);
        }
    }
    endDetec();
}

void Detec::treatOneFile(QString wavFile,QString dirPath)
{
    _logText << "Début de traitement du fichier " << wavFile << " : " << " dans le répertoire " << dirPath
                   << "   -   " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QString pathFile = dirPath + '/' + wavFile;
    if(_modeDirFile == FILESMODE)
    {
        if(!createTxtFile(dirPath)) return;
        _detecTreatment->SetDirParameters(dirPath,_txtPath,false,"","");
    }
    _detecTreatment->CallTreatmentsForOneFile(wavFile,pathFile);
    _logText << "Fin de traitement d'un fichier : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
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



