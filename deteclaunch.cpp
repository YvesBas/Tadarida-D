#include "deteclaunch.h"
#include "detec.h"


DetecLaunch::DetecLaunch(QObject *parent) : QObject(parent)
{
}

DetecLaunch::~DetecLaunch()
{
}

bool DetecLaunch::treat(int argc, char *argv[])
{
    _timeExpansion = 10;
    _nbThreads = 1;
    _nbProcess = 1;
    _nCalled = 0;
    _withTimeCsv = false;
    _paramVersion = 0;
    IDebug = false;
    // -----------------------------------------------------------------
    // 1) lecture des paramètres
    bool waitValue=false;
    int paramParam=PARAMNO;
    bool on_a_determine = false;
    QStringList   firstList;
    QString determine;
    for(int i=1;i<argc;i++)
    {
        QString alire = QString(argv[i]);
        if(alire.left(1)=="-")
        {
            paramParam = PARAMNO;
            waitValue = false;
            if(alire.length()==2)
            {
                if(alire.right(1) == "x") {paramParam = PARAMEXPANSION; waitValue=true;}
                if(alire.right(1) == "c") paramParam = PARAMCOMPRESS; // TODO...
                // if(alire.right(1)=="h") paramParam = PARAMHELP; // TODO...
                if(alire.right(1) == "t") {paramParam = PARAMNTHREADS; waitValue=true;}
                if(alire.right(1) == "p") {paramParam = PARAMNPROCESS; waitValue=true;}
                if(alire.right(1) == "s") _withTimeCsv = true;
                if(alire.right(1) == "v") {paramParam = PARAMNVERSION; waitValue=true;}
                if(alire.right(1) == "d") IDebug = true;
            }
            else
            {
                if(alire=="-called")  {paramParam = PARAMNCALLED; waitValue=true;}
            }
        }
        else
        {
            if(waitValue)
            {
                bool ok;
                int value = alire.toInt(&ok,10) ;
                if(paramParam==PARAMEXPANSION)
                {
                    if(ok==true && (value == 1 || value == 10)) _timeExpansion = value;
                }
                if(paramParam==PARAMNTHREADS)
                {
                    if(ok==true && (value > 1 || value <= 8)) _nbThreads = value;
                }
                if(paramParam==PARAMNPROCESS)
                {
                    if(ok==true && (value > 1 || value <= 4)) _nbProcess = value;
                }
                if(paramParam==PARAMNCALLED)
                {
                    if(ok==true && value > 0) _nCalled = value;
                }
                if(paramParam==PARAMNVERSION)
                {
                    if(ok==true && value > 0) _paramVersion = value;
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
                if(_modeDirFile == FILESMODE) firstList.append(alire);
            } // fin du else waitvalue = true
            waitValue=false;
        } // fin du cas non paramètre "-"
    }
    // -----------------------------------------------------------------
    // 2) initialisations

    if(_modeDirFile == DIRECTORYMODE)
    {
        QDir sdir(_wavPath);
        if(!sdir.exists())
        {
            //logText << "répertoire "<< _wavPath << " non trouvé !" << endl;
            //logFile.close();
            return(false);
        }
        _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
        if(_wavFileList.isEmpty())
        {
            //logText << "aucun fichier wav trouvé dans le répertoire "<< _wavPath << " !" << endl;
            //return(false);
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
        _nwf=0;
        foreach(QString wf,firstList)
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
                        _nwf++;
                    }
                }
            }
        }
        if(_nwf==0)
        {
            //logText << "aucun fichier wav trouvé !" << endl;
            //logFile.close();
            return(false);
        }
    }
    _nwf = _wavFileList.size();
    // -----------------------------------------------------------------
    //  3) Mise à jour de la liste pour le process en cours
    if(_modeDirFile!=DIRECTORYMODE || _nwf < 3)  _nbProcess = 1;
    else
    {
        if(_nwf <= _nbProcess ) _nbProcess = _nwf-1;
    }
    _wavFileListProcess.clear();
    if(_nbProcess>1)
    {
        int c=0;
        for(int j=0;j<_nwf;j++)
        {
            if(c==_nCalled)  _wavFileListProcess.append(_wavFileList.at(j));
            c++;
            if(c>=_nbProcess) c=0;
        }
    }
    else
    {
        _wavFileListProcess = _wavFileList;
    }
    _nwfp = _wavFileListProcess.size();
    QProcess **pProcess;
    bool *processRunning;
    // -------------------------------    QString txtPath = dirPath+"/txt";
    QString logDirPath = QDir::currentPath()+"/log";
    QDir logDir(logDirPath);
    if(!logDir.exists()) logDir.mkdir(logDirPath);
    // ----------------------------------
    QString processSuffixe = "";
    if(_nCalled>0 || _nbProcess > 1) processSuffixe = QString::number(_nCalled+1);


    QString logFilePath(logDirPath + QString("/tadaridaD")+processSuffixe+QString(".log"));
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    _logText.setDevice(&_logFile);
    _logText << "Idebug=" << IDebug << endl;
    _logText << "Lancement TadaridaD - " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    // -----------------------------------------------------------------
    // 4) Démarrage des autres process si nécessaire
    if(_nbProcess>1 && _nCalled == 0)
    {
        //
        _nbPec = _nbProcess - 1;
        pProcess = new QProcess*[_nbProcess];
        processRunning = new bool[_nbProcess];
        QString program = "TadaridaD.exe";
        QStringList  argumentsBase;
        for(int k=1;k<argc;k++) argumentsBase << QString(argv[k]);
        // arguments << "a" << "-tgzip" << compressedParametersPath <<  txtFilePath;
        for(int l=1;l<_nbProcess;l++)
        {
            pProcess[l] = new QProcess(this);
            connect(pProcess[l], SIGNAL(started()),this, SLOT(processStarted()));
            connect(pProcess[l], SIGNAL(error(QProcess::ProcessError)),this, SLOT(processError(QProcess::ProcessError)));
            connect(pProcess[l], SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(processFinished(int,QProcess::ExitStatus)));
            QStringList arguments = argumentsBase;
            arguments << "-called" << QString::number(l);
            _logText << "Avant lancement du processus " << l+1 << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            pProcess[l]->start(program,arguments);
            _logText << "Apres lancement du processus " << l+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            processRunning[l] = true;
            QString endFilePath(logDirPath + QString("/end")+QString::number(l)+QString(".log"));
            QFile endFile;
            endFile.setFileName(endFilePath);
            if(endFile.exists())
            {
                _logText << "Effacement initial du fichier " << endFilePath  << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                endFile.remove();
            }
        }
        _logText << "Fin du demarrage des autres processus "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

    }
    // -----------------------------------------------------------------
    // 5) Répartition entre les threads
    if(_modeDirFile!=DIRECTORYMODE || _nwfp <2) _nbThreads = 1;
    else
    {
        if(_nwfp < _nbThreads) _nbThreads = _nwfp;
    }
    _logText << "Nbthreads = " << _nbThreads << endl;
    Detec **pdetec = new Detec*[_nbThreads];
    QStringList   *pWavFileList = new QStringList[_nbThreads];
    bool *threadRunning = new bool[_nbThreads];
    if(_nbThreads>1)
    {
        int c=0;
        for(int j=0;j<_nwfp;j++)
        {
            pWavFileList[c].append(_wavFileListProcess.at(j));
            c++;
            if(c>=_nbThreads) c=0;
        }
    }
    // -----------------------------------------------------------------
    // 5) Lancement des threads
    QString threadSuffixe = "";
    for(int i=0;i<_nbThreads;i++)
    {
        _logText << "Creation du thread " << i+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

        if(_nbThreads>1) threadSuffixe = QString("_") + QString::number(i+1);
        if(_nbThreads==1) pdetec[i] = new Detec(processSuffixe,threadSuffixe,_modeDirFile,_wavPath,_wavFileListProcess,_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug);
        else  pdetec[i] = new Detec(processSuffixe,threadSuffixe,_modeDirFile,_wavPath,pWavFileList[i],_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug);
    // variables à initialiser

        _logText << "Lancement du thread " << i+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        pdetec[i]->start();
        threadRunning[i]=true;
    }
    // -----------------------------------------------------------------
    // 6) Boucle d'attente de fin des threads lancés par le processus en cours
    int nbtr = _nbThreads;
    while(nbtr>0)
    {
        for(int i=0;i<_nbThreads;i++)
        {
            if(threadRunning[i]) if(!pdetec[i]->isRunning())
            {
                threadRunning[i] = false;
                nbtr--;
                _logText << "Détecté fin du thread " << i+1 << "delete  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

                delete pdetec[i];
                _logText << "Détecté fin du thread " << i+1 << " après delete  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            }
        }
        SLEEP(100);
    }

    _logText << "Delete des tableaux pdetec etc...  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    delete[] pdetec;
    delete[] pWavFileList;
    delete[] threadRunning;
    _logText << "Après delete des tableaux  " << endl;
    // -----------------------------------------------------------------
    // 7) boucle d'attente de fin des autres processus
    // TODO : donner un temps d'attente limite proportionnel au nombre de fichiers
    if(_nCalled == 0 && _nbProcess > 1 && _nbPec>0)
    {
        while(_nbPec>0)
        {
            for(int i=1;i<_nbProcess;i++)
            {
                if(processRunning[i])
                {
                    QString endFilePath(logDirPath + QString("/end")+QString::number(i)+QString(".log"));
                    QFile endFile;
                    endFile.setFileName(endFilePath);
                    if(endFile.exists())
                    {
                        _logText << "Détecté fin du processus " << i+1 << "  -  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                        processRunning[i] = false;
                        _nbPec--;
                    }
                }
            }
            SLEEP(200);
        }
        _logText << "Avant delete des processus.  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        for(int i=1;i<_nbProcess;i++)
        {
            _logText << "Delete du processus " << i+1 << "  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

            delete pProcess[i];
            _logText << "Apres delete du processus " << i+1  << endl;
        }
        delete[] pProcess;
        _logText << "Après delete du tableau des processus  " << endl;
    }

    if(_nCalled > 0)
    {
        QString endFilePath(logDirPath + QString("/end")+QString::number(_nCalled)+QString(".log"));
        QFile endFile;
        QTextStream endText;
        endFile.setFileName(endFilePath);
        endFile.open(QIODevice::WriteOnly | QIODevice::Text);
        endText.setDevice(&endFile);
        endText << "Fin du processus TadaridaD" << _nCalled << " : " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        endFile.close();
    }

    // -----------------------------------------------------------------
    _logText << "Fin d'exécution TadaridaD" << processSuffixe << "  -  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    _logFile.close();
    return(true);

}

bool DetecLaunch::createTxtFile(QString dirPath)
{
    QString txtPath = dirPath+"/txt";
    QDir reptxt(txtPath);
    if(!reptxt.exists())
    {
        if(!reptxt.mkdir(txtPath)) return(false);
    }
    return(true);
}

void DetecLaunch::processFinished(int ec,QProcess::ExitStatus es)
{
    _nbPec--;
    if(es==0) _logText << "Fin de processus normale _ exitCode = " << ec << endl;
    else _logText << "Fin de processus avec crash" << endl;
}


void DetecLaunch::processStarted()
{
    _logText << "évt processStarted" << endl;
}

void DetecLaunch::processError(QProcess::ProcessError pe)
{
    int ipe = (int)pe;
    _logText << "évt processError erreur = " << ipe << endl;
}

