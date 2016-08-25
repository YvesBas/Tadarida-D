#include "deteclaunch.h"
#include "detec.h"


DetecLaunch::DetecLaunch(QObject *parent) : QObject(parent)
{
}

DetecLaunch::~DetecLaunch()
{
}

bool DetecLaunch::Treat(int argc, char *argv[])
{
    _timeExpansion = 10;
    _nbThreads = 1;
    _nbProcess = 1;
    _nCalled = 0;
    _withTimeCsv = false;
    _paramVersion = 2;
    IDebug = false;
    _mustCompress = false;
    _modeFreq = 1;
    // µµµµµ debut
    _launchRecord = false; _wavStock = false;
    _audioName = "";
    _modeDirFile = NODETERMINED;
    _recordSize = 5;
    _nRecords = 10;
    _nFilesPerTreatment=50;
    int compteur = 0;
    Detec **pdetec;
    QStringList   *pWavFileList;
    bool *threadRunning;
    QProcess **pProcess;
    bool *processRunning;
    QString processSuffixe;
    QString logDirPath;
    // µµµµµ fin
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
                if(alire.right(1) == "c") _mustCompress = true;
                // if(alire.right(1)=="h") paramParam = PARAMHELP; // TODO...
                if(alire.right(1) == "t") {paramParam = PARAMNTHREADS; waitValue=true;}
                if(alire.right(1) == "p") {paramParam = PARAMNPROCESS; waitValue=true;}
                if(alire.right(1) == "s") _withTimeCsv = true;
                if(alire.right(1) == "v") {paramParam = PARAMNVERSION; waitValue=true;}
                if(alire.right(1) == "d") IDebug = true;
                // µµµµµ debut

                if(alire.right(1) == "r") {_launchRecord = true; paramParam = PARAMRECORD; waitValue=true;}
                if(alire.right(1) == "a") {paramParam = PARAMAUDIO; waitValue=true;}

                if(alire.right(1) == "w") _wavStock = true;
                if(alire.right(1) == "f") {paramParam = PARAMFREQ; waitValue=true;}
                // µµµµµ fin
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
                    if(ok==true && (value > 1 && value <= 8)) _nbThreads = value;
                }
                if(paramParam==PARAMNPROCESS)
                {
                    if(ok==true && (value > 1 && value <= 4)) _nbProcess = value;
                }
                if(paramParam==PARAMNCALLED)
                {
                    if(ok==true && value > 0) _nCalled = value;
                }
                if(paramParam==PARAMNVERSION)
                {
                    if(ok==true && value > 0) _paramVersion = value;
                }
                // µµµµµ debut
                if(paramParam==PARAMRECORD)
                {
                    if(ok==true && value > 0) _nRecords = value;
                }
                if(paramParam==PARAMAUDIO)
                {
                    _audioName = alire;
                }
                if(paramParam==PARAMFREQ)
                {
                    if(ok==true && (value ==1 || value == 2)) _modeFreq = value;
                }
                // µµµµµ fin
            }
            else
            {
                // QString firstArgument = QString(argv[1]);
                if(!on_a_determine)
                {
                    // µµµµµ debut
                    determine = "";
                    // µµµµµ fin
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
    // µµµµµ debut
    if(_launchRecord)
    {
        // important : empêcher le multi-process si lancement de sox
        _nbProcess = 1;
        _nbThreads = 1; // nécessaire dans la version en cours en raison de la decoupe...
        _nSeries = (_nRecords + _nFilesPerTreatment - 1)/ _nFilesPerTreatment;
        _modeDirFile = DIRECTORYMODE;
    }
    else _nSeries = 0;
    // µµµµµ fin
    // -----------------------------------------------------------------
    // 2) initialisations
    if(_modeDirFile == NODETERMINED)
    {
        return(false);
    }
// -----------------------------------------------------------------
    // µµµµµ debut
    logDirPath = QDir::currentPath()+"/log";
    QDir logDir(logDirPath);
    if(!logDir.exists()) logDir.mkdir(logDirPath);
    if(_nCalled==0 && _nbProcess==1)
    {
        QString logFilePath(logDirPath + QString("/tadaridaD.log"));
        _logFile.setFileName(logFilePath);
        _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
        LogStream.setDevice(&_logFile);
        LogStream << "Lancement de TadaridaD - " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    }

for(int iSerie = 0;iSerie<_nSeries+1;iSerie++)
{
    if(_launchRecord && iSerie>0)
    {
        QString wavtravPath(QDir::current().path()+"/wavtrav");
        if((iSerie&1)==1) wavtravPath += "1"; else  wavtravPath += "2";
        _wavPath = wavtravPath;
    }

    if(!_launchRecord || iSerie>0)
    {
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
            if(_launchRecord)
            {
                LogStream << "iSerie = " << iSerie << " listsize=" << _wavFileList.size()  << " _wavPath=" << _wavPath
                         << " " <<  QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            }
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
        // µµµµµ debut
        // mise en commentaire
        // QProcess **pProcess;
        // bool *processRunning;
        processSuffixe = "";
        // µµµµµ debut
        // -------------------------
        // µµµµµ debut
        // logDirPath = QDir::currentPath()+"/log";
        // QDir logDir(logDirPath);
        //if(!logDir.exists()) logDir.mkdir(logDirPath);
        // -------------------------
        if(_nCalled>0 || _nbProcess > 1)
        {
            processSuffixe = QString::number(_nCalled+1);
            QString logFilePath(logDirPath + QString("/tadaridaD")+processSuffixe+QString(".log"));
            _logFile.setFileName(logFilePath);
            _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
            LogStream.setDevice(&_logFile);
            LogStream << "_nCalled=" << _nCalled << endl;
            LogStream << "_nbProcess=" << _nbProcess << endl;
            LogStream << "Lancement TadaridaD - " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        }
        // µµµµµ fin
        // -----------------------------------------------------------------
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
                LogStream << "Avant lancement du processus " << l+1 << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                pProcess[l]->start(program,arguments);
                LogStream << "Apres lancement du processus " << l+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                processRunning[l] = true;
                QString endFilePath(logDirPath + QString("/end")+QString::number(l)+QString(".log"));
                QFile endFile;
                endFile.setFileName(endFilePath);
                if(endFile.exists())
                {
                    LogStream << "Effacement initial du fichier " << endFilePath  << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                    endFile.remove();
                }
            }
            LogStream << "Fin du demarrage des autres processus "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

        }
        // -----------------------------------------------------------------
        // 5) Répartition entre les threads
        if(_modeDirFile!=DIRECTORYMODE || _nwfp <2) _nbThreads = 1;
        else
        {
            if(_nwfp < _nbThreads) _nbThreads = _nwfp;
        }
        LogStream << "Nbthreads = " << _nbThreads << endl;
        // µµµµµ debut
        // variables initialisées au début
        pdetec = new Detec*[_nbThreads];
        pWavFileList = new QStringList[_nbThreads];
        threadRunning = new bool[_nbThreads];
        // µµµµµ fin
        //
        // µµµµµ debut
        //
        if(_nbThreads>0)
        // µµµµµ fin
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
        // 6) Initialisation des variables fftw partagées
        //_plan = fftwf_plan_dft_1d( _fftHeight, _complexInput, _fftRes, FFTW_FORWARD, FFTW_ESTIMATE );
        int fh;
        for(int j=0;j<_nbThreads;j++)
        {
            FftRes[j] 		= ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * FFT_HEIGHT_MAX );
            ComplexInput[j]        = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * FFT_HEIGHT_MAX );
            for(int i=0;i<6;i++)
            {
                fh = pow(2,7+i);
                // _logText << "plan " << i << " = " << fh << endl;
                Plan[j][i] = fftwf_plan_dft_1d(fh, ComplexInput[j], FftRes[j], FFTW_FORWARD, FFTW_ESTIMATE );
            }
        }
        // -----------------------------------------------------------------
        // 7) Lancement des threads
        QString threadSuffixe = "";
        for(int i=0;i<_nbThreads;i++)
        {
            LogStream << "Creation du thread " << i+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            if(_nbThreads>1) threadSuffixe = QString("_") + QString::number(i+1);
            if(_nbThreads==1) pdetec[i] = new Detec(this,processSuffixe,i,threadSuffixe,_modeDirFile,_wavPath,_wavFileListProcess,_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug,_launchRecord,_mustCompress,_modeFreq);
            else  pdetec[i] = new Detec(this,processSuffixe,i,threadSuffixe,_modeDirFile,_wavPath,pWavFileList[i],_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug,_launchRecord,_mustCompress,_modeFreq);

        // variables à initialiser
            LogStream << "Lancement du thread " << i+1 << " " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            pdetec[i]->start();
            threadRunning[i]=true;
            //SLEEP(500); mise en rem pour version 11.1
        }
    // µµµµµ debut
    } // condition _launchRecord == false || iSerie>0
    // µµµµµ fin
    // -----------------------------------------------------------------
    // 7,5) Lancement de sox
    // µµµµµ debut
    if(_launchRecord)
    {
        if(iSerie < _nSeries)
        {
            int nfi = _nFilesPerTreatment;
            if(iSerie == _nSeries-1) nfi = _nRecords - (_nSeries - 1) * _nFilesPerTreatment;
            LogStream << "Avant lancement de sox - nfi = " << nfi << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            if(!lanceSox(iSerie,nfi,compteur)) return(false);
            compteur += nfi;
            LogStream << "Retour de sox - compteur = " << compteur << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            if(iSerie==0) continue;
        }
    }
    // µµµµµ fin
    // -----------------------------------------------------------------
    // 8) Boucle d'attente de fin des threads lancés par le processus en cours
    int nbtr = _nbThreads;
    while(nbtr>0)
    {
        for(int i=0;i<_nbThreads;i++)
        {
            if(threadRunning[i]) if(!pdetec[i]->isRunning())
            {
                threadRunning[i] = false;
                nbtr--;
                LogStream << "Détecté fin du thread " << i+1 << "delete  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

                delete pdetec[i];
                LogStream << "Détecté fin du thread " << i+1 << " après delete  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
            }
        }
        SLEEP(20); // descendu de 100 à 20 pour version 11.1
    }
    LogStream << "Delete des tableaux pdetec etc...  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    delete[] pdetec;
    delete[] pWavFileList;
    delete[] threadRunning;
    LogStream << "Après delete des tableaux  " << endl;
    // -----------------------------------------------------------------
    // 9) libération des variables des threads
    for(int j=0;j<_nbThreads;j++)
    {
        fftwf_free(FftRes[j]);
        fftwf_free(ComplexInput[j]);
    }
    // -----------------------------------------------------------------
    // 10) boucle d'attente de fin des autres processus
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
                        LogStream << "Détecté fin du processus " << i+1 << "  -  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
                        processRunning[i] = false;
                        _nbPec--;
                    }
                }
            }
            SLEEP(200);
        }
        LogStream << "Avant delete des processus.  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        for(int i=1;i<_nbProcess;i++)
        {
            LogStream << "Delete du processus " << i+1 << "  : "<< QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

            delete pProcess[i];
            LogStream << "Apres delete du processus " << i+1  << endl;
        }
        delete[] pProcess;
        LogStream << "Après delete du tableau des processus  " << endl;
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
    // µµµµµ debut
    if(_launchRecord && iSerie>0)
    {
        // TODO : récupérer les résultats de wavtrav... dans répertoire
        LogStream << "Avant copie de txt "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        QString taStockPath(QDir::current().path()+"/txt");
        QDir taStock(taStockPath);
        if(!taStock.exists()) if(!taStock.mkdir(taStockPath)) return(false);
        QString taDirPath =  _wavPath+"/txt";
        QDir taDir(taDirPath);
        if(!taDir.exists()) return(false);
        QStringList taList = taDir.entryList(QStringList("*.ta"), QDir::Files);
        foreach(QString taFileName, taList)
        {
            QFile taFile(taDirPath + "/" + taFileName);
            if(taFile.exists()) taFile.copy(taStockPath+"/"+taFileName);
        }
        LogStream << "Apres copie de txt "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        if(_wavStock)
        {
            QString wavStockPath(QDir::current().path()+"/wav");
            QDir wavStock(wavStockPath);
            if(!wavStock.exists()) if(!wavStock.mkdir(wavStockPath)) return(false);
            // on refait _wavFileList car découpe faite dans detec :
            QDir sdir(_wavPath);
            if(sdir.exists()) _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
            foreach(QString wavFileName,_wavFileList)
            {
                QFile wavFile(_wavPath + "/" + wavFileName);
                if(wavFile.exists()) wavFile.copy(wavStockPath+"/"+wavFileName);
            }
        }
    }
    // µµµµµ fin
    // µµµµµ debut
// fin de boucle iSerie
}
// µµµµµ fin
    // -----------------------------------------------------------------
    LogStream << "Fin d'exécution TadaridaD" << processSuffixe << "  -  " << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
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
    if(es==0) LogStream << "Fin de processus normale _ exitCode = " << ec << endl;
    else LogStream << "Fin de processus avec crash" << endl;
}


void DetecLaunch::processStarted()
{
    LogStream << "évt processStarted" << endl;
}

void DetecLaunch::processError(QProcess::ProcessError pe)
{
    int ipe = (int)pe;
    LogStream << "évt processError erreur = " << ipe << endl;
}

// µµµµµ debut
bool DetecLaunch::lanceSox(int iserie,int nfi,int compteur)
{
    QString wavtravPath(QDir::current().path()+"/wavtrav");
    if((iserie&1)==0) wavtravPath += "1"; else  wavtravPath += "2";
    QDir wavtrav(wavtravPath);
    if(!wavtrav.exists())
    {
        if(!wavtrav.mkdir(wavtravPath)) return(false);
    }
    else
    {
        // vider le répertoire
        LogStream << "Avant vidage de repertoire"  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        QStringList listBefore = wavtrav.entryList(QStringList("*.wav"), QDir::Files);
        foreach(QString f,listBefore) wavtrav.remove(f);
        QString taTravPath(wavtravPath+"/txt");
        QDir taTrav(taTravPath);
        if(taTrav.exists())
        {
            QString taDirPath =  wavtravPath+"/txt";
            QDir taDir(taDirPath);
            QStringList listTaBefore = taDir.entryList(QStringList("*.ta"), QDir::Files);
            foreach(QString f,listTaBefore) taDir.remove(f);
        }
        LogStream << "Apres vidage de repertoire"  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    }
    _wavTrav = wavtrav;
    QString program = "sox";
    /*
    for(int j=0;j<nfi;j++)
    {
        QString fichierWav = wavtravPath + "/f" + QString::number(compteur+j+1) + ".wav";
        QStringList  arguments;
        arguments << "-c" << "1" << "-d" << fichierWav << "trim" << "0" << QString::number(_recordSize) ;
        _logText << "Lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
        QProcess p;
        p.execute(program,arguments);
        // SLEEP(_recordSize * 1000);
        _logText << "Retour du lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    }
    */
    /*
    // remplacé par :
    // version 11.0
    QString fichierWav = wavtravPath + "/f" + QString::number(compteur+1) + ".wav";
    QStringList  arguments;
    arguments << "-c" << "1" << "-d" << fichierWav << "trim" << "0" << QString::number(_recordSize) ;
    for(int j=1;j<nfi;j++) arguments << ":" << "newfile" <<  ":"  << "trim" << "0" << QString::number(_recordSize) ;
    _logText << "Lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QProcess p;
    p.execute(program,arguments);
    // SLEEP(_recordSize * 1000);
    _logText << "Retour du lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    */
    // version 11.1 : un gros fichier à découper dans le thread !
    QString fichierWav = wavtravPath + "/f" + QString::number(compteur+1) + ".wav";
    QStringList  arguments;
    // arguments << "-c" << "1" << "-d" << fichierWav << "trim" << "0" << QString::number(_recordSize*nfi) ;
    QString paraudio = "hw:1,0";
    if(!_audioName.isEmpty()) paraudio = QString("hw:") + _audioName;
    arguments << "-c" << "1" << "-t" << "alsa" << paraudio << fichierWav << "trim" << "0" << QString::number(_recordSize*nfi) ;
    LogStream << "Lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;
    QProcess p;
    p.execute(program,arguments);
    // SLEEP(_recordSize * 1000);
    LogStream << "Retour du lancement reel de sox - fichierwav = " << fichierWav << " "  << QDateTime::currentDateTime().toString("hh:mm:ss:zzz") << endl;

  return(true);
}
// µµµµµ fin
