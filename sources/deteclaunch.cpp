#include "deteclaunch.h"
#include "detec.h"
using namespace std;
#include <iostream>

DetecLaunch::DetecLaunch(QObject *parent) : QObject(parent)
{
	// Deteclaunch : main class of tadaridaD (non graphic class)
}

DetecLaunch::~DetecLaunch()
{

}

bool DetecLaunch::Treat(int argc, char *argv[])
{
    // Treat : this method manages the whole treatment
	// Initializations
    _timeExpansion = 10;
    _nbThreads = 1;
    _withTimeCsv = false;
    _paramVersion = 1;  
    IDebug = false;
    _mustCompress = false;
    _modeFreq = 1;
    _launchRecord = false; 
	_wavStock = false;
    _helpShown = false;
    _audioName = "";
    _modeDirFile = NODETERMINED;
    _recordSize = 5;
    _nRecords = 10;
    _nFilesPerTreatment=50;
    int compteur = 0;
    Detec **pdetec;
    QStringList   *pWavFileList;
    bool *threadRunning;
    bool *processRunning;
    QString logDirPath;
    _nbTreatedFiles = 0;
    _nbErrorFiles = 0;
    // Use of parameters
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
                if(alire.right(1)=="h") showHelp();
                if(alire.right(1) == "t") {paramParam = PARAMNTHREADS; waitValue=true;}
                if(alire.right(1) == "s") _withTimeCsv = true;
                if(alire.right(1) == "v") {paramParam = PARAMNVERSION; waitValue=true;}
                if(alire.right(1) == "d") IDebug = true;
                if(alire.right(1) == "r") {_launchRecord = true; paramParam = PARAMRECORD; waitValue=true;}
                if(alire.right(1) == "a") {paramParam = PARAMAUDIO; waitValue=true;}
                if(alire.right(1) == "w") _wavStock = true;
                if(alire.right(1) == "f") {paramParam = PARAMFREQ; waitValue=true;}
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
                if(paramParam==PARAMNVERSION)
                {
                    if(ok==true && value > 0) _paramVersion = value;
                }
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
            }
            else
            {
                if(!on_a_determine)
                {
                    determine = "";
                    if(alire.length()>4) determine = alire.right(4).toLower();
                    if(determine == ".wav" ) _modeDirFile = FILESMODE;
                    else
                    {
                        _modeDirFile = DIRECTORYMODE;
                        _wavPath = alire;
                    }
                }
                if(_modeDirFile == FILESMODE) firstList.append(alire);
            }
            waitValue=false;
        }
    }
    // -----------------------------------------------------------------
    // _launchRecord = true : tadaridad launches sox to record sound files to process
    if(_launchRecord)
    {
        _nbThreads = 1;
        _nSeries = (_nRecords + _nFilesPerTreatment - 1)/ _nFilesPerTreatment;
        _modeDirFile = DIRECTORYMODE;
    }
    else _nSeries = 0;
    // -----------------------------------------------------------------
    // initializations
    logDirPath = QDir::currentPath()+"/log";
    QDir logDir(logDirPath);
    if(!logDir.exists()) logDir.mkdir(logDirPath);
    QString logFilePath(logDirPath + QString("/tadaridaD.log"));
    _logFile.setFileName(logFilePath);
    _logFile.open(QIODevice::WriteOnly | QIODevice::Text);
    LogStream.setDevice(&_logFile);
    if(!_helpShown) showInfo("Launch TadaridaD",true);
    // the files to be processed are not entered: program stopped
    if(_modeDirFile == NODETERMINED)
    {
        if(!_helpShown)
        showInfo(QString("\nNo file or folder to treat ! You must enter :\n   >TadaridaD [optional setttings] [folder name or .wav files list to treat]\n"));
        _logFile.close();
        return(false);
    }
    // -----------------------------------------------------------------
    // Main loop - single pass in standard mode (without using sox)
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
            // creation of the list of files to be processed
            // First case: a directory has been entered
            if(_modeDirFile == DIRECTORYMODE)
            {
                QDir sdir(_wavPath);
                if(!sdir.exists())
                {
                    showInfo(QString("Folder ")+_wavPath+" does not exist !");
                    _logFile.close();
                    return(false);
                }
                _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
                QString dirPath = sdir.absolutePath();
                if(!createTxtFile(dirPath))
                {
                    showInfo("Unable to create txt directory !");
                    _logFile.close();
                    return(false);
                }
            }
            else
            {
                // Second case: file names have been entered
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
            }
            _nwf = _wavFileList.size();
            if(_nwf==0)
            {
                showInfo("No wav file to treat ");
                _logFile.close();
                return(false);
            }
            // -----------------------------------------------------------------
            // Distribution of files between threads
            if(_modeDirFile!=DIRECTORYMODE || _nwf <2) _nbThreads = 1;
            else
            {
                if(_nwf < _nbThreads) _nbThreads = _nwf;
            }
            showInfo(QString::number(_nbThreads)+" thread(s) launched");
            showInfo(QString("See log files in log folder"));
            pdetec = new Detec*[_nbThreads];
            pWavFileList = new QStringList[_nbThreads];
            threadRunning = new bool[_nbThreads];
            if(_nbThreads>0)
            {
                int c=0;
                for(int j=0;j<_nwf;j++)
                {
                    pWavFileList[c].append(_wavFileList.at(j));
                    c++;
                    if(c>=_nbThreads) c=0;
                }
            }
            // -----------------------------------------------------------------
            // Initializing shared fftw variables
            int fh;
            for(int j=0;j<_nbThreads;j++)
            {
                FftRes[j] 		= ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * FFT_HEIGHT_MAX );
                ComplexInput[j]        = ( fftwf_complex* ) fftwf_malloc( sizeof( fftwf_complex ) * FFT_HEIGHT_MAX );
                for(int i=0;i<6;i++)
                {
                    fh = pow(2,7+i);
                    Plan[j][i] = fftwf_plan_dft_1d(fh, ComplexInput[j], FftRes[j], FFTW_FORWARD, FFTW_ESTIMATE );
                }
            }
            // -----------------------------------------------------------------
            // Launching threads
            QString threadSuffixe = "";
            for(int i=0;i<_nbThreads;i++)
            {
                // creating a thread
                if(_nbThreads>1) threadSuffixe = QString("_") + QString::number(i+1);
                if(_nbThreads==1) pdetec[i] = new Detec(this,i,threadSuffixe,_modeDirFile,_wavPath,_wavFileList,_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug,_launchRecord,_mustCompress,_modeFreq);
                else  pdetec[i] = new Detec(this,i,threadSuffixe,_modeDirFile,_wavPath,pWavFileList[i],_wavRepList,_timeExpansion,_withTimeCsv,_paramVersion,IDebug,_launchRecord,_mustCompress,_modeFreq);
                connect(pdetec[i], SIGNAL(info1(QString,int)),this, SLOT(detecInfoTreat(QString,int)));
                // launching this thread
                pdetec[i]->start();
                threadRunning[i]=true;
            }
        }
        // -----------------------------------------------------------------
        // Mode of operation with sound recording: launch of sox
        if(_launchRecord)
        {
            if(iSerie < _nSeries)
            {
                int nfi = _nFilesPerTreatment;
                if(iSerie == _nSeries-1) nfi = _nRecords - (_nSeries - 1) * _nFilesPerTreatment;
                if(!lanceSox(iSerie,nfi,compteur))
                {
                    showInfo("Unable to launch sox !");
                    _logFile.close();
                    return(false);
                }
                compteur += nfi;

                if(iSerie==0) continue;
            }
        }
        // -----------------------------------------------------------------
        // Waiting for the end of threads
        int nbtr = _nbThreads;
        while(nbtr>0)
        {
            for(int i=0;i<_nbThreads;i++)
            {
                if(threadRunning[i]) if(!pdetec[i]->isRunning())
                {
                    threadRunning[i] = false;
                    nbtr--;
                    delete pdetec[i];
                    }
            }
            SLEEP(20);
        }
        delete[] pdetec;
        delete[] pWavFileList;
        delete[] threadRunning;
        // -----------------------------------------------------------------
        // Release of thread memory allocations
        for(int j=0;j<_nbThreads;j++)
        {
            fftwf_free(FftRes[j]);
            fftwf_free(ComplexInput[j]);
        }
        // -----------------------------------------------------------------
        if(_launchRecord && iSerie>0)
        {
            // Copy of files in the case of live recordings
            QString taStockPath(QDir::current().path()+"/txt");
            QDir taStock(taStockPath);
            if(!taStock.exists()) if(!taStock.mkdir(taStockPath))
            {
                showInfo("Unable to create txt directory !");
                _logFile.close();
                return(false);
            }
            QString taDirPath =  _wavPath+"/txt";
            QDir taDir(taDirPath);
            if(!taDir.exists()) return(false);
            QStringList taList = taDir.entryList(QStringList("*.ta"), QDir::Files);
            foreach(QString taFileName, taList)
            {
                QFile taFile(taDirPath + "/" + taFileName);
                if(taFile.exists()) taFile.copy(taStockPath+"/"+taFileName);
            }
            // conservation of .wav files (-w option)
            if(_wavStock)
            {
                QString wavStockPath(QDir::current().path()+"/wav");
                QDir wavStock(wavStockPath);
                if(!wavStock.exists()) if(!wavStock.mkdir(wavStockPath))
                {
                    showInfo("Unable to create wav directory !");
                    _logFile.close();
                    return(false);
                }
                // list re-created after cutting of files by thread
                QDir sdir(_wavPath);
                if(sdir.exists()) _wavFileList = sdir.entryList(QStringList("*.wav"), QDir::Files);
                foreach(QString wavFileName,_wavFileList)
                {
                    QFile wavFile(_wavPath + "/" + wavFileName);
                    if(wavFile.exists()) wavFile.copy(wavStockPath+"/"+wavFileName);
                }
            }
        }
    }
    // -----------------------------------------------------------------
    showInfo(QString("\nEnd of TadaridaD"),true,true,false);
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

bool DetecLaunch::lanceSox(int iserie,int nfi,int compteur)
{
    // launching sox to record the sound files to be processed
    QString wavtravPath(QDir::current().path()+"/wavtrav");
    if((iserie&1)==0) wavtravPath += "1"; else  wavtravPath += "2";
    QDir wavtrav(wavtravPath);
    if(!wavtrav.exists())
    {
        if(!wavtrav.mkdir(wavtravPath)) return(false);
    }
    else
    {
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

		}
    _wavTrav = wavtrav;
    QString program = "sox";
    // a big file has to be cut off in the thread
    QString fichierWav = wavtravPath + "/f" + QString::number(compteur+1) + ".wav";
    QStringList  arguments;
    // arguments << "-c" << "1" << "-d" << fichierWav << "trim" << "0" << QString::number(_recordSize*nfi) ;
    QString paraudio = "hw:1,0";
    if(!_audioName.isEmpty()) paraudio = QString("hw:") + _audioName;
    arguments << "-c" << "1" << "-t" << "alsa" << paraudio << fichierWav << "trim" << "0" << QString::number(_recordSize*nfi) ;
    QProcess p;
    p.execute(program,arguments);
    return(true);
}

void DetecLaunch::showInfo(QString s,bool showTime,bool e,bool l)
{
    QString s2 = "";
    if(showTime) s2 = QString("  -   ")+QDateTime::currentDateTime().toString("hh:mm:ss:zzz");
    if(l) LogStream << s << s2 << endl;
    if(e) cout << s.toStdString() << s2.toStdString() << endl;
}

void DetecLaunch::detecInfoTreat(QString wavFile,int resu)
{
    //LogStream << "detecInfoTreat " << endl;
    if(resu)_nbTreatedFiles++; else _nbErrorFiles++;
    QString sresu;
    if(resu==1) sresu = "ok"; else sresu = "error";
    showInfo(QString("Treatment of ")+wavFile+" : "+sresu,true);
}

void DetecLaunch::showHelp()
{
    QString helpInfo = "\n-t [n] allows to execute n parallel threads (1 by default)\n";
    helpInfo += "-x [n] is the time expansion factor,\n   either 10 (default) for 10-times expanded .wav files\n   or 1 for direct recordings\n";
    helpInfo += "-v [n] sets the list of features to be extracted on each detected sound event\n   (2 by default)\n";
    helpInfo += "_f [n] sets the frequency bands to be used;\n   n = 2 allows to treat low frequencies (0.8 to 25 kHz)\n   whereas n=1 (default) treats high frequencies (8 to 250 kHz)\n";
    helpInfo += "-c gives compressed version of .ta output files\n";
    helpInfo += "\?nAfter optional settings, must be mentioned :\n   - either a directory path containing .wav files \n   - or a list of .wav files, to be processed.\n   Relative or absolute paths can be used.";
    showInfo(helpInfo,false,true,false);
    _helpShown = true;
}

