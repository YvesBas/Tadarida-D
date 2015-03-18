#ifndef DETEC_H
#define DETEC_H


#include <float.h>
#include <iostream>
#include <QDate>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QSettings>
#include <QtCore/qmath.h>
#include <QtGui/QMainWindow>
#include <QTextStream>
#include <QTimer>
#include <QThread>

#include <QVector>
//#include <QImage>
#include "sndfile.h"
#include "fftw3.h"
#define	FFT_HEIGHT_MAX 2048
#define FREQ_MAX 250 // KHz
#define MAXLARCRI 10000
#define MAXHAUCRI 500
#define PI	3.14159265358979
#define SONOGRAM_WIDTH_MAX 65535
#define SONOGRAM_WIDTH_MIN 64

#define MAXCRI 2000
#define NCRETES 5


enum NUMTAB {SH,CM,CS,CN,CO,CO2,NBTAB};
enum NUMPAR {StartTime,Dur,Prev,Fmax,Fmin,BW,Fdom,FdomDiff,Ldom,Slope,HiFc,LowFc,Lhi,Su,Sl,StartF,
     EndF,StartSlope,EndSlope,SlopeAtFc,FreqCtr,FBak5dB,FFwd5dB,Bndw5dB,DurOf5dB,Fc3,FreqMP,
     PosMP,FreqDomSum,FreqDomMean,PosPeakSum,PosPeakMean,FreqPeakSum,FreqPeakMean,PrevMP,
     PrevSmart1,NextMP,NextSmart1,Amp1,Amp2,Amp3,Amp4,NoiseLeft,NoiseRight,NoiseDown,NoiseUp,
     CVAmp,Hup_RFMP,Hup_PosMP,Hup_PosSt,Hup_PosEn,Hup_RAmp,Hup_RSlope,Hlo_RFMP,Hlo_PosMP,
     Hlo_PosSt,Hlo_PosEn,Hlo_RAmp,Hlo_RSlope,Ramp_2_1,Ramp_3_1,Ramp_3_2,Ramp_1_2,Ramp_4_3,
     Ramp_2_3,RAN_2_1,RAN_3_1,RAN_3_2,RAN_1_2,RAN_4_3,RAN_2_3,HetX,HetY,RPM8,Stab,
     NBPAR};

class ParamToSave
{
public:
    ParamToSave(int,int,QString);
    ParamToSave();
    ~ParamToSave();

    int NumTableau;
    int NumPar;
    QString ColumnTitle;
    QString InfoLabel;
};


class Detec : public QThread
{
    Q_OBJECT
public:
    Detec(QString);
    ~Detec();
    void run();
    bool IsRunning;

private:

    // methods
    bool computeFFT(QString &);
    void correctNoise();
    void detectsParameter2();
    void endDetec();
    bool InitializeDetec();
    void initVectorParams();
    bool openWavFile(QString &wavFile);
    void saveParameters(const QString&);
    void shapesDetects();
    void sortWaves();
    void treatOneFile(QString);

    // attributes
    float                        *_averagePerX;
    QDir                         _baseDayDir;
    float                        _callEnergyMax ;
    int                          _callEnergyMaxIndex;
    QVector< QVector<QPoint> >   _callsArray;
    QVector<QPoint>              _callMasterRidge;
    QVector< QVector<QPoint> >   _callMasterRidgeArray;
    QFile                        _callMatrixFile;
    QString                      _callMatrixName;
    QDataStream                  _callMatrixStream;
    QVector<QPoint>              _callNorthRidge;
    QVector< QVector<QPoint> >   _callNorthRidgeArray;
    QVector<QPoint>              _callSecondWestRidge;
    QVector< QVector<QPoint> >   _callSecondWestRidgeArray;
    QVector<QPoint>              _callSouthRidge;
    QVector< QVector<QPoint> >   _callSouthArray;
    QVector<QPoint>              _callWestRidge;
    QVector< QVector<QPoint> >   _callWestRidgeArray;
    int                          _callsNumber;
    char *                       _charParamsArray;
    char *                       _charTabX;
    char *                       _charTabYX;
    char *                       _charYEmaxPerX;
    QTimer                       *_clock;
    float                        *_coeff;
    fftwf_complex*	             _complexInput;
    float*			             _data;
    QString                      _datPath;
    int                          _detectionThreshold;
    float                        **_dpm;
    int                          **_dypm;
    float                        *_eMaxPerX;
    double			             _energyMax;
    double		                 _energyMin;
    float                        _energyShapeThreshold;
    float                        _energyStopThreshold;
    QFile                        _errorFile;
    QTextStream                  _errorStream;
    int			                 _fftHeight;
    int			                 _fftHeightHalf;
    fftwf_complex*	             _fftRes;
    QFile                        _fileInfo;
    QTextStream                  _fileStream;
    QTextStream                  _fileStream2;
    QTextStream                  _fIStream;
    float                        _freqCallMin;
    int                          _freqMin;
    int                          **_harmonic;
    QImage  	                 _image;
    bool                         _imageData;
    QString                      _imagePath;
    int                          _imaWidth;
    int                          *_inflexion1;
    int                          *_inflexion3;
    bool                         _initPassed;
    int				             _iOverlapMoving;
    float                        _khzPerY;
    bool                         _littleParams;
    QFile                        _logFile;
    QTextStream                  _logText;
    int                          *_lowSlope;
    QVector< QPoint >            _masterPoints;
    int                          _maxCallWidth;
    int                          _maxCallHeight;
    int                          _maxY;
    int			                 _medianNoise;
    int                          _minY;
    float                        _msPerX;
    int	                         _nbo;
    float**                      _noiseArray;
    int                          *_numberPixelsPerX;
    int                          *_numberPixelsPerY;
    QMainWindow                  *_parent;
    float***                     _paramsArray;
    int                          _patience;
    fftwf_plan		             _plan;
    char**                       _pointFlagsArray;
    bool                         _saveTitleLine;
    float                        *_slope;
    float**                      _sonogramArray;
    int			                 _sonogramWidth;
    SNDFILE*		             _soundFile;
    SF_INFO			             _soundFileInfo;
    int                          _stopThreshold;
    float                        **_tabX;
    float                        *_tabY;
    float                        **_tabYX;
    int	                         _timeExpansion;
    bool                         _treating;
    QFile                        _txtFile;
    QString                      _txtFilePath;
    QString                      _txtPath;
    QFile                        _txtFile2;
    QString                      _txtFilePath2;
    uint                         _tvaleur[2000];
    QVector< QPoint >            _vectorCallPoints;
    QVector < int >              _vectorXMin;
    QVector < int >              _vectorXMax;
    QVector< ParamToSave >       _vectPar;
    QStringList                  _wavFileList;
    QString                      _wavPath;
    int                          _xMax;
    int                          *_xMaxPerY;
    int                          _xMin;
    int                          *_xMinPerY;
    bool                         _xmoitie;
    int                          *_xSecondWestRidgePerY;
    int                          **_yEmaxPerX;
    int                          _yMax;
    int                          *_yMaxPerX;
    int                          _yMin;
    int                          *_yMinPerX;
    int                          **_ypm;
};

#endif
