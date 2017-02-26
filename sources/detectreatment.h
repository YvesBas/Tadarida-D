#ifndef DETECTREATMENT_H
#define DETECTREATMENT_H

#include <float.h>
#include <iostream>
#include <QDate>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QtCore/qmath.h>
#include <QTextStream>
#include <QProcess>
#include <QVector>
#include "sndfile.h"
#include "fftw3.h"

#include <QPoint>

#define	FFT_HEIGHT_MAX 4096
#define	FFT_HEIGHT_HALF_MAX 2048
#define FREQ_MAX 250 // KHz
#define MAXLARCRI 10000
#define MAXHEIGHT 500
#define PI	3.14159265358979
#define SONOGRAM_WIDTH_MAX 60008
#define  LD8 (SONOGRAM_WIDTH_MAX+15)/8
#define SONOGRAM_WIDTH_MIN 64
#define MAXCRI 1500
#define NCRETES 5
#define EMIN -100

#define WITHNEWPARAM true

class Detec;
class Recherche;

enum NUMTAB {SH,CM,CS,CN,CO,CO2,NBTAB};
enum NUMPAR {StTime,Dur,PrevSt,Fmax,Fmin,BW,FPk,FPkD,TPk,Slope,ISlope,HCF,FIF,THCF,UpSl,LoSl,StF,
     EnF,StSl,EnSl,FPSl,FISl,CeF,B5dBBF,B5dBAF,B5dBBW,B5dBDur,LCF,FreqMP,
     PosMP,FreqPkS,FreqPkM,PosPkS,PosPkM,FreqPkS2,FreqPkM2,PrevMP1,
     PrevMP2,NextMP1,NextMP2,Amp1,Amp2,Amp3,Amp4,NoisePrev,NoiseNext,NoiseDown,NoiseUp,
     CVAmp,Hup_RFMP,Hup_PosMP,Hup_PosSt,Hup_PosEn,Hup_AmpDif,Hup_RSlope,Hlo_RFMP,Hlo_PosMP,
     Hlo_PosSt,Hlo_PosEn,Hlo_AmpDif,Hlo_RSlope,Ramp_2_1,Ramp_3_1,Ramp_3_2,Ramp_1_2,Ramp_4_3,
     Ramp_2_3,RAN_2_1,RAN_3_1,RAN_3_2,RAN_1_2,RAN_4_3,RAN_2_3,HetX,HetY,Dbl8,Stab,Rythm1,Rythm2,Rythm3,
     HeiET,HeiEM,HeiRT,HeiRM,HeiETT,HeiEMT,HeiRTT,HeiRMT,MedInt,Int25,Int75,RInt1,IntDev,SmIntDev,LgIntDev,
     VarInt,VarSmInt,VarLgInt,RIntDev1,EnStabSm,EnStabLg,HetXr,HetYr,HetCMC,HetCMD,HetCTC,HetCTD,
     HetCMnP,HetCMfP,HetCTnP,HetCTfP,
     HetPicsMAD,HetPicsMALD,HetPicsMABD,HetPicsMRBLD,HetPicsTAD,HetPicsTALD,HetPicsTABD,HetPicsTRBLD,
     VDPicsM,VLDPicsM,VBDPicsM,VDPPicsM,VLDPPicsM,VBDPPicsM,
     VDPicsT,VLDPicsT,VBDPicsT,VDPPicsT,VLDPPicsT,VBDPPicsT,
     HetYr2,SDC,SDCR,SDCRY,SDCRXY,
     SDCL,SDCLR,SDCLRY,SDCLRXY,SDCLRXY2,
     SDCLOP,SDCLROP,SDCLRYOP,SDCLRXYOP,
     SDCLWB,SDCLRWB,SDCLRYWB,SDCLRXYWB,
     SDCLOPWB,SDCLROPWB,SDCLRYOPWB,SDCLRXYOPWB,
     SDCL_DNP,SDCLR_DNP,SDCLRY_DNP,SDCLRXY_DNP,SDCLRXY2_DNP,
     ELBPOS,ELBSB,ELB2POS,ELB2SB,
     RAF,RAE,RAFE,RAFP,RAFP2,RAFP3,
     SBMP,SAMP,SBAR,
     RAHP2,RAHP4,RAHP8,RAHP16,RAHE2,RAHE4,RAHE8,RAHE16,
     NBPAR};

enum NUMERROR {FNREC,MCNT,DTP,DTG,TNT,NTERRORS};


class ParamToSave
{
public:
    ParamToSave(int,int,QString);
    ParamToSave(int,int,QString,int);
    ParamToSave(int,int,QString,int,int);
    ParamToSave();
    ~ParamToSave();

    int                                                                    ArrayNumber;
    QString                                                           ColumnTitle;
    int                                                                    FromVersion;
    QString                                                           InfoLabel;
    int                                                                    ParameterNumber;
    int                                                                    ToVersion;
};


class DetecTreatment
{
public:
    DetecTreatment(Detec *);
    DetecTreatment();
    ~DetecTreatment();
    bool                                                                 CallTreatmentsForOneFile(QString&,QString &);
    void                                                                 EndDetecTreatment();
    void                                                                 InitializeDetecTreatment();
    void                                                                 SetDirParameters(QString,QString,bool,QString,QString);
    void                                                                 SetGlobalParameters(int,int,int,int,int,int,int,bool,int,int,int,int,int,int,int,int,int,int,bool);
    void                                                                 SortFloatIndArray(float *,int,int *);

    QVector< QVector<QPoint> >                    CallsArray;
    QVector< QVector<QPoint> >                    CallMasterRidgeArray;
    QVector< QVector<QPoint> >                    CallNorthRidgeArray;
    QVector< QVector<QPoint> >                    CallSecondWestRidgeArray;
    QVector< QVector<QPoint> >                    CallSouthArray;
    QVector< QVector<QPoint> >                    CallWestRidgeArray;
    float                                                                *EnergyColumAverage;
    double		                                                        EnergyMin;
    double			                                                    EnergyMax;
    double                                                              EnergyShapeThreshold;
    double                                                              EnergyStopThreshold;
    bool                                                                 *FlagGoodCol;
    bool                                                                 *FlagGoodColInitial;
    int			                                                        FftHeightHalf;
    int                                                                    *Inflexion1;
    int                                                                    *Inflexion3;
    double                                                                 KhzPerY;
    int                                                                    LimY;
    QVector< QPoint >                                       MasterPoints;
    double                                                                 MsPerX;
    int                                                                    NError;
    char**                                                              PointFlagsArray;
    qint16 **                                                         SonogramArray;
    int			                                                        SonogramWidth;
    int                                                                    TabErrors[NTERRORS];
    int	                                                                TimeExpansion;
    QVector< ParamToSave >                           VectPar;
    bool                                                                 WithSilence;
    int                                                                     *LowSlope;

private:
    void                                                                 aff(QString,qint64,int);
    void                                                                 clearVars();
    bool                                                                 computeFFT(QString &);
    void                                                                 correctNoise();
    void                                                                 detectsParameter2();
    bool                                                                 determineLeftOrRight(QString &);
    void                                                                 initVectorParams();
    bool                                                                 openWavFile(QString &);
    void                                                                 saveParameters(const QString&);
    void                                                                 saveCompressedParameters(const QString&);
    void                                                                 shapesDetects();
    void                                                                 sortWaves();
    void                                                                 sortDoubleArray(double *,int);
    void                                                                 sortIntArrays(int *,int,int *);

    float *                                                               _averagePerX;
    QDir                                                                _baseDayDir;
    double                                                               _callEnergyMax ;
    int                                                                    _callEnergyMaxIndex;
    QVector<QPoint>                                          _callMasterRidge;
    QVector<QPoint>                                          _callNorthRidge;
    QVector<QPoint>                                         _callSecondWestRidge;
    QVector<QPoint>                                         _callSouthRidge;
    QVector<QPoint>                                         _callWestRidge;
    int                                                                    _callsNumber;
    char *                                                              _charParamsArray;
    char*                                                               _charPointFlagsArray;
    char*                                                               _charSonogramArray;
    char *                                                              _charTabX;
    char *                                                              _charTabYX;
    char *                                                              _charYEmaxPerX;
    double  *                                                             _coeff;
    fftwf_complex*                                               _complexInput;
    bool                                                                _desactiveCorrectNoise;
    int                                                                    _compressedVersion;
    float *                                                              _data;
    QString                                                           _datPath;
    Detec                                                              *_detec;
    int                                                                    _detectionThreshold;
    float **                                                             _dpm;
    int **                                                                _dypm;
    float *                                                              _eMaxPerX;
    int                                                                   _fftHeight;
    fftwf_complex*	                                           _fftRes;
    float                                                                _freqCallMax;
    float                                                                _freqCallMin;
    int                                                                   _freqMin;
    int                                                                   **_harmonic;
    int                                                                   _highThresholdJB;
    int                                                                   _highThresholdC;
    int                                                                   _iH;
    bool                                                               _imageData;
    QString                                                          _imagePath;
    int*                                                                  _invMp;
    int                                                                   _iOverlapMoving;
    int                                                                   _jumpThreshold;
    bool                                                               _littleParams;
    int                                                                   _lowThresholdJB;
    int                                                                   _lowThresholdC;
    int                                                                   _maxCallWidth;
    int                                                                   _maxCallHeight;
    int                                                                   _maxY;
    int                                                                   _modeFreq;
    int                                                                   _minY;
    int                                                                   _nbo;
    int                                                                   _numberCallParameters;
    int                                                                  *_numberPixelsPerX;
    int                                                                  *_numberPixelsPerY;
    float***                                                          _paramsArray;
    int                                                                   _paramVersion;
    fftwf_plan                                                      *_pPlan;
    int                                                                  _qN;
    int                                                                  _qR;
    QString                                                         _resultCompressedSuffix;
    QString                                                         _resultSuffix;
    bool                                                               _saveTitleLine;
    float*                                                              _slope;
    float**                                                            _simpleParamsArray;
    int*                                                                _sortMp;
    SNDFILE*                                                    _soundFile;
    SF_INFO                                                      _soundFileInfo;
    int                                                                  _stopThreshold;
    float**                                                            _tabX;
    float*                                                              _tabY;
    float**                                                            _tabYX;
    int                                                                  _timeExpansionLeft;
    int                                                                  _timeExpansionRight;
    QString                                                        _txtPath;
    uint                                                                _tvalue[2000];
    bool                                                               _useValflag;
    float**                                                            _valuesToCompressArray;
    QVector< QPoint >                                     _vectorCallPoints;
    QVector < int >                                            _vectorXMin;
    QVector < int >                                            _vectorXMax;
    QString                                                        _wavFile;
    QStringList                                                  _wavFileList;
    QString                                                        _wavPath;
    int                                                                 _widthBigControl;
    int                                                                 _widthLittleControl;
    int                                                                  _xMax;
    int*                                                                _xMaxPerY;
    int                                                                 _xMin;
    quint16*                                                       _xMinPerY;
    int *                                                               _xMp;
    quint16*                                                       _xSecondWestRidgePerY;
    float*                                                             _totYEPerX;
    float*                                                            _yEbarPerX;
    quint16**                                                     _yEmaxPerX;
    int                                                                 _yMax;
    quint16*                                                       _yMaxPerX;
    int                                                                 _yMin;
    quint16*                                                       _yMinPerX;
    int**                                                              _ypm;
};

#endif // DETECTREATMENT_H
