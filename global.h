#include "define.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <direct.h>
#include <stdlib.h>	//exit�֐�


#ifdef GLOBAL_INSTANCE
#define GLOBAL 
//#define GNUPLOT_PATH "C:/PROGRA~2/gnuplot/bin/gnuplot.exe"

#else
#define GLOBAL extern
#endif


#include "function.h"


/*�ϐ��錾��*/
GLOBAL double FLandingTime;
GLOBAL double FMaxHeight;
GLOBAL double FMaxDistanceX;
GLOBAL double FMaxDistanceY;
GLOBAL double FMaxDistanceXt;
GLOBAL double FMaxDistanceYt;
GLOBAL double FMaxDistanceXf;
GLOBAL double FMaxDistanceYf;
GLOBAL double FMaxDistanceXr;
GLOBAL double FMaxDistanceYr;
GLOBAL double FMaxVelocity;
GLOBAL double FMaxTime;
GLOBAL double FWindDeg;
GLOBAL double FWindAbs;
GLOBAL double FPara1Time;

/*�\���̕ϐ��錾��*/
GLOBAL Wind_data FWind;

/*���͊֌W*/
//GLOBAL double data[266/*450*/];
GLOBAL double FThrustTimeStamp[512];
GLOBAL double FThrust[512];
GLOBAL int FThrustCount;
GLOBAL int FThrustIndex;