#define _CRT_SECURE_NO_WARNINGS
#include "global.h"

//インパルス計算巻数
double impulse(double t) {
    double Thrust;
    if (t < FThrustTimeStamp[0] || t >FThrustTimeStamp[FThrustCount -1]) {
        Thrust = 0;
        return Thrust;
    }
    while (t >= FThrustTimeStamp[0] && t <= FThrustTimeStamp[FThrustCount -1]) {
        if (t < FThrustTimeStamp[FThrustIndex + 1]) {
            Thrust = FThrust[FThrustIndex + 1] + (FThrust[FThrustIndex + 1] - FThrust[FThrustIndex]) * (t - FThrustTimeStamp[FThrustIndex + 1]) / (FThrustTimeStamp[FThrustIndex + 1] - FThrustTimeStamp[FThrustIndex]);
            return Thrust;
        }
        FThrustIndex += 1;
    }
    Thrust = 0;
    return Thrust;
}