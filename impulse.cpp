#define _CRT_SECURE_NO_WARNINGS
#include "global.h"

//インパルス計算巻数
double impulse(double t) {
    double T;
    if (t < ThrustTimeStamp[0] || t >ThrustTimeStamp[ThrustCount -1]) {
        T = 0;
        return T;
    }
    while (t >= ThrustTimeStamp[0] && t <= ThrustTimeStamp[ThrustCount -1]) {
        if (/*t >= ThrustTimeStamp[ThrustIndex] &&*/ t < ThrustTimeStamp[ThrustIndex + 1]/* && t <= ThrustTimeStamp[ThrustCount]*/) {
            T = Thrust[ThrustIndex + 1] + (Thrust[ThrustIndex + 1] - Thrust[ThrustIndex]) * (t - ThrustTimeStamp[ThrustIndex + 1]) / (ThrustTimeStamp[ThrustIndex + 1] - ThrustTimeStamp[ThrustIndex]);
            return T;
        }
        ThrustIndex += 1;
    }
    T = 0;
    return T;
}