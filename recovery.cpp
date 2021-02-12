#define _CRT_SECURE_NO_WARNINGS
#include "global.h"


//------------------------------ƒpƒ‰ƒVƒ…[ƒg“WŠJŽž----------------------------------//
void recovery(Rocket_data *Rocket, Wind_data *Wind, double *t) {
    //FPara1Time = PARA_1t;	//C³ˆÄ
    int flag = 0;
    //Rocket->v_y = (Rocket->v_abs) * sin(Rocket->gamma);		//³‚µ‚¢H

    FILE *fpXYZ, *fpV;
    if ((fpXYZ = fopen("recovery_XYZ_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_errorrecovery_XYZ_data");
        exit(1);
    }
    if ((fpV = fopen("recovery_V_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_errorrecovery_V_data");
        exit(1);
    }

    double initt = 0;
    initt = *t;

    while (Rocket->z > 0) {
        rungekutta_recovery(Rocket, Wind, (*t - initt));
        if ((int)(*t / DT) % (int)(SAMPLING_T / DT) == 0) {
            fprintf(fpXYZ, "%lf\t%lf\t%lf\t%lf\n", *t, Rocket->x, Rocket->y, Rocket->z);
            fprintf(fpV, "%lf\t%lf\n", *t, Rocket->v_y_abs);
        }
        *t += DT;
        //if ((*t > FPara1Time + PARA_2tIME + PARA_2_ERROR && (TYPE == 0 || TYPE == 3)) && flag == 0) {
        //	Rocket->v_abs = Rocket->v_z_abs;
        //	flag = 1;
        //}
        //—Ìˆæ(x•ûŒü)
        if (abs(Rocket->x) >= (FMaxDistanceXr)) {
            FMaxDistanceXr = MAX_DISTANCE(abs(Rocket->x), FMaxDistanceXr);
        }
        //—Ìˆæ(y•ûŒü)
        if (abs(Rocket->y) >= (FMaxDistanceYr)) {
            FMaxDistanceYr = MAX_DISTANCE(abs(Rocket->y), FMaxDistanceYr);
        }
    }
}

void rungekutta_recovery(Rocket_data *Rocket, Wind_data *Wind, double t) {
    double k1[6], k2[6], k3[6], k4[6];
    Rocket_data RocketTemp[4] =
    {
         *Rocket, *Rocket, *Rocket, *Rocket
    };

    Wind_data WindTemp[4] =
    {
         *Wind, *Wind, *Wind, *Wind
    };

    normalization_wind(Wind, *Rocket);

    //======================Step01=============================================//

    k1[V_X] = DT * recovery_v_x(RocketTemp[0], WindTemp[0], t);
    k1[V_Y] = DT * recovery_v_y(RocketTemp[0], WindTemp[0], t);
    k1[V_Z] = DT * recovery_v_z(RocketTemp[0], WindTemp[0], t);

    k1[X] = DT * recovery_x(RocketTemp[0], WindTemp[0], t);
    k1[Y] = DT * recovery_y(RocketTemp[0], WindTemp[0], t);
    k1[Z] = DT * recovery_z(RocketTemp[0], WindTemp[0], t);

    RocketTemp[1].v_x_abs += k1[V_X] / 2.0;
    RocketTemp[1].x += k1[X] / 2.0;

    RocketTemp[1].v_y_abs += k1[V_Y] / 2.0;
    RocketTemp[1].y += k1[Y] / 2.0;

    RocketTemp[1].v_z_abs += k1[V_Z] / 2.0;
    RocketTemp[1].z += k1[Z] / 2.0;

    //=========================================================================//

    //======================Step02=============================================//

    k2[V_X] = DT * recovery_v_x(RocketTemp[1], WindTemp[1], (t + DT / 2.0));
    k2[V_Y] = DT * recovery_v_y(RocketTemp[1], WindTemp[1], (t + DT / 2.0));
    k2[V_Z] = DT * recovery_v_z(RocketTemp[1], WindTemp[1], (t + DT / 2.0));

    k2[X] = DT * recovery_x(RocketTemp[1], WindTemp[1], (t + DT / 2.0));
    k2[Y] = DT * recovery_y(RocketTemp[1], WindTemp[1], (t + DT / 2.0));
    k2[Z] = DT * recovery_z(RocketTemp[1], WindTemp[1], (t + DT / 2.0));

    RocketTemp[2].v_x_abs += k2[V_X] / 2.0;
    RocketTemp[2].x += k2[X] / 2.0;

    RocketTemp[2].v_y_abs += k2[V_Y] / 2.0;
    RocketTemp[2].y += k2[Y] / 2.0;

    RocketTemp[2].v_z_abs += k2[V_Z] / 2.0;
    RocketTemp[2].z += k2[Z] / 2.0;

    //=========================================================================//

    //======================Step03=============================================//

    k3[V_X] = DT * recovery_v_x(RocketTemp[2], WindTemp[2], (t + DT / 2.0));
    k3[V_Y] = DT * recovery_v_y(RocketTemp[2], WindTemp[2], (t + DT / 2.0));
    k3[V_Z] = DT * recovery_v_z(RocketTemp[2], WindTemp[2], (t + DT / 2.0));

    k3[X] = DT * recovery_x(RocketTemp[2], WindTemp[2], (t + DT / 2.0));
    k3[Y] = DT * recovery_y(RocketTemp[2], WindTemp[2], (t + DT / 2.0));
    k3[Z] = DT * recovery_z(RocketTemp[2], WindTemp[2], (t + DT / 2.0));

    RocketTemp[3].v_x_abs += k3[V_X] / 2.0;
    RocketTemp[3].x += k3[X];

    RocketTemp[3].v_y_abs += k3[V_Y];
    RocketTemp[3].y += k3[Y];

    RocketTemp[3].v_z_abs += k3[V_Z];
    RocketTemp[3].z += k3[Z];

    //=========================================================================//

    //======================Step04=============================================//

    k4[V_X] = DT * recovery_v_x(RocketTemp[3], WindTemp[3], t + DT);
    k4[V_Y] = DT * recovery_v_y(RocketTemp[3], WindTemp[3], t + DT);
    k4[V_Z] = DT * recovery_v_z(RocketTemp[3], WindTemp[3], t + DT);

    k4[X] = DT * recovery_x(RocketTemp[3], WindTemp[3], t + DT);
    k4[Y] = DT * recovery_y(RocketTemp[3], WindTemp[3], t + DT);
    k4[Z] = DT * recovery_z(RocketTemp[3], WindTemp[3], t + DT);

    //=========================================================================//

    //======================Step05=============================================//

    Rocket->v_x_abs += (1.0 / 6.0) * (k1[V_X] + 2.0 * k2[V_X] + 2.0 * k3[V_X] + k4[V_X]);
    Rocket->x += (1.0 / 6.0) * (k1[X] + 2.0 * k2[X] + 2.0 * k3[X] + k4[X]);

    Rocket->v_y_abs += (1.0 / 6.0) * (k1[V_Y] + 2.0 * k2[V_Y] + 2.0 * k3[V_Y] + k4[V_Y]);
    Rocket->y += (1.0 / 6.0) * (k1[Y] + 2.0 * k2[Y] + 2.0 * k3[Y] + k4[Y]);

    Rocket->v_z_abs += (1.0 / 6.0) * (k1[V_Z] + 2.0 * k2[V_Z] + 2.0 * k3[V_Z] + k4[V_Z]);
    Rocket->z += (1.0 / 6.0) * (k1[Z] + 2.0 * k2[Z] + 2.0 * k3[Z] + k4[Z]);

    //=========================================================================//

    //======absolute_of_velocity============================================//
    Rocket->v_abs = sqrt(Rocket->v_x_abs * Rocket->v_x_abs + Rocket->v_y_abs * Rocket->v_y_abs + Rocket->v_z_abs * Rocket->v_z_abs);
    //======================================================================//

}

double recovery_v_x(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = (-1.0) * ((Cd_2 * S_2 * rho(Rocket)) / (2 * Rocket.m)) * (-Wind.w_x + Rocket.v_x_abs) * abs(-Wind.w_x + Rocket.v_x_abs);
    return ans;
}

double recovery_x(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_x_abs;
    return ans;
}

double recovery_v_z(Rocket_data Rocket, Wind_data Wind, double t) {
    FPara1Time = PARA_1t;	//C³ˆÄ
    double ans, para_s = PARA_S_1;
    if (((t > FPara1Time - PARA_1_ERROR + PARA_2tIME + PARA_2_ERROR) && TYPE == 0) || TYPE == 3) {
        para_s = PARA_S_2;
    }
    ans = (((Cd + PARA_CD) * (S + para_s) * rho(Rocket)) / (2 * Rocket.m)) * Rocket.v_z_abs *  Rocket.v_z_abs - G;
    return ans;
}

double recovery_z(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_z_abs;
    return ans;
}

double recovery_v_y(Rocket_data Rocket, Wind_data Wind, double t) {
    //FPara1Time = PARA_1t;	//C³ˆÄ
    //double ans, para_s = PARA_S_1;
    //if (t > FPara1Time + PARA_2tIME + PARA_2_ERROR && (TYPE == 0 || TYPE == 3)) {
    //	para_s = PARA_S_2;
    //}
    double ans;
    ans = (-1.0) * ((Cd_2 * S_2 * rho(Rocket)) / (2 * Rocket.m)) * (-Wind.w_y + Rocket.v_y_abs) * abs(-Wind.w_y + Rocket.v_y_abs);
    return ans;
}

double recovery_y(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_y_abs;
    return ans;
}
//-----------------------------------------------------------------------------//
