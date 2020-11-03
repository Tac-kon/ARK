#define _CRT_SECURE_NO_WARNINGS
#include "global.h"


//------------------------------ƒpƒ‰ƒVƒ…[ƒg“WŠJŽž----------------------------------//
void recovery(Rocket_data *Rocket, Wind_data *Wind, double *t) {
    //PARA_1_TIME = PARA_1_T;	//C³ˆÄ
    int flag = 0;
    //Rocket->v_y = (Rocket->v_abs) * sin(Rocket->gamma);		//³‚µ‚¢H

    FILE *fp1, *fp2;
    if ((fp1 = fopen("recovery_XYZ_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_error_recovery_XYZ_data");
        exit(1);
    }
    if ((fp2 = fopen("recovery_V_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_error_recovery_V_data");
        exit(1);
    }

    double init_t = 0;
    init_t = *t;

    while (Rocket->z > 0) {
        rungekutta_recovery(Rocket, Wind, (*t - init_t));
        if ((int)(*t / DT) % (int)(SAMPLING_T / DT) == 0) {
            fprintf(fp1, "%lf\t%lf\t%lf\t%lf\n", *t, Rocket->x, Rocket->y, Rocket->z);
            fprintf(fp2, "%lf\t%lf\n", *t, Rocket->v_y_abs);
        }
        *t += DT;
        //if ((*t > PARA_1_TIME + PARA_2_TIME + PARA_2_ERROR && (TYPE == 0 || TYPE == 3)) && flag == 0) {
        //	Rocket->v_abs = Rocket->v_z_abs;
        //	flag = 1;
        //}
        //—Ìˆæ(x•ûŒü)
        if (abs(Rocket->x) >= (max_distance_x_r)) {
            max_distance_x_r = MAX_DISTANCE(abs(Rocket->x), max_distance_x_r);
        }
        //—Ìˆæ(y•ûŒü)
        if (abs(Rocket->y) >= (max_distance_y_r)) {
            max_distance_y_r = MAX_DISTANCE(abs(Rocket->y), max_distance_y_r);
        }
    }
}

void rungekutta_recovery(Rocket_data *Rocket, Wind_data *Wind, double t) {
    double k1[6], k2[6], k3[6], k4[6];
    Rocket_data Rocket_temp[4] =
    {
         *Rocket, *Rocket, *Rocket, *Rocket
    };

    Wind_data Wind_temp[4] =
    {
         *Wind, *Wind, *Wind, *Wind
    };

    normalization_wind(Wind, *Rocket);

    //======================Step01=============================================//

    k1[V_X] = DT * recovery_v_x(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Y] = DT * recovery_v_y(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Z] = DT * recovery_v_z(Rocket_temp[0], Wind_temp[0], t);

    k1[X] = DT * recovery_x(Rocket_temp[0], Wind_temp[0], t);
    k1[Y] = DT * recovery_y(Rocket_temp[0], Wind_temp[0], t);
    k1[Z] = DT * recovery_z(Rocket_temp[0], Wind_temp[0], t);

    Rocket_temp[1].v_x_abs += k1[V_X] / 2.0;
    Rocket_temp[1].x += k1[X] / 2.0;

    Rocket_temp[1].v_y_abs += k1[V_Y] / 2.0;
    Rocket_temp[1].y += k1[Y] / 2.0;

    Rocket_temp[1].v_z_abs += k1[V_Z] / 2.0;
    Rocket_temp[1].z += k1[Z] / 2.0;

    //=========================================================================//

    //======================Step02=============================================//

    k2[V_X] = DT * recovery_v_x(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[V_Y] = DT * recovery_v_y(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[V_Z] = DT * recovery_v_z(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    k2[X] = DT * recovery_x(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[Y] = DT * recovery_y(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[Z] = DT * recovery_z(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    Rocket_temp[2].v_x_abs += k2[V_X] / 2.0;
    Rocket_temp[2].x += k2[X] / 2.0;

    Rocket_temp[2].v_y_abs += k2[V_Y] / 2.0;
    Rocket_temp[2].y += k2[Y] / 2.0;

    Rocket_temp[2].v_z_abs += k2[V_Z] / 2.0;
    Rocket_temp[2].z += k2[Z] / 2.0;

    //=========================================================================//

    //======================Step03=============================================//

    k3[V_X] = DT * recovery_v_x(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[V_Y] = DT * recovery_v_y(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[V_Z] = DT * recovery_v_z(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    k3[X] = DT * recovery_x(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[Y] = DT * recovery_y(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[Z] = DT * recovery_z(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    Rocket_temp[3].v_x_abs += k3[V_X] / 2.0;
    Rocket_temp[3].x += k3[X];

    Rocket_temp[3].v_y_abs += k3[V_Y];
    Rocket_temp[3].y += k3[Y];

    Rocket_temp[3].v_z_abs += k3[V_Z];
    Rocket_temp[3].z += k3[Z];

    //=========================================================================//

    //======================Step04=============================================//

    k4[V_X] = DT * recovery_v_x(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Y] = DT * recovery_v_y(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Z] = DT * recovery_v_z(Rocket_temp[3], Wind_temp[3], t + DT);

    k4[X] = DT * recovery_x(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Y] = DT * recovery_y(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Z] = DT * recovery_z(Rocket_temp[3], Wind_temp[3], t + DT);

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
    PARA_1_TIME = PARA_1_T;	//C³ˆÄ
    double ans, para_s = PARA_S_1;
    if (((t > PARA_1_TIME - PARA_1_ERROR + PARA_2_TIME + PARA_2_ERROR) && TYPE == 0) || TYPE == 3) {
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
    //PARA_1_TIME = PARA_1_T;	//C³ˆÄ
    //double ans, para_s = PARA_S_1;
    //if (t > PARA_1_TIME + PARA_2_TIME + PARA_2_ERROR && (TYPE == 0 || TYPE == 3)) {
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
