#define _CRT_SECURE_NO_WARNINGS
#include "global.h"


//--------------------------------Šµ«”òsŽž(‰º~Žž)-----------------------------//
void inertia_fall(Rocket_data *Rocket, Wind_data *Wind, double *t) {
    FILE *fp1, *fp2;
    double init_t = 0;
    if ((fp1 = fopen("inertia_fall_XYZ_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_error_inertia_fall_XYZ_data");
        exit(1);
    }
    if ((fp2 = fopen("inertia_fall_V_data.txt", "w")) == NULL) {
        fprintf(stderr, "File_error_inertia_fall_V_data");
        exit(1);
    }
    normalization_wind(Wind, *Rocket);
    init_t = *t;
    //Rocket->v_x = sin(Rocket->gamma) * cos(Rocket->psi) * Rocket->v_abs;
    //Rocket->v_z = sin(Rocket->gamma) * sin(Rocket->psi) * Rocket->v_abs;
    //Rocket->v_y = cos(Rocket->gamma) * Rocket->v_abs;
    while ((TYPE !=2 && (*t < max_time + PARA_1_ERROR)) || (TYPE == 2 && Rocket->z > 0))/*&& *t <= PARA_1_TIME + PARA_1_ERROR && (TYPE == 0 || TYPE == 1)) || (Rocket->z > 0 && TYPE == 2))*/ {
        rungekutta_inertia_fall(Rocket, Wind, (*t - init_t));
        //printf("wind_abs =%lf\n", Wind->w_abs);
        //printf("wind_x =%lf\n",Wind->w_x);
        //printf("a =%lf\n", Rocket->v_x_abs);
        if ((int)(*t / DT) % (int)(SAMPLING_T / DT) == 0) {
            fprintf(fp1, "%lf\t%lf\t%lf\t%lf\n", *t, Rocket->x, Rocket->y, Rocket->z);
            fprintf(fp2, "%lf\t%lf\t%lf\t%lf\n", *t, Rocket->v_abs, Wind->w_x, Wind->w_y);
        }
        *t += DT;
        //—Ìˆæ(x•ûŒü)
        if (abs(Rocket->x) >= (max_distance_x_f)) {
            max_distance_x_f = MAX_DISTANCE(abs(Rocket->x), max_distance_x_f);
        }
        //—Ìˆæ(y•ûŒü)
        if (abs(Rocket->y) >= (max_distance_y_f)) {
            max_distance_y_f = MAX_DISTANCE(abs(Rocket->y), max_distance_y_f);
        }
        landing_time = TIME_DEVELOP(*t, landing_time);
    }
    //Rocket->v_abs = sqrt(Rocket->v_x_abs * Rocket->v_x_abs + Rocket->v_y_abs * Rocket->v_y_abs + Rocket->v_z_abs * Rocket->v_z_abs);
    //Rocket->v_abs = Rocket->v_z_abs;
    fclose(fp1);
    fclose(fp2);
}

void rungekutta_inertia_fall(Rocket_data *Rocket, Wind_data *Wind, double t) {
    double k1[12], k2[12], k3[12], k4[12];
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

    //•Ài‰^“®
    k1[V_X] = DT * inertia_fall_equation_v_x(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Y] = DT * inertia_fall_equation_v_y(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Z] = DT * inertia_fall_equation_v_z(Rocket_temp[0], Wind_temp[0], t);

    k1[X] = DT * inertia_fall_equation_x(Rocket_temp[0], Wind_temp[0], t);
    k1[Y] = DT * inertia_fall_equation_y(Rocket_temp[0], Wind_temp[0], t);
    k1[Z] = DT * inertia_fall_equation_z(Rocket_temp[0], Wind_temp[0], t);

    //‰ñ“]‰^“®
    k1[R_V_X] = DT * fall_rotate_equation_x_2(Rocket_temp[0], Wind_temp[0], t);
    k1[R_V_Y] = DT * fall_rotate_equation_y_2(Rocket_temp[0], Wind_temp[0], t);
    k1[R_V_Z] = DT * fall_rotate_equation_z_2(Rocket_temp[0], Wind_temp[0], t);

    k1[R_X] = DT * fall_rotate_equation_x_1(Rocket_temp[0], Wind_temp[0], t);
    k1[R_Y] = DT * fall_rotate_equation_y_1(Rocket_temp[0], Wind_temp[0], t);
    k1[R_Z] = DT * fall_rotate_equation_z_1(Rocket_temp[0], Wind_temp[0], t);

    Rocket_temp[1].v_x_abs += k1[V_X] / 2.0;
    Rocket_temp[1].x += k1[X] / 2.0;

    Rocket_temp[1].v_y_abs += k1[V_Y] / 2.0;
    Rocket_temp[1].y += k1[Y] / 2.0;

    Rocket_temp[1].v_z_abs += k1[V_Z] / 2.0;
    Rocket_temp[1].z += k1[Z] / 2.0;

    Rocket_temp[1].r_v_x += k1[R_V_X] / 2.0;
    Rocket_temp[1].r_x += k1[R_X] / 2.0;

    Rocket_temp[1].r_v_y += k1[R_V_Y] / 2.0;
    Rocket_temp[1].r_y += k1[R_Y] / 2.0;

    Rocket_temp[1].r_v_z += k1[R_V_Z] / 2.0;
    Rocket_temp[1].r_z += k1[R_Z] / 2.0;

    //=========================================================================//

    //======================Step02=============================================//

    //•Ài‰^“®
    k2[V_X] = DT * inertia_fall_equation_v_x(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[V_Y] = DT * inertia_fall_equation_v_y(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[V_Z] = DT * inertia_fall_equation_v_z(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    k2[X] = DT * inertia_fall_equation_x(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[Y] = DT * inertia_fall_equation_y(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[Z] = DT * inertia_fall_equation_z(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    //‰ñ“]‰^“®
    k2[R_V_X] = DT * fall_rotate_equation_x_2(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[R_V_Y] = DT * fall_rotate_equation_y_2(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[R_V_Z] = DT * fall_rotate_equation_z_2(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    k2[R_X] = DT * fall_rotate_equation_x_1(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[R_Y] = DT * fall_rotate_equation_y_1(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));
    k2[R_Z] = DT * fall_rotate_equation_z_1(Rocket_temp[1], Wind_temp[1], (t + DT / 2.0));

    Rocket_temp[2].v_x_abs += k2[V_X] / 2.0;
    Rocket_temp[2].x += k2[X] / 2.0;

    Rocket_temp[2].v_y_abs += k2[V_Y] / 2.0;
    Rocket_temp[2].y += k2[Y] / 2.0;

    Rocket_temp[2].v_z_abs += k2[V_Z] / 2.0;
    Rocket_temp[2].z += k2[Z] / 2.0;
    
    Rocket_temp[2].r_v_x += k2[R_V_X] / 2.0;
    Rocket_temp[2].r_x += k2[R_X] / 2.0;

    Rocket_temp[2].r_v_y += k2[R_V_Y] / 2.0;
    Rocket_temp[2].r_y += k2[R_Y] / 2.0;

    Rocket_temp[2].r_v_z += k2[R_V_Z] / 2.0;
    Rocket_temp[2].r_z += k2[R_Z] / 2.0;

    //=========================================================================//

    //======================Step03=============================================//

    //•Ài‰^“®
    k3[V_X] = DT * inertia_fall_equation_v_x(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[V_Y] = DT * inertia_fall_equation_v_y(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[V_Z] = DT * inertia_fall_equation_v_z(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    k3[X] = DT * inertia_fall_equation_x(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[Y] = DT * inertia_fall_equation_y(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[Z] = DT * inertia_fall_equation_z(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    //‰ñ“]‰^“®
    k3[R_V_X] = DT * fall_rotate_equation_x_2(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[R_V_Y] = DT * fall_rotate_equation_y_2(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[R_V_Z] = DT * fall_rotate_equation_z_2(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    k3[R_X] = DT * fall_rotate_equation_x_1(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[R_Y] = DT * fall_rotate_equation_y_1(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));
    k3[R_Z] = DT * fall_rotate_equation_z_1(Rocket_temp[2], Wind_temp[2], (t + DT / 2.0));

    Rocket_temp[3].v_x_abs += k3[V_X] / 2.0;
    Rocket_temp[3].x += k3[X];

    Rocket_temp[3].v_y_abs += k3[V_Y];
    Rocket_temp[3].y += k3[Y];

    Rocket_temp[2].v_z_abs += k3[V_Z];
    Rocket_temp[3].z += k3[Z];

    Rocket_temp[3].r_v_x += k3[R_V_X] / 2.0;
    Rocket_temp[3].r_x += k3[R_X] / 2.0;

    Rocket_temp[3].r_v_y += k3[R_V_Y] / 2.0;
    Rocket_temp[3].r_y += k3[R_Y] / 2.0;

    Rocket_temp[2].r_v_z += k3[R_V_Z] / 2.0;
    Rocket_temp[3].r_z += k3[R_Z] / 2.0;

    //=========================================================================//

    //======================Step04=============================================//

    k4[V_X] = DT * inertia_fall_equation_v_x(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Y] = DT * inertia_fall_equation_v_y(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Z] = DT * inertia_fall_equation_v_z(Rocket_temp[3], Wind_temp[3], t + DT);

    k4[X] = DT * inertia_fall_equation_x(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Y] = DT * inertia_fall_equation_y(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Z] = DT * inertia_fall_equation_z(Rocket_temp[3], Wind_temp[3], t + DT);

    //‰ñ“]‰^“®
    k4[R_V_X] = DT * fall_rotate_equation_x_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_V_Y] = DT * fall_rotate_equation_y_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_V_Z] = DT * fall_rotate_equation_z_2(Rocket_temp[3], Wind_temp[3], t + DT);

    k4[R_X] = DT * fall_rotate_equation_x_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_Y] = DT * fall_rotate_equation_y_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_Z] = DT * fall_rotate_equation_z_1(Rocket_temp[3], Wind_temp[3], t + DT);

    //=========================================================================//

    //======================Step05=============================================//

    //•Ài‰^“®
    Rocket->v_x_abs += (1.0 / 6.0) * (k1[V_X] + 2.0 * k2[V_X] + 2.0 * k3[V_X] + k4[V_X]);
    Rocket->x += (1.0 / 6.0) * (k1[X] + 2.0 * k2[X] + 2.0 * k3[X] + k4[X]);

    Rocket->v_y_abs += (1.0 / 6.0) * (k1[V_Y] + 2.0 * k2[V_Y] + 2.0 * k3[V_Y] + k4[V_Y]);
    Rocket->y += (1.0 / 6.0) * (k1[Y] + 2.0 * k2[Y] + 2.0 * k3[Y] + k4[Y]);

    Rocket->v_z_abs += (1.0 / 6.0) * (k1[V_Z] + 2.0 * k2[V_Z] + 2.0 * k3[V_Z] + k4[V_Z]);
    Rocket->z += (1.0 / 6.0) * (k1[Z] + 2.0 * k2[Z] + 2.0 * k3[Z] + k4[Z]);

    //‰ñ“]‰^“®
    Rocket->r_v_x += (1.0 / 6.0) * (k1[R_V_X] + 2.0 * k2[R_V_X] + 2.0 * k3[R_V_X] + k4[R_V_X]);
    Rocket->r_x += (1.0 / 6.0) * (k1[R_X] + 2.0 * k2[R_X] + 2.0 * k3[R_X] + k4[R_X]);

    Rocket->r_v_y += (1.0 / 6.0) * (k1[R_V_Y] + 2.0 * k2[R_V_Y] + 2.0 * k3[R_V_Y] + k4[R_V_Y]);
    Rocket->r_y += (1.0 / 6.0) * (k1[R_Y] + 2.0 * k2[R_Y] + 2.0 * k3[R_Y] + k4[R_Y]);

    Rocket->r_v_z += (1.0 / 6.0) * (k1[R_V_Z] + 2.0 * k2[R_V_Z] + 2.0 * k3[R_V_Z] + k4[R_V_Z]);
    Rocket->r_z += (1.0 / 6.0) * (k1[R_Z] + 2.0 * k2[R_Z] + 2.0 * k3[R_Z] + k4[R_Z]);

    //=========================================================================//

    //======absolute_of_velocity============================================//
    Rocket->v_abs = sqrt(Rocket->v_x_abs * Rocket->v_x_abs + Rocket->v_y_abs * Rocket->v_y_abs + Rocket->v_z_abs * Rocket->v_z_abs);
    //======================================================================//

    //======absolute_of_velocity_0==========================================//
    Rocket->v_abs_0 = sqrt((Rocket->v_x_abs - Wind->w_x) * (Rocket->v_x_abs - Wind->w_x) + (Rocket->v_y_abs - Wind->w_y) * (Rocket->v_y_abs - Wind->w_y) + Rocket->v_z_abs * Rocket->v_z_abs);
    //======================================================================//

    //======gamma===========================================================//
    Rocket->gamma = atan(sqrt((Rocket->x * Rocket->x + Rocket->y * Rocket->y)) / Rocket->z);
    //======================================================================//

    //======psi=============================================================//
    Rocket->psi = asin(Rocket->y / sqrt(Rocket->x * Rocket->x + Rocket->y * Rocket->y));

    //======delta===========================================================//
    Rocket->delta = acos(Rocket->v_z_abs / Rocket->v_abs_0);

    //======alpha===========================================================//
    Rocket->alpha = Rocket->delta - (LAUNCHER_ANGLE + Rocket->r_x);

    //inertia_fall
}

double inertia_fall_equation_v_x(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 1.0 / Rocket.m_ * (-D_0(Rocket, Wind) /** cos(Rocket.alpha) */* sin(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x/* - Rocket7.r_y*/) * cos(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/));
    return ans;
}

double inertia_fall_equation_x(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_x_abs;
    return ans;
}

double inertia_fall_equation_v_y(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 1.0 / Rocket.m_ * (-D_0(Rocket, Wind) /* * cos(Rocket.alpha)*/ * sin(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x/* - Rocket7.r_y*/) * sin(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/));
    //if(ans<-10) printf("%lf %lf\n", 180.0 * (WIND_DEG) / M_PI, ans);
    return ans;
}

double inertia_fall_equation_y(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_y_abs;
    //if (ans<-20) printf("%lf %lf\n", 180.0 * (WIND_DEG) / M_PI, ans);
    return ans;
}

double inertia_fall_equation_v_z(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 1.0 / Rocket.m_ * (-D_0(Rocket, Wind) /** cos(Rocket.alpha)*/ * cos(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x)) - G;
    return ans;
}

double inertia_fall_equation_z(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_z_abs;
    return ans;
}

double fall_rotate_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_x;
    return ans;
}

double fall_rotate_equation_x_2(Rocket_data Rocket,Wind_data wind, double t) {
    double ans;
    ans = D_0(Rocket, Wind) * sin(Rocket.alpha) * Rocket.I_xx_inv * distane_between_Cp_and_Cd;// ((1.0 - (Rocket.alpha) * (Rocket.alpha) / 6.0) * Rocket.I_xx_inv /*+ (1.0 - Rocket7.r_y * Rocket7.r_y / 3.0) *Rocket.I_xy_inv*/);//Taylor“WŠJ
    return ans;
}

double fall_rotate_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_y;
    return ans;
}

double fall_rotate_equation_y_2(Rocket_data Rocket, Wind_data wind, double t) {
    double ans;
    ans = D_0(Rocket, Wind) * sin(Rocket.r_y) * Rocket.I_yy_inv * distane_between_Cp_and_Cd;// (-0.5) * rho(Rocket) * Rocket.v_abs * Rocket.v_abs * S * Cn *  distane_between_Cp_and_Cd * ((1.0 - (Rocket.alpha) * (Rocket.alpha) / 6.0) * Rocket.I_xx_inv)/*+(1.0 - Rocket7.r_y * Rocket7.r_y / 3.0) *Rocket.I_yy_inv)*/;//Taylor“WŠJ
    return ans;
}

double fall_rotate_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_z;
    return ans;
}

double fall_rotate_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 0.0;
    return ans;
}

//-----------------------------------------------------------------------------//
