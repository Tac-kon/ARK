#define _CRT_SECURE_NO_WARNINGS
#include "global.h"

//--------------------------„—Í‹y‚ÑŠµ«”òsŽž(ã¸Žž)-----------------------------//
void thrust_inertia(Rocket_data *Rocket, Wind_data *Wind, double *t) {
	FILE *fp1, *fp2/*, *fp3*/;
	if((fp1 = fopen("thrust_inertia_XYZ_data.txt", "w")) == NULL) {
		fprintf(stderr, "File_error_thrust_inertia_XYZ_data");
		exit(1);
	}
	if((fp2 = fopen("thrust_inertia_V_data.txt", "w")) == NULL) {
		fprintf(stderr, "File_error_thrust_inertia_V_data");
		exit(1);
	}

    normalization_wind(Wind, *Rocket);

	while ( Rocket->z  >= (g_max_height)) {
		//printf("max_distance_x_t = %lf\n", max_distance_x_t);
		//printf("gamma = %lf\n", Rocket->gamma * 180.0 / M_PI);
		//printf("r_x = %lf\n", Rocket->r_x * 180.0 / M_PI);
		//printf("alpha = %lf\n", Rocket->alpha * 180.0 / M_PI);
		rungekutta_thrust_inertia_rotate(Rocket, Wind, *t);
		if((int)(*t / DT) % (int)(SAMPLING_T / DT) == 0) {
			fprintf(fp1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", *t, Rocket->x, Rocket->y, Rocket->z,Rocket->psi,Rocket->alpha, Rocket->r);
			fprintf(fp2, "%lf\t%lf\t%lf\t%lf\t%lf\tf\n", *t, Rocket->v_abs, Rocket->v_x_abs, Rocket->v_y_abs, Rocket->v_z_abs);
		}
		
		*t += DT; 
		//printf("t=%lf,z=%lf,v=%lf\n", *t,Rocket->z, Rocket->v_abs);
		g_max_height = MAX_HEIGHT(Rocket->z, g_max_height);
		//‘¬‚³‚ÌÅ‘å’l
		if (Rocket->v_abs >= (max_velocity)) {
			max_velocity = MAX_VELOCITY(Rocket->v_abs, max_velocity);
		}
		//—Ìˆæ(x•ûŒü)
		if (abs(Rocket->x) >= (max_distance_x_t)) {
			max_distance_x_t = MAX_DISTANCE(abs(Rocket->x), max_distance_x_t);
		}
		//—Ìˆæ(y•ûŒü)
		if (abs(Rocket->y) >= (max_distance_y_t)) {
			max_distance_y_t = MAX_DISTANCE(abs(Rocket->y), max_distance_y_t);
		}

	} 

	
	fclose(fp1);
	fclose(fp2);
	
}



void rungekutta_thrust_inertia_rotate(Rocket_data *Rocket, Wind_data *Wind, double t) {
    double k1[12], k2[12], k3[12], k4[12];
    Rocket_data Rocket_temp[4] =
    {
        *Rocket, *Rocket, *Rocket, *Rocket
    };

    Wind_data Wind_temp[4] =
    {
        *Wind, *Wind, *Wind, *Wind
    };

    //===================Step_01=========================================//

    //•Ài‰^“®
    k1[V_X] = DT * thrust_inertia_equation_x_2(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Y] = DT * thrust_inertia_equation_y_2(Rocket_temp[0], Wind_temp[0], t);
    k1[V_Z] = DT * thrust_inertia_equation_z_2(Rocket_temp[0], Wind_temp[0], t);

    k1[X] = DT * thrust_inertia_equation_x_1(Rocket_temp[0], Wind_temp[0], t);
    k1[Y] = DT * thrust_inertia_equation_y_1(Rocket_temp[0], Wind_temp[0], t);
    k1[Z] = DT * thrust_inertia_equation_z_1(Rocket_temp[0], Wind_temp[0], t);


    //‰ñ“]‰^“®
    k1[R_V_X] = DT * thrust_rotate_equation_x_2(Rocket_temp[0], Wind_temp[0], t);
    k1[R_V_Y] = DT * thrust_rotate_equation_y_2(Rocket_temp[0], Wind_temp[0], t);
    k1[R_V_Z] = DT * thrust_rotate_equation_z_2(Rocket_temp[0], Wind_temp[0], t);

    k1[R_X] = DT * thrust_rotate_equation_x_1(Rocket_temp[0], Wind_temp[0], t);
    k1[R_Y] = DT * thrust_rotate_equation_y_1(Rocket_temp[0], Wind_temp[0], t);
    k1[R_Z] = DT * thrust_rotate_equation_z_1(Rocket_temp[0], Wind_temp[0], t);



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


    //normalization(&Rocket_temp[0][1]);//’¼Œð‰»
    //normalization(&Rocket_temp[1][1]);
    //normalization(&Rocket_temp[2][1]);


    //===================================================================//

    //===================Step_02=========================================//

    //•Ài‰^“®
    k2[V_X] = DT * thrust_inertia_equation_x_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[V_Y] = DT * thrust_inertia_equation_y_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[V_Z] = DT * thrust_inertia_equation_z_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));

    k2[X] = DT * thrust_inertia_equation_x_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[Y] = DT * thrust_inertia_equation_y_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[Z] = DT * thrust_inertia_equation_z_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));


    //‰ñ“]‰^“®
    k2[R_V_X] = DT * thrust_rotate_equation_x_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[R_V_Y] = DT * thrust_rotate_equation_y_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[R_V_Z] = DT * thrust_rotate_equation_z_2(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));

    k2[R_X] = DT * thrust_rotate_equation_x_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[R_Y] = DT * thrust_rotate_equation_y_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));
    k2[R_Z] = DT * thrust_rotate_equation_z_1(Rocket_temp[1], Wind_temp[1], t + (DT / 2.0));

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

    //normalization(&Rocket_temp[0][2]);//’¼Œð‰»
    //normalization(&Rocket_temp[1][2]);
    //normalization(&Rocket_temp[2][2]);

    //===================================================================//

    //===================Step_03=========================================//

    //•Ài‰^“®
    k3[V_X] = DT * thrust_inertia_equation_x_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[V_Y] = DT * thrust_inertia_equation_y_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[V_Z] = DT * thrust_inertia_equation_z_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));

    k3[X] = DT * thrust_inertia_equation_x_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[Y] = DT * thrust_inertia_equation_y_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[Z] = DT * thrust_inertia_equation_z_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));


    //‰ñ“]‰^“®
    k3[R_V_X] = DT * thrust_rotate_equation_x_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[R_V_Y] = DT * thrust_rotate_equation_y_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[R_V_Z] = DT * thrust_rotate_equation_z_2(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));

    k3[R_X] = DT * thrust_rotate_equation_x_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[R_Y] = DT * thrust_rotate_equation_y_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));
    k3[R_Z] = DT * thrust_rotate_equation_z_1(Rocket_temp[2], Wind_temp[2], t + (DT / 2.0));

    Rocket_temp[3].v_x_abs += k3[V_X];
    Rocket_temp[3].x += k3[X];

    Rocket_temp[3].v_y_abs += k3[V_Y] / 2.0;
    Rocket_temp[3].y += k3[Y];

    Rocket_temp[2].v_z_abs += k3[V_Z] / 2.0;
    Rocket_temp[3].z += k3[Z];

    Rocket_temp[3].r_v_x += k3[V_X];
    Rocket_temp[3].r_x += k3[X];

    Rocket_temp[3].r_v_y += k3[V_Y] / 2.0;
    Rocket_temp[3].r_y += k3[Y];

    Rocket_temp[2].r_v_z += k3[V_Z] / 2.0;
    Rocket_temp[3].r_z += k3[Z];


    //normalization(&Rocket_temp[0][3]);//’¼Œð‰»
    //normalization(&Rocket_temp[1][3]);
    //normalization(&Rocket_temp[2][3]);

    //===================================================================//

    //===================Step_04=========================================//

    //•Ài‰^“®
    k4[V_X] = DT * thrust_inertia_equation_x_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Y] = DT * thrust_inertia_equation_y_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[V_Z] = DT * thrust_inertia_equation_z_2(Rocket_temp[3], Wind_temp[3], t + DT);

    k4[X] = DT * thrust_inertia_equation_x_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Y] = DT * thrust_inertia_equation_y_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[Z] = DT * thrust_inertia_equation_z_1(Rocket_temp[3], Wind_temp[3], t + DT);


    //‰ñ“]‰^“®
    k4[R_V_X] = DT * thrust_rotate_equation_x_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_V_Y] = DT * thrust_rotate_equation_y_2(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_V_Z] = DT * thrust_rotate_equation_z_2(Rocket_temp[3], Wind_temp[3], t + DT);

    k4[R_X] = DT * thrust_rotate_equation_x_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_Y] = DT * thrust_rotate_equation_y_1(Rocket_temp[3], Wind_temp[3], t + DT);
    k4[R_Z] = DT * thrust_rotate_equation_z_1(Rocket_temp[3], Wind_temp[3], t + DT);


    //===================================================================//

    //===================Step_05=========================================//

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


    //======================================================================//

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
    Rocket->alpha = -Rocket->delta + (LAUNCHER_ANGLE + Rocket->r_x);

    //normalization(Rocket);

    //======alpha===========================================================//
    //Rocket->alpha = ALPHA_ERROR * (atan2(Rocket->v_abs * cos(Rocket->gamma), Rocket->v_abs * sin(Rocket->gamma) - Wind.w_gamma) - Rocket->gamma);
    //======================================================================//

    //=======beta===========================================================//
    //Rocket->beta = ALPHA_ERROR * (atan2(-1.0 * Wind.w_phi, sqrt(((Rocket->v_abs * sin(Rocket->gamma)) - Wind.w_gamma) * ((Rocket->v_abs * sin(Rocket->gamma)) - Wind.w_gamma) + (Rocket->v_abs * cos(Rocket->gamma) * Rocket->v_abs * cos(Rocket->gamma)))));
    //======================================================================//
}


//‰^“®•û’öŽ®‚ð•ÏX(8/19)
double thrust_inertia_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_x_abs;
    return ans;
}

double thrust_inertia_equation_x_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 1.0 / Rocket.m * (impulse(t) * sin(LAUNCHER_ANGLE + Rocket.r_x) * cos(AZIMUTH_ANGLE + Rocket.r_y) - D_0(Rocket, Wind) /** cos(Rocket.alpha) */* sin(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x/* - Rocket7.r_y*/) * cos(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/));// -N_0(Rocket, Wind) * sin(LAUNCHER_ANGLE + Rocket.r_x) / Rocket.m * cos(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/);
    return ans;
}

double thrust_inertia_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_y_abs;
    return ans;
}

double thrust_inertia_equation_y_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 1.0 / Rocket.m * (impulse(t) * sin(LAUNCHER_ANGLE + Rocket.r_x) * sin(AZIMUTH_ANGLE + Rocket.r_y) - D_0(Rocket, Wind) /* * cos(Rocket.alpha)*/ * sin(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x/* - Rocket7.r_y*/) * sin(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/));// -N_0(Rocket, Wind) * sin(LAUNCHER_ANGLE + Rocket.r_x) / Rocket.m * sin(AZIMUTH_ANGLE + Rocket.r_y /*+ Rocket7.r_y*/);
    return ans;
}

double thrust_inertia_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.v_z_abs;
    return ans;
}

double thrust_inertia_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    //if(impulse(t)>0&&t<1) printf("%lf %lf\n",t, impulse(t));
    ans = 1.0 / Rocket.m * (impulse(t) * cos(LAUNCHER_ANGLE + Rocket.r_x) - D_0(Rocket, Wind) /** cos(Rocket.alpha)*/ * cos(Rocket.alpha + LAUNCHER_ANGLE + Rocket.r_x)) - G;// +N_0(Rocket, Wind) * sin(LAUNCHER_ANGLE + Rocket.r_x) / Rocket.m;
    return ans;
}

double thrust_rotate_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_x;
    return ans;

}

double thrust_rotate_equation_x_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = D_0(Rocket, Wind) * sin(Rocket.alpha) * Rocket.I_xx_inv * distane_between_Cp_and_Cd;// ((1.0 - (Rocket.alpha) * (Rocket.alpha) / 6.0) * Rocket.I_xx_inv /*+ (1.0 - Rocket7.r_y * Rocket7.r_y / 3.0) *Rocket.I_xy_inv*/);//Taylor“WŠJ
    return ans;
}

double thrust_rotate_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_y;
    return ans;
}

double thrust_rotate_equation_y_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    //ans =  (-0.5) * rho(Rocket) * Rocket.v_abs * Rocket.v_abs * S * Cn * (/*(1.0 - Rocket6.r_x * Rocket6.r_x / 3.0) * Rocket.I_xy_inv*/ +(1.0 - Rocket.r_y * Rocket.r_y / 3.0) *Rocket.I_yy_inv) * distane_between_Cp_and_Cd;//Taylor“WŠJ
    ans = D_0(Rocket, Wind) * sin(Rocket.r_y) * Rocket.I_yy_inv * distane_between_Cp_and_Cd;// (-0.5) * rho(Rocket) * Rocket.v_abs * Rocket.v_abs * S * Cn *  distane_between_Cp_and_Cd * ((1.0 - (Rocket.alpha) * (Rocket.alpha) / 6.0) * Rocket.I_xx_inv)/*+(1.0 - Rocket7.r_y * Rocket7.r_y / 3.0) *Rocket.I_yy_inv)*/;//Taylor“WŠJ
    return ans;
}

double thrust_rotate_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = Rocket.r_v_z;
    return ans;
}

double thrust_rotate_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t) {
    double ans;
    ans = 0.0;
    return ans;
}
//-----------------------------------------------------------------------------//

