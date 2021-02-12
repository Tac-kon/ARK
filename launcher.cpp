#define _CRT_SECURE_NO_WARNINGS
#include "global.h"

//--------------------------------ƒ‰ƒ“ƒ`ƒƒ[ŠŠ‘–Žž---------------------------------//
void launcher(Rocket_data *Rocket, double *t) {
	FILE *fpXYZ, *fpV;
	if ((fpXYZ = fopen("launcher_XYZ_data.txt", "w")) == NULL) {
		fprintf(stderr, "File_error_launcher_XYZ_data");
		exit(1);
	}
	if ((fpV = fopen("launcher_V_data.txt", "w")) == NULL) {
		fprintf(stderr, "File_error_launcher_V_data");
		exit(1);
	}
	
	while (Rocket->z <= (LAUNCHER_L * cos(LAUNCHER_ANGLE))) {
		rungekutta_launcher(Rocket, *t);
		if ((int)(*t / DT) % (int)(SAMPLING_T / DT) == 0) {
			fprintf(fpXYZ, "%lf\t%lf\t%lf\t%lf\t%lf\n", *t, Rocket->x, Rocket->y, Rocket->z,Rocket->gamma);
			fprintf(fpV, "%lf\t%d\t%lf\t%lf\t%lf\t%lf\n", *t,FThrustIndex,FThrustTimeStamp[FThrustIndex],FThrust[FThrustIndex],impulse(*t), Rocket->v_abs);
		}
		*t += DT;
	}
	fclose(fpXYZ);
	fclose(fpV);
}

void rungekutta_launcher(Rocket_data *Rocket, double t) {
	double k1[6], k2[6], k3[6], k4[6];
	Rocket_data RocketTemp[4] = 
	{
		*Rocket, *Rocket, *Rocket, *Rocket	
	};


	//================Step_01=========================================//
	
	k1[V_X] = DT * launcher_equation_x_2(RocketTemp[0], t);
	k1[V_Y] = DT * launcher_equation_y_2(RocketTemp[0], t);
	k1[V_Z] = DT * launcher_equation_z_2(RocketTemp[0], t);

	k1[X] = DT * launcher_equation_x_1(RocketTemp[0], t);
	k1[Y] = DT * launcher_equation_y_1(RocketTemp[0], t);
	k1[Z] = DT * launcher_equation_z_1(RocketTemp[0], t);
	
	RocketTemp[1].v_x_abs += k1[V_X] / 2.0;
	RocketTemp[1].x += k1[X] / 2.0;

	RocketTemp[1].v_y_abs += k1[V_Y] / 2.0;
	RocketTemp[1].y += k1[Y] / 2.0;

	RocketTemp[1].v_z_abs += k1[V_Z] / 2.0;
	RocketTemp[1].z += k1[Z] / 2.0;

	/*
	k1[V_X] = DT * launcher_equation_x_2( v_x_abs, t);
	k1[V_Y] = DT * launcher_equation_y_2( v_y_abs, t);
	k1[V_Z] = DT * launcher_equation_z_2( v_z_abs, t);

	k1[X] = DT * launcher_equation_x_1( x, t);
	k1[Y] = DT * launcher_equation_y_1( y, t);
	k1[Z] = DT * launcher_equation_z_1( z, t);
	*/

	//===================================================================//

	//===================Step_02=========================================//

	
	k2[V_X] = DT * launcher_equation_x_2(RocketTemp[1], t + (DT / 2.0));
	k2[V_Y] = DT * launcher_equation_y_2(RocketTemp[1], t + (DT / 2.0));
	k2[V_Z] = DT * launcher_equation_z_2(RocketTemp[1], t + (DT / 2.0));

	k2[X] = DT * launcher_equation_x_1(RocketTemp[1], t + (DT / 2.0));
	k2[Y] = DT * launcher_equation_y_1(RocketTemp[1], t + (DT / 2.0));
	k2[Z] = DT * launcher_equation_z_1(RocketTemp[1], t + (DT / 2.0));

	RocketTemp[2].v_x_abs += k2[V_X] / 2.0;
	RocketTemp[2].x += k2[X] / 2.0;

	RocketTemp[2].v_y_abs += k2[V_Y] / 2.0;
	RocketTemp[2].y += k2[Y] / 2.0;

	RocketTemp[2].v_z_abs += k2[V_Z] / 2.0;
	RocketTemp[2].z += k2[Z] / 2.0;
    
	/*
	k2[V_X] = DT * launcher_equation_x_2( v_x_abs + k1[V_X] / 2.0, t + (DT / 2.0));
	k2[V_Y] = DT * launcher_equation_y_2( v_y_abs + k1[V_Y] / 2.0, t + (DT / 2.0));
	k2[V_Z] = DT * launcher_equation_z_2( v_z_abs + k1[V_Z] / 2.0, t + (DT / 2.0));

	k2[X] = DT * launcher_equation_x_1( x + k1[X] / 2.0, t + (DT / 2.0));
	k2[Y] = DT * launcher_equation_y_1( y + k1[Y] / 2.0, t + (DT / 2.0));
	k2[Z] = DT * launcher_equation_z_1( z + k1[Z] / 2.0, t + (DT / 2.0));
	*/
	//===================================================================//

	//===================Step_03=========================================//

	
	k3[V_X] = DT * launcher_equation_x_2(RocketTemp[2], t + (DT / 2.0));
	k3[V_Y] = DT * launcher_equation_y_2(RocketTemp[2], t + (DT / 2.0));
	k3[V_Z] = DT * launcher_equation_z_2(RocketTemp[2], t + (DT / 2.0));

	k3[X] = DT * launcher_equation_x_1(RocketTemp[2], t + (DT / 2.0));
	k3[Y] = DT * launcher_equation_y_1(RocketTemp[2], t + (DT / 2.0));
	k3[Z] = DT * launcher_equation_z_1(RocketTemp[2], t + (DT / 2.0));

	RocketTemp[3].v_x_abs += k3[V_X];
	RocketTemp[3].x += k3[X];

	RocketTemp[3].v_y_abs += k3[V_Y] ;
	RocketTemp[3].y += k3[Y];

	RocketTemp[3].v_z_abs += k3[V_Z] ;
	RocketTemp[3].z += k3[Z];
	
	/*
	k3[V_X] = DT * launcher_equation_x_2( v_x_abs + k2[V_X] / 2.0, t + (DT / 2.0));
	k3[V_Y] = DT * launcher_equation_y_2( v_y_abs + k2[V_Y] / 2.0, t + (DT / 2.0));
	k3[V_Z] = DT * launcher_equation_z_2( v_z_abs + k2[V_Z] / 2.0, t + (DT / 2.0));

	k3[X] = DT * launcher_equation_x_1( x + k2[X] / 2.0, t + (DT / 2.0));
	k3[Y] = DT * launcher_equation_y_1( y + k2[Y] / 2.0, t + (DT / 2.0));
	k3[Z] = DT * launcher_equation_z_1( z + k2[Z] / 2.0, t + (DT / 2.0));
	*/
	//===================================================================//

	//===================Step_04=========================================//

	
	k4[V_X] = DT * launcher_equation_x_2(RocketTemp[3], t + DT);
	k4[V_Y] = DT * launcher_equation_y_2(RocketTemp[3], t + DT);
	k4[V_Z] = DT * launcher_equation_z_2(RocketTemp[3], t + DT);

	k4[X] = DT * launcher_equation_x_1(RocketTemp[3], t+DT);
	k4[Y] = DT * launcher_equation_y_1(RocketTemp[3], t+DT);
	k4[Z] = DT * launcher_equation_z_1(RocketTemp[3], t+DT);
	
	/*
	k4[V_X] = DT * launcher_equation_x_2( v_x_abs + k3[X] / 2.0, t + (DT / 2.0));
	k4[V_Y] = DT * launcher_equation_y_2( v_y_abs + k3[Y] / 2.0, t + (DT / 2.0));
	k4[V_Z] = DT * launcher_equation_z_2( v_z_abs + k3[Z] / 2.0, t + (DT / 2.0));

	k4[X] = DT * launcher_equation_x_1( x + k3[X] / 2.0, t + (DT / 2.0));
	k4[Y] = DT * launcher_equation_y_1( y + k3[Y] / 2.0, t + (DT / 2.0));
	k4[Z] = DT * launcher_equation_z_1( z + k3[Z] / 2.0, t + (DT / 2.0));
	*/
	//===================================================================//

	//===================Step_05=========================================//

	
	Rocket->v_x_abs += (1.0 / 6.0) * (k1[V_X] + 2.0 * k2[V_X] + 2.0 * k3[V_X] + k4[V_X]);

	Rocket->x += (1.0 / 6.0) * (k1[X] + 2.0 * k2[X] + 2.0 * k3[X] + k4[X]);

	Rocket->v_y_abs += (1.0 / 6.0) * (k1[V_Y] + 2.0 * k2[V_Y] + 2.0 * k3[V_Y] + k4[V_Y]);

	Rocket->y += (1.0 / 6.0) * (k1[Y] + 2.0 * k2[Y] + 2.0 * k3[Y] + k4[Y]);

	Rocket->v_z_abs += (1.0 / 6.0) * (k1[V_Z] + 2.0 * k2[V_Z] + 2.0 * k3[V_Z] + k4[V_Z]);

	Rocket->z += (1.0 / 6.0) * (k1[Z] + 2.0 * k2[Z] + 2.0 * k3[Z] + k4[Z]);
	
	/*
	v_x_abs += (1.0 / 6.0) * (k1[V_X] + 2.0 * k2[V_X] + 2.0 * k3[V_X] + k4[V_X]);

	x += (1.0 / 6.0) * (k1[X] + 2.0 * k2[X] + 2.0 * k3[X] + k4[X]);

	v_y_abs += (1.0 / 6.0) * (k1[V_Y] + 2.0 * k2[V_Y] + 2.0 * k3[V_Y] + k4[V_Y]);

	y += (1.0 / 6.0) * (k1[Y] + 2.0 * k2[Y] + 2.0 * k3[Y] + k4[Y]);

	v_z_abs += (1.0 / 6.0) * (k1[V_Z] + 2.0 * k2[V_Z] + 2.0 * k3[V_Z] + k4[V_Z]);

	z += (1.0 / 6.0) * (k1[Z] + 2.0 * k2[Z] + 2.0 * k3[Z] + k4[Z]);
	*/
	//======================================================================//

    //======absolute_of_velocity============================================//
    Rocket->v_abs = sqrt(Rocket->v_x_abs * Rocket->v_x_abs + Rocket->v_y_abs * Rocket->v_y_abs + Rocket->v_z_abs * Rocket->v_z_abs);
    //======================================================================//

    Rocket->gamma = atan(sqrt((Rocket->x * Rocket->x + Rocket->y * Rocket->y)) / Rocket->z);
}




double launcher_equation_x_1(Rocket_data Rocket, double t) {
	double ans;
	ans = Rocket.v_x_abs;
	return ans;
}

double launcher_equation_x_2(Rocket_data Rocket, double t) {
	double ans;
	ans = (1.0 / Rocket.m) * (impulse(t) - D(Rocket)) * sin(LAUNCHER_ANGLE) * cos(AZIMUTH_ANGLE) - G * sin(LAUNCHER_ANGLE) * cos(LAUNCHER_ANGLE) * cos(AZIMUTH_ANGLE);
	return ans;
}

double launcher_equation_y_1(Rocket_data Rocket, double t) {
	double ans;
	ans = Rocket.v_y_abs;
	return ans;
}

double launcher_equation_y_2(Rocket_data Rocket, double t) {
	double ans;
    ans = (1.0 / Rocket.m) * (impulse(t) - D(Rocket)) * sin(LAUNCHER_ANGLE) * sin(AZIMUTH_ANGLE) - G * sin(LAUNCHER_ANGLE) * cos(LAUNCHER_ANGLE) * sin(AZIMUTH_ANGLE);
	return ans;
}

double launcher_equation_z_1(Rocket_data Rocket, double t) {
	double ans;
	ans = Rocket.v_z_abs ;
	return ans;
}

double launcher_equation_z_2(Rocket_data Rocket, double t) {
	double ans;
    ans = (1.0 / Rocket.m) * (impulse(t) - D(Rocket)) * cos(LAUNCHER_ANGLE) - G + G * sin(LAUNCHER_ANGLE) * sin(LAUNCHER_ANGLE);
	return ans;
}
//-----------------------------------------------------------------------------//
