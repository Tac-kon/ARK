

/*�����`���[�������@->�@launcher.cpp�@*/
GLOBAL void launcher(Rocket_data *Rocket, double *t);
GLOBAL void rungekutta_launcher(Rocket_data *Rocket, double t);
GLOBAL double launcher_equation_x_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_x_2(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_y_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_y_2(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_z_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_z_2(Rocket_data Rocket, double t);

/*���͋y�ъ�����s���@->�@thrust_inertia.cpp�@*/
GLOBAL void thrust_inertia(Rocket_data *Rocket, Wind_data *Wind, double *t);
GLOBAL void rungekutta_thrust_inertia_rotate(Rocket_data *Rocket, Wind_data *Wind, double t);
GLOBAL double thrust_inertia_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_inertia_equation_x_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_inertia_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_inertia_equation_y_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_inertia_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_inertia_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t);

GLOBAL double thrust_rotate_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_rotate_equation_x_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_rotate_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_rotate_equation_y_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_rotate_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double thrust_rotate_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t);


/*������s��(���~��)�@->�@inertia_fall.cpp�@*/
GLOBAL void inertia_fall(Rocket_data *Rocket, Wind_data *Wind, double *t);
GLOBAL void rungekutta_inertia_fall(Rocket_data *Rocket, Wind_data *Wind, double t);
GLOBAL double inertia_fall_equation_v_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double inertia_fall_equation_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double inertia_fall_equation_v_y(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double inertia_fall_equation_y(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double inertia_fall_equation_v_z(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double inertia_fall_equation_z(Rocket_data Rocket, Wind_data Wind, double t);

GLOBAL double fall_rotate_equation_x_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double fall_rotate_equation_x_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double fall_rotate_equation_y_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double fall_rotate_equation_y_2(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double fall_rotate_equation_z_1(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double fall_rotate_equation_z_2(Rocket_data Rocket, Wind_data Wind, double t);

/*�p���V���[�g�W�J���@->�@recovery.cpp�@*/
GLOBAL void recovery(Rocket_data *Rocket,Wind_data *Wind, double *t);
GLOBAL void rungekutta_recovery(Rocket_data *Rocket,Wind_data *Wind, double t);
GLOBAL double recovery_v_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_v_z(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_z(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_v_y(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_y(Rocket_data Rocket, Wind_data Wind, double t);


/*������&�ϊ��@->�@init.cpp�@*/
GLOBAL void init(Rocket_data *Rocket);
GLOBAL void normalization(Rocket_data *Rocket);
GLOBAL void normalization_wind(Wind_data *Wind, Rocket_data Rocket);


/*�C���p���X�v�Z�֐�*/
GLOBAL double impulse(double t);


/*��C��R�v�Z�֐�*/
GLOBAL double D(Rocket_data Rocket);
GLOBAL double D_0(Rocket_data Rocket, Wind_data Wind);

/*�@���͌v�Z�֐�*/
GLOBAL double N_0(Rocket_data Rocket, Wind_data Wind);


/*��C���x�v�Z�֐�*/
GLOBAL double rho(Rocket_data Rocket);


/*���x-�����֐�*/
GLOBAL double wind_height(Rocket_data Rocket);

/*�Ώ̍s��̍s��*/
GLOBAL double determinant_of_I();


/*��C��R�W��*/
GLOBAL double Mach_Cd(Rocket_data Rocket);