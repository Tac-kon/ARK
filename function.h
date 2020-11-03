

/*ランチャー滑走時　->　launcher.cpp　*/
GLOBAL void launcher(Rocket_data *Rocket, double *t);
GLOBAL void rungekutta_launcher(Rocket_data *Rocket, double t);
GLOBAL double launcher_equation_x_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_x_2(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_y_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_y_2(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_z_1(Rocket_data Rocket, double t);
GLOBAL double launcher_equation_z_2(Rocket_data Rocket, double t);

/*推力及び慣性飛行時　->　thrust_inertia.cpp　*/
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


/*慣性飛行時(下降時)　->　inertia_fall.cpp　*/
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

/*パラシュート展開時　->　recovery.cpp　*/
GLOBAL void recovery(Rocket_data *Rocket,Wind_data *Wind, double *t);
GLOBAL void rungekutta_recovery(Rocket_data *Rocket,Wind_data *Wind, double t);
GLOBAL double recovery_v_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_x(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_v_z(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_z(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_v_y(Rocket_data Rocket, Wind_data Wind, double t);
GLOBAL double recovery_y(Rocket_data Rocket, Wind_data Wind, double t);


/*初期化&変換　->　init.cpp　*/
GLOBAL void init(Rocket_data *Rocket);
GLOBAL void normalization(Rocket_data *Rocket);
GLOBAL void normalization_wind(Wind_data *Wind, Rocket_data Rocket);


/*インパルス計算関数*/
GLOBAL double impulse(double t);


/*空気抵抗計算関数*/
GLOBAL double D(Rocket_data Rocket);
GLOBAL double D_0(Rocket_data Rocket, Wind_data Wind);

/*法線力計算関数*/
GLOBAL double N_0(Rocket_data Rocket, Wind_data Wind);


/*空気密度計算関数*/
GLOBAL double rho(Rocket_data Rocket);


/*高度-風速関数*/
GLOBAL double wind_height(Rocket_data Rocket);

/*対称行列の行列式*/
GLOBAL double determinant_of_I();


/*空気抵抗係数*/
GLOBAL double Mach_Cd(Rocket_data Rocket);