#include "global.h"


//----------------------------------その他---------------------------------------//

//初期化
void init(Rocket_data *Rocket) {
    /*
    x,y,zはシミュレーションを行う空間にはった3次元デカルト座標
    特にx-y平面を地面にとり、z成分は高度と一致するようにとった

    r_x,r_y.r_zは、ロケットにはった回転を表す座標
    特に、r_z方向はロケットの機軸方向に一致する
    r_x,r_yはピッチ・ヨー
    r_zはロール
    */


    Rocket->m = M;
    Rocket->m_ = M_;
    Rocket->x = 0.0;
    Rocket->y = 0.0;
    Rocket->z = 0.0;

    Rocket->r = 0.0;
    Rocket->gamma = LAUNCHER_ANGLE;
    Rocket->psi = AZIMUTH_ANGLE;

    Rocket->alpha = 0.0;
    //Rocket->beta = 0.0;
    Rocket->v_abs = 0.0;
    Rocket->v_abs_0 = 0.0;
    Rocket->v_r_abs = 0.0;
    Rocket->v_gamma_abs = 0.0;
    Rocket->v_psi_abs = 0.0;

    Rocket->v_y_abs = 0.0;
    Rocket->v_x_abs = 0.0;
    Rocket->v_z_abs = 0.0;

    Rocket->r_x = 0.0;
    Rocket->r_y = 0.0;
    Rocket->r_z = 0.0;
    Rocket->r_v_x = 0.0;
    Rocket->r_v_y = 0.0;
    Rocket->r_v_z = 0.0;

    Rocket->I_xx_inv = (I_yy * I_zz - I_yz * I_yz) / determinant_of_I();
    Rocket->I_xy_inv = (I_xz * I_yz - I_xy * I_zz) / determinant_of_I();
    Rocket->I_xz_inv = (I_xy * I_yz - I_xy * I_zz) / determinant_of_I();
    Rocket->I_yy_inv = (I_xx * I_zz - I_xz * I_xz) / determinant_of_I();
    Rocket->I_yz_inv = (I_xz * I_xy - I_xx * I_xz) / determinant_of_I();
    Rocket->I_zz_inv = (I_xx * I_yy - I_xy * I_xy) / determinant_of_I();


}

void normalization_wind(Wind_data *Wind, Rocket_data Rocket) {
    //Wind->phi = WIND_DEG;
    Wind->w_abs = wind_height(Rocket);
    Wind->w_x = -Wind->w_abs * cos(WIND_DEG);
    Wind->w_y = -Wind->w_abs * sin(WIND_DEG);
    //Wind->w_gamma = Wind->w_x * cos(Rocket.psi) + Wind->w_y * sin(Rocket.psi);
    //Wind->w_phi = -1.0 * Wind->w_x * sin(Rocket.psi) + Wind->w_y * cos(Rocket.psi);
}

//高度-風速関数
double wind_height(Rocket_data Rocket) {
    double ans;
    ans = (WIND_ABS * pow((Rocket.z / 10.0), 1.0 / 7.0));
    return ans;
}


//空気抵抗計算関数(風の影響なし)
double D(Rocket_data Rocket) {
    double ans, ans_1;
    ans_1 = 0.5 * rho(Rocket) * S * (Rocket.v_abs * Rocket.v_abs  * Mach_Cd(Rocket) + Rocket.v_abs * 68);//慣性抵抗
    ans = ans_1;
    return ans;
}

//空気抵抗計算関数(風の影響あり)
double D_0(Rocket_data Rocket, Wind_data Wind) {
    double ans, ans_1;
    ans_1 = 0.5 * rho(Rocket) * S * (Rocket.v_abs_0 * Rocket.v_abs_0  * Mach_Cd(Rocket) + Rocket.v_abs_0 * 68);//慣性抵抗
    ans = ans_1;
    return ans;
}

//法線力計算関数(風の影響あり)
double N_0(Rocket_data Rocket, Wind_data Wind) {
    double ans;
    ans = 0.5 * rho(Rocket) * Rocket.v_abs_0 * Rocket.v_abs_0 * S * Cn * abs(Rocket.alpha);
    return ans;
}

//Mach_Cd_curve(
double Mach_Cd(Rocket_data Rocket) {
    double ans;
    ans = Cd;
    /*0.429068
    + (-0.0092415) * Rocket.v_abs
    + (0.000318424) * Rocket.v_abs * Rocket.v_abs
    + (-0.00000498557) * Rocket.v_abs * Rocket.v_abs * Rocket.v_abs
    + (0.0000000288) * Rocket.v_abs * Rocket.v_abs * Rocket.v_abs * Rocket.v_abs;*/
    return ans;
}

//空気密度計算関数
double rho(Rocket_data Rocket) {
    double ta = T0 - (0.6 / 100.0) * Rocket.z + 273.15;
    return ((Pa_AVE / (Ra * ta)) * exp((-1 * G * Rocket.z) / (Ra * Ta_AVE)));
}

//慣性テンソルの行列式
double determinant_of_I() {
    double ans;
    ans = I_xx * I_yy * I_zz + I_xy * I_yz * I_xz + I_xz * I_xy * I_yz - I_xx * I_yz * I_yz - I_xz * I_yy * I_xz - I_xy * I_xy * I_zz;
    return ans;
}

//r,gamma,psi → x,y,z に変換
void normalization(Rocket_data *Rocket) {
    Rocket->z = Rocket->r * cos(Rocket->gamma);
    Rocket->y = Rocket->r * sin(Rocket->gamma) * sin(Rocket->psi);
    Rocket->x = Rocket->r * sin(Rocket->gamma) * cos(Rocket->psi);
}
//-----------------------------------------------------------------------------//
