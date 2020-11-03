/*内部パラメータ*/


/*ロケットに関する構造体*/
typedef struct Rocket_data {
	double m;
	double m_;

	double x;
	double y;			
	double z;	//高度と一致

	double r;	//一時的な値
	double gamma;
	double psi;

	double alpha;
	//double beta;
	double delta;


	double v_abs;		//機体の最大速さ
	double v_abs_0;		//最大対気速さx
	double v_r_abs;
	double v_gamma_abs;
	double v_psi_abs;
	double v_x_abs;
	double v_y_abs;
	double v_z_abs;
	double v_y;
	double v_x;
	double v_z;

	double I_xx_inv;
	double I_xy_inv;
	double I_xz_inv;
	double I_yy_inv;
	double I_yz_inv;
	double I_zz_inv;

	double r_x;
	double r_y;
	double r_z;
	double r_v_x;
	double r_v_y;
	double r_v_z;



	double v0_x;	//弾道飛行下降時の初期速度
	double v0_z;	//弾道飛行下降時の初期速度


	double ans;

} Rocket_data;


/*風に関する構造体*/
typedef struct wind_data{
	double w_abs;
	double w_x;
	double w_y;
	double w_gamma;
	double w_phi;
	double phi;
}Wind_data;