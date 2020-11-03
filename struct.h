/*�����p�����[�^*/


/*���P�b�g�Ɋւ���\����*/
typedef struct Rocket_data {
	double m;
	double m_;

	double x;
	double y;			
	double z;	//���x�ƈ�v

	double r;	//�ꎞ�I�Ȓl
	double gamma;
	double psi;

	double alpha;
	//double beta;
	double delta;


	double v_abs;		//�@�̂̍ő呬��
	double v_abs_0;		//�ő�΋C����x
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



	double v0_x;	//�e����s���~���̏������x
	double v0_z;	//�e����s���~���̏������x


	double ans;

} Rocket_data;


/*���Ɋւ���\����*/
typedef struct wind_data{
	double w_abs;
	double w_x;
	double w_y;
	double w_gamma;
	double w_phi;
	double phi;
}Wind_data;