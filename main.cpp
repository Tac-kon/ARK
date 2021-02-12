#define _CRT_SECURE_NO_WARNINGS
#define GLOBAL_INSTANCE
#define GNUPLOT_PATH "gnuplot.exe"
#define GNU_DIR "gnuplot"
#define RESULT_DIR "result"

#include <direct.h>
#include <string.h>
#include <sys/stat.h>
#include "global.h"


int main() {
    int i, j, k;
    int PrevProg;
    double progress;
    char progress_bar[11] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };

    double stablity;
    double launch_clear;

    Rocket_data FTE_06;
    FILE *fp, *gp, *fpResult, *fpT, *fpXYData;
    struct stat statBuf;
    char fname[50];
    
    // ����䂪�K���͈͓��ɂ��邩�`�F�b�N
    stablity = distane_between_Cp_and_Cd * 100 / Length;
    printf("����� : %.1lf\n", stablity);
    if (stablity < 10 || stablity > 20) printf("Caution!! : �����𒲐����Ă�������\n");

    // gnuplot�̃p�X���쐬
    if (stat(GNU_DIR, &statBuf) != 0) {
        if (_mkdir(GNU_DIR) != 0 ) {
            printf("�t�H���_(gnuplot)�̍쐬�Ɏ��s���܂���\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // �V�~�����[�g���ʕۑ��p�p�X���쐬
    if (stat(RESULT_DIR, &statBuf) != 0) {
        if (_mkdir("result") != 0) {
            printf("�t�H���_(result)�̍쐬�Ɏ��s���܂���\n");
            exit(EXIT_FAILURE);
        }
    }

    // ���̓f�[�^�擾����
    // �e�L�X�g���J���邩�ǂ����`�F�b�N
    fpT = fopen("ThrustK.txt", "r");
    if (fpT == NULL) {
        printf("Error opening Thrust data file.\n");
        exit(EXIT_FAILURE);
    }

    i = 0;
    FThrustCount = 0;
    while (fscanf(fpT, "%lf %lf", &FThrustTimeStamp[i], &FThrust[i]) != EOF) {
        i += 1;
    }
    FThrustCount = i;
    fclose(fpT);

    // gnuplot�̃p�X�����������`�F�b�N
    if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {
        fprintf(stderr, "�t�@�C����������܂��� %s.", GNUPLOT_PATH);
        exit(EXIT_FAILURE);
    }
    fflush(gp); // �o�b�t�@�Ɋi�[����Ă���f�[�^��f���o���i�K�{�j
    fprintf(gp, "exit\n");
    _pclose(gp);

    // ���̓f�[�^�������ݏ���
    fpXYData = fopen("data.txt", "w");
    if (fpXYData == NULL) {
        fprintf(stderr, "File_error_XY_data");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < FThrustCount -1 ; i++) fprintf(fpXYData, "%d %lf %lf\n", i + 1, FThrustTimeStamp[i], FThrust[i]);
    fclose(fpXYData);

    // �V�~�����[�g�J�n
    printf("Simulate Start!\n");

    PrevProg = 0;

    // �����̃��[�v
    for (j = 1; j < (MAX_WIND + 1); j++) {
        FWindAbs = j;   // ����
        // �����̃��[�v
        for (k = 0; k < DEV; k++) {
            FWindDeg = 2 * M_PI * k / DEV;    // ����(�p�x)

            init(&FTE_06);

            double t = 0;
            FThrustIndex = 0;

            launcher(&FTE_06, &t);
            launch_clear = FTE_06.v_abs;
            FMaxHeight = FTE_06.z;
            thrust_inertia(&FTE_06, &FWind, &t);
            FMaxTime = t;

            // �����^�C�v�ɉ����ď�����ύX
            switch (TYPE) {

                // �^�C�v0 : 2�i�p���V���[�g�W�J������
            case 0:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                break;

                // �^�C�v1 : 1�i�ڂ̂ݓW�J������
            case 1:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                break;

                // �^�C�v2 : �e������������
            case 2:
                inertia_fall(&FTE_06, &FWind, &t);
                //�̈�(x����)
                FMaxDistanceX = MAX_DISTANCE(FMaxDistanceXt, FMaxDistanceXf);
                //�̈�(y����)
                FMaxDistanceY = MAX_DISTANCE(FMaxDistanceYt, FMaxDistanceYf);
                break;

                // �^�C�v3 : �J�P������
            case 3:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                //�̈�(x����)
                FMaxDistanceX = MAX_DISTANCE(FMaxDistanceXt, FMaxDistanceXr);
                //�̈�(y����)
                FMaxDistanceY = MAX_DISTANCE(FMaxDistanceYt, FMaxDistanceYr);
                break;

            default:
                printf("Please Select fall Type. \n");
                exit(EXIT_FAILURE);
            }

            if ((fp = fopen("XY_data.txt", "a")) == NULL) {
                fprintf(stderr, "File_error_XY_data");
                exit(EXIT_FAILURE);
            }
            fprintf(fp, "%5.1lf, %5.1lf\n", FTE_06.x, FTE_06.y);
            fclose(fp);

            /*result*/
            sprintf(fname, "%s/result [%d, %d].txt", RESULT_DIR, int(FWindAbs), int(180.0 * (FWindDeg) / M_PI));
            fpResult = fopen(fname, "w");
            if (fpResult == NULL) {
                printf("error opening %s\n", fname);
                exit(EXIT_FAILURE);
            }
            fprintf(fpResult, "Max_height    : %7.2lf [m]\n", FMaxHeight);
            fprintf(fpResult, "FMaxTime      : %7.2lf [s]\n", FMaxTime);
            fprintf(fpResult, "Launch_clear  : %7.2lf [m/s]\n", launch_clear);
            fprintf(fpResult, "Max_vercity   : %7.2lf [m/s]\n", FMaxVelocity);
            fprintf(fpResult, "PARA1tIME    : %6.2lf [m/s]\n", FPara1Time);
            fprintf(fpResult, "PARA1_2_V     : %5.1lf [m/s]\n", FTE_06.v_abs);
            fprintf(fpResult, "PARA2_V       : %5.1lf [m/s]\n", FTE_06.v_z_abs);
            fprintf(fpResult, "distance_y    : %5.1lf [m]\n", FTE_06.y);
            fprintf(fpResult, "distance_x    : %5.1lf [m]\n\n", FTE_06.x);
            fclose(fpResult);

            /*gnuplot*/
            if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {
                fprintf(stderr, "�t�@�C����������܂��� %s.", GNUPLOT_PATH);
                exit(EXIT_FAILURE);
            }
            fprintf(gp, "set grid\n");
            fprintf(gp, "set xl'distance_x (m)'\n");
            fprintf(gp, "set yl'distance_y (m)'\n");
            fprintf(gp, "set zl'height (m)'\n");
            fprintf(gp, "set zlabel rotate by 90\n");
            fprintf(gp, "set xrange[-abs(%5.1lf)-50.:abs(%5.1lf)+50]\n", FMaxDistanceX, FMaxDistanceX);
            fprintf(gp, "set yrange[-abs(%5.1lf)-50.:abs(%5.1lf)+50]\n", FMaxDistanceY, FMaxDistanceY);
            fprintf(gp, "set zrange[0.:%7.2lf+50]\n", FMaxHeight);
            fprintf(gp, "set ticslevel 0\n");
            fprintf(gp, "set title'Flight-Trajectory [%d, %d]'\n", int(FWindAbs), int(180.0 * (FWindDeg) / M_PI));
            fprintf(gp, "set terminal png\n");
            fprintf(gp, "set output '%s/Flight-Trajectory [%d, %d].jpg'\n", GNU_DIR, int(FWindAbs), int(180.0 * (FWindDeg) / M_PI));
            if (TYPE == 2) {
                fprintf(gp, "splot 'launcher_XYZ_data.txt' using($2):($3):($4) , 'thrust_inertia_XYZ_data.txt' using($2):($3):($4), 'inertiafall_XYZ_data.txt' using($2):($3):($4) \n");
            }
            else {
                fprintf(gp, "splot 'launcher_XYZ_data.txt' using($2):($3):($4) , 'thrust_inertia_XYZ_data.txt' using($2):($3):($4), 'inertiafall_XYZ_data.txt' using($2):($3):($4), 'recovery_XYZ_data.txt' using($2):($3):($4) \n");
            }
            fflush(gp); // �o�b�t�@�Ɋi�[����Ă���f�[�^��f���o���i�K�{�j
            fprintf(gp, "exit\n");
            _pclose(gp);

            // �i�s�󋵕\��
            progress = 100 * (DEV * (int(FWindAbs) - 1) + k + 1) / (MAX_WIND * DEV);
            if (int(progress / 10) - PrevProg > 0)
            {
                PrevProg += 1;
                for (i = 0; i < PrevProg; i++) {
                    progress_bar[i] = '#';
                }
                for (i = PrevProg; i < 10; i++) {
                    progress_bar[i] = '\0';
                }
            }
            printf("\r|%s|  �i�s�� : %.0lf ��", progress_bar, progress);
            fflush(stdout);
        }
    }
    printf("\nSimulate finished!\n");

    return 0;
}