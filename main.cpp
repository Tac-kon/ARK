#define _CRT_SECURE_NO_WARNINGS
#define GLOBAL_INSTANCE
#define GNUPLOT_PATH "C:/PROGRA~1/gnuplot/bin/gnuplot.exe"

#include<direct.h>
#include<string.h>
#include "global.h"


int main() {
    int i, j, k;
    int PrevProg;
    double stablity;
    double launch_clear;
    Rocket_data FTE_06;
    FILE* fp, * gp, * hp, *fpT, *fpXYData;
    char fname[50];
    double progress;
    char progress_bar[11] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };

    // 安定比が適正範囲内にあるかチェック
    stablity = distane_between_Cp_and_Cd * 100 / Length;
    printf("安定比 : %.1lf\n", stablity);
    if (stablity < 10 || stablity > 20) printf("Caution!! : 安定比を調整してください\n");

    // gnuplotのパスを作成
    _mkdir("gnuplot");
    //if (_mkdir("gnuplot") > 0 ) {
    //    printf("フォルダ(gnuplot)の作成に失敗しました\n");
    //    exit(EXIT_FAILURE);
    //}

    // シミュレート結果保存用パスを作成
    _mkdir("result");
    //if (_mkdir("result") > 0) {
    //    printf("フォルダ(result)の作成に失敗しました\n");
    //    exit(EXIT_FAILURE);
    //}

    // 推力データ取得処理
    // テキストが開けるかどうかチェック
    fpT = fopen("ThrustK.txt", "r");
    if (fpT == NULL) {
        printf("Error opening Thrust data file.\n");
        exit(EXIT_FAILURE);
    }

    i = 0;
    ThrustCount = 0;
    while (fscanf(fpT, "%lf %lf", &ThrustTimeStamp[i], &Thrust[i]) != EOF) {
        i += 1;
    }
    ThrustCount = i;
    fclose(fpT);

    // gnuplotのパスが正しいかチェック
    if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {
        fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
        exit(EXIT_FAILURE);
    }
    fflush(gp); // バッファに格納されているデータを吐き出す（必須）
    fprintf(gp, "exit\n");
    _pclose(gp);

    // 推力データ書き込み処理
    fpXYData = fopen("data.txt", "w");
    if (fpXYData == NULL) {
        fprintf(stderr, "File_error_XY_data");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < ThrustCount -1 ; i++) fprintf(fpXYData, "%d %lf %lf\n", i + 1, ThrustTimeStamp[i], Thrust[i]);
    fclose(fpXYData);

    // シミュレート開始
    printf("Simulate Start!\n");

    PrevProg = 0;

    // 風速のループ
    for (j = 1; j < (MAX_WIND + 1); j++) {
        WIND_ABS = j;   // 風速
        // 風向のループ
        for (k = 0; k < DEV; k++) {
            WIND_DEG = 2 * M_PI * k / DEV;    // 風向(角度)

            init(&FTE_06);

            double t = 0;
            ThrustIndex = 0;

            launcher(&FTE_06, &t);
            launch_clear = FTE_06.v_abs;
            g_max_height = FTE_06.z;
            thrust_inertia(&FTE_06, &Wind, &t);
            max_time = t;

            // 落下タイプに応じて処理を変更
            switch (TYPE) {

                // タイプ0 : 2段パラシュート展開時処理
            case 0:
                inertia_fall(&FTE_06, &Wind, &t);
                PARA_1_TIME = t;
                recovery(&FTE_06, &Wind, &t);
                break;

                // タイプ1 : 1段目のみ展開時処理
            case 1:
                inertia_fall(&FTE_06, &Wind, &t);
                PARA_1_TIME = t;
                recovery(&FTE_06, &Wind, &t);
                break;

                // タイプ2 : 弾道落下時処理
            case 2:
                inertia_fall(&FTE_06, &Wind, &t);
                //領域(x方向)
                max_distance_x = MAX_DISTANCE(max_distance_x_t, max_distance_x_f);
                //領域(y方向)
                max_distance_y = MAX_DISTANCE(max_distance_y_t, max_distance_y_f);
                break;

                // タイプ3 : 開傘時処理
            case 3:
                inertia_fall(&FTE_06, &Wind, &t);
                PARA_1_TIME = t;
                recovery(&FTE_06, &Wind, &t);
                //領域(x方向)
                max_distance_x = MAX_DISTANCE(max_distance_x_t, max_distance_x_r);
                //領域(y方向)
                max_distance_y = MAX_DISTANCE(max_distance_y_t, max_distance_y_r);
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
            sprintf(fname, "./result/result [%d, %d].txt", int(WIND_ABS), int(180.0 * (WIND_DEG) / M_PI));
            hp = fopen(fname, "w");
            if (hp == NULL) {
                printf("error opening %s\n", fname);
                exit(EXIT_FAILURE);
            }
            fprintf(hp, "Max_height    : %7.2lf [m]\n", g_max_height);
            fprintf(hp, "Max_time      : %7.2lf [s]\n", max_time);
            fprintf(hp, "Launch_clear  : %7.2lf [m/s]\n", launch_clear);
            fprintf(hp, "Max_vercity   : %7.2lf [m/s]\n", max_velocity);
            fprintf(hp, "PARA1_TIME    : %6.2lf [m/s]\n", PARA_1_TIME);
            fprintf(hp, "PARA1_2_V     : %5.1lf [m/s]\n", FTE_06.v_abs);
            fprintf(hp, "PARA2_V       : %5.1lf [m/s]\n", FTE_06.v_z_abs);
            fprintf(hp, "distance_y    : %5.1lf [m]\n", FTE_06.y);
            fprintf(hp, "distance_x    : %5.1lf [m]\n\n", FTE_06.x);
            fclose(hp);

            /*gnuplot*/
            if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {
                fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
                exit(EXIT_FAILURE);
            }
            fprintf(gp, "set grid\n");
            fprintf(gp, "set xl'distance_x (m)'\n");
            fprintf(gp, "set yl'distance_y (m)'\n");
            fprintf(gp, "set zl'height (m)'\n");
            fprintf(gp, "set zlabel rotate by 90\n");
            fprintf(gp, "set xrange[-abs(%5.1lf)-50.:abs(%5.1lf)+50]\n", max_distance_x, max_distance_x);
            fprintf(gp, "set yrange[-abs(%5.1lf)-50.:abs(%5.1lf)+50]\n", max_distance_y, max_distance_y);
            fprintf(gp, "set zrange[0.:%7.2lf+50]\n", g_max_height);
            fprintf(gp, "set ticslevel 0\n");
            fprintf(gp, "set title'Flight-Trajectory [%d, %d]'\n", int(WIND_ABS), int(180.0 * (WIND_DEG) / M_PI));
            fprintf(gp, "set terminal png\n");
            fprintf(gp, "set output './gnuplot/Flight-Trajectory [%d, %d].jpg'\n", int(WIND_ABS), int(180.0 * (WIND_DEG) / M_PI));
            if (TYPE == 2) {
                fprintf(gp, "splot 'launcher_XYZ_data.txt' using($2):($3):($4) , 'thrust_inertia_XYZ_data.txt' using($2):($3):($4), 'inertia_fall_XYZ_data.txt' using($2):($3):($4) \n");
            }
            else {
                fprintf(gp, "splot 'launcher_XYZ_data.txt' using($2):($3):($4) , 'thrust_inertia_XYZ_data.txt' using($2):($3):($4), 'inertia_fall_XYZ_data.txt' using($2):($3):($4), 'recovery_XYZ_data.txt' using($2):($3):($4) \n");
            }
            fflush(gp); // バッファに格納されているデータを吐き出す（必須）
            fprintf(gp, "exit\n");
            _pclose(gp);

            // 進行状況表示
            progress = 100 * (DEV * (int(WIND_ABS) - 1) + k + 1) / (MAX_WIND * DEV);
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
            printf("\r|%s|  進行状況 : %.0lf ％", progress_bar, progress);
            fflush(stdout);
        }
    }
    printf("\nSimulate finished!\n");

    return 0;
}