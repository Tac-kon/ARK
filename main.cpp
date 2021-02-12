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
    
    // 安定比が適正範囲内にあるかチェック
    stablity = distane_between_Cp_and_Cd * 100 / Length;
    printf("安定比 : %.1lf\n", stablity);
    if (stablity < 10 || stablity > 20) printf("Caution!! : 安定比を調整してください\n");

    // gnuplotのパスを作成
    if (stat(GNU_DIR, &statBuf) != 0) {
        if (_mkdir(GNU_DIR) != 0 ) {
            printf("フォルダ(gnuplot)の作成に失敗しました\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // シミュレート結果保存用パスを作成
    if (stat(RESULT_DIR, &statBuf) != 0) {
        if (_mkdir("result") != 0) {
            printf("フォルダ(result)の作成に失敗しました\n");
            exit(EXIT_FAILURE);
        }
    }

    // 推力データ取得処理
    // テキストが開けるかどうかチェック
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
    for (i = 0; i < FThrustCount -1 ; i++) fprintf(fpXYData, "%d %lf %lf\n", i + 1, FThrustTimeStamp[i], FThrust[i]);
    fclose(fpXYData);

    // シミュレート開始
    printf("Simulate Start!\n");

    PrevProg = 0;

    // 風速のループ
    for (j = 1; j < (MAX_WIND + 1); j++) {
        FWindAbs = j;   // 風速
        // 風向のループ
        for (k = 0; k < DEV; k++) {
            FWindDeg = 2 * M_PI * k / DEV;    // 風向(角度)

            init(&FTE_06);

            double t = 0;
            FThrustIndex = 0;

            launcher(&FTE_06, &t);
            launch_clear = FTE_06.v_abs;
            FMaxHeight = FTE_06.z;
            thrust_inertia(&FTE_06, &FWind, &t);
            FMaxTime = t;

            // 落下タイプに応じて処理を変更
            switch (TYPE) {

                // タイプ0 : 2段パラシュート展開時処理
            case 0:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                break;

                // タイプ1 : 1段目のみ展開時処理
            case 1:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                break;

                // タイプ2 : 弾道落下時処理
            case 2:
                inertia_fall(&FTE_06, &FWind, &t);
                //領域(x方向)
                FMaxDistanceX = MAX_DISTANCE(FMaxDistanceXt, FMaxDistanceXf);
                //領域(y方向)
                FMaxDistanceY = MAX_DISTANCE(FMaxDistanceYt, FMaxDistanceYf);
                break;

                // タイプ3 : 開傘時処理
            case 3:
                inertia_fall(&FTE_06, &FWind, &t);
                FPara1Time = t;
                recovery(&FTE_06, &FWind, &t);
                //領域(x方向)
                FMaxDistanceX = MAX_DISTANCE(FMaxDistanceXt, FMaxDistanceXr);
                //領域(y方向)
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
                fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
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
            fflush(gp); // バッファに格納されているデータを吐き出す（必須）
            fprintf(gp, "exit\n");
            _pclose(gp);

            // 進行状況表示
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
            printf("\r|%s|  進行状況 : %.0lf ％", progress_bar, progress);
            fflush(stdout);
        }
    }
    printf("\nSimulate finished!\n");

    return 0;
}