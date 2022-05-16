#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    int N;
    double *data;
} CoordinateData;

typedef struct
{
    double aveX;
    double aveY;
} Average;

typedef struct
{
    double *deviX;
    double *deviY;
} Deviation;

typedef struct
{
    double stdDeviX;
    double stdDeviY;
} StdDevi;

//プロトタイプ宣言
double CalAverage(CoordinateData *);
double *CalDeviation(CoordinateData *);
double CalStandardDeviation(CoordinateData *);
double CalCovariance(CoordinateData *, CoordinateData *);
double CalCorrectionCoefficient(CoordinateData *, CoordinateData *);
void readData(char *, CoordinateData *, CoordinateData *);

double CalAverage(CoordinateData *specimen)
{
    double sum = 0;
    for (int i = 0; i < specimen->N; i++)
    {
        sum += specimen->data[i];
    }
    double ave = sum / specimen->N;
    return ave;
}

double *CalDeviation(CoordinateData *specimen)
{
    double ave = CalAverage(specimen);
    double *devi = (double *)malloc(sizeof(double) * specimen->N);
    for (int i = 0; i < specimen->N; i++)
    {
        devi[i] = specimen->data[i] - ave;
    }
    return devi;
}

double CalStandardDeviation(CoordinateData *specimen)
{
    double *deviation = CalDeviation(specimen);
    double sum = 0;
    for (int i = 0; i < specimen->N; i++)
    {
        sum += deviation[i] * deviation[i];
    }
    double variance = sum / specimen->N;
    double stdDevi = sqrt(variance);
    return stdDevi;
}

double CalCovariance(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double *deviationX = (double *)malloc(sizeof(double) * specimenX->N);
    double *deviationY = (double *)malloc(sizeof(double) * specimenY->N);
    deviationX = CalDeviation(specimenX);
    deviationY = CalDeviation(specimenY);
    // Nが同数であるのは間違いない
    int wholeNum = specimenX->N;
    double sum = 0;
    for (int i = 0; i < wholeNum; i++)
    {
        sum += deviationX[i] * deviationY[i];
    }

    double covarience = sum / wholeNum;

    return covarience;
}

double CalCorrectionCoefficient(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double covarience = CalCovariance(specimenX, specimenY);
    double stdDeviX = CalStandardDeviation(specimenX);
    double stdDeviY = CalStandardDeviation(specimenY);
    double CorrectionCoefficient = covarience / (stdDeviX * stdDeviY);

    return CorrectionCoefficient;
}

void readData(char *filename, CoordinateData *dataX, CoordinateData *dataY)
{
    FILE *fp;
    int i;
    if ((fp = fopen(filename, "r")) == NULL)
    {
        printf("ファイル読み込みerror");
        exit(1);
    }
    int N;
    fscanf(fp, "%d", &N);
    dataX->N = N;
    dataY->N = N;
    if (dataX->data == NULL)
    {
        printf("NULL");
    }
    dataX->data = (double *)malloc(sizeof(double) * N);
    if (dataX->data == NULL)
    {
        printf("NULL");
    }
    dataY->data = (double *)malloc(sizeof(double) * N);
    for (i = 0; i < N; i++)
    {
        fscanf(fp, "%lf %lf\n", &dataX->data[i], &dataY->data[i]);
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    CoordinateData *dataX, *dataY;
    dataX = (CoordinateData *)malloc(sizeof(CoordinateData));
    dataY = (CoordinateData *)malloc(sizeof(CoordinateData));

    readData(argv[1], dataX, dataY);

    Average *ave;
    ave = (Average *)malloc(sizeof(Average));
    ave->aveX = CalAverage(dataX);
    printf("aveX : %lf\n", ave->aveX);
    ave->aveY = CalAverage(dataY);
    printf("aveY : %lf\n", ave->aveY);

    //標準偏差
    StdDevi *stdDevi;
    stdDevi = (StdDevi *)malloc(sizeof(StdDevi));
    stdDevi->stdDeviX = CalStandardDeviation(dataX);
    stdDevi->stdDeviY = CalStandardDeviation(dataY);
    printf("stdDeviX : %lf\nstdDeviY : %lf\n", stdDevi->stdDeviX, stdDevi->stdDeviY);
    double covarience = CalCovariance(dataX, dataY);
    double correctionCoefficient = CalCorrectionCoefficient(dataX, dataY);
    printf("correctionCoefficient : %lf\n", correctionCoefficient);
}
