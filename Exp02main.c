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

typedef struct
{
    int N;
    double *data;
} EstimatedValueData;

double CalAverage(CoordinateData *);
double *CalDeviation(CoordinateData *);
double CalStandardDeviation(CoordinateData *);
double CalCovariance(CoordinateData *, CoordinateData *);
double CalCorrectionCoefficient(CoordinateData *, CoordinateData *);
void readData(char *, CoordinateData *, CoordinateData *);
double CalRegressionCoefficient(CoordinateData *, CoordinateData *);
double CalIntercept(CoordinateData *, CoordinateData *);
double *CalEstimatedValueData(CoordinateData *, CoordinateData *);
double CalRegressionSumOfSquares(CoordinateData *, CoordinateData *);
double CalTotalVarience(CoordinateData *);
double CalCoefficientOfDetermination(CoordinateData *, CoordinateData *);

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

double CalRegressionCoefficient(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double regressionCoefficient = 0;
    double correctionCoefficient = CalCorrectionCoefficient(specimenX, specimenY);
    double standardDeviationX = CalStandardDeviation(specimenX);
    double standardDeviationY = CalStandardDeviation(specimenY);
    regressionCoefficient = correctionCoefficient * (standardDeviationX / standardDeviationY);
    return regressionCoefficient;
}

double CalIntercept(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double intercept = 0;
    Average *ave;
    ave = (Average *)malloc(sizeof(Average));
    ave->aveX = CalAverage(specimenX);
    ave->aveY = CalAverage(specimenY);

    double regressionCoefficient = CalRegressionCoefficient(specimenX, specimenY);
    intercept = ave->aveY - ave->aveX * regressionCoefficient;
    return intercept;
}

double *CalEstimatedValueData(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double intercept = CalIntercept(specimenX, specimenY);
    double regressionCoefficient = CalRegressionCoefficient(specimenX, specimenY);

    double *estimatedData = (double *)malloc(sizeof(double) * specimenX->N);

    for (int i = 0; i < specimenY->N; i++)
    {
        estimatedData[i] = regressionCoefficient * specimenX->data[i] + intercept;
    }
    return estimatedData;
}

double CalRegressionSumOfSquares(CoordinateData *specimenX, CoordinateData *specimenY)
{
    EstimatedValueData *estimatedData;
    estimatedData = (EstimatedValueData *)malloc(sizeof(EstimatedValueData));

    double regressionSum = 0;
    double aveY = CalAverage(specimenY);

    estimatedData->N = specimenX->N;
    estimatedData->data = CalEstimatedValueData(specimenX, specimenY);
    for (int i = 0; i < estimatedData->N; i++)
    {
        regressionSum += (estimatedData->data[i] - aveY) * (estimatedData->data[i] - aveY);
    }

    return regressionSum;
}

double CalTotalVarience(CoordinateData *specimenY)
{
    double *deviationY = (double *)malloc(sizeof(double) * specimenY->N);
    deviationY = CalDeviation(specimenY);

    double totalVarience = 0;
    for (int i = 0; i < specimenY->N; i++)
    {
        totalVarience += deviationY[i] * deviationY[i];
    }

    return totalVarience;
}

double CalCoefficientOfDetermination(CoordinateData *specimenX, CoordinateData *specimenY)
{
    double regressionSum = 0, totalVarience = 0;
    regressionSum = CalRegressionSumOfSquares(specimenX, specimenY);
    totalVarience = CalTotalVarience(specimenY);
    double coefficientOfDetermination = regressionSum / totalVarience;
    return coefficientOfDetermination;
}

int main(int argc, char **argv)
{
    CoordinateData *dataX, *dataY;
    dataX = (CoordinateData *)malloc(sizeof(CoordinateData));
    dataY = (CoordinateData *)malloc(sizeof(CoordinateData));
    readData(argv[1], dataX, dataY);
    double regressionCoefficient = CalRegressionCoefficient(dataX, dataY);
    double intercept = CalIntercept(dataX, dataY);
    printf("regreesionCoefficient intercept  *%lf+%lf\n", regressionCoefficient, intercept);
    double coefficientOfDetermination = CalCoefficientOfDetermination(dataX, dataY);
    printf("coefficientOfDetermination : %lf\n", coefficientOfDetermination);
}