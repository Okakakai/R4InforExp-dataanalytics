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

typedef struct
{
    int N;
    double **data;
} DataMatrix;

typedef struct
{
    int N;
    double *data;
} Solution;

//プロトタイプ宣言
double CalAverage(CoordinateData *);
double *CalDeviation(CoordinateData *);
double CalStandardDeviation(CoordinateData *);
double CalCovariance(CoordinateData *, CoordinateData *);
double CalCorrectionCoefficient(CoordinateData *, CoordinateData *);
void readData(char *, CoordinateData *, CoordinateData *);
double CalRegressionCoefficient(CoordinateData *, CoordinateData *);
double CalIntercept(CoordinateData *, CoordinateData *);
double *CalEstimatedValueData(CoordinateData *, CoordinateData *, DataMatrix *, Solution *);
double CalRegressionSumOfSquares(CoordinateData *, CoordinateData *, DataMatrix *, Solution *);
double CalTotalVarience(CoordinateData *);
double CalCoefficientOfDetermination(CoordinateData *, CoordinateData *, DataMatrix *, Solution *);
void GaussianElimination(DataMatrix *, Solution *);
double NormalEquation(CoordinateData *, CoordinateData *, DataMatrix *, Solution *);

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

double *CalEstimatedValueData(CoordinateData *specimenX, CoordinateData *specimenY, DataMatrix *dataMatrix, Solution *solution)
{
    NormalEquation(specimenX, specimenY, dataMatrix, solution);

    double *estimatedData = (double *)malloc(sizeof(double) * specimenX->N);
    double sum;
    for (int i = 0; i < specimenY->N; i++)
    {
        sum = 0;
        for (int j = 0; j < solution->N; j++)
        {
            sum += pow(specimenX->data[i], j) * solution->data[j];
        }
        sum += solution->data[solution->N];
        estimatedData[i] = sum;
    }
    return estimatedData;
}

double CalRegressionSumOfSquares(CoordinateData *specimenX, CoordinateData *specimenY, DataMatrix *dataMatrix, Solution *solution)
{
    EstimatedValueData *estimatedData;
    estimatedData = (EstimatedValueData *)malloc(sizeof(EstimatedValueData));

    double regressionSum = 0;
    double aveY = CalAverage(specimenY);

    estimatedData->N = specimenX->N;
    estimatedData->data = CalEstimatedValueData(specimenX, specimenY, dataMatrix, solution);
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

double CalCoefficientOfDetermination(CoordinateData *specimenX, CoordinateData *specimenY, DataMatrix *dataMatrix, Solution *solution)
{
    double regressionSum = 0, totalVarience = 0;
    regressionSum = CalRegressionSumOfSquares(specimenX, specimenY, dataMatrix, solution);
    totalVarience = CalTotalVarience(specimenY);
    double coefficientOfDetermination = regressionSum / totalVarience;
    return coefficientOfDetermination;
}

void GaussianElimination(DataMatrix *dataMatrix, Solution *solution)
{
    double x[solution->N];
    double tmpColumn[1][dataMatrix->N + 1];
    int Max, pivot;
    double diagonalElement;
    for (int i = 0; i < dataMatrix->N; i++)
    {
        Max = 0;
        pivot = i;

        for (int l = i; l < dataMatrix->N; l++)
        {
            if (fabs(dataMatrix->data[l][i]) > Max)
            { // i列の中で一番値が大きい行を選ぶ
                Max = fabs(dataMatrix->data[l][i]);
                pivot = l;
            }
        }

        if (pivot != i)
        { // pivotがiと違えば、行の入れ替え
            for (int j = 0; j < dataMatrix->N + 1; j++)
            {
                tmpColumn[0][j] = dataMatrix->data[i][j];
                dataMatrix->data[i][j] = dataMatrix->data[pivot][j];
                dataMatrix->data[pivot][j] = tmpColumn[0][j];
            }
        }
    }

    for (int k = 0; k < dataMatrix->N; k++)
    {
        diagonalElement = dataMatrix->data[k][k]; //対角要素を保存
        dataMatrix->data[k][k] = 1;               //対角要素は１になることがわかっているから

        for (int j = k + 1; j < dataMatrix->N + 1; j++)
        {
            dataMatrix->data[k][j] /= diagonalElement;
        }

        for (int i = k + 1; i < dataMatrix->N; i++)
        {
            diagonalElement = dataMatrix->data[i][k];

            for (int j = k + 1; j < dataMatrix->N + 1; j++)
            {
                dataMatrix->data[i][j] -= diagonalElement * dataMatrix->data[k][j];
            }
            dataMatrix->data[i][k] = 0; //０となることがわかっているところ
        }
    }

    //解の計算
    for (int i = dataMatrix->N - 1; i >= 0; i--)
    {
        x[i] = dataMatrix->data[i][dataMatrix->N];
        for (int j = dataMatrix->N - 1; j > i; j--)
        {
            x[i] -= dataMatrix->data[i][j] * x[j];
        }
    }

    printf("解は\n");
    for (int i = 0; i < dataMatrix->N; i++)
    {
        printf("%f\n", x[i]);
        solution->data[i] = x[i];
    }
}

double NormalEquation(CoordinateData *specimenX, CoordinateData *specimenY, DataMatrix *dataMatrix, Solution *solution)
{
    dataMatrix->data = (double **)malloc(sizeof(double *) * dataMatrix->N);
    for (int i = 0; i < dataMatrix->N; i++)
    {
        dataMatrix->data[i] = (double *)malloc(sizeof(double) * dataMatrix->N + 1);
    }
    solution->data = (double *)malloc(sizeof(double) * solution->N);

    /*初期化*/
    for (int i = 0; i < dataMatrix->N; i++)
    {
        for (int j = 0; j < dataMatrix->N + 1; j++)
        {
            dataMatrix->data[i][j] = 0.0;
        }
    }

    /*ガウスの消去法で解く行列の作成*/
    for (int i = 0; i < dataMatrix->N; i++)
    {
        for (int j = 0; j < dataMatrix->N; j++)
        {
            for (int k = 0; k < specimenX->N; k++)
            {
                // シグマの実行
                dataMatrix->data[i][j] += pow(specimenX->data[k], i + j);
            }
        }
    }
    for (int i = 0; i < dataMatrix->N; i++)
    {
        for (int k = 0; k < specimenX->N; k++)
        {
            // 一番右のやつ
            dataMatrix->data[i][dataMatrix->N] += pow(specimenX->data[k], i) * specimenY->data[k];
        }
    }
    GaussianElimination(dataMatrix, solution);
}

int main(int argc, char **argv)
{
    CoordinateData *dataX, *dataY;
    dataX = (CoordinateData *)malloc(sizeof(CoordinateData));
    dataY = (CoordinateData *)malloc(sizeof(CoordinateData));
    readData(argv[1], dataX, dataY);

    DataMatrix *dataMatrix;
    dataMatrix = (DataMatrix *)malloc(sizeof(DataMatrix));

    Solution *solution;
    solution = (Solution *)malloc(sizeof(Solution));

    int degree = 3;
    // printf("n次を計算するプログラム\n\n nを入力してください : ");
    // scanf("%d", &degree);
    dataMatrix->N = degree + 1;
    solution->N = degree + 1;
    double coefficientOfDetermination = CalCoefficientOfDetermination(dataX, dataY, dataMatrix, solution);
    printf("CoefficientOfDetermination : %lf\n", coefficientOfDetermination);
}