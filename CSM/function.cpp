#include "function.h"


//存储整体刚度矩阵
FILE* fp1;
//存储结构节点位移
FILE* fp2;
//存储单元内力
FILE* fp3;
//存储节点力
FILE* fp4;
//显示输入数据
FILE* fp5;
//原始数据
FILE* fp6;
//单元的节点总数
int ND;
//单个节点的自由度数
int NF;
//节点总数
int NP;
//单元总数
int NE;
//受约束的自由度总数
int NR;
//单元类型总数，单元的类别数
int NM, NMN;
//N=NP*NF
int N;
//一维存储AK的总容量
int NN;
//各个节点的三维坐标
double X[50], Y[50], Z[50];
double X2, X1, Y2, Y1, Z2, Z1, b;
//每个单元的节点号，ME[1][i]存放第i个单元第一个节点坐标号，ME[2][i]存放第i个单元第二个节点坐标号
int ME[3][50];
//约束的位移号
int NRR[30];
int LD[50];
//每个单元的类别
int NAE[100];
//AE[1][i]存放第i种类型的单元的杨氏模量，AE[2][i]存放第i种类型的单元的横截面积
double AE[3][100];
//节点载荷
double P[100], P1[100];
//整体结构节点力
double PP[100];
//单元横截面
double A;
//单元杨氏模量
double E;
//单元刚度矩阵
double TK[3][3];
//坐标转换矩阵
double T[3][7];
//坐标转换矩阵的转置
double TT[7][3];
//整体刚度矩阵
double AK[100];
//整体坐标系下的单元刚度阵
double AKEE[7][7];
//作矩阵乘法时的中间矩阵
double s[7][3];
int IS[7];
//杆单元的长度
double L;
//单元应力
double SG;
//结构位移矩阵
double d[100];
//单元位移矩阵
double ue[100];
//局部坐标下的单元位移矩阵
double dee[50][3];
//单元内力（局部坐标系下的单元节点力）
double Fee[50][7];
//单元整体坐标系中的节点力
double F[50][100];
//解方程时用到的L,Y矩阵
double l[100][100], y[100];


void scan()
{
	int i;
	fp5 = fopen("显示输入数据1.txt", "w");
	fp6 = fopen("原始数据1.txt", "r");

	fscanf(fp6, "%d", &ND);
	fprintf(fp5, "ND(单元节点数)=%d\n", ND);

	fscanf(fp6, "%d", &NF);
	fprintf(fp5, "NF(单元自由度数)=%d\n", NF);

	fscanf(fp6, "%d", &NP);
	fprintf(fp5, "NP(节点总数)=%d\n", NP);

	fscanf(fp6, "%d", &NE);
	fprintf(fp5, "NE(单元类别总数)=%d\n", NE);

	fscanf(fp6, "%d", &NM);
	fprintf(fp5, "NM(受约束的自由度总数)=%d\n", NM);

	fscanf(fp6, "%d", &NR);
	fprintf(fp5, "NR(受约束的自由度总数)=%d\n", NR);

	//for (i = 1; i <= NE; ++i) {
	//	fscanf(fp6, "%d", &NAE[i]);

	//}
}