
#ifndef INCFILE_H_
#define INCFILE_H_


#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996) 

#include "stdio.h"

//存储整体刚度矩阵
extern FILE* fp1;
//存储结构节点位移
extern FILE* fp2;
//存储单元内力
extern FILE* fp3;
//存储节点力
extern FILE* fp4;
//显示输入数据
extern FILE* fp5;
//原始数据
extern FILE* fp6;
//单元的节点总数
extern int ND;
//单个节点的自由度数
extern int NF;
//节点总数
extern int NP;
//单元总数
extern int NE;
//受约束的自由度总数
extern int NR;
//单元类型总数，单元的类别数
extern int NM, NMN;
//N=NP*NF
extern int N;
//一维存储AK的总容量
extern int NN;
//各个节点的三维坐标
extern double X[50], Y[50], Z[50];
extern double X2, X1, Y2, Y1, Z2, Z1, b;
//每个单元的节点号，ME[1][i]存放第i个单元第一个节点坐标号，ME[2][i]存放第i个单元第二个节点坐标号
extern int ME[3][50];
//约束的位移号
extern int NRR[30];
extern int LD[50];
//每个单元的类别
extern int NAE[100];
//AE[1][i]存放第i种类型的单元的杨氏模量，AE[2][i]存放第i种类型的单元的横截面积
extern double AE[3][100];
//节点载荷
extern double P[100], P1[100];
//整体结构节点力
extern double PP[100];
//单元横截面
extern double A;
//单元杨氏模量
extern double E;
//单元刚度矩阵
extern double TK[3][3];
//坐标转换矩阵
extern double T[3][7];
//坐标转换矩阵的转置
extern double TT[7][3];
//整体刚度矩阵
extern double AK[100];
//整体坐标系下的单元刚度阵
extern double AKEE[7][7];
//作矩阵乘法时的中间矩阵
extern double s[7][3];
extern int IS[7];
//杆单元的长度
extern double L;
//单元应力
extern double SG;
//结构位移矩阵
extern double d[100];
//单元位移矩阵
extern double ue[100];
//局部坐标下的单元位移矩阵
extern double dee[50][3];
//单元内力（局部坐标系下的单元节点力）
extern double Fee[50][7];
//单元整体坐标系中的节点力
extern double F[50][100];
//解方程时用到的L,Y矩阵
extern double l[100][100], y[100];


//输入输出仿真所需要的参数
void scan();

#endif  /* INCFILE_H_ */

