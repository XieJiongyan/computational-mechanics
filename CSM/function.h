
#ifndef INCFILE_H_
#define INCFILE_H_


#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996) 

#include "stdio.h"

//�洢����նȾ���
extern FILE* fp1;
//�洢�ṹ�ڵ�λ��
extern FILE* fp2;
//�洢��Ԫ����
extern FILE* fp3;
//�洢�ڵ���
extern FILE* fp4;
//��ʾ��������
extern FILE* fp5;
//ԭʼ����
extern FILE* fp6;
//��Ԫ�Ľڵ�����
extern int ND;
//�����ڵ�����ɶ���
extern int NF;
//�ڵ�����
extern int NP;
//��Ԫ����
extern int NE;
//��Լ�������ɶ�����
extern int NR;
//��Ԫ������������Ԫ�������
extern int NM, NMN;
//N=NP*NF
extern int N;
//һά�洢AK��������
extern int NN;
//�����ڵ����ά����
extern double X[50], Y[50], Z[50];
extern double X2, X1, Y2, Y1, Z2, Z1, b;
//ÿ����Ԫ�Ľڵ�ţ�ME[1][i]��ŵ�i����Ԫ��һ���ڵ�����ţ�ME[2][i]��ŵ�i����Ԫ�ڶ����ڵ������
extern int ME[3][50];
//Լ����λ�ƺ�
extern int NRR[30];
extern int LD[50];
//ÿ����Ԫ�����
extern int NAE[100];
//AE[1][i]��ŵ�i�����͵ĵ�Ԫ������ģ����AE[2][i]��ŵ�i�����͵ĵ�Ԫ�ĺ�����
extern double AE[3][100];
//�ڵ��غ�
extern double P[100], P1[100];
//����ṹ�ڵ���
extern double PP[100];
//��Ԫ�����
extern double A;
//��Ԫ����ģ��
extern double E;
//��Ԫ�նȾ���
extern double TK[3][3];
//����ת������
extern double T[3][7];
//����ת�������ת��
extern double TT[7][3];
//����նȾ���
extern double AK[100];
//��������ϵ�µĵ�Ԫ�ն���
extern double AKEE[7][7];
//������˷�ʱ���м����
extern double s[7][3];
extern int IS[7];
//�˵�Ԫ�ĳ���
extern double L;
//��ԪӦ��
extern double SG;
//�ṹλ�ƾ���
extern double d[100];
//��Ԫλ�ƾ���
extern double ue[100];
//�ֲ������µĵ�Ԫλ�ƾ���
extern double dee[50][3];
//��Ԫ�������ֲ�����ϵ�µĵ�Ԫ�ڵ�����
extern double Fee[50][7];
//��Ԫ��������ϵ�еĽڵ���
extern double F[50][100];
//�ⷽ��ʱ�õ���L,Y����
extern double l[100][100], y[100];


//���������������Ҫ�Ĳ���
void scan();

#endif  /* INCFILE_H_ */

