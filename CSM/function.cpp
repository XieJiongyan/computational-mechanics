#include "function.h"


//�洢����նȾ���
FILE* fp1;
//�洢�ṹ�ڵ�λ��
FILE* fp2;
//�洢��Ԫ����
FILE* fp3;
//�洢�ڵ���
FILE* fp4;
//��ʾ��������
FILE* fp5;
//ԭʼ����
FILE* fp6;
//��Ԫ�Ľڵ�����
int ND;
//�����ڵ�����ɶ���
int NF;
//�ڵ�����
int NP;
//��Ԫ����
int NE;
//��Լ�������ɶ�����
int NR;
//��Ԫ������������Ԫ�������
int NM, NMN;
//N=NP*NF
int N;
//һά�洢AK��������
int NN;
//�����ڵ����ά����
double X[50], Y[50], Z[50];
double X2, X1, Y2, Y1, Z2, Z1, b;
//ÿ����Ԫ�Ľڵ�ţ�ME[1][i]��ŵ�i����Ԫ��һ���ڵ�����ţ�ME[2][i]��ŵ�i����Ԫ�ڶ����ڵ������
int ME[3][50];
//Լ����λ�ƺ�
int NRR[30];
int LD[50];
//ÿ����Ԫ�����
int NAE[100];
//AE[1][i]��ŵ�i�����͵ĵ�Ԫ������ģ����AE[2][i]��ŵ�i�����͵ĵ�Ԫ�ĺ�����
double AE[3][100];
//�ڵ��غ�
double P[100], P1[100];
//����ṹ�ڵ���
double PP[100];
//��Ԫ�����
double A;
//��Ԫ����ģ��
double E;
//��Ԫ�նȾ���
double TK[3][3];
//����ת������
double T[3][7];
//����ת�������ת��
double TT[7][3];
//����նȾ���
double AK[100];
//��������ϵ�µĵ�Ԫ�ն���
double AKEE[7][7];
//������˷�ʱ���м����
double s[7][3];
int IS[7];
//�˵�Ԫ�ĳ���
double L;
//��ԪӦ��
double SG;
//�ṹλ�ƾ���
double d[100];
//��Ԫλ�ƾ���
double ue[100];
//�ֲ������µĵ�Ԫλ�ƾ���
double dee[50][3];
//��Ԫ�������ֲ�����ϵ�µĵ�Ԫ�ڵ�����
double Fee[50][7];
//��Ԫ��������ϵ�еĽڵ���
double F[50][100];
//�ⷽ��ʱ�õ���L,Y����
double l[100][100], y[100];


void scan()
{
	int i;
	fp5 = fopen("��ʾ��������1.txt", "w");
	fp6 = fopen("ԭʼ����1.txt", "r");

	fscanf(fp6, "%d", &ND);
	fprintf(fp5, "ND(��Ԫ�ڵ���)=%d\n", ND);

	fscanf(fp6, "%d", &NF);
	fprintf(fp5, "NF(��Ԫ���ɶ���)=%d\n", NF);

	fscanf(fp6, "%d", &NP);
	fprintf(fp5, "NP(�ڵ�����)=%d\n", NP);

	fscanf(fp6, "%d", &NE);
	fprintf(fp5, "NE(��Ԫ�������)=%d\n", NE);

	fscanf(fp6, "%d", &NM);
	fprintf(fp5, "NM(��Լ�������ɶ�����)=%d\n", NM);

	fscanf(fp6, "%d", &NR);
	fprintf(fp5, "NR(��Լ�������ɶ�����)=%d\n", NR);

	//for (i = 1; i <= NE; ++i) {
	//	fscanf(fp6, "%d", &NAE[i]);

	//}
}