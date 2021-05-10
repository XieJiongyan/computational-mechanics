#include <math.h>
#include <stdio.h>
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
    if(fp5 = NULL) {
        printf("ERR");
    }
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

	for (i = 1; i <= NE; i++) {
		fscanf(fp6, "%d", &NAE[i]);}
        for(i=1;i<=NE;i++)
        {fprintf(fp5,"第%d个单元的类别=%d/n",i,NAE[i]);}
    
    for (i = 1; i <= NM; i++) {
		fscanf(fp6, "%1e%1e", &AE[1][i],&AE[2][i]);
        }
    for(i=1;i<=NM;i++)
    {  fprintf(fp5,"第%d个单元类别的E和A值：%e %e\n",i,AE[1][i],AE[2][i]); }
    
    for (i = 1; i <= NP; i++) {
		fscanf(fp6, "%lf%lf%lf", &X[i],&Y[i],&Z[i]);}
    for(i=1;i<=NP;i++)
    {  fprintf(fp5,"第%d个节点的坐标值：%f %f%f\n",i,X[i],Y[i],Z[i]); }
    
    for (i = 1; i <= NE; i++) {
		fscanf(fp6, "%d%d", &ME[1][i],&ME[2][i]);}
    for(i=1;i<=NE;i++)
    {  fprintf(fp5,"第%d个单元的节点号：%d %d\n",i,ME[1][i],ME[2][i]); }
    
    for (i = 1; i <= NR; i++) {
		fscanf(fp6, "%d", &NRR[i]);}
    for(i=1;i<=NR;i++)
    {  fprintf(fp5,"第%d个单元的位移号：%d\n",i,NRR[i]); }

    for (i = 1; i <= NP*NF; i++) {
		fscanf(fp6, "%lf", &P[i]);}
    for(i=1;i<=NP*NF;i++)
    {   P1[i]=P[i];
        fprintf(fp5,"第%d个位移方向上的外加载荷：%f %f\n",i,P[i],P1[i]); }
    
}
//求每个杆的长度（i）表示单元号
void Length(int i)
{
    X2=X[ME[2][i]]; X1=X[ME[1][i]];
    Y2=Y[ME[2][i]]; Y1=Y[ME[1][i]];
    Z2=Y[ME[2][i]]; Z1=Z[ME[1][i]];
    L=sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1));
}
//杆单元的单元刚度阵（单元坐标系下）（i表示单元号）
void StiffnessMatrix_unit(int i)
{
    Length(i);
    NMN=NAE[i];
    E=AE[1][NMN];
    A=AE[2][NMN];
    TK[1][1]=E*A/L; TK[1][2]=-E*A/L;
    TK[2][1]=-E*A/L; TK[2][2]=E*A/L;
}
//坐标转换矩阵（i表示单元号）
void TransformMatrix(int i)
{
    int m, n;
    Length(i);
    T[1][1]=(X2-X1)/L; T[2][4]=(X2-X1)/L;
    T[1][2]=(Y2-Y1)/L; T[2][5]=(Y2-Y1)/L;
    T[1][3]=(Z2-Z1)/L; T[2][6]=(X2-X1)/L;//其他元素已在声明时赋值
    for(m=1;m<=2;m++)
    for(n=1;n<=(NF*ND);n++)
    TT[n][m]=T[m][n];//求出坐标转换矩阵[T]的转置

}
//总体坐标下的单元刚度阵（矩阵乘法）（i表示单元号）
void MultiplyMatrix(int i)
{
    int j,m,n;
    double b;
    StiffnessMatrix_unit(i);
    for(n=1;n<=(NF*ND);n++)
    {
        for(m=1;m<=2;m++)
        {
            b=0.0;
        for (j=1;j<=2;j++)
        b+=TT[n][j]*TK[j][m];
        s[n][m]=b;}}
    for(m=1;m<=(NF*ND);m++)
    {
        for(n=1;n<=(NF*ND);N++){
            b=0.0;
            for(j=1;j<=2;j++)
            b+=s[m][j]*T[j][n];
            AKEE[m][n]=b;
        }
    }
            
}
//形成LD数组
void FLd()
{
    int k,i,j,L,IG,J,NN,N;
    LD[1]=1;
    //按节点循环，确定与其有关的的最小节点号和该点所占行的LD数组
    for(k=1;k<=NP;k++)
    {IG=100000;
    //按单元循环，确定其最小节点号
    for(i=1;i<=NE;i++)
    for(j=1;j<=ND;j++)
    //判别单元是否含有K点
    {if(ME[j][i]!=k)continue;
    //寻找与k点有关的最小节点号放入IG
    for(L=1;L<=ND;L++)
    {if(ME[L][i]>=IG)continue;
    IG=ME[L][i];
    }}
    //确定K点所含的LD数组
    for(i=1;i<=NF;i++)
    {
    J=NF*(k-1)+i;//K点所对应的刚阵的行号
    if(J==1)continue;
    LD[J]=LD[J-1]+NF*(k-IG)+i;
    }
    }
    //确定一维数组的总容量
    N=NP*NF;
    NN=LD[N];

}

//形成IS数组
void FIS(int i)
{
    //对于杆单元，单元位移为I（I=1，2，3.....),他在结构位移中的编号为IS[I]
    IS[1]=(ME[1][i]-1)*NF+1;
    IS[2]=(ME[1][i]-1)*NF+2;
    IS[3]=(ME[1][i]-1)*NF+3;
    IS[4]=(ME[2][i]-1)*NF+1;
    IS[5]=(ME[2][i]-1)*NF+2;
    IS[6]=(ME[2][i]-1)*NF+3;
}
//组集结构刚度阵
void StructureMatrix()
{
    int i, j,m,ISS,NI,NJ,IJ;
    FLd();
    for(m=1;m<=NE;m++)
    {
        MultiplyMatrix(m);
        FIS(m);
        for(i=1;i<=(NF*ND);i++)
        for(j=1;j<=(NF*ND);j++)
        {
            ISS=IS[i]-IS[j];
            if(ISS>=0)
            {
                NI=IS[i];IJ=LD[NI]-(NI-IS[j]);
                AK[IJ]+=AKEE[i][j];
            }
        }
    }
    //约束处理(置大数法)
    for(i=1;i<=NR;i++)
    {
        NI=NRR[i];NJ=LD[NI];
        AK[NJ]=1e25;
    }
}

void cholesky(int n,double a[100],double x[100])
{
    int i,j,k,ij,kj,ii,jj,ik,jk,kk,iig,ig,igp,jgp,mi,mj,mij;
    for(i=1;i<=n;i++)
    {
        if(i!=1)
        {
            mi=i-(LD[i]-LD[i-1])+1;
            if(mi!=i)
            {
                iig=LD[i]-i;
                for(j=mi;j<=i-1;j++)
                {
                    if(j!=mi)
                    {
                        mj=j-(LD[j]-LD[j-1])+1;
                        igp=LD[i]-(i-j);
                        if(mj<mi) mij=mi;
                        else mij=mj;
                        jgp=LD[j]-j;
                        if(mij<=j-1)
                        {
                            for(k=mij;k<=j-1;k++)
                            {
                                ik=iig+k;jk=jgp+k;kk=LD[k];
                                a[igp]-=a[ik]*a[kk]*a[jk];
                            }
                        }
                    }

                    if(j==mi)igp=LD[i-1]+1;
                    ii=LD[j];
                    a[igp]=a[igp]/a[ii];
                    x[i]-=a[igp]*a[ii]*x[j];
                }
                ij=LD[i];
                for(k=mi;k<=i-1;k++)
                { 
                    ii=iig+k;jj=LD[k];
                    a[ij]-=a[ii]*a[ii]*a[jj];
                }
            }
        }
        ij=LD[i];
        x[i]=x[i]/a[ij];


    }
    for(i=n;i>=2;i--)
    {
        mi=i-(LD[i]-LD[i-1])+1;
        if(mi==i)continue;
        iig=LD[i]-i;
        for(k=mi;k<=i-1;k++)
        {
            ij=iig+k;
            x[k]-=a[ij]*x[i];
        }
    }
    fp2=fopen("整体节点位移.txt","w");
    for(i=1;i<=n;i++)    {
        fprintf(fp2,"%d %1e%1e\n",i,x[i],a[i]);
    }
    fclose(fp2);
}

void InternalForce() {
    int n, m, i, k, j, ii, mm, jj;
    fp3 = fopen("单元内力.txt", "w");
    fp4 = fopen("约束反力.txt", "w");
    for (i = 1; i <= NE; i++) {
        TransformMatrix(i);
        StiffnessMatrix_unit(i);
        FIS(i);
        for (k = 1; k <= (NF * ND); k++) {
            ue[k] = P[IS[k]];
        }
        fprintf(fp3, "第%d个单元整体坐标系下的节点位移\n", i);
        for (ii = 1; ii <= ND * NF; ii++) {
            fprintf(fp3, " %d %le\n", ii, ue[ii]);
        }
        for (j = 1; j <= 2; j++) {
            for (k = 1; k <= 6; k++) {
                dee[i][j] += T[j][k] * ue[k];
            }
        }
        fprintf(fp3, "第%d个单元局部坐标系下的节点位移\n", i);
        fprintf(fp3, " %1e %1e \n", dee[i][1], dee[i][2]);
        for (m = 1; m<= ND; m++) {
            for (j = 1; j <= ND; j++) {
                Fee[i][m] += TK[m][j] * dee[i][j];
            }
        }
        fprintf(fp4, "第%d个单元整体坐标系下的节点力\n", i);
        for (ii = 1; ii <= ND * NF; ii++) {
            fprintf(fp4, " %d %1e \n", ii, F[i][ii]);
        }

        for (jj = 1; jj <= 6; jj++) {
            PP[IS[jj]] = PP[IS[jj]] + F[i][jj];
        }
    }
    for (k = 1; k <= NF * NP; k++) 
        P1[k] = PP[k] - P1[k];
    for (m = 1; m <= NP; m++) {
        for (n = 1; n <= NF; n++) 
            fprintf(fp4, "第%d个节点的第%d个约束力%f\n", m, n, P1[(m - 1) * NF + n]);
        fprintf(fp3, "               \n");
    }
    fclose(fp3);
    fclose(fp4);
        
}

int main()
{
	scan();
    StructureMatrix();//计算刚度矩阵
    N=NP*NF;
    cholesky(N,AK,P);//计算节点位移
    InternalForce();//计算单元内力
}
