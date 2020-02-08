//      program Adaptive Natural Density Partitioning AdNDP
//!
//!     version A.6 May 07, 2009
//!

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define PNAt 128
#define PNVl 9
#define PNTot 15
#define PNMax 4096
#define PBSz 1024
/************************************************************************
! NAt - number of atoms
! NVal - number of valence electronic pairs
! NTot - number of core and valence electronic pairs
! NMax - the maximum number of NBO vectors to be found;
!        defines the size of NBOVec, NBOOcc, NMax
! BSz - total number of basis functions;
!       defines the size of the density matrix
!************************************************************************/

char CMOFile[256], AdNDPFile[256];
int NAt, NVal, NTot, BSz;
int AtBs[PNAt], AtBsRng[PNAt][2];
double DMNAO[PBSz][PBSz];
double Thr[PNAt];

void AdNBO(int *MnSrch, double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int *NBOAmnt, double *DResid);
void SortPrel(double PrelOcc[PNMax], double PrelVec[PBSz][PNMax], int PrelCtr[PNAt][PNMax],
			  int PP, double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int NCtr, int *IndF);
void DepleteDMNAO(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int IndS, int IndF);
void TraceDMNAO(int BSz, double *DResid);
void BlockDMNAO(int CBl[PNAt], int NCtr, double DUMMY[PBSz][PBSz]);
void EigenSystem(double a[PBSz][PBSz], int n, double d[PBSz], double v[PBSz][PBSz]);
void EigenSrt(double d[PBSz], double v[PBSz][PBSz], int n);
void Subsets(int AtBl[100000][PNAt], int *NLn, int NClmn, int NAt);
void Input(double NAOAO[PBSz][PBSz]);
void Output(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int NBOAmnt, double DResid, double NAOAO[PBSz][PBSz]);
void BasisChange(double Old[PBSz][PNMax], double New[PBSz][PNMax], double Trans[PBSz][PBSz], int Row, int Col);
void NBOPlotMolden(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOAmnt);

void main()
{
	int MnSrch;

	double NBOOcc[PNMax], NBOVec[PBSz][PNMax];
	double NBOVecAO[PBSz][PNMax];
	double DResid;
	int NBOCtr[PNAt][PNMax];
	int NBOAmnt;
	double NAOAO[PBSz][PBSz];

	int i, ii, j, jj, k, l, m, n;

	MnSrch = 0;

	Input(NAOAO);
	AdNBO(&MnSrch, NBOOcc, NBOVec, NBOCtr, &NBOAmnt, &DResid);
	Output(NBOOcc, NBOVec, NBOCtr, NBOAmnt, DResid, NAOAO);
	BasisChange(NBOVec, NBOVecAO, NAOAO, BSz, NBOAmnt);
	NBOPlotMolden(NBOOcc, NBOVecAO, NBOAmnt);
}

//************************************************************************
void AdNBO(int *MnSrch, double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int *NBOAmnt, double *DResid)
{
	//prel初步的
	double PrelOcc[PNMax], PrelVec[PBSz][PNMax];
	int PrelCtr[PNAt][PNMax];

	int i, j, k, l, m, imax;
	int Cnt, PP;
	int IndS, IndF;
	//NCtr counter
	int NCtr, AtBlQnt;
	int AtBl[100000][PNAt], CBl[PNAt];
	double threshold, vmax;
	double DUMMY[PBSz][PBSz];
	double EiVal[PBSz], EiVec[PBSz][PBSz];
	FILE *fp;

	fp = fopen("out_debug", "w");

	IndS = 0;
	IndF = 0;
	*NBOAmnt = 0;

	for (i = 0; i < PBSz; i++)
		for (j = 0; j < PBSz; j++)
			DUMMY[i][j] = 0.0;

	for (NCtr = 1; NCtr <= NAt; NCtr++)
	{ //NAt=20
		threshold = Thr[NCtr - 1];
		if (threshold == 0)
		{
			printf("skippint NCtr=%d\n", NCtr);
			continue;
		}
		if (*MnSrch == 0)
			Subsets(AtBl, &AtBlQnt, NCtr, NAt); //对Atbl赋值
		printf("NCtr = %3d, AtBlQnt= %6d\n", NCtr, AtBlQnt);
		PP = 0;

		Cnt = 1;
		while (1)
		{
			for (k = 0; k < AtBlQnt; k++)
			{
				for (j = 0; j < NCtr; j++)
					CBl[j] = AtBl[k][j];
				BlockDMNAO(CBl, NCtr, DUMMY);
				EigenSystem(DUMMY, BSz, EiVal, EiVec);
				EigenSrt(EiVal, EiVec, BSz);

				for (i = 0, vmax = 0.0; i < BSz; i++)
				{
					if (vmax < EiVal[i])
					{
						vmax = EiVal[i];
						imax = i;
					}
				}

				i = imax;
				if (fabs(vmax - 2.0) <= threshold)
				{ //occupied
					PrelOcc[PP] = EiVal[i];
					for (j = 0; j < BSz; j++)
						PrelVec[j][PP] = EiVec[j][i];
					for (j = 0; j < NCtr; j++)
						PrelCtr[j][PP] = CBl[j];
					PP += 1;
					fprintf(fp, "AtBl ");
					for (j = 0; j < NCtr; j++)
						fprintf(fp, "%3d", CBl[j]);
					fprintf(fp, "\nEiVal \n\n", EiVal[i]);
				}
			} //end for
			Cnt++;
			if (PP > 0)
			{
				for (i = 0; i < PP; i++)
				{
					for (j = 0; j < NCtr; j++)
						AtBl[i][j] = PrelCtr[j][i];
				} //end for
				AtBlQnt = PP;
				printf("***AtBlQnt %5d\n", AtBlQnt);

				SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NBOOcc, NBOVec, NBOCtr, NCtr, &IndF);
				DepleteDMNAO(NBOOcc, NBOVec, IndS, IndF);
				TraceDMNAO(BSz, DResid);
				IndS = IndF;
				*NBOAmnt = IndF;
				//				printf("PP= %4d, IndS= %4d, IndF= %4d, NBOAmnt= %4d\n",PP,IndS,IndF,*NBOAmnt);
				PP = 0;
			}
			else
				break;
			if (Cnt > 20)
				break;
		} //while
	}	 //end for
	fclose(fp);
}

//!************************************************************************
void SortPrel(double PrelOcc[PNMax], double PrelVec[PBSz][PNMax], int PrelCtr[PNAt][PNMax],
			  int PP, double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int NCtr, int *IndF)
{

	int i, ii, j, jj, k, l, Cnt1;
	double ThrOcc;
	int PrelAccpt[PNMax];

	ThrOcc = 0.0;

	printf(" SortPrel OK\n");
	for (i = 0; i < PP; i++)
	{
		if ((int)(PrelOcc[i] * 100000) > (int)(ThrOcc * 100000))
		{
			ThrOcc = (int)(PrelOcc[i] * 100000) / 100000.0;
			Cnt1 = 0;
			PrelAccpt[Cnt1++] = i;
			printf("A Ctr=");
			for (jj = 0; jj < NCtr; jj++)
				printf("%4d", PrelCtr[jj][i] + 1);
			printf(" ***ThrOcc= %.5f\n", ThrOcc);
		}
		else if ((int)(PrelOcc[i] * 100000) == (int)(ThrOcc * 100000))
		{
			PrelAccpt[Cnt1++] = i;
			printf("A Ctr=");
			for (jj = 0; jj < NCtr; jj++)
				printf("%4d", PrelCtr[jj][i] + 1);
			printf(" ***ThrOcc= %.5f\n", ThrOcc);
		}
	} // end for

	for (j = 0; j < Cnt1; j++)
	{
		i = PrelAccpt[j];
		NBOOcc[*IndF] = PrelOcc[i];
		for (jj = 0; jj < BSz; jj++)
			NBOVec[jj][*IndF] = PrelVec[jj][i];
		for (jj = 0; jj < NCtr; jj++)
			NBOCtr[jj][*IndF] = PrelCtr[jj][i];
		NBOCtr[jj][*IndF] = -1;
		*IndF += 1;
	}
}

//!************************************************************************
void DepleteDMNAO(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int IndS, int IndF)
{
	int i, j, k, l;

	printf("DeleteDMNAO OK\n");
	for (l = IndS; l < IndF; l++)
	{
		for (i = 0; i < BSz; i++)
		{
			for (j = 0; j < BSz; j++)
			{
				DMNAO[i][j] -= NBOOcc[l] * NBOVec[j][l] * NBOVec[i][l];
			}
		}
	}
}

//!************************************************************************
void TraceDMNAO(int BSz, double *DResid)
{
	int i, j;

	*DResid = 0.0;
	for (i = 0; i < BSz; i++)
		*DResid += DMNAO[i][i];
}

//!************************************************************************
void BlockDMNAO(int CBl[PNAt], int NCtr, double DUMMY[PBSz][PBSz])
{
	int i, j, l, n1, n2;
	int flag[PBSz];

	for (i = 0; i < PBSz; i++)
		flag[i] = 0;
	for (l = 0; l < NCtr; l++)
	{
		n1 = AtBsRng[CBl[l]][0];
		n2 = AtBsRng[CBl[l]][1];
		for (i = n1; i < n2; i++)
			flag[i] = 1;
	}
	for (i = 0; i < BSz; i++)
	{
		for (j = 0; j < BSz; j++)
		{
			if (flag[i] && flag[j])
				DUMMY[i][j] = DMNAO[i][j];
			else
				DUMMY[i][j] = 0.0;
		}
	}
}

//!************************************************************************
void EigenSystem(double a[PBSz][PBSz], int n, double d[PBSz], double v[PBSz][PBSz])
{
	//! a  - real symmetric matrix 实对称矩阵
	//! n  - size of a	矩阵大小
	//! np - size of the physical array storing a - replaced by PBSz
	//! d  - returns evalues of a in first n elements  特征值
	//! v  - contains normalized evectors of a in columns */ 归一化的特征向量

	int i, ip, iq, j;
	double c, g, h, s, sm, t, tau, theta, tresh, b[PBSz], z[PBSz];

	//把v定义为一个单位矩阵
	for (ip = 0; ip < n; ip++)
	{
		for (iq = 0; iq < n; iq++)
		{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}

	//b为a的对角的值 d也为对角的值
  	for (ip = 0; ip < n; ip++)
	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
	}


	for (i = 1; i <= 100; i++)
	{
		sm = 0.0;
		//求矩阵的绝对值的和
		for (ip = 0; ip < n - 1; ip++)
			for (iq = ip + 1; iq < n; iq++)
				sm += fabs(a[ip][iq]);
		//如果为0矩阵
		if (sm == 0.0)
			return;
		//		printf("EigenSystem: loop %4d, sum= %lf\n",i,sm);

		//前3次
		if (i < 4)
			//平均值的0.2
			tresh = 0.2 * sm / (n * n);
		else
			tresh = 0.0;

		//上半矩阵
		for (ip = 0; ip < n - 1; ip++)
		{
			for (iq = ip + 1; iq < n; iq++)
			{
				g = 100.0 * fabs(a[ip][iq]);
				if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh)
				{
					h = d[iq] - d[ip];

					if (fabs(h) + g == fabs(h))
					{
						t = a[ip][iq] / h;
					}
					else
					{
						//t为tan(4)
						theta = 0.5 * h / a[ip][iq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);	//cos(4)
					s = t * c;		//sin(4)
					tau = s / (1.0 + c); //tan(4/2)
					h = t * a[ip][iq];	//tan(4)*a(pq)
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;

					//
					for (j = 0; j <= ip - 1; j++)
					{
						g = a[j][ip];
						h = a[j][iq];
						a[j][ip] = g - s * (h + g * tau);
						a[j][iq] = h + s * (g - h * tau);
					}
					for (j = ip + 1; j <= iq - 1; j++)
					{
						g = a[ip][j];
						h = a[j][iq];
						a[ip][j] = g - s * (h + g * tau);
						a[j][iq] = h + s * (g - h * tau);
					}
					for (j = iq + 1; j < n; j++)
					{
						g = a[ip][j];
						h = a[iq][j];
						a[ip][j] = g - s * (h + g * tau);
						a[iq][j] = h + s * (g - h * tau);
					}
					for (j = 0; j < n; j++)
					{
						g = v[j][ip];
						h = v[j][iq];
						v[j][ip] = g - s * (h + g * tau);
						v[j][iq] = h + s * (g - h * tau);
					}
				} //endif
			}
		}

		for (ip = 0; ip < n; ip++)
		{
			b[ip] = b[ip] + z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	} //end for
	printf("\n******************\nToo many iterations in jacobi\n*******************\n");
	//     pause 'too many iterations in jacobi'
}

//!************************************************************************
void EigenSrt(double d[PBSz], double v[PBSz][PBSz], int n)
{
	int i, j, k;
	double p;

	for (i = 0; i < n - 1; i++)
	{
		k = i;
		p = d[i];
		for (j = i + 1; j < n; j++)
		{
			if (d[j] >= p)
			{
				k = j;
				p = d[j];
			}
		}

		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for (j = 0; j < n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	} //end for
} //return

//!************************************************************************
void Subsets(int AtBl[100000][PNAt], int *NLn, int NClmn, int NAt)
{
	//AtBlQnt->NLn , NCtr->NClmn
	int i, j, k, n, Done;

	//到最后一个原子
	if (NClmn == NAt)
	{
		for (j = 0; j < NAt; j++)
			AtBl[0][j] = j;
		AtBl[0][j] = -1;
		//AtBlQnt = 1
		*NLn = 1;
		return;
	}
	//第一个原子
	if (NClmn == 1)
	{
		for (i = 0; i < NAt; i++)
		{
			AtBl[i][0] = i;
			AtBl[i][1] = -1;
		}
		*NLn = NAt;
		return;
	}

	Done = 0;
	*NLn = 1;
	for (j = 0; j < NClmn; j++)
		AtBl[*NLn - 1][j] = j;

	while (1)
	{
		for (k = 0; k < NClmn; k++)
		{
			if (AtBl[*NLn - 1][k] == NAt - NClmn + k)
			{
				Done = k + 1;
				break;
			}
		}
		if (Done == 1)
			return;
		if (Done == 0)
		{
			*NLn += 1;
			for (k = 0; k < NClmn - 1; k++)
				AtBl[*NLn - 1][k] = AtBl[*NLn - 2][k];
			if (AtBl[*NLn - 2][NClmn - 1] < NAt)
				AtBl[*NLn - 1][NClmn - 1] = AtBl[*NLn - 2][NClmn - 1] + 1;
		}
		else
		{
			*NLn += 1;
			for (k = 0; k < Done - 2; k++)
				AtBl[*NLn - 1][k] = AtBl[*NLn - 2][k];
			AtBl[*NLn - 1][Done - 2] = AtBl[*NLn - 2][Done - 2] + 1;
			for (k = Done - 1; k < NClmn; k++)
				AtBl[*NLn - 1][k] = AtBl[*NLn - 1][k - 1] + 1;
			Done = 0;
		}
	} //while
}

//!************************************************************************
void Input(double NAOAO[PBSz][PBSz])
{

	char dtln[256], dtln1[256], nbofile[256];
	int ioflg;
	int Resid, Cnt;
	int i, j, k;
	FILE *fp;

	if ((fp = fopen("AdNDP.in", "r")) == 0)
	{
		printf("File AdNDP.in not found!\n");
		exit(0);
	}
	printf("AdNDP.in file is OK\n");
	fgets(dtln, 256, fp);		 //NBO filename
	fscanf(fp, "%s\n", nbofile); //nbo.log
	fgets(dtln, 256, fp);		 //Number of atoms
	fscanf(fp, "%d\n", &NAt);	//20
	fgets(dtln, 256, fp);		 //Amount of valence electronic pairs	价电子对数
	fscanf(fp, "%d\n", &NVal);   //454
	fgets(dtln, 256, fp);		 //Total amount of electronic pairs	电子总数
	fscanf(fp, "%d\n", &NTot);   //15
	fgets(dtln, 256, fp);		 //Total amount of basis functions 基函数总数
	fscanf(fp, "%d\n", &BSz);	//21
	fgets(dtln, 256, fp);		 //Amount of basis functions on each atom 每个原子的基函数总数
	for (i = 0; i < NAt; i++)
	{
		fscanf(fp, "%d\n", &AtBs[i]); //19 2
	}
	fgets(dtln, 256, fp); //Occupation number thresholds 阈值
	for (i = 0; i < NAt; i++)
	{
		fscanf(fp, "%lf\n", &Thr[i]); //0.1 0.1
	}
	fgets(dtln, 256, fp);		   //CMO filename
	fscanf(fp, "%s\n", CMOFile);   //mo-03.log
	fscanf(fp, "%s\n", AdNDPFile); //adndp.log
	fclose(fp);

	j = 0;
	for (i = 0; i < NAt; i++)
	{
		AtBsRng[i][0] = j;
		j += AtBs[i];
		AtBsRng[i][1] = j;
	}
	if (BSz != j)
	{ //BSz=21
		printf("Size of basis set error!! %5d != %5d\n", BSz, j);
		exit(0);
	}

	fp = fopen(nbofile, "r"); //nob.log
	if (fp == 0)
	{
		printf("Cannot open the NBO file: %s\n", nbofile);
		exit(0);
	}

	while (
		(dtln, " NAO density matrix:", 11))
	{ //定位到nbo.log的366行
		fgets(dtln, 256, fp);
	}
	fgets(dtln, 256, fp); //
	fgets(dtln, 256, fp); // NAO  1  2  3
	fgets(dtln, 256, fp); // -----------

	Resid = BSz; //21
	k = 0;
	while (Resid > 0)
	{
		if (Resid < 8)
			Cnt = Resid;
		else
			Cnt = 8;
		for (i = 0; i < BSz; i++) //行号
		{
			fgets(dtln, 16, fp);	  //   1. Na 1 (S)  留一个空格
			for (j = 0; j < Cnt; j++) //列号
			{
				//DMNAO 1024*1024
				fscanf(fp, "%lf", &DMNAO[i][8 * k + j]);
				//				printf("%10.6f",DMNAO[i][8*k+j]);
			}
			fgets(dtln, 16, fp);
			//			putchar('\n');
		}
		fgets(dtln, 256, fp);
		fgets(dtln, 256, fp);
		fgets(dtln, 256, fp);
		Resid -= Cnt;
		k++;
	}
	rewind(fp); //指向文件开头
	while (strncmp(dtln, " NAOs in the AO basis:", 12))
	{ //291行
		fgets(dtln, 256, fp);
	}
	fgets(dtln, 256, fp);
	fgets(dtln, 256, fp);
	fgets(dtln, 256, fp);

	Resid = BSz; //21
	k = 0;
	while (Resid > 0)
	{
		if (Resid < 8)
			Cnt = Resid;
		else
			Cnt = 8;

		for (i = 0; i < BSz; i++)
		{
			fgets(dtln, 17, fp);
			for (j = 0; j < Cnt; j++)
			{
				fscanf(fp, "%lf", &NAOAO[i][8 * k + j]);
				//				printf("%10.6f",NAOAO[i][8*k+j]);
			}
			fgets(dtln, 17, fp);
			//			putchar('\n');
		}
		fgets(dtln, 256, fp);
		fgets(dtln, 256, fp);
		fgets(dtln, 256, fp);
		Resid -= Cnt;
		k++;
	}
	printf("Read DMNAO NAOAO matrixs in nbofile OK!\n");
} //over

typedef struct
{
	int nc;
	int ctr[PNAt];
	double occ;
	double vec[PBSz];
	double vocc[PNAt];
} BOND;

void calbond(BOND *b)
{
	int i, j, m, n1, n2;
	double v;

	for (i = 0; i <= NAt; i++)
	{
		if (b->ctr[i] == -1)
		{
			b->nc = i;
			break;
		}
	}

	for (i = 0; i < b->nc; i++)
	{
		m = b->ctr[i];
		n1 = AtBsRng[m][0];
		n2 = AtBsRng[m][1];
		for (j = n1, v = 0.0; j < n2; j++)
			v += b->occ * b->vec[j] * b->vec[j];
		b->vocc[i] = v;
	}
}

//!************************************************************************
void Output(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOCtr[PNAt][PNMax], int NBOAmnt, double DResid, double NAOAO[PBSz][PBSz])
{
	char outfile[256];
	int i, j, k, l, m;
	FILE *fp;
	BOND b;

	strcpy(outfile, "NDP_nbo.log");

	fp = fopen(outfile, "w");
	fprintf(fp, "NBOAmnt=%4d\n", NBOAmnt);

	for (i = 0; i < NBOAmnt; i++)
	{
		b.occ = NBOOcc[i];
		for (j = 0; j <= NAt; j++)
			b.ctr[j] = NBOCtr[j][i];
		for (j = 0; j < BSz; j++)
			b.vec[j] = NBOVec[j][i];

		calbond(&b);

		//		fprintf(fp,"[%4d] Occ= %.4f, Atom ID:",i+1,NBOOcc[i]);
		fprintf(fp, "bond%3d (%2d center, %.3f|e|) ***:", i + 1, b.nc, b.occ);
		for (j = 0; j < b.nc; j++)
		{
			//			k = NBOCtr[j][i];
			//			if(k == -1) break;
			fprintf(fp, "%3d(%.3f)  ", b.ctr[j] + 1, b.vocc[j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\nResidual Density: %10.6f\n", DResid);
	fclose(fp);

	fp = fopen("AdNDP_Resid.log", "wb");

	for (i = 0; i < PBSz; i++)
	{
		fwrite(DMNAO[i], sizeof(double), PBSz, fp);
		fwrite(NAOAO[i], sizeof(double), PBSz, fp);
	}
	for (i = 0; i < PNAt; i++)
		fwrite(AtBsRng[i], sizeof(int), 2, fp);

	fwrite(&BSz, sizeof(int), 1, fp);
	fclose(fp);
}

//!************************************************************************
void BasisChange(double Old[PBSz][PNMax], double New[PBSz][PNMax], double Trans[PBSz][PBSz], int Row, int Col)
{
	int i, j, k;

	for (k = 0; k < Col; k++)
		for (i = 0; i < Row; i++)
			for (j = 0; j < Row; j++)
				New[j][k] += Old[i][k] * Trans[j][i];
}

//!************************************************************************
void NBOPlotMolden(double NBOOcc[PNMax], double NBOVec[PBSz][PNMax], int NBOAmnt)
{
	int i, j, k;
	int Resid, Cnt;
	char dtln[256];
	FILE *fp1, *fp2;

	if ((fp1 = fopen(CMOFile, "r")) == NULL)
	{
		printf("Cannot open the fine %s\n", CMOFile);
		return;
	}
	if ((fp2 =

			 (AdNDPFile, "w")) == NULL)
	{
		printf("Cannot open the fine %s\n", AdNDPFile);
		return;
	}

	while (1)
	{
		fgets(dtln, 256, fp1);
		fputs(dtln, fp2);
		if (strncmp(dtln, "     Molecular Orbital Coefficients", 35) == 0)
			break;
	}

	Cnt = 0;
	Resid = NBOAmnt;
	k = 0;
	while (Resid > 0)
	{
		fgets(dtln, 256, fp1);
		fputs(dtln, fp2);
		fgets(dtln, 256, fp1);
		fputs(dtln, fp2);

		fgets(dtln, 256, fp1);
		dtln[21] = '\0';
		fprintf(fp2, "%s%10.5f%10.5f%10.5f%10.5f%10.5f\n", dtln,
				NBOOcc[5 * k], NBOOcc[5 * k + 1], NBOOcc[5 * k + 2], NBOOcc[5 * k + 3], NBOOcc[5 * k + 4]);

		for (i = 0; i < BSz; i++)
		{
			fgets(dtln, 256, fp1);
			dtln[21] = '\0';
			fprintf(fp2, "%s%10.5f%10.5f%10.5f%10.5f%10.5f\n", dtln,
					NBOVec[i][5 * k], NBOVec[i][5 * k + 1], NBOVec[i][5 * k + 2], NBOVec[i][5 * k + 3], NBOVec[i][5 * k + 4]);
		}
		k++;
		Resid -= 5;
	}
	while (fgets(dtln, 256, fp1))
		fputs(dtln, fp2);
	fclose(fp1);
	fclose(fp2);
}
