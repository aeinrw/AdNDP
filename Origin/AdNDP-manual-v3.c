//      program Adaptive Natural Density Partitioning AdNDP
//!
//!     version A.6 May 07, 2009
//!

#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <stdlib.h>
#define NAt 128
#define PNVl 9
#define PNTot 15
#define PNMax 1024
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
int AtBs[NAt], AtBsRng[NAt][2];
double DMNAO[PBSz][PBSz];
double Thr[NAt];

typedef struct
{
	int nc;
	int ctr[NAt];
	double occ;
	double vec[PBSz];
	double vocc[NAt];
} BOND;

void calbond(BOND *b)
{
	int i, j, m, n1, n2;
	double v;

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

void AdNBO(int *MnSrch, BOND bond[PNMax], int *NBOAmnt, double *DResid);
void SortPrel(double PrelOcc[PNMax], double PrelVec[PBSz][PNMax], int PrelCtr[NAt][PNMax], int PP, int NCtr, BOND bond[PNMax]);
void DepleteDMNAO(BOND b);
void TraceDMNAO(int BSz, double *DResid);
void BlockDMNAO(int CBl[NAt], int NCtr, double DUMMY[PBSz][PBSz]);
void EigenSystem(double a[PBSz][PBSz], int n, double d[PBSz], double v[PBSz][PBSz]);
void EigenSrt(double d[PBSz], double v[PBSz][PBSz], int n);
void Subsets(int AtBl[100000][NAt], int *NLn, int NClmn, int NAt);
void Input(double NAOAO[PBSz][PBSz]);
void Output(BOND bond[PNMax], int NBOAmnt, double DResid, double NAOAO[PBSz][PBSz]);
void BasisChange(BOND bond[PNMax], double New[PBSz][PNMax], double Trans[PBSz][PBSz], int Row, int Col);
void NBOPlotMolden(BOND bond[PNMax], double NBOVec[PBSz][PNMax], int NBOAmnt);

void main()
{
	int MnSrch;

	double NBOOcc[PNMax], NBOVec[PBSz][PNMax];
	double NBOVecAO[PBSz][PNMax];
	double DResid;
	int NBOCtr[NAt][PNMax];
	int NBOAmnt;
	double NAOAO[PBSz][PBSz];
	BOND bond[PNMax];

	int i, ii, j, jj, k, l, m, n;

	MnSrch = 0;

	Input(NAOAO);
	//    AdNBO(&MnSrch,bond,&NBOAmnt,&DResid);
	//    Output(bond,NBOAmnt,DResid,NAOAO);
	//    BasisChange(bond,NBOVecAO,NAOAO,BSz,NBOAmnt);
	//    NBOPlotMolden(bond,NBOVecAO,NBOAmnt);
}

//************************************************************************
void AdNBO(int *MnSrch, BOND bond[PNMax], int *NBOAmnt, double *DResid)
{
	double PrelOcc[PNMax], PrelVec[PBSz][PNMax];
	int PrelCtr[NAt][PNMax];

	int i, j, k, l, m, imax, smode;
	int Cnt, PP;
	int IndS, IndF;
	int NCtr, AtBlQnt;
	int AtBl[100000][NAt], CBl[NAt];
	BOND prebond[PNMax], b;
	double threshold, vmax;
	double DUMMY[PBSz][PBSz];
	double EiVal[PBSz], EiVec[PBSz][PBSz];
	FILE *fp;

	fp = fopen("out_debug", "w");

	IndS = 0;
	IndF = 0;

	for (i = 0; i < PBSz; i++)
		for (j = 0; j < PBSz; j++)
			DUMMY[i][j] = 0.0;

	printf("Welcom to AdNDP-manual program v2.0, revised by Longjiu-cheng\n");

	*NBOAmnt = 0;
	while (1)
	{
		NCtr = 0;
		while (1)
		{
			printf("Input the number of centers of the bond to be found (nc = 1,2,3...):");
			scanf("%d", &NCtr);
			if (NCtr < 1 || NCtr > NAt)
				printf("\nInvalid nc number nc=%4d! Reinit please!\n", NCtr);
			else
				break;
		}
		threshold = -1.0;
		while (1)
		{
			printf("Input the threshold value for the %d-center bonds(0.0~2.0): ", NCtr);
			scanf("%lf", &threshold);
			if (threshold < 0 || threshold > 2)
				printf("Invalid threshold value: %lf! Reinit please!\n", threshold);
			else
				break;
		}
		while (1)
		{
			printf("Input the searching mode: 1--direct searching; 2--given centers: ");
			scanf("%d", &smode);
			if (smode != 1 && smode != 2)
		}
		if (smode == 1)
		{
			printf("Searching for all possible %d-center bonds ......\n", NCtr);
			Subsets(AtBl, &AtBlQnt, NCtr, NAt);
		}
		if (smode == 2)
		{
			AtBlQnt = 1;
			printf("\nInput the given %3d centers: ", NCtr);
			for (i = 0; i < NCtr; i++)
			{
				scanf("%d", &j);
				AtBl[0][i] = j - 1;
			}
		}

		printf("NCtr = %3d, AtBlQnt= %6d\n", NCtr, AtBlQnt);

		PP = 0;
		for (k = 0; k < AtBlQnt; k++)
		{
			for (j = 0; j < NCtr; j++)
				CBl[j] = AtBl[k][j];
			BlockDMNAO(CBl, NCtr, DUMMY);
			EigenSystem(DUMMY, BSz, EiVal, EiVec);
			EigenSrt(EiVal, EiVec, BSz);
			for (i = 0; i < BSz; i++)
			{
				if (EiVal[i] >= threshold)
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
			}
		}
		if (PP <= 0)
		{
			printf("No bonds found with occ >= %lf\n", threshold);
		}
		else
		{
			SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NCtr, prebond);
			printf("********\n\t%5d %3d-center bonds found with occ = %lf:\n", PP, NCtr, threshold);
			for (i = 0; i < PP; i++)
			{
				b = prebond[i];
				printf("[%4d]**** occ=%.4f):**  ", i, b.occ);
				for (j = 0; j < b.nc; j++)
				{
					printf("%3d(%.3f)  ", b.ctr[j] + 1, b.vocc[j]);
				}
				printf("\n");
			}
		}
		printf("Please input the bond indexes to be selected(end with -1):");
		while (1)
		{
			scanf("%d", &i);
			if (i < 0)
				break;
			b = prebond[i];
			DepleteDMNAO(b);
			bond[*NBOAmnt] = b;
			*NBOAmnt += 1;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!
		TraceDMNAO(BSz, DResid);
		printf("\n*** %5d bonds found, Density residure = %lf***\n", *NBOAmnt, *DResid);
		printf("0 --- end the AdNDP program\n");
		printf("1 --- redo AdNDP search\n");
		printf("*** Input you selection(0 or 1): ");
		while (1)
		{
			scanf("%d", &i);
			if (i != 0 && i != 1)
			{
				printf("\n*** invalid input, reinput your selection please (0 or 1): ");
				continue;
			}
			else
				break;
		}
		if (i == 1)
			continue;
		else
		{
			close(fp);
			return;
		}
	}
}

//!************************************************************************
void SortPrel(double PrelOcc[PNMax], double PrelVec[PBSz][PNMax], int PrelCtr[NAt][PNMax], int PP, int NCtr, BOND bond[PNMax])
{

	int i, ii, j, jj, k, l, Cnt1;
	double ThrOcc;
	BOND b;
	int PrelAccpt[PNMax];

	printf(" SortPrel OK\n");
	for (j = 0; j < PP; j++)
	{
		for (i = 0, ThrOcc = 0.0; i < PP; i++)
		{
			if (PrelOcc[i] > ThrOcc)
			{
				ThrOcc = PrelOcc[i];
				Cnt1 = i;
			}
		}
		b.nc = NCtr;
		i = Cnt1;
		b.occ = PrelOcc[i];
		PrelOcc[i] = -1.0;
		for (jj = 0; jj < BSz; jj++)
			b.vec[jj] = PrelVec[jj][i];
		for (jj = 0; jj < NCtr; jj++)
			b.ctr[jj] = PrelCtr[jj][i];
		calbond(&b);
		bond[j] = b;
	} // end for
}

//!************************************************************************
void DepleteDMNAO(BOND b)
{
	int i, j;

	for (i = 0; i < BSz; i++)
	{
		for (j = 0; j < BSz; j++)
		{
			DMNAO[i][j] -= b.occ * b.vec[j] * b.vec[i];
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
void BlockDMNAO(int CBl[NAt], int NCtr, double DUMMY[PBSz][PBSz])
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
	//! a  - real symmetric matrix
	//! n  - size of a
	//! np - size of the physical array storing a - replaced by PBSz
	//! d  - returns evalues of a in first n elements
	//! v  - contains normalized evectors of a in columns */

	int i, ip, iq, j;
	double c, g, h, s, sm, t, tau, theta, tresh, b[PBSz], z[PBSz];

	for (ip = 0; ip < n; ip++)
	{
		for (iq = 0; iq < n; iq++)
		{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}

	for (ip = 0; ip < n; ip++)
	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
	}

	for (i = 1; i <= 100; i++)
	{
		sm = 0.0;
		for (ip = 0; ip < n - 1; ip++)
			for (iq = ip + 1; iq < n; iq++)
				sm += fabs(a[ip][iq]);
		if (sm == 0.0)
			return;
		//		printf("EigenSystem: loop %4d, sum= %lf\n",i,sm);

		if (i < 4)
			tresh = 0.2 * sm / (n * n);
		else
			tresh = 0.0;

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
						theta = 0.5 * h / a[ip][iq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;

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
void Subsets(int AtBl[100000][NAt], int *NLn, int NClmn, int NAt)
{
	int i, j, k, n, Done;

	if (NClmn == NAt)
	{
		for (j = 0; j < NAt; j++)
			AtBl[0][j] = j;
		AtBl[0][j] = -1;
		*NLn = 1;
		return;
	}
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
	fgets(dtln, 256, fp);
	fscanf(fp, "%s\n", nbofile);
	fgets(dtln, 256, fp);
	fscanf(fp, "%d\n", &NAt);
	fgets(dtln, 256, fp);
	fscanf(fp, "%d\n", &NVal);
	fgets(dtln, 256, fp);
	fscanf(fp, "%d\n", &NTot);
	fgets(dtln, 256, fp);
	fscanf(fp, "%d\n", &BSz);
	fgets(dtln, 256, fp);
	for (i = 0; i < NAt; i++)
	{
		fscanf(fp, "%d\n", &AtBs[i]);
	}
	fgets(dtln, 256, fp);
	for (i = 0; i < NAt; i++)
	{
		fscanf(fp, "%lf\n", &Thr[i]);
	}
	fgets(dtln, 256, fp);
	fscanf(fp, "%s\n", CMOFile);
	fscanf(fp, "%s\n", AdNDPFile);
	fclose(fp);
	/*
	j=0;
	for(i=0;i<NAt;i++){
		AtBsRng[i][0] = j;
		j += AtBs[i];
		AtBsRng[i][1] = j;
	}
	if(BSz != j){
		printf("Size of basis set error!! %5d != %5d\n",BSz,j);
		exit(0);
	}

	fp=fopen(nbofile,"r");
	if( fp == 0){
		printf("Cannot open the NBO file: %s\n",nbofile);
		exit(0);
	}

	while(strncmp(dtln," NAO density matrix:",11)){
		fgets(dtln,256,fp);
	}
	fgets(dtln,256,fp);
	fgets(dtln,256,fp);
	fgets(dtln,256,fp);

	Resid=BSz;
	k=0;
	while(Resid>0)
	{
		if(Resid<8)
			Cnt=Resid;
		else
			Cnt=8;
		for(i=0;i<BSz;i++)
		{
			fgets(dtln,16,fp);
			for(j=0;j<Cnt;j++)
			{
				fscanf(fp,"%lf",&DMNAO[i][8*k+j]);
//				printf("%10.6f",DMNAO[i][8*k+j]);
			}
			fgets(dtln,16,fp);
//			putchar('\n');
		}
		fgets(dtln,256,fp);
		fgets(dtln,256,fp);
		fgets(dtln,256,fp);
		Resid -= Cnt;
		k++;
	}
	rewind(fp);
	while(strncmp(dtln," NAOs in the AO basis:",12)){
		fgets(dtln,256,fp);
	}
	fgets(dtln,256,fp);
	fgets(dtln,256,fp);
	fgets(dtln,256,fp);

	Resid=BSz;
	k=0;
	while(Resid>0){
		if(Resid<8)
			Cnt=Resid;
		else
			Cnt=8;

		for(i=0;i<BSz;i++){
			fgets(dtln,17,fp);
			for(j=0;j<Cnt;j++) {
				fscanf(fp,"%lf",&NAOAO[i][8*k+j]);
//				printf("%10.6f",NAOAO[i][8*k+j]);
			}
			fgets(dtln,17,fp);
//			putchar('\n');
		}
		fgets(dtln,256,fp);
		fgets(dtln,256,fp);
		fgets(dtln,256,fp);
		Resid -= Cnt;
		k++;
	}*/
	printf("Read DMNAO NAOAO matrixs in nbofile OK!\n");
} //over

//!************************************************************************
void Output(BOND bond[PNMax], int NBOAmnt, double DResid, double NAOAO[PBSz][PBSz])
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
		b = bond[i];

		fprintf(fp, "bond%3d (%2d center, %.3f|e|) ***:", i + 1, b.nc, b.occ);
		for (j = 0; j < b.nc; j++)
		{
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
	for (i = 0; i < NAt; i++)
		fwrite(AtBsRng[i], sizeof(int), 2, fp);

	fwrite(&BSz, sizeof(int), 1, fp);
	fclose(fp);
}

//!************************************************************************
void BasisChange(BOND bond[PNMax], double New[PBSz][PNMax], double Trans[PBSz][PBSz], int Row, int Col)
{
	int i, j, k;

	for (k = 0; k < Col; k++)
		for (i = 0; i < Row; i++)
			New[i][k] = 0.0;
	for (k = 0; k < Col; k++)
		for (i = 0; i < Row; i++)
			for (j = 0; j < Row; j++)
				New[j][k] += bond[k].vec[i] * Trans[j][i];
}

//!************************************************************************
void NBOPlotMolden(BOND bond[PNMax], double NBOVec[PBSz][PNMax], int NBOAmnt)
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
	if ((fp2 = fopen(AdNDPFile, "w")) == NULL)
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
				bond[5 * k].occ, bond[5 * k + 1].occ, bond[5 * k + 2].occ, bond[5 * k + 3].occ, bond[5 * k + 4].occ);

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
