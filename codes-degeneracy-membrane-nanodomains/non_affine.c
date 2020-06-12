/*********************************************************************************************************/
/* Calculate the local strain tensor epsilon_ij at every particle/node of a system using configurations
   of the system at two subsequent time steps

   Following the algorithm by Falk and Langer (PRE, 1998)

   Author: Madhusmita Tripathy (madhu.cfl@gmail.com)
   Date: 26th October 2016
   Modified: 14th November 2016
             6th December 2016 (PBC correction for a periodic system) */
/*********************************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define max_nodes     144	        /* Maximum number of nodes at which deformation has to be calculated */

# define dim		2		/* dimension */

# define r_cut		14.0		/* radius of the coarse-graining volume */

//# define box		101		/* box dimension for the periodic system */


void dgetrf_(int*, int *, double*, int*, int*, int*);			/* LU decomoposition of a general matrix */

void dgetri_(int*, double*, int*, int*, double*, int*, int*);		/* Generate inverse of a matrix given its LU decomposition */

void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 	/* Multiply two general matrices */


int N_nodes;								/* Number of neighboring nodes for a given one */
double r_o[max_nodes+1][dim+1], r_n[max_nodes+1][dim+1];		/* Node coordinates: old and new */
double dr_o[dim+1], mod_dr_o;
double rij_old[max_nodes+1][dim+1], rij_new[max_nodes+1][dim+1];	/* Displacement vectors: old and new */

double X[dim+1][dim+1], Y[dim+1][dim+1], eps[dim+1][dim+1], sym_eps[dim+1][dim+1];	/* Matrices for epsilon calculation */
double A[dim*dim], A1[dim*dim], Y_inv[dim+1][dim+1], del[dim+1][dim+1];

double y[dim*dim], y_inv[dim*dim], C[dim*dim];				/* Matrices for inverse calculation */
double v[dim+1];

double bx, by, box[2726][3], dummy;

double sum, chi_square;

int IPIV[dim];
int LWORK = dim*dim;
double WORK[dim*dim];
int INFO;


int main()
{

 int i, j, k, p, q, dum, N = dim, N1 = 1;
 int n, num_nbr;
 double a, b, c;

 char no = 'N';
 double alpha= 1.0, beta= 0.0;

 FILE *fp, *fp1, *fp2, *fp3;
 char file[150], file1[150], file2[150], file3[150];
 int run, maxrun;

run = 2;
maxrun = 2728;        // number of data files
    fp = fopen("box-dimensions.dat","r");
    i = 1;
    while(i <= 2726)
    {
        fscanf(fp, "%lf %lf %lf", &bx, &by, &dummy);
        box[i][1] = bx;
        box[i][2] = by;
        i++;
    }
    fclose(fp);
 while ( run < maxrun )
 {
    /* Open the data file to read in bead positions */
     sprintf(file,"bottom-%d.txt",run);
    fp = fopen(file,"r");    
    
    if(fp == NULL)
    {
        printf("cannot open file %s\n",file);
    }
    printf("Reading file %s\n",file);

 i = 0;
 while(!feof(fp))
 {
	i++;
	fscanf(fp, "%lf %lf %lf", &a, &b, &dummy);
//	fscanf(fp, "%d %lf %lf",  &dum, &a, &b);
	r_o[i][1] = a;
	r_o[i][2] = b;
 }
 fclose(fp);
 N_nodes = i-1;					/* Number of nodes */
    
 




    /* Open the data file to read in bead positions */
     sprintf(file1,"bottom-%d.txt",run+1);
    fp1 = fopen(file1,"r");    
    
    if(fp1 == NULL)
    {
        printf("cannot open file %s\n",file1);
    }
    printf("Reading file %s\n",file1);
 	i = 0;

    while(!feof(fp1))
    { 
	i++;
	fscanf(fp1, "%lf %lf %lf", &a, &b, &dummy);
//	fscanf(fp1, "%d %lf %lf", &dum, &a, &b);
	r_n[i][1] = a;
	r_n[i][2] = b;
     }
 fclose(fp1);
 N_nodes = i-1;

   sprintf(file3,"chi_sq_frames-%d-%d.dat",run, run+1);
   fp3 = fopen(file3,"w");

/* Calculate the local deformation at all nodes using each one of them as reference */
 for(i=1; i<=N_nodes; i++)
 {
	n = 0;						/* Number of neighbors of ith node within cutoff radius */
	for(j=1; j<i; j++)
	{
		mod_dr_o = 0;
		for(k=1; k<=dim; k++)
		{
			dr_o[k] = r_o[j][k] - r_o[i][k];
			dr_o[k] = dr_o[k] - (rint(dr_o[k]/box[run][k])) * box[run][k];

			mod_dr_o = mod_dr_o + (dr_o[k]*dr_o[k]);
		}

		if(sqrt(mod_dr_o) <= r_cut)		/* Pair of nodes within cutoff */
		{
			n++;
			for(k=1; k<=dim; k++)
			{
				rij_old[n][k] = dr_o[k];
				rij_new[n][k] = r_n[j][k] - r_n[i][k];
				rij_new[n][k] = rij_new[n][k] - (rint(rij_new[n][k]/box[run][k])) * box[run][k];	/* PBC correction */
			}
		}
	}
	for(j=i+1; j<=N_nodes; j++)
	{
		mod_dr_o = 0;
		for(k=1; k<=dim; k++)
		{
			dr_o[k] = r_o[j][k] - r_o[i][k];
			dr_o[k] = dr_o[k] - (rint(dr_o[k]/box[run][k])) * box[run][k];

			mod_dr_o = mod_dr_o + (dr_o[k]*dr_o[k]);
		}

		if(sqrt(mod_dr_o) <= r_cut)		/* Pair of nodes within cutoff */
		{
			n++;
			for(k=1; k<=dim; k++)
			{
				rij_old[n][k] = dr_o[k];
				rij_new[n][k] = r_n[j][k] - r_n[i][k];
				rij_new[n][k] = rij_new[n][k] - (rint(rij_new[n][k]/box[run][k])) * box[run][k];	/* PBC correction */
			}
		}
	}
	num_nbr = n;
//	printf("Node: %d\t with number of neighbors: %d\n", i, num_nbr);

	/* Initialise matrices */
 	for(p=1; p<=dim; p++)
 	{
		for(q=1; q<=dim; q++)
		{
			X[p][q] = 0.0;
			Y[p][q] = 0.0;
			eps[p][q] = 0.0;
			del[p][q] = 0.0;
		}

		del[p][p] = 1.0;
 	}

	/* Compute the elements of X and Y using dr_o and dr_n (reference: Falk and Langer PRE 1998)*/
 	for(p=1; p<=dim; p++)
 	{
		for(q=1; q<=dim; q++)
		{
			for(k=1; k<=num_nbr; k++)
			{
				X[p][q] = X[p][q] + rij_new[k][p]*rij_old[k][q];
				Y[p][q] = Y[p][q] + rij_old[k][p]*rij_old[k][q];
			}
		}
 	}

	/* Storing Y matrix as a 1-D column major array for calculating inverse */
	/* To call a Fortran routine from C, we have to transform the matrix */
 	for (p=0; p<dim; p++)
 	{
		for(q=0; q<dim; q++)
		{
			A[q + dim*p] = Y[q+1][p+1];
		}
 	}

	/* Invert the matrix A(=Y) using LAPACK routines */
 	dgetrf_(&N, &N, A, &N, IPIV, &INFO);			/* LAPACK routine for matrix inversion, */
	dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);		/* using LU decomposition */


	/* Transform the column major 1-D array A to dim-D matrix Y_inv */
 	for (p=0; p<dim; p++)
 	{
		for(q=0; q<dim; q++)
		{
			Y_inv[q+1][p+1] = A[q + dim*p];
		}
 	}


	/* Calculate local deformation eps using "eps[p][q] = \sum_k X[p][k]Y_inv[q][k] - del[p][q]" */
	for(p=1; p<=dim; p++)
 	{
		for(q=1; q<=dim; q++)
		{
			for(k=1; k<=dim; k++)
			{
				eps[p][q] = eps[p][q] + X[p][k]*Y_inv[q][k];
			}
			eps[p][q] = eps[p][q] - del[p][q];
			//printf("%lf\t", eps[p][q]);
		}
		//printf("\n");
 	}


	/* Make the local deformation matrix (eps) symmetric (if needed) */
	//printf("\nLocal deformation=\n");
 	for(p=1; p<=dim; p++)
 	{
		for(q=1; q<=dim; q++)
		{
			sym_eps[p][q] = 0.5*(eps[p][q] + eps[q][p]);
			// printf("%lf\t", sym_eps[p][q]);
		}
		// printf("\n");
	}

	/* Calculate the minimum chi-square: local deviation from affineness */
	chi_square = 0.0;
 	for(k=1; k<=num_nbr; k++)
 	{
		for(p=1; p<=dim; p++)
		{
			sum = 0.0;
			for(q=1; q<=dim; q++)
			{
				sum = sum + (del[p][q] + eps[p][q])*rij_old[k][q];
			}
			chi_square =  chi_square + (rij_new[k][p] - sum)*(rij_new[k][p] - sum);
		}
 	}
 
	//fprintf(fp3,"%d\t %lf\t %lf\t %lf\n", i, r_n[i][1], r_n[i][2], chi_square);
    fprintf(fp3,"%lf\t %lf\t %lf\n", r_n[i][1], r_n[i][2], chi_square);
 }
 fclose(fp3);
run++;
}



}
