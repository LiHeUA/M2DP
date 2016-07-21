/*======================================================================*/ 
/*	Count number of points in bins					*/ 	  
/*	Takes 2D points (in spherical coordinates) and the bin 		*/ 	  
/*	parameters as input; output the number of points in bins	*/
/*	An analogue of Matlab built-in function histcounts2		*/
/*									*/
/*	Author: Li He, Dept. of Computing Science, University of Alberta*/
/*	lhe2@ualberta.ca						*/
/*======================================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*
// Count points in bins
// Input:
//			pTheta			pointer, indicate theta of one data
//			pRho			pointer, indicate rho of one data
//			pTList			pointer, Theta (circle) list
//			pRList			pointer, Radius list
//			numT			number of Theta (circle)
//			numR			number of Radius
//			numData			number of data
//			pVoteCnt		output
//			pIndexPoint		output
// Output:
//			pVoteCnt		pointer, number of points in bins
//			pIndexPoint		pointer, bin index of data belonging to
*/
int CountVote2D(double *pTheta, double *pRho, double *pTList, double *pRList, int numT, int numR, int numData, long int *pVoteCnt, double *pIndexPoint)
{
	double *pT, *pR;
	int i,n;
	int idxR, idxT;

	pT = pTheta;
	pR = pRho;

	/* loop on data */
	for (n=0; n<numData; n++)
	{
		/* find which theta bin this point belonging to */ 
		idxT = numT-1;
		for (i=0;i<numT-1;i++)
			if(*pT<pTList[i+1])
			{
				idxT = i;
				break;
			}
		/* find which radius bin this point belonging to */ 
		idxR = numR-1;
		for (i=0;i<numR-1;i++)
			if(*pR<pRList[i+1])
			{
				idxR = i;
				break;
			}
		
		/* count = count + 1 */
		(pVoteCnt[idxR*numT+idxT])++;
		/* record which bin this point in */	
		pIndexPoint[n] = idxR*numT+idxT;

		/* move to next point */
		pT++;
		pR++;
	}
	return 1;
}


/*
// Matlab interface
// Input:
//			prhs[0]			n*1			theta of points, azimuth
//			prhs[1]			n*1			distances of points to origin, rho
//			prhs[2]			(dt+1)*1		Theta bins, dt: number of bins along theta (Azimuth)
//			prhs[3]			(dr+1)*1		radius bins, dr: number of bins along radius
// Output:
//			plhs[0]			(dr*dt)*1		votes in bins, size = #Radius * #Theta
//			plhs[1]			n*1			bin index of each point belonging to 
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	int numData = mxGetM(prhs[0]); /* number of data */
	double *pIndexPoint; /* bin index of data */
	double *pVote; /* select one bin */
	double *pTList, *pRList; /* Theta (circle) list and Radii list */
	int numR, numT; /* number of Thetas and Rhos */
	long int *pVoteCnt; /* number of points in one bin */
	int i;
	double *pTheta, *pRho; /* data, Theta (data) and Rho */

	/* get data, Theta (data) and Rho */
	pTheta = (double *)mxGetPr(prhs[0]);
	pRho = (double *)mxGetPr(prhs[1]);

	/* get bin parameters, Theta (circle) and Radius */
	pTList = (double *)mxGetPr(prhs[2]);
	pRList = (double *)mxGetPr(prhs[3]);
	
	/* get number of Theta (circle) and number of Radius */
	numT = mxGetM(prhs[2])-1;
	numR = mxGetM(prhs[3])-1;

	/* creat output */
	plhs[0] = mxCreateDoubleMatrix(numR*numT, 1, mxREAL);
	pVote = (double*)mxGetPr(plhs[0]); /* get the start of output */
	plhs[1] = mxCreateDoubleMatrix(numData, 1, mxREAL);
	pIndexPoint = (double*)mxGetPr(plhs[1]); /* get the start of bin index of data */
	
	/* malloc vote count */
	pVoteCnt = (long int*)malloc(numR*numT*sizeof(long int));
	if (pVoteCnt==NULL)
	{
		printf("Not enough memory for pVoteCnt!\n");
		exit(-1);
	}
	memset(pVoteCnt,0,numR*numT*sizeof(long int));	/* initialize vote by zero */

	/* main function, count vote */
	CountVote2D(pTheta, pRho, pTList, pRList, numT, numR, numData, pVoteCnt, pIndexPoint);

	/* normalize counts */
	for (i=0;i<numR*numT;i++)
		pVote[i] = (double)(pVoteCnt[i])/numData;

	free(pVoteCnt);
}
