/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GLinear.h"
#include "GTwt.h"
#include "GTransform.h"
#include "GRand.h"
#include "GVec.h"

using namespace GClasses;

GLinearRegressor::GLinearRegressor(GRand* pRand)
: m_pRand(pRand), m_pPCA(NULL)
{
}

GLinearRegressor::GLinearRegressor(GTwtNode* pNode, GRand* pRand)
: m_pRand(pRand)
{
	m_pPCA = new GPCA(pNode->field("pca"), pRand);
}

// virtual
GLinearRegressor::~GLinearRegressor()
{
	delete(m_pPCA);
}

// virtual
GTwtNode* GLinearRegressor::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GLinearRegressor");
	pNode->addField(pDoc, "pca", m_pPCA->toTwt(pDoc));
	return pNode;
}

// virtual
int GLinearRegressor::featureDims()
{
	if(!m_pPCA)
		ThrowError("Not trained yet");
	return m_pPCA->targetDims();
}

// virtual
int GLinearRegressor::labelDims()
{
	if(!m_pPCA)
		ThrowError("Not trained yet");
	return 1;
}

// virtual
void GLinearRegressor::train(GData* pData, int labelDims)
{
	if(labelDims != 1)
		ThrowError("GLinearRegressor only supports 1 label dim");
	clear();
	m_pPCA = new GPCA(pData->cols() - 1, m_pRand);
	m_pPCA->train(pData);
}

// virtual
void GLinearRegressor::predictDistribution(const double* pIn, GPrediction* pOut)
{
	ThrowError("Sorry, predictDistribution is not implemented for this class. Use \"predict\" instead.");
}
#include <math.h>
// virtual
void GLinearRegressor::predict(const double* pIn, double* pOut)
{
	int inputs = featureDims();
	*pOut = m_pPCA->mean()[inputs];
	for(int i = 0; i < inputs; i++)
	{
		double dp = GVec::dotProduct(m_pPCA->mean(), pIn, m_pPCA->basis(i), inputs);
		double mag = GVec::squaredMagnitude(m_pPCA->basis(i), inputs);
		if(mag > 1e-12)
			*pOut += (dp / mag) * m_pPCA->basis(i)[inputs];
	}
}

// virtual
void GLinearRegressor::clear()
{
	delete(m_pPCA);
	m_pPCA = NULL;
}

#ifndef NO_TEST_CODE
// static
void GLinearRegressor::test()
{
	// Train
	GRand prng(0);
	GData set1(4);
	for(int i = 0; i < 1000; i++)
	{
		double* pVec = set1.newRow();
		pVec[0] = prng.uniform();
		pVec[1] = prng.uniform(); // irrelevant attribute
		pVec[2] = prng.uniform();
		pVec[3] = 7.0 * pVec[0] + 3.0 * pVec[2] + 5.0;
	}
	GLinearRegressor lr(&prng);
	lr.train(&set1, 1);

	// Test
	GData set2(4);
	for(int i = 0; i < 1000; i++)
	{
		double* pVec = set2.newRow();
		pVec[0] = prng.uniform();
		pVec[1] = prng.uniform(); // irrelevant attribute
		pVec[2] = prng.uniform();
		pVec[3] = 7.0 * pVec[0] + 3.0 * pVec[2] + 5.0;
	}
	double results;
	lr.accuracy(&set2, &results);
	if(results > 0.01)
		ThrowError("failed");
}
#endif







/*************************
* The following code was derived (with permission) from code by Jean-Pierre Moreau,
* which was posted at: http://jean-pierre.moreau.pagesperso-orange.fr/
* Here is the text of an email exchange with Dr. Moreau:

On 07/28/2010 06:58 AM, Jean-Pierre Moreau wrote:
> I suppose you can use freely the sources given in my website
> on a "as is" basis (with the given reference).
>
> Best regards.
>
>Jean-Pierre Moreau, Paris.
>
-----Message d'origine-----
De : Mike Gashler [mailto:mikegashler@gmail.com] 
Envoyé : mercredi 14 juillet 2010 20:42
À : jpmoreau@wanadoo.fr
Objet : numerical analysis code licensing

Dr. Moreau,

    Thank you for the useful code that you have shared on your web site. 
I would like to know if I may use some of your code in my open source
projects. Your site does not mention any particular license. Since you
posted it on the web, I assume you would like other people to use your code,
but if there is no license, I cannot safely use your code without fear of
legal repercussions. (I would recommend that you might choose the Creative
Commons CC0 license, which essentially puts the code into the public domain
for anyone to use, but any open source license would be welcome.)

I very much appreciate your excellent work,
-Mike Gashler
*/

void simp1(GData& a, int mm, int* ll, int nll, int iabf, int* kp, double* bmax);
void simp2(GData& a, int m, int n, int* l2, int nl2, int* ip, int kp, double* q1);
void simp3(GData& a, int i1, int k1, int ip, int kp);

void simplx(GData& a, int m, int n, int m1, int m2, int* icase, int* izrov, int* iposv)
{
	int m3 = m - (m1 + m2);
	int i, ip, ir, is, k, kh, kp, m12, nl1, nl2;
	GTEMPBUF(int, l1, n + 2);
	GTEMPBUF(int, l2, m + 2);
	GTEMPBUF(int, l3, m + 2);
	double bmax, q1, EPS = 1e-6; // EPS is the absolute precision, which should be adjusted to the scale of your variables.
	nl1 = n;
	for(k = 1; k<=n; k++)
	{
		l1[k] = k;     //Initialize index list of columns admissible for exchange.
		izrov[k] = k;  //Initially make all variables right-hand.
	}
	nl2 = m;
	for(i = 1; i <= m; i++)
	{
		if (a[i + 1][1] < 0.0)
		{
			printf(" Bad input tableau in simplx, Constants bi must be nonnegative.\n");
			return;
		}
		l2[i] = i;
		iposv[i] = n + i;
		// ------------------------------------------------------------------------------------------------
		// Initial left-hand variables. m1 type constraints are represented by having their slackv ariable
		// initially left-hand, with no artificial variable. m2 type constraints have their slack
		// variable initially left-hand, with a minus sign, and their artificial variable handled implicitly
		// during their first exchange. m3 type constraints have their artificial variable initially
		// left-hand.
		// ------------------------------------------------------------------------------------------------
	}
	for(i = 1; i <= m2; i++)
		l3[i] = 1;
	ir = 0;
	if(m2 + m3 == 0)
		goto e30; //The origin is a feasible starting solution. Go to phase two.
	ir = 1;
	for(k = 1; k <= n + 1; k++)
	{ //Compute the auxiliary objective function.
		q1 = 0.0;
		for(i = m1 + 1; i <= m; i++)
			q1 += a[i + 1][k];
		a[m + 2][k] = -q1;
	}
e10:
	simp1(a,m+1,l1,nl1,0,&kp,&bmax);    //Find max. coeff. of auxiliary objective fn
	if(bmax <= EPS && a[m + 2][1] < -EPS)
	{
		*icase=-1;    //Auxiliary objective function is still negative and cant be improved,
		return;       //hence no feasible solution exists.
	}
	else if (bmax <= EPS && a[m + 2][1] <= EPS)
	{
		//Auxiliary objective function is zero and cant be improved; we have a feasible starting vector.
		//Clean out the artificial variables corresponding to any remaining equality constraints by
		//goto 1s and then move on to phase two by goto 30.
		m12 = m1 + m2 + 1;
		if(m12 <= m)
			for (ip=m12; ip<=m; ip++)
				if(iposv[ip] == ip+n)
				{   //Found an artificial variable for an equalityconstraint.
					simp1(a,ip,l1,nl1,1,&kp,&bmax);
					if(bmax > EPS)
						goto e1; //Exchange with column corresponding to maximum
				} //pivot element in row.
		ir=0;
		m12=m12-1;
		if (m1+1 > m12) goto e30; 
		for (i=m1+1; i<=m1+m2; i++)     //Change sign of row for any m2 constraints
		if(l3[i-m1] == 1)             //still present from the initial basis.
		for (k=1; k<=n+1; k++) 
			a[i + 1][k] *= -1.0; 
		goto e30;                       //Go to phase two.
	}
	simp2(a,m,n,l2,nl2,&ip,kp,&q1); //Locate a pivot element (phase one).

	if(ip == 0)
	{                         //Maximum of auxiliary objective function is
		*icase=-1;                          //unbounded, so no feasible solution exists
		return; 
	} 
e1:
	simp3(a, m + 1, n, ip, kp);
	//Exchange a left- and a right-hand variable (phase one), then update lists. 
	if(iposv[ip] >= n + m1 + m2 + 1)
	{ //Exchanged out an artificial variable for an equality constraint. Make sure it stays out by removing it from the l1 list.
		for (k=1; k<=nl1; k++)
			if(l1[k] == kp)
				break;
		nl1 = nl1 - 1;
		for(is = k; is <= nl1; is++)
			l1[is] = l1[is+1];
	}
	else
	{
		if(iposv[ip] < n+m1+1) goto e20;
		kh=iposv[ip]-m1-n; 
		if(l3[kh] == 0) goto e20;  //Exchanged out an m2 type constraint. 
		l3[kh]=0;                  //If it's the first time, correct the pivot column 
					//or the minus sign and the implicit 
					//artificial variable. 
	}
	a[m + 2][kp + 1] += 1.0;
	for (i=1; i<=m+2; i++)
		a[i][kp + 1] *= -1.0;
e20:
	is=izrov[kp];             //Update lists of left- and right-hand variables. 
	izrov[kp]=iposv[ip];
	iposv[ip]=is;
	if (ir != 0) goto e10;       //if still in phase one, go back to 10. 

	//End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 
e30:
	simp1(a,0,l1,nl1,0,&kp,&bmax); //Test the z-row for doneness.
	if(bmax <= EPS)
	{          //Done. Solution found. Return with the good news.
		*icase=0;
		return;
	}
	simp2(a,m,n,l2,nl2,&ip,kp,&q1);   //Locate a pivot element (phase two).
	if(ip == 0)
	{             //Objective function is unbounded. Report and return.
		*icase=1;
		return;
	} 
	simp3(a,m,n,ip,kp);       //Exchange a left- and a right-hand variable (phase two),
	goto e20;                 //update lists of left- and right-hand variables and
}                           //return for another iteration.


void simp1(GData& a, int mm, int* ll, int nll, int iabf, int* kp, double* bmax)
{
//Determines the maximum of those elements whose index is contained in the supplied list
//ll, either with or without taking the absolute value, as flagged by iabf. 
	int k;
	double test;
	*kp=ll[1];
	*bmax=a[mm + 1][*kp + 1];
	if (nll < 2) return;
	for (k=2; k<=nll; k++)
	{
		if(iabf == 0)
			test=a[mm + 1][ll[k]+1]-(*bmax);
		else
			test=fabs(a[mm + 1][ll[k] + 1]) - fabs(*bmax);
		if(test > 0.0)
		{ 
			*bmax=a[mm + 1][ll[k]+1];
			*kp=ll[k];
		}
	}
	return;
}

void simp2(GData& a, int m, int n, int* l2, int nl2, int* ip, int kp, double* q1)
{
	double EPS=1e-6;
	//Locate a pivot element, taking degeneracy into account.
	int i,ii,k;
	double q;
	double q0 = 0.0;
	double qp = 0.0;
	*ip=0;
	if(nl2 < 1)
		return;
	for (i=1; i<=nl2; i++)
		if (a[i+1][kp+1] < -EPS)
			goto e2;
	return;  //No possible pivots. Return with message.
e2:
	*q1 = -a[l2[i] + 1][1] / a[l2[i] + 1][kp + 1];
	*ip=l2[i];
	if (i+1 > nl2)
		return;
	for (i=i+1; i<=nl2; i++)
	{
		ii=l2[i];
		if(a[ii+1][kp+1] < -EPS)
		{
			q=-a[ii+1][1]/a[ii+1][kp+1];
			if (q <  *q1)
			{
				*ip=ii;
				*q1=q;
			}
			else if (q == *q1)
			{  //We have a degeneracy.
				for (k=1; k<=n; k++)
				{
					qp=-a[*ip+1][k+1]/a[*ip+1][kp+1];
					q0=-a[ii+1][k+1]/a[ii+1][kp+1];
					if (q0 != qp)
						goto e6;
				}
e6:
				if (q0 < qp)
					*ip=ii;
			}
		}
	}
	return;
}

void simp3(GData& a, int i1, int k1, int ip, int kp)
{
	int ii,kk;
	double piv;
	piv=1.0/a[ip+1][kp+1];
	if (i1 >= 0)
		for(ii = 1; ii <= i1 + 1; ii++)
			if (ii-1 != ip)
			{
				a[ii][kp+1] *= piv;
				for (kk=1; kk<=k1+1; kk++)
					if (kk-1 != kp)
						a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
			}
	for (kk=1; kk<=k1+1; kk++)
		if(kk-1 !=  kp) a[ip+1][kk] = -a[ip+1][kk]*piv;
			a[ip+1][kp+1]=piv;
	return;
}

/*
* End of code derived from Dr. Moreau's numerical analysis code
**************************/




// static
bool GLinearProgramming::simplexMethod(GData* pA, const double* pB, int leConstraints, int geConstraints, const double* pC, double* pOutX)
{
	// Set up the matrix in the expected form
	if((size_t)leConstraints + (size_t)geConstraints > pA->rows())
		ThrowError("The number of constraints must be >= leConstraints + geConstraints");
	GData aa(pA->cols() + 2);
	aa.newRows(pA->rows() + 3);
	aa.setAll(0.0);
	aa[1][1] = 0.0;
	GVec::copy(aa.row(1) + 2, pC, pA->cols());
	for(size_t i = 1; i <= pA->rows(); i++)
	{
		GVec::copy(aa.row(i + 1) + 2, pA->row(i - 1), pA->cols());
		GVec::multiply(aa.row(i + 1) + 2, -1.0, pA->cols());
		aa[i + 1][1] = pB[i - 1];
	}

	// Solve it
	int icase;
	GTEMPBUF(int, iposv, aa.rows());
	GTEMPBUF(int, izrov, aa.cols());
	simplx(aa, pA->rows(), pA->cols(), leConstraints, geConstraints, &icase, izrov, iposv);

	// Extract the results
	if(icase)
		return false; // No solution. (icase gives an error code)
	GVec::setAll(pOutX, 0.0, pA->cols());
	for(size_t i = 1; i <= pA->rows(); i++)
	{
		int index = iposv[i];
		if(index >= 1 && index <= pA->cols())
			pOutX[index - 1] = aa[i + 1][1];
	}
	return true;
}


#ifndef NO_TEST_CODE
// static
void GLinearProgramming::test()
{
	GData a(4);
	a.newRows(4);
	a[0][0] = 1.0; a[0][1] = 0.0; a[0][2] = 2.0; a[0][3] = 0.0;
	a[1][0] = 0.0; a[1][1] = 2.0; a[1][2] = 0.0; a[1][3] = -7.0;
	a[2][0] = 0.0; a[2][1] = 1.0; a[2][2] = -1.0; a[2][3] = 2.0;
	a[3][0] = 1.0; a[3][1] = 1.0; a[3][2] = 1.0; a[3][3] = 1.0;
	double b[4]; b[0] = 740.0; b[1] = 0.0; b[2] = 0.5; b[3] = 9.0;
	double c[4]; c[0] = 1.0; c[1] = 1.0; c[2] = 3.0; c[3] = -0.5;
	double x[4];
	if(!GLinearProgramming::simplexMethod(&a, b, 2, 1, c, x))
		ThrowError("failed to find a solution");
	if(ABS(0.0 - x[0]) > 1e-6)
		ThrowError("failed");
	if(ABS(3.325 - x[1]) > 1e-6)
		ThrowError("failed");
	if(ABS(4.725 - x[2]) > 1e-6)
		ThrowError("failed");
	if(ABS(0.95 - x[3]) > 1e-6)
		ThrowError("failed");
}
#endif
