/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GMacros.h"
#include "GEvolutionary.h"
#include "GBits.h"
#include "GRand.h"
#include "GVec.h"
#include <math.h>

using namespace GClasses;

namespace GClasses {
class GEvolutionaryOptimizerNode
{
protected:
	double* m_pVector;
	double m_error;

public:
	GEvolutionaryOptimizerNode(int dims)
	{
		m_pVector = new double[dims];
	}

	virtual ~GEvolutionaryOptimizerNode()
	{
		delete[] m_pVector;
	}

	double* GetVector()
	{
		return m_pVector;
	}

	void SetError(double d)
	{
		m_error = d;
	}

	double GetError()
	{
		return m_error;
	}
};
}

GEvolutionaryOptimizer::GEvolutionaryOptimizer(GTargetFunction* pCritic, int nPopulation, GRand* pRand, double moreFitSurvivalRate)
: GOptimizer(pCritic), m_pRand(pRand)
{
	GAssert(nPopulation >= 2); // can't have a tournament without two members of the population

	// Compute the probability that we need to do a tournament, such that it works out
	// to be equivalent to always performing the tournament and picking the winner with
	// the specified probability. (This saves time by performing fewer tournaments.)
	m_tournamentProbability = 2 * moreFitSurvivalRate - 1;

	int dims = pCritic->relation()->size();
	double* pVec;
	for(int i = 0; i < nPopulation; i++)
	{
		GEvolutionaryOptimizerNode* pNode = new GEvolutionaryOptimizerNode(dims);
		pVec = pNode->GetVector();
		m_pCritic->initVector(pVec);
		m_pCritic->constrain(pVec);
		if(m_pCritic->isStable())
		{
			double error = pCritic->computeError(pVec);
			pNode->SetError(error);
		}
		m_population.push_back(pNode);
	}

	m_bestIndex = 0;
	m_bestErr = 1e308;
}

// virtual
GEvolutionaryOptimizer::~GEvolutionaryOptimizer()
{
	for(size_t i = 0; i < m_population.size(); i++)
		delete(node(i));
}

GEvolutionaryOptimizerNode* GEvolutionaryOptimizer::node(size_t index)
{
	return m_population[index];
}

int GEvolutionaryOptimizer::doTournament()
{
	size_t popSize = m_population.size();
	size_t a = (size_t)m_pRand->next(popSize);
	if(m_pRand->uniform() >= m_tournamentProbability)
		return a;
	size_t b = (size_t)m_pRand->next(popSize - 1);
	if(b >= a)
		b++;
	GEvolutionaryOptimizerNode* pA = node(a);
	GEvolutionaryOptimizerNode* pB = node(b);
	if(!m_pCritic->isStable())
	{
		if(m_pRand->next(popSize) == 0)
			m_bestErr = 1e308;
		recomputeError(a, pA, pA->GetVector());
		recomputeError(b, pB, pB->GetVector());
	}
	if(pA->GetError() >= pB->GetError())
		return a;
	else
		return b;
}

void GEvolutionaryOptimizer::recomputeError(size_t index, GEvolutionaryOptimizerNode* pNode, const double* pVec)
{
	double err = m_pCritic->computeError(pVec);
	pNode->SetError(err);
	if(err < m_bestErr)
	{
		m_bestErr = err;
		m_bestIndex = index;
	}
}

// virtual
double GEvolutionaryOptimizer::iterate()
{
	int dims = m_pCritic->relation()->size();
	int target = doTournament();
	GEvolutionaryOptimizerNode* pNode = node(target);
	int popSize = (int)m_population.size();
	double* pVec = pNode->GetVector();
	int technique = (int)m_pRand->next(8);
	switch(technique)
	{
		case 0: // clone and mutate in all dimensions
			{
				GEvolutionaryOptimizerNode* pParent = node((size_t)m_pRand->next(popSize));
				GVec::copy(pVec, pParent->GetVector(), dims);
				double dev = exp(m_pRand->uniform() * 16.0 - 8.0);
				for(int i = 0; i < dims; i++)
				{
					int vals = m_pCritic->relation()->valueCount(i);
					if(vals == 0)
						pVec[i] += m_pRand->normal() * dev;
					else
					{
						if(m_pRand->next(4) == 0)
							pVec[i] = (double)m_pRand->next(vals);
					}
				}
			}
			break;
		case 1: // random mix
			{
				double* pPar1 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				double* pPar2 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				int i;
				for(i = 0; i < dims; i++)
					pVec[i] = (m_pRand->next(2) == 0 ? pPar1[i] : pPar2[i]);
			}
			break;
		case 2: // single-point cross-over
			{
				double* pPar1 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				double* pPar2 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				unsigned int pivot = (unsigned int)m_pRand->next(dims);
				unsigned int i;
				for(i = 0; i < pivot; i++)
					pVec[i] = pPar1[i];
				for( ; i < (unsigned int)dims; i++)
					pVec[i] = pPar2[i];
			}
			break;
		case 3: // interpolate
			{
				double* pPar1 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				double* pPar2 = (node((size_t)m_pRand->next(popSize)))->GetVector();
				double t = m_pRand->uniform() * 3; // values > 1 will catapult beyond the parent
				int i;
				for(i = 0; i < dims; i++)
				{
					int vals = m_pCritic->relation()->valueCount(i);
					if(vals == 0)
						pVec[i] = t * pPar1[i] + (1.0 - t) * pPar2[i];
					else
					{
						if(t < 0.5)
							pVec[i] = pPar1[i];
						else
							pVec[i] = pPar2[i];
					}
				}
			}
			break;
		default: // clone and mutate in one dimension
			{
				GEvolutionaryOptimizerNode* pParent = node((size_t)m_pRand->next(popSize));
				GVec::copy(pVec, pParent->GetVector(), dims);
				int mutateDim = (int)m_pRand->next(dims);
				int vals = m_pCritic->relation()->valueCount(mutateDim);
				if(vals == 0)
					pVec[mutateDim] = m_pRand->normal() * exp(m_pRand->uniform() * 16.0 - 8.0);
				else
					pVec[mutateDim] = (int)m_pRand->next(vals);
			}
			break;
	}
	m_pCritic->constrain(pVec);
	if(m_pCritic->isStable())
		recomputeError(target, pNode, pVec);
	return m_bestErr;
}

// virtual
double* GEvolutionaryOptimizer::currentVector()
{
	return (node(m_bestIndex))->GetVector();
}
