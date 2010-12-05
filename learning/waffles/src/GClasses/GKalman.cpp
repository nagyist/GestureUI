/*
	Copyright (C) 2009, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GKalman.h"
#include "GData.h"
#include "GVec.h"

using namespace GClasses;

GExtendedKalmanFilter::GExtendedKalmanFilter(int stateDims, int observationDims, int controlDims)
: m_stateDims(stateDims), m_obsDims(observationDims)
{
	m_x = new double[stateDims + observationDims + stateDims];
	m_z = m_x + stateDims;
	m_zz = m_z + observationDims;
	GVec::setAll(m_x, 0.0, stateDims);

	m_pP = new GData(stateDims);
	m_pP->newRows(stateDims);
	m_pP->makeIdentity();
	m_pP->multiply(1000.0); // todo: inefficient since non-diagonal values are known to be 0
}

GExtendedKalmanFilter::~GExtendedKalmanFilter()
{
	delete[] m_x;
	delete(m_pP);
}

void GExtendedKalmanFilter::advance(const double* pControl, GData* pA)
{
	// Check values
	GAssert(pA->rows() == (size_t)m_stateDims && pA->cols() == m_stateDims); // transition Jacobian wrong size

	// Compute uncorrected next estimated state
	transition(m_x, pControl);

	// Compute uncorrected next estimated covariance of state
	GData* pTemp = GData::multiply(*m_pP, *pA, false, true);
	Holder<GData> hTemp(pTemp);
	delete(m_pP);
	m_pP = GData::multiply(*pA, *pTemp, false, false);
	addTransitionNoise(m_pP);
}

void GExtendedKalmanFilter::correct(const double* pObservation, GData* pH)
{
	// Check values
	GAssert(pH->rows() == (size_t)m_obsDims && pH->cols() == m_stateDims); // observation Jacobian wrong size

	// Compute the Kalman gain
	GData* pK;
	{
		GData* pTemp = GData::multiply(*m_pP, *pH, false, true);
		Holder<GData> hTemp(pTemp);
		GData* pTemp2 = GData::multiply(*pH, *pTemp, false, false);
		Holder<GData> hTemp2(pTemp2);
		addObservationNoise(pTemp2);
		GData* pTemp3 = pTemp2->pseudoInverse();
		Holder<GData> hTemp3(pTemp3);
		pK = GData::multiply(*pTemp, *pTemp3, false, false);
	}
	Holder<GData> hK(pK);

	// Correct the estimated state
	observation(m_z, m_x);
	GVec::multiply(m_z, -1.0, m_obsDims);
	GVec::add(m_z, pObservation, m_obsDims);
	pK->multiply(m_z, m_zz, false);
	GVec::add(m_x, m_zz, m_stateDims);

	// Correct the estimated covariance of state
	{
		GData* pTemp = GData::multiply(*pK, *pH, false, false);
		Holder<GData> hTemp(pTemp);
		pTemp->multiply(-1.0);
		for(int i = 0; i < m_stateDims; i++)
			pTemp->row(i)[i] += 1.0;
		GData* pTemp2 = GData::multiply(*pTemp, *m_pP, false, false);
		delete(m_pP);
		m_pP = pTemp2;
	}
}
