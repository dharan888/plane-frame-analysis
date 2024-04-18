/*********************************************
Derived class from CXSType base class
to store rectangular solid data and properties
*********************************************/
#include <cmath>
#include <iostream>
#include "circsolid.h"

CCircSolid::CCircSolid (const CVector<float>& fV) 
                     : CXSType (m_numCircDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with circular solid dimension
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_numCircDimensions; i++)
        m_fVDimensions(i) = fV(i);
        ComputeProperties ();
}

CCircSolid::~CCircSolid()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CCircSolid::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the circular solid properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    const float PI = 3.141593f;

    // radius
    float fR = m_fVDimensions(1); 

    // cross-sectional area
    m_fArea = PI * pow(fR, 2.0f);
    // MOI y-axis
    m_fIyy  = PI * pow(fR, 4.0f) / 4.0f;
    // MOI z-axis
    m_fIzz  = m_fIyy;
    // Section Modulus
    m_fSzz = m_fSyy = m_fIzz / fR;
    // Shear Factor y-axis
    m_fSFy = (3.0f / 4.0f) * m_fArea;
    // Shear Factor z-axis
    m_fSFz = (3.0f / 4.0f) * m_fArea;
}