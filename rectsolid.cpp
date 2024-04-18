/*********************************************
Derived class from CXSType base class
to store rectangular solid data and properties
*********************************************/
#include <cmath>
#include <iostream>
#include "rectsolid.h"

CRectSolid::CRectSolid (const CVector<float>& fV) 
                     : CXSType (m_numRectDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with rectangular solid dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_numRectDimensions; i++)
        m_fVDimensions(i) = fV(i);
    ComputeProperties ();
}

CRectSolid::~CRectSolid ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CRectSolid::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the rectangular solid properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // height
    float fH = m_fVDimensions(1); 
    // width
    float fW = m_fVDimensions(2); 

    // cross-sectional area
    m_fArea = fH * fW;
    // MOI y-axis
    m_fIyy = fH*pow(fW, 3.0f)/12.0f;
    // MOI z-axis
    m_fIzz = fW*pow(fH, 3.0f)/12.0f;
    // Section Modulus y-axis
    m_fSyy = m_fIyy / (0.50f*fW);
    // Section Modulus z-axis
    m_fSzz = m_fIzz / (0.50f*fH);
    // Shear Factor y-axis
    m_fSFy = (2.0f / 3.0f) * m_fArea;
    // Shear Factor z-axis
    m_fSFz = (2.0f / 3.0f) * m_fArea;
}