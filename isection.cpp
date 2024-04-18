/*********************************************
Derived class from CXSType base class
to store I-section data and properties
*********************************************/
#include <cmath>
#include <iostream>
#include "isection.h"

CISection::CISection (const CVector<float>& fV) : 
                      CXSType (m_numISDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with I-section dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_numISDimensions; i++)
        m_fVDimensions(i) = fV(i);
    ComputeProperties ();
}

CISection::~CISection ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CISection::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the I-section properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // web height
    float wH = m_fVDimensions(1); 
    // web thickness
    float wT = m_fVDimensions(2); 
    // flange width
    float fW = m_fVDimensions(3); 
    // flange thickness
    float fT = m_fVDimensions(4); 

    // cross-sectional area
    m_fArea = wH*wT + 2.0f*fW*fT;
    // MOI z-axis
    float fOI, fII;
    fOI = (pow(wH+2.0f*fT, 3.0f) * fW)/12.0f;
    fII = 2.0f*(pow(wH, 3.0f) * (0.5f*(fW-wT)))/12.0f;
    m_fIzz = fOI - fII;
    // MOI y-axis
    float fIyy = 2.0f*(pow(fW,3.0f)*fT)/12.0f + 
                 pow(wT,3.0f)*wH/12.0f;
    m_fIyy = fIyy;

    //Section Modulus 
    m_fSyy = m_fIyy / (0.50f*fW);
    m_fSzz = m_fIzz / (0.50f*wH + fT);

    //Shear factors
    float fQy = (fW * fW * fT / 4.0f) + (wH * wT * wT / 4.0f);
    float fQz = fW * fT * (wH + fT) * 0.50f + (wH * wT * 0.5f) * (wH * 0.25f);

    m_fSFy = wT * m_fIzz / fQz;
    m_fSFz = fT * m_fIyy / fQy;
}