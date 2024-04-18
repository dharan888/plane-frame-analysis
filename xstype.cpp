/*********************************************
Base class for element cross-sectional properties
*********************************************/
#include "xstype.h"

CXSType::CXSType ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    Initialize ();
}

CXSType::CXSType (int numDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    # of cross-sectional dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    Initialize ();
    m_numDimensions = numDimensions;
    m_fVDimensions.SetSize (m_numDimensions);
}

void CXSType::Initialize ()
// ---------------------------------------------------------------------------
// Function: initializes all the member variables with default values
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_numDimensions = 0;
    m_fArea = 0.0f;
    m_fIyy = 0.0f;
    m_fIzz = 0.0f;
    m_fSyy = 0.0f;
    m_fSzz = 0.0f;
    m_fSFy = 0.0f;
    m_fSFz = 0.0f;
}

CXSType::~CXSType ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CXSType::GetProperties (float& fArea, float& fIyy, float& fIzz)
// ---------------------------------------------------------------------------
// Function: gets the computed cross-sectional properties
// Input:    variables to hold area, Iyy, Izz
// Output:   area, Iyy, Izz values
// ---------------------------------------------------------------------------
{
    ComputeProperties ();
    fArea = m_fArea;
    fIyy = m_fIyy;
    fIzz = m_fIzz;
}

void CXSType::GetSectionModuli (float& fSyy, float& fSzz)
// ---------------------------------------------------------------------------
// Function: gets the computed cross-sectional properties
// Input:    variables to hold area, Iyy, Izz
// Output:   area, Iyy, Izz values
// ---------------------------------------------------------------------------
{
    ComputeProperties ();
    fSyy = m_fSyy;
    fSzz = m_fSzz;
}

void CXSType::GetShearFactors (float& fSFy, float& fSFz)
// ---------------------------------------------------------------------------
// Function: gets the computed cross-sectional properties
// Input:    variables to hold area, Iyy, Izz
// Output:   area, Iyy, Izz values
// ---------------------------------------------------------------------------
{
    ComputeProperties ();
    fSFy = m_fSFy;
    fSFz = m_fSFz;
}

void CXSType::GetDimensions (CVector<float>& fV) const
// ---------------------------------------------------------------------------
// Function: gets the cross-sectional dimensions
// Input:    vector to hold x/s dimensions
// Output:   x/s dimensions
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_numDimensions; i++)
        fV(i) = m_fVDimensions(i);
}

