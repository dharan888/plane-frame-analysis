/*********************************************
Implementation of the CElement class.
*********************************************/
#include "element.h"

CElement::CElement ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nSN = m_nEN = m_nMPGroup = m_nEPGroup = 0;
    m_fL = 0.0f;
    m_dl = m_dm = 0.0f;
    m_pEPGroup = NULL;
}

CElement::CElement (const CElement& EO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nSN = EO.m_nSN;
    m_nEN = EO.m_nEN;
    m_nMPGroup = EO.m_nMPGroup;
    m_pEPGroup = EO.m_pEPGroup;
}

CElement::~CElement ()
// ---------------------------------------------------------------------------
// Function: dtor 
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElement::GetENodes (int& nSN, int& nEN) const
// ---------------------------------------------------------------------------
// Function: obtains the element start and end node numbers
// Input:    variables to hold the values
// Output:   the start and end node numbers
// ---------------------------------------------------------------------------
{
    nSN = m_nSN;
    nEN = m_nEN;
}

void CElement::SetENodes (const int nSN, const int nEN) 
// ---------------------------------------------------------------------------
// Function: sets the element start and end node numbers
// Input:    the start and end node numbers
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nSN = nSN;
    m_nEN = nEN;
}

void CElement::SetEPropertyGroup (const int nEPG, CXSType* pEPG)
// ---------------------------------------------------------------------------
// Function: sets the pointer to the element property group associated
//           with the element along with group no
// Input:    element property pointer 
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nEPGroup = nEPG;
    m_pEPGroup = pEPG;
}

CXSType* CElement::GetEPropertyGroup () const
// ---------------------------------------------------------------------------
// Function: gets the pointer to the element property group associated
//           with the element
// Input:    none
// Output:   returns the element property pointer
// ---------------------------------------------------------------------------
{
    return m_pEPGroup;
}

void CElement::SetMatPropertyGroup (const int nMPG)
// ---------------------------------------------------------------------------
// Function: sets the material group group # associated with the element
// Input:    material group #
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nMPGroup = nMPG;
}

int CElement::GetMatPropertyGroup () const
// ---------------------------------------------------------------------------
// Function: gets the material group group # associated with the element
// Input:    none
// Output:   returns the material group #
// ---------------------------------------------------------------------------
{
    return m_nMPGroup;
}

float CElement::GetLength() const
// ---------------------------------------------------------------------------
// Function: obtains the element length
// Input:    none
// Output:   element length
// ---------------------------------------------------------------------------
{
    return m_fL;
}

void CElement::SetLength(const float fL)
// ---------------------------------------------------------------------------
// Function: sets the element length
// Input:    element length
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fL = fL;
}

void CElement::GetDirectionCosines (double& dl, double& dm) const
// ---------------------------------------------------------------------------
// Function: obtains the direction cosines of the element
// Input:    variables to hold direction cosines
// Output:   variables holding direction cosines
// ---------------------------------------------------------------------------
{
    dl = m_dl;
    dm = m_dm;
}

void CElement::SetDirectionCosines (const double dl, const double dm)
// ---------------------------------------------------------------------------
// Function: sets the direction cosines of the element
// Input:    variables holding direction cosines
// Output:   none
// ---------------------------------------------------------------------------
{
    m_dl = dl;
    m_dm = dm;
}

int CElement::GetEPropertyGroupNo() const
// ---------------------------------------------------------------------------
// Function: gets the material group group # associated with the element
// Input:    none
// Output:   returns the material group #
// ---------------------------------------------------------------------------
{
    return m_nEPGroup;
}