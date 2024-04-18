/*********************************************
Implementation of the CNodalResponse class.
*********************************************/
#include "nodalresponse.h"

CNodalResponse::CNodalResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVD.SetSize (DOFPN, 0.0f);
    m_fVR.SetSize (DOFPN, 0.0f);
}


CNodalResponse::~CNodalResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}


void CNodalResponse::GetDisplacements (CVector<float>& fVD) const
// ---------------------------------------------------------------------------
// Function: gets the nodal displacements at the current node
// Input:    vector to hold the nodal displacement values
// Output:   nodal displacement values
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        fVD(i) = m_fVD(i);
}


void CNodalResponse::GetReactions (CVector<float>& fVR) const
// ---------------------------------------------------------------------------
// Function: gets the nodal reactions at the current node
// Input:    vector to hold the nodal reactions values
// Output:   nodal reactions values
// ---------------------------------------------------------------------------
{
    for (int i = 1; i <= DOFPN; i++)
        fVR(i) = m_fVR(i);
}

void CNodalResponse::SetDisplacements (const CVector<float>& fVD)
// ---------------------------------------------------------------------------
// Function: sets the nodal displacements at the current node
// Input:    vector holding the nodal displacement values
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        m_fVD(i) = fVD(i);
}

void CNodalResponse::SetReactions (const CVector<float>& fVR)
// ---------------------------------------------------------------------------
// Function: sets the nodal reactions at the current node
// Input:    vector holding the nodal reactions values
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i = 1; i <= DOFPN; i++)
        m_fVR(i) = fVR(i);
}