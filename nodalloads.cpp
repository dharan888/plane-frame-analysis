/*********************************************
Implementation of the CNodalLoads class.
*********************************************/
#include "NodalLoads.h"

CNodalLoads::CNodalLoads ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nNode = 0;
    m_fVLoads.SetSize (DOFPN); m_fVLoads.Set(0.0f);
}

CNodalLoads::CNodalLoads (const CNodalLoads& NLO)
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVLoads.SetSize (DOFPN); 
    m_fVLoads = NLO.m_fVLoads;
}

CNodalLoads::~CNodalLoads ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CNodalLoads::GetLoadValues (int& nN, CVector<float>& fV) const
// ---------------------------------------------------------------------------
// Function: gets the nodal loads at the current node
// Input:    vector to hold the nodal load values and the integer to node no
// Output:   nodal load values
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        fV(i) = m_fVLoads(i);
    nN = m_nNode;
}

void CNodalLoads::SetValues (const int nN, const CVector<float>& fV)
// ---------------------------------------------------------------------------
// Function: sets the nodal loads at the current node and node no
// Input:    vector with the nodal load values
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        m_fVLoads(i) = fV(i);
    m_nNode = nN;
}
