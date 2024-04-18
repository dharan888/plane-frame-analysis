/*********************************************
Implementation of the CElementResponse class.
*********************************************/
#include "ElementResponse.h"

CElementResponse::CElementResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVFStartNode.SetSize (DOFPN);  m_fVFStartNode.Set(0.0f);
    m_fVFEndNode.SetSize (DOFPN);    m_fVFEndNode.Set(0.0f);
    m_fMaxCS = m_fMaxTS = m_fMaxSS = 0.0f;
}

CElementResponse::~CElementResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementResponse::GetForces (CVector<float>& fVFSN,
                                  CVector<float>& fVFEN) const
// ---------------------------------------------------------------------------
// Function: gets the element force response
// Input:    vectors to hold element forces at start and end nodes
// Output:   the element force values
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
    {
        fVFSN(i) = m_fVFStartNode(i);
        fVFEN(i) = m_fVFEndNode(i);
    }
}

void CElementResponse::SetForces(const CVector<float>& fVFSN,
                                  const CVector<float>& fVFEN)
// ---------------------------------------------------------------------------
// Function: gets the element force response
// Input:    vectors to hold element forces at start and end nodes
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
    {
        m_fVFStartNode(i) = fVFSN(i);
        m_fVFEndNode (i) = fVFEN(i);
    }
}

void CElementResponse::GetMaxStresses (float& fMaxCS,float& fMaxTS, float& fMaxSS) const
    // ---------------------------------------------------------------------------
    // Function: gets the maximum element compressive, tensile and shear stresses
    // Input:    variables to hold compressive, tensile and shear stresses
    // Output:   variables holding compressive(fMaxCS), tensile(fMaxTS) and shear stresses (fMaxSS)
    // ---------------------------------------------------------------------------
{
    fMaxCS = m_fMaxCS;
    fMaxTS = m_fMaxTS;
    fMaxSS = m_fMaxSS;
}

void CElementResponse::SetMaxStresses (const float fMaxCS, const float fMaxTS, const float fMaxSS)
// ---------------------------------------------------------------------------
// Function: Sets the maximum element compressive, tensile and shear stresses
// Input:    variables holding compressive(fMaxCS), tensile(fMaxTS) and shear stresses(fMaxSS)
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fMaxCS = fMaxCS;
    m_fMaxTS = fMaxTS;
    m_fMaxSS = fMaxSS;
}

void CElementResponse::GetMaxForces (float& fMaxP, float& fMaxV, float& fMaxM) const
// ---------------------------------------------------------------------------
// Function: Gets the maximum element axial, shear and moment
// Input:    variables to hold  max axial(fMaxP), shear(fMaxV) and moment(fMaxM)
// Output:   variables holding  max axial(fMaxP), shear(fMaxV) and moment(fMaxM)
// ---------------------------------------------------------------------------
{
    fMaxP = m_fMaxP;
    fMaxV = m_fMaxV;
    fMaxM = m_fMaxM;
}


void CElementResponse::SetMaxForces (const float fMaxP, const float fMaxV, const float fMaxM)
// ---------------------------------------------------------------------------
// Function: Sets the maximum element axial, shear and moment
// Input:    variables holding  max axial(fMaxP), shear(fMaxV) and moment(fMaxM)
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fMaxP = fMaxP;
    m_fMaxV = fMaxV;
    m_fMaxM = fMaxM;
}