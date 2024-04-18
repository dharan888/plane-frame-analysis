/*********************************************
Implementation of Node Class
*********************************************/
#include "Node.h"

CNode::CNode ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVCoor.SetSize (NDIM);   m_fVCoor.Set (0.0f);
    m_VFC.SetSize (DOFPN);     m_VFC.Set (Fixity::FREE);
    m_fVDisp.SetSize (DOFPN);  m_fVDisp.Set (0.0f);
    m_fdT = 0.0f;
}

CNode::CNode (const CNode& NO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVCoor = NO.m_fVCoor;
    m_VFC = NO.m_VFC;   
    m_fVDisp = NO.m_fVDisp;
}

CNode::~CNode ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CNode::GetCoords (CVector<float>& fV) const
// ---------------------------------------------------------------------------
// Function: gets the nodal coordinates at the current node
// Input:    vector to hold the nodal coordinates
// Output:   nodal coordinates
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= NDIM; i++)
        fV(i) = m_fVCoor(i);
}

void CNode::SetCoords (const CVector<float>& fV)
// ---------------------------------------------------------------------------
// Function: sets the nodal coordinates at the current node
// Input:    vector holding the nodal coordinates
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= NDIM; i++)
        m_fVCoor(i) = fV(i);
}

void CNode::GetFixity (CVector<Fixity>& VFC) const
// ---------------------------------------------------------------------------
// Function: gets the nodal fixity conditions at the current node
// Input:    vector to hold the nodal fixities
// Output:   nodal fixities
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        VFC(i) = m_VFC(i);
}

void CNode::SetFixity (const CVector<Fixity>& VFC)
// ---------------------------------------------------------------------------
// Function: sets the nodal fixity conditions at the current node
// Input:    vector holding the nodal fixities
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        m_VFC(i) = VFC(i);
}

std::ostream &operator<< (std::ostream& ofs, const CNode::Fixity& FC)
// ---------------------------------------------------------------------------
// Function: non-class function to display the nodal fixity as a text string
// Input:    output stream, fixity condition
// Output:   none
// ---------------------------------------------------------------------------
{
    if (FC == CNode::Fixity::FREE)
        ofs << "Free";
    else if (FC == CNode::Fixity::SPECIFIED)
        ofs << "Specified";
    else
        ofs << "Hinged";

    return ofs;
}

float CNode::ToDistance(const CNode& Nd) const
// ---------------------------------------------------------------------------
// Function: sets the nodal fixity conditions at the current node
// Input:    vector holding the nodal fixities
// Output:   none
// ---------------------------------------------------------------------------
{
    float fD = 0.0;
    for (int i = 1; i <= NDIM; i++)
    {
        fD += pow(Nd.m_fVCoor(i) - m_fVCoor(i), 2.0f);
    }
    fD = sqrt(fD);

    return fD;

}

void CNode::DirectionCosines(const CNode& Nd, CVector<double>& dVDC) const
// ---------------------------------------------------------------------------
// Function: sets the nodal fixity conditions at the current node
// Input:    vector holding the nodal fixities
// Output:   none
// ---------------------------------------------------------------------------
{
    double dL = double(ToDistance(Nd));
    for (int i = 1; i <= NDIM; i++)
    {
        dVDC(i) = (double(Nd.m_fVCoor(i) - m_fVCoor(i))/dL);
    }

}


float CNode::GetNodalTemperature() const
// ---------------------------------------------------------------------------
// Function: gets the nodal Temperature Change
// Input:    none
// Output:   nodal temperature change
// ---------------------------------------------------------------------------
{
    return m_fdT;
}

void CNode::SetNodalTemperature (const float fdT)
// ---------------------------------------------------------------------------
// Function: sets the nodal Temperature Change
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fdT = fdT;
}