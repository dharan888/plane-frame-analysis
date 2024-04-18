/*********************************************
Implementation of the CElementLoads class.
*********************************************/
#include "ElementLoads.h"

CElementLoads::CElementLoads ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = 0;
    m_Type = ELType::CONCLOADX;
    m_fValue1 = m_fValue2 = 0.0f;
}

CElementLoads::CElementLoads (const CElementLoads& ELO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = ELO.m_nElement;
    m_Type = ELO.m_Type;
    m_fValue1 = ELO.m_fValue1;
    m_fValue2 = ELO.m_fValue2;
}

CElementLoads::~CElementLoads ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementLoads::GetValues (int& nE, ELType& Type,
                               float& f1, float& f2) const
// ---------------------------------------------------------------------------
// Function: gets the element load information
// Input:    element #, load type, 2 load associated values
// Output:   values for all these variables
// ---------------------------------------------------------------------------
{
    nE = m_nElement;
    Type = m_Type;
    f1 = m_fValue1;
    f2 = m_fValue2;
}

void CElementLoads::SetValues (const int nE, const ELType Type,
                               const float f1, const float f2)
// ---------------------------------------------------------------------------
// Function: sets the element load information
// Input:    element #, load type, 2 load associated values
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = nE;
    m_Type = Type;
    m_fValue1 = f1;
    m_fValue2 = f2;
}

void CElementLoads::GetENF (CVector<float>& fVENFL1, CVector<float>& fVENFL2,
                            CVector<float>& fVENFG1, CVector<float>& fVENFG2,
                            const float L, const double l, const double m,
                            bool bLH, bool bRH) const
// ---------------------------------------------------------------------------
// Function: gets the equivalent nodal forces
// Input:    vectors to store ENF at start node and end node, element length
// Output:   the two vectors suitably populated
// ---------------------------------------------------------------------------
{
    fVENFL1.Set(0.0);
    fVENFL2.Set(0.0);
    fVENFG1.Set(0.0);
    fVENFG2.Set(0.0);


    if (m_Type == ELType::DISTLOAD)
    {
        float wL = m_fValue1;
        float wR = m_fValue2;

        fVENFL1(2) =  (L / 20.0f) * (7.0f * wL + 3.0f * wR);
        fVENFL2(2) =  (L / 20.0f) * (3.0f * wL + 7.0f * wR);

        fVENFL1(3) =  (L * L / 60.0f) * (3.0f * wL + 2.0f * wR);
        fVENFL2(3) = -(L * L / 60.0f) * (2.0f * wL + 3.0f * wR);

    }
    else
    {
        float P = m_fValue2;
        float a = m_fValue1;
        float b = L - a;

        if (m_Type == ELType::CONCLOADX)
        {
            fVENFL1(1) = P * b / L;
            fVENFL2(1) = P * a / L;
        }
        else if (m_Type == ELType::CONCLOADY)
        {
            fVENFL1(2) = P * b * b * (L + 2.0f * a) / pow(L, 3.0f);
            fVENFL2(2) = P * a * a * (L + 2.0f * b) / pow(L, 3.0f);

            fVENFL1(3) =  P * a * b * b / (L * L);
            fVENFL2(3) = -P * a * a * b / (L * L);
        }
        else if (m_Type == ELType::CONCMOMENT)
        {
            double M = P;

            fVENFL1(2) = -6.0f * M * a * b / pow(L, 3.0f);
            fVENFL2(2) =  6.0f * M * a * b / pow(L, 3.0f);

            fVENFL1(3) = -M * b * (2.0f * a - b) / (L * L);
            fVENFL2(3) =  M * a * (a - 2.0f * b) / (L * L);
        }

    }

    if (bLH && bRH)
    {
        fVENFL1(2) -= ((fVENFL1(3) + fVENFL2(3)) / L);
        fVENFL1(2) += ((fVENFL1(3) + fVENFL2(3)) / L);
        fVENFL1(3) = 0.0f;
        fVENFL2(3) = 0.0f;
    }
    else if (bLH)
    {
        fVENFL1(2) -= (1.5f * fVENFL1(3) / L);
        fVENFL2(2) += (1.5f * fVENFL1(3) / L);
        fVENFL2(3) -= 0.5f * fVENFL1(3);
        fVENFL1(3) = 0.0f;
    }
    else if (bRH)
    {
        fVENFL1(2) -= (1.5f * fVENFL2(3) / L);
        fVENFL2(2) += (1.5f * fVENFL2(3) / L);
        fVENFL1(3) -= 0.5f * fVENFL2(3);
        fVENFL2(3) = 0.0f;
    }

    float fl = float(l);
    float fm = float(m);

    fVENFG1(1) = fl * fVENFL1(1) - fm * fVENFL1(2);
    fVENFG1(2) = fm * fVENFL1(1) + fl * fVENFL1(2);
    fVENFG1(3) = fVENFL1(3);
    fVENFG2(1) = fl * fVENFL2(1) - fm * fVENFL2(2);
    fVENFG2(2) = fm * fVENFL2(1) + fl * fVENFL2(2);
    fVENFG2(3) = fVENFL2(3);

}

