/*********************************************
Implementation of the CMaterial class.
*********************************************/
#include "material.h"

CMaterial::CMaterial ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fYM = 0.0f;
    m_fCTE = 0.0f;
}

CMaterial::~CMaterial ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

float CMaterial::GetYM () const
// ---------------------------------------------------------------------------
// Function: gets the young's modulus value
// Input:    none
// Output:   returns the YM value
// ---------------------------------------------------------------------------
{
    return m_fYM;
}

float CMaterial::GetCTE () const
// ---------------------------------------------------------------------------
// Function: gets the coef of thermal expansion value
// Input:    none
// Output:   returns the CTE value
// ---------------------------------------------------------------------------
{
    return m_fCTE;
}

void CMaterial::SetYM (const float fYM)
// ---------------------------------------------------------------------------
// Function: sets the young's modulus value
// Input:    YM value
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fYM = fYM;
}

void CMaterial::SetCTE (const float fCTE)
// ---------------------------------------------------------------------------
// Function: sets the coef of thermal expansion value
// Input:    CTE value
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fCTE = fCTE;
}

