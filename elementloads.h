/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "constants.h"

class CElementLoads
{
    public:
        enum class ELType {CONCLOADX, CONCLOADY, CONCMOMENT, DISTLOAD};

        CElementLoads ();                         // default ctor
        CElementLoads (const CElementLoads& ELO); // copy ctor
        ~CElementLoads ();                        // dtor

        // accessor functions
        void GetValues (int& nE, ELType&, 
                        float& fLeft, float& fRight) const;

        void CElementLoads::GetENF(CVector<float>& fVENFL1, CVector<float>& fVENFL2,
                                   CVector<float>& fVENFG1, CVector<float>& fVENFG2,
                                   const float L, const double l, const double m,
                                   bool bLHinged, bool bRHinged) const;

        // modifier functions
        void SetValues (const int nE, const ELType Type,
                        const float fLeft, const float fRight);
    private:
        int    m_nElement;	// element number
        ELType m_Type;	    // load type: concentrated or distributed
        float  m_fValue1;	// distance from start node
                            // or load intensity at start node
        float  m_fValue2;	// load value or
                            // load intensity at end node
};		
