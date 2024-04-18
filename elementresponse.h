/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "constants.h"

class CElementResponse
{
    public:
        CElementResponse ();    // default ctor
        ~CElementResponse ();   // dtor

        // accessor functions
        void GetForces (CVector<float>& fVFStartNode,
                        CVector<float>& fVFEndNode) const;
        void GetMaxStresses (float& fMaxCS, float& fMaxTS, 
                             float& fMaxSS) const;
        void GetMaxForces (float& fMaxP, float& fMaxV,
                           float& fMaxM) const;

        // modifier functions
        void SetForces (const CVector<float>& fVFStartNode,
                        const CVector<float>& fVFEndNode);
        void SetMaxStresses (const float fMaxCS, const float fMaxTS,
                            const float fMaxSS);
        void SetMaxForces (const float fMaxP,const float fMaxV,
                           const float fMaxM);

    private:
        float m_fMaxCS, m_fMaxTS, m_fMaxSS;
        float m_fMaxP, m_fMaxV, m_fMaxM;
        CVector<float> m_fVFStartNode;
                // nodal forces at start node
        CVector<float> m_fVFEndNode;
                // nodal forces at end node
};
