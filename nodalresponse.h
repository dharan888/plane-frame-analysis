/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "constants.h"

class CNodalResponse
{
    public:
        CNodalResponse ();   // ctor
        ~CNodalResponse ();  // dtor

        // accessor functions
        void GetDisplacements (CVector<float>& fVD) const;

        void GetReactions (CVector<float>& fVR) const;

        // modifier functions
        void SetDisplacements (const CVector<float>& fVD);

        void SetReactions (const CVector<float>& fVR);

    private:
        CVector<float> m_fVD; // nodal displacements	
        CVector<float> m_fVR; // nodal reactions
};
