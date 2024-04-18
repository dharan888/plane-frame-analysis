/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "constants.h"

class CNodalLoads
{
    public:
        CNodalLoads ();                   // default ctor 
        CNodalLoads (const CNodalLoads&); //copy ctor
        ~CNodalLoads ();

        // accessor functions
        void GetLoadValues (int& nN, CVector<float>&) const;

        // modifier functions
        void SetValues (const int nN, const CVector<float>&);

    private:
        CVector<float> m_fVLoads;	// nodal loads
        int m_nNode; //node number
};
