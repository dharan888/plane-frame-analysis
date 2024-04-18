/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "xstype.h"

class CRectSolid: public CXSType
{
    static const int m_numRectDimensions = 2;
    public:
        CRectSolid (const CVector<float>& fV);
        CRectSolid (const CRectSolid&);
        ~CRectSolid ();

        // helper functions
        virtual void ComputeProperties ();

    private:
};
