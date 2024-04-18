/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "xstype.h"

class CCircSolid: public CXSType
{
    static const int m_numCircDimensions = 1;
    public:
        CCircSolid(const CVector<float>& fV);
        CCircSolid(const CCircSolid&);
        ~CCircSolid();

        // helper functions
        virtual void ComputeProperties ();

    private:
};
