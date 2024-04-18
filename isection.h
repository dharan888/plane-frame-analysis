/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "arraycontainersEXH.h"
#include "xstype.h"

class CISection: public CXSType
{
    static const int m_numISDimensions = 4;
    public:
        CISection (const CVector<float>& fV);
        ~CISection ();

        // helper functions
        virtual void ComputeProperties ();

    private:
};
