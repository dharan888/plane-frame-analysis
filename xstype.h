/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include <string>
#include "arraycontainersEXH.h"

class CXSType
{
    public:
        CXSType ();
        CXSType (int);
        virtual ~CXSType ();

        // helper function
        //void DisplayProperties ();

        // accessor functions
        void GetProperties (float& fA, float& fIyy, float& fIzz);
        void GetDimensions (CVector<float>&) const;
        void GetSectionModuli (float& fSyy, float& fSzz);
        void GetShearFactors (float& fSFy, float& fSFz);

        virtual void ComputeProperties () = 0;

    private:
        void Initialize ();

    protected:
        std::string m_strID;           // identification tag
        float m_fArea;                 // x/s area
        float m_fIyy;                  // MOI y-axis
        float m_fIzz;                  // MOI z-axis
        float m_fSyy;                  // SM y-axis
        float m_fSzz;                  // SM z-axis
        float m_fSFy;                  // shear factor Y
        float m_fSFz;                  // shear factor Z
        int   m_numDimensions;         // number of dimensions
        CVector<float> m_fVDimensions; // the dimensions
};
