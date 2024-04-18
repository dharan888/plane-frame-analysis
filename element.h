/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include "xstype.h"

class CElement
{
    public:
        CElement ();                  // default ctor
        CElement(const CElement& EO); // copy ctor
        ~CElement ();                 // dtor

        // accessor functions
        void GetENodes (int& nSN, int& nEN) const;
        int  GetMatPropertyGroup () const;
        CXSType* GetEPropertyGroup () const;
        int  GetEPropertyGroupNo() const;
        float GetLength() const;
        void GetDirectionCosines (double& dl, double& dm) const;
        
        // modifier functions
        void SetENodes (const int nSN, const int nEN);
        void SetEPropertyGroup (const int nEPG, CXSType*);
        void SetMatPropertyGroup (const int);
        void SetLength (const float);
        void SetDirectionCosines (const double dl, const double dm);

    private:
        int m_nSN;		     // start node
        int m_nEN;		     // end node
        float m_fL;          // length of the member
        double m_dl;          // direction cosine along x
        double m_dm;          // direction cosine along y
        int m_nMPGroup;      // material property group
        int m_nEPGroup;          // element group number
        CXSType* m_pEPGroup; // element property group
};
