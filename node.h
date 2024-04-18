/*********************************************
Program Planar Frame
*********************************************/
#pragma once
#include <iostream>
#include "arraycontainersEXH.h"
#include "constants.h"

class CNode
{
    public:
        enum class Fixity {FREE, SPECIFIED, HINGED};
        friend std::ostream &operator << (std::ostream&, 
                                         const CNode::Fixity&); //comment: what is happening here?
        CNode ();             // ctor
        CNode (const CNode&); // copy ctor
        ~CNode ();            // dtor

        // accessor functions
        void GetCoords (CVector<float>& fVC) const;
        void GetFixity (CVector<Fixity>& VFC) const;
        float GetNodalTemperature() const;

        // modifier functions
        void SetCoords (const CVector<float>& fVC);
        void SetFixity (const CVector<Fixity>& VFC);
        void SetNodalTemperature(const float fdT);

        // helper functions
        float ToDistance(const CNode& Nd) const;
        void DirectionCosines(const CNode& Nd, CVector<double>& dVDC) const;

    private:
        CVector<float>	m_fVCoor;	// coordinates
        CVector<Fixity>	m_VFC;		// fixity conditions
        CVector<float>	m_fVDisp;	// known displacements
        float m_fdT; //nodal temperature change;
};
