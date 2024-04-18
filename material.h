/*********************************************
Program Planar Frame
*********************************************/
#pragma once

class CMaterial
{
    public:
        CMaterial ();   // ctor
        ~CMaterial ();  // dtor

        // accessor functions
        float GetYM () const;
        float GetCTE () const;

        // modifier functions
        void SetYM (const float);
        void SetCTE (const float);

    private:
        float m_fYM;   // young's modulus
        float m_fCTE;  // coef of thermal expansion
};
