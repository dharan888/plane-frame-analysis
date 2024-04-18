/*********************************************
Implementation of CFrame Class
*********************************************/
#include "frame.h"
#include "MatToolBox.h"
#include <algorithm>

/* ==================================================================
   ======================= CFrame class =============================
   ================================================================== */

CFrame::CFrame ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nNodes = m_nElements = m_nEPGroups = m_nMatGroups
             = m_nDOF = m_nSDOF = m_nEDOF = m_nLineNumber 
              = m_nElementLoads = m_nNodalLoads = 0;
    m_nDebugLevel = 0;
    m_strDelimiters = "\t, "; // tab space and comma delimited file
    m_strComment    = "**";
}

CFrame::~CFrame ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate ... deallocate ... deallocate
    for (int i=1; i <= m_nEPGroups; i++)
    {
        delete m_EPData(i); 
    }
}

void CFrame::Analyze ()
// ---------------------------------------------------------------------------
// Function: Reads the input data and analyzes the frame
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	// read the problem size
	ReadProblemSize ();

	// set problem size
	SetSize ();

	// read nodal and element data
	ReadFrameModel ();

	// construct system equations
	ConstructSystemKandF ();

	// solve for the nodal displacements
	Solve ();

	// compute secondary unknowns
	Response ();

    // extract support reactions and find error norms
    ObtainSupportReactions();

	// create output file
	CreateOutput ();
}

void CFrame::SetSize ()
// ---------------------------------------------------------------------------
// Function: memory allocation for all major arrays in the program
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // allocate space for nodal data
    m_NodalData.SetSize (m_nNodes);
    // allocate space for nodal loads data
    m_MaterialData.SetSize (m_nMatGroups);
    // allocate space for x/s data
    m_EPData.SetSize (m_nEPGroups);
    // allocate space for nodal response data
    m_NodalResponseData.SetSize (m_nNodes);
    // allocate space for element data
    m_ElementData.SetSize (m_nElements);
    // allocate space for element response data
    m_ElementResponseData.SetSize (m_nElements);
    // allocate space for element loads, if required
    if (m_nElementLoads > 0)
        m_ElementLoadData.SetSize (m_nElementLoads);
    m_ELL.SetSize (DOFPE, m_nElements, 0.0f); 
    m_ELG.SetSize (DOFPE, m_nElements, 0.0f);
    // allocate space for nodal loads, if required
    if (m_nNodalLoads > 0)
        m_NodalLoadData.SetSize (m_nNodes);
    // allocate and initialize major matrices and vectors
    m_nDOF = DOFPN*m_nNodes;
    m_SADOF.SetSize(m_nDOF, 0);
    m_SCDOF.SetSize(m_nDOF, false);
    m_fVR.SetSize(m_nDOF, 0.0f);
}


void CFrame::SetSystemSize()
// ---------------------------------------------------------------------------
// Function: memory allocation for all major arrays in the program
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // allocate and initialize major matrices and vectors
    m_SSM.SetSize(m_nEDOF, m_nEDOF,0.0);
    m_SNF.SetSize(m_nEDOF, 0.0);
    m_SND.SetSize(m_nEDOF, 0.0);

}



void CFrame::ConstructSystemKandF()
// ---------------------------------------------------------------------------
// Function: constructs the system stiffness matrix and force vector
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int nSN, nEN, nMatGrp;
    CXSType* pXSGrp;
    float fE, fA, fL, fCTE, fI, fIyy, fIzz;
    double dl, dm;
    CVector<float> fVC1(NDIM), fVC2(NDIM);
    CVector<CNode::Fixity> VFC1(DOFPN), VFC2(DOFPN);
    CMatrix<double> dMKe(DOFPE, DOFPE);
    CMatToolBox<double> MTB;

    // loop throught the elements
    for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
    {

        // -----------------------------------------------------------------
        //                  GETTING THE ELEMENT PROPERTIES
        // -----------------------------------------------------------------
       
        //Getting the numbers of the element nodes in nSN and nEN
        m_ElementData(ElNo).GetENodes(nSN, nEN);

        //Getting the fixities of the start and end node in VFC1 and VFC2 resp.
        m_NodalData(nSN).GetFixity(VFC1);
        m_NodalData(nEN).GetFixity(VFC2);

        //obtaining the pointer of the element c/s group
        pXSGrp = m_ElementData(ElNo).GetEPropertyGroup();
        pXSGrp->GetProperties(fA, fIyy, fIzz);

        //obtaining the members E and CTE based on  material group
        nMatGrp = m_ElementData(ElNo).GetMatPropertyGroup();
        fE = m_MaterialData(nMatGrp).GetYM();
        fCTE = m_MaterialData(nMatGrp).GetCTE();

        //finding the element length
        fL = m_ElementData(ElNo).GetLength();

        //Finding MI for the bending
        fI = fIzz;
        if (fIyy > fIzz) fI = fIyy;

        //finding the direction cosines
        m_ElementData(ElNo).GetDirectionCosines(dl, dm);

        //finding the boolean value indicating the left and right node is hinged or not
        bool bLH = (VFC1(3) == CNode::Fixity::HINGED ? true : false);
        bool bRH = (VFC2(3) == CNode::Fixity::HINGED ? true : false);


        // -----------------------------------------------------------------
        //       ASSEMBLING THE Effective K matrix and F vector
        // -----------------------------------------------------------------

        //finding the element stiffness matrix
        ElementStiffness(fE, fA, fI, fL, dl, dm, bLH, bRH, dMKe);
        
        CVector<float> fVD1(DOFPN), fVD2(DOFPN), fVD(DOFPN);
        m_NodalResponseData(nSN).GetDisplacements(fVD1);
        m_NodalResponseData(nEN).GetDisplacements(fVD2);

        int nNi, nNj, nAi, nAj, DOFi, DOFj;
        float fD;

        //loop through the element degrees of freedom
        for (int i = 1; i <= DOFPE; i++)
        {
            // find the degrees of freedom
            nNi     = (i <= 3 ? nSN - 1 : nEN - 2);
            DOFi    = nNi * DOFPN + i;
            nAi     = m_SADOF(DOFi); // find the active dof no

            if (nAi > 0) // if the node active
            {
                // loop through the elements
                for (int j = 1; j <= DOFPE; j++)
                {
                    nNj     = (j <= 3 ? nSN - 1 : nEN - 2);
                    DOFj    = nNj * DOFPN + j;
                    nAj     = m_SADOF(DOFj);

                    // if the dof active
                    if (nAj > 0)
                    {
                        m_SSM(nAi, nAj) += dMKe(i, j); //store the structural stiffness matrix
                    }
                    else
                    {
                        fD = (j <= 3 ? fVD1(j) : fVD2(j-DOFPN));
                        if (abs(fD) > 0.0f)
                        {
                            m_SNF(nAi) -= double(fD) * dMKe(i, j);  //store the structural nodal forces
                        }
                    }
                }
                m_SNF(nAi) += double (m_ELG(i, ElNo)); //store the structural nodal forces
            }

        }

        std::string strElem = std::to_string(ElNo);
        if (m_nDebugLevel == 1)
        {
            MTB.PrintMatrixRowWise(dMKe, "Element Stiffness " + strElem, m_FileOutput);
        }

    }


    // debug?
    if (m_nDebugLevel == 1)
    {
        MTB.PrintMatrixRowWise(m_SSM, "Structural Stiffness (After BCs)", m_FileOutput);
        MTB.PrintVector(m_SNF, "Structural Nodal Forces (After BCs)", m_FileOutput);
    }

}


void CFrame::ElementStiffness(const float fE, const float fA, const float fI, 
                              const float fL, const double l, const double m, 
                              const bool bLHinged, const bool bRHinged,
                              CMatrix<double>& Ke)
    // ---------------------------------------------------------------------------
    // Function: Obtain the element stiffness matrix local and global coordinates
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
{
    
    //initiating the Ke matrix
    Ke.Set(0.0);

    //converting the float to double upfront for clarity
    double L = double(fL);
    double E = double(fE);
    double A = double(fA);
    double I = double(fI);

    // all the variables used in teh calculation here are double 
    // so prefix d that normally indicates double variable 
    // is removed to more clearly relate it with actual expressions

    // axial stiffness term
    double ka = A * E / L;

    // shear, moment stiffness terms
    double ks = 0.0, km = 0.0, kms = 0.0;

    if (!bLHinged && !bRHinged)
    {
        //Stiffness Terms
        ks = 12.0 * E * I / pow(L, 3.0);
        km = 4.0 * E * I / L;
        kms = 6.0 * E * I / pow(L, 2.0);
    }
    else if (bLHinged || bRHinged)
    {
        //Stiffness Terms
        ks = 3.0 * E * I / pow(L, 3.0);
        km = 3.0 * E * I / L;
        kms = 3.0 * E * I / pow(L, 2.0);
    }
    //building the element K matrix - global coords
    Ke(1, 1) = ka * l * l + ks * m * m;
    Ke(1, 2) = Ke(2, 1) = (ka - ks) * l * m;
    Ke(1, 4) = Ke(4, 1) = -Ke(1, 1);
    Ke(1, 5) = Ke(5, 1) = -Ke(1, 2);
    Ke(2, 2) = ka * m * m + ks * l * l;
    Ke(2, 4) = Ke(4, 2) = -Ke(2, 1);
    Ke(2, 5) = Ke(5, 2) = -Ke(2, 2);
    Ke(4, 4) = Ke(1, 1);
    Ke(4, 5) = Ke(5, 4) = Ke(1, 2);
    Ke(5, 5) = Ke(2, 2);

    if (!bLHinged)
    {
        Ke(1, 3) = Ke(3, 1) = -kms * m;
        Ke(2, 3) = Ke(3, 2) =  kms *l;
        Ke(3, 3) = km;
        Ke(3, 4) = Ke(4, 3) = -Ke(3, 1);
        Ke(3, 5) = Ke(5, 3) = -Ke(3, 2);
        Ke(3, 6) = Ke(6, 3) =  km / 2.0;
    }

    if (!bRHinged)
    {
        Ke(1, 6) = Ke(6, 1) = -kms * m;
        Ke(2, 6) = Ke(6, 2) =  kms * l;
        Ke(4, 6) = Ke(6, 4) = -Ke(1, 6);
        Ke(5, 6) = Ke(6, 5) = -Ke(2, 6);
        Ke(6, 6) = km;
    }
}


void CFrame::Solve ()
// ---------------------------------------------------------------------------
// Function: solves for the nodal displacements
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CMatToolBox<double> MTB;

    double TOL = 1.0e-6;

    m_SND.SetSize(m_nEDOF, 1);

    MTB.LDLTFactorization(m_SSM, TOL);
    MTB.LDLTSolve(m_SSM, m_SND, m_SNF);

    if (m_nDebugLevel == 1)
    {
        MTB.PrintVector(m_SND, "Structural Nodal Displacements (After BCs)", m_FileOutput);
    }

    CVector<float> fVD(DOFPN);
    
    int DOF, ADOF;
    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalResponseData(i).GetDisplacements(fVD);

        for (int j = 1; j <= DOFPN; j++)
        {
            DOF = (i - 1) * DOFPN + j;
            ADOF = m_SADOF(DOF);
            if (ADOF > 0)
                fVD(j) = float(m_SND(ADOF));
        }

        m_NodalResponseData(i).SetDisplacements(fVD);

    }
}


void CFrame::Response()
// ---------------------------------------------------------------------------
// Function: computes the element nodal forces and support reactions
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int nSN, nEN, nMatGrp;
    CXSType* pXSGrp;
    float fE, fA, fL, fCTE;
    float fIyy, fIzz, fI, fSyy, fSzz, fS;
    float fSFy, fSFz, fSF;

    CVector<float> fVD1(DOFPN), fVD2(DOFPN);
    CVector<float> fVDL1(DOFPN), fVDL2(DOFPN);
    CVector<CNode::Fixity> VFC1(DOFPN), VFC2(DOFPN);
    CVector<float> fVEF1(DOFPN), fVEF2(DOFPN);
    CVector<double> fVEFG1(DOFPN), fVEFG2(DOFPN);


    for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
    {
        // -----------------------------------------------------------------
        //                  GETTING THE ELEMENT PROPERTIES
        // -----------------------------------------------------------------

        //Getting the numbers of the element nodes in nSN and nEN
        m_ElementData(ElNo).GetENodes(nSN, nEN);

        //obtaining the pointer of the element c/s group
        pXSGrp = m_ElementData(ElNo).GetEPropertyGroup();
        pXSGrp->GetProperties(fA, fIyy, fIzz);

        //obtaining the members E and CTE based on  material group
        nMatGrp = m_ElementData(ElNo).GetMatPropertyGroup();
        fE = m_MaterialData(nMatGrp).GetYM();
        fCTE = m_MaterialData(nMatGrp).GetCTE();

        //finding the element length
        fL = m_ElementData(ElNo).GetLength();

        //Finding MI for the bending
        fI = fIzz;
        if (fIyy > fIzz) fI = fIyy;

        //finding the direction cosines
        double dl, dm;
        m_ElementData(ElNo).GetDirectionCosines(dl, dm);

        m_NodalResponseData(nSN).GetDisplacements(fVD1);
        m_NodalResponseData(nEN).GetDisplacements(fVD2);
        float fl = float(dl), fm = float(dm);


        m_NodalData(nSN).GetFixity(VFC1);
        m_NodalData(nEN).GetFixity(VFC2);

        //finding the boolean value indicating the left and right node is hinged or not
        bool bLH = (VFC1(3) == CNode::Fixity::HINGED ? true : false);
        bool bRH = (VFC2(3) == CNode::Fixity::HINGED ? true : false);

        // element displacements in local direction
        fVDL1(1) = fl * fVD1(1) + fm * fVD1(2);
        fVDL1(2) = -fm * fVD1(1) + fl * fVD1(2);
        fVDL2(1) = fl * fVD2(1) + fm * fVD2(2);
        fVDL2(2) = -fm * fVD2(1) + fl * fVD2(2);
        
        // modify the element forces for the hinge
        if (!bLH && !bRH)
        {
            fVDL1(3) = fVD1(3);
            fVDL2(3) = fVD2(3);
        }
        else if (bLH && !bRH)
        {
            fVDL2(3) = fVD2(3);
            fVDL1(3) = (1.5f/fL)*(-fVDL1(2) + fVDL2(2)) - 0.50f*fVDL2(3)
                     + (fL/(4.0f*fE*fI))* m_ELL(3, ElNo);
        }
        else if (bRH && !bLH)
        {
            fVDL1(3) = fVD1(3);
            fVDL2(3) = (1.5f / fL) * (-fVDL1(2) + fVDL2(2)) - 0.50f * fVDL1(3)
                     + (fL / (4.0f * fE * fI)) * m_ELL(6, ElNo);
        }
        else
        {
            fVDL1(3) = (1 / fL) * (-fVDL1(2) + fVDL2(2)) + (fL / (6.0f * fE * fI))
                       * (2 * m_ELL(3, ElNo) - m_ELL(6, ElNo));
            fVDL2(3) = (1 / fL) * (-fVDL1(2) + fVDL2(2)) + (fL / (6.0f * fE * fI))
                       * (2 * m_ELL(6, ElNo) - m_ELL(3, ElNo)); 
        }

        // stiffness terms
        float ka = fA * fE / fL;         
        float ks = 12.0f * fE * fI / pow(fL, 3.0f);
        float km = 4.0f * fE * fI / fL;
        float kms = 6.0f * fE * fI / pow(fL, 2.0f);


        // element forces
        fVEF1(1) = ka * (fVDL1(1) - fVDL2(1));
        fVEF1(2) = ks * (fVDL1(2) - fVDL2(2)) + kms * (fVDL1(3) + fVDL2(3));
        fVEF1(3) = kms * (fVDL1(2) - fVDL2(2)) + km * (fVDL1(3) + 0.5f * fVDL2(3));
        fVEF2(1) = -fVEF1(1);
        fVEF2(2) = -fVEF1(2);
        fVEF2(3) = kms * (fVDL1(2) - fVDL2(2)) + km * (0.50f * fVDL1(3) + fVDL2(3));

        for (int i = 1; i <= DOFPN; i++)
            fVEF1(i) -= m_ELL(i, ElNo);
        for (int i = 1; i <= DOFPN; i++)
            fVEF2(i) -= m_ELL(i + DOFPN, ElNo);

        m_ElementResponseData(ElNo).SetForces(fVEF1, fVEF2);

        //finding the element forces in global direction
        fVEFG1(1) = fl * fVEF1(1) - fm * fVEF1(2);
        fVEFG1(2) = fm * fVEF1(1) + fl * fVEF1(2);
        fVEFG1(3) = fVEF1(3);
        fVEFG2(1) = fl * fVEF2(1) - fm * fVEF2(2);
        fVEFG2(2) = fm * fVEF2(1) + fl * fVEF2(2);
        fVEFG2(3) = fVEF2(3);

        // find the residual vector
        int DOF;
        for (int i = 1; i <= DOFPN; i++)
        {
            DOF = (nSN - 1) * DOFPN + i;
            m_fVR(DOF) += fVEFG1(i);
        }
        for (int i = 1; i <= DOFPN; i++)
        {
            DOF = (nEN - 1) * DOFPN + i;
            m_fVR(DOF) += fVEFG2(i);
        }

        pXSGrp->GetSectionModuli(fSyy, fSzz);
        pXSGrp->GetShearFactors(fSFy, fSFz);

        //Finding MI for the bending
        fI = fIzz; fS = fSzz; fSF = fSFy;
        if (fIyy > fIzz)
        {
            fI = fIyy;
            fS = fSyy;
            fSF = fSFz;
        }

        int nSections = 50;
        //Finding the peak stresses
        ComputePeakStresses(nSections, ElNo, fL, fA, fS, fSF);

    }


}

void CFrame::ComputePeakStresses (const int nSections, const int ElNo, 
                                  const float fL, const float fA, 
                                  const float fS, const float fSF)
// ---------------------------------------------------------------------------
// Function: computes the peak stress values in a member
// Input:    nSections - number of sections to analyze, ElNo - Element No, fL - length
//           fA - Area, fI - MOI, fS - Section Modulus, fSF - shear factor
// Output:   none
// ---------------------------------------------------------------------------
{
    float fP, fV, fM;
    float fCS, fTS, fSS;
    float fMaxCS, fMaxTS, fMaxSS;
    float fMaxP, fMaxV, fMaxM;
    float fSpacing;

    CVector<float> fVFL(DOFPN), fVFR(DOFPN);

    m_ElementResponseData(ElNo).GetForces(fVFL, fVFR);

    fP = fVFL(1), fV = fVFL(2), fM = fVFL(3);

    MaxStresses(fA, fS, fSF, fP, fV, fM, fCS, fTS, fSS);

    fMaxCS = fCS, fMaxTS = fTS, fMaxSS = fSS;
    fMaxP = abs(fP); fMaxV = abs(fV); fMaxM = abs(fM);


    fSpacing = fL / float(nSections);

    float fX;
    for (int SNo = 2; SNo <= nSections + 1; SNo++)
    {
        fX = float(SNo - 1) * fSpacing;

        ComputeSectionalLoads(ElNo, fL, fX, fP, fV, fM);

        fP += fVFL(1), fV += fVFL(2); fM += fVFL(3);
        fM -= fVFL(2) * fX;

        if (abs(fP) > fMaxP)
            fMaxP = abs(fP);
        if (abs(fV) > fMaxV)
            fMaxV = abs(fV);
        if (abs(fM) > fMaxM)
            fMaxM = abs(fM);

        MaxStresses(fA, fS, fSF, fP, fV, fM, fCS, fTS, fSS);

        if (fCS > fMaxCS)
            fMaxCS = fCS;
        if (fTS > fMaxTS)
            fMaxTS = fTS;
        if (fSS > fMaxSS)
            fMaxSS = fSS;
    }


    m_ElementResponseData(ElNo).SetMaxStresses(fMaxCS, fMaxTS, fMaxSS);
    m_ElementResponseData(ElNo).SetMaxForces(fMaxP, fMaxV, fMaxM);
}


void CFrame::MaxStresses (const float fA, const float fS, const float fSF,
                          const float fP, const float fV, const float fM,
                          float& fCS, float& fTS, float& fSS)
//---------------------------------------------------------------------------
//Function: computes maximum stress values at a section
//Input:    fA - Area, fS - Section Modulus, fP - Axial force
//          fV - Shear Force, fM - Moment
//Output:   none
//---------------------------------------------------------------------------
{
    fCS = abs(std::min((-fP / fA) - (fabs(fM) / fS), 0.0f));
    fTS = std::max((-fP / fA) + (fabs(fM) / fS), 0.0f);
    fSS = fabs(fV) / fSF;
}


void CFrame::ObtainSupportReactions ()
// ---------------------------------------------------------------------------
// Function: Finds the absolute and relative error norms
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // --------------------------------------------------------
// Finding Error Norms and Extracting Support Reactions
// --------------------------------------------------------

    float fFNorm = 0.0f;
    m_fAbsErrNorm = 0.0f;
    CVector<float> fVSR(DOFPN);
    fVSR.Set(0.0);

    for (int NdNo = 1; NdNo <= m_nNodes; NdNo++)
    {
        for (int i = 1; i <= DOFPN; i++)
        {
            int DOF = (NdNo - 1) * DOFPN + i;
            int ADOF = m_SADOF(DOF);
            if (ADOF == 0)
            {
                fVSR(i) = m_fVR(DOF);
            }
            else
            {
                fFNorm += float(m_SNF(ADOF) * m_SNF(ADOF)); //adding the forces square
                m_fAbsErrNorm += (m_fVR(DOF) * m_fVR(DOF)); //adding the residual square
            }
        }
        m_NodalResponseData(NdNo).SetReactions(fVSR);
    }

    fFNorm = sqrt(fFNorm);
    m_fAbsErrNorm = sqrt(m_fAbsErrNorm);
    m_fRelErrNorm = m_fAbsErrNorm / fFNorm;
}


void CFrame::ComputeSectionalLoads (const int ElNo, const float fL, const float fX,
                                    float& fAF, float& fSF, float& fBM)
// ---------------------------------------------------------------------------
// Function: computes the element nodal forces due to loads at a section
// Input:    ElNo - Element Number, fX - section distance from left
// Output:   fAF - Axial force, fSF - Shear Force, fBM - Bending Moment
// ---------------------------------------------------------------------------
{
    CElementLoads::ELType LType;
    int nETag;
    float fVal1, fVal2;
    
    fAF = 0.0f; fSF = 0.0f; fBM = 0.0f;

    for (int i = 1; i <= m_nElementLoads; i++)
    {
        m_ElementLoadData(i).GetValues(nETag, LType, fVal1, fVal2);

        if (ElNo == nETag)
        {
            if (LType == CElementLoads::ELType::DISTLOAD)
            {
                float wL = fVal1;
                float wR = fVal2;
                float wX = wL + ((wR - wL) * fL / fX);
                fSF += 0.50f * (wL + wX) * fX;
                fBM -= ((wL / 2.0f) + ((wX - wL) / 3.0f)) * fX * fX;
            }
            if (LType == CElementLoads::ELType::CONCLOADX)
            {
                float a = fVal1;
                float P = fVal2;
                if (a < fX)
                    fAF += P;
            }
            if (LType == CElementLoads::ELType::CONCLOADY)
            {
                float a = fVal1;
                float P = fVal2;
                if (a < fX)
                {
                    fSF += P;
                    fBM -= P * (fX - a);
                } 
            }
            if (LType == CElementLoads::ELType::CONCMOMENT)
            {
                float a = fVal1;
                float M = fVal2;
                if (a < fX)
                    fBM += M;
            }
        }
    }
}


void CFrame::ErrorHandler (CLocalErrorHandler::ERRORCODE ErrorCode)
{
    throw ErrorCode;
}

void CFrame::ErrorHandler (CGlobalErrorHandler::ERRORCODE ErrorCode) const
{
    throw ErrorCode;
}

void CFrame::DisplayErrorMessage (CLocalErrorHandler::ERRORCODE err)
{
    m_LEH.ErrorHandler (err, m_nLineNumber);
}
