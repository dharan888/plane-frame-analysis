/*********************************************
Planar Frame Analysis Program
*********************************************/
#pragma once
#include <fstream>
#include <iostream>
#include <functional>
#include <algorithm>
#include "constants.h"
#include "arraycontainersEXH.h"
#include "parserEXH.h"
#include "GlobalErrorHandler.h"
#include "fileioEXH.h"
#include "MatToolBox.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "xstype.h"
#include "rectsolid.h"
#include "isection.h"
#include "circsolid.h"
#include "nodalloads.h"
#include "elementloads.h"
#include "nodalresponse.h"
#include "elementresponse.h"
#include "LocalErrorHandler.h"
#include "printtableEXH.h"
#include "clockEXH.h"

class CFrame
{
    public:
        CFrame ();    // ctor
        ~CFrame ();   // dtor
        enum class ERRORCODE {NUMNODES, NUMELEMENTS, DEBUGCODE,
                              NODENUMBER, ELEMENTNUMBER, XSAREA,
                              YOUNGSMODULUS, UNSTABLE, INVALIDINPUT,
                              INVALIDLASTLINE, NODALFIXITY, NUMMATGROUPS,
                              NUMEPROPGROUPS, EPROPNUMBER, XSDIMENSION,
                              MATGROUPNUMBER, MATPROPERTY, ELOADTYPE,
                              XSTYPE, ELEMEQUILIBRIUM, INVALIDCOMMANDLINE,
		                      CANNOTOPENIFILE, CANNOTOPENOFILE};

        // helper functions
        void Banner (std::ostream& OF);
        void PrepareIO (int argc, char *argv[]);
        void Analyze ();
        void TerminateProgram ();
        void DisplayErrorMessage (CLocalErrorHandler::ERRORCODE);

    private:
        int m_nNodes;		   // number of nodes
        int m_nElements;       // number of elements
        int m_nEPGroups;       // number of x/s properties
        int m_nMatGroups;      // number of material groups
        int m_nElementLoads;   // number of element loads
        int m_nNodalLoads;     // number of nodal loads

        int m_nDOF;			   // total degrees-of-freedom
        int m_nSDOF;           // total number of suppressed degrees of freedom
        int m_nEDOF;           // total number of effective degrees of freedom
        int m_nDebugLevel;	   // debugging level

        CVector<CNode>				m_NodalData;               // nodal data
        CVector<CElement>			m_ElementData;             // element data
        CVector<CMaterial>			m_MaterialData;            // material data
        CVector<CXSType*>			m_EPData;                  // element property data
        CVector<CNodalLoads>		m_NodalLoadData;           // nodal load data
        CVector<CElementLoads>		m_ElementLoadData;         // element load data
        CVector<CNodalResponse>		m_NodalResponseData;       // nodal response data
        CVector<CElementResponse>	m_ElementResponseData;     // element response data
        

        CMatrix<double> m_SSM;	 // structural stiffness matrix
        CVector<double> m_SND;	 // structural nodal displacements
        CVector<double> m_SNF;	 // structural nodal forces

        CMatrix<float> m_ELL;	 // element loads (local coordinate system)
        CMatrix<float> m_ELG;	 // element loads (global coordinate system)

        CVector<int>  m_SADOF;    // structural active degrees of freedom
        CVector<bool> m_SCDOF;    // constrained degrees of freedom

        CVector<float> m_fVR;    // residual vector
        float m_fAbsErrNorm;     // Absolute Residual Norm
        float m_fRelErrNorm;     // Relative Residual Norm

        // input and output
        std::ifstream m_FileInput;	   // File Input
        int m_nLineNumber;	           // current line number in input file
        CParser m_Parse;               // parser for free format read
        int     m_nV;                  // integer value that is read in
        float   m_fV;                  // float value that is read in
        std::string m_strDelimiters;   // delimiters used in input file
        std::string m_strComment;      // characters to signify comment line
        std::vector<std::string> m_strVTokens; // vector to store tokens read
        int     m_nTokens;             // number of tokens read
        std::ofstream m_FileOutput;	   // File Output
        std::string m_strOutFileName;  // Output File Name
        CFileIO m_FIO;                 // for file operations
        CLocalErrorHandler m_LEH;      // for handling errors detected by program
        CClock m_Clock;                // for finding the elapsed time

        bool IsInStringsList(const std::string& strQuery, const std::string strArray[]);


        // FEA-related functions
        void ReadProblemSize ();
        void ReadFrameModel ();
        void ElementStiffness(const float fE, const float fA, const float fI,
                              const float fL, const double l, const double m,
                              const bool bLHinged, const bool bRHinged, CMatrix<double>& dMKe);
        void ConstructSystemKandF ();
        void Solve ();
        void Response ();
        void ComputeElementForces();
        void ComputePeakStresses (const int nSections, const int ElNo,
                                  const float fL, const float fA,
                                  const float fS, const float fSF);
        void ComputeSectionalLoads (const int ElNo, const float fL, const float fX,
                                    float& fAF, float& fSF, float& fBM);
        void MaxStresses (const float fA, const float fS, const float fSF, 
                          const float fP, const float fV, const float fM,
                          float& fCS, float& fTS, float& fSS);
        void CFrame::ObtainSupportReactions();
        void CreateOutput ();
        void PrintPropertyData();
        void PrintNodalData();
        void PrintElementData();
        void PrintResponseData();
        void PrintStatistics (CPrintTable& LoadTable, int nColumns, std::string strTag);
        void PrintStatistics(CPrintTable& LoadTable, CVector<int> nVTags,
                             int nColumns, std::string strTag);

        // modifier functions
        void SetSize ();
        void SetSystemSize();

        // error handlers
        void ErrorHandler (CLocalErrorHandler::ERRORCODE);        // gateway to local error handler
        void ErrorHandler (CGlobalErrorHandler::ERRORCODE) const; // gateway to global error handler
        void IOErrorHandler (ERRORCODE nCode) const;
};
