/*********************************************
Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "frame.h"
#include "getinteractiveEXH.h"
#include "parserEXH.h"
#include "fileioEXH.h"
#include "clockEXH.h"
#include "printtableEXH.h"


void CFrame::Banner (std::ostream& OF)
// ---------------------------------------------------------------------------
// Function: Prints program banner
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    OF << '\n';
    OF << "\t\t------------------------------------------------" << '\n';
    OF << "\t\t           Planar Frame Analysis Program        " << '\n';
    OF << "\t\t  Introduction to Structural Analysis & Design  " << '\n';
    OF << "\t\t             (c) 2000-22, S. D. Rajan           " << '\n';
    OF << "\t\t     Enhanced By: Dharanidharan Arumugam        " << '\n';
    OF << "\t\t------------------------------------------------" << '\n';
}

void CFrame::PrepareIO (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
    if (argc == 1)
    {
        // open the input file
        m_FIO.OpenInputFileByName ("Complete input file name: ", m_FileInput,
                                   std::ios::in);

        // open the output file
        m_FIO.OpenOutputFileByName ("Complete output file name: ", m_FileOutput,
                                    std::ios::out);
    }
    else if (argc == 3) // planar frame input_file output_file
    {
        m_FileInput.open (argv[1], std::ios::in);
        if (!m_FileInput)
            ErrorHandler (CGlobalErrorHandler::ERRORCODE::CANNOTOPENIFILE);
        m_FileOutput.open (argv[2], std::ios::out);
        if (!m_FileOutput)
            ErrorHandler (CGlobalErrorHandler::ERRORCODE::CANNOTOPENOFILE);
        std::cout << "\n";
        std::cout << argv[1] << " opened as input file.\n";
        std::cout << argv[2] << " opened as output file.\n";
        m_strOutFileName = argv[2];
    }
	else
    {
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDCOMMANDLINE);
    }

    // print banner
    Banner (m_FileOutput);


    std::string strDateTime, strNextLine;
    m_Clock.GetDateTime(strDateTime);
    strNextLine = "\t\tReport created on " + strDateTime;

    m_FileOutput << strNextLine << "\n";
}

void CFrame::ReadProblemSize ()
// ---------------------------------------------------------------------------
// Function: reads the size of the problem from input file but does not
//           store data
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    bool bEOF = false;
    try
    {
        // read the problem description
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);
        if (m_strVTokens[0] != "*heading")
            IOErrorHandler (ERRORCODE::INVALIDINPUT);

        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);

        // nodal coordinates
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);
        if (m_strVTokens[0] != "*nodal" && m_strVTokens[1] != "coordinates")
            IOErrorHandler (ERRORCODE::INVALIDINPUT);
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*nodal" && m_strVTokens[1] == "fixity")
                break;
            ++m_nNodes;
        }

        // nodal fixity
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*nodal" && m_strVTokens[1] == "loads")
                break;
        }

        // nodal loads
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*material" && m_strVTokens[1] == "data")
                break;

            m_nNodalLoads++;
        }

        // material data
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*cross-sectional" && m_strVTokens[1] == "data")
                break;
            ++m_nMatGroups;
        }

        // cross-sectional data
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*element" && m_strVTokens[1] == "data")
                break;
            ++m_nEPGroups;
        }

        // element data
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*element" && m_strVTokens[1] == "loads")
                break;
            ++m_nElements;
        }

        // element loads
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_strVTokens[0] == "*end")
                break;
            ++m_nElementLoads;
        }
    }
    // trap all input file errors here
    catch (CLocalErrorHandler::ERRORCODE &err)
    {
        m_LEH.ErrorHandler (err);
    }

    catch (CGlobalErrorHandler::ERRORCODE &err)
    {
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDINPUT);
    }

    catch (...)
    {
        std::cout << "Sorry, could not catch the error whatever it is.\n";
    }

    // check data for validity
    if (m_nNodes <= 1) 
        IOErrorHandler (ERRORCODE::NUMNODES);
    if (m_nElements <= 0) 
        IOErrorHandler (ERRORCODE::NUMELEMENTS);
    if (m_nMatGroups <= 0) 
        IOErrorHandler (ERRORCODE::NUMMATGROUPS); 
    if (m_nEPGroups <= 0) 
        IOErrorHandler (ERRORCODE::NUMEPROPGROUPS); 
    if (m_nDebugLevel < 0 || m_nDebugLevel > 1) 
        IOErrorHandler (ERRORCODE::DEBUGCODE);
}

void CFrame::ReadFrameModel ()
// ---------------------------------------------------------------------------
// Function: reads the rest of the frame data from input file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int nTag;                // stores tag
    CVector<float> fVC(3);   // stores float input
    bool bEOF = false;       // end of file flag

    try
    {
        // rewind the file to read the input file again 
        m_FIO.Rewind (m_FileInput);
        m_nLineNumber = 0;


        // --------------------------------------------------------
        //   READING HEADER BLOCK 
        // --------------------------------------------------------          
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);


        // --------------------------------------------------------
        // READING NODAL COORDINATES
        // --------------------------------------------------------
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment, 
                           bEOF);
        for (int i = 1; i <= m_nNodes; i++)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_nTokens != 3)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);
            if (!m_Parse.GetIntValue(m_strVTokens[0], nTag) ||
                !m_Parse.GetFloatValue(m_strVTokens[1], fVC(1)) ||
                !m_Parse.GetFloatValue(m_strVTokens[2], fVC(2)))
                IOErrorHandler (ERRORCODE::INVALIDINPUT);
            if (nTag <= 0 || nTag > m_nNodes)
                IOErrorHandler (ERRORCODE::NODENUMBER);
            m_NodalData(nTag).SetCoords (fVC);
        }


        // --------------------------------------------------------
        // READING NODAL FIXITY CONDITIONS
        // --------------------------------------------------------
        
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);

        CVector<CNode::Fixity> VFC(DOFPN);
        CVector<int> nVFC;
        nVFC.SetSize(3, 0);
        int nHC = 0;

        std::string strArray[3] = {  "free", "specified", "hinged" };
        std::string strNFType;
        bool bInList;

        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);

            if (m_strVTokens[0] == "*nodal" && m_strVTokens[1] == "loads")
                break;

            if (m_nTokens != 7)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);

            if (!m_Parse.GetIntValue(m_strVTokens[0], nTag) ||
                !m_Parse.GetFloatValue(m_strVTokens[4], fVC(1)) ||
                !m_Parse.GetFloatValue(m_strVTokens[5], fVC(2)) ||
                !m_Parse.GetFloatValue(m_strVTokens[6], fVC(3)))
                IOErrorHandler (ERRORCODE::INVALIDINPUT);

            if (nTag <= 0 || nTag > m_nNodes)
                IOErrorHandler (ERRORCODE::NODENUMBER);

            // looping through degrees of freedom at node
            for (int i = 1; i <= 3; i++)
            {
                int DOF = (nTag - 1) * DOFPN + i;
                strNFType = m_strVTokens[i];
                //check whether the fixity type is available
                bInList = IsInStringsList(strNFType, strArray);
                // if not raise error
                if (!bInList)
                    IOErrorHandler (ERRORCODE::NODALFIXITY);

                if (strNFType == strArray[0])
                { 
                    VFC(i) = CNode::Fixity::FREE;
                }
                else if (strNFType == strArray[1])
                {
                    VFC(i) = CNode::Fixity::SPECIFIED;
                    m_SCDOF(DOF) = true; // true if there is a constraint
                    nVFC(i)++; // increment the number of constraints
                }
                else
                { 
                    if (i != 3)
                        IOErrorHandler(ERRORCODE::NODALFIXITY);
                    else
                    {
                        VFC(i) = CNode::Fixity::HINGED;
                        m_SCDOF(DOF) = true; //hinged condition is considered as a constraint
                        nHC++;
                    }
                }
            }
            m_NodalData(nTag).SetFixity (VFC);
            m_NodalResponseData(nTag).SetDisplacements(fVC);
        }

        // check for the minimum constraint conditions
        int nTFC = nVFC(1) + nVFC(2) + nVFC(3);
        if (nVFC(1) < 1 || nVFC(2) < 1 || nTFC < 3)
            IOErrorHandler(ERRORCODE::UNSTABLE);

        // effective degrees of freedom
        m_nEDOF = m_nDOF - nTFC - nHC;

        // number the degrees of freedom
        int ADOF = 1;
        for (int DOF = 1; DOF <= m_nDOF; DOF++)
        {
            if (!m_SCDOF(DOF))
            {
                m_SADOF(DOF) = ADOF;
                ADOF++;
            }
        }

        // set the sytem size
        SetSystemSize();

        // ---------------------------------------------------------------
        // READING NODAL LOADS
        // ----------------------------------------------------------------
        
        float fdT; // declaring temperature change at node
        int K = 1;
        for (;;)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);

            if (m_strVTokens[0] == "*material" && m_strVTokens[1] == "data")
                break;

            if (m_nTokens != 5)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);

            m_Parse.GetIntValue(m_strVTokens[0], nTag);
            m_Parse.GetFloatValue(m_strVTokens[1], fVC(1));
            m_Parse.GetFloatValue(m_strVTokens[2], fVC(2));
            m_Parse.GetFloatValue(m_strVTokens[3], fVC(3));
            m_Parse.GetFloatValue(m_strVTokens[4], fdT);

            if (nTag <= 0 || nTag > m_nNodes)
                IOErrorHandler (ERRORCODE::NODENUMBER);

            m_NodalLoadData(K).SetValues(nTag, fVC);
            m_NodalData(nTag).SetNodalTemperature(fdT);

            // looping through nodal degrees of freedom
            for (int i = 1; i <= DOFPN; i++)
            {
                int DOF = (nTag - 1) * DOFPN + i;
                int ADOF = m_SADOF(DOF);
                
                if (ADOF > 0) // if the dof is active, store the nodal forces
                {
                    m_SNF(ADOF) =  fVC(i);
                    m_fVR(DOF)  = -fVC(i); //residual vector
                }
            }

            K++;
        }

        // ----------------------------------------------------------------
        // READING MATERIAL DATA
        // ----------------------------------------------------------------
        for (int i = 1; i <= m_nMatGroups; i++)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_nTokens != 4)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);
            m_Parse.GetIntValue(m_strVTokens[0], nTag);
            m_Parse.GetFloatValue(m_strVTokens[1], fVC(1));
            m_Parse.GetFloatValue(m_strVTokens[2], fVC(2));
            m_Parse.GetFloatValue(m_strVTokens[3], fVC(3));
            if (fVC(1) <= 0.0f)
                IOErrorHandler (ERRORCODE::MATPROPERTY);
            if (nTag <= 0 || nTag > m_nMatGroups)
                IOErrorHandler (ERRORCODE::MATGROUPNUMBER);
            m_MaterialData(nTag).SetYM  (fVC(1)); 
            m_MaterialData(nTag).SetCTE (fVC(3));
        }


        // ---------------------------------------------------------------
        // READING CROSS SECTIONAL DATA
        // ---------------------------------------------------------------
        std::string strTag;
        CVector<float> fVXSDims(MAXEPDIM);
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);
        for (int i = 1; i <= m_nEPGroups; i++)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            m_Parse.GetIntValue(m_strVTokens[0], nTag);
            if (nTag < 0 || nTag > m_nEPGroups)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);
            if (m_strVTokens[1] == "rects")
            {
                if (m_nTokens != 4)
                    IOErrorHandler (ERRORCODE::INVALIDINPUT);
                m_Parse.GetFloatValue(m_strVTokens[2], fVXSDims(1));
                m_Parse.GetFloatValue(m_strVTokens[3], fVXSDims(2));
                if (fVXSDims(1) <= 0.0f || fVXSDims(2) <= 0.0f)
                    IOErrorHandler (ERRORCODE::XSDIMENSION);
                m_EPData(nTag) = new CRectSolid (fVXSDims);
            }
            else if (m_strVTokens[1] == "circs")
            {
                if (m_nTokens != 3)
                    IOErrorHandler(ERRORCODE::INVALIDINPUT);
                m_Parse.GetFloatValue(m_strVTokens[2], fVXSDims(1));
                if (fVXSDims(1) <= 0.0f)
                    IOErrorHandler(ERRORCODE::XSDIMENSION);
                m_EPData(nTag) = new CCircSolid (fVXSDims);
            }
            else if (m_strVTokens[1] == "isection")
            {
                if (m_nTokens != 6)
                    IOErrorHandler (ERRORCODE::INVALIDINPUT);
                m_Parse.GetFloatValue(m_strVTokens[2], fVXSDims(1));
                m_Parse.GetFloatValue(m_strVTokens[3], fVXSDims(2));
                m_Parse.GetFloatValue(m_strVTokens[4], fVXSDims(3));
                m_Parse.GetFloatValue(m_strVTokens[5], fVXSDims(4));
                if (fVXSDims(1) <= 0.0f || fVXSDims(2) <= 0.0f ||
                    fVXSDims(3) <= 0.0f || fVXSDims(4) <= 0.0f)
                    IOErrorHandler (ERRORCODE::XSDIMENSION);
                m_EPData(nTag) = new CISection (fVXSDims);
            }
            else
                IOErrorHandler (ERRORCODE::XSTYPE);
        }


        // ----------------------------------------------------------------
        // READING ELEMENT DATA
        // -----------------------------------------------------------------
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);

        int nSN, nEN, nMatGrp, nEPGrp;
        for (int i = 1; i <= m_nElements; i++)
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_nTokens != 5)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);
            m_Parse.GetIntValue(m_strVTokens[0], nTag);
            if (nTag <= 0 || nTag > m_nElements)
                IOErrorHandler (ERRORCODE::ELEMENTNUMBER);
            m_Parse.GetIntValue(m_strVTokens[1], nSN);
            m_Parse.GetIntValue(m_strVTokens[2], nEN);
            m_Parse.GetIntValue(m_strVTokens[3], nMatGrp);
            m_Parse.GetIntValue(m_strVTokens[4], nEPGrp);
            if (nMatGrp <= 0 || nMatGrp > m_nMatGroups)
                IOErrorHandler (ERRORCODE::MATGROUPNUMBER);
            if (nEPGrp <= 0 || nEPGrp > m_nEPGroups)
                IOErrorHandler (ERRORCODE::EPROPNUMBER);
            if (nSN <= 0 || nSN > m_nNodes)
                IOErrorHandler (ERRORCODE::NODENUMBER);
            if (nEN <= 0 || nEN > m_nNodes)
                IOErrorHandler (ERRORCODE::NODENUMBER);
            m_ElementData(nTag).SetMatPropertyGroup (nMatGrp);
            m_ElementData(nTag).SetEPropertyGroup (nEPGrp, m_EPData(nEPGrp));
            m_ElementData(nTag).SetENodes (nSN, nEN);
            //finding the element length
            float fL = m_NodalData(nEN).ToDistance(m_NodalData(nSN));
            m_ElementData(nTag).SetLength(fL);
            //finding the direction cosines
            CVector<double> dVDC(NDIM);
            m_NodalData(nSN).DirectionCosines(m_NodalData(nEN), dVDC);
            double dl = dVDC(1), dm = dVDC(2);
            m_ElementData(nTag).SetDirectionCosines(dl, dm);
        }


        // -------------------------------------------------------------
        // READ ELEMENT LOADS
        // -------------------------------------------------------------
        m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);
        std::string strEType;
        CVector<float> fVELoads(2);

        K = 1;
        for (; ; )
        {
            m_Parse.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);

            if (m_strVTokens[0] == "*end")
                break;

            if (m_nTokens != 4)
                IOErrorHandler (ERRORCODE::INVALIDINPUT);

            m_Parse.GetIntValue(m_strVTokens[0], nTag);

            if (nTag <= 0 || nTag > m_nElements)
                IOErrorHandler (ERRORCODE::ELEMENTNUMBER);

            strEType = m_strVTokens[1];
            std::string strETypeList[4] = {"dly'", "ploadx'", "ploady'", "cmoment"};

            bInList = IsInStringsList(strEType, strETypeList);
            if (!bInList)
                IOErrorHandler(ERRORCODE::ELOADTYPE);

            m_Parse.GetFloatValue(m_strVTokens[2], fVELoads(1));
            m_Parse.GetFloatValue(m_strVTokens[3], fVELoads(2));
                
            if (strEType == "dly'")
                m_ElementLoadData(K).SetValues(nTag, CElementLoads::ELType::DISTLOAD,
                                               fVELoads(1), fVELoads(2));
            else if (strEType == "ploadx'")
                m_ElementLoadData(K).SetValues(nTag, CElementLoads::ELType::CONCLOADX,
                                               fVELoads(1), fVELoads(2));
            else if (strEType == "ploady'")
                m_ElementLoadData(K).SetValues(nTag, CElementLoads::ELType::CONCLOADY,
                                               fVELoads(1), fVELoads(2));
            else
                m_ElementLoadData(K).SetValues(nTag, CElementLoads::ELType::CONCMOMENT,
                                               fVELoads(1), fVELoads(2));

            // finding the element length
            float fL = m_ElementData(nTag).GetLength();
            // finding the direction cosines
            double dl, dm;
            m_ElementData(nTag).GetDirectionCosines(dl, dm);

            //Getting the fixities of the start and end node in VFC1 and VFC2 resp.
            CVector<CNode::Fixity> VFC1(DOFPN), VFC2(DOFPN);
            m_ElementData(nTag).GetENodes(nSN, nEN);
            m_NodalData(nSN).GetFixity(VFC1);
            m_NodalData(nEN).GetFixity(VFC2);

            //finding the boolean value indicating the left and right node is hinged or not
            bool bLH = (VFC1(3) == CNode::Fixity::HINGED ? true : false);
            bool bRH = (VFC2(3) == CNode::Fixity::HINGED ? true : false);

            CVector<float> fVENFL1(3), fVENFL2(3), fVENFG1(3), fVENFG2(3);
            m_ElementLoadData(K).GetENF(fVENFL1, fVENFL2, fVENFG1, fVENFG2,
                                        fL, dl, dm, bLH, bRH); // getting the equivalent nodal forces

            //storing the equivalent nodal forces in a matrix
            for (int i = 1; i <= 3; i++)
            {
                m_ELL(i, nTag) += fVENFL1(i);
                m_ELG(i, nTag) += fVENFG1(i);
            }
            for (int i = 1; i <= 3; i++)
            {
                m_ELL(i+DOFPN, nTag) += fVENFL2(i);
                m_ELG(i+DOFPN, nTag) += fVENFG2(i);
            }

            K++;
        }

        if (m_nDebugLevel == 1)
        {
            CMatToolBox<float> MTB;
            MTB.PrintMatrixRowWise(m_ELL, "Element Local ENFs", m_FileOutput);
            MTB.PrintMatrixRowWise(m_ELG, "Element Global ENFs", m_FileOutput);
        }

    }

    // trap all input file errors here
    catch (CLocalErrorHandler::ERRORCODE &err)
    {
        m_LEH.ErrorHandler (err);
    }

    catch (CGlobalErrorHandler::ERRORCODE &err)
    {
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDINPUT);
    }

    catch (...)
    {
        std::cout << "Sorry, could not catch the error whatever it is.\n";
    }
}

void CFrame::CreateOutput()
// ---------------------------------------------------------------------------
// Function: creates the output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_FileOutput.close();


    CPrintTable LoadPrintProblemSize(1, m_strOutFileName);
    int TFORMAT = CPrintTable::TFORMAT::CENTERTEXT;
    std::string strNextLine;

    // Print Problem Size data
    LoadPrintProblemSize.PrintNextLine(TFORMAT, "---------------------------------------------");
    LoadPrintProblemSize.PrintNextLine(TFORMAT, "                PROBLEM SIZE                 ");
    LoadPrintProblemSize.PrintNextLine(TFORMAT, "---------------------------------------------");
    LoadPrintProblemSize.SkipLines(1);

    strNextLine = "Number of Nodes                      : " + std::to_string(m_nNodes);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    strNextLine = "Number of Elements                   : " + std::to_string(m_nElements);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    strNextLine = "Number of Material Groups            : " + std::to_string(m_nMatGroups);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    strNextLine = "Number of Element Property Groups    : " + std::to_string(m_nEPGroups);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    strNextLine = "Number of DOF                        : " + std::to_string(m_nDOF);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    strNextLine = "Number of Effective DOF              : " + std::to_string(m_nEDOF);
    LoadPrintProblemSize.PrintNextLine(TFORMAT, strNextLine);
    LoadPrintProblemSize.SkipLines(1);
    LoadPrintProblemSize.AllDone();

    PrintPropertyData();

    PrintNodalData();

    PrintElementData();

    PrintResponseData();

}

void CFrame::PrintPropertyData()
// ---------------------------------------------------------------------------
// Function: Prints property data to output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    {
        // Print Material Property Data
        int NUMCOLUMNS = 3;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Group #", "Young's Modulus", "CTE" };
        int nWidth1[] = { 10, 20, 20 };
        CPrintTable LoadMatPropTable(NUMCOLUMNS, m_strOutFileName);
        LoadMatPropTable.PrintTableHeading("MATERIAL PROPERTIES", strColumnHeadings,
                                            nWidth1, nAlignment);

        for (int MatGrpNo = 1; MatGrpNo <= m_nMatGroups; MatGrpNo++)
        {
            float YM = m_MaterialData(MatGrpNo).GetYM();
            float CTE = m_MaterialData(MatGrpNo).GetCTE();

            LoadMatPropTable.PrintNext(MatGrpNo);
            LoadMatPropTable.PrintNext(YM);
            LoadMatPropTable.PrintNext(CTE);

        }
        LoadMatPropTable.SkipLines(1);
        LoadMatPropTable.AllDone();
    }

    {
    // Print Element Property Data
    int NUMCOLUMNS = 3;
    int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
    std::string strColumnHeadings[] = { "Group #", "Type", "Dimensions" };
    int nWidth[] = { 10, 15, 20 };
    CPrintTable LoadElemPropTable(NUMCOLUMNS, m_strOutFileName);
    LoadElemPropTable.PrintTableHeading("ELEMENT SECTIONS", strColumnHeadings, nWidth, nAlignment);

    std::string strType;
    CXSType* pXSGrp;
    for (int ElGrpNo = 1; ElGrpNo <= m_nEPGroups; ElGrpNo++)
    {

        LoadElemPropTable.PrintNext(ElGrpNo);
        pXSGrp = m_EPData(ElGrpNo);
        CVector<float> fVDims(MAXEPDIM);
        m_EPData(ElGrpNo)->GetDimensions(fVDims);

        if (CRectSolid* pRS = dynamic_cast<CRectSolid*> (pXSGrp))
        {
            LoadElemPropTable.PrintNext("Rect. Solid");
            LoadElemPropTable.PrintNext("Height - " + std::to_string(fVDims(1)));
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("Width  - " + std::to_string(fVDims(2)));
        }
        else if (CISection* pIS = dynamic_cast<CISection*> (pXSGrp))
        {
            LoadElemPropTable.PrintNext("I - Section");
            LoadElemPropTable.PrintNext("Web Ht.    - " + std::to_string(fVDims(1)));
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("Web Tk.    - " + std::to_string(fVDims(2)));
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("Flange Ht. - " + std::to_string(fVDims(3)));
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("");
            LoadElemPropTable.PrintNext("Flange Tk. - " + std::to_string(fVDims(4)));
        }
        else if (CCircSolid* pCS = dynamic_cast<CCircSolid*> (pXSGrp))
        {
            LoadElemPropTable.PrintNext("Circ. Solid");
            LoadElemPropTable.PrintNext("Radius - " + std::to_string(fVDims(1)));
        }
    }
    LoadElemPropTable.SkipLines(1);
    LoadElemPropTable.AllDone();
}

{
        // Print Cross Sectional Property Data
        int NUMCOLUMNS = 5;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Group #", "Sectional Area", "Moment of Inertia",
                                            "Shear Factor", "Section Modulus"};
        int nWidth[] = { 10, 20, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("CROSS SECTIONAL PROPETIES", strColumnHeadings,
            nWidth, nAlignment);
        for (int ElGrpNo = 1; ElGrpNo <= m_nEPGroups; ElGrpNo++)
        {
            std::string strType;
            CXSType* pXSGrp;
            float fIyy, fIzz, fSyy, fSzz, fSFy, fSFz;
            float fA, fI, fS, fSF;

            pXSGrp = m_EPData(ElGrpNo);
            m_EPData(ElGrpNo)->GetProperties(fA, fIyy, fIzz);
            m_EPData(ElGrpNo)->GetSectionModuli(fSyy, fSzz);
            m_EPData(ElGrpNo)->GetShearFactors(fSFy, fSFz);

            //Finding MI, SM, SF for the bending and shear
            fI = fIzz; fS = fSzz; fSF = fSFy;
            if (fIyy > fIzz)
            {
                fI = fIyy;
                fS = fSyy;
                fSF = fSFz;
            }

            LoadTable.PrintNext(ElGrpNo);
            LoadTable.PrintNext(fA);
            LoadTable.PrintNext(fI);
            LoadTable.PrintNext(fSF);
            LoadTable.PrintNext(fS);
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

}



void CFrame::PrintNodalData()

{
    {
        // Print Nodal Coordinates
        int NUMCOLUMNS = 3;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Node", "X Coordinate", "Y Coordinate" };
        int nWidth[] = { 10, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("NODAL COORDINATES", strColumnHeadings,
            nWidth, nAlignment);
        for (int NdNo = 1; NdNo <= m_nNodes; NdNo++)
        {
            CVector<float> fC(NDIM);
            m_NodalData(NdNo).GetCoords(fC);

            LoadTable.PrintNext(NdNo);
            LoadTable.PrintNext(fC(1));
            LoadTable.PrintNext(fC(2));
        }
        PrintStatistics(LoadTable, NUMCOLUMNS, "Node");
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Nodal Fixities
        int NUMCOLUMNS = 7;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Node", "X-Fixity", "X-Disp.","Y-Fixity", "Y-Disp.",
                                            "RZ-Fixity", "Z-Rot." };
        int nWidth[] = { 10, 15, 15, 15, 15, 15, 15 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("NODAL FIXITIES (default is free)", strColumnHeadings,
            nWidth, nAlignment);
        for (int NdNo = 1; NdNo <= m_nNodes; NdNo++)
        {
            CVector<CNode::Fixity> VFC(DOFPN);
            m_NodalData(NdNo).GetFixity(VFC);
            CVector<float> fVD(DOFPN);
            m_NodalResponseData(NdNo).GetDisplacements(fVD);

            LoadTable.PrintNext(NdNo);

            for (int i = 1; i <= DOFPN; i++)
            {
                if (VFC(i) == CNode::Fixity::FREE)
                {
                    LoadTable.PrintNext("free");
                    LoadTable.PrintNext("");
                }
                if (VFC(i) == CNode::Fixity::HINGED)
                {
                    LoadTable.PrintNext("hinged");
                    LoadTable.PrintNext("");
                }
                if (VFC(i) == CNode::Fixity::SPECIFIED)
                {
                    LoadTable.PrintNext("specified");
                    LoadTable.PrintNext(fVD(i));
                }
            }

        }
        LoadTable.SkipLines(2);
        LoadTable.AllDone();
    }

    {
        // Print Nodal Loads
        int NUMCOLUMNS = 5;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Node", "X - Force", "Y - Force", "Z-Force", "Delta T" };
        int nWidth[] = { 10, 15, 15, 15, 15 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("NODAL LOADS (Default is no force)", strColumnHeadings,
            nWidth, nAlignment);

        CVector<float> fNL(DOFPN); int NdNo; float fdT;
        for (int i = 1; i <= m_nNodalLoads; i++)
        {

            m_NodalLoadData(i).GetLoadValues(NdNo, fNL);
            fdT = m_NodalData(NdNo).GetNodalTemperature();

            LoadTable.PrintNext(NdNo);
            LoadTable.PrintNext(fNL(1));
            LoadTable.PrintNext(fNL(2));
            LoadTable.PrintNext(fNL(3));
            LoadTable.PrintNext(fdT);
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }
}

void CFrame::PrintElementData()
// ---------------------------------------------------------------------------
// Function: Prints element data to output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{

    {
        // Print Element Data
        int NUMCOLUMNS = 5;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Start Node", "End Node", "Prop Group",
                                            "Mat Group" };
        int nWidth[] = { 10, 15, 15, 15, 15 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("ELEMENT DATA", strColumnHeadings,
            nWidth, nAlignment);

        int nSN, nEN, nEPG, nMPG;
        for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
        {
            m_ElementData(ElNo).GetENodes(nSN, nEN);
            nEPG = m_ElementData(ElNo).GetEPropertyGroupNo();
            nMPG = m_ElementData(ElNo).GetMatPropertyGroup();
            LoadTable.PrintNext(ElNo);
            LoadTable.PrintNext(nSN);
            LoadTable.PrintNext(nEN);
            LoadTable.PrintNext(nEPG);
            LoadTable.PrintNext(nMPG);
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }


    {
        // Print Element Loads
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Type", "Dist. from left",
                                            "Load Magnitude" };
        int nWidth[] = { 10, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("ELEMENT CONCENTRATED LOAD DATA",
            strColumnHeadings, nWidth, nAlignment);
        int ElNo;
        float fDistance, fMagnitude;
        CElementLoads::ELType ElemType;

        for (int nEL = 1; nEL <= m_nElementLoads; nEL++)
        {
            m_ElementLoadData(nEL).GetValues(ElNo, ElemType, fDistance, fMagnitude);

            if (ElemType != CElementLoads::ELType::DISTLOAD)
            {
                LoadTable.PrintNext(ElNo);

                if (ElemType == CElementLoads::ELType::CONCLOADX)
                    LoadTable.PrintNext("Local x conc.");
                if (ElemType == CElementLoads::ELType::CONCLOADY)
                    LoadTable.PrintNext("Local y conc.");
                if (ElemType == CElementLoads::ELType::CONCMOMENT)
                    LoadTable.PrintNext("Conc. moment");

                LoadTable.PrintNext(fDistance);
                LoadTable.PrintNext(fMagnitude);

            }
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Element Loads
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Type", "Intensity1",
                                            "Intensity2" };
        int nWidth[] = { 10, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("ELEMENT DISTRIBUTED LOAD DATA", strColumnHeadings,
            nWidth, nAlignment);

        int ElNo;
        float fVal1, fVal2;
        CElementLoads::ELType ElemType;

        for (int nEL = 1; nEL <= m_nElementLoads; nEL++)
        {
            m_ElementLoadData(nEL).GetValues(ElNo, ElemType, fVal1, fVal2);

            if (ElemType == CElementLoads::ELType::DISTLOAD)
            {
                LoadTable.PrintNext(ElNo);
                LoadTable.PrintNext("Local y dist.");
                LoadTable.PrintNext(fVal1);
                LoadTable.PrintNext(fVal2);

            }
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();

    }
}


void CFrame::PrintResponseData()
// ---------------------------------------------------------------------------
// Function: Prints element data to output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    {
        CPrintTable LoadTable(1, m_strOutFileName);
        int TFORMAT = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strNextLine;

        strNextLine = "==============================================================================";
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        LoadTable.PrintNextLine(TFORMAT, "FEA RESULTS");
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        LoadTable.SkipLines(1);
    }
    {
        // Print Nodal Displacements
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Node", "X-Displacement", "Y-Displacement","Z-Rotation" };
        int nWidth[] = { 10, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("NODAL DISPLACEMENTS", strColumnHeadings,
            nWidth, nAlignment);

        CVector<int> Nodes(m_nNodes);
        int N = 1;
        for (int NdNo = 1; NdNo <= m_nNodes; NdNo++)
        {
            CVector<CNode::Fixity> VFC(DOFPN);
            m_NodalData(NdNo).GetFixity(VFC);
            int nFree = 0;
            for (int i = 1; i <= DOFPN; i++)
            {
                if (VFC(i) == CNode::Fixity::FREE)
                    nFree++;
            }
            if (nFree > 0)
            {
                LoadTable.PrintNext(NdNo);
                CVector<float> fVD(DOFPN);
                m_NodalResponseData(NdNo).GetDisplacements(fVD);
                for (int i = 1; i <= DOFPN; i++)
                {
                    LoadTable.PrintNext(fVD(i));
                }
                Nodes(N) = NdNo;
                N++;
            }
        }
        PrintStatistics(LoadTable, Nodes, NUMCOLUMNS, "Node");
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Element Nodal Forces
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Axial", "Shear","Moment" };
        int nWidth[] = { 15, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("ELEMENT NODAL FORCES", strColumnHeadings,
            nWidth, nAlignment);

        CVector<float> fVFL(DOFPN), fVFR(DOFPN);
        for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
        {
            m_ElementResponseData(ElNo).GetForces(fVFL, fVFR);
            LoadTable.PrintNext(ElNo);
            for (int i = 1; i <= DOFPN; i++)
                LoadTable.PrintNext(fVFL(i));
            LoadTable.PrintNext("");
            for (int i = 1; i <= DOFPN; i++)
                LoadTable.PrintNext(fVFR(i));
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Maximum Element Forces
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Axial", "Shear","Moment" };
        int nWidth[] = { 15, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("MAX.(ABS) ELEMENT INTERNAL FORCES", strColumnHeadings,
            nWidth, nAlignment);

        float fMaxP, fMaxV, fMaxM;
        for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
        {
            m_ElementResponseData(ElNo).GetMaxForces(fMaxP, fMaxV, fMaxM);
            LoadTable.PrintNext(ElNo);
            LoadTable.PrintNext(fMaxP);
            LoadTable.PrintNext(fMaxV);
            LoadTable.PrintNext(fMaxM);
        }
        PrintStatistics(LoadTable, NUMCOLUMNS, "Element");
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Maximum Element Stresses
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Element", "Compression", "Tension","Shear" };
        int nWidth[] = { 15, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("MAX. ELEMENT STRESSES", strColumnHeadings,
            nWidth, nAlignment);

        float fMaxCS, fMaxTS, fMaxSS;
        for (int ElNo = 1; ElNo <= m_nElements; ElNo++)
        {
            m_ElementResponseData(ElNo).GetMaxStresses(fMaxCS, fMaxTS, fMaxSS);
            LoadTable.PrintNext(ElNo);
            LoadTable.PrintNext(fMaxCS);
            LoadTable.PrintNext(fMaxTS);
            LoadTable.PrintNext(fMaxSS);
        }
        PrintStatistics(LoadTable, NUMCOLUMNS, "Element");
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }

    {
        // Print Nodal Reactions
        int NUMCOLUMNS = 4;
        int nAlignment = CPrintTable::TFORMAT::CENTERTEXT;
        std::string strColumnHeadings[] = { "Node", "X-Reaction", "Y-Reaction","RZ-Reaction" };
        int nWidth[] = { 10, 20, 20, 20 };
        CPrintTable LoadTable(NUMCOLUMNS, m_strOutFileName);
        LoadTable.PrintTableHeading("SUPPORT(NODAL) REACTIONS", strColumnHeadings,
            nWidth, nAlignment);

        CVector<int> Nodes(m_nNodes);
        int N = 1;
        for (int NdNo = 1; NdNo <= m_nNodes; NdNo++)
        {
            CVector<CNode::Fixity> VFC(DOFPN);
            m_NodalData(NdNo).GetFixity(VFC);
            int nSpecified = 0;
            for (int i = 1; i <= DOFPN; i++)
            {
                if (VFC(i) == CNode::Fixity::SPECIFIED)
                    nSpecified++;
            }
            if (nSpecified > 0)
            {
                LoadTable.PrintNext(NdNo);
                CVector<float> fVRF(DOFPN);
                m_NodalResponseData(NdNo).GetReactions(fVRF);
                for (int i = 1; i <= DOFPN; i++)
                {
                    if (VFC(i) == CNode::Fixity::SPECIFIED)
                        LoadTable.PrintNext(fVRF(i));
                    else
                        LoadTable.PrintNext("");
                }
                Nodes(N) = NdNo;
                N++;
            }
        }
        LoadTable.SkipLines(1);
        LoadTable.AllDone();
    }
    {
        std::string strNextLine;
        int TFORMAT = CPrintTable::TFORMAT::CENTERTEXT;
        CPrintTable LoadTable(1, m_strOutFileName);
        strNextLine = "-------------------------------------------------------------------------";
        LoadTable.PrintNextLine(TFORMAT, "SOLVER PERFORMANCE");
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        strNextLine = "Solver Type          : Cholesky Solver    ";
        LoadTable.PrintNextLine(TFORMAT, strNextLine);;
        strNextLine = "Absolute Error Norm  : " + std::to_string(m_fAbsErrNorm) + "          ";
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        strNextLine = "Relative Error Norm  : " + std::to_string(m_fRelErrNorm) + "          ";
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        LoadTable.SkipLines(2);

        TFORMAT = CPrintTable::TFORMAT::ASIS;
        std::string strDateTime;
        m_Clock.GetDateTime(strDateTime);
        strNextLine = "\tEnding at " + strDateTime;
        LoadTable.PrintNextLine(TFORMAT, strNextLine);
        strNextLine = "\tElapsed time : " + std::to_string(m_Clock.DiffTime()) + "s";
        LoadTable.PrintNextLine(TFORMAT, strNextLine);

        LoadTable.AllDone();
    }

}


void CFrame::PrintStatistics (CPrintTable& LoadTable, int nColumns, std::string strTag)
// ---------------------------------------------------------------------------
// Function: Prints statistics to the output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CVector<double> dVMin(nColumns),        // min values one per column
                    dVMax(nColumns),        // max values one per column
                    dVSum(nColumns),        // sum values one per column
                    dVAbsSum(nColumns);     // sum of absolute values one per column
    CVector<int>    nVMinIndex(nColumns),   // location of min values one per column
                    nVMaxIndex(nColumns);   // location of max values one per column
    int nHit;                                 // if zero means no rows in the table 

    LoadTable.GetStatistics(dVMin, dVMax, dVSum, dVAbsSum, nVMinIndex,
                                nVMaxIndex, nHit);
    LoadTable.SkipLinesAndSpaces(1);
    LoadTable.PrintNext("Min");
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(dVMin(i));
    }
    LoadTable.PrintNext(strTag);
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(nVMinIndex(i));
    }
    LoadTable.PrintNext("Max");
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(dVMax(i));
    }
    LoadTable.PrintNext(strTag);
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(nVMaxIndex(i));
    }
}


void CFrame::PrintStatistics(CPrintTable& LoadTable, CVector<int> nVTags, 
                             int nColumns, std::string strTag)
// ---------------------------------------------------------------------------
// Function: Prints statistics to the output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CVector<double> dVMin(nColumns),        // min values one per column
                    dVMax(nColumns),        // max values one per column
                    dVSum(nColumns),        // sum values one per column
                    dVAbsSum(nColumns);     // sum of absolute values one per column
    CVector<int>    nVMinIndex(nColumns),   // location of min values one per column
                    nVMaxIndex(nColumns);   // location of max values one per column
    int nHit;                                 // if zero means no rows in the table 

    LoadTable.GetStatistics(dVMin, dVMax, dVSum, dVAbsSum, nVMinIndex,
                            nVMaxIndex, nHit);
    LoadTable.SkipLinesAndSpaces(1);
    LoadTable.PrintNext("Min");
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(dVMin(i));
    }
    LoadTable.PrintNext(strTag);
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(nVTags(nVMinIndex(i)));
    }
    LoadTable.PrintNext("Max");
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(dVMax(i));
    }
    LoadTable.PrintNext(strTag);
    for (int i = 2; i <= nColumns; i++)
    {
        LoadTable.PrintNext(nVTags(nVMaxIndex(i)));
    }
}

void CFrame::IOErrorHandler (ERRORCODE ECode) const
// ---------------------------------------------------------------------------
// Function: displays error messages related to input data
// Input:    error code
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';

    if (ECode == ERRORCODE::NUMNODES) // invalid number of nodes
        std::cerr << "Number of nodes must be >= 2.";
    else if (ECode == ERRORCODE::NUMELEMENTS) // invalid number of elements
        std::cerr << "Number of elements must be >= 1.";
    else if (ECode == ERRORCODE::DEBUGCODE) // invalid debug level
        std::cerr << "Debug level must be 0 or 1.";
    else if (ECode == ERRORCODE::NODENUMBER) // invalid node number
        std::cerr << "Invalid node number.";
    else if (ECode == ERRORCODE::ELEMENTNUMBER) // invalid element number
        std::cerr << "Invalid element number.";
    else if (ECode == ERRORCODE::XSAREA) // invalid x/s area
        std::cerr << "Area must be positive.";
    else if (ECode == ERRORCODE::YOUNGSMODULUS) // invalid E
        std::cerr << "Modulus of elasticity must be positive.";
    else if (ECode == ERRORCODE::INVALIDINPUT) // invalid input
        std::cerr << "Invalid input.";
    else if (ECode == ERRORCODE::INVALIDLASTLINE) // invalid input
        std::cerr << "Input file needs *end as last line.";
    else if (ECode == ERRORCODE::NODALFIXITY) // invalid fixity code
        std::cerr << "Nodal fixity code must be 'free' or 'specified'.";
    else if (ECode == ERRORCODE::UNSTABLE) // invalid fixity code
        std::cerr << "Unstable Structure. Not sufficient fixity conditions";
    else if (ECode == ERRORCODE::NUMMATGROUPS) // invalid number of material groups
        std::cerr << "Number of material groups must be >= 1.";
    else if (ECode == ERRORCODE::NUMEPROPGROUPS) // invalid number of property groups
        std::cerr << "Number of element property groups must be >= 1.";
    else if (ECode == ERRORCODE::EPROPNUMBER) // invalid element property group
        std::cerr << "Invalid element property group number.";
    else if (ECode == ERRORCODE::XSDIMENSION) // invalid x/s dimension
        std::cerr << "Invalid cross-section dimension.";
    else if (ECode == ERRORCODE::MATGROUPNUMBER) // invalid material group
        std::cerr << "Invalid material property group number.";
    else if (ECode == ERRORCODE::MATPROPERTY) // invalid material property
        std::cerr << "Invalid material property.";
    else if (ECode == ERRORCODE::ELOADTYPE) // invalid element load type
        std::cerr << "Invalid element load type.";
    else if (ECode == ERRORCODE::XSTYPE) // invalid x/s type
        std::cerr << "Invalid cross-section type.";
    else
        std::cerr << "Unknown error ...?";

    std::cerr << '\n' << "Error in input file line : " << m_nLineNumber;
    std::cerr << std::endl;

    exit (1);
}

void CFrame::TerminateProgram ()
// ---------------------------------------------------------------------------
// Function: terminates the program steps by closing the input/output files
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // close the input and output files
    m_FileInput.close ();
    //m_FileOutput.close ();

    std::cout << "\nExecution completed successfully." << std::endl;
}

bool CFrame::IsInStringsList(const std::string& strQuery, const std::string strList[])
    // ---------------------------------------------------------------------------
    // Function: Check whether the given string is in the string array
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
{
    int nStrings = sizeof(strList);
    bool bInList  = false;

    for (int i = 0; i < nStrings; i++)
    {
        if (strQuery == strList[i])
        {
            bInList = true;
            break;
        }      
    }

    return bInList;
}
