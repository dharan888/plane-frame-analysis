/* Matrix Toolbox

  NOTES
   (1) Illustrates the development of a matrix toolbox
       based on the CVector and CMatrix classes
       developed earlier.
   (2) Exceptions should be thrown for invalid operations

*/
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <exception>
#include "ArrayContainersEXH.h"
#include "LocalErrorHandler.h"

template <class T>
class CMatToolBox
{
    const int NUM_ELEMENTS_PER_LINE = 6;        // # of vector/matrix elements per line
    const int FW = 16;                          // field width
    const T ZERO = static_cast<T>(0);           // defining zero in data type T              
    

    public:

        CMatToolBox ();
        ~CMatToolBox ();

        // vector-related functions
        void Display (const std::string& strMessage, const CVector<T>& A) const;
        void Add (const CVector<T>& A, const CVector<T>& B, CVector<T>& C);
        void Subtract (const CVector<T>& A, const CVector<T>& B, CVector<T>& C);
        T DotProduct (const CVector<T>& A, const CVector<T>& B);
        void Normalize (CVector<T>& A);
        void Scale (CVector<T>& A, const T factor);
        T    MaxValue (const CVector<T>& A) const;
        T    MinValue (const CVector<T>& A) const;
        T    TwoNorm (const CVector<T>& A);
        T    MaxNorm (const CVector<T>& A) const;
        void CrossProduct (const CVector<T>& A, const CVector<T>& B, CVector<T>& C);

        // matrix-related functions
        void Display (const std::string& strMessage, const CMatrix<T>& A) const;
        void Add (const CMatrix<T>& A, const CMatrix<T>& B, CMatrix<T>& C);
        void Subtract (const CMatrix<T>& A, const CMatrix<T>& B, CMatrix<T>& C);
        void Multiply (const CMatrix<T>& A, const CMatrix<T>& B, CMatrix<T>& C);
        void Transpose(const CMatrix<T>& A, CMatrix<T>& B);
        void Scale (CMatrix<T>& A, const T scale);
        T    MaxNorm (const CMatrix<T>& A) const;

        void Determinant(const CMatrix<T>& A, T& det);
        void Inverse (const CMatrix<T>& A, CMatrix<T>& B, const T TOL);
        void MatMultVec (const CMatrix<T>& A, const CVector<T>& x, CVector<T>& b);

        void LUFactorization (CMatrix<T>& A, const T TOL);
        void LUSolve (const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b);
        void GaussElimination (CMatrix<T>& A, CVector<T>& x, CVector<T>& b, const T TOL);
        void LDLTFactorization (CMatrix<T>& A, const T TOL);
        void LDLTSolve (const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b);
        void LDLTFactorizationBanded(CMatrix<T>& A, T TOL);
        void LDLTSolveBanded(const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b);

        // helper function
        void MatrixAfromLU(CMatrix<T>& A);
        void MatrixAfromLDLT(CMatrix<T>& A);
        void ResidualVector (const CMatrix<T>& A, const CVector<T>& x, const CVector<T>& b,
                             CVector<T>& R, T& AbsError, T& RelError);           
        bool IsEqual (const T d1, const T d2, T TOL = TOLDEF) const;

        bool IsEqual (const CMatrix<T>& dMA, const CMatrix<T>& dMB, T TOL = TOLDEF) const;
        bool IsEqual (const CVector<T>& dVA, const CVector<T>& dVB, T TOL = TOLDEF) const;
        int Min(const int val1, const int val2);
        int Max(const int val1, const int val2);

        void GetFLOPStats (double& dAS, double& dM, double& dD) const;
        void PrintVector(const CVector<T>& A, const std::string& strMessage, std::ostream& Out) const;
        void PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading, std::ostream& Out) const;
        void PrintMatrixColumn (CMatrix<T>& A, const std::string& heading, int i, std::ostream& Out) const;
        void PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading, std::ostream& Out) const;

    private:
        // these are declared double to avoid integer overflow
        double m_dASOP; // # of floating point additions and subtractions
        double m_dMOP;  // # of floating point multiplications
        double m_dDOP;  // # of floating point divisions
        void ErrorHandler (CLocalErrorHandler::ERRORCODE);

    protected: 
       inline static const T TOLDEF = static_cast<T>(1.0e-6);      
};

// ctor
template <class T>
CMatToolBox<T>::CMatToolBox ()
// ==================================================================
// Function: default constructor
//    Input: none
//   Output: none
// ==================================================================
{
    m_dASOP = m_dMOP = m_dDOP = 0.0;
}

// dtor
template <class T>
CMatToolBox<T>::~CMatToolBox ()
// ==================================================================
// Function: destructor
//    Input: none
//   Output: none
// ==================================================================
{
}

// ---------------------------------------------------------------
// ------------------------ vector functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& strMessage,
                              const CVector<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: message string and the vector 
//   Output: None
// ==================================================================
{
    std::cout << '\n' << strMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetSize(); i++)
    {
        std::cout << "(" << i << ") "
                  << std::setw(FW) << A(i) << " ";
        if ((i % NUM_ELEMENTS_PER_LINE) == 0)
            std::cout << '\n';
    }
}

template <class T>
void CMatToolBox<T>::Add (const CVector<T>& A, const CVector<T>& B,
                          CVector<T>& C)
// ==================================================================
// Function: adds two vectors and stores the result in the
//           third vector C = A + B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_VECTORADDERROR);

    // add
    for (int i=1; i <= n; i++)
        C(i) = A(i) + B(i);
    m_dASOP += static_cast<double>(n);
}

template <class T>
void CMatToolBox<T>::Subtract (const CVector<T>& A,
                               const CVector<T>& B, CVector<T>& C)
// ==================================================================
// Function: subtracts one vector from another and stores the result
//           in the third vector C = A - B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_VECTORSUBTRACTERROR);

    // subtract
    for (int i=1; i <= n; i++)
        C(i) = A(i) - B(i);
    m_dASOP += static_cast<double>(n);
}

template <class T>
T CMatToolBox<T>::DotProduct (const CVector<T>& A,
                              const CVector<T>& B)
// ==================================================================
// Function: computes the dot product of two vectors such that
//           product = A dot B
//    Input: vectors A and B 
//   Output: product 
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize())
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_VECTORDOTPRODUCTTERROR);

    // dot product
    T product = ZERO; //initializing the product value as zero
    for (int i=1; i <= n; i++)
        product += (A(i) * B(i));

    m_dASOP += static_cast<double>(n); //updating the no. of multiplication operations
    m_dMOP  += static_cast<double>(n); //updating the no. of multiplication operations

    return product;
}

template <class T>
void CMatToolBox<T>::Normalize (CVector<T>& A)
// ==================================================================
// Function: normalizes a vector
//    Input: vector A 
//   Output: normalized vector A 
// ==================================================================
{
    // get the size of the vector
    T norm = TwoNorm(A);  

    if (IsEqual(norm, ZERO)) // if the norm is near zero raise error
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_VECTORNORMALIZEERROR);
    else // if the norm is not zero compute the normalized vector
    {
        int n = A.GetSize(); // no. of elements in the vector
        for (int i = 1; i <= n; i++)
            A(i) /= norm;
        m_dDOP += static_cast<double>(n);
    }
}

template <class T>
void CMatToolBox<T>::Scale (CVector<T>& A, T c)
// ==================================================================
// Function: scales a vector by a constant c such that A = c A
//    Input: vector A and constant c 
//   Output: scaled vector A
// ==================================================================
{
   
    // scaled vector
    int n = A.GetSize(); // find the size of A
    for (int i = 1; i <= n; i++)
        A(i) = c*A(i);

    m_dMOP += static_cast<double>(n); //updating the no. of multiplication operations
}

template <class T>
T CMatToolBox<T>::MaxValue (const CVector<T>& A) const
// ==================================================================
// Function: finds the largest value among all the elements in A
//    Input: vector A 
//   Output: return value is the largest element in A
// ==================================================================
{
    int n = A.GetSize(); // find the size of A
    T MaxVal = A(1); // initially assuming the maximum value as the first vector element

    // finding the maximum value
    for (int i = 2; i <= n; i++)
    {
        if (A(i) > MaxVal) //if the current vector element value is more than the max val
            MaxVal = A(i); //then Max val is the current vector element value
    }

    return MaxVal;
}

template <class T>
T CMatToolBox<T>::MinValue (const CVector<T>& A) const
// ==================================================================
// Function: finds the smallest value among all the elements in A
//    Input: vector A 
//   Output: return value is the smallest element in A
// ==================================================================
{
    int n = A.GetSize(); // find the size of A
    T MinVal = A(1); // initially assuming the minimum value as the first vector element

    // finding the minimum value
    for (int i = 2; i <= n; i++)
    {
        if (A(i) < MinVal) //if the current vector element value is less than the min val
            MinVal = A(i); //then Min val is the current vector element value
    }

    return MinVal;
}

template <class T>
T CMatToolBox<T>::TwoNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{
    // get the size of the vector
    int n = A.GetSize();

    // finding the norm
    T norm = ZERO; //initializing the norm as zero

    // initially finding the squared sum
    for (int i=1; i <= n; i++)
        norm += (A(i) * A(i));

    // finding the square root and casting the data type back to T
    norm = static_cast<T>(sqrt(norm));

    m_dASOP += static_cast<double>(n);
    m_dMOP += static_cast<double>(n);

    return norm;
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CVector<T>& A) const
// ==================================================================
// Function: computes the max norm of vector A
//    Input: vector A 
//   Output: return value is the max norm
// ==================================================================
{
    int n = A.GetSize(); // find the size of A
    T MaxNormVal = abs(A(1)); // initially assuming the maximum norm value as the first vector element

    // finding the max norm value
    for (int i = 2; i <= n; i++)
    {
        if (abs(A(i)) > MaxNormVal) //if the current vector element value is more than the absolute max val
            MaxNormVal = abs(A(i)); //then absolute max val is the current vector element value
    }

    return MaxNormVal;
}

template <class T>
void CMatToolBox<T>::CrossProduct (const CVector<T>& A,
                                   const CVector<T>& B,
                                   CVector<T>& C)
// ==================================================================
// Function: computes the cross-product of two vectors and stores the
//           result in the third vector such that C = A x B
//           (3-dimensional space)
//    Input: vectors A, B and C
//   Output: vector C
// ==================================================================
{
    int nA = A.GetSize();
    int nB = B.GetSize();
    int nC = C.GetSize();

    // check all the input vectors are three dimensional
    if ((nA != 3) || (nB != 3) || (nC!=3))
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_VECTORCROSSPRODUCTERROR);

    // find the cross product and assign the values to the vector C
    C(1) = A(2) * B(3) - B(2) * A(3);
    C(2) = A(3) * B(1) - B(3) * A(1);
    C(3) = A(1) * B(2) - B(1) * A(2);

    m_dASOP += static_cast<double>(3);
    m_dMOP  += static_cast<double>(6);

}

// ---------------------------------------------------------------
// ------------------------ matrix functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& strMessage,
                              const CMatrix<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a matrix
//           row wise
//    Input: message string and the matrix
//   Output: None
// ==================================================================
{
    std::cout << '\n' << strMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetRows(); i++)
    {
        int nC = 0;
        for (int j=1; j <= A.GetColumns(); j++)
        {
            ++nC;
            std::cout << "(" << i << "," << j << ") "
                      << std::setw(FW) << A(i,j) << " ";
            if ((nC % NUM_ELEMENTS_PER_LINE) == 0)
                std::cout << '\n';
        }
        std::cout << '\n';
    }
}

template <class T>
void CMatToolBox<T>::Add (const CMatrix<T>& A, const CMatrix<T>& B,
                          CMatrix<T>& C)
// ==================================================================
// Function: adds two matrices and stores the result in the
//           third matrix C = A + B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
    // Getting the dimensions of input and output matrices
    int m  = C.GetRows(); int n  = C.GetColumns(); //Getting the row and the column dimension of matrix C
    int mA = A.GetRows(); int nA = A.GetColumns(); //Getting the row and the column dimension of matrix A
    int mB = B.GetRows(); int nB = B.GetColumns(); //Getting the row and the column dimension of matrix B

    // check for incompatible vectors
    bool browschk = ((m == mA) && (m == mB));
    bool bcolschk = ((n == nA) && (n == nB));
    if (!browschk || !bcolschk)     // if column or row dimensions are not consistent
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXADDERROR); // then raise error

    // add
    for (int i = 1; i <= m; i++) // looping along row elements 
    {
        for (int j = 1; j <= n; j++) // looping along column elements
        {
            C(i, j) = A(i, j) + B(i, j);
        }
    }

    m_dASOP += static_cast<double>(m*n);

}

template <class T>
void CMatToolBox<T>::Subtract (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: subtracts one matrix from another and stores the result
//           in the third matrix C = A - B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
    // Getting the dimensions of input and output matrices
    int m = C.GetRows(); int n = C.GetColumns(); //Getting the row and the column dimension of matrix C
    int mA = A.GetRows(); int nA = A.GetColumns(); //Getting the row and the column dimension of matrix A
    int mB = B.GetRows(); int nB = B.GetColumns(); //Getting the row and the column dimension of matrix B

    // check for incompatible vectors
    bool browschk = ((m == mA) && (m == mB));
    bool bcolschk = ((n == nA) && (n == nB));
    if (!browschk || !bcolschk)     // if column or row dimensions are not consistent
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXSUBTRACTERROR);// then raise error

    // subtract
    for (int i = 1; i <= m; i++) // looping along row elements 
    {
        for (int j = 1; j <= n; j++) // looping along column elements
        {
            C(i, j) = A(i, j) - B(i, j);
        }
    }

    m_dASOP += static_cast<double>(m * n);
}

template <class T>
void CMatToolBox<T>::Multiply (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: multiplies two matrices and stores the result
//           in the third matrix C = A * B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
    // Getting the dimensions of input and output matrices
    int m  = C.GetRows(); int n  = C.GetColumns(); //Getting the row and the column dimension of matrix C
    int mA = A.GetRows(); int nA = A.GetColumns(); //Getting the row and the column dimension of matrix A
    int mB = B.GetRows(); int nB = B.GetColumns(); //Getting the row and the column dimension of matrix B

    // check for incompatible vectors
    if ((m != mA) || (n != nB) || (nA != mB))     // if column or row dimensions are not consistent
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXMULTIPLYERROR);  // then raise error

    // multiply
    T sum; //variable used to find sum of the products
    for (int i = 1; i <= m; i++) // looping along row elements of C
    {
        for (int j = 1; j <= n; j++) // looping along column elements of C
        {
            sum = ZERO; // initiazling the sum of the products is zero
            for (int k = 1; k <= nA; k++) // looping through row-i elements of A and column-j elements of B
            {
                sum += A(i, k) * B(k, j); // summation of the element products
            }
            C(i, j) = sum;
        }
    }

    m_dASOP += static_cast<double>(m * n * nA); //number of addition operations
    m_dMOP  += static_cast<double>(m * n * nA); //number of multiplication operations

}

// not needed to be coded for this current project
template <class T>
void CMatToolBox<T>::Determinant (const CMatrix<T>& A, T& det)
// ==================================================================
// Function: computes the determinant of matrix A
//    Input: matrix A and variable to hold the determinant
//   Output: determinant
// ==================================================================
{
    ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_UNSUPPORTEDOPERATION);
}


// not needed to be coded for this current project
template <class T>
void  CMatToolBox<T>::Inverse (const CMatrix<T>& A, CMatrix<T>& B, const T TOL)
// ==================================================================
// Function: computes the inverse of matrix A
//    Input: matrix A and TOl for tolerance value
//   Output: matrix B containing the A inverse
// ==================================================================
{
    ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_UNSUPPORTEDOPERATION);
}


template <class T>
void CMatToolBox<T>::Scale (CMatrix<T>& A, T c)
// ==================================================================
// Function: scales all the elements of a matrix by a constant c
//           such that A = c A
//    Input: matrix A and constant c
//   Output: scaled matrix A
// ==================================================================
{
    //Getting the row and the column dimensions of matrix A
    int m = A.GetRows(); int n = A.GetColumns(); 

    //Scaling the matrix
    for (int i = 1; i <= m; i++) // looping along row elements of A
    {
        for (int j = 1; j <= n; j++) // looping along column elements of A 
        {
            A(i, j) = c * A(i, j); //scaling the matrix elements
        }
    }


    m_dMOP += static_cast<double>(m * n); //number of addition operations
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CMatrix<T>& A) const
// ==================================================================
// Function: computes the max norm of matrix A
//    Input: matrix A 
//   Output: return value is the max norm
// ==================================================================
{
    //Getting the row and the column dimensions of matrix A
    int m = A.GetRows(); int n = A.GetColumns();

    //initializing the max norm value
    T MaxNormVal = static_cast<T>(-1.0);

    //finding the maximum norm of matrix A
    for (int i = 1; i <= m; i++) // looping along row elements of A
    {
        for (int j = 1; j <= n; j++) // looping along column elements of A 
        {
            if (abs(A(i, j)) > MaxNormVal)
                MaxNormVal = abs(A(i, j)); // Update the Max norm value to Aij since it is greater than existing norm value
        }

    }

    return MaxNormVal;
}

template <class T>
void CMatToolBox<T>::Transpose (const CMatrix<T>& A,
                                CMatrix<T>& B)
// ==================================================================
// Function: computes the transpose of a matrix and stores the result
//           in another matrix B = A(T)
//    Input: matrices A and B
//   Output: matrix B
// ==================================================================
{
    int mA = A.GetRows(); int nA = A.GetColumns(); //Getting the row and the column dimensions of matrix A
    int mB = B.GetRows(); int nB = B.GetColumns(); //Getting the row and the column dimensions of matrix B

    // check for incompatible vectors
    if ((mA != nB) || (nA != mB))     // if column or row dimensions are not consistent
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXTRANSPOSEERROR);// then raise error

    for (int i = 1; i <= mA; i++) // looping through row elements of matrix A
        for (int j = 1; j <= mB; j++) // looping through column elements of matrix A
            B(j, i) = A(i, j);
}

template <class T>
void CMatToolBox<T>::MatMultVec (const CMatrix<T>& A,
                                 const CVector<T>& x,
                                 CVector<T>& b)
// ==================================================================
// Function: multiplies a matrix and a vector and stores the result
//           in a vector b = A * x
//    Input: vectors A and x 
//   Output: vector b
// ==================================================================
{
    // Getting the dimensions of input and output matrices
    int m  = A.GetRows(); int n  = A.GetColumns(); //Getting the row and the column dimension of matrix A
    int nx = x.GetSize(); int nb = b.GetSize(); //Getting the number of elements of vector x and b, resptectively.


    // check for incompatible vectors
    if ((n != nx) || (m != nb))     // if dimensions of A, x and b are not consistent
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXMULTVECERROR);// then raise error

    // multiply
    for (int i = 1; i <= m; i++) // looping along row elements of A and C
    {
        b(i) = ZERO;
        for (int j = 1; j <= n; j++) // looping along column elements of A and C
            b(i) += A(i, j) * x(j);
    }

    m_dASOP += static_cast<double>(m * n); //number of addition operations
    m_dMOP += static_cast<double>(m * n); //number of multiplication operations
}

template <class T>
void CMatToolBox<T>::GaussElimination (CMatrix<T>& A, CVector<T>& x,
                                       CVector<T>& b, T TOL)
// ==================================================================
// Function: solves A x = b using Gaussian Elimination Technique
//           (this version only for one rhs vector)
//    Input: Matrices A and b
//   Output: Matrix x
// ==================================================================
{
    int i, j, k;   // loop indices
    T c;           // multiplier (Step 4)

    // number of equations to solve
    int n = A.GetRows();

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)              // Step 1
    {
        if (fabs(A(k,k)) <= TOL)          // Step 2
            ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXGAUSSELIMINATIONERROR);
        for (i=k+1; i <= n; i++)          // Step 3
        {
            c = A(i,k)/A(k,k);            // Step 4
            for (j=k+1; j <= n; j++)      // Step 5
                A(i,j) -= c * A(k,j);     // Step 6
            x(i) -= c * x(k);             // Step 8
        }                                 // Step 9
        int nC = n-k;
        if (nC > 0)
        {
            m_dDOP += static_cast<double>(nC);
            m_dMOP += static_cast<double>(nC*nC);
            m_dASOP += static_cast<double>(nC+nC*nC);
        }
    }                                     // Step 10 

    // back substitution
    if (fabs(A(n,n)) <= TOL)              
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXGAUSSELIMINATIONERROR);
    x(n) /= A(n,n);                       // Step 11

    for (i=n-1; i >= 1; i--)              // Step 12
    {
        T sum = T(0);
        for (j=i+1; j <= n; j++)
            sum += A(i,j) * x(j);         // Step 13
        if ((n-i) > 0)
        {
            m_dASOP += static_cast<double>(n-i);
            m_dMOP += static_cast<double>(n-i);
        }
        x(i) = (x(i) - sum)/A(i,i);       // Step 14
    }                                     // Step 15
    m_dASOP += static_cast<double>(n);
    m_dDOP += static_cast<double>(n+1);

}


template <class T>
void CMatToolBox<T>::LUFactorization (CMatrix<T>& A, T TOL)
// ==================================================================
// Function: carries out LU factorization (DoLittle method) 
//           of matrix A 
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A (A is replaced with L and U)
// ==================================================================
{
    // Getting the column dimension of the matrix A
    int n = A.GetColumns(); 

    // Checking whether the A matrix is square
    if (n != A.GetRows())
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_LUFACTDIMERROR);

    // check if the matrix is trivial 1 by 1
    if (n == 1)
        return; //if it is return the original matrix

    // Step 2
    for (int j = 1; j <= n; j++)
    {
        // Step 3
        for (int i = 1; i <= j; i++)
        {
            for (int k = 1; k < i; k++)
            {
                A(i, j) -= A(i, k) * A(k, j);
            }
            //Step 4
            if (abs(A(j, j)) <= TOL)
                ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_LUFACTSINGULARERROR);
        }
        // Step 5
        for (int i = j+1; i <= n; i++)
        {
            for (int k = 1; k < j; k++)
            {
                A(i, j) -= A(i, k) * A(k, j);
            }
            A(i, j) /= A(j, j);
        }
    }

    double DivOps = static_cast<double>(n * (n - 1)/2);
    double MulOps = static_cast<double>(n * (n - 1) * (2 * n - 1)/6);

    m_dDOP += DivOps;  // number of divisional operations
    m_dMOP += MulOps;  // number of multiplication operations
    m_dASOP += MulOps; // number of subtraction operations (same as multiplication operations)

}


template <class T>
void CMatToolBox<T>::LUSolve (const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b)            
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and U terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    // Getting the column dimension of the matrix A
    int n = A.GetColumns();

    // Checking whether the A matrix is square
    if (n != A.GetRows() || n!= x.GetSize() || n!= b.GetSize())
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXLUSOLVEERROR);

    // check if the matrix is trivial 1 by 1
    if (n == 1)
    {
        x(1) = b(1) / A(1, 1);
        m_dDOP += static_cast<double>(1);
        return;
    }

    // x initially contains b
    x = b;

    // forward substitution
    for (int i = 2; i <= n; i++)
    {
        for (int j = 1; j <= i - 1; j++)
        {
            x(i) -= A(i, j) * x(j);
        }
    }

    // backward substitution
    for (int i = n; i >= 1; i--)
    {
        for (int j = i+1; j <= n; j++)
        {
            x(i) -= A(i, j) * x(j);
        }

        x(i) /= A(i, i); 
    }

    m_dDOP  += static_cast<double>(n);  // number of divisional operations
    m_dMOP  += static_cast<double>(n * (n - 1));  // number of multiplication operations
    m_dASOP += static_cast<double>(n * (n - 1)); // number of subtraction operations (same as multiplication operations)
}

template <class T>
void CMatToolBox<T>::LDLTFactorization (CMatrix<T>& A,
                                        T TOL)
// ==================================================================
// Function: carries out LDL(T) factorization of matrix A
//           A is replaced with L and D. A is a symmetric matrix.
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
    // Getting the column dimension of the matrix A
    int n = A.GetColumns();

    // Checking whether the A matrix is square
    if (n != A.GetRows())
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXLDLTFACTERROR);

    // check if the matrix is trivial 1 by 1
    if (n == 1)
        return; //if it is return the original matrix

    // Step 2
    for (int i = 1; i <= n; i++)
    {
        // Step 3
        for (int k = 1; k <= i - 1; k++)
        {
            A(i,i) -= A(i, k) * A(i, k) * A(k, k);
        }

        if (A(i, i) <= TOL)
            ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXLDLTFACTERROR);

        // Step 3
        for (int j = i + 1; j <= n; j++)
        {
            for (int k = 1; k <= i - 1; k++)
            {
                A(j, i) -= A(j, k) * A(k, k) * A(i, k);
            }
            A(j, i) /= A(i, i);
        } 
    }

    double Ops = static_cast<double> ((pow(n,3)-n)/6);

    m_dDOP += static_cast<double>(n*(n-1)/2);  // number of divisional operations
    m_dMOP += 2*Ops;  // number of multiplication operations
    m_dASOP += Ops; // number of subtraction operations (half of multiplication operations)

}


template <class T>
void CMatToolBox<T>::LDLTSolve (const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and D terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    // Getting the column dimension of the matrix A
    int n = A.GetColumns();

    // Checking whether the A matrix is square
    if (n != A.GetRows() || n != x.GetSize() || n != b.GetSize())
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXLDLTSOLVEERROR);

    //  if the matrix is 1 by 1
    if (n == 1)
    {
        x(1) = b(1) / A(1, 1);
        return;
    }

    // x initially contains b
    x = b;

    // forward substitution
    for (int i = 2; i <= n; i++)
    {
        for (int j = 1; j <= i - 1; j++)
        {
            x(i) -= A(i, j) * x(j);
        }
    }

    // backward substitution
    for (int i = n; i >= 1; i--)
    {
        x(i) /= A(i, i);

        for (int j = i + 1; j <= n; j++)
        {
            x(i) -= A(j, i) * x(j);
        }
    }

    m_dDOP += static_cast<double>(n);  // number of divisional operations
    m_dMOP += static_cast<double>(n * (n - 1));  // number of multiplication operations
    m_dASOP += static_cast<double>(n * (n - 1)); // number of subtraction operations (same as multiplication operations)
}

template <class T>
void CMatToolBox<T>::LDLTFactorizationBanded (CMatrix<T>& A, T TOL)
    // ==================================================================
    // Function: carries out LDL(T) factorization of matrix A
    //           A is replaced with L and D. A is a symmetric matrix.
    //    Input: matrix A and tolerance value to detect singular A
    //   Output: matrix A 
    // ==================================================================
{
    // Getting the column dimension of the matrix A
    int w = A.GetColumns(); //half band width
    int n = A.GetRows(); //no. of rows

    int jL, jU, kL, kU; //upper and lower limits of indices

    // Step 2
    for (int i = 1; i <= n; i++)
    {
        kL = Max(i - w + 1, 1); kU = i - 1;

        // Step 3
        for (int k = kL; k <= kU; k++)
        {
            A(i, 1) -= A(k, i - k + 1) * A(k, 1) * A(k, i - k + 1);
        }

        if (A(i, 1) <= TOL)
            ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXLDLTFACTERROR);

        jL = i + 1; jU = Min(i + w - 1, n);

        // Step 3
        for (int j = jL; j <= jU; j++)
        {

            kL = Max(j - w + 1, 1); kU = i -1;

            for (int k = kL; k <= kU; k++)
            {
                A(i, j - i + 1) -= A(k, j - k + 1) * A(k, 1) * A(k, i - k + 1);
            }
            A(i, j - i + 1) = A(i, j - i + 1) / A(i, 1);

        }
    }

}

template <class T>
void CMatToolBox<T>::LDLTSolveBanded(const CMatrix<T>& A, CVector<T>& x, const CVector<T>& b)
    // ==================================================================
    // Function: carries out forward and backward substitution so as to
    //           solve A x = b. A contains L and D terms.
    //    Input: matrix A, vectors x and b
    //   Output: vector x 
    // ==================================================================
{
    // Getting the column dimension of the matrix A
    int n = A.GetRows();
    int w = A.GetColumns();

    // Checking whether the dimensions are consistent
    if (n != x.GetSize() || n != b.GetSize())
        ErrorHandler(CLocalErrorHandler::ERRORCODE::MTB_MATRIXLDLTDIMERROR);

    //  if the matrix is 1 by 1
    if (n == 1)
    {
        x(1) = b(1) / A(1, 1);
        m_dDOP += static_cast<double>(1);
        return;
    }

    // x initially contains b
    x = b;

    // forward substitution
    for (int i = 2; i <= n; i++)
    {
        int jL = Max(i - w + 1, 1); 

        for (int j = jL; j <= i - 1; j++)
        {
            x(i) -= A(j, i-j+1) * x(j);
        }
    }

    // backward substitution
    for (int i = n; i >= 1; i--)
    {
        x(i) /= A(i, 1);
        int jU = Min(i + w - 1, n);

        for (int j = i + 1; j <= jU; j++)
        {
            x(i) -= A(i, j-i+1) * x(j);
        }
    }

}

template <class T>
void CMatToolBox<T>::MatrixAfromLU(CMatrix<T>& A)
// ==================================================================
// Function: computes the original matrix A from L and U 
//    Input: matrix A, containing the decomposed L and U
//   Output: matrix A
// ==================================================================
{
    int n = A.GetRows();

    T sum;
    // Loop through rows of A
    for (int i = n; i >= 1; i--)
    {
        // Loop through columns of A
        for (int j = n; j >= 1; j--)
        {
            sum = ZERO;
            for (int k = 1; ;k++)
            {
                if (i == k) //check the Lik = 1
                    sum += A(k, j); 
                else
                    sum += A(i, k) * A(k, j); 

                if (k >= i || k >= j) // check if the index exceeds the lower and upper triangle max index.
                    break;
            }
            A(i, j) = sum;
        }
    }
}


template <class T>
void CMatToolBox<T>::MatrixAfromLDLT(CMatrix<T>& A)
// ==================================================================
// Function: computes the original matrix A from L and D 
//    Input: matrix A, containing the decomposed L and D
//   Output: matrix A
// ==================================================================
{
    int n = A.GetRows();

    T sum;

    // Loop through rows of A
    for (int i = n; i >= 1; i--) 
    {
        // Loop through columns of A
        for (int j = i; j >= 1; j--) 
        {
            sum = ZERO;
            for (int k = 1; ; k++) 
            {
                // check if the index exceeds the lower and upper triangle max index.
                if (k > i || k > j)
                    break;

                if ((i == k) && (j == k))
                    sum += A(k, k); // check if Lik and Lki = 1
                else if ((i == k))
                    sum += A(k, k) * A(j, k); // check if Lik = 1
                else if ((j == k))
                    sum += A(i, k) * A(k, k); // check if Ljk = 1
                else
                    sum += A(i, k) * A(k, k) * A(j, k); 
            }

            A(i, j) = sum;
      
        }

    }
}


template <class T>
void CMatToolBox<T>::ResidualVector (const CMatrix<T>& A, 
                                     const CVector<T>& x,
                                     const CVector<T>& b,
                                     CVector<T>& R,
                                     T& AbsError, T& RelError)
// ==================================================================
// Function: computes the residual vector arising from the solution
//           of A x = b, i.e., R = A x - b
//    Input: matrix A, vectors x and b, 
//   Output: vector R, abs. and rel. error values via TwoNorm function
// ==================================================================
{
    // check for incompatible sizes
    int m  = A.GetRows();
    int n  = A.GetColumns();
    int nx = x.GetSize();
    int nb = b.GetSize();

    if (n != nx || m != nb)
        ErrorHandler (CLocalErrorHandler::ERRORCODE::MTB_MATRIXRESIDUALERROR);

    MatMultVec (A, x, R);
    for (int i = 1; i <= n; i++)
    {
        R(i) -= b(i);
    }
    m_dASOP += static_cast<double>(n);

    AbsError = TwoNorm(R);

    T bnorm = TwoNorm(b);
    if (IsEqual(bnorm, ZERO))
        RelError = AbsError; 
    else 
        RelError = AbsError / bnorm; 
}


template <class T>
bool CMatToolBox<T>::IsEqual (const T d1, const T d2, T TOL) const
// ==================================================================
// Function: checks if d1 and d2 are 'nearly' equal
//    Input: d1, d2, tolerance to use
//   Output: true if they are
// ==================================================================
{
    return (abs(d1-d2) <= TOL);
}

template <class T>
bool CMatToolBox<T>::IsEqual (const CMatrix<T>& dMA,
                              const CMatrix<T>& dMB, T TOL) const
// ==================================================================
// Function: checks if matrices A and B are 'nearly' equal
//    Input: A, B, tolerance to use
//   Output: true if they are
// ==================================================================
{
    for (int i=1; i <= dMA.GetRows(); i++)
    {
        for (int j=1; j <= dMA.GetColumns(); j++)
        {
            if (!IsEqual(dMA(i,j), dMB(i,j), TOL))
                return false;
        }
    }

    return true;
}

template <class T>
bool CMatToolBox<T>::IsEqual (const CVector<T>& dVA,
                              const CVector<T>& dVB, T TOL) const
// ==================================================================
// Function: checks if vectors A and B are 'nearly' equal
//    Input: A, B, tolerance to use
//   Output: true if they are
// ==================================================================
{
    for (int i=1; i <= dVA.GetSize(); i++)
    {
        if (!IsEqual(dVA(i), dVB(i), TOL))
            return false;
    }

    return true;
}

template <class T>
int CMatToolBox<T>::Min(const int val1, const int val2)
// ---------------------------------------------------------------------------
// Function: find the minimum of two numbers
// Input:    value 1 and value 2
// Output:   min values
// ---------------------------------------------------------------------------
{
    if (val1 < val2)
        return val1;
    else
        return val2;
}

template <class T>
int CMatToolBox<T>::Max(const int val1, const int val2)
// ---------------------------------------------------------------------------
// Function: find the maximum of two numbers
// Input:    value 1 and value 2
// Output:   min values
// ---------------------------------------------------------------------------
{
    if (val1 > val2)
        return val1;
    else
        return val2;
}

template <class T>
void CMatToolBox<T>::GetFLOPStats (double& dAS, double& dM, double& dD) const
// ==================================================================
// Function: retrieves floating point operations
//    Input: variables to store +-, * and / operations
//   Output: variables with their values
// ==================================================================
{
    dAS = m_dASOP;
    dM  = m_dMOP;
    dD  = m_dDOP;
}

template <class T>
void CMatToolBox<T>::PrintVector (const CVector<T>& A, const std::string& strMessage,
                                  std::ostream& Out) const
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: vector, a heading, output stream object 
//   Output: None
// ==================================================================
{
    Out << '\n' << strMessage << '\n';
    Out.setf(std::ios::left);
    for (int i=1; i <= A.GetSize(); i+=NUM_ELEMENTS_PER_LINE)
    {
        for (int k=i; k <= std::min(i+NUM_ELEMENTS_PER_LINE-1, A.GetSize());
             k++)
        {
            Out << "[" << std::setw (4) << k << "]";
            Out << std::setw (15) << A(k) << " ";
        }
        Out << '\n';
    }
}

template <class T>
void CMatToolBox<T>::PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading,
                                         std::ostream& Out) const
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    Out << '\n' << heading << '\n';
    Out << std::setiosflags (std::ios::left);

    int nRows = A.GetRows();
    int nColumns = A.GetColumns();
    for (int i=1; i <= nRows; i++)
    {
        Out << "Row No: " << i << '\n';
        for (int j=1; j <= nColumns; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (int k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nColumns); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(i,k) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void CMatToolBox<T>::PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading,
                                            std::ostream& Out) const
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    Out << '\n' << heading << '\n';
    Out << std::setiosflags (std::ios::left);

    int nRows = A.GetRows();
    int nColumns = A.GetColumns();
    for (int i=1; i <= nColumns; i++)
    {
        Out << "Column No: " << i << '\n';
        for (int j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (int k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(k,i) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void CMatToolBox<T>::PrintMatrixColumn (CMatrix<T>& A, const std::string& heading, int i,
                                        std::ostream& Out) const
// ---------------------------------------------------------------------------
// Function: outputs a column of the matrix into stream Out
// Input:    matrix, a heading, column index, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    Out << '\n' << heading << '\n';
    Out << std::setiosflags (std::ios::left);

    int nRows = A.GetRows();
    for (int j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
    {
        for (int k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
        {
            Out << "[" << std::setw (4) << k << "]";
            Out << std::setw (15) << A(k,i) << " ";
        }
        Out << '\n';
    }
}



template <class T>
void CMatToolBox<T>::ErrorHandler (CLocalErrorHandler::ERRORCODE err)
// ---------------------------------------------------------------------------
// Function: gateway to error handling. useful for setting breakpoint in
//           the debugger
// Input:    error code
// Output:   none
// ---------------------------------------------------------------------------
{
    throw err;
}