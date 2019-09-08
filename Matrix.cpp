#include <stdio.h>
#include "ARmodel.h"
using namespace std;

// constructor from vector<vector>
Matrix::Matrix(vector<vector<double>> Mat)
{
    matrix = Mat;
    size = matrix.size();
}

// constructor for a lagged correlation matrix from a correlation vector
// if correlation vector is [a, b, c]
// correlation matrix is [a, b, c]
//                       [b, a, b]
//                       [c, b, a]
Matrix::Matrix(vector<double> seed)
{
    // set size of matrix
    size = seed.size();
    // initialize with matrix of 0s
    vector<vector<double>> mat(size, vector<double> (size, 0));
    
    // loop by row
    for (int row=0; row<size; row++)
    {
        // loop from index row to end (forward filling: seed[0]->seed[size-row])
        for (int i=0; i<size-row; i++)
        {
            mat[row][row+i] = seed[i];
        }
        // loop from index row-1 to 0 (backward filling: seed[0]->seed[row])
        for (int i=0; i<row; i++)
        {
            mat[row][row-i-1] = seed[i+1];
        }
    }
    // save
    matrix = mat;
}

void Matrix::setRow(double* data, int row)
{
    for (int i=0; i<size; i++)
    {
        matrix[row][i] = *(data+i);
    }
}

void Matrix::print()
{
    int n = matrix.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            cout << matrix[i][j] << "\t\t";
        }
        cout << "\n";
    }
    cout << endl;
}


// Gaussian Elimination
vector<double> Matrix::Solve(vector<double> B)
{
    //vector<vector<double>> mat = matrix;
    int n = matrix.size();
    
    // rewrite the system as a triangular matrix mat (gaussian elimination) and a corresponding modified matrix B
    for (int i=0; i<n; i++)
    {
        // search for maximum in the column
        double maxCol = abs(matrix[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++)
        {
            if (abs(matrix[k][i]) > maxCol)
            {
                maxCol = abs(matrix[k][i]);
                maxRow = k;
            }
        }
        
        // swap maximum row with current row
        // swap matrix (column by column)
        for (int k=i; k<n;k++)
        {
            double tmp = matrix[maxRow][k];
            matrix[maxRow][k] = matrix[i][k];
            matrix[i][k] = tmp;
        }
        // same swap in B
        double tmp = B[maxRow];
        B[maxRow] = B[i];
        B[i] = tmp;
        
        // make all rows below the max row 0 in current column
        for (int k=i+1; k<n; k++)
        {
            double c = -matrix[k][i]/matrix[i][i];
            for (int j=i; j<n; j++)
            {
                if (i==j)
                {
                    matrix[k][j] = 0;
                }
                else
                {
                    matrix[k][j] += c * matrix[i][j];
                }
            }
            // for element in B
            if (i==n)
            {
                B[k] = 0;
            }
            else
            {
                B[k] += c * B[i];
            }
            
        }
    }  // end of for loop
    
    //this->print();
    // now mat is an upper triangular matrix
    // we can solve easily an upper triangular matrix mat*x = B
    // Example:
    //[ a11, a12]   [x1]   [B1]
    //            *      =
    //[ 0  , a22]   [x2]   [B2]
    // -> we get directly that x2 = B2/a22
    // -> then x1 = B1/a11 - a12*x2/a11
    // So we always go from bottom to top
    
    vector<double> x(n);
    for (int i=n-1; i>=0; i--)
    {
        // initial set just like x2 = B2/a22 (diagonal)
        x[i] = B[i]/matrix[i][i];
        // modify mat to directly obtain x1 = B1/a11 - a12*x2/a11 during next diagonal operation
        for (int k=i-1;k>=0; k--)
        {
            B[k] -= matrix[k][i] * x[i];
        }
    }
    
    return x;
}

// for testing the class
void testMatrix()
{
    cout << "Create and print matrix" << endl;
    vector<double> seed {1.0, 2.0, 3.0, 4.0};
    Matrix mat = Matrix(seed);
    mat.print();
    
    cout << "Create and print matrix" << endl;
    vector<double> seed2 {1.0, 0.0, 0.0, 0.0};
    Matrix mat2 = Matrix(seed2);
    mat2.print();
    vector<double> sol = mat2.Solve(seed);
    
    for (int i=0; i<4; i++)
    {
        cout << sol[i] << endl;
    }
}
