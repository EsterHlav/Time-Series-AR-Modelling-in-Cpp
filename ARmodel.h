#ifndef ARmodel_h
#define ARmodel_h

#include <stdio.h>
#include "fstream"
#include "sstream"
#include <stdio.h>
#include <iostream>
#include <list>
#include <vector>
#include <cmath>

// utility for generating data
double generateGaussianNoise(const double, const double);
double* generate(double*, double*, const int, const int, const double);


// class Matrix

class Matrix
{
public:
    std::vector< std::vector<double>> matrix;
    int size;
    
    // constructor from vector<vector>
    Matrix(std::vector<std::vector<double>>);
    
    // constructor for a lagged correlation matrix from a correlation vector
    Matrix(std::vector<double>);
    
    // set a row of the vector
    void setRow(double*, int);
    
    // methods useful for ARmodel
    void print();
    std::vector<double> Solve(std::vector<double>);
    
};

void testMatrix();

// class Data

class Data
{
public:
    double* array;
    std::string name;
    int sizeArray;
    std::vector<double> autocorrelations;
    
    Data();    // constructor
    Data(const Data&);  // copy constructor
    Data(std::vector<double> init, std::string name_input="");   // constructor from vector
    Data(double* d, int size, std::string name_input=""); // constructor from array of double
    ~Data();  // destructor
    
    // methods
    void copyToArray( double* newarray, int start, int end, int start_new);
    void addItems(double value=0.0, int howMany=1);
    void removeItems(int howMany, int start=-1);
    double at(int index);
    void put(int index, double data);
    void put(int index, double* data, int howMany);
    void insert(int index, double data);
    void insert(int index, double* data, int howMany);
    void copyTo(std::vector<double>& v);
    void print();
    int size();
    void clear();
    
    // methods specific time series
    double mean();
    double var();
    double autocorrelation(int);
    void setAutocorrelations();
};

void testData();

// class ARmodel

class ARmodel
{
public:
    Data data;
    double* phis;
    int order;
    
    // constructor
    ARmodel();
    ARmodel(Data);
    ARmodel(Data, int);
    ARmodel(Data, int, double*);
    
    // copy constructor
    ARmodel(ARmodel&);
    
    // destructor
    ~ARmodel();
    
    //methods
    void setData(Data);
    double* ACF();
    double PACF(int k);
    double* PACF();
    void determineOrder();
    void fit();
    void fit(int);
    void print();
    Data predict();
    Data predictOOS(const int);
    double test();
    double test(Data);
    
};

void testAR2();
void testAR5();

#endif /* ARmodel_h */
