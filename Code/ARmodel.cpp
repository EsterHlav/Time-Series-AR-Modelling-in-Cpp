#include "ARmodel.h"
using namespace std;

// constructor
ARmodel::ARmodel() {}

ARmodel::ARmodel(Data data_input)
{
    // save data
    data = data_input;
}
ARmodel::ARmodel(Data data_input, int order_input)
{
    // save data
    data = data_input;
    
    // set input order
    order = order_input;
}

ARmodel::ARmodel(Data data_input, int order_input, double phis_input[])
{
    // save data
    data = data_input;
    
    // set input order
    order = order_input;
    
    // copy phis
    double* phis = new double[order];
    for (int i=0; i<order; i++)
    {
        *(phis+i) = *(phis_input+i);
    }
    
}

// copy constructor
ARmodel::ARmodel(ARmodel& model)
{
    // save data
    data = model.data;
    
    // set input order
    order = model.order;
    
    // copy phis
    double* phis = new double[order];
    for (int i=0; i<order; i++)
    {
        *(phis+i) = *(model.phis+i);
    }
}

// destructor
ARmodel::~ARmodel()
{
    delete &data;
    delete[] phis;
}

//methods
void ARmodel::setData(Data data_input)
{
    data = data_input;
}

// determine the coefficient of the autocorrelation function
// c.f.: https://en.wikipedia.org/wiki/Autocorrelation_function
double* ARmodel::ACF()
{
    // create array of autocorrelation coefficient
    double* autocorr = new double[data.sizeArray];
    for (int i=0; i<data.sizeArray; i++)
    {
        *(autocorr+i) = data.autocorrelation(i);
    }
    
    return autocorr;
}

// PACF(k) requires to solve a linear system
// c.f.: https://stats.stackexchange.com/questions/129052/acf-and-pacf-formula/129374
double ARmodel::PACF(int k)
{
    // by definition is 1 in k=0
    if (k==0)
    {
        return 1.0;
    }
    
    // create sub autocorrelation vector [autocorr[0], autocorr[1], ..., autocorr[k-1]]
    vector<double>::const_iterator first = data.autocorrelations.begin();
    vector<double>::const_iterator last = data.autocorrelations.begin()+k;
    vector<double> veck(first, last);
    
    //    cout << "Vector veck:" << endl;
    //    cout << veck.size() << endl;
    //    for(auto &i: veck)
    //        cout << i << endl;
    
    // create the autocorrelation matrix system of order k from order k vector autocorrelation
    Matrix matk = Matrix(veck);
    //matk.print();
    
    // create sub autocorrelation vector [autocorr[1], autocorr[2], ..., autocorr[k]]
    vector<double>::const_iterator first2 = data.autocorrelations.begin()+1;
    vector<double>::const_iterator last2 = data.autocorrelations.begin()+k+1;
    vector<double> veck2(first2, last2);
    
    // solve the system matk * x = veck+1 for x
    vector<double> sol = matk.Solve(veck2);
    
    // cout << "Solution for k = " + to_string(k) + " is " << sol[k-1]  << endl ;
    
    // the PACF(k) is the last element in that sol vector (index k-1 for vector of size k)
    return sol[k-1];
    
}
// PACF requires to solve a linear system for each k=1..size
// c.f.: https://stats.stackexchange.com/questions/129052/acf-and-pacf-formula/129374
double* ARmodel::PACF()
{
    double* pacf = new double[data.sizeArray];
    
    // cut-off after 50 lags
    int maxLags;
    if (data.sizeArray>50)
    {
        maxLags = 50;
    }
    else
    {
        maxLags = data.sizeArray;
    }
    
    // loop over all lags
    for (int i=0; i<maxLags; i++)
    {
        *(pacf+i) = this->PACF(i);
    }
    
    return pacf;
}


// determine the number of order that are significant
// c.f.: https://en.wikipedia.org/wiki/Partial_autocorrelation_function
void ARmodel::determineOrder()
{
    double* pacf = this->PACF();
    
    // compute the limit for the statistical significance at 95%
    // formula is 0 \pm 1.96/{\sqrt{n}}
    double limit = 1.96/sqrt(data.sizeArray);
    
    cout << "Limit significance: |gamma|>" << limit << endl;
    
    // cut-off after 50 lags
    int maxLags;
    if (data.sizeArray>50)
    {
        maxLags = 50;
    }
    else
    {
        maxLags = data.sizeArray;
    }
    
    bool* significance = new bool[maxLags];
    // create array of significance (boolean)
    for (int i=0; i<maxLags; i++)
    {
        *(significance+i) = (*(pacf+i)>limit) || (*(pacf+i)<(-limit));
    }
    
    // print significance
    for (int i=0; i<maxLags; i++)
    {
        cout << "Significance " << i << " = " << *(significance+i) << " (" << *(pacf+i) << ")." << endl;
    }
    
    // determine the coefficient by looking at the length of the sequence of initial true
    // for example: [true, true, true, false, false, true] should give 3
    int count = 0;
    while (*(significance+count) && count<data.sizeArray) {
        count++;
    };
    
    // create a warning when higher lags are significant
    // for example: [true, true, true, false, false, true] should give a warning
    // for example: [true, true, true, false, false, false] should NOT give a warning
    
    for (int i=count; i<maxLags; i++)
    {
        if (*(significance+i))
        {
            cout << "WARNING: Lag " + to_string(i) + " is significant and order is " + to_string(count) + "! Maybe an AR model is not suitable." << endl;
        }
    }
    
    // if nothing is significant, then order by default is 1.
    if (count==0)
    {
        // save the default order
        order = 1;
    }
    else
    {
        // save the order
        order = count-1;
    }
    
    
    // save the order
}
void ARmodel::fit()
{
    // determine the order
    this->determineOrder();
    
    cout << "Order of model was determined to be " << this->order << endl;
    
    // fit the model
    this->fit(order);
    
}

// the coefficient are determined using the Yule Walker equations
// c.f.: https://en.wikipedia.org/wiki/Autoregressive_model#Yule%E2%80%93Walker_equations and http://paulbourke.net/miscellaneous/ar/
// AutocorrelationMatrix * \Phi = AutocorrelationVector, solve for \Phi
void ARmodel::fit(int order2fit)
{
    order = order2fit;
    // create the matrix AutocorrelationMatrix //
    
    // create sub autocorrelation vector [autocorr[0], autocorr[1], ..., autocorr[order-1]]
    vector<double>::const_iterator first = data.autocorrelations.begin();
    vector<double>::const_iterator last = data.autocorrelations.begin()+order;
    vector<double> vecOrder(first, last);
    
    // create the autocorrelation matrix system of order k
    Matrix matOrder = Matrix(vecOrder);
    cout << "Matrix to be inversed or solved:" << endl;
    matOrder.print();
    
    // create sub autocorrelation vector [autocorr[1], autocorr[2], ..., autocorr[order]]
    vector<double>::const_iterator first2 = data.autocorrelations.begin()+1;
    vector<double>::const_iterator last2 = data.autocorrelations.begin()+order+1;
    vector<double> veckOrder2(first2, last2);
    
    // solve the system matOrder * x = veckOrder2 for x
    vector<double> sol = matOrder.Solve(veckOrder2);
    
    cout << "Size of solution: " << sol.size() << endl;
    
    // fill-in the phis
    double* phisnew = new double[order];
    
    for (int i=0; i<order; i++)
    {
        *(phisnew+i) = sol[i];
    }
    phis = phisnew;
}

void ARmodel::print()
{
    cout << "Order of model: " << to_string(order) << endl;
    for (int i=0; i<order; i++)
    {
        cout << "Coeff phi_" + to_string(i+1) << " = " << *(phis+i) << endl;
    }
    cout << "End of model" << endl;
}

// predict on the training set
Data ARmodel::predict()
{
    Data Predictions = Data();
    
    // computing X of N+j, i.e. all Xs given by N
    for (int t = order; t < data.sizeArray; t++)
    {
        double sum = 0;
        
        // here computing Xhat of j
        // sum(phi_i * X_t-i for i: 0->order)
        for (int i=0; i < order; i++)
        {
            sum += ( *(phis+i) * data.at(t-i-1) );
        }
        
        Predictions.addItems(sum);
    }
    
    return Predictions;
}
// predict out of sample (predict after the end of the provided data)
Data ARmodel::predictOOS(const int N)
{
    Data Predictions = Data();
    
    // add first values to be able to apply recurrence relationship
    for (int t = data.sizeArray-order; t < data.sizeArray; t++)
    {
        Predictions.addItems(data.at(t));
    }
    
    for (int t = data.sizeArray; t < data.sizeArray+N; t++)
    {
        double sum = 0;
        
        // here computing Xhat of j for j>data.sizeArray
        // sum(phi_i * X_t-i for i: 0->order)
        for (int i=0; i < order; i++)
        {
            sum += ( *(phis+i) * Predictions.at(t-i-1) );
        }
        
        Predictions.addItems(sum);
    }
    
    // remove the first order items since there were not predicted values
    Predictions.removeItems(order, 0);
    
    return Predictions;
}
double ARmodel::test()
{
    // predict on training data
    Data predictions = this->predict();
    
    // sum up errors to compute MSE
    double errors = 0;
    
    for (int i=order; i<data.sizeArray; i++)
    {
        errors += pow(predictions.at(i-order)-data.at(i), 2);
    }
    
    return errors/static_cast<double>(data.sizeArray-order);
    
}

// given a dataset, compute the MSE
double ARmodel::test(Data inputs)
{
    // predict
    
    Data Predictions = Data();
    
    // computing X of N+j, i.e. all Xs given by N
    for (int t = order; t < inputs.sizeArray; t++)
    {
        double sum = 0;
        
        // here computing Xhat of j
        // sum(phi_i * X_t-i for i: 0->order)
        for (int i=0; i < order; i++)
        {
            sum += ( *(phis+i) * inputs.at(t-i-1) );
        }
        
        Predictions.addItems(sum);
    }
    
    // sum up errors to compute MSE
    double errors = 0;
    
    for (int i=order; i<inputs.sizeArray; i++)
    {
        errors += pow(Predictions.at(i-order)-inputs.at(i), 2);
    }
    
    return errors/static_cast<double>(inputs.sizeArray-order);
}

void testAR2()
{
    // Example 1
    // linear case with \phi=2 -->   $X_t = 1/2 * X_t-1,
    double phis_ex1[] = {0.3, 0.7};
    double x[] = {1.0, 1.0};
    const int ORDER = 2; // AR(1) --> autoregressive model of order 1
    const int N = 500; // generating 3 Xs, X_1, X_2, X_3
    const double SIGMA = 0.01; // normal distribtuion with standard deviation of 0.001
    
    double* gen = generate(phis_ex1, x, ORDER, N, SIGMA);
    
    Data d = Data(gen, N, "AR(2) train");
    d.print();
    
    ARmodel* model = new ARmodel(d);
    
    double* acf = model->ACF();
    for (int i=0; i<d.sizeArray; i++)
    {
        cout << "ACF(" << i << "):" << endl;
        cout << *(acf+i) << endl;
    }
    
    double* pacf = model->PACF();
    for (int i=0; i<d.sizeArray; i++)
    {
        cout << "PACF(" << i << "):" << endl;
        cout << *(pacf+i) << endl;
    }
    
    // estimate the best order and the coefficients for that order
    model->fit();
    model->print();
    
    double MSEtrain = model->test();
    cout << "MSE training: " << MSEtrain << endl;
    
    // create test set with same coefficients and small noise
    double x2[] = {4.0, 2.0};
    const int N2 = 30;
    const double SIGMA2 = 0.1; // normal distribtuion with standard deviation of 0.1
    
    double* gen_test = generate(phis_ex1, x2, ORDER, N2, SIGMA2);
    Data d2 = Data(gen_test, N2, "AR(2) test1");
    // d2.print();
    
    double MSEtest1 = model->test(d2);
    cout << "MSE test AR(2) 1: " << MSEtest1 << endl;
    
    // create test set with same coefficients and small noise
    double x3[] = {3.0, 1.5};
    const int N3 = 30;
    const double SIGMA3 = 1; // normal distribtuion with standard deviation of 1
    
    double* gen_test2 = generate(phis_ex1, x3, ORDER, N3, SIGMA3);
    Data d3 = Data(gen_test2, N3, "AR(2) test2");
    //d3.print();
    
    double MSEtest2 = model->test(d3);
    cout << "MSE test AR(2) 2: " << MSEtest2 << endl;
    
    // test case with AR(5) process data
    double phis_ex4[] = {0.3, 0.7, -0.8, 0.2, -0.1};
    double x4[] = {1.0, 2.1, -1.0, 1.3, 1.4};
    const int ORDER4 = 5; // AR(5) --> autoregressive model of order 5
    const int N4 = 50;
    const double SIGMA4 = 0.5;
    
    double* gen_test3 = generate(phis_ex4, x4, ORDER4, N4, SIGMA4);
    Data d4 = Data(gen_test3, N3, "AR(5) test3");
    //d4.print();
    
    double MSEtest3 = model->test(d4);
    cout << "MSE test AR(5) 3: " << MSEtest3 << endl;
}

void testAR5()
{
    // Example 2
    // linear case with \phi=5 -->   $X_t = 0.3 * X_t-1 + 0.6 * X_t-2 - 0.9 * X_t-3 + 0.3 * X_t-4 - 0.1 * X_t-5,
    double phis_ex1[] = {0.4, 0.2, -0.2, 0.3, -0.1};
    double x[] = {1.0, 2.1, -1.0, 1.3, 1.4};
    const int ORDER = 5; // AR(5) --> autoregressive model of order 5
    const int N = 3000;
    const double SIGMA = 0.1; // normal distribution with standard deviation of 0.1
    
    double* gen = generate(phis_ex1, x, ORDER, N, SIGMA);
    
    Data d = Data(gen, N, "AR(1) train");
    d.print();
    
    ARmodel* model = new ARmodel(d);
    
    double* acf = model->ACF();
    for (int i=0; i<d.sizeArray; i++)
    {
        cout << "ACF(" << i << "):" << endl;
        cout << *(acf+i) << endl;
    }
    
    double* pacf = model->PACF();
    for (int i=0; i<d.sizeArray; i++)
    {
        cout << "PACF(" << i << "):" << endl;
        cout << *(pacf+i) << endl;
    }
    
    // estimate the best order and the coefficients for that order
    model->fit();
    model->print();
    
    double MSEtrain = model->test();
    cout << "MSE training: " << MSEtrain << endl;
    
    // create test set with same coefficients and small noise
    double x2[] = {4.0, 2.0, 1.2, 3.5, 2.5};
    const int N2 = 30;
    const double SIGMA2 = 0.05; // normal distribtuion with standard deviation of 0.1
    
    double* gen_test = generate(phis_ex1, x2, ORDER, N2, SIGMA2);
    Data d2 = Data(gen_test, N2, "AR(5) test1");
    // d2.print();
    
    double MSEtest1 = model->test(d2);
    cout << "MSE test AR(5) 1: " << MSEtest1 << endl;
    
    // create test set with same coefficients and small noise
    double x3[] = {3.0, 1.5, -1.0, 4.3, 0.5};
    const int N3 = 30;
    const double SIGMA3 = 1; // normal distribtuion with standard deviation of 1
    
    double* gen_test2 = generate(phis_ex1, x3, ORDER, N3, SIGMA3);
    Data d3 = Data(gen_test2, N3, "AR(5) test2");
    //d3.print();
    
    double MSEtest2 = model->test(d3);
    cout << "MSE test AR(5) 2: " << MSEtest2 << endl;
    
    // test case with AR(2) process data
    double phis_ex4[] = {0.3, 0.7};
    double x4[] = {1.0, 2.1};
    const int ORDER4 = 2; // AR(2) --> autoregressive model of order 2
    const int N4 = 50;
    const double SIGMA4 = 0.5;
    
    double* gen_test3 = generate(phis_ex4, x4, ORDER4, N4, SIGMA4);
    Data d4 = Data(gen_test3, N3, "AR(5) test3");
    //d4.print();
    
    double MSEtest3 = model->test(d4);
    cout << "MSE test AR(2) 3: " << MSEtest3 << endl;
}

