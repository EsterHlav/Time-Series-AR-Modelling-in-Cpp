#include "ARmodel.h"
using namespace std;

// class Data

Data::Data()    // constructor
{
    array = new double;
    string name = "";
    sizeArray = 0;
}

Data::Data(const Data& a)  // copy constructor
{
    sizeArray = a.sizeArray;
    name = a.name;
    array = new double[sizeArray];
    for (int i=0; i<sizeArray; i++)
    {
        *(array+i) = *(a.array+i);
    }
    // compute autocorrelations to have it available
    this->setAutocorrelations();
}

Data::Data(vector<double> vec, string name_input)   // constructor from vector
{
    sizeArray = vec.size();
    name = name_input;
    array = new double[sizeArray];
    for (int i=0; i<sizeArray; i++)
    {
        *(array+i) = vec.at(i);
    }
    // compute autocorrelations to have it available
    this->setAutocorrelations();
}

Data::Data(double* d, int size, string name_input)   // constructor from array of double
{
    // copy d in a new array array
    sizeArray = size;
    name = name_input;
    
    double* newarray = new double[sizeArray];
    for (int i=0; i<sizeArray; i++)
    {
        *(newarray+i) = *(d+i);
    }
    array = newarray;
    // compute autocorrelations to have it available
    this->setAutocorrelations();
}

Data::~Data()
{
    delete[] array;
}

// methods
void Data::copyToArray(double* newarray, int start, int end, int start_new)
{
    if (end==-1) // convention -1 means end at end of array
    {
        end = sizeArray;
    }
    
    // only if copy makes sense
    for (int i=0; i<end-start; i++)
    {
        *(newarray+start_new+i) = *(array+i+start);
    }
    
}

void Data::addItems(double value, int howMany)
{
    double* newarray = new double[sizeArray+howMany];
    
    // copy old values
    copyToArray(newarray, 0, -1, 0);
    
    // add new values
    for (int i=sizeArray; i<sizeArray+howMany; i++)
    {
        *(newarray+i) = value;
    }
    
    // change attribute of the class
    array = newarray;
    sizeArray += howMany;
}

void Data::removeItems(int howMany, int start)
{
    // if nothing is specified for start
    if (start == -1)
    {
        start = sizeArray-howMany;
    }
    
    // not to allow out of bonds
    if (start+howMany>sizeArray)
    {
        howMany = sizeArray-start+1;
    }
    
    double* newarray = new double[sizeArray-howMany];
    
    // copy beggining of old values
    copyToArray(newarray, 0, start, 0);
    
    // copy end of old values
    copyToArray(newarray, start+howMany, -1, start);
    
    // change attributes of the class
    array = newarray;
    sizeArray -= howMany;
    
}
double Data::at(int index)
{
    // verify it is within the bounds
    if (index>=sizeArray)
    {
        string error = "Out of range for index " + to_string(index) + " array of size " + to_string(sizeArray) +".";
        throw std::out_of_range (error);
    }
    else{
        return *(array+index);
    }
}

void Data::put(int index, double data)
{
    // verify it is within the bounds
    if (index>=sizeArray)
    {
        string error = "Out of range for index " + to_string(index) + " array of size " + to_string(sizeArray) +".";
        throw std::out_of_range (error);
    }
    
    *(array+index) = data;
}

void Data::put(int index, double *data, int howMany)
{
    // verify it is within the bounds
    if (index>=sizeArray)
    {
        string error = "Out of range for index " + to_string(index) + " array of size " + to_string(sizeArray) +".";
        throw std::out_of_range (error);
    }
    
    // add the values from index
    for (int i=0; i<howMany; i++)
    {
        *(array+index+i) = *(data+i);
    }
    
}

void Data::insert(int index, double data)
{
    // verify it is within the bounds
    if (index>=sizeArray)
    {
        string error = "Out of range for index " + to_string(index) + " array of size " + to_string(sizeArray) +".";
        throw std::out_of_range (error);
    }
    
    double* newarray = new double[sizeArray+1];
    
    // copy beggining of old values
    copyToArray(newarray, 0, index, 0);
    
    // copy end of old values
    copyToArray(newarray, index, -1, index+1);
    
    // add the value at index
    *(newarray+index) = data;
    
    // change attributes of the class
    array = newarray;
    sizeArray++;
}

void Data::insert(int index, double *data, int howMany)
{
    // verify it is within the bounds
    if (index>=sizeArray)
    {
        string error = "Out of range for index " + to_string(index) + " array of size " + to_string(sizeArray) +".";
        throw std::out_of_range (error);
    }
    
    double* newarray = new double[sizeArray+howMany];
    
    // copy beggining of old values
    copyToArray(newarray, 0, index, 0);
    
    // add the values from index
    for (int i=0; i<howMany; i++)
    {
        *(newarray+index+i) = *(data+i);
    }
    
    // copy end of old values
    copyToArray(newarray, index, -1, index+howMany);
    
    // change attributes of the class
    array = newarray;
    sizeArray+=howMany;
    
}

void Data::copyTo(vector<double>& v)
{
    // allocate the size for the vector
    v.resize(sizeArray);
    
    // copy the values to vector
    for (int i=0; i<sizeArray; i++)
    {
        v.at(i) = *(array+i);
    }
}

void Data::print()
{
    cout << "Content of array '" + name + "' of size " + to_string(sizeArray) + ":" << endl;
    for (int i=0; i<sizeArray; i++)
    {
        cout << *(array+i) << endl;
    }
    cout << "--- end of array ---" << endl;
}

int Data::size()
{
    return this->sizeArray;
}

void Data::clear()
{
    array = new double;
    sizeArray = 0;
}

// compute the mean of the array
double Data::mean()
{
    double sum = 0;
    for (int i=0; i<sizeArray; i++)
    {
        sum += *(array+i);
    }
    return sum/static_cast<double>(sizeArray);
}

// compute variance of array
double Data::var()
{
    double mean = this->mean();
    double sum = 0;
    for (int i=0; i<sizeArray; i++)
    {
        sum += pow(*(array+i)-mean, 2);
    }
    return sum/static_cast<double>(sizeArray);
}

// compute the autocorrelation of lag
// cf: https://en.wikipedia.org/wiki/Autocorrelation
// Formula is \hat{R}(k)=\frac{1}{(n-k) \sigma^2} \sum_{t=1}^{n-k} (X_t-\mu)(X_{t+k}-\mu)
// this is then used to compute the PACF function
double Data::autocorrelation(int lag)
{
    // compute \frac{1}{(n-k) \sigma^2}
    double coeff = 1.0/( (static_cast<double>(sizeArray)-static_cast<double>(lag))*this->var() );
    double mean = this->mean();
    double sum = 0;
    
    // compute \sum_{t=1}^{n-k} (X_t-\mu)(X_{t+k}-\mu)
    for (int i=0; i<sizeArray-lag; i++)
    {
        sum += (*(array+i) - mean) * (*(array+i+lag)-mean);
    }
    
    return sum*coeff;
}

void Data::setAutocorrelations()
{
    vector<double> autos(sizeArray, 0);
    for (int i=0; i<sizeArray; i++)
    {
        autos[i] = this->autocorrelation(i);
    }
    autocorrelations = autos;
}

void testData()
{
    vector<double> vec { 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 };
    Data d = Data(vec, "myarray");
    d.print();
    cout << d.mean() << endl;
    cout << d.var() << endl;
    for (int i=0; i<vec.size(); i++)
    {
        cout << d.autocorrelation(i) << endl;
        cout << d.autocorrelations[i] << endl;
    }
}
