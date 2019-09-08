# Fast Autoregressive Modeling for Time Series in C++

#### Ester Hlav, 2019

The goal of this project is to analyze time series data with an Autoregressive (AR) model. To fit such a model, we construct a set of classes that store data as a dataframe/matrix; for example, one class method was developed to perform Gaussian Elimination on matrices. The fitting of AR models can be computationally expensive as they require the solving of multiple linear systems. Time series analysis is commonly performed in R or Python, which include high-level libraries for experimenting with time series models such as AR, ARMA, ARIMA, GARCH, etc. As for the power and speed of C++, we demonstrate how fast an AR model can perform in a low level language such as C++.
