//
// Created by UTENTE on 19/12/2024.
//

#include <vector>
#include <iostream>

#ifndef PROJET_HELPERS_H
#define PROJET_HELPERS_H

#endif //PROJET_HELPERS_H

std::vector<double> operator * (std::vector<double> vect, double mu);
std::vector<double> operator ^ (std::vector<double> vect, int n);
std::vector<double> operator / (std::vector<double> vect, std::vector<double> vect2);
std::vector<double> operator * (double mu, std::vector<double> vect);
std::vector<double> operator + (std::vector<double> vect1, std::vector<double> vect);
std::vector<double> operator - (std::vector<double> vect1, std::vector<double> vect);
std::vector<double> operator + (std::vector<double> vect1, std::vector<double> vect);
double superBee(double r);
std::vector<double> superBee(std::vector<double> & r);
void print(std::vector<double> vect);
double max(double a, double b);
double min(double a, double b);
double polynomial(double x);
