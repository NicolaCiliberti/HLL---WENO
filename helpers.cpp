//
// Created by UTENTE on 19/12/2024.
//
#include <valarray>
#include "helpers.h"

std::vector<double> operator * (std::vector<double> vect, double mu){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect[i]*mu;
    }
    return out;
}

std::vector<double> operator * (double mu, std::vector<double> vect){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect[i]*mu;
    }
    return out;
}

std::vector<double> operator + (std::vector<double> vect1, std::vector<double> vect){
    if(vect1.size() == vect.size()){
        std::vector<double> out(vect.size());
        for(int i = 0; i<vect.size(); i++){
            out[i] = vect[i]+vect1[i];
        }
        return out;
    } else {
        std::vector<double> ciao;
        std::cerr << "Errore: qualcosa è andato storto!" << std::endl;
        return ciao;
    }

}

std::vector<double> operator - (std::vector<double> vect1, std::vector<double> vect){
    std::vector<double> out(vect.size());
    for(int i = 0; i<vect.size(); i++){
        out[i] = vect1[i]-vect[i];
    }
    return out;
}

void print(std::vector<double> vect){
    for(int i=0; i<vect.size(); i++){
        std::cout << vect[i] << "    ";
    };
    std::cout << std::endl;
}

double max(double a, double b) {
    if(a>b){
        return a;
    } else {
        return b;
    }
}

double min(double a, double b) {
    if(a<b){
        return a;
    } else {
        return b;
    }
}

std::vector<double> operator / (std::vector<double> vect, std::vector<double> vect2){
    std::vector<double> out (vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = vect[i]/vect2[i];
    }
}



std::vector<double> operator ^ (std::vector<double> vect, int n){
    std::vector<double> out(vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = pow(vect[i],n);
    }
    return out;
}

std::vector<double> operator * (std::vector<double> vect, std::vector<double> vect2){
    std::vector<double> out (vect.size());
    for(int i=0; i<vect.size(); i++){
        out[i] = vect[i]*vect2[i];
    }
}

double superBee(double r){
    return max(0,r);
}

std::vector<double> superBee(std::vector<double> & r){
    std::vector<double> out(r.size());
    for(int i=0; i<r.size(); i++){
        out[i] = superBee(r[i]);
    }
    return out;
}

double polynomial(double x){
        return /*x*x*x*x*x + x*x*x*x +*/  x*x*x + x * x  + 2 * x + 1;
}
