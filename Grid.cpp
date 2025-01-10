#include "Grid.h"
#include <iostream>
#include <valarray>


// Costruttore che inizializza la griglia con le dimensioni e la spaziatura
Grid::Grid(int nx, int ny, double dx, double dy, std::string & prob)
        : nx(nx), ny(ny), dx(dx), dy(dy), prob(prob) {
    // Inizializzazione dei vettori per densità, velocità e pressione
    rho_.resize(nx * ny, 0.0);
    u_.resize(nx * ny, 0.0);
    v_.resize(nx * ny, 0.0);
    p_.resize(nx * ny, 0.0);
}

// Metodo per inizializzare i valori iniziali (qui puoi aggiungere logica per valori specifici)
void Grid::initialize() {
    double x = 0;
    // Esempio di inizializzazione: settiamo la densità a 1, e velocità a 0
    if(prob=="Sod1"){
        for (int i = 0; i < nx; ++i) {
            double y = 0;
            x += dx;
            for (int j = 0; j < ny; ++j) {
                y += dy;
                if (x < 0.5) {
                    rho(i, j) = 1.0;
                    p(i, j) = 1.0;
                } else {
                    rho(i, j) = 0.125;
                    p(i, j) = 0.1;
                }
                u(i, j) = 0.0;
                v(i, j) = 0.0;
            }
        }
    } else if (prob=="RTI") {
        gamma = 1.666667;
        for (int i = 0; i < nx; ++i) {
            double y = 0;
            x += dx;
            for (int j = 0; j < ny; ++j) {
                y += dy;
                if (y > 0.5 + 0.1*sin(4*M_PI*(x-dx/2))) {
                    rho(i, j) = 2.0;
                    p(i, j) = 2 * ( 1 - y );
                } else {
                    rho(i, j) = 1.0;
                    p(i, j) = 1.5 - y;
                }
                u(i, j) = 0.0;
                v(i, j) = 0; // -0.1 * sqrt(gamma * p(i, j) / rho(i, j)) * cos(8 * M_PI * (i * dx)); // Perturbazion*/
            }
        }
    } else if (prob=="Case3") {
        for (int i = 0; i < nx; ++i) {
            double y = 0;
            x += dx;
            for (int j = 0; j < ny; ++j) {
                y += dy;
                if (y <= 0.5 && x <= 0.5) {
                    p(i, j) = 0.029;
                    rho(i, j) = 0.138;
                    u(i, j) = 1.206;
                    v(i, j) = 1.206;
                } else if (y > 0.5 && x <= 0.5) {
                    p(i, j) = 0.3;
                    rho(i, j) = 0.5323;
                    u(i, j) = 1.206;
                    v(i, j) = 0;
                } else if (y > 0.5 && x > 0.5) {
                    p(i, j) = 1.5;
                    rho(i, j) = 1.5;
                    u(i, j) = 0;
                    v(i, j) = 0;
                } else {
                    p(i, j) = 0.3;
                    rho(i, j) = 0.5323;
                    u(i, j) = 0;
                    v(i, j) = 1.206;
                }
            }
        }
    } else if (prob=="Case12") {
        for (int i = 0; i < nx; ++i) {
            double y = 0;
            x += dx;
            for (int j = 0; j < ny; ++j) {
                y += dy;
                if (y <= 0.5 && x <= 0.5) {
                    p(i, j) = 1.0;
                    rho(i, j) = 0.8;
                    u(i, j) = 0.0;
                    v(i, j) = 0.0;
                } else if (y > 0.5 && x <= 0.5) {
                    p(i, j) = 1.0;
                    rho(i, j) = 1.0;
                    u(i, j) = 0.7276;
                    v(i, j) = 0.0;
                } else if (y > 0.5 && x > 0.5) {
                    p(i, j) = 0.4;
                    rho(i, j) = 0.5313;
                    u(i, j) = 0;
                    v(i, j) = 0;
                } else {
                    p(i, j) = 1.0;
                    rho(i, j) = 1.0;
                    u(i, j) = 0;
                    v(i, j) = 0.7276;
                }
            }
        }
    }  else if (prob=="Case7") {
        for (int i = 0; i < nx; ++i) {
            double y = 0;
            x += dx;
            for (int j = 0; j < ny; ++j) {
                y += dy;
                if (y <= 0.5 && x <= 0.5) {
                    p(i, j) = 1.0;
                    rho(i, j) = 1.0;
                    u(i, j) = -0.75;
                    v(i, j) = 0.5;
                } else if (y > 0.5 && x <= 0.5) {
                    p(i, j) = 1.0;
                    rho(i, j) = 2.0;
                    u(i, j) = 0.75;
                    v(i, j) = 0.5;
                } else if (y > 0.5 && x > 0.5) {
                    p(i, j) = 1.0;
                    rho(i, j) = 1.0;
                    u(i, j) = 0.75;
                    v(i, j) = -0.5;
                } else {
                    p(i, j) = 1.0;
                    rho(i, j) = 3.0;
                    u(i, j) = -0.75;
                    v(i, j) = -0.5;
                }
            }
        }
    }
    applyBound();
    intConserved = integrate();
    std::cout << "Initial conserved: " << std::endl;
    print(intConserved);
}

double Grid::getDx() const {
    return dx;
}

double Grid::getDy() const {
    return dy;
}

std::vector<double> Grid::asCoords(int i, int j) {
    std::vector<double> coords(4);
    coords[0] = rho(i,j);
    coords[1] = rho(i,j)*u(i,j);
    coords[2] = rho(i,j)*v(i,j);
    coords[3] = E(i,j);
    return coords;
}

double Grid::E(int i, int j) {
    return p(i,j)/(gamma-1) + 0.5*rho(i,j)*(u(i,j)*u(i,j) + v(i,j)*v(i,j));
}

double Grid::H(int i, int j) {
    return (E(i,j) + p(i,j))/rho(i,j);
}

void Grid::updateRho(std::vector<double> &rhonew) {
    rho_ = rhonew;
}

void Grid::updateU(std::vector<double> & unew) {
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            u_[i*ny + j] = (1/rho(i,j))*unew[i*ny + j];
        }
    }
}

void Grid::updateV(std::vector<double> &vnew) {
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            v_[i*ny + j] = (1/rho(i,j))*vnew[i*ny + j];
        }
    }
}

void Grid::updateEP(std::vector<double> &Enew) {
    E_ = Enew;
    double tempP;
    for(int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            tempP = (gamma-1)*(E_[i*ny +  j] - 0.5*rho(i,j)*( u(i,j) * u(i,j) + v(i,j) * v(i,j) ));
            p_[i*ny + j] = max(0,tempP); //max(tempP,tempP);
        }
    }
}

void Grid::update(std::vector<std::vector<double>> &unew) {
    updateRho(unew[0]);
    updateU(unew[1]);
    updateV(unew[2]);
    updateEP(unew[3]);
}

double &Grid::rho(int i, int j) {
    if (i <= -1 || i >= getNx() || j <= -1 || j >= getNy()) {
        return ghostCell(i, j, 0);
    } else {
        return rho_[i * ny + j];
    }
}

double &Grid::u(int i, int j) {
    if (i <= -1 || i >= getNx() || j <= -1 || j >= getNy()) {
        return ghostCell(i, j, 1);
    } else {
        return u_[i * ny + j];
    }
}

double &Grid::v(int i, int j) {
    if (i <= -1 || i >= getNx() || j <= -1 || j >= getNy()) {
        return ghostCell(i, j, 2);
    } else {
        return v_[i * ny + j];
    }
}

double &Grid::p(int i, int j) {
    if (i <= -1 || i >= getNx() || j <= -1 || j >= getNy()) {
        return ghostCell(i, j, 3);
    } else {
        return p_[i * ny + j];
    }
}

double &Grid::ghostCell(int i, int j, int var) {
    if(i <= -1){
        return boundLeft[j][var];
    } else if (i >= getNx()) {
        return boundRight[j][var];
    } else if (j <= -1) {
        if(prob=="RTI"){
            return boundBottom[0][var];
        }
        return boundBottom[i][var];
    } else if (j >= getNy()) {
        if(prob=="RTI"){
            return boundTop[0][var];
        }
        return boundTop[i][var];
    } else {
        std::cerr << "Chiamata errata Ghost Cell" << std::endl;
        return boundLeft[100][var];
    }
}

void Grid::applyBound() {
    //Dirichlet Top
    //boundTop.resize(1,std::vector<double> (4,0));
    //boundTop[0] = {.125,0,0,.1};
    if (prob=="Sod1"){
        boundBottom.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundBottom[i] = {rho(i, 0), u(i, 0), v(i, 0), p(i, 0)};
        }
        //Reflective left
        boundLeft.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundLeft[i] = {rho(0, i), u(0, i), v(0, i), p(0, i)};
        }
        boundTop.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundTop[i] = {rho(i, ny-1), u(i, ny-1), v(i, ny-1), p(i, ny-1)};
        }

        boundRight.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundRight[i] = {rho(nx-1, i), u(nx-1, i), v(nx-1, i), p(nx-1, i)};
        }
    } else if (prob=="RTI"){
        //Reflective left
        boundLeft.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundLeft[i] = {rho(0, i), -u(0, i), v(0, i), p(0, i)};
        }
        /*
        boundBottom.resize(1,std::vector<double> (4,0));
        boundBottom[0] = {2.0,0,0,1.5};
        */

        boundBottom.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundBottom[i] = {rho(i, 0), u(i, 0), -v(i, 0), p(i, 0)};
        }
        /*
        boundTop.resize(1,std::vector<double> (4,0));
        boundTop[0] = {1.0,0,0,0.};
        */

        boundTop.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundTop[i] = {rho(i, ny-1), u(i, ny-1), -v(i, ny-1), p(i, ny-1)};
        }

        boundRight.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundRight[i] = {rho(nx-1, i), -u(nx-1, i), v(nx-1, i), p(nx-1, i)};
        }
    } else if (prob=="Case3" || prob=="Case12" || prob=="Case16"){
        boundBottom.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundBottom[i] = {rho(i, 0), u(i, 0), v(i, 0), p(i, 0)};
        }
        //Reflective left
        boundLeft.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundLeft[i] = {rho(0, i), u(0, i), v(0, i), p(0, i)};
        }
        boundTop.resize(nx, std::vector<double>(4, 0.0));
        for (int i = 0; i < nx; i++) {
            boundTop[i] = {rho(i, ny-1), u(i, ny-1), v(i, ny-1), p(i, ny-1)};
        }

        boundRight.resize(ny, std::vector<double>(4, 0.0));
        for (int i = 0; i < ny; i++) {
            boundRight[i] = {rho(nx-1, i), u(nx-1, i), v(nx-1, i), p(nx-1, i)};
        }
    }

}

std::vector<double> Grid::integrate() {
    std::vector<double> sums(4,0);
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            sums[0] += dx*dy*rho(i,j);
            sums[1] += dx*dy*rho(i,j)*u(i,j);
            sums[2] += dx*dy*rho(i,j)*v(i,j);
            sums[3] += dx*dy*E(i,j);
        }
    }
    return sums;
}

const std::vector<double> &Grid::getIntConserved() const {
    return intConserved;
}

double Grid::sMax() {
    double max = 0;
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            if(std::abs(u(i,j)) + sqrt(gamma*p(i,j)/rho(i,j)) > max) {
                max = std::abs(u(i,j)) + sqrt(gamma*p(i,j)/rho(i,j));
            }
        }
    }
    return max;
}

const std::string &Grid::getProb() const {
    return prob;
}

double Grid::pCoords(std::vector<double> &coords) {
    return max(0, (gamma-1)*(coords[3] - 0.5*(coords[1]*coords[1] + coords[2]*coords[2])/coords[0]));
}

double Grid::getGamma() const {
    return gamma;
}











