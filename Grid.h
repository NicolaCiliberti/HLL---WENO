// grid.h
#ifndef GRID_H
#define GRID_H

#include <vector>
#include <iostream>
#include "helpers.h"

class Grid {
public:
    const std::string &getProb() const;

public:
    // Costruttore
    Grid(int nx, int ny, double dx, double dy, std::string &prob);

    // Inizializza i valori iniziali
    void initialize();

    // Getter per le dimensioni
    int getNx() const { return nx; }
    int getNy() const { return ny; }

    double getDx() const;

    double getDy() const;

    // Applica boundary conditions
    void applyBoundaryConditions();

    // Accesso ai dati
    double& rho(int i, int j);
    double& u(int i, int j);
    double& v(int i, int j);
    double& p(int i, int j);
    double E(int i, int j);
    double H(int i, int j);
    std::vector<double> asCoords(int i, int j);
    void updateRho(std::vector<double> & rhonew);
    void updateU(std::vector<double> & unew);
    void updateV(std::vector<double> & vnew);
    void updateEP(std::vector<double> & Enew);
    void update(std::vector<std::vector<double>> & unew);
    double& ghostCell(int i, int j, int var);
    double getGamma() const;
    void applyBound();
    void applyBoundCond(std::vector<std::vector<double>> & boundVect, int cond, int side);
    std::vector<double> integrate();
    double pCoords(std::vector<double> & coords);
    double sMax();
private:
    int nx, ny; // Numero di celle
    double dx, dy; // Spaziatura del reticolo

    // Variabili fisiche
    std::vector<double> rho_; // Densità
    std::vector<double> u_;   // Velocità x
    std::vector<double> v_;   // Velocità y
    std::vector<double> p_;   // Pressione
    std::vector<double> E_;
    std::vector<std::vector<double>> boundTop;
    std::vector<std::vector<double>> boundBottom;
    std::vector<std::vector<double>> boundLeft;
    std::vector<std::vector<double>> boundRight;
    std::vector<double> intConserved;
    std::vector<int> bound; //0 Outflow, 1 Reflective, 2 Periodic, 3 Dirichlet 0
    std::string prob;
    double gamma = 1.4;
public:
    const std::vector<double> &getIntConserved() const;
};


#endif // GRID_H
