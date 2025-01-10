//
// Created by UTENTE on 09/12/2024.
//

#include <valarray>
#include "Flux.h"
#include "Time.h"
#include "Print.h"
#include "Grid1D.h"

#ifndef PROJET_TEST_H
#define PROJET_TEST_H

#endif //PROJET_TEST_H

void test1(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test2(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test3(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1DWENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test4(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    double T = 0.25;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test5(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    double T = 0.2;
    Grid1D grid(nx, dx, {{1,0.75,1},{0.125,0,0.1}});
    grid.initialize("Sod2"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test6(){
    int nx = 100; // Numero di celle in x
    double dx = 1.0/nx; // Spaziatura in x
    double dt = 0.001; // Passo temporale
    double T = 0.15;
    Grid1D grid(nx, dx, {{1,-2,0.4},{1,2,0.4}});
    grid.initialize("123"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<10; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.advance1D(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    for (int i = 10; i<T/dt; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles1D(grid);
    //grid.applyBoundaryConditions();
}

void test7(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit = 300;
    std::string prob = "Case3";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test8(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.001; // Passo temporale
    int nit =  250;
    std::string prob = "Case12";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
        //print.exportToVTK(grid, "vtk/results" + std::to_string(i) + ".vtk");
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test9(){
    int nx = 20; // Numero di celle in x
    int ny = 20; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.00125; // Passo temporale
    int nit = 100;
    std::string prob = "Case3";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,ny,nx,dt);
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        grid = time.advance(flux,grid);
        time.isConserved(grid.integrate(),residue0);
    }
    Print print;
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test10(){
    int nx = 64; // Numero di celle in x
    int ny = 8*nx; // Numero di celle in y
    double dx = 0.25/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.0001; // Passo temporale
    int nit = 20000;
    std::string prob = "RTI";
    Grid grid(nx, ny, dx, dy,prob);
    grid.initialize(); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    Time time(dt,residue0);
    Flux flux(dx,dy,nx,ny,dt);
    Print print;
    for (int i = 0; i<nit; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);

        if(std::fmod(i,500)==0){
            print.exportToVTK(grid, "vtk/results" + std::to_string(i/500) + ".vtk");
        }
    }
    print.printGridToFiles(grid);
    print.exportToVTK(grid, "results.vtk");
    std::cout<< "Final T: " << nit*dt << std::endl;
}

void test11(){
    int nx = 200; // Numero di celle in x
    int ny = 1; // Numero di celle in y
    double dx = 1.0/nx; // Spaziatura in x
    double dy = 1.0/ny; // Spaziatura in y
    double dt = 0.002; // Passo temporale
    int nit = 125;
    Grid1D grid(nx, dx, {{1,0,1},{0.125,0,0.1}});
    grid.initialize("Sod1"); // Configura i valori iniziali
    std::vector<double> residue0 = grid.getIntConserved();
    std::pair<std::vector<double>,std::vector<double>> us;
    Time time(dt,residue0);
    Flux flux(dx,nx,dt);
    for (int i = 0; i<nit/*T/dt*/; i++) {
        std::cout << "Iteration: " << i+1 << std::endl;
        grid = time.evolveRK3WENO(grid,flux);
        time.isConserved(grid.integrate(),residue0);
    }
    std::cout << "Final t: " << nit*dt <<std::endl;
    //us = flux.WENO1D(grid,5);
    //print(us.first);
    //print(us.second);
    Print print;
    print.printGridToFiles1D(grid);
}
