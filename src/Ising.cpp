#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include <string>
#include <cmath>

#include "Ising.h"

std::random_device rd;  // a seed source for the random number engine
std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> prob(0.0, 1.0);

Ising2d::Ising2d(int size, float temp, long iter)
: size(size), temp(temp), iter(iter){
    grid = new int*[size];

    // initialize model with random spins
    for(int i = 0; i < size; i++) {

        grid[i] = new int[size];
        for(int j = 0; j < size; j++){

            if(prob(gen) < 0.5){
                grid[i][j] = 1;
                totalSpin++;
            }else{
                grid[i][j] = -1;
                totalSpin--;
            }
        }
    }
    totalU = Utotal();
}

Ising2d::~Ising2d(){
    for (int i = 0; i < size; i++){
        delete[] grid[i];
    }
    delete[] grid;
}

float Ising2d::du(int i, int j){ 
    int top = (j == 0) ? grid[i][size-1] : grid[i][j-1];
    int bottom = (j == size-1) ? grid[i][0] : grid[i][j+1];
    int left = (i == 0) ? grid[size-1][j] : grid[i-1][j];
    int right = (i == size-1) ? grid[0][j] : grid[i+1][j];

    return 2*grid[i][j]*(top+bottom+left+right);
}

//  returns total energy of a given grid state
int Ising2d::Utotal(){
    int U = 0;
    // loop over all dipoles in grid
    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < size; j++){
            // calc energy for this dipole and add to sum
            int top = (j == 0) ? grid[i][size-1] : grid[i][j-1];
            int bottom = (j == size-1) ? grid[i][0] : grid[i][j+1];
            int left = (i == 0) ? grid[size-1][j] : grid[i-1][j];
            int right = (i == size-1) ? grid[0][j] : grid[i+1][j];
            U += grid[i][j]*(top+bottom+left+right);
        }
    }
    //return half U because we count each interaction twice
    return U/2;
}

int Ising2d::SpinTotal() {
    int spin = 0;
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            spin += grid[i][j];
        }
    }
    return spin;
}

void Ising2d::Run(){

    for(long long i = 0; i < iter;i++){
        int x = prob(gen)*size;
        int y = prob(gen)*size; 
        float E = du(x, y);

        if(E <= 0){
            grid[x][y] *= -1;
            // totalSpin += 2*grid[x][y];
        }else if(prob(gen) < exp(-E/temp)){
            grid[x][y] *= -1;
            // totalSpin += 2*grid[x][y];
        }

        // if (i%10'000==0) {
        //     DrawToImage(std::to_string(i).c_str());
        // }
    }

}

void Ising2d::timeSeries(std::ofstream& file) {

    for(long long i = 0; i < iter;i++){
        int x = prob(gen)*size;
        int y = prob(gen)*size; 
        float E = du(x, y);

        //flip if energy reduced
        if(E <= 0){
            grid[x][y] *= -1;
            totalSpin += 2*grid[x][y];
            totalU += E;
        }//flip with probability if energy increase
        else if(prob(gen) < exp(-E/temp)){
            grid[x][y] *= -1;
            totalSpin += 2*grid[x][y];
            totalU += E;
        }
        if (i % (iter/100) == 0) {
            file << temp << "," << int(i/(iter/100)) << "," << totalU << "," << totalSpin << "\n";
        }
    }
}

void Ising2d::tempRun(std::ofstream& file)  {

    for(long long i = 0; i < iter;i++){
        int x = prob(gen)*size;
        int y = prob(gen)*size; 
        float E = du(x, y);

        if(E <= 0){
            grid[x][y] *= -1;
            totalSpin += 2*grid[x][y];
        }else if(prob(gen) < exp(-E/temp)){
            grid[x][y] *= -1;
            totalSpin += 2*grid[x][y];
        }

        if(i%10000 == 0){
            file << i/10000 << "," << totalSpin <<"\n";
        }
    }
    file.close();
}

void Ising2d::smoothTemp(int runs, float maxTemp, std::ofstream& file) {

    grid = new int*[size];
    // initialize model with random spins
    for(int i = 0; i < size; i++) {

        grid[i] = new int[size];
        for(int j = 0; j < size; j++){
            grid[i][j] = 1;
            totalSpin++;
        }
    }
    for (int i = 0; i < runs; i++) {
        temp = (float(i)/runs) * maxTemp;
        Run();
        file << temp << "," << Utotal() << "," << SpinTotal() << "\n";
        std::cout << i << std::endl;
    }
}

void Ising2d::DrawToImage(const char* name){
    std::string n = std::string("../video/") + std::string(name) + std::string(".ppm");
    std::ofstream file(n);
    file << "P3\n" << size  << " " << size << "\n" << "255\n";

    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < size; j++){
            if(grid[i][j] == 1) file << "255 255 255 ";
            else file << "0 0 0 ";
        }
    }

    file.close();
    std::string conv = std::string("convert ../video/")
        + std::string(name)
        + std::string(".ppm ../video/")
        + std::string(name)
        + std::string(".png");

    std::string rm = std::string("rm ../video/")
        + std::string(name)
        + std::string(".ppm");

    int i = system(conv.c_str());
    i = system(rm.c_str());

}

int Ising2d::corrFunc(std::ofstream& file) {

    float corr[size/2];   // array to store summed correlations

    for (int x = 0; x < size/2; x++) {
        //init correlation to 0
        corr[x] = 0;
    }

    // loop over each dipole
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {

            // loop over lengths between 1 and size/2
            for (int r = 1; r <= size/2; r++) {

                // calc vetical correlation
                int top = (j - r >= 0) ? grid[i][(j - r)] : grid[i][(size + (j - r))];
                int bottom = (j + r <= size-1) ? grid[i][(j + r)] : grid[i][(j + r - size)];
                // calc horizontal correlation
                int left = (i - r >= 0) ? grid[(i - r)][j] : grid[(size + (i - r))][j];
                int right = (i + r <= size-1) ? grid[(i + r)][j] : grid[(i + r - size)][j];

                corr[r-1] += grid[i][j] * (top+bottom+left+right);
            }
        }
    }

    float max = 0.0f;
    float min = 1.0f;
    for (int i = 1; i <= size/2; i++) {

        std::cout << corr[i-1] << std::endl;
        // take average over all pairs 
        corr[i-1] = float(corr[i-1])/float(size*size*4);

        // subtract mean square magnetization
        // corr[i-1] -= float(Utotal())/float(size*size) * float(Utotal())/float(size*size);

        // track max and min
        if (corr[i-1] < min) min = corr[i-1];
        if (corr[i-1] > max) max = corr[i-1];

        file << i-1 << "," << corr[i-1] << "\n"; 
    }

    // find correlation length
    int corrlength = 0;
    float thresh = min + (max - min)/M_Ef;
    for (int i = 1; i < size/2; i++) {
        if (corr[i-1] <= thresh) {
            corrlength = i;
            break;
        }
    }
        
    file.close();
    return corrlength;
}

int Ising2d::corrLength() {

    float corr[size/2];   // array to store summed correlations

    for (int x = 0; x < size/2; x++) {
        //init correlation to 0
        corr[x] = 0;
    }

    // loop over each dipole
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {

            // loop over lengths between 1 and size/2
            for (int r = 1; r <= size/2; r++) {

                // calc vetical correlation
                int top = (j - r >= 0) ? grid[i][(j - r)] : grid[i][(size + (j - r))];
                int bottom = (j + r <= size-1) ? grid[i][(j + r)] : grid[i][(j + r - size)];
                // calc horizontal correlation
                int left = (i - r >= 0) ? grid[(i - r)][j] : grid[(size + (i - r))][j];
                int right = (i + r <= size-1) ? grid[(i + r)][j] : grid[(i + r - size)][j];

                corr[r-1] += grid[i][j] * (top+bottom+left+right);
            }
        }
    }

    float max = 0.0f;
    float min = 1.0f;
    for (int i = 1; i <= size/2; i++) {

        // std::cout << corr[i-1] << std::endl;
        // take average over all pairs 
        corr[i-1] = float(corr[i-1])/float(size*size*4);

        // subtract mean square magnetization
        // corr[i-1] -= float(Utotal())/float(size*size) * float(Utotal())/float(size*size);

        // track max and min
        if (corr[i-1] < min) min = corr[i-1];
        if (corr[i-1] > max) max = corr[i-1];

        // file << i-1 << "," << corr[i-1] << "\n"; 
    }

    // find correlation length
    int corrlength = 0;
    float thresh = min + (max - min)/M_Ef;
    for (int i = 1; i < size/2; i++) {
        if (corr[i-1] <= thresh) {
            corrlength = i;
            break;
        }
    }
        
    // file.close();
    return corrlength;
}

void Ising2d::smoothcorrLength(int runs, std::ofstream& file) {
    grid = new int*[size];
    // initialize model with random spins
    for(int i = 0; i < size; i++) {

        grid[i] = new int[size];
        for(int j = 0; j < size; j++){
            grid[i][j] = 1;
            totalSpin++;
        }
    }
    for (int i = 0; i < runs; i++) {
        temp = (float(i)/runs) * 5.0;
        Run();
        file << temp << "," << corrLength() << "\n";
        std::cout << i << std::endl;
    }
    file.close();
}

void Ising2d::blockSpin() {

    // declare block (uninitialized)
    int** block = new int*[size/3];

    for (int i = 0; i < size/3; i++) {
        block[i] = new int[size/3];

        // populate with average of grid spins
        for (int j = 0; j < size/3; j++) {
            //access 3*3 grid block centered at [i,j]
            // i,j      i+1,j       i+2,j
            // i,j+1    i+1,j+1     i+2,j+1
            // i,j+2    i+1,j+2     i+2,j+2
            int x = i*3;
            int y = j*3;
            int sum = 
                grid[x][y] +    grid[x+1][y] +  grid[x+2][y]
            +   grid[x][y+1]+   grid[x+1][y+1]+ grid[x+2][y+1]
            +   grid[x][y+2]+   grid[x+1][y+2]+ grid[x+2][y+2];
            block[i][j] = (sum < 0) ? -1 : 1;
            std::cout << block[i][j];
        }
        std::cout << std::endl;
    }
    // Draw block to image
    std::ofstream file("../video/block.ppm");

    file << "P3\n" << size/3  << " " << size/3 << "\n" << "255\n";

    for(size_t i = 0; i < size/3; i++){
        for(size_t j = 0; j < size/3; j++){
            if(block[i][j] == 1) file << "255 255 255 ";
            else file << "0 0 0 ";
        }
    }

    std::string name("block");
    file.close();
    std::string conv = std::string("convert ../video/")
        + std::string(name)
        + std::string(".ppm ../video/")
        + std::string(name)
        + std::string(".png");

    std::string rm = std::string("rm ../video/")
        + std::string(name)
        + std::string(".ppm");

    int i = system(conv.c_str());
    i = system(rm.c_str());


}
