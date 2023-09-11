#pragma once
#include <cstdint>
#include <random>
#include <fstream>

class Ising2d{
public:
    int size;
    int totalU = 0;
    int totalSpin = 0;
    int spinUp = 0;
    int spinDown = 0;
    float temp;
    int** grid;
    long iter;
    

    Ising2d(int size, float temp, long iter);
    ~Ising2d();

    int Utotal();
    int SpinTotal();

    void Run();
    void tempRun(std::ofstream& file);
    void smoothTemp(int runs, float maxTemp, std::ofstream& file);
    void timeSeries(std::ofstream& file);
    int corrFunc(std::ofstream& file);
    void smoothcorrLength(int runs, std::ofstream& file);
    int corrLength();
    void DrawToImage(const char* name);
    float du(int i, int j);

    // Delete the array created
};