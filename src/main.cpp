#include "Ising.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <cstdlib>
#include <cmath>

void timeSeries(int size, int iter, int runs) {
  std::ofstream f1("../data/" + std::to_string(size) + "TimeSeries01.txt");
  // run 50 for each size
  for (size_t i = 0; i < runs; i++) {
    // size 20
    Ising2d x1(size, 0.1, iter);
	x1.timeSeries(f1);
	Ising2d x2(size, 1.0, iter);
	x2.timeSeries(f1);
	Ising2d x3(size, 2.5, iter);
	x3.timeSeries(f1);
	Ising2d x4(size, 5.0, iter);
	x4.timeSeries(f1);

	std::cout << i << std::endl;

  }
  f1.close();
}

void meanField(int runs, int nn ,float thresh, float maxT) {

	// loop over temp [0,5]
	std::ofstream file("../data/meanfield" + std::to_string(nn) + ".txt");
	for (int i = 0; i < runs; i++) {
		double t = (double(i)/runs) * maxT;
		// loop over M values
		for (int j = 0; j < 10000; j++) {
			// test if func = 0
			double x = double(j)/10000.0f;
			if (std::abs(std::tanh((nn/t)*x)-x) <= thresh) {
				// add value to list
				file << t << "," << x << "\n";
			}
		}
		std::cout << t << std::endl;
	}
	file.close();
}

void temprun(int size, int iter, int runs) {
	std::ofstream file("../data/" + std::to_string(size) +"Temp.txt");
  // loop over temperature range
  for (int i = 0; i < runs; i++) {
	float temp = 5.0*(float(i)/runs);
	Ising2d x(size, temp, iter);
	x.Run();
	file << temp << "," << x.Utotal() << "," << x.SpinTotal() << "\n";
	printf("%i\n", i);
  }
}

void Correlation(int size, int iter) {
	// float temps[] = {0.1, 1.0, 1.5, 2.0, 2.27, 2.5, 3.0};
	float temps[] = {1.5};
	
	// loop over temps
	for (float t : temps) {
		Ising2d x(size, t, iter);
		x.Run();
		x.DrawToImage((std::to_string(size) + std::to_string(t)).c_str());
		std::ofstream file("../data/corr" + std::to_string(size) + std::to_string(t) + ".txt");
		x.corrFunc(file);
	}
}

void CorrLength(int size, int iter, int runs) {

	std::ofstream clfile("../data/correlations/corrlen" + std::to_string(size) + ".txt");
	for (int i = 0; i < runs; i++) {
		std::cout << i << std::endl;
		float temp = 5.0f * float(i)/float(runs);
		Ising2d x(size, temp, iter);
		x.Run();
		x.DrawToImage((std::to_string(size) + std::to_string(temp)).c_str());
		// std::ofstream file("../data/correlations/corr" + std::to_string(size) + std::to_string(t) + ".txt");
		int cl = x.corrLength();
		clfile << temp << "," << cl << "\n";
	}
	clfile.close();
}

int main() {
  auto begin = std::chrono::high_resolution_clock::now();
  int nruns = 500;

  // std::vector<std::thread> threads;
  // for(size_t i = 0; i < nruns; i++){
  // threads.push_back(std::thread(RunIsing,i, nruns)); } for(size_t i = 0; i <
  // nruns; i++){ threads[i].join(); }

	// Ising2d x(400, 0, 1'000'000);
	// std::ofstream file("../data/400smooth.txt");
	// x.smoothTemp(1000, 5, file);

	// meanField(1000, 12, float(1)/10000, 14.0);
	Ising2d x(100, 2.0, 1'000'000);
	x.Run();

	// Ising2d x(100, 1.5, 100'000'000);
	// x.Run();
	// std::ofstream file("../data/corr10001.txt");
	// x.corrFunc(file);
	// x.DrawToImage("corr");
	// CorrLength(40, 10'000'000, 50);
	// std::ofstream file("../data/correlations/corrlen20.txt");
	// Ising2d x(20, 0.0f, 1'000'000);
	// x.smoothcorrLength(100, file);


  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time Elapsed: " << (end - begin).count() / 1000000.0
            << "milliseconds" << std::endl;

  return 0;
}