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
	// x.Run();
	file << temp << "," << x.Utotal() << "," << x.SpinTotal() << "\n";
	printf("%i\n", i);
  }
}

void Correlation(int size, int iter, float temp, int runs) {

	float** corr = new float*[size/2];
	// initialize correlation
	for (size_t i = 0; i < size/2; i++) {
		// each corr length will have runs elements
		corr[i] = new float[runs];
	}

	// number of runs to average over
	for (size_t i = 0; i < runs; i++) {
		
		Ising2d x(size, temp, iter);
		// x.Run();
		float* current = x.corrFunc();
		for (size_t j = 0; j < size/2; j++) {
			corr[j][i] = current[j];
		}
		delete[] current;
		std::cout << "temp: " << temp << "\trun: " << i << "/" << runs << std::endl;
	}

	// find average
	std::ofstream file("../data/corrfunc/" 
		+ std::to_string(size) 
		+ "/S" + std::to_string(size) 
		+ "T" + std::to_string(temp) + ".txt");

	for (size_t i = 0; i < size/2; i++) {
	float sum = 0.0f;
		for (size_t j = 0; j < runs; j++) {
			sum += corr[i][j];
		}
		file << i << "," << sum/float(runs) << '\n';
	}

	file.close();

	for (int i = 0; i < size/2; i++) {
		delete[] corr[i];
	}
	delete[] corr;
}

void CorrLength(int size, int iter, int runs) {

	std::ofstream clfile("../data/correlations/corrlen" + std::to_string(size) + ".txt");
	for (int i = 0; i < runs; i++) {
		std::cout << i << std::endl;
		float temp = 5.0f * float(i)/float(runs);
		Ising2d x(size, temp, iter);
		// x.Run();
		x.DrawToImage((std::to_string(size) + std::to_string(temp)).c_str());
		// std::ofstream file("../data/correlations/corr" + std::to_string(size) + std::to_string(t) + ".txt");
		int cl = x.corrLength();
		clfile << temp << "," << cl << "\n";
	}
	clfile.close();
}

void BlockRun(int size , int iter, int tempruns, int avgruns) {
	
	std::ofstream file("../data/block/" + std::to_string(size) + ".txt");
	for (int i = 0; i < tempruns; i++) {
		float temp = 10.0f*float(i)/float(tempruns);
		Ising2d x(size, temp, iter);
		float bu = 0.0f;
		float bu2 = 0.0f;
		float u = 0.0f;
		float bs = 0.0f;
		float bs2 = 0.0f;
		float s = 0.0f;

		float buv[avgruns];
		float bu2v[avgruns];
		float uv[avgruns];
		float bsv[avgruns];
		float bs2v[avgruns];
		float sv[avgruns];
		// run n times and take average
		for (int j = 0; j < avgruns; j++) {
			float* vals = x.Runblock();
			bu += vals[0]/avgruns;
			buv[j] = vals[0];
			bu2 += vals[1]/avgruns;
			bu2v[j] = vals[1];
			u += vals[2]/avgruns;
			uv[j] = vals[2];
			bs += vals[3]/avgruns;
			bsv[j] = vals[3];
			bs2 += vals[4]/avgruns;
			bs2v[j] = vals[4];
			s += vals[5]/avgruns;
			sv[j] = vals[5];
		}
		// calculate variance
		float vbu = 0.0f;
		float vbu2 = 0.0f;
		float vu = 0.0f;
		float vbs = 0.0f;
		float vbs2 = 0.0f;
		float vs = 0.0f;
		for (int j = 0; j < avgruns; j++) {
			vbu += (bu-buv[j])*(bu-buv[j]);
			vbu2 +=(bu2-bu2v[j])*(bu2-bu2v[j]);
			vu +=(u-uv[j])*(u-uv[j]);
			vbs +=(bs-bsv[j])*(bs-bsv[j]);
			vbs2 +=(bs2-bs2v[j])*(bs2-bs2v[j]);
			vs +=(s-sv[j])*(s-sv[j]);
		}
		vbu =vbu/(avgruns-1);
		vbu2 = vbu2/(avgruns-1);
		vu = vu/(avgruns-1);
		vbs = vbs/(avgruns-1);
		vbs2 = vbs2/(avgruns-1);
		vs = vs/(avgruns-1);
		
		// print to file
		file << temp << "," << bu << "," << vbu 
			<< "," << bu2 << "," << vbu2 
			<< "," << u << "," << vu 
			<< "," << bs << "," << vbs
			<< "," << bs2 << "," << vbs2
			<< "," << s << "," << vs
			<< "\n";
		std::cout << "run: " << i << std::endl;
	}
	file.close();
}

int main() {
  auto begin = std::chrono::high_resolution_clock::now();
  int nruns = 500;

  // std::vector<std::thread> threads;
  // for(size_t i = 0; i < nruns; i++){
  // threads.push_back(std::thread(RunIsing,i, nruns)); } for(size_t i = 0; i <
  // nruns; i++){ threads[i].join(); }

	// float temps[] = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.27, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};

	// for (auto t : temps) {
	// 	Correlation(100, 1'000'000, t, 20);
	// }
	BlockRun(360, 10'000'000,50, 10);




  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time Elapsed: " << (end - begin).count() / 1000000.0
            << "milliseconds" << std::endl;

  return 0;
}