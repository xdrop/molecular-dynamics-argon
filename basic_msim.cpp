// The C++ standard used is C++11
// g++ -O4 -std=c++11 simple_msim.cpp

#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <thread>

const int N = 2;
double pos[N][7];

#define EMPTY -1

double v[N][3];

const double a = 0.00001;

const double s = a;

double timeStepSize;

inline double periodicBoundaryDelta(int i, int j, int dim) {

  double dt = pos[j][dim] - pos[i][dim];

  if (dt > 0.5) dt = dt - 1.0;
  if (dt <= -0.5) dt = dt + 1.0;

  return -dt;

}

inline const double calcDistance(double dx, double dy, double dz);

inline double force_potential(double r) {

  return 4 * a * ((12 * (std::pow(s, 12) / std::pow(r, 13))) - (6 * (std::pow(s, 6) / std::pow(r, 7))));

}

const double calcDistance(double dx, double dy, double dz) {
  const double distance = sqrt(
      dx * dx +
      dy * dy +
      dz * dz
  );
  return distance;
}

float rand_FloatRange(float a, float b) {
  return ((b - a) * ((float) rand() / RAND_MAX)) + a;
}

void setUp() {
  pos[0][0] = 0.4;
  pos[0][1] = 0.5;
  pos[0][2] = 0.5;

  pos[1][0] = 0.6;
  pos[1][1] = 0.5;
  pos[1][2] = 0.5;

//  pos[2][0] = 0.5;
//  pos[2][1] = 0.98;
//  pos[2][2] = 0.5;
//
//  pos[3][0] = 0.5;
//  pos[3][1] = 0.02;
//  pos[3][2] = 0.5;
//
//  pos[4][0] = 0.02;
//  pos[4][1] = 0.02;
//  pos[4][2] = 0.5;
//
//
//  pos[5][0] = 0.98;
//  pos[5][1] = 0.98;
//  pos[5][2] = 0.5;


  v[0][0] = 0.0000005;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = -0.0000005;
  v[1][1] = 0.0;
  v[1][2] = 0.0;


}

double randUnit() {
  return rand_FloatRange(0.0, 1.0);
}

void randSetup() {

  for (int i = 0; i < N; i++) {
    pos[i][0] = randUnit();
    pos[i][1] = randUnit();
    pos[i][2] = randUnit();
    v[i][0] = 0;
    v[i][1] = 0;
    v[i][2] = 0;
  }

}

void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter << ".csv";
  std::ofstream out(filename.str().c_str());

  out << "x, y, z, force, dist, wrap, vel" << std::endl;

  for (int i = 0; i < N; i++) {
    out << pos[i][0]
        << ","
        << pos[i][1]
        << ","
        << pos[i][2]
        << ","
        << pos[i][3]
        << ","
        << pos[i][4]
        << ","
        << pos[i][5]
        << ","
        << pos[i][6]
        << std::endl;
  }
}

double wrap(double pos) {
  double wrapped = fmod(pos, 1.0);
  double d = wrapped > 0.0 ? wrapped : wrapped + 1;
  if(!(d < 1.0))
    std::cout << "break";
  return d;
}

void updateBody(int begin, int end) {

  double minDistance = 100.0;

  for (int i = begin; i < end; i++) {
    double force[3];

    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    // for every neighbour j in the vertlet list
    for (int j = 0; j < N; j++) {

      if(i==j) continue;

      // calculate difference (with periodic boundary) in x,y,z
      double dx = periodicBoundaryDelta(i, j, 0);
      double dy = periodicBoundaryDelta(i, j, 1);
      double dz = periodicBoundaryDelta(i, j, 2);

      const double distance = sqrt(dx * dx + dy * dy + dz * dz);

      if (distance < minDistance)
        minDistance = distance;

      pos[i][4] = distance;

      double potential = force_potential(distance);

      force[0] += dx * potential;
      force[1] += dy * potential;
      force[2] += dz * potential;
    }

    if (minDistance < 0.0002) {
      timeStepSize = (1E6 * (minDistance/0.0002)) * 0.0000001;
    }
    else {
      timeStepSize = 1.0E2;
    }

    double displacement[3];
    displacement[0] = timeStepSize * v[i][0];
    displacement[1] = timeStepSize * v[i][1];
    displacement[2] = timeStepSize * v[i][2];

    assert(pos[i][0] > 0.0);
    assert(pos[i][0] < 1.0);

    // perform the actual movement of the particles while taking into account the
    // wrap around force
    pos[i][0] = wrap(pos[i][0] + displacement[0]);
    pos[i][1] = wrap(pos[i][1] + displacement[1]);
    pos[i][2] = wrap(pos[i][2] + displacement[2]);

    assert(pos[i][0] > 0.0);
    assert(pos[i][0] < 1.0);


    pos[i][3] = force[0];
    pos[i][6] = fabs(v[i][0]) + fabs(v[i][1]) + fabs(v[i][2]);


    v[i][0] += timeStepSize * force[0];
    v[i][1] += timeStepSize * force[1];
    v[i][2] += timeStepSize * force[2];


  }


}

int main() {

  srand(time(NULL));

//  randSetup();
  setUp();
  printCSVFile(0);

  const int timeSteps = 16000;
  const int plotEveryKthStep = 8;


  // create a few threads to parallelize the simulation
  int noOfThreads = 2;
  std::vector<std::thread> threads(noOfThreads);
  // the "piece" of work that each thread should do
  int grainSize = N / noOfThreads;


  for (int i = 0; i < timeSteps; i++) {
//    int iteration = 0;
//
//    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
//      *it = std::thread(updateBody, iteration, iteration + grainSize);
//      iteration += grainSize;
//    }
//    threads.back() = std::thread(updateBody, iteration, N);
//
//    for (auto &&i: threads) {
//      i.join();
//    }

    updateBody(0,N);

    if (i % plotEveryKthStep == 0) {
      int i1 = i / plotEveryKthStep + 1;
      std::cout << "Round " << i1 << std::endl;
      if (i1 % 100 == 0) {
        std::cout << "Round " << i1 << std::endl;
      }
      printCSVFile(i1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}


