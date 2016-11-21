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

#define VARIABLE_TIME_STEP false

double v[N][3];

#if VARIABLE_TIME_STEP
const double a = 0.00001;
#else
const double a = 0.1;
#endif

const double s = a;

double timeStepSize = 0.0001;

bool lock = false;

inline double delta(double x1, double x2) {
  return x1 - x2;
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
  pos[0][0] = 0.25;
  pos[0][1] = 0.5;
  pos[0][2] = 0.5;

  pos[1][0] = 0.75;
  pos[1][1] = 0.5;
  pos[1][2] = 0.5;

  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;


}


void randSetup() {

  for (int i = 0; i < N; i++) {
    pos[i][0] = rand_FloatRange(0.0, 1.0);
    pos[i][1] = rand_FloatRange(0.0, 1.0);
    pos[i][2] = rand_FloatRange(0.0, 1.0);
    v[i][0] = rand_FloatRange(0.0, 0.000001);
    v[i][1] = rand_FloatRange(0.0, 0.000001);
    v[i][2] = rand_FloatRange(0.0, 0.000001);
  }

}

void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter << ".csv";
  std::ofstream out(filename.str().c_str());

  out << "x, y, z, force, dist, timestep, vel" << std::endl;

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
  double d = wrapped > 0.0 || wrapped == 0 ? wrapped : wrapped + 1;
  if (!(d < 1.0) || !(d > 0.0))
    std::cout << "break";
  return d;
}

void updateBody(int begin, int end) {

  double forces[N][3];

  for (int i = 0; i < N; i++) {
    forces[i][0] = 0.0;
    forces[i][1] = 0.0;
    forces[i][2] = 0.0;

  }

  double minDistance = 10000;
#if VARIABLE_TIME_STEP
  double maxForce = 0;
#endif

  for (int i = begin; i < end; i++) {

    for (int j = i + 1; j < N; j++) {


      for (int cX = -1; cX < 2; cX++) {
        for (int cY = -1; cY < 2; cY++) {
          for (int cZ = -1; cZ < 2; cZ++) {


            double nX = pos[j][0] + cX;
            double nY = pos[j][1] + cY;
            double nZ = pos[j][2] + cZ;


            // calculate difference (with periodic boundary) in x,y,z
            double dx = delta(pos[i][0], nX);
            double dy = delta(pos[i][1], nY);
            double dz = delta(pos[i][2], nZ);

            const double distance = sqrt(dx * dx + dy * dy + dz * dz);

            pos[i][4] = minDistance;

            double potential = force_potential(distance);

            double fDx = dx * potential;
            double fDy = dy * potential;
            double fDz = dz * potential;

            #if VARIABLE_TIME_STEP
            if (distance < minDistance)
              minDistance = distance;

            if(maxForce < fabs(fDx))
              maxForce = fabs(fDx);
            if(maxForce < fabs(fDy))
              maxForce = fabs(fDy);
            if(maxForce < fabs(fDz))
              maxForce = fabs(fDz);

            #endif

            // compute the forces on i
            forces[i][0] += fDx;
            forces[i][1] += fDy;
            forces[i][2] += fDz;

            // compute the forces on j
            forces[j][0] -= fDx;
            forces[j][1] -= fDy;
            forces[j][2] -= fDz;


          }
        }
      }


    }





    assert(pos[i][0] >= 0.0);
    assert(pos[i][0] < 1.0);

    // perform the actual movement of the particles while taking into account the
    // wrap around force
    double displacement[3];
    displacement[0] = timeStepSize * v[i][0];
    displacement[1] = timeStepSize * v[i][1];
    displacement[2] = timeStepSize * v[i][2];

#if VARIABLE_TIME_STEP

    if(minDistance < 0.00003){
      timeStepSize = 0.0001;
    } else if(minDistance + displacement[0] ){

    }
#endif

    pos[i][0] = wrap(pos[i][0] + displacement[0]);
    pos[i][1] = wrap(pos[i][1] + displacement[1]);
    pos[i][2] = wrap(pos[i][2] + displacement[2]);

    assert(pos[i][0] >= 0.0);
    assert(pos[i][0] < 1.0);


    pos[i][3] = forces[i][0];
    pos[i][5] = timeStepSize;
    pos[i][6] = fabs(v[i][0]) + fabs(v[i][1]) + fabs(v[i][2]);


    v[i][0] += timeStepSize * forces[i][0];
    v[i][1] += timeStepSize * forces[i][1];
    v[i][2] += timeStepSize * forces[i][2];


  }


}

int main() {

  srand(time(NULL));

//  randSetup();
  setUp();
  printCSVFile(0);

  const int timeSteps = 2000000;
  const int plotEveryKthStep = 1000;


  // create a few threads to parallelize the simulation
  int noOfThreads = 1;
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


    if(i == 4097){
      std::cout << "";
    }
    updateBody(0, N);

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


