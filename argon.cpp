// Translate this file with
//
// g++ -O3 spaceboddies.c -o spaceboddies
//
// Run it with
//
// ./spaceboddies
//
// Open Paraview (www.paraview.org) and do the following:
// - Select File/open and select all the results files. Press the Apply button.
// - Click into the left visualisation window (usually titled Layout #1).
// - Click the result-* item in the window Pipeline Browser. Ensure that your Layout #1 and the item result-* is marked.
// - Select Filters/Alphabetical/TableToPoints. Press Apply button.
// - Switch the representation (on top) from Surface into Points.
// - Press the play button and adopt colours and point sizes.
// - For some Paraview versions, you have to mark your TableToPoints item (usually called TableToPoints1) and explicitly select that X Column is x, Y Column is y, and Z Column is z.
// - What is pretty cool is the Filter TemporalParticlesToPathlines. If you set Mask Points to 1, you see a part of the trajactory.
//
// (C) 2015 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <chrono>
#include <thread>

const int N = 10;

const double a = 0.00001;
const double s = a;

double pos[3][4];
double v[3][3];

void setUp() {
  pos[0][0] = 0.4;
  pos[0][1] = 0.5;
  pos[0][2] = 0.5;

  pos[1][0] = 0.6;
  pos[1][1] = 0.5;
  pos[1][2] = 0.5;

  pos[2][0] = 0.7;
  pos[2][1] = 0.0;
  pos[2][2] = 0.0;

  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;

  v[2][0] = 0.0;
  v[2][1] = 0.0;
  v[2][2] = 0.0;

}

auto leonard_differential = [](double r) -> double {
  return 4*a * ((12 * (std::pow(s, 12) / std::pow(r, 13))) -
                (6 * (std::pow(s, 6) / std::pow(r, 7))));
};

void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".csv";
  std::ostream& out = std::cout;
  // std::ofstream out( filename.str().c_str() );

  out << "x, y, z, d" << std::endl;

  for (int i=0; i<3; i++) {
    out << pos[i][0]
        << ","
        << pos[i][1]
        << ","
        << pos[i][2]
        << ","
        << pos[i][3]
        << std::endl;
  }
}

double euler(std::function<double(double)> dydx, double dist, int iters) {

  double x = 1;
  double y = 0;
  double dx = (dist - x) / (iters - 1);

  for(int i = 0; i < iters; i++) {
    double dy = dx * dydx(x);
    x += dx;
    y += dy;
  }

  return y;
}


double potential(double r) {

  return 4*a * (std::pow((s/r), 12) - std::pow((s/r),6));

}


double force_potential(double r){

  return 4*a * (12 * (std::pow(s,12) / std::pow(r,13)) - 6 * (std::pow(s, 6)/std::pow(r, 7)));

}

void updateBody() {

  double force[3];

  force[0] = 0.0;
  force[1] = 0.0;
  force[2] = 0.0;


  const double timeStepSize = 0.1;

  for(int i = 0; i < N -1; i++){
    for (int j = i + 1; j < N; j++){

      double dx = pos[i][0] - pos[j][0];
      double dy = pos[i][1] - pos[j][1];
      double dz = pos[i][2] - pos[j][2];

      const double distance = sqrt(
          dx * dx +
          dy * dy +
          dz * dz
      );

      pos[i][3] = distance;

      force[0]  += dx * force_potential(distance);
      force[1]  += dy * force_potential(distance);
      force[2]  += dz * force_potential(distance);

    }

    pos[i][0] = timeStepSize * v[2][0];
    pos[i][1] = timeStepSize * v[2][1];
    pos[i][2] = timeStepSize * v[2][2];

    v[i][0] += timeStepSize * force[0];
    v[i][1] += timeStepSize * force[1];
    v[i][2] += timeStepSize * force[2];

  }


}


int main() {

  setUp();
  printCSVFile(0);

  const int timeSteps        = 200000;
  const int plotEveryKthStep = 10;
  for (int i=0; i<timeSteps; i++) {

    updateBody();
    if (i%plotEveryKthStep==0) {
      printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}