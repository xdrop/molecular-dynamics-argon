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

#define NO_OF_PARTICLES 2

#define LEOJ_CONST = 0.00001

double pos[NO_OF_PARTICLES][NO_OF_PARTICLES];
double v[NO_OF_PARTICLES][NO_OF_PARTICLES];


void setUp() {
  pos[0][0] = 0.4;
  pos[0][1] = 0.0;

  pos[1][0] = 2.0;
  pos[1][1] = 0.0;

  pos[2][0] = 3.0;
  pos[2][1] = 0.0;

  v[0][0] = 0.0;
  v[0][1] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;

  v[2][0] = 0.0;
  v[2][1] = 1.0;

}


void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter << ".csv";
  std::ofstream out(filename.str().c_str());

  out << "x, y, z" << std::endl;

  for (int i = 0; i < 3; i++) {
    out << pos[i][0]
        << ","
        << pos[i][1]
        << ","
        << pos[i][2]
        << std::endl;
  }
}


void updateBody() {
  double force[NO_OF_PARTICLES];
  force[0] = 0.0;
  force[1] = 0.0;

  auto leonard_differential = [](int r) -> double {
    return 4*LEOJ_CONST * ((12 * (std::pow(LEOJ_CONST, 12) / std::pow(r, 13))) -
        (6 * (std::pow(LEOJ_CONST, 6) / std::pow(r, 7))));
  };

  for (int i=0; i<2; i++) {
    const double distance = sqrt(
        (pos[2][0]-pos[i][0]) * (pos[2][0]-pos[i][0]) +
        (pos[2][1]-pos[i][1]) * (pos[2][1]-pos[i][1]) +
        (pos[2][2]-pos[i][2]) * (pos[2][2]-pos[i][2])
    );
    force[0] += (pos[i][0]-pos[2][0]) * mass[i]*mass[2] / distance / distance / distance ;
    force[1] += (pos[i][1]-pos[2][1]) * mass[i]*mass[2] / distance / distance / distance ;
    force[2] += (pos[i][2]-pos[2][2]) * mass[i]*mass[2] / distance / distance / distance ;
  }
  double dy = leonard_differential()

  const double timeStepSize = 0.0001;

  pos[2][0] = pos[2][0] + timeStepSize * v[2][0];
  pos[2][1] = pos[2][1] + timeStepSize * v[2][1];

  v[2][0] = v[2][0] + timeStepSize * force[0];
  v[2][1] = v[2][1] + timeStepSize * force[1];
}


int main() {

  setUp();
  printCSVFile(0);

  const int timeSteps = 2000000;
  const int plotEveryKthStep = 100;
  for (int i = 0; i < timeSteps; i++) {
    updateBody();
    if (i % plotEveryKthStep == 0) {
      printCSVFile(i / plotEveryKthStep + 1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}