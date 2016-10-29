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


const int N = 20;
double pos[N][3];


double v[N][3];

const double a = 0.00001;
const double s = a;

double boundedIncrease(double a, double b);

void setUp() {
  pos[0][0] = 0.4;
  pos[0][1] = 0.0;
  pos[0][2] = 0.0;

  pos[1][0] = 2.0;
  pos[1][1] = 0.0;
  pos[1][2] = 0.0;

  pos[2][0] = 3.0;
  pos[2][1] = 0.0;
  pos[2][2] = 0.0;

  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;

  v[2][0] = 0.0;
  v[2][1] = 1.0;
  v[2][2] = 0.0;

}

float rand_FloatRange(float a, float b)
{
  return ((b-a)*((float)rand()/RAND_MAX))+a;
}


double randUnit(){
  return rand_FloatRange(0.0, 1.0);
}


void randSetup() {

  for(int i = 0; i < N; i++) {
     pos[i][0] = randUnit();
     pos[i][1] = randUnit();
     pos[i][2] = randUnit();
     v[i][0] = randUnit();
     v[i][1] = randUnit();
     v[i][2] = randUnit();
  }

}


void printCSVFile(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;

  for (int i=0; i<N; i++) {
    out << pos[i][0]
        << ","
        << pos[i][1]
        << ","
        << pos[i][2]
        << std::endl;
  }
}

double force_potential(double r){

  return 4*a * (12 * (std::pow(s,12) / std::pow(r,13)) - 6 * (std::pow(s, 6)/std::pow(r, 7)));

}



void updateBody() {
  double force[3];
  force[0] = 0.0;
  force[1] = 0.0;
  force[2] = 0.0;

  for(int i = 0; i < N -1; i++) {
    for (int j = i + 1; j < N; j++) {

      double dx = pos[i][0] - pos[j][0];
      double dy = pos[i][1] - pos[j][1];
      double dz = pos[i][2] - pos[j][2];

      const double distance = sqrt(
          dx * dx +
          dy * dy +
          dz * dz
      );

      force[0] += dx * force_potential(distance);
      force[1] += dy * force_potential(distance);
      force[2] += dz * force_potential(distance);

    }
    const double timeStepSize = 0.0001;

    pos[i][0] = boundedIncrease(pos[i][0], timeStepSize * v[2][0]);
    pos[i][1] = boundedIncrease(pos[i][1], timeStepSize * v[2][1]);
    pos[i][2] = boundedIncrease(pos[i][2], timeStepSize * v[2][2]);

    v[i][0] = v[i][0] + timeStepSize * force[0];
    v[i][1] = v[i][1] + timeStepSize * force[1];
    v[i][2] = v[i][2] + timeStepSize * force[2];
  }

}

double boundedIncrease(double cur, double change){
  double d = cur + change;

  if (d > 1.0){
    double offset = d - ((long)d);
    double distToAxis = 1.0 - cur;

    return offset - distToAxis;
  } else if (d < 0) {
    d = std::abs(d);
    double offset = d - ((long)d);
    double interim = offset - cur;

    return 1.0 - interim;
  }

  return d;
}


int main() {

  randSetup();
  printCSVFile(0);

  const int timeSteps        = 20000;
  const int plotEveryKthStep = 100;
  for (int i=0; i<timeSteps; i++) {
    updateBody();
    if (i%plotEveryKthStep==0) {
      printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
    }
  }

  return 0;
}