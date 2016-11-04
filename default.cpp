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
#include <cassert>


const int N = 2;
double pos[N][7];

#define EMPTY -1

double v[N][3];

const double a = 0.01;

const double s = a;
const double cutoff = 10 * s;
const double rSqrd = cutoff * cutoff;
const double cellSize = 0.2;
const int M = 25;

int cells[M][M][M];
int *neighbourCells[M][M][M];
int next[N];

const double timeStepSize = 0.0003;

inline double getDelta(int i, int j, int dim) {

  double dt = pos[i][dim] - pos[j][dim];

  if (1 - fabs(dt) < fabs(dt)) {
    if (dt > 0)
      dt = 1 - fabs(dt);
    else
      dt = -(1 - fabs(dt));
  }

  return dt;

}

double boundedIncrease(double a, double b);

inline const double calcDistance(double dx, double dy, double dz);

inline void calculateForces(int i, int j, double *in);

void setUp() {
  pos[0][0] = 0.46;
  pos[0][1] = 0.0;
  pos[0][2] = 0.0;

  pos[1][0] = 0.53;
  pos[1][1] = 0.0;
  pos[1][2] = 0.0;


  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
  v[1][1] = 0.0;
  v[1][2] = 0.0;


}

float rand_FloatRange(float a, float b) {
  return ((b - a) * ((float) rand() / RAND_MAX)) + a;
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

double boundedIncrease(double cur, double change) {

  double d = cur + change;

  if (d > 1.0) {
    double offset = d - ((long) d);
    double distToAxis = 1.0 - cur;

    return (offset - distToAxis);
  } else if (d < 0.0) {
    d = std::abs(d);
    double offset = d - ((long) d);
    double interim = offset - cur;

    return 1.0 - interim;
  }

  return d;
}

double force_potential(double r) {

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

void cell_update() {

  // all heads of all cells point to empty
  for (int x = 0; x < M; x++) {
    for (int y = 0; y < M; y++) {
      for (int z = 0; z < M; z++) {
        cells[x][y][z] = EMPTY;
      }
    }
  }


  for (int i = 0; i < N; i++) {

    int indices[3];

    // for each direction, calculate the index of the cell which should receive it
    for (int dir = 0; dir < 3; dir++) indices[dir] = (int) (pos[i][dir] / cellSize);

    // the current cell now points to the (soon to be previous) head of the cell
    next[i] = cells[indices[0]][indices[1]][indices[2]];

    // the head of the cell is now this point
    cells[indices[0]][indices[1]][indices[2]] = i;

    //std::cout << "Particle" << i << " is in cell " << indices[0] << std::endl;
  }

}

void writeNeighbourCells() {


  for (int x = 0; x < M; x++) {
    for (int y = 0; y < M; y++) {
      for (int z = 0; z < M; z++) {

        int *neighbours = new int[27 * 3];
        // generate neighbours in 3 dimensions

        int i = 0;

        for (int xx = x - 1; xx <= x + 1; xx++) {
          for (int yy = y - 1; yy <= y + 1; yy++) {
            for (int zz = z - 1; zz <= z + 1; zz++) {


              if (xx > 0 && xx < M && yy > 0 && yy < M && zz > 0 && zz < M) {
                // TODO: Memory leak
                ((int (*)[3]) neighbours)[i][0] = xx;
                ((int (*)[3]) neighbours)[i][1] = yy;
                ((int (*)[3]) neighbours)[i][2] = zz;
              } else{
                ((int (*)[3]) neighbours)[i][0] = EMPTY;
              }

              i++;

            }
          }
        }

        neighbourCells[x][y][z] = neighbours;
      }
    }
  }

  printf("%d", ((int (*)[3]) neighbourCells[5][5][5])[0][0]);
}

void updateBody() {


  double distA = 0;


  for (int x = 0; x < M; x++) {
    for (int y = 0; y < M; y++) {
      for (int z = 0; z < M; z++) {


        for (int di = 0; di < 27; di++) {
          int (*neighbours)[3] = (int (*)[3]) neighbourCells[x][y][z];

          int xx = neighbours[di][0];
          if(xx == EMPTY) continue;

          int yy = neighbours[di][1];
          int zz = neighbours[di][2];


          int i = cells[x][y][z];
          // current cell
          while (i != EMPTY) {

            int j = cells[xx][yy][zz];
            double force[3];
            bool dirty = false;

            force[0] = 0.0;
            force[1] = 0.0;
            force[2] = 0.0;


            while (j != EMPTY) {
              if (i != j) {

                double dx = getDelta(i, j, 0);
                double dy = getDelta(i, j, 1);
                double dz = getDelta(i, j, 2);

                double dSqrd = dx * dx + dy * dy + dz * dz;

                if (dSqrd <= rSqrd) {

                  dirty = true;

                  const double distance = sqrt(dSqrd);
                  distA = distance;

                  force[0] += dx * force_potential(distance);
                  force[1] += dy * force_potential(distance);
                  force[2] += dz * force_potential(distance);
                }
              }

              j = next[j];

            }

            if (dirty) {
              pos[i][0] = boundedIncrease(pos[i][0], timeStepSize * v[i][0]);
              pos[i][1] = boundedIncrease(pos[i][1], timeStepSize * v[i][1]);
              pos[i][2] = boundedIncrease(pos[i][2], timeStepSize * v[i][2]);


              pos[i][3] = force[0];
              pos[i][4] = distA;
              pos[i][5] = xx;
              pos[i][6] = v[i][0];

              v[i][0] += timeStepSize * force[0];
              v[i][1] += timeStepSize * force[1];
              v[i][2] += timeStepSize * force[2];
            }

            i = next[i];

          }

        }
      }
    }

  }

}


int main() {

//  randSetup();
  setUp();
  printCSVFile(0);
  writeNeighbourCells();

  const int timeSteps = 20000000;
  const int plotEveryKthStep = 10000;
  for (int i = 0; i < timeSteps; i++) {
    cell_update();
    updateBody();
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


