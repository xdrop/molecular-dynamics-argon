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
#include <vector>

const int N = 100;
double pos[N][7];

#define MOD(a, b) ((((a)%(b))+(b))%(b))
#define EMPTY -1

double v[N][3];

const double a = 0.01;

const double s = a;
const double cutoff = 10 * s;
const double rSqrd = cutoff * cutoff;


const double skinDiameter = cutoff + 0.5 * cutoff;
const double cellSize = skinDiameter;
const int M = 7;

double maximumDisplacementPerParticle[N][3];
double maximumDistanceDisplaced[2];

std::vector<int> *vertletList;

int cells[M][M][M];
int neighbourCells[M][M][M][27][3];
int next[N];

const double timeStepSize = 0.0001;

inline double getDelta(int i, int j, int dim) {

  double dt = pos[i][dim] - pos[j][dim];
  double absolute = fabs(dt);

  if (1 - absolute < absolute) {
    if (dt < 0)
      dt = 1 - absolute;
    else
      dt = -(1 - absolute);
  }

  return dt;

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
  pos[0][0] = 0.44;
  pos[0][1] = 0.5;
  pos[0][2] = 0.5;

  pos[1][0] = 0.56;
  pos[1][1] = 0.5;
  pos[1][2] = 0.5;

//  pos[2][0] = 0.5;
//  pos[2][1] = 0.53;
//  pos[2][2] = 0.5;
//
//  pos[3][0] = 0.5;
//  pos[3][1] = 0.46;
//  pos[3][2] = 0.5;
//
//  pos[4][0] = 0.0;
//  pos[4][1] = 0.0;
//  pos[4][2] = 0.0;



  v[0][0] = 0.0;
  v[0][1] = 0.0;
  v[0][2] = 0.0;

  v[1][0] = 0.0;
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
  return wrapped > 0.0 ? wrapped : wrapped + 1;
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

    // for each direction, calculate the index of the cell which should receive it
    int x = (int) (pos[i][0] / cellSize);
    int y = (int) (pos[i][1] / cellSize);
    int z = (int) (pos[i][2] / cellSize);

    // the current cell now points to the (soon to be previous) head of the cell
    next[i] = cells[x][y][z];

    // the head of the cell is now this point
    cells[x][y][z] = i;

  }

}

void writeNeighbourCells() {
  for (int x = 0; x < M; x++) {
    for (int y = 0; y < M; y++) {
      for (int z = 0; z < M; z++) {

        // generate neighbours in 3 dimensions

        int i = 0;

        for (int xx = x - 1; xx <= x + 1; xx++) {
          for (int yy = y - 1; yy <= y + 1; yy++) {
            for (int zz = z - 1; zz <= z + 1; zz++) {


              neighbourCells[x][y][z][i][0] = MOD(xx, M);
              neighbourCells[x][y][z][i][1] = MOD(yy, M);
              neighbourCells[x][y][z][i][2] = MOD(zz, M);

              i++;

            }
          }
        }
      }
    }
  }
}

void updateBody() {

  double distA = 0;

  for (int i = 0; i < N; i++) {
    double force[3];

    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    // for every neighbour j in the vertlet list
    for (auto &j : vertletList[i]) {

      // calculate difference in x,y,z
      double dx = getDelta(i, j, 0);
      double dy = getDelta(i, j, 1);
      double dz = getDelta(i, j, 2);

      // calculate the distance squared
      double dSqrd = dx * dx + dy * dy + dz * dz;

      // compare the distance with the cutoff squared
      if (dSqrd <= rSqrd) {

        const double distance = sqrt(dSqrd);

        distA = distance;

        double potential = force_potential(distance);

        force[0] += dx * potential;
        force[1] += dy * potential;
        force[2] += dz * potential;
      }
    }

    // perform the actual movement of the particles while taking into
    // account the wrap around force
    double displacement[3];
    displacement[0] = timeStepSize * v[i][0];
    displacement[1] = timeStepSize * v[i][1];
    displacement[2] = timeStepSize * v[i][2];

    maximumDisplacementPerParticle[i][0] += displacement[0];
    maximumDisplacementPerParticle[i][1] += displacement[1];
    maximumDisplacementPerParticle[i][2] += displacement[2];

    // we keep track of the two maximum displacement's between
    // particles and update them accordingly
    double absDist = maximumDisplacementPerParticle[i][0] * maximumDisplacementPerParticle[i][0] +
                     maximumDisplacementPerParticle[i][1] * maximumDisplacementPerParticle[i][1] +
                     maximumDisplacementPerParticle[i][2] * maximumDisplacementPerParticle[i][2];

    if (absDist > maximumDistanceDisplaced[0]) {
      maximumDistanceDisplaced[1] = maximumDistanceDisplaced[0];
      maximumDistanceDisplaced[0] = absDist;
    } else {
      if (absDist > maximumDistanceDisplaced[1]) {
        maximumDistanceDisplaced[1] = absDist;
      }
    }

    pos[i][0] = wrap(pos[i][0] + displacement[0]);
    pos[i][1] = wrap(pos[i][1] + displacement[1]);
    pos[i][2] = wrap(pos[i][2] + displacement[2]);

    assert(pos[i][0] > 0.0);


    pos[i][3] = force[0];
    pos[i][4] = distA;
    pos[i][6] = v[i][0];


    v[i][0] += timeStepSize * force[0];
    v[i][1] += timeStepSize * force[1];
    v[i][2] += timeStepSize * force[2];


  }


}

void createVertletLists() {
  vertletList = new std::vector<int>[N];
  for (int i = 0; i < N; i++) {
    vertletList[i] = std::vector<int>(N / 8, EMPTY);
  }
}

void updateNeighbourList() {

  for (int i = 0; i < N; i++) {

    // we use a vector to store the neighbours of the current cell
    // we don't know how many there are in advance so vector is
    // needed for this.
    std::vector<int> &neighbours = vertletList[i];

    // clear the previous neighbours
    neighbours.clear();

  }

  // we iterate through all the cells in three dimensions
  for (int x = 0; x < M; x++) {
    for (int y = 0; y < M; y++) {
      for (int z = 0; z < M; z++) {

        // for each of the 27 neighbours of this cell
        for (int neighbourIndex = 0; neighbourIndex < 27; neighbourIndex++) {

          // fetch the neighbours x,y,z coordinates which
          // are precalculated
          int xx = neighbourCells[x][y][z][neighbourIndex][0];
          int yy = neighbourCells[x][y][z][neighbourIndex][1];
          int zz = neighbourCells[x][y][z][neighbourIndex][2];


          // fetch the "head" particle of the current cell
          int i = cells[x][y][z];


          // while the "head" continues to point at a particle
          while (i != EMPTY) {

            // get the neighbours of particle i
            std::vector<int> &neighbours = vertletList[i];

            // fetch the "head" of the neighbouring cell
            int j = cells[xx][yy][zz];

            // while the neighbour cell contains more particles
            while (j != EMPTY) {

              if (i != j) {
                // j is a neighbour of i
                neighbours.push_back(j);
              }

              j = next[j];

            }

            i = next[i];

          }

        }
      }
    }
  }

}

bool shouldUpdateNeighbourList() {
  return maximumDistanceDisplaced[0] + maximumDistanceDisplaced[1] > skinDiameter - cutoff;
}

int main() {

  srand(time(NULL));

//randSetup();
  setUp();
  printCSVFile(0);
  createVertletLists();
  writeNeighbourCells();
  cell_update();
  updateNeighbourList();

  for (int i = 0; i < N; i++) {
    maximumDisplacementPerParticle[i][0] = 0.0;
    maximumDisplacementPerParticle[i][1] = 0.0;
    maximumDisplacementPerParticle[i][2] = 0.0;
  }

  int vertletSaved = 0;

  const int timeSteps = 20000000;
  const int plotEveryKthStep = 10000;
  for (int i = 0; i < timeSteps; i++) {
    if (shouldUpdateNeighbourList()) {
      cell_update();
      updateNeighbourList();
      for (int i = 0; i < N; i++) {
        maximumDisplacementPerParticle[i][0] = 0.0;
        maximumDisplacementPerParticle[i][1] = 0.0;
        maximumDisplacementPerParticle[i][2] = 0.0;
      }
    } else {
      vertletSaved++;
    }

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

  std::cout << "Vertlet saved: " << vertletSaved << " timesteps" << std::endl;

  return 0;
}


