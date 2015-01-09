#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>

#define nc       768        // number of cells per dimension in one node
#define num      (nc*nc*nc) // number of cells in one node
#define box      250.       // box size of the cube from one node
#define corr_num 9          // bins of correlation

char   den_[] = {"/home/zhm/test/den0"};
char _dmdat[] = {"_dm.dat"};
char _nudat[] = {"_nu.dat"};
char  velx_[] = {"/home/zhm/test/vx0"};
char  vely_[] = {"/home/zhm/test/vy0"};
char  velz_[] = {"/home/zhm/test/vz0"};
char _dmbin[] = {"_dm.bin"};
char _nubin[] = {"_nu.bin"};

char out_[] = {"/home/zhm/test/020_3fast/output_"};
char _d[]   = {".dat"};

int index_f (int i, int j, int k);
double costheta (double v_dir[3], double vec[3]);
