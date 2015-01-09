#include"./head.h"
int main(int argc, char *argv[])
{
    int i, j, k, n;
    int grid_num, l_index, index_0;
    int index_1, index_2, index_3, index_4, index_5, index_6, index_7, index_8;
    int index_px, index_nx, index_py, index_ny, index_pz, index_nz;
    int index_upx, index_unx, index_upy, index_uny;
    int index_dpx, index_dnx, index_dpy, index_dny;
    int index_m1,  index_m2,  index_m3,  index_m4 ;
    int index_u1,  index_u2,  index_u3,  index_u4 ;
    int index_d1,  index_d2,  index_d3,  index_d4 ;

    float *del_dm = (float *)malloc(sizeof(float)*num);
    float *del_nu = (float *)malloc(sizeof(float)*num);

    float *vx_bg = (float *)malloc(sizeof(float)*num);
    float *vy_bg = (float *)malloc(sizeof(float)*num);
    float *vz_bg = (float *)malloc(sizeof(float)*num);
    float *vx_dm = (float *)malloc(sizeof(float)*num);
    float *vy_dm = (float *)malloc(sizeof(float)*num);
    float *vz_dm = (float *)malloc(sizeof(float)*num);

    double vec_px[3] = { 1., 0, 0};  //Unit vector, arranges as "x y z"
    double vec_nx[3] = {-1., 0, 0};  
    double vec_py[3] = {0,  1., 0};
    double vec_ny[3] = {0, -1., 0};
    double vec_pz[3] = {0, 0,  1.};
    double vec_nz[3] = {0, 0, -1.};

    double vec_upx[3] = { sqrt(2.)/2., 0, sqrt(2.)/2.};
    double vec_unx[3] = {-sqrt(2.)/2., 0, sqrt(2.)/2.};
    double vec_upy[3] = {0,  sqrt(2.)/2., sqrt(2.)/2.};
    double vec_uny[3] = {0, -sqrt(2.)/2., sqrt(2.)/2.}; 
    double vec_dpx[3] = { sqrt(2.)/2., 0, -sqrt(2.)/2.};
    double vec_dnx[3] = {-sqrt(2.)/2., 0, -sqrt(2.)/2.};
    double vec_dpy[3] = {0,  sqrt(2.)/2., -sqrt(2.)/2.};
    double vec_dny[3] = {0, -sqrt(2.)/2., -sqrt(2.)/2.};
    double vec_m1[3] = { sqrt(2.)/2., -sqrt(2.)/2., 0};
    double vec_m2[3] = { sqrt(2.)/2.,  sqrt(2.)/2., 0};
    double vec_m3[3] = {-sqrt(2.)/2.,  sqrt(2.)/2., 0};
    double vec_m4[3] = {-sqrt(2.)/2., -sqrt(2.)/2., 0};

    double vec_u1[3] = { sqrt(3.)/3., -sqrt(3.)/3., sqrt(3.)/3.};    
    double vec_u2[3] = { sqrt(3.)/3.,  sqrt(3.)/3., sqrt(3.)/3.};
    double vec_u3[3] = {-sqrt(3.)/3.,  sqrt(3.)/3., sqrt(3.)/3.};
    double vec_u4[3] = {-sqrt(3.)/3., -sqrt(3.)/3., sqrt(3.)/3.};
    double vec_d1[3] = { sqrt(3.)/3., -sqrt(3.)/3., -sqrt(3.)/3.};
    double vec_d2[3] = { sqrt(3.)/3.,  sqrt(3.)/3., -sqrt(3.)/3.};
    double vec_d3[3] = {-sqrt(3.)/3.,  sqrt(3.)/3., -sqrt(3.)/3.};
    double vec_d4[3] = {-sqrt(3.)/3., -sqrt(3.)/3., -sqrt(3.)/3.};

    double vx_med, vy_med, vz_med, v_mag;
    double v_dir[3] = {0.};
    double mu;

    double corr_px, corr_nx, corr_py, corr_ny, corr_pz, corr_nz;
    double corr_upx, corr_unx, corr_upy, corr_uny;
    double corr_dpx, corr_dnx, corr_dpy, corr_dny;
    double corr_m1, corr_m2, corr_m3, corr_m4;
    double corr_u1, corr_u2, corr_u3, corr_u4;
    double corr_d1, corr_d2, corr_d3, corr_d4;

    double corr_sum1, corr_sum2, corr_sum3;

    double     dipole_1[corr_num] = {0.};
    int    dipole_num_1[corr_num] = {0.};
    double     dipole_2[corr_num] = {0.};
    int    dipole_num_2[corr_num] = {0.};
    double     dipole_3[corr_num] = {0.};
    int    dipole_num_3[corr_num] = {0.};

    double L1 = box/nc;
    double L2 = L1*sqrt(2.);
    double L3 = L1*sqrt(3.);

    int id;
    int p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    char name_den_dm[80];
    char name_den_nu[80];
    char name_vlx_dm[80];
    char name_vly_dm[80];
    char name_vlz_dm[80];
    char name_vlx_nu[80];
    char name_vly_nu[80];
    char name_vlz_nu[80];

    sprintf(name_den_dm, "%s%d%s", den_, id, _dmdat);
    sprintf(name_den_nu, "%s%d%s", den_, id, _nudat);
    sprintf(name_vlx_nu, "%s%d%s", velx_, id, _nubin);
    sprintf(name_vly_nu, "%s%d%s", vely_, id, _nubin);
    sprintf(name_vlz_nu, "%s%d%s", velz_, id, _nubin);
    sprintf(name_vlx_dm, "%s%d%s", velx_, id, _dmbin);
    sprintf(name_vly_dm, "%s%d%s", vely_, id, _dmbin);
    sprintf(name_vlz_dm, "%s%d%s", velz_, id, _dmbin);

    if(id == 0) printf("%s\n", name_den_dm);
    if(id == 0) printf("%s\n", name_vlx_dm);

    FILE *den_dm = fopen(name_den_dm, "rb");
    FILE *den_nu = fopen(name_den_nu, "rb");

    FILE *velx_bg = fopen(name_vlx_nu, "rb");
    FILE *vely_bg = fopen(name_vly_nu, "rb");
    FILE *velz_bg = fopen(name_vlz_nu, "rb");

    FILE *velx_dm = fopen(name_vlx_dm, "rb");
    FILE *vely_dm = fopen(name_vly_dm, "rb");
    FILE *velz_dm = fopen(name_vlz_dm, "rb");

    if(id == 0) printf("Begin to readin data!\n");
    fread(del_dm, sizeof(float)*num, 1, den_dm);
    fread(del_nu, sizeof(float)*num, 1, den_nu);

    fread(vx_bg, sizeof(float)*num, 1, velx_bg);
    fread(vy_bg, sizeof(float)*num, 1, vely_bg);
    fread(vz_bg, sizeof(float)*num, 1, velz_bg);

    fread(vx_dm, sizeof(float)*num, 1, velx_dm);
    fread(vy_dm, sizeof(float)*num, 1, vely_dm);
    fread(vz_dm, sizeof(float)*num, 1, velz_dm);
    if(id == 0) printf("Finish reading data!\n");

    for(i = 0; i < num; i++)
    {
        del_dm[i] = del_dm[i] - 1.;
        del_nu[i] = del_nu[i] - 1.;

        vx_bg[i] = vx_bg[i] - vx_dm[i];
        vy_bg[i] = vy_bg[i] - vy_dm[i];
        vz_bg[i] = vz_bg[i] - vz_dm[i];

        if(((i+1)%(nc*nc) == 0) && id ==0) 
            printf("Processing initial data..%d\n", ((i+1)/nc/nc));
    }

    grid_num = nc;

    for (n = 0; n < corr_num; n++)
    {
        l_index = grid_num - 1;
     
        if(id == 0) printf("Processing...n=%d\n", n);

        for(k = 1; k < l_index; k++)  // Increasing z index
        {
            for(j = 1; j < l_index; j++)  // Increasing y index
            {
                for(i = 1; i < l_index; i++) // Increasing x index
                {
                    index_0  = index_f (k, j, i  ); //Origin point  
                    index_px = index_f (k, j, i+1); //Points at 1 away
                    index_nx = index_f (k, j, i-1);
                    index_py = index_f (k, j+1, i);
                    index_ny = index_f (k, j-1, i);
                    index_pz = index_f (k+1, j, i);
                    index_nz = index_f (k-1, j, i);

                    index_upx = index_f (k+1, j, i+1); //Points at sqrt2 away
                    index_unx = index_f (k+1, j, i-1);
                    index_upy = index_f (k+1, j+1, i);
                    index_uny = index_f (k+1, j-1, i);
                    index_dpx = index_f (k-1, j, i+1);
                    index_dnx = index_f (k-1, j, i-1);
                    index_dpy = index_f (k-1, j+1, i);
                    index_dny = index_f (k-1, j-1, i);
                    index_m1 = index_f (k, j-1, i+1);
                    index_m2 = index_f (k, j+1, i+1);
                    index_m3 = index_f (k, j+1, i-1);
                    index_m4 = index_f (k, j-1, i-1);

                    index_u1 = index_f (k+1, j-1, i+1); //Points at sqrt3 away
                    index_u2 = index_f (k+1, j+1, i+1);
                    index_u3 = index_f (k+1, j+1, i-1);
                    index_u4 = index_f (k+1, j-1, i-1);
                    index_d1 = index_f (k-1, j-1, i+1);
                    index_d2 = index_f (k-1, j+1, i+1);
                    index_d3 = index_f (k-1, j+1, i-1);
                    index_d4 = index_f (k-1, j-1, i-1);

                    //Processing points at 1 away
                    vx_med = (vx_bg[index_0] + vx_bg[index_px])/2.;//px 1
                    vy_med = (vy_bg[index_0] + vy_bg[index_px])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_px])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_px);
                    corr_px = del_dm[index_0]*del_nu[index_px]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_nx])/2.;//nx 2
                    vy_med = (vy_bg[index_0] + vy_bg[index_nx])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_nx])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_nx);
                    corr_nx = del_dm[index_0]*del_nu[index_nx]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_py])/2.;//py 3
                    vy_med = (vy_bg[index_0] + vy_bg[index_py])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_py])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_py);
                    corr_py = del_dm[index_0]*del_nu[index_py]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_ny])/2.;//ny 4
                    vy_med = (vy_bg[index_0] + vy_bg[index_ny])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_ny])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_ny);
                    corr_ny = del_dm[index_0]*del_nu[index_ny]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_pz])/2.;//pz 5
                    vy_med = (vy_bg[index_0] + vy_bg[index_pz])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_pz])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_pz);
                    corr_pz = del_dm[index_0]*del_nu[index_pz]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_nz])/2.;//nz 6
                    vy_med = (vy_bg[index_0] + vy_bg[index_nz])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_nz])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_nz);
                    corr_nz = del_dm[index_0]*del_nu[index_nz]*mu;

                    corr_sum1 = corr_px + corr_nx + corr_py + corr_ny +
                                corr_pz + corr_nz;
                    corr_sum1 = corr_sum1/6.;

                    dipole_1[n]     = dipole_1[n] + corr_sum1;
                    dipole_num_1[n] = dipole_num_1[n] + 1; 

                    //Processing points at sqrt2 away
                    vx_med = (vx_bg[index_0] + vx_bg[index_upx])/2.;//upx 1
                    vy_med = (vy_bg[index_0] + vy_bg[index_upx])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_upx])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_upx);
                    corr_upx = del_dm[index_0]*del_nu[index_upx]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_unx])/2.;//unx 2
                    vy_med = (vy_bg[index_0] + vy_bg[index_unx])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_unx])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_unx);
                    corr_unx = del_dm[index_0]*del_nu[index_unx]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_upy])/2.;//upy 3
                    vy_med = (vy_bg[index_0] + vy_bg[index_upy])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_upy])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_upy);
                    corr_upy = del_dm[index_0]*del_nu[index_upy]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_uny])/2.;//uny 4
                    vy_med = (vy_bg[index_0] + vy_bg[index_uny])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_uny])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_uny);
                    corr_uny = del_dm[index_0]*del_nu[index_uny]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_dpx])/2.;//dpx 5
                    vy_med = (vy_bg[index_0] + vy_bg[index_dpx])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_dpx])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_dpx);
                    corr_dpx = del_dm[index_0]*del_nu[index_dpx]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_dnx])/2.;//dnx 6
                    vy_med = (vy_bg[index_0] + vy_bg[index_dnx])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_dnx])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_dnx);
                    corr_dnx = del_dm[index_0]*del_nu[index_dnx]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_dpy])/2.;//dpy 7
                    vy_med = (vy_bg[index_0] + vy_bg[index_dpy])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_dpy])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_dpy);
                    corr_dpy = del_dm[index_0]*del_nu[index_dpy]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_dny])/2.;//dny 8
                    vy_med = (vy_bg[index_0] + vy_bg[index_dny])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_dny])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_dny);
                    corr_dny = del_dm[index_0]*del_nu[index_dny]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_m1])/2.;//m1  9
                    vy_med = (vy_bg[index_0] + vy_bg[index_m1])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_m1])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_m1);
                    corr_m1 = del_dm[index_0]*del_nu[index_m1]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_m2])/2.;//m2  10
                    vy_med = (vy_bg[index_0] + vy_bg[index_m2])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_m2])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_m2);
                    corr_m2 = del_dm[index_0]*del_nu[index_m2]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_m3])/2.;//m3  11
                    vy_med = (vy_bg[index_0] + vy_bg[index_m3])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_m3])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_m3);
                    corr_m3 = del_dm[index_0]*del_nu[index_m3]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_m4])/2.;//m4  12
                    vy_med = (vy_bg[index_0] + vy_bg[index_m4])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_m4])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_m4);
                    corr_m4 = del_dm[index_0]*del_nu[index_m4]*mu;

                    corr_sum2 = corr_upx + corr_unx + corr_upy + corr_uny +
                                corr_dpx + corr_dnx + corr_dpy + corr_dny +
                                corr_m1  + corr_m2  + corr_m3  + corr_m4;
                    corr_sum2 = corr_sum2/12.;

                    dipole_2[n]     = dipole_2[n] + corr_sum2;
                    dipole_num_2[n] = dipole_num_2[n] + 1;
   
                    //Processing points sqrt3 away 
                    vx_med = (vx_bg[index_0] + vx_bg[index_u1])/2.;//u1  1
                    vy_med = (vy_bg[index_0] + vy_bg[index_u1])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_u1])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_u1);
                    corr_u1 = del_dm[index_0]*del_nu[index_u1]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_u2])/2.;//u2  2
                    vy_med = (vy_bg[index_0] + vy_bg[index_u2])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_u2])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_u2);
                    corr_u2 = del_dm[index_0]*del_nu[index_u2]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_u3])/2.;//u3  3
                    vy_med = (vy_bg[index_0] + vy_bg[index_u3])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_u3])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_u3);
                    corr_u3 = del_dm[index_0]*del_nu[index_u3]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_u4])/2.;//u4  4
                    vy_med = (vy_bg[index_0] + vy_bg[index_u4])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_u4])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_u4);
                    corr_u4 = del_dm[index_0]*del_nu[index_u4]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_d1])/2.;//d1  5
                    vy_med = (vy_bg[index_0] + vy_bg[index_d1])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_d1])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_d1);
                    corr_d1 = del_dm[index_0]*del_nu[index_d1]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_d2])/2.;//d2  6
                    vy_med = (vy_bg[index_0] + vy_bg[index_d2])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_d2])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_d2);
                    corr_d2 = del_dm[index_0]*del_nu[index_d2]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_d3])/2.;//d3  7
                    vy_med = (vy_bg[index_0] + vy_bg[index_d3])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_d3])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_d3);
                    corr_d3 = del_dm[index_0]*del_nu[index_d3]*mu;

                    vx_med = (vx_bg[index_0] + vx_bg[index_d4])/2.;//d4  8
                    vy_med = (vy_bg[index_0] + vy_bg[index_d4])/2.;
                    vz_med = (vz_bg[index_0] + vz_bg[index_d4])/2.;
                    v_mag = sqrt(vx_med*vx_med + vy_med*vy_med + vz_med*vz_med);
                    v_dir[0] = vx_med/v_mag;
                    v_dir[1] = vy_med/v_mag;
                    v_dir[2] = vz_med/v_mag;

                    mu = costheta(v_dir, vec_d4);
                    corr_d4 = del_dm[index_0]*del_nu[index_d4]*mu;

                    corr_sum3 = corr_u1 + corr_u2 + corr_u3 + corr_u4 +
                                corr_d1 + corr_d2 + corr_d3 + corr_d4;
                    corr_sum3 = corr_sum3/8.;

                    dipole_3[n]     = dipole_3[n] + corr_sum3;
                    dipole_num_3[n] = dipole_num_3[n] + 1;
                }
            }
        }

        dipole_1[n] = dipole_1[n]/dipole_num_1[n];
        dipole_2[n] = dipole_2[n]/dipole_num_2[n];
        dipole_3[n] = dipole_3[n]/dipole_num_3[n];

        if(n == (corr_num-1)) continue;

        grid_num = nc/pow(2, (n+1));

        for (i = 0; i < grid_num; i++)
        {
            if(id == 0) printf("Combining...%d\n", i);
            for(j = 0; j < grid_num; j++)
            {
                for(k = 0; k < grid_num; k++)
                {
                    index_0 = index_f (i, j, k);
                    index_1 = index_f (2*i  , 2*j  , 2*k  );
                    index_2 = index_f (2*i+1, 2*j  , 2*k  );
                    index_3 = index_f (2*i  , 2*j+1, 2*k  );
                    index_4 = index_f (2*i  , 2*j  , 2*k+1);
                    index_5 = index_f (2*i+1, 2*j+1, 2*k  );
                    index_6 = index_f (2*i  , 2*j+1, 2*k+1);
                    index_7 = index_f (2*i+1, 2*j  , 2*k+1);
                    index_8 = index_f (2*i+1, 2*j+1, 2*k+1);

                    del_dm[index_0]
                  = del_dm[index_1] + del_dm[index_2] +
                    del_dm[index_3] + del_dm[index_4] +
                    del_dm[index_5] + del_dm[index_6] +
                    del_dm[index_7] + del_dm[index_8];

                    del_nu[index_0]
                  = del_nu[index_1] + del_nu[index_2] +
                    del_nu[index_3] + del_nu[index_4] +
                    del_nu[index_5] + del_nu[index_6] +
                    del_nu[index_7] + del_nu[index_8];

                    vx_bg[index_0]
                  = vx_bg[index_1] + vx_bg[index_2] + 
                    vx_bg[index_3] + vx_bg[index_4] +
                    vx_bg[index_5] + vx_bg[index_6] +
                    vx_bg[index_7] + vx_bg[index_8];

                    vy_bg[index_0]
                  = vy_bg[index_1] + vy_bg[index_2] + 
                    vy_bg[index_3] + vy_bg[index_4] +
                    vy_bg[index_5] + vy_bg[index_6] +
                    vy_bg[index_7] + vy_bg[index_8];

                    vz_bg[index_0]
                  = vz_bg[index_1] + vz_bg[index_2] + 
                    vz_bg[index_3] + vz_bg[index_4] +
                    vz_bg[index_5] + vz_bg[index_6] +
                    vz_bg[index_7] + vz_bg[index_8];

                    del_dm[index_0] = del_dm[index_0]/8.;
                    del_nu[index_0] = del_nu[index_0]/8.;
                    vx_bg[index_0] = vx_bg[index_0]/8.;
                    vy_bg[index_0] = vy_bg[index_0]/8.;
                    vz_bg[index_0] = vz_bg[index_0]/8.;
                }
            }
        }    
    }

    char name_out[80];
    sprintf(name_out, "%s%d%s", out_, id, _d);
    FILE *out = fopen(name_out, "w");

    for(n = 0; n < corr_num; n++)
    {
        fprintf(out,"%f\t%f\t%d\n",L1*pow(2, n), dipole_1[n], dipole_num_1[n]);
        fprintf(out,"%f\t%f\t%d\n",L2*pow(2, n), dipole_2[n], dipole_num_2[n]);
        fprintf(out,"%f\t%f\t%d\n",L3*pow(2, n), dipole_3[n], dipole_num_3[n]);
        
        if(id == 0) printf("Outputing...%d\n", n);
    }

/*-----------------*/
    MPI_Finalize();
    return 0;
}

int index_f (int i, int j, int k)
{
    int n = (i*nc + j)*nc + k;
    return n;
}

double costheta (double v_dir[3], double vec[3])
{
    double product;
    product = vec[0]*v_dir[0] + vec[1]*v_dir[1] + vec[2]*v_dir[2];
    return product;
}























