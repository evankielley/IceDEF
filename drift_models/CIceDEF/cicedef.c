/*Compile with command: gcc -g cicedef.c -lnetcdf -lm -o runcicedef */

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>

/* Define physical constants  */
#define PI 3.14159265358979323846

/* This is the name of the data file we will read. */
#define FILE_NAME "../Inputs/test_data.nc"

/* We are reading 3D data, a 340 x 360 x 122 grid. */
#define NX 340
#define NY 360
#define NZ 122

double UAF[NX][NY][NZ];
double VAF[NX][NY][NZ];
double UWF[NX][NY][NZ];
double VWF[NX][NY][NZ];
double SST[NX][NY][NZ];
int MASK[NX][NY];


/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main()
{

    /********************* GET INPUT DATA FROM NETCDF FILE ************************/


    /* This will be the netCDF ID for the file and data variable. */
    int ncid, varid;

    /* Loop indexes, and error handling. */
    int xind, yind, zind, retval;

    /* Open the file. NC_NOWRITE gives read-only access to the file.*/
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    printf("file open");


    /* uaF */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "uaF", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, varid, &UAF[0][0][0])))
        ERR(retval);

    printf("\nuaF data read\n");

    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            for (zind = 0; zind < 3; zind++)
                printf("%f ", UAF[xind][yind][zind]);


    /* vaF */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "vaF", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, varid, &VAF[0][0][0])))
        ERR(retval);

    printf("\nvaF data read\n");

    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            for (zind = 0; zind < 3; zind++)
                printf("%f ", VAF[xind][yind][zind]);


    /* uwF */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "uwF", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, varid, &UWF[0][0][0])))
        ERR(retval);

    printf("\nuwF data read\n");

    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            for (zind = 0; zind < 3; zind++)
                printf("%f ", UWF[xind][yind][zind]);


    /* vwF */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "vwF", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, varid, &VWF[0][0][0])))
        ERR(retval);

    printf("\nvwF data read\n");

    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            for (zind = 0; zind < 3; zind++)
                printf("%f ", VWF[xind][yind][zind]);


    /* sst */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "sst", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, varid, &SST[0][0][0])))
        ERR(retval);

    printf("\nsst data read\n");

    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            for (zind = 0; zind < 3; zind++)
                printf("%f ", SST[xind][yind][zind]);


    /* mask */

    /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "mask", &varid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = nc_get_var_int(ncid, varid, &MASK[0][0])))
        ERR(retval);

    printf("\nmask data read\n");

    /* Check the data. */
    for (xind = 0; xind < 3; xind++)
        for (yind = 0; yind < 3; yind++)
            printf("%d ", MASK[xind][yind]);


    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("\nfile closed\n");

    printf("\n*** SUCCESS reading example file %s!\n", FILE_NAME);

    
    /********************* MODEL PREAMBLE *************************/

    const double R = 6.378e6;
    const double om = 7.2921e-5;
    const double rhow = 1.027e3;
    const double rhoa = 1.2;
    const double rhoi = 0.850e3;
    const double drho = rhow - rhoi;
    const double Cw = 0.9;
    const double Ca = 1.3;
    const double gam = sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw));
    const double sst0 = -4;
    const double Cs1 = 1.5;
    const double Cs2 = 0.5;
    const double Cs3 = 0.1;
    const double CMv1 = 7.62e-3;
    const double CMv2 = 1.29e-3;
    const double CMe1 = 0.5;
    const double CMb1 = 0.58;
    const double CMb2 = 0.8;
    const double CMb3 = 0.2;

    const int t0 = 0;
    const int tn = 122;
    const int dt = 1*24*3600;

    int t_all[tn];
    for(int count = 0; count < tn; count++)
        t_all[count] = count;
    printf("last t_all: %d\n", t_all[tn-1]);

    double x_all[NX];
    for(int count1 = 0; count1 < NX; count1++)
        x_all[count1] = 275.125 + 0.250*count1; 
    printf("last x_all: %f\n", x_all[NX-1]);

    double y_all[NY];
    for(int count2 = 0; count2 < NY; count2++)
        y_all[count2] = 0.125 + 0.250*count2; 
    printf("last y_all: %f\n", y_all[NY-1]);

    double x[tn], y[tn], l[tn], w[tn], h[tn];
    x[0] = 310;
    y[0] = 50;
    l[0] = 600;
    w[0] = 500;
    h[0] = 400;

    int t = 0;

    while(t <= t_all[0])
    {
        printf("time: %d\n", t);

        double smallest_diffx = 0.250;
        double smallest_diffy = 0.250;
        double nearest_x, diffx, nearest_y, diffy;
        int nearest_x_ind, nearest_y_ind;

        /* Find nearest x */
        for(int jx = 0; jx < NX; jx++)
        {   
            //printf("x_all[jx]: %f\n", x_all[jx]);
            diffx = fabs(x[t] - x_all[jx]);
            //printf("diffx %f\n", diffx);
            if(diffx < smallest_diffx)
            {
                smallest_diffx = diffx;
                nearest_x = x_all[jx];
                nearest_x_ind = jx;
            }
        }
        printf("smallest diffx: %f\n", smallest_diffx);
        printf("nearest x: %f\n", nearest_x);
        printf("associated x index: %d\n", nearest_x_ind);

        /* Find nearest y */
        for(int jy = 0; jy < NY; jy++)
        {   
            //printf("y_all[jy]: %f\n", y_all[jy]);
            diffy = fabs(y[t] - y_all[jy]);
            //printf("diffy %f\n", diffy);
            if(diffy < smallest_diffy)
            {
                smallest_diffy = diffy;
                nearest_y = y_all[jy];
                nearest_y_ind = jy;
            }
        }
        printf("smallest diffy: %f\n", smallest_diffy);
        printf("nearest y: %f\n", nearest_y);
        printf("associated y index: %d\n", nearest_y_ind);

        double uaf, vaf, uwf, vwf, sst;
         
        uaf = UAF[nearest_x_ind][nearest_y_ind][t];
        vaf = VAF[nearest_x_ind][nearest_y_ind][t];
        uwf = UWF[nearest_x_ind][nearest_y_ind][t];
        vwf = VWF[nearest_x_ind][nearest_y_ind][t];
        sst = SST[nearest_x_ind][nearest_y_ind][t];

        
        /* DRIFT */

        double S, ff, lam;
        S = PI*((l[t]*w[t])/(l[t]+w[t]));
        ff = 2*om*sin((fabs(y[t])*PI)/180);
        lam = sqrt(2)*Cw*(gam*sqrt(pow(uaf, 2) + pow(vaf,2)))/(ff*S);

        double alpha, beta;

        if(lam < 0.1)
        {
            printf("Taylor approx used for alpha\n");
            alpha = lam*(pow(lam, 4)*(pow(lam, 4)*(pow(lam, 4)*(-0.0386699020961393*pow(lam, 4) + 
                    0.055242717280199) - 0.0883883476483184) + 
                    0.176776695296637) - 0.707106781186548);        
        }
        else
            alpha = (sqrt(2) / pow(lam, 3))*(1 - sqrt(1 + pow(lam, 4)));        

        if(lam < 0.6)
        {
            printf("Taylor approx used for beta\n");
            beta = pow(lam, 3) * (pow(lam, 4) * (pow(lam, 4) * 
                    (pow(lam, 4) * (pow(lam, 4) * (pow(lam, 4) * 
                    (pow(lam, 4) * (pow(lam, 4) * (pow(lam, 4) *
                    (0.0153268598203613 * pow(lam, 4) - 0.0151656272365985) + 
                    0.0180267866272764) + 0.0219176256311202) - 
                    0.0274446790511418) + 0.0357675015202851) - 
                    0.0493731785691779) + 0.0745776683282687) - 
                    0.132582521472478) + 0.353553390593274);
        }
        else
            beta = ((1./pow(lam, 3.)))*((sqrt((4. + pow(lam, 4.))) *
                    (sqrt(1. + pow(lam, 4.))) - 3.*pow(lam, 4.) - 4.));

        double viu, viv;
        viu = uwf + gam*(-alpha*vaf + beta*uaf);
        viv = vwf + gam*(alpha*uaf + beta*vaf);

        y[t+1] = y[t] + (viv*dt)*(180/(PI*R));
        x[t+1] = x[t] + (viu*dt)/(cos((((y[t] + y[t+1])/2)*PI)/180))*(180/(PI*R));

        printf("y[t+1]: %f\n", y[t+1]);
        printf("x[t+1]: %f\n", x[t+1]);

        t += 1;
    }

    return 0;
}
