/* INCLUDE ESSENTIAL LIBRARIES */

#include <stdio.h>
#include <math.h>                                                     
#include <time.h>
#include <stdlib.h>
#include "stringlib.h"  

///////////////////////////////////////////////////////////////////////////////////////////////////////

/* DECLARATION OF VARIABLES */

void readConfigFile(char *, char *, char *, char *, double *, int *, int *, int *, double*, double *);

void initPositions(double **, int, double *);
void write_xyz(double **, int, int, double, FILE *);
void PrevCoord(double, double, double **, double **, double **);
void initVelocities(int, int, double, double, double **);
void write_Vel(int, int, double **, FILE *);
void computeForce_Energy(int, int, int, double, double, double, double, double, double **, double **); 
void write_Force(int, int, double **, FILE *);
void CoordUpdate(int, double, double, double, double, double **, double **, double **, double **);

///////////////////////////////////////////////////////////////////////////////////////////////////////

/*STATIC VARIABLES */

const double kB = 8.314;		 // Boltzman Constat times Avogadro # (Gass constant) (changed after talkiing with Russell)
const double mass = 0.0399;		 // Mass of Ar or any Lennard-Jones Liquid,Units in kg per mole (changed after talkiing with Russell)
const double eps = 0.210849;		 // Units of kcal/mol, it is a meassure of strength
const double sigma = 3.345;		 // Sigma value for Ar or any Lennard-Jones Liquid, it is a measure of range of potential, Units of Angstrums

///////////////////////////////////////////////////////////////////////////////////////////////////////

/*MAIN CODE */

int main(){
	int n_atoms;				// Number of atoms
	double box;				// Size of the cubic box
	double length;				// Length of cubic box
	int n_iter;				// Number of MD iterations
	double temp;				// Temperature
	int delta_write;			// How often to write coordinates and log file in MD iterations
	double cutoff;				// cutoff distance (Angstrom)
	double cutoff_squ;			// nonbonding interaction cutoff distance squared

	double kBT;				// Boltzman constant*Temperature Units of energy (J)
	//double mass;				// Mass of a Lennard-Jones Liquid 
	// double eps;				// Epsilon value in Lennard-Jones Potential expression
	double **atom_vel;			// Velocity array	
	
	double **coord;				// Particle coordinates array
	double **O_coord;	
	int i,j,k;				// Generic indeces
	int seed=1;				// Random seed for velocity initialization
	int iter;

	double sigma;				// sigma value
	double sigma6;				// LJL sigma^6 value

	
	double **forces_on_atom;		// Force acting on atom array fx, fy, fz
	
	double Tff;				// Total LJ potential energy 
	
	double Kt_en;				// Total Kinetic Energy
	double Tot_en;                  	// Total Energy
	double dt;				// delta t value	
	
	double dt2;				// delta t value squared
	double im;				// Argon's inverse mass


	char log_FileName[1024];         	// Log file name
	char vel_FileName[1024];         	// Velocity file name
	char traj_FileName[1024];         	// Trajectoruy file name
	char force_FileName[1024];         	// Force file name



	FILE *logOut;
	FILE *xyzOut;
	FILE *velOut;
	FILE *forceOut;


/* READ CONFIG FILE FROM STANDARD IN */

        readConfigFile(log_FileName, vel_FileName, traj_FileName, force_FileName, &temp, &n_atoms, &n_iter, &delta_write, &cutoff, &dt);
	
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
	dt2=dt*dt;
	cutoff_squ = cutoff*cutoff;
	kBT = kB*temp;
	im = 1/mass;
	printf("kB*T=%3.6E\n",kBT);

	// allocate coordinate  and velocity arrays
	coord = (double**) malloc(n_atoms*sizeof(double*));
	O_coord = (double**) malloc(n_atoms*sizeof(double*));
	// allocate velocity array memory
	atom_vel = (double**) calloc(n_atoms,sizeof(double*));   
	
	// allocate force array memory
	forces_on_atom = (double**) calloc(n_atoms,sizeof(double*));        
	  
	for (i=0;i<n_atoms;i++) {
		coord[i] = (double*) malloc(3*sizeof(double));
		O_coord[i] = (double*) malloc(3*sizeof(double));	
		atom_vel[i] = (double*) malloc(3*sizeof(double));	// Does not complie wiht calloc???
		forces_on_atom[i] = (double*) malloc(3*sizeof(double));
	}
	// open log files
	logOut = fopen(log_FileName,"w");
	xyzOut = fopen(traj_FileName,"w");
	velOut = fopen(vel_FileName,"w");
	forceOut = fopen(force_FileName, "w");

	// inintialize positions and velocities+ compute forces//

	n_iter=0;		                 			// Number of iteration before MD loop
	initPositions(coord, n_atoms, &box);
	write_xyz(coord, n_atoms, n_iter, box, xyzOut);
	

	initVelocities(n_atoms, seed, mass, kBT, atom_vel);
	write_Vel(n_atoms, n_iter, atom_vel, velOut);
	
	computeForce_Energy(n_atoms, n_iter, delta_write, box, cutoff_squ, eps, sigma6, Tff, coord, forces_on_atom);	
	write_Force(n_atoms, n_iter, forces_on_atom, forceOut);

	fflush(xyzOut);
	fflush(velOut);
	fflush(forceOut);
	

////////////////////////////////////////////////////

/* Run MD itterations  */
// Use Verlet integration
	
	
	for(iter=1;iter>n_iter;iter++) {
		CoordUpdate(n_atoms, im, dt, dt2, box, coord, atom_vel, forces_on_atom, O_coord);

		computeForce_Energy(n_atoms, iter, delta_write, box, cutoff_squ, eps, sigma6, Tot_en, coord, forces_on_atom);
	
		

		if(iter%delta_write==0) {

		
			
			write_xyz(coord, n_atoms, iter, box, xyzOut);
			
			
			
			write_Force(n_atoms, iter, forces_on_atom, forceOut);

			fflush(xyzOut);
			
			fflush(forceOut);

		}	

	}



	fclose(xyzOut);
	fclose(forceOut);

}

	



///////////////////////////////////////////////////////////////////////////////////////////////////////

/* SUBROUTINS AND FUNCTIONS */

/* READ CONFIG FILE */
	
void readConfigFile(char *log_FileName, char *vel_FileName, char *traj_FileName, char *force_FileName, double *temp, int *n_atoms, int *n_iter, int *delta_write, double *cutoff, double *dt)
	{

	char buffer[1024];
        char tempBuffer[1024];
        char check[15];
        char *firstWord;	
	
	while (fgets(buffer,1024,stdin) != NULL) {
		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
                if (strncmp(firstWord,"log_File",8)==0) {
		       strcpy(log_FileName,string_secondword(buffer));
               } else if (strncmp(firstWord,"vel_File",8)==0) {
		       strcpy(vel_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"traj_File",9)==0) {
                       strcpy(traj_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"force_File",10)==0) {
		       strcpy(force_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"temperature",11)==0) {
	               *temp = atof(string_secondword(buffer));
	       } else if (strncmp(firstWord,"n_atoms",7)==0) {
	               *n_atoms = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"n_iter",6)==0) {
	               *n_iter = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"delta_write",11)==0) {
	               *delta_write = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"cutoff",6)==0) {
		       *cutoff = atof(string_secondword(buffer));
	       } else if (strncmp(firstWord,"dt",2)==0) {
		       *cutoff = atof(string_secondword(buffer));
	       }
	}
}



////////////////////////////////////////////////////


/* Initilazie Positions */

void initPositions(double **coord, int n_atoms, double *box) {

	double cbrt(double x);          // cube root function
	int iBoxD;                      // integer box dimension
	double fBoxD;                   // float box dimension
	
	int x, y, z;
	double xPos, yPos, zPos;
	int atomCount;

	// determine how many bins to divide the box into
	iBoxD = (int) cbrt((double) n_atoms);	// number of atoms should be an integer to the third power
	if (iBoxD*iBoxD*iBoxD < n_atoms) {
		iBoxD++;
	}

	// determine the size of the bins
	fBoxD = 3.55;        	// seems unitless since multipleied by (x+0.5) gives xPos but I think it gives dimention to particles 
	*box = iBoxD*fBoxD;     // calculate dimension of the box

	// add an atom to each created bin 
	atomCount = 0;
	for(x=0;x<iBoxD;x++) {
		xPos = (x+0.5)*fBoxD;
		for(y=0;y<iBoxD;y++) {
			yPos = (y+0.5)*fBoxD;
			for(z=0; z<iBoxD; z++) {
				if (atomCount < n_atoms) {
					zPos = (z+0.5)*fBoxD;
					coord[atomCount][0]=xPos;
					coord[atomCount][1]=yPos;
					coord[atomCount][2]=zPos;
					atomCount++;
				} else {
					break;
				}
			}

			if (atomCount>=n_atoms) {
				break;
			}
		}
	}
}

/////////////////////////////////////////////////////////

// Compute pervious step coordinates
	
void PrevCoord(double n_atoms, double dt, double **coord, double **atom_vel, double **O_coord) {

	int i, j;
	
	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {
		O_coord[i][j] = coord[i][j] - atom_vel[i][j]*dt;	
		}
	}
}	
//////////////////////////////////////////////////////////////

/* Initialize Velocities */

void initVelocities(int n_atoms, int seed, double mass, double kBT, double **atom_vel) {
	 
	double sumv[3] = {0.0, 0.0, 0.0};
	double msv = 0;
	int i, j;
	double scale_factor;
	double total_number_of_element = n_atoms*3;

	double sqrt(double x);

	srand((unsigned) seed);

	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {		
			// randomly assigned velocities in th range of -0.5 - 0.5
			atom_vel[i][j] = (rand()/((double) RAND_MAX))-0.5;
			//printf("%f  ", atom_vel[i][j]);
			// sum in each coordinate
			sumv[j] = sumv[j] + atom_vel[i][j]/n_atoms;
			// mean square velocity
			msv = msv + (atom_vel[i][j]*atom_vel[i][j])/total_number_of_element;
		}
                //printf("\n");
	}
	scale_factor = sqrt(3.0*kBT/(mass*msv)); 
	// printf("scale factor %f mass %f", kBT, mass);
	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {
			atom_vel[i][j] = (atom_vel[i][j] - sumv[j])*scale_factor;
			//printf("%f  ", atom_vel[i][j]);
		}

	}

}
////////////////////////////////////////////////////

// Sub: Compute Energy and forces

void computeForce_Energy(int n_atoms, int n_iter, int delta_write, double box, double cutoff_squ, double eps, double sigma6, double Tff, double **coord, double **forces_on_atom) {

	 
	int atom1;
        int atom2;
        double r2;
        double r2rev;
        double r6rev;
        double rev6sig6;
        int j;
        double element[3] = {0.0, 0.0, 0.0};
        double ff;

	Tff = 0.0;

	for(atom1=0; atom1<n_atoms; atom1++) {
		for(j=0;j<3;j++) {
			forces_on_atom[atom1][j] = 0.0;
		}
	}
	

	for(atom1=0; atom1<n_atoms-1; atom1++){

		for(atom2=atom1+1;atom2<n_atoms;atom2++){
			r2=0.0;
			
			for(j=0;j<3;j++) {

				element[j]= coord[atom1][j]-coord[atom2][j];
				if(element[j]<-box/2.0){
					element[j]=element[j]+box;
				} else if (element[j]>box/2.0) {
					element[j]=element[j]-box;
				}
				r2+=element[j]*element[j];

			}
//		printf("r2=%f", r2);
		
			if(r2<cutoff_squ) {
				r2rev=1.0/r2;
				r6rev=r2rev*r2rev*r2rev;
				rev6sig6=sigma6*r6rev;
				ff=48*eps*r2rev*rev6sig6*(rev6sig6-0.5);

				for(j=0; j<3; j++) {
					forces_on_atom[atom1][j]+= ff*element[j];	  // forces acting on atoms
					forces_on_atom[atom2][j]+= -ff*element[j];
				}
			}

			
		}
	}

//printf("ff=%f", ff);
}

////////////////////////////////////////////////////

// Sub: Integrations (Propagations)

/* VERLET INTEGRATION STEP */



// compute new coordinates


void CoordUpdate(int n_atoms, double im, double dt, double dt2, double box, double **coord, double **atom_vel, double **forces_on_atom, double **O_coord) {
	
	int i, j;
	double element;
	
	for (i=0; i<n_atoms; i++) {
		for (j=0; j<3; j++) {
			element = 2*coord[i][j] - O_coord[i][j] + im*dt2*forces_on_atom[i][j];
			// wrapping the particles within the box...
			if(element<0) {
				element += box;
			} else if(element > box) {
				element -= box;
			}
			coord[i][j] = element;
			element = 0.0;
		}
	}
}

////////////////////////////////////////////////////

// Sub: Write Output data

// WRITE 3D COORDINATES

void write_xyz(double **coord, int n_atoms, int n_iter, double box, FILE *xyzOut) {

	int i,j,k;
	int atom;
	
	fprintf(xyzOut, "%d \n", n_atoms);
	fprintf(xyzOut, "Step %d , box: %8.3f  %8.3f  %8.3f\n", n_iter, box, box, box);
	for (atom=0;atom<n_atoms;atom++) {
		fprintf(xyzOut, "Ar%12.6f%12.6f%12.6f\n",coord[atom][0],coord[atom][1],coord[atom][2]);

	}
}



// Write Velocities 


void write_Vel(int n_atoms, int n_iter, double **atom_vel, FILE *velOut) {
	int atom;

	fprintf(velOut, "Step %d\n", n_iter);
	for(atom=0; atom<n_atoms; atom++) {
		fprintf(velOut, "Ar%f   %f   %f\n", atom_vel[atom][0], atom_vel[atom][1], atom_vel[atom][2]);
	}
}

// Write Forces 


void write_Force(int n_atoms, int n_iter, double **forces_on_atom, FILE *forceOut) {

	int atom;

	fprintf(forceOut, "Step %d\n", n_iter);
	for(atom=0; atom<n_atoms; atom++) {
		fprintf(forceOut, "Ar%16.6E%16.6E%16.6E\n", forces_on_atom[atom][0], forces_on_atom[atom][1], forces_on_atom[atom][2]);
	}
}

////////////////////////////////////////////////////
//
//
//
//
//
//
