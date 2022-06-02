#include <iostream>
#include "mex.h"
#include <cmath>
#include "voro++_2d.hh"
//#include "matrix.h"

using namespace std;
using namespace voro;

/*
  This is the signature of the function that we will write our mex file to interface with
*/

int powerfn(double, double, double, double, const int&, double*, double*, double*, double*, double*, bool, bool);

/*
  Function definition
*/

int powerfn(double bx_xmin,double bx_xmax,double bx_ymin,double bx_ymax,const int& N,double* X,double* w,double* vol,double* xc, double* t, bool periodx, bool periody){

  /* Inputs are x_min,x_max,y_min,y_max which are coordinates of the box
     N is the number of seeds/generators
     X is a pointer to the coordinates of the seeds/generators
     w is a pointer to the weights
     vol is a pointer to the volumes
     t is a pointer to the transport costs
     xc is a pointer to the cell centroids
     periodx is a boolean controlling the periodicity in the x-direction
     periody is a boolean controlling the periodicity in the y-direction */
  
  // The numbers n_x,n_y are related to efficiency of the calculations in making a periodic cell

  const int n_x=10,n_y=10;
  const double rr=1.0;

  /*
    Initially the box is set to be [bx_xmin,bx_xmax] x [bx_ymin,bx_ymax]
    We allow these positions to be modified in case seeds lie outside the box
  */
  
  double xmin=bx_xmin,xmax=bx_xmax;
  double ymin=bx_ymin,ymax=bx_ymax;

  /* Construct planes that are the actual locations of the box walls */
  
  wall_plane_2d* left;
  wall_plane_2d* right;
  wall_plane_2d* top;
  wall_plane_2d* bottom;

  // Create the walls, these are not necessarily added to the container unless the domain is not periodic in the relevant direction
  
  left=new wall_plane_2d(-1.,0.,bx_xmin,-99);
  right=new wall_plane_2d(1.,0.,bx_xmax,-97);
  top=new wall_plane_2d(0.,1.,bx_ymax,-98);
  bottom=new wall_plane_2d(0.,-1.,bx_ymin,-96);

  /* Pointer to the container */
  container_poly_2d* con;
  
  /* Find the minimum weight, this is used to define the 'radius' of the particles in voro++ */

  double wmin=w[0];
  for(int i=1;i<N;i++){
    if(w[i]<wmin){wmin=w[i];}
  }

  /*
    If we are not periodic then we must allow for the fact that seeds can lie outside the domain
    To deal with this we first find the maximum and minimum values of the the x and y coordinates of the seeds
    then we specify a box that is large enough to contain these seeds, we will later clip the cells generated
    to the original box
  */
  
  if(!periodx){
    // Find the maximum and minimum x coordinates of all the points
    for(int i=1;i<N;i++){
      if(X[i]<xmin){xmin=X[i];}
      else if(X[i]>xmax){xmax=X[i];}
    }

    // Dilate the computational box so that the seed locations are not exactly on the perimeter
    xmin=0.5*((xmin+xmax)-(xmax-xmin)*(1.+rr));
    xmax=0.5*((xmin+xmax)+(xmax-xmin)*(1.+rr));
  }

  if(!periody){
    // Find the maximum and minimum y coordinates of all the points    
    for(int i=1;i<N;i++){
      if(X[i+N]<ymin){ymin=X[i+N];}
      else if(X[i+N]>ymax){ymax=X[i+N];}
    }

    // Dilate the box so the seed locations are not on the perimeter
    ymin=0.5*((ymin+ymax)-(ymax-ymin)*(1.+rr));
    ymax=0.5*((ymin+ymax)+(ymax-ymin)*(1.+rr));
  }    
      
  // Specify the container that is large enough to fit all the seeds
  con=new container_poly_2d(xmin,xmax,ymin,ymax,n_x,n_y,periodx,periody,16);
  
  // Depending on periodicity add the walls
  if(!periodx){
    con->add_wall(left);
    con->add_wall(right);
  }

  if(!periody){
    con->add_wall(bottom);
    con->add_wall(top);
  }
  
  /* Add particles into the container
     the particle radius is related to the weights after we subtract the minimum weight r=sqrt(w-w_min) */
  for(int i=0;i<N;i++){
    double r=sqrt(w[i]-wmin);
    con->put(i,X[i],X[i+N],r);
  }

  // Calculate the cells
  c_loop_all_2d cla(*con);
  voronoicell_2d c;

  int Np=0;
  if(cla.start()) do if (con->compute_cell(c,cla)) {
	int id;
	
	// Get the position and ID information for the particle
	// currently being considered by the loop.

	// Obtain the centroid in local coordinates (relative to the seed location)
	double cx;double cy;
	c.centroid(cx,cy);

	// Get the location of the seed (in the periodic case if seeds lie outside the container, they are first mapped to an equivalent location inside the container, in the non-periodic case the seed location is the same as the input)
	
	double x;double y;double r;
	cla.pos(id,x,y,r);

	// Get the cell area
	double A=c.area();
	// Get the second moment (relative to the seed location) of the cell
	double T=c.transport_cost();

	// Populate the arrays for return variables
	vol[id]=A;
	t[id]=T;
	
	// Return centroid relative to the (remapped) seed location)
	xc[id]=cx+x;xc[id+N]=cy+y;

	// Increment the number of particles
	Np++;
      } while (cla.inc());
        
  // Make sure we delete any 'new' variables, to free up memory and avoid memory leaks
  delete left;
  delete right;
  delete top;
  delete bottom;

  delete con;
  return Np;
}

/*
  The main program
*/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){

  // The default box size is [0,1]x[0,1]
  double xmin=0;double ymin=0;
  double xmax=1;double ymax=1;

  double* box;
  double* X;
  double* W;

  // Number of seeds/generators
  int N;

  // Periodicity flags
  bool periodx=false;
  bool periody=false;
  
  bool Wflag=false;

  if(nrhs==0){
    // If there are no inputs then user warned that they must at least specify the locations of the generators
    mexErrMsgIdAndTxt("mexPD:InputArguments","You must specify at least the locations of the generators.");
  }
  else if(nrhs==1){
    /*
      In this case we have the following inputs
      1. generator locations - x (Nx2)
    */
    
    mxClassID cat1;
    cat1=mxGetClassID(prhs[0]);
    if(cat1!=mxDOUBLE_CLASS){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying one argument, it must be a numeric double array consisting of generator locations");
    }
    X=mxGetPr(prhs[0]);
    int Xrows = mxGetM(prhs[0]);
    int Xcols = mxGetN(prhs[0]);
    if(Xcols!=2){
      mexErrMsgIdAndTxt("mexPD:Xshape","The locations of generators must be an Nx2 array");
    }

    N=Xrows;
    // Given only a list of generator locations we must make the weights all zero. Periodicity is assumed to be false and the box size is the default
    W=new double[N];Wflag=true;
    for(int i=0;i<N;i++){
      W[i]=0.0;
    }
  }
  else if(nrhs==2){
    /*
      In this case we have either of the following inputs
      1. box corners and generator locations - bx (1x4) and x (Nx2)
      2. generator locations and generator weights - x (Nx2) and w (Nx1)
    */

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    // Check if the first argument is a double array
    if(cat1!=mxDOUBLE_CLASS){ // Check if the first argument is NOT a double array
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, the possible choices are box corners (1x4 array) and generator locations (Nx2) array or generator locations (Nx2 array) and weights (Nx1 array)");
    }
    else { // The first argument is a double array
      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);

      if(rows1==1 && cols1==4){ // Check that the first argument is a 1x4 array
	// In this case the first input is the box size
	box=mxGetPr(prhs[0]);
	xmin=box[0];ymin=box[1];
	xmax=box[2];ymax=box[3];

	if(xmax<xmin || ymax<ymin){ // Now check that the box vertices make sense
	  mexErrMsgIdAndTxt("mexPD:InputArguments", "Specify the box dimensions in the form [xmin ymin xmax ymax] with xmax>xmin and ymax>ymin");
	}
	else{ // At this stage we can be sure that the first argument is a valid specification of the box
	  if(cat2==mxDOUBLE_CLASS){  // Check the second argument
	    // Now we now the second argument is a double array, it should be Nx2 to be a valid second input (the generator locations)
	    int rows2=mxGetM(prhs[1]);
	    int cols2=mxGetN(prhs[1]);
	    N=rows2;
	  
	    if(cols2!=2){
	      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying the inputs as box plus generator locations the 2nd argument must be an Nx2 array");
	    }
	    else{
	      X=mxGetPr(prhs[1]);	  
	      N=rows2;

	      // Assume zero weights
	      W=new double[N];Wflag=true;
	      for(int i=0;i<N;i++){
		W[i]=0.0;
	      }
	    }
	  }
	  else{
	    mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments with the first as box vertices, the second must be an Nx2 double array"); 
	  }
	}
      }
      else if(cols1==2){
	// Now the first argument is an Nx2 array which must be the generator locations
	// Now check the second argument
	
	if(cat2==mxDOUBLE_CLASS){
	  int rows2=mxGetM(prhs[1]);
	  int cols2=mxGetN(prhs[1]);
	
	  if(cols2==1){
	    // Here we are in the case where the second argument is the weights
	    if(rows1!=rows2){
	      // The number of generator locations and the number of weights must be the same
	      mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the generator locations and weights must be the same size");
	    }
	    else{
	      X=mxGetPr(prhs[0]);
	      N=rows1;
	      W=mxGetPr(prhs[1]);
	    }
	  }
	  else{
	    mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, when the first is an Nx2 array of generator locations, the second must be an Nx1 array of weights");
	  }
	}
	else{
	  // In this case the second argument is not a double array
	  mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, they should be one of the following combinations (box,x), (x,w)");
	}
      }
      else{
	mexErrMsgIdAndTxt("mexPD:InputArguments", "When specifying two arguments they should be one of the following combinations (box,x) where box is a 1x4 double array and x is an Nx2 double array, or (x,w) where x is an Nx2 double array and w is an Nx1 double array");
      }
    }
  }
  else if(nrhs==3){

    /*
      In this case we have the following possible inputs
      1. box corners, generator locations, weights - bx (1x4), x (Nx2) and w (Nx1)
      2. box corners, generator locations and a periodic flag - bx (1x4), x (Nx2) and periodic (bool) 
      3. generator locations, generator weights and a periodic flag - x (Nx2), w (Nx1) and periodic (bool)
    */

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    mxClassID cat3;
    cat3=mxGetClassID(prhs[2]);

    if(cat3==mxLOGICAL_CLASS){
      /*
	In this case we have provided the third argument as a logical and so the periodicity in both directions should be
	set equal to the value provided
      */
      periodx=mxIsLogicalScalarTrue(prhs[2]);
      periody=periodx;
      
      // We must also check that the first two arguments are double arrays of the right size

      if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS)){
	mexErrMsgIdAndTxt("mexPD:InputArguments","If the third argument is specified as a logical variable then the first two must be double arrays");
      }
      else{
	int rows1=mxGetM(prhs[0]);
	int cols1=mxGetN(prhs[0]);
	
	int rows2=mxGetM(prhs[1]);
	int cols2=mxGetN(prhs[1]);
	
	if(rows1==1 && cols1==4&& cols2==2){
	  // In this case we have specified a box size and some number of generator locations and must set the weights to be zero
	  box=mxGetPr(prhs[0]);
	  xmin=box[0];ymin=box[1];
	  xmax=box[2];ymax=box[3];

	  if(xmax<xmin || ymax<ymin){
	    mexErrMsgIdAndTxt("mexPD:InputArguments", "Specify the box dimensions in the form [xmin ymin xmax ymax] with xmax>xmin and ymax>ymin");
	  }

	  
	  X=mxGetPr(prhs[1]);
	  N=rows2;

	  // Set the weights to be zero
	  W=new double[N];Wflag=true;
	  for(int i=0;i<N;i++){
	    W[i]=0.0;
	  }
	}
	else if(cols1==2 && cols2==1){
	  // In this case we have specified generators and weights and so must check that the number of generators and weights is the same
	  if(rows1!=rows2){
	    mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the generator locations and weights must be the same size");
	  }
	  else{
	    X=mxGetPr(prhs[0]);
	    N=rows1;

	    W=mxGetPr(prhs[1]);
	  }
	} else{
	  // In this case the size of the double arrays is incompatible and cannot be decided what the data is
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The first two arguments must either be box coordinates (1x4 array) and generator locations (Nx2 array), or generator locations (Nx2 array) and weights (Nx1 array) ");
	}
      }
    } // This is the case where the third argument is not a logical type, so we must check that all the other arguments are double arrays
    else if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS && cat3==mxDOUBLE_CLASS)){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifiying three arguments these can be one of (box,x,periodic) or (x,w,periodic) or (box,x,w) where box is a 1x2 double array, x is an Nx2 double array, w is an Nx1 double array and periodic is a boolean variable");
    }
    else{
      // In this case all the three inputs are double arrays and we must check the sizes

      // Box
      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);

      // X
      int rows2=mxGetM(prhs[1]);
      int cols2=mxGetN(prhs[1]);

      // W
      int rows3=mxGetM(prhs[2]);
      int cols3=mxGetN(prhs[2]);

      if(!(rows1==1 && cols1==4 && cols2==2 && cols3==1)){
	mexErrMsgIdAndTxt("mexPD:ArraySize","The three arguments must be box coordinates (1x4 array), generator locations (Nx2 array), and weights (Nx1 array) ");
      }
      else{
	box=mxGetPr(prhs[0]);
	xmin=box[0];ymin=box[1];
	xmax=box[2];ymax=box[3];

	if(xmax<xmin || ymax<ymin){
	  mexErrMsgIdAndTxt("mexPD:InputArguments", "Specify the box dimensions in the form [xmin ymin xmax ymax] with xmax>xmin and ymax>ymin");
	}
	  
	X=mxGetPr(prhs[1]);
	int Xrows = mxGetM(prhs[1]);
	    
	W=mxGetPr(prhs[2]);
	int Wrows = mxGetM(prhs[2]);
	    
	// Error check is Xrows=Wrows?
	if(Xrows!=Wrows){
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the seed locations and weights must be the same size");
	}
	N=Xrows;
      }
    }
  } else if(nrhs==4){
    // In this case we have specified four arguments, they must be box size, generator locations, weights and a periodicity flag for both directions

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    mxClassID cat3;
    cat3=mxGetClassID(prhs[2]);

    mxClassID cat4;
    cat4=mxGetClassID(prhs[3]);

    if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS && cat3==mxDOUBLE_CLASS && cat4==mxLOGICAL_CLASS)){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifiying four arguments these should be (box,x,w,periodic) where box is a 1x2 double array, x is an Nx2 double array, w is an Nx1 double array and periodic is a boolean variable");
    }
    else{

      // Decide whether periodic or not
      periodx=mxIsLogicalScalarTrue(prhs[3]);
      periody=periodx;
      
      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);
	
      int rows2=mxGetM(prhs[1]);
      int cols2=mxGetN(prhs[1]);

      int rows3=mxGetM(prhs[2]);
      int cols3=mxGetN(prhs[2]);

      if(!(rows1==1 && cols1==2 && cols2==2 && cols3==1)){
	mexErrMsgIdAndTxt("mexPD:ArraySize","The first three arguments must be box dimensions (1x2 array), generator locations (Nx2 array), and weights (Nx1 array) ");
      }
      else{

	box=mxGetPr(prhs[0]);
	xmax=box[0];
	ymax=box[1];
	  
	X=mxGetPr(prhs[1]);
	int Xrows = mxGetM(prhs[1]);
	  
	W=mxGetPr(prhs[2]);
	int Wrows = mxGetM(prhs[2]);
	 
	// Error check is Xrows=Wrows?
	if(Xrows!=Wrows){
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the seed locations and weights must be the same size");
	}
	N=Xrows;
      }
    }
  }  else if(nrhs==5){
    // In this case we have specified five arguments, they must be box size, generator locations, weights and two periodicity flags

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    mxClassID cat3;
    cat3=mxGetClassID(prhs[2]);

    mxClassID cat4;
    cat4=mxGetClassID(prhs[3]);

    mxClassID cat5;
    cat5=mxGetClassID(prhs[4]);
    
    if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS && cat3==mxDOUBLE_CLASS && cat4==mxLOGICAL_CLASS && cat5==mxLOGICAL_CLASS)){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifiying five arguments these should be (box,x,w,periodic_x,periodic_y) where box is a 1x2 double array, x is an Nx2 double array, w is an Nx1 double array and periodic_x and periodic_y are boolean variables");
    }
    else{

      // Decide whether periodic or not
      periodx=mxIsLogicalScalarTrue(prhs[3]);
      periody=mxIsLogicalScalarTrue(prhs[4]);
      
      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);
	
      int rows2=mxGetM(prhs[1]);
      int cols2=mxGetN(prhs[1]);

      int rows3=mxGetM(prhs[2]);
      int cols3=mxGetN(prhs[2]);

      if(!(rows1==1 && cols1==2 && cols2==2 && cols3==1)){
	mexErrMsgIdAndTxt("mexPD:ArraySize","The first three arguments must be box dimensions (1x2 array), generator locations (Nx2 array), and weights (Nx1 array) ");
      }
      else{

	box=mxGetPr(prhs[0]);
	xmax=box[0];
	ymax=box[1];
	  
	X=mxGetPr(prhs[1]);
	int Xrows = mxGetM(prhs[1]);
	  
	W=mxGetPr(prhs[2]);
	int Wrows = mxGetM(prhs[2]);
	 
	// Error check is Xrows=Wrows?
	if(Xrows!=Wrows){
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the seed locations and weights must be the same size");
	}
	N=Xrows;
      }
    }
  }

  // At this stage we have parsed all the arguments and now can sensibly process the arguments to make the output data

  // Output variables

  double* XC;
  double* V;
  double* T;

  plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
  V=mxGetPr(plhs[0]);

  plhs[1]=mxCreateDoubleMatrix(N,1,mxREAL);
  T=mxGetPr(plhs[1]);

  plhs[2]=mxCreateDoubleMatrix(N,2,mxREAL);
  XC=mxGetPr(plhs[2]);

  int NP;
  NP=powerfn(xmin,xmax,ymin,ymax,N,X,W,V,XC,T,periodx,periody);

  // Clean up

  if(Wflag){
    delete [] W;
  }

}

