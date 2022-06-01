#include <iostream>
#include "mex.h"
#include <cmath>
#include "voro++.hh"
//#include "matrix.h"

using namespace std;
using namespace voro;

/*
  This is the function that we will write our mex file to interface with
*/

int powerfn(double, double, double, double, double, double, const int&, double*, double*, double*, double*, double*,mxArray*,bool);

/*
  Function definition
*/

int powerfn(double x_min,double x_max,double y_min,double y_max,double z_min,double z_max,const int& N,double* X,double* w, double* vol, double* trans,double* xc,mxArray* VFN,bool period){
  
  /* Inputs are x_min,x_max,y_min,y_max,z_min,z_max which are coordinates of the box
  N is the number of seeds/generators
  X is a pointer to the coordinates of the seeds/generators
  w is a pointer to the weights
  vol is a pointer to the volumes
  trans is a pointer to the 2nd moments
  xc is a pointer to the cell centroids
  VFN is a pointer to a mex array
  period is a boolean controlling the periodicity */
  
  // The numbers n_x,n_y,n_z are related to efficiency of the calculations in making a periodic cell
  
  const int n_x=6,n_y=6,n_z=6;
  container_poly* con;
  con=new container_poly(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,period,period,period,16);
  
  // Find the minimum weight, this is because voro++ uses particle radius rather than the weight
  double wmin=w[0];
  for(int i=1;i<N;i++){
    if(w[i]<wmin){wmin=w[i];}
  }

  // Add particles into the container, the particle radius is related to the weights after we subtract the minimum weight r=sqrt(w-w_min)
  for(int i=0;i<N;i++){
    double r=sqrt(w[i]-wmin);
    con->put(i,X[i],X[i+N],X[i+2*N],r);
  }

  // Calculate the cells
  c_loop_all cla(*con);
  voronoicell_neighbor c;

  int Np=0;
  if(cla.start()) do if (con->compute_cell(c,cla)) {
	int id;
       	// Get the position and ID information for the particle
	// currently being considered by the loop. Ignore the radius
	// information.

	// Obtain the centroid in local coordinates (relative to the seed location)
	double cx;double cy;double cz;
	c.centroid(cx,cy,cz);

	// Get the location of the seed (in the periodic case if seeds lie outside the container, they are first mapped to an equivalent location inside the container, in the non-periodic case the seed location is the same as the input)
	
	double x;double y;double z;double r;
	cla.pos(id,x,y,z,r);

	// Get the cell volume
	double V=c.volume();
	// Get the second moment (relative to the seed location) of the cell
	double T=c.transportcost();

	// Populate the arrays for return variables
	vol[id]=V;
	trans[id]=T;
	
	// Return the centroid relative to the (remapped) seed location
	xc[id]=cx+x;xc[id+N]=cy+y;xc[id+2*N]=cz+z;

	// Increment the number of particles
	Np++;

	// Create a vector to store the vertices of the cell
	std::vector<double> vs;
	c.vertices(x,y,z,vs);
	// Calculate the number of vertices of the cell
	int Nv=vs.size()/3;

	// Create a double matrix to enter into the cell array to return to MatLab, to contain the vertices
	mxArray* tmp;double* tmpV;
	tmp=mxCreateDoubleMatrix(Nv,3,mxREAL);
	tmpV=mxGetPr(tmp);

	// Transfer the data from vs (the vector of vertices) to the output
	for(int j=0;j<Nv;j++){
	  tmpV[j]=vs[3*j];tmpV[j+Nv]=vs[3*j+1];tmpV[j+2*Nv]=vs[3*j+2];
	}

	// Obtain the single index to each cell entry and assign the data
	mwIndex subs[2];
	// The row is the seed/generator/particle id
	subs[0]=id;
	// The column for the vertex data is the 1st column (zero indexing so equal to 0)
	subs[1]=0;
	mwIndex ind=mxCalcSingleSubscript(VFN,(mwSize)2,subs);
	// Assign the vertex data for cell 'id'
	mxSetCell(VFN,ind,tmp);

	/* Face data
	   The column for the face data is the 2nd column (zero indexing so equal to 1) */
	subs[1]=1;
	ind=mxCalcSingleSubscript(VFN,(mwSize)2,subs);
	
	// All the vertices of all the faces as an integer list
	std::vector<int> fv;
	c.face_vertices(fv);
	int Nfv=fv.size();

	// First work out how many faces there are
	int nf=c.number_of_faces();

	// Now we extract the individual faces
	mxArray* face_cell_array;
	face_cell_array=mxCreateCellMatrix((mwSize)nf,(mwSize)1);
	int j=0;
	
	for(int fi=0;fi<nf;fi++){
	  /*
	    For each face we extract the face indices (remember that in C they are zero indexed
	    but in MatLab we want indexing from 1 */

	  // Make appropriate pointers
	  mxArray* tmpj;
	  int* tmpjFV;

	  // Find the number of vertices for the face
	  int l=fv[j];

	  // Create an integer array to contain the face vertices of the current face
	  tmpj=mxCreateNumericMatrix(1,l,mxINT32_CLASS,mxREAL);
	  tmpjFV=(int*) mxGetData(tmpj);

	  // Write the face vertex indices
	  for(int k=0;k<l;k++){
	    tmpjFV[k]=fv[j+k+1]+1;
	  }
	  j=j+l+1;
	  // Assign the face vertices to the cell array of face vertices
	  mxSetCell(face_cell_array,fi,tmpj);
	}

	// Assign the set of face vertices for each face to the output
	mxSetCell(VFN,ind,face_cell_array);
	
	/* Neighbour data
	   The column for the neighbour data is the 3rd column (zero indexing so equal to 2) */
	
	subs[1]=2;
	ind=mxCalcSingleSubscript(VFN,(mwSize)2,subs);

	std::vector<int> cell_neighbors;
	c.neighbors(cell_neighbors);
	int Ncn=cell_neighbors.size();

	mxArray* tmpic;
	int* tmpicn;
	tmpic=mxCreateNumericMatrix(Ncn,1,mxINT32_CLASS,mxREAL);
	tmpicn=(int*) mxGetData(tmpic);

	//	mexPrintf("Cell %d has %d neighbours\n",id,Ncn);
	for(int j=0;j<Ncn;j++){
	  // MatLab uses indexing from 1, so we add 1 to the cell id to correctly record the neighbour information
	  tmpicn[j]=cell_neighbors[j]+1;
	}

	// Assign the neighbour information to the output
	mxSetCell(VFN,ind,tmpic);
	
      } while (cla.inc());
     
  // Make sure we delete any 'new' variables, to free up memory
    
  delete con;
  return Np;
}

/*
  The main program
*/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){

  // Cell array variables
  
  mxArray* verts_faces_neighs; // This is the name of the cell array

  mwSize ndim=2; // This is the number of dimensions of the cell array
  mwSize* dims; // This is an array of the dimensions of the cell array

  // Default box size [0,1]x[0,1]x[0,1]
  double xmin=0;double ymin=0;double zmin=0;
  double xmax=1;double ymax=1;double zmax=1;

  double* box;
  double* X;
  double* W;

  // Number of seeds/generators
  int N;
  bool period=false;

  bool Wflag=false;

  if(nrhs==0){
    // If there are no inputs then user warned that they must at least specify the locations of the generators
    mexErrMsgIdAndTxt("mexPD:InputArguments","You must specify at least the locations of the generators.");
  }
  else if(nrhs==1){
    /*
      In the case of only one input argument it should be a double array containing generator locations. Size and shape should be checked to see that valid list of generator locations has been provided.
    */
    mxClassID cat1;
    cat1=mxGetClassID(prhs[0]);
    if(cat1!=mxDOUBLE_CLASS){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying one argument, it must be a numeric double array consisting of generator locations");
    }
    X=mxGetPr(prhs[0]);
    int Xrows = mxGetM(prhs[0]);
    int Xcols = mxGetN(prhs[0]);
    if(Xcols!=3){
      mexErrMsgIdAndTxt("mexPD:Xshape","The locations of generators must be an Nx3 array");
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
      In this case we could miss out the box and weights, in which case it should be assumed a unit cube and zero weights, or we could miss out the box and periodic flag in which case it should be assumed a unit cube and not periodic, or we could miss out the weights and periodic flag in which case the weights should be assumed zero and it should be assumed not periodic
    */

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    // Check if the first argument is a double array
    if(cat1!=mxDOUBLE_CLASS){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, the first must be a double numeric array, either the box size or generator locations");
    }
    else {
      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);

      // Check that the first argument has the correct number of columns
      if(cols1!=3){
	mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, the first must be a double numeric array, either the box size or generator locations. The box size is a 1x3 array and the generator locations is an Nx3 array");
      }
      else{
	// Here we can be sure that the first argument is a double array, and is either a 1x3 or Nx3 array.
	if(cat2==mxLOGICAL_CLASS){
	  // If the second argument is a logical then we assume that the first argument contains the generator locations and take the box size to be the default, we also set the weights to be zero
	  period=mxIsLogicalScalarTrue(prhs[1]);
	  X=mxGetPr(prhs[0]);
	  int Xrows=mxGetM(prhs[0]);
	  N=Xrows;
	  
	  W=new double[N];Wflag=true;
	  for(int i=0;i<N;i++){
	    W[i]=0.0;
	  }
	} else if(cat2==mxDOUBLE_CLASS){
	  // If the second argument is a double array then we are now in one of the cases of box size, generators being specified or generators and weights being specified
	  int rows2=mxGetM(prhs[1]);
	  int cols2=mxGetN(prhs[1]);
	
	  if(rows1==1 && cols2==3){
	    // Here we have box size and generators being specified, we assume no periodicity and that the weights are all zero
	    box=mxGetPr(prhs[0]);
	    xmax=box[0];
	    ymax=box[1];
	    zmax=box[2];
	  
	    X=mxGetPr(prhs[1]);
	    int Xrows = mxGetM(prhs[1]);
	    N=Xrows;
	    W=new double[N];Wflag=true;
	    for(int i=0;i<N;i++){
	      W[i]=0.0;
	    }
	  } else if(cols2==1){
	    // Here we are in the case where the second argument is the weights
	    if(rows1!=rows2){
	      // The number of generator locations and the number of weights must be the same
	      mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the generator locations and weights must be the same size");
	    } else
	      {
		X=mxGetPr(prhs[0]);
		int Xrows = mxGetM(prhs[0]);
		N=Xrows;
		W=mxGetPr(prhs[1]);
		int Wrows = mxGetM(prhs[1]);
	      }
	  }
	  else {
	    mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, the first must be a double numeric array, either the box size or generator locations. The box size is a 1x3 array and the generator locations is an Nx3 array");
	  }
	}
	else{
	  // In this case the second argument is neither a logical or a double array
	  mexErrMsgIdAndTxt("mexPD:InputArguments","When specifying two arguments, they should be one of the following combinations (box,x), (x,w), (x,periodic)");
	}
      }
    }
  }
  else if(nrhs==3){
    /*
      In this case we can miss out the box size, in which the code should assume a unit cube, or we could miss out the weights, in which case the code should assume zero weights, or we could miss out the periodic flag, in which case it should be assumed not periodic.
    */

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    mxClassID cat3;
    cat3=mxGetClassID(prhs[2]);

    if(cat3==mxLOGICAL_CLASS){
      // In this case we have provided the third argument as a logical and so the periodicity should be set equal to the value provided
      period=mxIsLogicalScalarTrue(prhs[2]);

      // We must also check that the first two arguments are double arrays of the right size

      if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS)){
	mexErrMsgIdAndTxt("mexPD:InputArguments","If the third argument is specified as a logical variable then the first two must be double arrays");
      }
      else{
	int rows1=mxGetM(prhs[0]);
	int cols1=mxGetN(prhs[0]);
	
	int rows2=mxGetM(prhs[1]);
	int cols2=mxGetN(prhs[1]);
	
	if(rows1==1 && cols1==3 && cols2==3){
	  // In this case we have specified a box size and some number of generator locations and must set the weights to be zero
	  box=mxGetPr(prhs[0]);
	  xmax=box[0];
	  ymax=box[1];
	  zmax=box[2];
	  
	  X=mxGetPr(prhs[1]);
	  int Xrows=mxGetM(prhs[1]);
	  N=Xrows;
	  W=new double[N];Wflag=true;
	  for(int i=0;i<N;i++){
	    W[i]=0.0;
	  }
	}
	else if(cols1==3 && cols2==1){
	  // In this case we have specified generators and weights and so must check that the number of generators and weights is the same
	  if(rows1!=rows2){
	    mexErrMsgIdAndTxt("mexPD:ArraySize","The arrays containing the generator locations and weights must be the same size");
	  }
	  else{
	    X=mxGetPr(prhs[0]);
	    int Xrows = mxGetM(prhs[0]);
	    N=Xrows;
	    W=mxGetPr(prhs[1]);
	    int Wrows = mxGetM(prhs[1]);
	  }
	} else{
	  // In this case the size of the double arrays is incompatible and cannot be decided what the data is
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The first two arguments must either be box dimensions (1x3 array) and generator locations (Nx3 array), or generator locations (Nx3 array) and weights (Nx1 array) ");
	}
      }
    } else{
      // This is the case where the third argument is not a logical type, so we must check that all the other arguments are double arrays

      if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS && cat3==mxDOUBLE_CLASS)){
	mexErrMsgIdAndTxt("mexPD:InputArguments","When specifiying three arguments these can be one of (box,x,periodic) or (x,w,periodic) or (box,x,w) where box is a 1x3 double array, x is an Nx3 double array, w is an Nx1 double array and periodic is a boolean variable");
      }
      else{
	// In this case all the three inputs are double arrays and we must check the sizes

	int rows1=mxGetM(prhs[0]);
	int cols1=mxGetN(prhs[0]);
	
	int rows2=mxGetM(prhs[1]);
	int cols2=mxGetN(prhs[1]);

	int rows3=mxGetM(prhs[2]);
	int cols3=mxGetN(prhs[2]);

	if(!(rows1==1 && cols1==3 && cols2==3 && cols3==1)){
	  mexErrMsgIdAndTxt("mexPD:ArraySize","The three arguments must be box dimensions (1x3 array), generator locations (Nx3 array), and weights (Nx1 array) ");
	}
	else{
	  box=mxGetPr(prhs[0]);
	  xmax=box[0];
	  ymax=box[1];
	  zmax=box[2];
	    
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
  } else if(nrhs==4){
    // In this case we have specified four arguments, they must be box size, generator locations, weights and a periodicity flag

    mxClassID cat1;
    mxClassID cat2;
    
    cat1=mxGetClassID(prhs[0]);
    cat2=mxGetClassID(prhs[1]); 

    mxClassID cat3;
    cat3=mxGetClassID(prhs[2]);

    mxClassID cat4;
    cat4=mxGetClassID(prhs[3]);

    if(!(cat1==mxDOUBLE_CLASS && cat2==mxDOUBLE_CLASS && cat3==mxDOUBLE_CLASS && cat4==mxLOGICAL_CLASS)){
      mexErrMsgIdAndTxt("mexPD:InputArguments","When specifiying four arguments these should be (box,x,w,periodic) where box is a 1x3 double array, x is an Nx3 double array, w is an Nx1 double array and periodic is a boolean variable");
    }
    else{

      // Decide whether periodic or not
      period=mxIsLogicalScalarTrue(prhs[3]);

      int rows1=mxGetM(prhs[0]);
      int cols1=mxGetN(prhs[0]);
	
      int rows2=mxGetM(prhs[1]);
      int cols2=mxGetN(prhs[1]);

      int rows3=mxGetM(prhs[2]);
      int cols3=mxGetN(prhs[2]);

      if(!(rows1==1 && cols1==3 && cols2==3 && cols3==1)){
	mexErrMsgIdAndTxt("mexPD:ArraySize","The three arguments must be box dimensions (1x3 array), generator locations (Nx3 array), and weights (Nx1 array) ");
      }
      else{

	box=mxGetPr(prhs[0]);
	xmax=box[0];
	ymax=box[1];
	zmax=box[2];
	  
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
  
  /*
    Here we create the cell array that contains the vertices, face information and neighbour information of each Laguerre cell
    The cell array will be an Nc x 3 array, the first entry is the vertices of the cell
  */
  
  verts_faces_neighs=mxCreateCellMatrix((mwSize)N,(mwSize)3);
  
  // Output variable pointers
  double* XC;
  double* V;
  double* T;

  // Set the first return value to be the volumes of the cells
  plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
  V=mxGetPr(plhs[0]);

  // Set the second return value to be the transport costs of the cells
  plhs[1]=mxCreateDoubleMatrix(N,1,mxREAL);
  T=mxGetPr(plhs[1]);

  // Set the third reutrn value to be the centroids of the cells
  plhs[2]=mxCreateDoubleMatrix(N,3,mxREAL);
  XC=mxGetPr(plhs[2]);

  // Set the fourth return value to be a cell array containing vertices, faces and neighbours
  plhs[3]=verts_faces_neighs;

  int NP;
  NP=powerfn(xmin,xmax,ymin,ymax,zmin,zmax,N,X,W,V,T,XC,verts_faces_neighs,period);

  // Clean up

  if(Wflag){
    delete [] W;
  }
}

