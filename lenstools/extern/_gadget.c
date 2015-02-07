/*Python wrapper module for gadget snapshot reading;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _gadget and it defines the methods below (see docstrings)
*/

#include <stdio.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "gadget.h"
#include "grid.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for reading Gadget2 snapshots";
static char getHeader_docstring[] = "Reads the header of a Gadget2 snapshot";
static char getPosVel_docstring[] = "Gets the positions or velocities of the particles in a Gadget2 snapshot";
static char getID_docstring[] = "Gets the 4 byte int particles IDs from the Gadget2 snapshot";
static char write_docstring[] = "Writes the particles information to a Gadget snapshot, with a proper header";
static char grid3d_docstring[] = "Put the snapshot particles on a regularly spaced grid";
static char adaptive_docstring[] = "Put the snapshot particles on a regularly spaced grid using adaptive smoothing";

//Method declarations
static PyObject *_gadget_getHeader(PyObject *self,PyObject *args);
static PyObject *_gadget_getPosVel(PyObject *self,PyObject *args);
static PyObject *_gadget_getID(PyObject *self,PyObject *args);
static PyObject * _gadget_write(PyObject *self,PyObject *args);
static PyObject * _gadget_grid3d(PyObject *self,PyObject *args);
static PyObject * _gadget_adaptive(PyObject *self,PyObject *args);

//_gadget method definitions
static PyMethodDef module_methods[] = {

	{"getHeader",_gadget_getHeader,METH_VARARGS,getHeader_docstring},
	{"getPosVel",_gadget_getPosVel,METH_VARARGS,getPosVel_docstring},
	{"getID",_gadget_getID,METH_VARARGS,getID_docstring},
	{"write",_gadget_write,METH_VARARGS,write_docstring},
	{"grid3d",_gadget_grid3d,METH_VARARGS,grid3d_docstring},
	{"adaptive",_gadget_adaptive,METH_VARARGS,adaptive_docstring},
	{NULL,NULL,0,NULL}

} ;

//_gadget constructor
PyMODINIT_FUNC init_gadget(void){

	PyObject *m = Py_InitModule3("_gadget",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

//getHeader() implementation
static PyObject *_gadget_getHeader(PyObject *self,PyObject *args){

	struct io_header_1 header;
	PyObject *file_obj;
	int k;

	//Header contents
	int NumPart=0,NumPartFile=0,Ngas=0,Ngas_file=0,Nwithmass=0,Nwithmass_file=0;
	
	//Interpret the tuple of arguments (there should be only one: the file descriptor)
	if(!PyArg_ParseTuple(args,"O",&file_obj)){
		return NULL;
	}

	//Get a file pointer out of the file object
	FILE *fp = PyFile_AsFile(file_obj);
	PyFile_IncUseCount((PyFileObject *)file_obj);

	//Read in the header
	int endianness = getHeader(fp,&header);

	//Release the file object
	PyFile_DecUseCount((PyFileObject *)file_obj);

	//Exit if there was a problem (endianness check shall return 1 or 0, no exceptions; if not, the snapshot has the wrong format)
	if(endianness==-1){
		PyErr_SetString(PyExc_IOError,"Not a valid Gadget2 snapshot!");
		return NULL;
	}

	//Count the number of particles of each kind and total

	//Allocate resources
	npy_intp Nkinds[] = {(npy_intp) 6};
	PyObject *NumPart_array = PyArray_ZEROS(1,Nkinds,NPY_INT32,0);
	PyObject *NumPart_array_file = PyArray_ZEROS(1,Nkinds,NPY_INT32,0);
	PyObject *npartHighWord_array = PyArray_ZEROS(1,Nkinds,NPY_UINT32,0);
	PyObject *Mass_array = PyArray_ZEROS(1,Nkinds,NPY_DOUBLE,0);

	if(NumPart_array==NULL || Mass_array==NULL || npartHighWord_array==NULL || NumPart_array_file==NULL){
		Py_XDECREF(NumPart_array);
		Py_XDECREF(Mass_array);
		Py_XDECREF(npartHighWord_array);
		Py_XDECREF(NumPart_array_file);
		return NULL;
	}

	//Get pointers to the array elements
	int *NumPart_data = (int *)PyArray_DATA(NumPart_array);
	int *NumPart_file_data = (int *)PyArray_DATA(NumPart_array_file);
	unsigned int *npartHighWord_data = (unsigned int *)PyArray_DATA(npartHighWord_array);
	double *Mass_data = (double *)PyArray_DATA(Mass_array);

	//Fill in the values
	Ngas = header.npartTotal[0];
	Ngas_file = header.npart[0];

	for(k=0;k<6;k++){
			
		NumPart_file_data[k] = header.npart[k];
		NumPartFile += header.npart[k];
		NumPart_data[k] = header.npartTotal[k];
		npartHighWord_data[k] = header.npartTotalHighWord[k];
		NumPart += header.npartTotal[k];
		Mass_data[k] = header.mass[k];
			
		if(header.mass[k]==0){ 
			Nwithmass+=header.npartTotal[k];
			Nwithmass_file+=header.npart[k];
		} 
		
	}

	//Build a dictionary with the header information
	PyObject *header_dict = PyDict_New();
	if(header_dict==NULL){
		return NULL;
	}

	//Fill in dictionary values
	if(PyDict_SetItemString(header_dict,"endianness",Py_BuildValue("i",endianness))) return NULL;
	if(PyDict_SetItemString(header_dict,"scale_factor",Py_BuildValue("d",header.time))) return NULL;
	if(PyDict_SetItemString(header_dict,"redshift",Py_BuildValue("d",header.redshift))) return NULL;
	if(PyDict_SetItemString(header_dict,"Om0",Py_BuildValue("d",header.Omega0))) return NULL;
	if(PyDict_SetItemString(header_dict,"Ode0",Py_BuildValue("d",header.OmegaLambda))) return NULL;
	if(PyDict_SetItemString(header_dict,"w0",Py_BuildValue("d",header.w0))) return NULL;
	if(PyDict_SetItemString(header_dict,"wa",Py_BuildValue("d",header.wa))) return NULL;
	if(PyDict_SetItemString(header_dict,"comoving_distance",Py_BuildValue("d",header.comoving_distance))) return NULL;
	if(PyDict_SetItemString(header_dict,"h",Py_BuildValue("d",header.HubbleParam))) return NULL;
	if(PyDict_SetItemString(header_dict,"box_size",Py_BuildValue("d",header.BoxSize))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_files",Py_BuildValue("i",header.num_files))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_total",Py_BuildValue("i",NumPart))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_file",Py_BuildValue("i",NumPartFile))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_total_gas",Py_BuildValue("i",Ngas))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_file_gas",Py_BuildValue("i",Ngas_file))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_total_with_mass",Py_BuildValue("i",Nwithmass))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_file_with_mass",Py_BuildValue("i",Nwithmass_file))) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_total_of_type",NumPart_array)) return NULL;
	if(PyDict_SetItemString(header_dict,"num_particles_file_of_type",NumPart_array_file)) return NULL;
	if(PyDict_SetItemString(header_dict,"npartTotalHighWord",npartHighWord_array)) return NULL;
	if(PyDict_SetItemString(header_dict,"masses",Mass_array)) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_cooling",Py_BuildValue("i",header.flag_cooling))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_feedback",Py_BuildValue("i",header.flag_feedback))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_sfr",Py_BuildValue("i",header.flag_sfr))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_stellarage",Py_BuildValue("i",header.flag_stellarage))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_metals",Py_BuildValue("i",header.flag_metals))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_entropy_instead_u",Py_BuildValue("i",header.flag_entropy_instead_u))) return NULL;

	//return
	return header_dict;

}


//getPosVel() implementation
static PyObject *_gadget_getPosVel(PyObject *self,PyObject *args){

	PyObject *file_obj;
	long offset;
	int NumPart;

	//Interpret the tuple of arguments
	if(!PyArg_ParseTuple(args,"Oli",&file_obj,&offset,&NumPart)){
		return NULL;
	}

	//Build the numpy array that will hold the particles positions, or velocities
	npy_intp dims[] = { (npy_intp) NumPart, (npy_intp) 3 };
	PyObject *particle_data_array = PyArray_ZEROS(2,dims,NPY_FLOAT32,0);

	if(particle_data_array==NULL){
		return NULL;
	}

	//Get a data pointer out of the array
	float *particle_data = (float *)PyArray_DATA(particle_data_array);

	//Get a file pointer out of the file object
	FILE *fp = PyFile_AsFile(file_obj); 
	PyFile_IncUseCount((PyFileObject *)file_obj);

	//Read in the positions of the partcles
	if(getPosVel(fp,offset,particle_data,NumPart)==-1){

		PyFile_DecUseCount((PyFileObject *)file_obj);
		Py_DECREF(particle_data_array);
		PyErr_SetString(PyExc_IOError,"End of file reached, the information requested is not available!");
		return NULL;
	
	}

	//Release the file pointer
	PyFile_DecUseCount((PyFileObject *)file_obj);

	//Return the array with the particle data (positions or velocities)
	return particle_data_array;

}


//getID() implementation
static PyObject *_gadget_getID(PyObject *self,PyObject *args){

	PyObject *file_obj;
	long offset;
	int NumPart;

	//Interpret the tuple of arguments
	if(!PyArg_ParseTuple(args,"Oli",&file_obj,&offset,&NumPart)){
		return NULL;
	}

	//Build the numpy array that will hold the particles IDs
	npy_intp dims[] = { (npy_intp) NumPart };
	PyObject *id_data_array = PyArray_ZEROS(1,dims,NPY_INT32,0);

	if(id_data_array==NULL){
		return NULL;
	}

	//Get a data pointer out of the array
	int *id_data = (int *)PyArray_DATA(id_data_array);

	//Get a file pointer out of the file object
	FILE *fp = PyFile_AsFile(file_obj); 
	PyFile_IncUseCount((PyFileObject *)file_obj);

	//Read in the IDs of the particles
	if(getID(fp,offset,id_data,NumPart)==-1){


		PyFile_DecUseCount((PyFileObject *)file_obj);
		Py_DECREF(id_data_array);
		PyErr_SetString(PyExc_IOError,"End of file reached, the information requested is not available!");
		return NULL;

	}

	//Release the file pointer
	PyFile_DecUseCount((PyFileObject *)file_obj);

	//Return the array with the particle data (positions or velocities)
	return id_data_array;

}

//write() implementation
static PyObject *_gadget_write(PyObject *self,PyObject *args){

	PyObject *header_obj,*positions_obj,*velocities_obj;
	const char *filename;
	int k,NumPart,firstID,writeVel;

	struct io_header_1 header;

	//interpret input tuple
	if(!PyArg_ParseTuple(args,"OOOisi",&header_obj,&positions_obj,&velocities_obj,&firstID,&filename,&writeVel)){
		return NULL;
	}

	//open the file to which the snapshot will be written, quit if cannot open the file
	FILE *fp = fopen(filename,"wb");
	if(fp==NULL){
		PyErr_SetString(PyExc_IOError,"Couldn't open snapshot file!");
		return NULL;
	}

	//interpret arrays
	PyObject *mass_array = PyArray_FROM_OTF(PyDict_GetItemString(header_obj,"masses"),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *NumPart_array = PyArray_FROM_OTF(PyDict_GetItemString(header_obj,"num_particles_total_of_type"),NPY_INT32,NPY_IN_ARRAY);
	PyObject *NumPart_file_array = PyArray_FROM_OTF(PyDict_GetItemString(header_obj,"num_particles_file_of_type"),NPY_INT32,NPY_IN_ARRAY);
	PyObject *npartHighWord_array = PyArray_FROM_OTF(PyDict_GetItemString(header_obj,"npartTotalHighWord"),NPY_UINT32,NPY_IN_ARRAY);

	if(mass_array==NULL || NumPart_array==NULL || NumPart_file_array==NULL || npartHighWord_array==NULL){
		Py_XDECREF(mass_array);
		Py_XDECREF(NumPart_array);
		Py_XDECREF(NumPart_file_array);
		Py_XDECREF(npartHighWord_array);
		fclose(fp);
		return NULL;
	}

	//get pointers
	double *mass_data = (double *)PyArray_DATA(mass_array);
	int *NumPart_data = (int *)PyArray_DATA(NumPart_array);
	int *NumPart_file_data = (int *)PyArray_DATA(NumPart_file_array);
	unsigned int *npartHighWord_data = (unsigned int *)PyArray_DATA(npartHighWord_array);

	//Fill in the header values

	//simple doubles
	header.Omega0 = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"Om0"));
	header.OmegaLambda = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"Ode0"));
	header.w0 = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"w0"));
	header.wa = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"wa"));
	header.time = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"scale_factor"));
	header.redshift = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"redshift"));
	header.BoxSize = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"box_size"));
	header.HubbleParam = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"h"));
	header.comoving_distance = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"comoving_distance"));

	//simple ints
	header.flag_cooling = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_cooling"));
	header.flag_sfr = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_sfr"));
	header.flag_feedback = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_feedback"));
	header.num_files = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"num_files"));
	header.flag_stellarage = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_stellarage"));
	header.flag_metals = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_metals"));
	header.flag_entropy_instead_u = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_entropy_instead_u"));


	//double and int arrays
	for(k=0;k<6;k++){
		header.mass[k] = mass_data[k];
		header.npart[k] = NumPart_file_data[k];
		header.npartTotal[k] = NumPart_data[k];
		header.npartTotalHighWord[k] = npartHighWord_data[k];
	}

	//release resources
	Py_DECREF(mass_array);
	Py_DECREF(NumPart_array);
	Py_DECREF(NumPart_file_array);
	Py_DECREF(npartHighWord_array);

	//now interpret positions and velocities as numpy arrays
	PyObject *positions_array = PyArray_FROM_OTF(positions_obj,NPY_FLOAT32,NPY_IN_ARRAY);
	PyObject *velocities_array = PyArray_FROM_OTF(velocities_obj,NPY_FLOAT32,NPY_IN_ARRAY);

	if(positions_array==NULL || velocities_array==NULL){
		Py_XDECREF(positions_array);
		Py_XDECREF(velocities_array);
		fclose(fp);
		return NULL;
	}

	//get the number of particles
	NumPart = (int)PyArray_DIM(positions_array,0);
	//get data pointers
	float *positions_data = (float *)PyArray_DATA(positions_array);
	float *velocities_data = (float *)PyArray_DATA(velocities_array);

	//ready to write Gadget snapshot, do it!
	if(writeSnapshot(fp,&header,positions_data,velocities_data,firstID,NumPart,writeVel)==-1){
		
		fclose(fp);
		PyErr_SetString(PyExc_IOError,"Couldn't write snapshot!");
		return NULL;

	}

	//release resources and close the snapshot file
	Py_DECREF(positions_array);
	Py_DECREF(velocities_array);
	fclose(fp);	

	Py_RETURN_NONE;

}

//grid3d() implementation
static PyObject *_gadget_grid3d(PyObject *self,PyObject *args){

	PyObject *positions_obj,*bins_obj;

	//parse input tuple
	if(!PyArg_ParseTuple(args,"OO",&positions_obj,&bins_obj)){
		return NULL;
	}

	//interpret parsed objects as arrays
	PyObject *positions_array = PyArray_FROM_OTF(positions_obj,NPY_FLOAT32,NPY_IN_ARRAY);
	PyObject *binsX_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,0),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binsY_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,1),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binsZ_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,2),NPY_DOUBLE,NPY_IN_ARRAY);

	//check if anything failed
	if(positions_array==NULL || binsX_array==NULL || binsY_array==NULL || binsZ_array==NULL){
		
		Py_XDECREF(positions_array);
		Py_XDECREF(binsX_array);
		Py_XDECREF(binsY_array);
		Py_XDECREF(binsZ_array);

		return NULL;
	}

	//Get data pointers
	float *positions_data = (float *)PyArray_DATA(positions_array);
	double *binsX_data = (double *)PyArray_DATA(binsX_array);
	double *binsY_data = (double *)PyArray_DATA(binsY_array);
	double *binsZ_data = (double *)PyArray_DATA(binsZ_array);

	//Get info about the number of bins
	int NumPart = (int)PyArray_DIM(positions_array,0);
	int nx = (int)PyArray_DIM(binsX_array,0) - 1;
	int ny = (int)PyArray_DIM(binsY_array,0) - 1;
	int nz = (int)PyArray_DIM(binsZ_array,0) - 1;

	//Allocate the new array for the grid
	npy_intp gridDims[] = {(npy_intp) nx,(npy_intp) ny,(npy_intp) nz};
	PyObject *grid_array = PyArray_ZEROS(3,gridDims,NPY_FLOAT32,0);

	if(grid_array==NULL){

		Py_DECREF(positions_array);
		Py_DECREF(binsX_array);
		Py_DECREF(binsY_array);
		Py_DECREF(binsZ_array);

		return NULL;

	}

	//Get a data pointer
	float *grid_data = (float *)PyArray_DATA(grid_array);

	//Snap the particles on the grid
	grid3d(positions_data,NumPart,binsX_data[0],binsY_data[0],binsZ_data[0],binsX_data[1] - binsX_data[0],binsY_data[1] - binsY_data[0],binsZ_data[1] - binsZ_data[0],nx,ny,nz,grid_data);

	//return the grid
	Py_DECREF(positions_array);
	Py_DECREF(binsX_array);
	Py_DECREF(binsY_array);
	Py_DECREF(binsZ_array);

	return grid_array;

}

//adaptive() implementation
static PyObject * _gadget_adaptive(PyObject *self,PyObject *args){

	PyObject *positions_obj,*rp_obj,*binning_obj;
	double center;
	int direction0,direction1,normal;

	//Parse argument tuple
	if(!PyArg_ParseTuple(args,"OOOdiii",&positions_obj,&rp_obj,&binning_obj,&center,&direction0,&direction1,&normal)){
		return NULL;
	}

	//Parse arrays
	PyObject *positions_array = PyArray_FROM_OTF(positions_obj,NPY_FLOAT32,NPY_IN_ARRAY);
	PyObject *rp_array = PyArray_FROM_OTF(rp_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binning0_array = PyArray_FROM_OTF(PyList_GetItem(binning_obj,0),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binning1_array = PyArray_FROM_OTF(PyList_GetItem(binning_obj,1),NPY_DOUBLE,NPY_IN_ARRAY);

	//Check if anything went wrong
	if(positions_array==NULL || rp_array==NULL || binning0_array==NULL || binning1_array==NULL){
		
		Py_XDECREF(positions_array);
		Py_XDECREF(rp_array);
		Py_XDECREF(binning0_array);
		Py_XDECREF(binning1_array);

		return NULL;
	}

	//Compute the number of particles
	int NumPart = (int)PyArray_DIM(positions_array,0);

	//Allocate space for lensing plane
	npy_intp dims[] =  {PyArray_DIM(binning0_array,0)-1,PyArray_DIM(binning1_array,0)-1};
	int size0 = (int)dims[0];
	int size1 = (int)dims[1];
	
	PyObject *lensingPlane_array = PyArray_ZEROS(2,dims,NPY_DOUBLE,0);
	
	if(lensingPlane_array==NULL){

		Py_DECREF(positions_array);
		Py_DECREF(rp_array);
		Py_DECREF(binning0_array);
		Py_DECREF(binning1_array);

		return NULL;

	}

	//Get data pointers
	float *positions = (float *)PyArray_DATA(positions_array);
	double *rp = (double *)PyArray_DATA(rp_array);
	double *binning0 = (double *)PyArray_DATA(binning0_array);
	double *binning1 = (double *)PyArray_DATA(binning1_array);
	double *lensingPlane = (double *)PyArray_DATA(lensingPlane_array);

	//Compute the adaptive smoothing using C backend
	adaptiveSmoothing(NumPart,positions,rp,binning0,binning1,center,direction0,direction1,normal,size0,size1,lensingPlane);

	//Cleanup
	Py_DECREF(positions_array);
	Py_DECREF(rp_array);
	Py_DECREF(binning0_array);
	Py_DECREF(binning1_array);

	//Return
	return lensingPlane_array;


}