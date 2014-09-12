/*Python wrapper module for gadget snapshot reading;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _gadget and it defines the methods below (see docstrings)
*/

#include <stdio.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "read_gadget.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for reading Gadget2 snapshots";
static char getHeader_docstring[] = "Reads the header of a Gadget2 snapshot";
static char getPosVel_docstring[] = "Gets the positions or velocities of the particles in a Gadget2 snapshot";
static char getID_docstring[] = "Gets the 4 byte int particles IDs from the Gadget2 snapshot";
static char write_docstring[] = "Writes the particles information to a Gadget snapshot, with a proper header";

//Method declarations
static PyObject *_gadget_getHeader(PyObject *self,PyObject *args);
static PyObject *_gadget_getPosVel(PyObject *self,PyObject *args);
static PyObject *_gadget_getID(PyObject *self,PyObject *args);
static PyObject * _gadget_write(PyObject *self,PyObject *args);

//_gadget method definitions
static PyMethodDef module_methods[] = {

	{"getHeader",_gadget_getHeader,METH_VARARGS,getHeader_docstring},
	{"getPosVel",_gadget_getPosVel,METH_VARARGS,getPosVel_docstring},
	{"getID",_gadget_getID,METH_VARARGS,getID_docstring},
	{"write",_gadget_write,METH_VARARGS,write_docstring},
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

	//Exit if there was a problem (endianness check shall return 1 or 0, no exceptions; if not, the snapshot is messed up)
	if(endianness==-1){
		PyErr_SetString(PyExc_ValueError,"Couldn't determine snapshot endianness!");
		return NULL;
	}

	//Count the number of particles of each kind and total

	//Allocate resources
	npy_intp Nkinds[] = {(npy_intp) 6};
	PyObject *NumPart_array = PyArray_ZEROS(1,Nkinds,NPY_INT32,0);
	PyObject *NumPart_array_file = PyArray_ZEROS(1,Nkinds,NPY_INT32,0);
	PyObject *Mass_array = PyArray_ZEROS(1,Nkinds,NPY_DOUBLE,0);

	if(NumPart_array==NULL || Mass_array==NULL || NumPart_array_file==NULL){
		Py_XDECREF(NumPart_array);
		Py_XDECREF(Mass_array);
		Py_XDECREF(NumPart_array_file);
		return NULL;
	}

	//Get pointers to the array elements
	int *NumPart_data = (int *)PyArray_DATA(NumPart_array);
	int *NumPart_file_data = (int *)PyArray_DATA(NumPart_array_file);
	double *Mass_data = (double *)PyArray_DATA(Mass_array);

	//Fill in the values
	Ngas = header.npartTotal[0];
	Ngas_file = header.npart[0];

	for(k=0;k<6;k++){
			
		NumPart_file_data[k] = header.npart[k];
		NumPartFile += header.npart[k];
		NumPart_data[k] = header.npartTotal[k];
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
	if(PyDict_SetItemString(header_dict,"masses",Mass_array)) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_cooling",Py_BuildValue("i",header.flag_cooling))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_feedback",Py_BuildValue("i",header.flag_feedback))) return NULL;
	if(PyDict_SetItemString(header_dict,"flag_sfr",Py_BuildValue("i",header.flag_sfr))) return NULL;

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
		return NULL;

	}

	//Release the file pointer
	PyFile_DecUseCount((PyFileObject *)file_obj);

	//Return the array with the particle data (positions or velocities)
	return id_data_array;

}

//write() implementation
static PyObject *_gadget_write(PyObject *self,PyObject *args){

	PyObject *header_obj;
	const char *filename;
	int k;

	struct io_header_1 header;

	//interpret input tuple
	if(!PyArg_ParseTuple(args,"Os",&header_obj,&filename)){
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

	if(mass_array==NULL || NumPart_array==NULL || NumPart_file_array==NULL){
		Py_XDECREF(mass_array);
		Py_XDECREF(NumPart_array);
		Py_XDECREF(NumPart_file_array);
		fclose(fp);
		return NULL;
	}

	//get pointers
	double *mass_data = (double *)PyArray_DATA(mass_array);
	int *NumPart_data = (int *)PyArray_DATA(NumPart_array);
	int *NumPart_file_data = (int *)PyArray_DATA(NumPart_file_array);

	//Fill in the header values

	//simple doubles
	header.Omega0 = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"Om0"));
	header.OmegaLambda = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"Ode0"));
	header.time = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"scale_factor"));
	header.redshift = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"redshift"));
	header.BoxSize = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"box_size"));
	header.HubbleParam = PyFloat_AsDouble(PyDict_GetItemString(header_obj,"h"));

	//simple ints
	header.flag_cooling = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_cooling"));
	header.flag_sfr = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_sfr"));
	header.flag_feedback = (int)PyInt_AsLong(PyDict_GetItemString(header_obj,"flag_feedback"));
	header.num_files = 1;

	//double and int arrays
	for(k=0;k<6;k++){
		header.mass[k] = mass_data[k];
		header.npart[k] = NumPart_file_data[k];
		header.npartTotal[k] = NumPart_data[k];
	}

	//release resources
	Py_DECREF(mass_array);
	Py_DECREF(NumPart_array);
	Py_DECREF(NumPart_file_array);

	//ready to write Gadget snapshot, do it!
	if(writeSnapshot(fp,&header,NULL,NULL,1,1)==-1){
		
		fclose(fp);
		PyErr_SetString(PyExc_IOError,"Coulnd't write snapshot!");
		return NULL;

	}

	//close the snapshot file
	fclose(fp);	

	Py_RETURN_NONE;

}