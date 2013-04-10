/* File: Source_functionsmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Wed Mar  6 10:46:58 2013
 * $Revision:$
 * $Date:$
 * Do not edit this file directly unless you know what you are doing!!!
 */
#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
/*need_includes0*/

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *Source_functions_error;
static PyObject *Source_functions_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (((PyArrayObject *)(capi_ ## var ## _tmp))->nd)
#define old_shape(var,dim) (((PyArrayObject *)(capi_ ## var ## _tmp))->dimensions[dim])
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyFloat_Check(obj)) {
#ifdef __sgi
    *v = PyFloat_AsDouble(obj);
#else
    *v = PyFloat_AS_DOUBLE(obj);
#endif
    return 1;
  }
  tmp = PyNumber_Float(obj);
  if (tmp) {
#ifdef __sgi
    *v = PyFloat_AsDouble(tmp);
#else
    *v = PyFloat_AS_DOUBLE(tmp);
#endif
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = Source_functions_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyInt_Check(obj)) {
    *v = (int)PyInt_AS_LONG(obj);
    return 1;
  }
  tmp = PyNumber_Int(obj);
  if (tmp) {
    *v = PyInt_AS_LONG(tmp);
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = Source_functions_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** gravity **********************************/
static char doc_f2py_rout_Source_functions_source_functions_gravity[] = "\
Function signature:\n\
  gravity = gravity(in,t_in,g,n)\n\
Required arguments:\n"
"  in : input rank-1 array('d') with bounds (n)\n"
"  t_in : input float\n"
"  g : input rank-1 array('d') with bounds (3)\n"
"  n : input int\n"
"Return objects:\n"
"  gravity : rank-1 array('d') with bounds (n)";
/* #declfortranroutine# */
static PyObject *f2py_rout_Source_functions_source_functions_gravity(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double*,double*,double*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double *gravity = NULL;
  npy_intp gravity_Dims[1] = {-1};
  const int gravity_Rank = 1;
  PyArrayObject *capi_gravity_tmp = NULL;
  int capi_gravity_intent = 0;
  double *in = NULL;
  npy_intp in_Dims[1] = {-1};
  const int in_Rank = 1;
  PyArrayObject *capi_in_tmp = NULL;
  int capi_in_intent = 0;
  PyObject *in_capi = Py_None;
  double t_in = 0;
  PyObject *t_in_capi = Py_None;
  double *g = NULL;
  npy_intp g_Dims[1] = {-1};
  const int g_Rank = 1;
  PyArrayObject *capi_g_tmp = NULL;
  int capi_g_intent = 0;
  PyObject *g_capi = Py_None;
  int n = 0;
  PyObject *n_capi = Py_None;
  static char *capi_kwlist[] = {"in","t_in","g","n",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOO:Source_functions.source_functions.gravity",\
    capi_kwlist,&in_capi,&t_in_capi,&g_capi,&n_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable t_in */
    f2py_success = double_from_pyobj(&t_in,t_in_capi,"Source_functions.source_functions.gravity() 2nd argument (t_in) can't be converted to double");
  if (f2py_success) {
  /* Processing variable g */
  g_Dims[0]=3;
  capi_g_intent |= F2PY_INTENT_IN;
  capi_g_tmp = array_from_pyobj(PyArray_DOUBLE,g_Dims,g_Rank,capi_g_intent,g_capi);
  if (capi_g_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(Source_functions_error,"failed in converting 3rd argument `g' of Source_functions.source_functions.gravity to C/Fortran array" );
  } else {
    g = (double *)(capi_g_tmp->data);

  /* Processing variable n */
    f2py_success = int_from_pyobj(&n,n_capi,"Source_functions.source_functions.gravity() 4th argument (n) can't be converted to int");
  if (f2py_success) {
  /* Processing variable in */
  in_Dims[0]=n;
  capi_in_intent |= F2PY_INTENT_IN;
  capi_in_tmp = array_from_pyobj(PyArray_DOUBLE,in_Dims,in_Rank,capi_in_intent,in_capi);
  if (capi_in_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(Source_functions_error,"failed in converting 1st argument `in' of Source_functions.source_functions.gravity to C/Fortran array" );
  } else {
    in = (double *)(capi_in_tmp->data);

  /* Processing variable gravity */
  gravity_Dims[0]=n;
  capi_gravity_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_gravity_tmp = array_from_pyobj(PyArray_DOUBLE,gravity_Dims,gravity_Rank,capi_gravity_intent,Py_None);
  if (capi_gravity_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(Source_functions_error,"failed in converting hidden `gravity' of Source_functions.source_functions.gravity to C/Fortran array" );
  } else {
    gravity = (double *)(capi_gravity_tmp->data);

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(gravity,in,&t_in,g,&n);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("N",capi_gravity_tmp);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_gravity_tmp == NULL) ... else of gravity*/
  /* End of cleaning variable gravity */
  if((PyObject *)capi_in_tmp!=in_capi) {
    Py_XDECREF(capi_in_tmp); }
  }  /*if (capi_in_tmp == NULL) ... else of in*/
  /* End of cleaning variable in */
  } /*if (f2py_success) of n*/
  /* End of cleaning variable n */
  if((PyObject *)capi_g_tmp!=g_capi) {
    Py_XDECREF(capi_g_tmp); }
  }  /*if (capi_g_tmp == NULL) ... else of g*/
  /* End of cleaning variable g */
  } /*if (f2py_success) of t_in*/
  /* End of cleaning variable t_in */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************* end of gravity *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/

static FortranDataDef f2py_source_functions_def[] = {
  {"gravity",-1,{{-1}},0,NULL,(void *)f2py_rout_Source_functions_source_functions_gravity,doc_f2py_rout_Source_functions_source_functions_gravity},
  {NULL}
};

static void f2py_setup_source_functions(char *gravity) {
  int i_f2py=0;
  f2py_source_functions_def[i_f2py++].data = gravity;
}
extern void F_FUNC_US(f2pyinitsource_functions,F2PYINITSOURCE_FUNCTIONS)(void (*)(char *));
static void f2py_init_source_functions(void) {
  F_FUNC_US(f2pyinitsource_functions,F2PYINITSOURCE_FUNCTIONS)(f2py_setup_source_functions);
}

/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "Source_functions",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyObject *PyInit_Source_functions(void) {
#else
#define RETVAL
PyMODINIT_FUNC initSource_functions(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = Source_functions_module = PyModule_Create(&moduledef);
#else
  m = Source_functions_module = Py_InitModule("Source_functions", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module Source_functions (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'Source_functions' is auto-generated with f2py (version:2).\nFunctions:\n"
"Fortran 90/95 modules:\n""  source_functions --- gravity()"".");
  PyDict_SetItemString(d, "__doc__", s);
  Source_functions_error = PyErr_NewException ("Source_functions.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));

/*eof initf2pywraphooks*/
  PyDict_SetItemString(d, "source_functions", PyFortranObject_New(f2py_source_functions_def,f2py_init_source_functions));
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"Source_functions");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
