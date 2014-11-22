/* ============================================================ *
 * config.h							*
 * Martin Kilbinger 2007					*
 * ============================================================ */

#ifndef __CONFIG_H
#define __CONFIG_H

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "errorlist.h"
#include "io.h"


#define conf_base  -400
#define conf_undef -1 + conf_base
#define conf_io    -2 + conf_base
#define conf_eof   -3 + conf_base


/* Maximal length of a short string config entry */
#define CSLENS 32


typedef enum {c_i, c_d, c_s} config_element_t;
#define sconfig_element_t(i) ( \
  i==c_i ? "integer" : \
  i==c_d ? "double"  : \
  i==c_s ? "string"  : \
  "" \
)

typedef struct {
  int i;
  double d;
  char s[512];
} config_element;


/* ============================================================ *
 * The following macros read entries from a configuration file. *
 * In the calling function the variables c (of type             *
 * config_element) and s (char[]) have to be defined.           *
 * ============================================================ */


/* Reads a config entry of the form "key    value" */
#define CONFIG_READ(config, element, type, F, c, err) \
 (config)->element = read_element(F, #element, c, c_##type, 1, err).type; \
 forwardError(*err, __LINE__,)



#define CONFIG_READ_S(config, element, s, F, c, err) \
  c = read_element(F, #element, c, c_s, 1, err); \
  strcpy((config)->element, c.s); \
  forwardError(*err, __LINE__,)


/* Reads an n-dimensional array of config entries of the form "key   value_1 value_2 value_3 ... value_n" */
#define CONFIG_READ_ARR(config, element, type, i, n, stmp, F, c, err)	\
  STRING_READ(#element, stmp, F, c, err);			\
  for (i=0; i<n; i++) {							\
     (config)->element[i] = read_element(F, "", c, c_##type, i==n-1, err).type; \
     forwardError(*err, __LINE__,);					\
  }


/* Reads an array of config entries of type 'string' of the form "key   string_1 string_2 string_3 ... string_n" */
#define CONFIG_READ_S_ARR(config, element, stmp, i, n, F, c, err)	\
  STRING_READ(#element, stmp, F, c, err);				\
  for (i=0; i<n; i++) {							\
  c = read_element(F, "", c, c_s, i==n-1, err);                         \
     strcpy((config)->element[i], c.s);                                 \
     forwardError(*err, __LINE__,);					\
  }


/* Assigns a value to element according to the string selement */
#define STRING2ENUM(element, selement, type, stype, j, N, err)	\
  element = -1;							\
    for (j=0; j<N; j++) { \
       if (strcmp(selement, stype((type)j))==0) element = (type)j;	\
    } \
    testErrorRetVA((int)element==-1, conf_undef, "Parametrization '%s' not found in type %s, cannot assign value of type %s", \
		   *err, __LINE__,, selement, #stype, #type)


config_element read_element(FILE *F, char *key, config_element c,
			    config_element_t t, int all, error **err);


/* ============================================================ *
 * The following macros are called internally by the above      *
 * defined macros.						*
 * ============================================================ */

/* If key="", only one string is read */
#define ENTRY_READ(F, key, c, type, err)				\
  (c = read_element(F, key, c, c_##type, 0, err), c.type)

#define STRING_READ(str, stmp, F, c, err)	\
  strcpy(stmp, ENTRY_READ(F, "", c, s, err));	\
  forwardError(*err, __LINE__,);			\
  testErrorRetVA(strcmp(stmp, str)!=0, conf_undef,	\
		 "String \"%s\" expected but \"%s\" found instead in config file", *err, __LINE__,, str, stmp);



#endif

