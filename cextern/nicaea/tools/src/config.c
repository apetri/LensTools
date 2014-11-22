#include "config.h"


#define this_testErrorRetVA(test, error_type, message, err, LINE, ...) { \
  if (test) {								\
    err = addErrorVA(error_type, message, err, LINE, __VA_ARGS__);	\
    free(line);							\
    c.i = -1;								\
    c.d = -1.0;								\
    strcpy(c.s, "");							\
    return c;								\
  } \
}



/* ============================================================ *
 * Reads a line (all=1) or string (all=0) from the (config)     *
 * file F."key value" of type t. Returns config_element with    *
 * value on success. all=0 is used if reading only the key or   *
 * elements of an array.					*
 * ============================================================ */
config_element read_element(FILE *F, char *key, config_element c,
			    config_element_t t, int all, error **err)
{
   int i, n, commentline;
   double d;
   char stmp[1024], formatkey[1024], *sptr;
   char *line=NULL;

   /* Empty config element */
   c.i = 0;
   c.d = 0.0;
   strcpy(c.s, "");

   strcpy(formatkey, key);

   line = malloc(1024*sizeof(char));

   /* Read from file while discarding comment lines (starting with '#') and empty lines */

   /* Get next string from file */
   if (all==0) {

      /* Read only one string from file */
      do {
	 n = fscanf(F, "%s", line);
	 this_testErrorRetVA(n != 1, conf_io, "%d elements read from file, expected %d", *err, __LINE__, n, 1);
	 if (line[0]=='#' || line[0]=='\n') {
	    sptr = fgets(line, 1024, F);
	    this_testErrorRetVA(sptr==NULL, conf_eof, "Premature end of file reached while scanning for key '%s'",
			      *err, __LINE__, key);
	    commentline = 1;
	 } else {
	    /* Check for whitespace-only line */
	    commentline = 1;
	    for (i=0; i<strlen(line); i++) {
	       //printf("Check %d '%c'\n", i, line[i]);
	       if (line[i]!=' ') {
		  commentline = 0;
		  break;
	       }
	    }
	 }
      } while (commentline==1);

   } else {

      /* Get next line and discard rest (e.g. comment) after value */
      do {
	 sptr = fgets(line, 1024, F);
	 //fprintf(stderr, "(%s) (%d) %d\n", line, sptr, 1);
	 this_testErrorRetVA(sptr==NULL, conf_eof, "Premature end of file reached while scanning for key '%s'",
			      *err, __LINE__, key);

	 /* Check for whitespace-only line */
	 commentline = 1;
	 for (i=0; i<strlen(line)-1; i++) { /* Last element is newline */
	    if (line[i]!=' ') {
	       commentline = 0;
	       break;
	    }
	 }

      } while (line[0]=='#' || line[0]=='\n' || commentline==1);

   }
   //fprintf(stderr, "all=%d,line=(%s)\n", all, line);

   /* Read value according to type (integer, double, or string) */
   switch (t) {
      case c_i :
	 strcat(formatkey, " %d ");
	 n = sscanf(line, formatkey, &i);
	 c.i = i;
	 break;

      case c_d :
	 strcat(formatkey, " %lf ");
	 n = sscanf(line, formatkey, &d);
	 c.d = d;
	 break;

      case c_s :
	 strcat(formatkey, " %s ");
	 n = sscanf(line, formatkey, stmp);
	 strcpy(c.s, stmp);
	 break;

      default : 
	 this_testErrorRetVA(0, conf_undef, "config_element type %d unknown", *err,
			     __LINE__, t);
   }


   /* Key not matched or end-of-file */
   if (n==0 || n==EOF) {
      if (n==0) {
	 chomp(line);
	 *err = addErrorVA(conf_io, "Error while reading file in line '%s': Expected key '%s' with %s value",
			   *err, __LINE__, line, key, sconfig_element_t(t));

      } else {
	 *err = addErrorVA(conf_io, "Unexpected eof in config file, line '%s'", *err, __LINE__, line);
      }

      free(line);
      c.i = -1; 
      c.d = -1.0; 
      strcpy(c.s, "");
      return c;
   }

   free(line);

   /* Success: return config element */
   return c;
}

