/*************************************************************/
/*                                                           */
/*  Copyright (c) 1990                                       */
/*  Darrell L. Whitley                                       */
/*  Computer Science Department                              */
/*  Colorado State University                                */
/*                                                           */
/*  Permission is hereby granted to copy all or any part of  */
/*  this program for free distribution.   The author's name  */
/*  and this copyright notice must be included in any copy.  */
/*                                                           */
/*************************************************************/
#ifndef _GA_PARAMS_H_
#define _GA_PARAMS_H_

#ifdef ANSI
void 
parse_command_line (int argc, 
		    char *argv[]);

int  
parse_config_file (FILE *fp);

void 
print_params (FILE *fp);

int  
set_parameter (char tag, 
	       char value[]);

void usage (void);

#else

void parse_command_line ();
int  parse_config_file ();
void print_params ();
int  set_parameter ();
void usage ();

#endif



#endif
