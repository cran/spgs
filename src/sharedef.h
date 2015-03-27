/*
This file is part of the source code for
SPGS: an R package for identifying statistical patterns in genomic sequences.
Copyright (C) 2015  Universidad de Chile and INRIA-Chile

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of Version 2 of the GNU Public License is available in the 
share/licenses/gpl-2 file in the R installation directory or from 
http://www.R-project.org/Licenses/GPL-2.
*/



/* Shared object defines */
/* Define DLL_EXPORT, DLL_IMPORT and DLL_LOCAL based on compiler. */

#if defined(_WIN32) || defined(__CYGWIN__)
  #define DLL_IMPORT __declspec(dllimport)
  #define DLL_EXPORT __declspec(dllexport)
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define DLL_IMPORT __attribute__ ((visibility ("default")))
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_IMPORT
    #define DLL_EXPORT
    #define DLL_LOCAL
  #endif
#endif

