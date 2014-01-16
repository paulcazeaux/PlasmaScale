/* HEADER XGrafixTools ====================================================================================
 * author: Paul Cazeaux
 * date: 2013-11-28
 * 
 * description:
 *  -> Link with the Xgrafix library
 * 
 * ============================================================================================== */

#ifndef DEF_PLASMASCALE_XGRAFIX_INTERFACE
#define DEF_PLASMASCALE_XGRAFIX_INTERFACE

/* includes ===================================================================================== */
#include "PlasmaScale.h"


extern "C"
void Quit();

extern "C"
void XGMainLoop();

extern "C"
void Dump(char *filename);

#endif