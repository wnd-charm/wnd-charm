/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2007 Open Microscopy Environment                             */
/*         Massachusetts Institue of Technology,                                 */
/*         National Institutes of Health,                                        */
/*         University of Dundee                                                  */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:  Lior Shamir <shamirl [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef FuzzyCalcH
#define FuzzyCalcH
//---------------------------------------------------------------------------


#ifndef WIN32
typedef enum {clMin=-0x7fffffff-1, clMax=0x7fffffff} TColor;
#else
#include <vcl.h>
#endif

#define MIN_ORANGE_VAL 0.5

/* hue */
#define HUE_RED_START -0.1
#define HUE_RED_MAX 0
#define HUE_RED_END 11

#define HUE_DARK_ORANGE_START 0
#define HUE_DARK_ORANGE_MAX1 17   /* 15 */
#define HUE_DARK_ORANGE_MAX2 17   /* 18 */
#define HUE_DARK_ORANGE_END  27

#define HUE_LIGHT_ORANGE_START 20
#define HUE_LIGHT_ORANGE_MAX1 27  /* 26 */
#define HUE_LIGHT_ORANGE_MAX2 27
#define HUE_LIGHT_ORANGE_END 37

#define HUE_YELLOW_MAX 39
#define HUE_YELLOW_START 27
#define HUE_YELLOW_END 39.1

/* saturation */
#define SATURATION_GREY        1
#define SATURATION_ALMOST_GREY 2
#define SATURATION_TEND_GREY   3
#define SATURATION_MEDIUM_GREY 4
#define SATURATION_TEND_CLEAR  5
#define SATURATION_CLEAR       6

#define MIN_GREY 0 //20

#define GREY_START -0.1   // to avoid div by zero 
#define GREY_MAX 0
#define GREY_END 40

#define ALMOST_GREY_START 10
#define ALMOST_GREY_MAX   30
#define ALMOST_GREY_END   70

#define TEND_GREY_START 30
#define TEND_GREY_MAX 70
#define TEND_GREY_END 120

#define MEDIUM_GREY_START 60
#define MEDIUM_GREY_MAX 120
#define MEDIUM_GREY_END 180

#define TEND_CLEAR_START 120
#define TEND_CLEAR_MAX 180
#define TEND_CLEAR_END 240

#define CLEAR_START 180
#define CLEAR_END 240.1  // to avoid div by zero 
#define CLEAR_MAX 240

/* value */

#define VALUE_DARK        1
#define VALUE_ALMOST_DARK 2
#define VALUE_TEND_DARK   3
#define VALUE_TEND_LIGHT  4
#define VALUE_LIGHT       5

#define MIN_DARK 0 //40

#define DARK_START -0.1
#define DARK_MAX 0
#define DARK_END 60

#define ALMOST_DARK_START 0
#define ALMOST_DARK_MAX  60
#define ALMOST_DARK_END 120

#define TEND_DARK_START 60
#define TEND_DARK_MAX 120
#define TEND_DARK_END 180

#define TEND_LIGHT_START 120
#define TEND_LIGHT_MAX 180
#define TEND_LIGHT_END 240

#define LIGHT_START 200 //180
#define LIGHT_MAX 240
#define LIGHT_END 240.1

#define COLOR_WHITE        1
#define COLOR_LIGHT_GREY   2
#define COLOR_DARK_GREY    3
#define COLOR_BLACK        4

#define COLOR_RED          5
#define COLOR_PINK         6
#define COLOR_DARK_BROWN   7
#define COLOR_LIGHT_BROWN  8
#define COLOR_DARK_ORANGE  9
#define COLOR_LIGHT_ORANGE 10
#define COLOR_YELLOW       11
#define COLOR_OLIVE        12
#define COLOR_LIGHT_GREEN  13
#define COLOR_DARK_GREEN   14
#define COLOR_TEAL         15
#define COLOR_AQUA         16
#define COLOR_BLUE         17
#define COLOR_DARK_FUCIA   18
#define COLOR_LIGHT_FUCIA  19

#define COLORS_NUM 19


typedef struct COLOR_TYPE
{  char name[30];     /* the name of the color (e.g. "red") */
   TColor color;      /* the typycal color value (e.g. for the color red, color=0,0,255 */
}color_type;

typedef struct FRULE
{  int hue;
   int saturation;
   int value;
   int color;
} fuzzy_rule;

/* data structure that stores the 2 minimum points and 2 maximum point
   for a triangle-like membership function
   this function represent a trapez */
typedef struct trapez_FUNCTION
{  char name[30];
   int color;
   double start;
   double maximum1;
   double maximum2;
   double end;
}trapez_function;

/* membership functions */
double hue_red_value(double value);
double hue_dark_orange_value(double value);
double hue_light_orange_value(double value);
double hue_yellow_value(double value);

double grey_value(double value);
double almost_grey_value(double value);
double tend_grey_value(double value);
double medium_grey_value(double value);
double tend_clear_value(double value);
double clear_value(double value);

double dark_value(double value);
double almost_dark_value(double value);
double tend_dark_value(double value);
double tend_light_value(double value);
double light_value(double value);

double CalculateRules2(double hue,double saturation,double value);
long FindColor(short hue, short saturation, short value, double *certainties);
int color2num(char *color);
int saturation2num(char *saturation);
int value2num(char *value);
void num2color(int color,char *color_name);
TColor num2tcolor(int color);
void SetColors();
int LoadRules(char *filename);
int LoadColorsFunctions(char *filename);
double color_value(int color_index, double color);
//---------------------------------------------------------------------------
#endif



