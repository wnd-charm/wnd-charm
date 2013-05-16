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


#include <stdio.h>
#include <math.h>
#include <string.h>

#include <stdlib.h>

#include "FuzzyCalc.h"

fuzzy_rule fuzzy_rules[1000];
trapez_function color_functions[50];

color_type colors[COLORS_NUM+1];

int rules_loaded=0;

char rulesfile[]="\ncolor_functions:\n\
\n\
red            0   0   0    11\n\
dark_orange    0   17  17  27\n\
light_orange   20  27  27  37\n\
yellow         27  39  39  47\n\
light_green    40  50  50  80\n\
dark_green     50  80  80  120\n\
aqua           80   120  120  160\n\
blue           120  160  160  200\n\
dark_fucia     160   200  200  220\n\
light_fucia    200   220  220  230\n\
red            220   240  240  240\n\
\n\
\n\
rules:\n\
\n\
\n\
\n\
//  ***  red ***\n\
\n    \
red         grey         dark             black\n\
red         grey         almost_dark      black\n\
red         grey         tend_dark        dark_grey\n\
red         grey         tend_light       light_grey\n\
red         grey         light            white\n\
\n\
red         almost_grey    dark             black\n\
red         almost_grey    almost_dark      dark_grey\n\
red         almost_grey    tend_dark        dark_grey\n\
red         almost_grey    tend_light       light_fucia\n\
red         almost_grey    light            pink\n\
\n\
red         tend_grey      dark             black\n\
red         tend_grey      almost_dark      dark_brown\n\
red         tend_grey      tend_dark        dark_fucia\n\
red         tend_grey      tend_light       pink\n\
red         tend_grey      light            pink\n\
\n\
red         medium_grey    dark             black\n\
red         medium_grey    almost_dark      dark_brown\n\
red         medium_grey    tend_dark        light_brown\n\
red         medium_grey    tend_light       light_fucia\n\
red         medium_grey    light            pink\n\
\n\
red         tend_clear   dark             black\n\
red         tend_clear   almost_dark      dark_brown\n\
red         tend_clear   tend_dark        dark_brown\n\
red         tend_clear   tend_light       red\n\
red         tend_clear   light            red\n\
\n\
red         clear        dark             black\n\
red         clear        almost_dark      black\n\
red         clear        tend_dark        red\n\
red         clear        tend_light       red\n\
red         clear        light            red\n\
\n\
//  ***  dark_orange ***    \n\
\n\
dark_orange   grey        dark            black\n\
dark_orange   grey        almost_dark     dark_grey\n\
dark_orange   grey        tend_dark       dark_grey\n\
dark_orange   grey        tend_light      light_grey\n\
dark_orange   grey        light           white\n\
\n\
dark_orange   almost_grey   dark            black\n\
dark_orange   almost_grey   almost_dark     dark_grey\n\
dark_orange   almost_grey   tend_dark       dark_grey\n\
dark_orange   almost_grey   tend_light      light_grey\n\
dark_orange   almost_grey   light           white\n\
\n\
dark_orange   tend_grey   dark            black\n\
dark_orange   tend_grey   almost_dark     dark_brown\n\
dark_orange   tend_grey   tend_dark       dark_grey\n\
dark_orange   tend_grey   tend_light      light_brown\n\
dark_orange   tend_grey   light           pink\n\
\n\
dark_orange   medium_grey   dark            black\n\
dark_orange   medium_grey   almost_dark     dark_brown\n\
dark_orange   medium_grey   tend_dark       dark_brown\n\
dark_orange   medium_grey   tend_light      light_brown\n\
dark_orange   medium_grey   light           light_orange\n\
\n\
dark_orange   tend_clear  dark            black\n\
dark_orange   tend_clear  almost_dark     dark_brown\n\
dark_orange   tend_clear  tend_dark       dark_brown\n\
dark_orange   tend_clear  tend_light      light_brown\n\
dark_orange   tend_clear  light           dark_orange\n\
\n\
dark_orange   clear       dark           black\n\
dark_orange   clear       almost_dark    dark_brown\n\
dark_orange   clear       tend_dark      dark_brown\n\
dark_orange   clear       tend_light     light_brown\n\
dark_orange   clear       light          dark_orange\n\
\n\
// *** light_orange ***\n\
\n\
light_orange   grey       dark           black\n\
light_orange   grey       almost_dark    dark_grey\n\
light_orange   grey       tend_dark      dark_grey\n\
light_orange   grey       tend_light     light_grey\n\
light_orange   grey       light          white\n\
\n\
light_orange   almost_grey    dark           black\n\
light_orange   almost_grey    almost_dark    dark_grey\n\
light_orange   almost_grey    tend_dark      dark_grey\n\
light_orange   almost_grey    tend_light     light_grey\n\
light_orange   almost_grey    light          white\n\
\n\
light_orange   tend_grey    dark           black\n\
light_orange   tend_grey    almost_dark    dark_grey\n\
light_orange   tend_grey    tend_dark      dark_grey\n\
light_orange   tend_grey    tend_light     light_grey\n\
light_orange   tend_grey    light          light_orange\n\
\n\
light_orange   medium_grey   dark         black\n\
light_orange   medium_grey   almost_dark  dark_brown\n\
light_orange   medium_grey   tend_dark    light_brown\n\
light_orange   medium_grey   tend_light   light_brown\n\
light_orange   medium_grey   light        light_orange\n\
\n\
light_orange   tend_clear  dark         black\n\
light_orange   tend_clear  almost_dark  dark_brown\n\
light_orange   tend_clear  tend_dark    light_brown\n\
light_orange   tend_clear  tend_light   light_brown\n\
light_orange   tend_clear  light        light_orange\n\
\n\
light_orange   clear       dark         black\n\
light_orange   clear       almost_dark  dark_brown\n\
light_orange   clear       tend_dark    light_brown\n\
light_orange   clear       tend_light   light_brown\n\
light_orange   clear       light        light_orange\n\
\n\
// *** yellow ***\n\
\n\
yellow         grey        dark         black\n\
yellow         grey        almost_dark  dark_grey\n\
yellow         grey        tend_dark    dark_grey\n\
yellow         grey        tend_light   light_grey\n\
yellow         grey        light        white\n\
\n\
yellow         almost_grey    dark         black\n\
yellow         almost_grey    almost_dark  dark_grey\n\
yellow         almost_grey    tend_dark    dark_grey\n\
yellow         almost_grey    tend_light   light_grey\n\
yellow         almost_grey    light        white\n\
\n\
yellow         tend_grey    dark         black\n\
yellow         tend_grey    almost_dark  dark_grey\n\
yellow         tend_grey    tend_dark    olive\n\
yellow         tend_grey    tend_light   olive\n\
yellow         tend_grey    light        yellow\n\
\n\
yellow         medium_grey  dark          black\n\
yellow         medium_grey  almost_dark   olive\n\
yellow         medium_grey  tend_dark     olive\n\
yellow         medium_grey  tend_light    olive\n\
yellow         medium_grey  light         yellow\n\
\n\
yellow         tend_clear  dark          black\n\
yellow         tend_clear  almost_dark   olive\n\
yellow         tend_clear  tend_dark     olive\n\
yellow         tend_clear  tend_light    olive\n\
yellow         tend_clear  light         yellow\n\
\n\
yellow         clear        dark          black\n\
yellow         clear        almost_dark   olive\n\
yellow         clear        tend_dark     olive\n\
yellow         clear        tend_light    olive\n\
yellow         clear        light         yellow\n\
\n\
// *** light_green ***\n\
\n\
light_green    grey        dark         black\n\
light_green    grey        almost_dark  dark_grey\n\
light_green    grey        tend_dark    dark_grey\n\
light_green    grey        tend_light   light_grey\n\
light_green    grey        light        white\n\
\n\
light_green    almost_grey     dark         black\n\
light_green    almost_grey     almost_dark  dark_grey\n\
light_green    almost_grey     tend_dark    dark_grey\n\
light_green    almost_grey     tend_light   light_grey\n\
light_green    almost_grey     light        white\n\
\n\
light_green    tend_grey     dark         black\n\
light_green    tend_grey     almost_dark  olive\n\
light_green    tend_grey     tend_dark    olive\n\
light_green    tend_grey     tend_light   light_green\n\
light_green    tend_grey     light        light_green\n\
\n\
light_green    medium_grey  dark          black\n\
light_green    medium_grey  almost_dark   dark_green\n\
light_green    medium_grey  tend_dark     dark_green\n\
light_green    medium_grey  tend_light    light_green\n\
light_green    medium_grey  light         light_green\n\
\n\
light_green    tend_clear  dark          black\n\
light_green    tend_clear  almost_dark   dark_green\n\
light_green    tend_clear  tend_dark     dark_green\n\
light_green    tend_clear  tend_light    light_green\n\
light_green    tend_clear  light         light_green\n\
\n\
light_green    clear        dark           black\n\
light_green    clear        almost_dark    dark_green\n\
light_green    clear        tend_dark      dark_green\n\
light_green    clear        tend_light     light_green\n\
light_green    clear        light          light_green\n\
\n\
// *** dark_green ***\n\
\n\
dark_green    grey         dark          black\n\
dark_green    grey         almost_dark   dark_grey\n\
dark_green    grey         tend_dark     dark_grey\n\
dark_green    grey         tend_light    light_grey\n\
dark_green    grey         light         white\n\
\n\
dark_green    almost_grey         dark          black\n\
dark_green    almost_grey         almost_dark   dark_grey\n\
dark_green    almost_grey         tend_dark     dark_grey\n\
dark_green    almost_grey         tend_light    light_green\n\
dark_green    almost_grey         light         light_green\n\
\n\
dark_green    tend_grey         dark          black\n\
dark_green    tend_grey         almost_dark   dark_green\n\
dark_green    tend_grey         tend_dark     dark_green\n\
dark_green    tend_grey         tend_light    light_green\n\
dark_green    tend_grey         light         light_green\n\
\n\
dark_green    medium_grey   dark          black\n\
dark_green    medium_grey   almost_dark   dark_green\n\
dark_green    medium_grey   tend_dark     dark_green\n\
dark_green    medium_grey   tend_light    light_green\n\
dark_green    medium_grey   light         light_green\n\
\n\
dark_green    tend_clear   dark          black\n\
dark_green    tend_clear   almost_dark   dark_green\n\
dark_green    tend_clear   tend_dark     dark_green\n\
dark_green    tend_clear   tend_light    light_green\n\
dark_green    tend_clear   light         light_green\n\
\n\
dark_green    clear         dark          black\n\
dark_green    clear         almost_dark   dark_green\n\
dark_green    clear         tend_dark     dark_green\n\
dark_green    clear         tend_light    light_green\n\
dark_green    clear         light         light_green\n\
\n\
// *** aqua ***\n\
\n\
aqua           grey          dark          black\n\
aqua           grey          almost_dark   dark_grey\n\
aqua           grey          tend_dark     dark_grey\n\
aqua           grey          tend_light    light_grey\n\
aqua           grey          light         white\n\
\n\
aqua           almost_grey     dark          black\n\
aqua           almost_grey     almost_dark   dark_grey\n\
aqua           almost_grey     tend_dark     dark_grey\n\
aqua           almost_grey     tend_light    teal\n\
aqua           almost_grey     light         aqua\n\
\n\
aqua           tend_grey     dark          black\n\
aqua           tend_grey     almost_dark   teal\n\
aqua           tend_grey     tend_dark     blue\n\
aqua           tend_grey     tend_light    blue\n\
aqua           tend_grey     light         aqua\n\
\n\
aqua           medium_grey    dark          black\n\
aqua           medium_grey    almost_dark   teal\n\
aqua           medium_grey    tend_dark     blue\n\
aqua           medium_grey    tend_light    aqua\n\
aqua           medium_grey    light         aqua\n\
\n\
aqua           tend_clear    dark         black\n\
aqua           tend_clear    almost_dark  teal\n\
aqua           tend_clear    tend_dark    blue\n\
aqua           tend_clear    tend_light   aqua\n\
aqua           tend_clear    light        aqua\n\
\n\
aqua           clear          dark          black\n\
aqua           clear          almost_dark   teal\n\
aqua           clear          tend_dark     teal\n\
aqua           clear          tend_light    aqua\n\
aqua           clear          light         aqua\n\
\n\
// *** blue ***\n\
\n\
blue           grey           dark         black\n\
blue           grey           almost_dark  dark_grey\n\
blue           grey           tend_dark    dark_grey\n\
blue           grey           tend_light   light_grey\n\
blue           grey           light        white\n\
\n\
\n\
blue           almost_grey    dark         black\n\
blue           almost_grey    almost_dark  dark_grey\n\
blue           almost_grey    tend_dark    light_gray\n\
blue           almost_grey    tend_light   aqua\n\
blue           almost_grey    light        aqua\n\
\n\
blue           tend_grey    dark         black\n\
blue           tend_grey    almost_dark  dark_fucia\n\
blue           tend_grey    tend_dark    ligh_fucia\n\
blue           tend_grey    tend_light   light_fucia\n\
blue           tend_grey    light        aqua\n\
\n\
blue           medium_grey     dark         black\n\
blue           medium_grey     almost_dark  blue\n\
blue           medium_grey     tend_dark    blue\n\
blue           medium_grey     tend_light   blue\n\
blue           medium_grey     light        blue\n\
\n\
blue           tend_clear     dark         black\n\
blue           tend_clear     almost_dark  blue\n\
blue           tend_clear     tend_dark    blue\n\
blue           tend_clear     tend_light   blue\n\
blue           tend_clear     light        blue\n\
\n\
blue           clear          dark         black\n\
blue           clear          almost_dark  blue\n\
blue           clear          tend_dark    blue\n\
blue           clear          tend_light   blue\n\
blue           clear          light        blue\n\
\n\
// *** dark_fucia ***\n\
\n\
dark_fucia    grey          dark        black\n\
dark_fucia    grey          almost_dark dark_grey\n\
dark_fucia    grey          tend_dark   dark_grey\n\
dark_fucia    grey          tend_light  light_grey\n\
dark_fucia    grey          light       white\n\
\n\
dark_fucia    almost_grey   dark        black\n\
dark_fucia    almost_grey   almost_dark dark_grey\n\
dark_fucia    almost_grey   tend_dark   dark_fucia\n\
dark_fucia    almost_grey   tend_light  light_fucia\n\
dark_fucia    almost_grey   light       light_fucia\n\
\n\
dark_fucia    tend_grey   dark        black\n\
dark_fucia    tend_grey   almost_dark dark_fucia\n\
dark_fucia    tend_grey   tend_dark   dark_fucia\n\
dark_fucia    tend_grey   tend_light  light_fucia\n\
dark_fucia    tend_grey   light       light_fucia\n\
\n\
dark_fucia    medium_grey     dark        black\n\
dark_fucia    medium_grey     almost_dark dark_fucia\n\
dark_fucia    medium_grey     tend_dark   dark_fucia\n\
dark_fucia    medium_grey     tend_light  light_fucia\n\
dark_fucia    medium_grey     light       light_fucia\n\
\n\
dark_fucia    tend_clear    dark         black\n\
dark_fucia    tend_clear    almost_dark  dark_fucia\n\
dark_fucia    tend_clear    tend_dark    dark_fucia\n\
dark_fucia    tend_clear    tend_light   light_fucia\n\
dark_fucia    tend_clear    light        light_fucia\n\
\n\
dark_fucia    clear          dark         black\n\
dark_fucia    clear          almost_dark  dark_brown\n\
dark_fucia    clear          tend_dark    dark_fucia\n\
dark_fucia    clear          tend_light   dark_fucia\n\
dark_fucia    clear          light        light_fucia\n\
\n\
// *** light fucia ***\n\
\n\
light_fucia    grey          dark          black\n\
light_fucia    grey          almost_dark   dark_grey\n\
light_fucia    grey          tend_dark     dark_grey\n\
light_fucia    grey          tend_light    light_grey\n\
light_fucia    grey          light         white\n\
\n\
light_fucia    almost_grey   dark          black\n\
light_fucia    almost_grey   almost_dark   dark_grey\n\
light_fucia    almost_grey   tend_dark     light_grey\n\
light_fucia    almost_grey   tend_light    light_fucia\n\
light_fucia    almost_grey   light         light_fucia\n\
\n\
light_fucia    tend_grey   dark          black\n\
light_fucia    tend_grey   almost_dark   dark_fucia\n\
light_fucia    tend_grey   tend_dark     dark_fucia\n\
light_fucia    tend_grey   tend_light    light_fucia\n\
light_fucia    tend_grey   light         light_fucia\n\
\n\
light_fucia    medium_grey    dark          black\n\
light_fucia    medium_grey    almost_dark   dark_fucia\n\
light_fucia    medium_grey    tend_dark     dark_fucia\n\
light_fucia    medium_grey    tend_light    light_fucia\n\
light_fucia    medium_grey    light         light_fucia\n\
\n\
light_fucia    tend_clear   dark         black\n\
light_fucia    tend_clear   almost_dark  dark_brown\n\
light_fucia    tend_clear   tend_dark    dark_fucia\n\
light_fucia    tend_clear   tend_light   light_fucia\n\
light_fucia    tend_clear   light        light_fucia\n\
\n\
light_fucia    clear         dark           black\n\
light_fucia    clear         almost_dark    dark_brown\n\
light_fucia    clear         tend_dark      dark_fucia\n\
light_fucia    clear         tend_light     dark_fucia\n\
light_fucia    clear         light          light_fucia\n\
\n\
\n\
";


TColor num2tcolor(int color) {
	return(colors[color].color);
}

void num2color(int color,char *color_name) {
	if (color>COLORS_NUM) strcpy(color_name,"");
	else strcpy(color_name,colors[color].name);
}

int color2num(char *color) {
	int color_index;
	for (color_index=COLOR_WHITE;color_index<=COLORS_NUM;color_index++)
	if (colors[color_index].color>=0) {
		if (strcmp(colors[color_index].name,color)==0)
		return(color_index);
	}
	return(0); /* color was not found */
}


int saturation2num(char *saturation) {
	if (strcmp(saturation,"grey")==0) return(SATURATION_GREY);
	if (strcmp(saturation,"almost_grey")==0) return(SATURATION_ALMOST_GREY);
	if (strcmp(saturation,"tend_grey")==0) return(SATURATION_TEND_GREY);
	if (strcmp(saturation,"medium_grey")==0) return(SATURATION_MEDIUM_GREY);
	if (strcmp(saturation,"tend_clear")==0) return(SATURATION_TEND_CLEAR);
	if (strcmp(saturation,"clear")==0) return(SATURATION_CLEAR);
	return(-1);
}

int value2num(char *value) {
	if (strcmp(value,"dark")==0) return(VALUE_DARK);
	if (strcmp(value,"almost_dark")==0) return(VALUE_ALMOST_DARK);
	if (strcmp(value,"tend_dark")==0) return(VALUE_TEND_DARK);
	if (strcmp(value,"tend_light")==0) return(VALUE_TEND_LIGHT);
	if (strcmp(value,"light")==0) return(VALUE_LIGHT);
	return(-1);
}

/* load the colors to the colors structures */
void SetColors() {
	int color_index;

	/* initialize */
	for (color_index=COLOR_WHITE;color_index<=COLORS_NUM;color_index++)
		colors[color_index].color=(TColor)-1;

	/* set the colors */
	strcpy(colors[COLOR_WHITE].name,"white");
	colors[COLOR_WHITE].color=(TColor)0x00FFFFFF;
	strcpy(colors[COLOR_LIGHT_GREY].name,"light_grey");
	colors[COLOR_LIGHT_GREY].color=(TColor)0x00BDBEBD;
	strcpy(colors[COLOR_DARK_GREY].name,"dark_grey");
	colors[COLOR_DARK_GREY].color=(TColor)0x007B7D7B;
	strcpy(colors[COLOR_BLACK].name,"black");
	colors[COLOR_BLACK].color=(TColor)0x00000000;
	strcpy(colors[COLOR_RED].name,"red");
	colors[COLOR_RED].color=(TColor)0x0000FF;
	strcpy(colors[COLOR_PINK].name,"pink");
	colors[COLOR_PINK].color=(TColor)0x008080FF;
	strcpy(colors[COLOR_DARK_BROWN].name,"dark_brown");
	colors[COLOR_DARK_BROWN].color=(TColor)0x0018497B;
	strcpy(colors[COLOR_LIGHT_BROWN].name,"light_brown");
	colors[COLOR_LIGHT_BROWN].color=(TColor)0x005A8ABD;
	strcpy(colors[COLOR_DARK_ORANGE].name,"dark_orange");
	colors[COLOR_DARK_ORANGE].color=(TColor)0x006DFF;
	strcpy(colors[COLOR_LIGHT_ORANGE].name,"light_orange");
	colors[COLOR_LIGHT_ORANGE].color=(TColor)0x0000AEFF;
	strcpy(colors[COLOR_YELLOW].name,"yellow");
	colors[COLOR_YELLOW].color=(TColor)0x0000FFFF;
	strcpy(colors[COLOR_OLIVE].name,"olive");
	colors[COLOR_OLIVE].color=(TColor)0x00008080;
	strcpy(colors[COLOR_LIGHT_GREEN].name,"light_green");
	colors[COLOR_LIGHT_GREEN].color=(TColor)0x0000FF00;
	strcpy(colors[COLOR_DARK_GREEN].name,"dark_green");
	colors[COLOR_DARK_GREEN].color=(TColor)0x00008000;
	strcpy(colors[COLOR_TEAL].name,"teal");
	colors[COLOR_TEAL].color=(TColor)0x00808000;
	strcpy(colors[COLOR_AQUA].name,"aqua");
	colors[COLOR_AQUA].color=(TColor)0x00FFFF00;
	strcpy(colors[COLOR_BLUE].name,"blue");
	colors[COLOR_BLUE].color=(TColor)0x00FF0000;
	strcpy(colors[COLOR_DARK_FUCIA].name,"dark_fucia");
	colors[COLOR_DARK_FUCIA].color=(TColor)0x00800080;
	strcpy(colors[COLOR_LIGHT_FUCIA].name,"light_fucia");
	colors[COLOR_LIGHT_FUCIA].color=(TColor)0x00FF00FF;
}

char *getline(char *buffer) {
	char *p,*base;
	base=p=buffer;
	if (*p=='\0') {p++;base++;}
	while (*p!='\n' && *p!='\0') p++;
	if (*p=='\0') return(NULL);
	*p='\0';
	return(base);
}

int LoadRules(char *rulesfile) {
	char *p_line;
	int RulesCounter=0,rules=0;
	p_line=getline(rulesfile);

	/* read the file and create the rules */
	while (p_line) {
		if (strstr(p_line,"rules:")) {
			rules=1;
			p_line=getline(&(p_line[strlen(p_line)]));
			continue;
		}

		if (rules && strlen(p_line) > 4 && p_line[0]!='/') {
			p_line=strtok(p_line," \n\t");
			fuzzy_rules[RulesCounter].hue=color2num(p_line);
			p_line=strtok(NULL," \n\t");
			fuzzy_rules[RulesCounter].saturation=saturation2num(p_line);
			p_line=strtok(NULL," \n\t");
			fuzzy_rules[RulesCounter].value=value2num(p_line);
			p_line=strtok(NULL," \n\t");
			fuzzy_rules[RulesCounter].color=color2num(p_line);
			RulesCounter++;
		}
		p_line=getline(&(p_line[strlen(p_line)]));
	}
	/* add the last rule which is an empty rule */
	fuzzy_rules[RulesCounter].hue = -1;
	fuzzy_rules[RulesCounter].saturation = -1;
	fuzzy_rules[RulesCounter].value = -1;
	fuzzy_rules[RulesCounter].color = -1;
	return(1);
}

int LoadColorsFunctions(char *filename) {
	char *p_line;
	int FunctionCounter=0,colorfunctions=0;

	p_line=getline(rulesfile);
	/* read the file and create the rules */
	while (p_line) {
		if (strstr(p_line,"color_functions:")) {
			colorfunctions=1;
			p_line=getline(&(p_line[strlen(p_line)]));
			continue;
		}
		if (strstr(p_line,"rules:")) {
			*(strchr(p_line,'\0'))='\n';
			return(1);
		}
		if (colorfunctions && strlen(p_line) > 4) {
			p_line=strtok(p_line," \n\t");
			strcpy(color_functions[FunctionCounter].name,p_line);
			color_functions[FunctionCounter].color=color2num(p_line);
			p_line=strtok(NULL," \n\t");
			color_functions[FunctionCounter].start=atof(p_line);
			p_line=strtok(NULL," \n\t");
			color_functions[FunctionCounter].maximum1=atof(p_line);
			p_line=strtok(NULL," \n\t");
			color_functions[FunctionCounter].maximum2=atof(p_line);
			p_line=strtok(NULL," \n\t");
			color_functions[FunctionCounter].end=atof(p_line);
			FunctionCounter++;
		}
		p_line=getline(&(p_line[strlen(p_line)]));
	}
	/* add the last membership function which is an empty function */
	strcpy(color_functions[FunctionCounter].name,"");
	color_functions[FunctionCounter].color=-1;
	color_functions[FunctionCounter].start=0;
	color_functions[FunctionCounter].maximum1=0;
	color_functions[FunctionCounter].maximum2=0;
	color_functions[FunctionCounter].end=0;

	return(1);
}

double color_value(int color_index, double color) {
	int FunctionCounter=0;
	double res=0;
	while (color_functions[FunctionCounter].color>=0) {
		if (color_index==color_functions[FunctionCounter].color) {
			if (!res
				&& (color<color_functions[FunctionCounter].start
				|| color>color_functions[FunctionCounter].end)
			)
				res = 0;           /* out of the range */
			else if (!res
				&& color>=color_functions[FunctionCounter].maximum1
				&& color<=color_functions[FunctionCounter].maximum2
			)
				res = 1; /* the top of the trapez maximum */
			else if (!res && color<color_functions[FunctionCounter].maximum1)
				res = (color-color_functions[FunctionCounter].start)/(color_functions[FunctionCounter].maximum1-color_functions[FunctionCounter].start);
			else if (!res && color>color_functions[FunctionCounter].maximum2)
				res = 1-(color-color_functions[FunctionCounter].maximum2)/(color_functions[FunctionCounter].end-color_functions[FunctionCounter].maximum2);
			if (res>0) return(res);
		}
		FunctionCounter++;
	}
	return(res);
}

/***************************************/
/* membership functions for saturation */
/***************************************/
/* membership function for "grey" */
double grey_value(double value) {
	if (value<MIN_GREY) return(0);   /* minimum */
	if (value<GREY_START || value>GREY_END) return(0);   /* out of the range      */
	if (value<GREY_MAX) return((value-GREY_START)/(GREY_MAX-GREY_START));
	else return(1-((value-GREY_MAX)/(GREY_END-GREY_MAX)));
}
/* membership function for "almost_grey" */
double almost_grey_value(double value) {
	if (value<ALMOST_GREY_START || value>ALMOST_GREY_END) return(0);   /* out of the range      */
	if (value<ALMOST_GREY_MAX) return((value-ALMOST_GREY_START)/(ALMOST_GREY_MAX-ALMOST_GREY_START));
	else return(1-((value-ALMOST_GREY_MAX)/(ALMOST_GREY_END-ALMOST_GREY_MAX)));
}
/* membership function for "tend_grey" */
double tend_grey_value(double value) {
	if (value<TEND_GREY_START || value>TEND_GREY_END) return(0);   /* out of the range      */
	if (value<TEND_GREY_MAX) return((value-TEND_GREY_START)/(TEND_GREY_MAX-TEND_GREY_START));
	else return(1-((value-TEND_GREY_MAX)/(TEND_GREY_END-TEND_GREY_MAX)));
}
/* membership function for "medium_grey" */
double medium_grey_value(double value) {
	if (value<MEDIUM_GREY_START || value>MEDIUM_GREY_END) return(0);   /* out of the range      */
	if (value<MEDIUM_GREY_MAX) return((value-MEDIUM_GREY_START)/(MEDIUM_GREY_MAX-MEDIUM_GREY_START));
	else return(1-((value-MEDIUM_GREY_MAX)/(MEDIUM_GREY_END-MEDIUM_GREY_MAX)));
}
/* membership function for "tend_clear" */
double tend_clear_value(double value) {
	if (value<TEND_CLEAR_START || value>TEND_CLEAR_END) return(0);   /* out of the range      */
	if (value<TEND_CLEAR_MAX) return((value-TEND_CLEAR_START)/(TEND_CLEAR_MAX-TEND_CLEAR_START));
	else return(1-((value-TEND_CLEAR_MAX)/(TEND_CLEAR_END-TEND_CLEAR_MAX)));
}
/* membership function for "clear" */
double clear_value(double value) {
	if (value<CLEAR_START || value>CLEAR_END) return(0);   /* out of the range      */
	if (value<=CLEAR_MAX) return((value-CLEAR_START)/(CLEAR_MAX-CLEAR_START));
	else return(1-((value-CLEAR_MAX)/(CLEAR_END-CLEAR_MAX)));
}

/**********************************/
/* memnership functions for value */
/**********************************/
/* membership function for "dark" */
double dark_value(double value) {
	if (value<MIN_DARK) return(0);  /* below the minimum */
	if (value<DARK_START || value>DARK_END) return(0);   /* out of the range      */
	if (value<DARK_MAX) return((value-DARK_START)/(DARK_MAX-DARK_START));
	else return(1-((value-DARK_MAX)/(DARK_END-DARK_MAX)));
}
/* membership function for "almost_dark" */
double almost_dark_value(double value) {
	if (value<ALMOST_DARK_START || value>ALMOST_DARK_END) return(0);   /* out of the range      */
	if (value<ALMOST_DARK_MAX) return((value-ALMOST_DARK_START)/(ALMOST_DARK_MAX-ALMOST_DARK_START));
	else return(1-((value-ALMOST_DARK_MAX)/(ALMOST_DARK_END-ALMOST_DARK_MAX)));
}
/* membership function for "tend_dark" */
double tend_dark_value(double value) {
	if (value<TEND_DARK_START || value>TEND_DARK_END) return(0);   /* out of the range      */
	if (value<TEND_DARK_MAX) return((value-TEND_DARK_START)/(TEND_DARK_MAX-TEND_DARK_START));
	else return(1-((value-TEND_DARK_MAX)/(TEND_DARK_END-TEND_DARK_MAX)));
}
/* membership function for "tend_light" */
double tend_light_value(double value) {
	if (value<TEND_LIGHT_START || value>TEND_LIGHT_END) return(0);   /* out of the range      */
	if (value<TEND_LIGHT_MAX) return((value-TEND_LIGHT_START)/(TEND_LIGHT_MAX-TEND_LIGHT_START));
	else return(1-((value-TEND_LIGHT_MAX)/(TEND_LIGHT_END-TEND_LIGHT_MAX)));
}
/* membership function for "light" */
double light_value(double value) {
	if (value<LIGHT_START || value>LIGHT_END) return(0);   /* out of the range      */
	if (value<=LIGHT_MAX) return((value-LIGHT_START)/(LIGHT_MAX-LIGHT_START));
	else return(1-((value-LIGHT_MAX)/(LIGHT_END-LIGHT_MAX)));
}

/*********************************/
/* colors membership functions   */
/*********************************/
/* membership function for "red" */
double hue_red_value(double value) {
	if (value<HUE_RED_START || value>HUE_RED_END) return(0); /* out of the range */
	if (value<HUE_RED_MAX) return((value-HUE_RED_START)/(HUE_RED_MAX-HUE_RED_START));
	else return(1-((value-HUE_RED_MAX)/(HUE_RED_END-HUE_RED_MAX)));
}
/* membership function for "dark orange" */
double hue_dark_orange_value(double value) {
	if (value<HUE_DARK_ORANGE_START || value>HUE_DARK_ORANGE_END) return(0);   /* out of the range      */
	if (value>=HUE_DARK_ORANGE_MAX1 && value<=HUE_DARK_ORANGE_MAX2) return(1); /* the top of the trapez */
	if (value<HUE_DARK_ORANGE_MAX1) return((value-HUE_DARK_ORANGE_START)/(HUE_DARK_ORANGE_MAX1-HUE_DARK_ORANGE_START));
	else return(1-((value-HUE_DARK_ORANGE_MAX2)/(HUE_DARK_ORANGE_END-HUE_DARK_ORANGE_MAX2)));
}
/* membership function for "light orange" */
double hue_light_orange_value(double value) {
	if (value<HUE_LIGHT_ORANGE_START || value>HUE_LIGHT_ORANGE_END) return(0);   /* out of the range      */
	if (value>=HUE_LIGHT_ORANGE_MAX1 && value<=HUE_LIGHT_ORANGE_MAX2) return(1); /* the top of the trapez */
	if (value<HUE_LIGHT_ORANGE_MAX1) return((value-HUE_LIGHT_ORANGE_START)/(HUE_LIGHT_ORANGE_MAX1-HUE_LIGHT_ORANGE_START));
	else return(1-((value-HUE_LIGHT_ORANGE_MAX2)/(HUE_LIGHT_ORANGE_END-HUE_LIGHT_ORANGE_MAX2)));
}
/* membership function for "yellow" */
double hue_yellow_value(double value) {
	if (value<HUE_YELLOW_START || value>HUE_YELLOW_END) return(0);   /* out of the range      */
	if (value<=HUE_YELLOW_MAX) return((value-HUE_YELLOW_START)/(HUE_YELLOW_MAX-HUE_YELLOW_START));
	else return(1-((value-HUE_YELLOW_MAX)/(HUE_YELLOW_END-HUE_YELLOW_MAX)));
}
//---------------------------------------------------------------------------
double CalculateRules2(double hue, double saturation, double value, int color) {
	double ret_val,lower_sum=0;
	int RulesCounter=0;
	while (fuzzy_rules[RulesCounter].hue>=0) {
		/*
		if (do_grey_scale || (strcmp(fuzzy_rules[RulesCounter].color,"white")
		&& strcmp(fuzzy_rules[RulesCounter].color,"light_grey")
		&& strcmp(fuzzy_rules[RulesCounter].color,"dark_grey")
		&& strcmp(fuzzy_rules[RulesCounter].color,"black")))
		{
		*/
		
		// calculate only rules of the same color
		if (fuzzy_rules[RulesCounter].color == color) {
			/* hue functions */
			ret_val=color_value(fuzzy_rules[RulesCounter].hue,hue);
			if (ret_val==0) {
				RulesCounter++;
				continue;
			}
			/* saturation function */
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_GREY) ret_val=ret_val*grey_value(saturation);
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_ALMOST_GREY) ret_val=ret_val*almost_grey_value(saturation);
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_TEND_GREY) ret_val=ret_val*tend_grey_value(saturation);
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_MEDIUM_GREY) ret_val=ret_val*medium_grey_value(saturation);
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_TEND_CLEAR) ret_val=ret_val*tend_clear_value(saturation);
			if (fuzzy_rules[RulesCounter].saturation==SATURATION_CLEAR) ret_val=ret_val*clear_value(saturation);
			if (ret_val==0) {
				RulesCounter++;
				continue;
			}
			/* value functions */
			if (fuzzy_rules[RulesCounter].value==VALUE_DARK) ret_val=ret_val*dark_value(value);
			if (fuzzy_rules[RulesCounter].value==VALUE_ALMOST_DARK) ret_val=ret_val*almost_dark_value(value);
			if (fuzzy_rules[RulesCounter].value==VALUE_TEND_DARK) ret_val=ret_val*tend_dark_value(value);
			if (fuzzy_rules[RulesCounter].value==VALUE_TEND_LIGHT) ret_val=ret_val*tend_light_value(value);
			if (fuzzy_rules[RulesCounter].value==VALUE_LIGHT) ret_val=ret_val*light_value(value);
			/* add the rule values */
			lower_sum=lower_sum+ret_val;
		}
		RulesCounter++;
	}
	return(lower_sum);
}
//---------------------------------------------------------------------------
long FindColor(short hue, short saturation, short value, double *color_certainties) {
	double max_membership,membership;
	int color_index,res;
	if (!rules_loaded) {
		char *ColorFunctionsStart,*RulesStart;

		SetColors();
		ColorFunctionsStart=strstr(rulesfile,"color_functions:");
		RulesStart=strstr(rulesfile,"rules:");
		if (!LoadColorsFunctions(ColorFunctionsStart) || !LoadRules(RulesStart)) {
			printf("Could not load rules \n");
			return(-1);
		}
		rules_loaded=1;
	}

	max_membership = 0;
	res = COLOR_LIGHT_GREY;
	for (color_index=COLOR_WHITE;color_index<=COLOR_LIGHT_FUCIA;color_index++) {
		if (colors[color_index].color>=0) {
			membership=CalculateRules2(hue,saturation,value,color_index);
			if (color_certainties) color_certainties[color_index]=membership;
			if (membership>max_membership) {
				max_membership=membership;
				res=color_index;
			}
		}
	}
	return(res);
}
//---------------------------------------------------------------------------


