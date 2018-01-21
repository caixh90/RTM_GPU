/* Operations with SEGY standard files. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simtab.h"

#include "segy.h"


typedef struct Segy {
    const char *name;
    unsigned int size;
} segy;

static segy bheadkey[] = {
    {"jobid",4},  /* job identification number */
    {"lino",4},	  /* line number (only one line per reel) */
    {"reno",4},	  /* reel number */
    {"ntrpr",2},  /* number of data traces per record */
    {"nart",2},	  /* number of auxiliary traces per record */
    {"hdt",2},    /* sample interval in micro secs for this reel */
    {"dto",2},    /* same for original field recording */
    {"hns",2},    /* number of samples per trace for this reel */
    {"nso",2},    /* same for original field recording */
    {"format",2}, /* data sample format code:
		     1 = floating point, 4 byte (32 bits)
		     2 = fixed point, 4 byte (32 bits)
		     3 = fixed point, 2 byte (16 bits)
		     4 = fixed point w/gain code, 4 byte (32 bits)
		     5 = IEEE floating point, 4 byte (32 bits)
		     8 = two's complement integer, 1 byte (8 bits)
		  */
    {"fold",2},	  /* CDP fold expected per CDP ensemble */
    {"tsort",2},  /* trace sorting code: 
		     1 = as recorded (no sorting)
		     2 = CDP ensemble
		     3 = single fold continuous profile
		     4 = horizontally stacked */
    {"vscode",2}, /* vertical sum code:
		     1 = no sum
		     2 = two sum ...
		     N = N sum (N = 32,767) */
    
    {"hsfs",2},	  /* sweep frequency at start */
    
    {"hsfe",2},	  /* sweep frequency at end */

    {"hslen",2},  /* sweep length (ms) */

    {"hstyp",2},  /* sweep type code:
		     1 = linear
		     2 = parabolic
		     3 = exponential
		     4 = other */
    
    {"schn",2},	  /* trace number of sweep channel */

    {"hstas",2},  /* sweep trace taper length at start if tapered (the
		     taper starts at zero time and is effective for
		     this length) */
    
    {"hstae",2},  /* sweep trace taper length at end (the ending taper
		     starts at sweep length minus the taper length at
		     end) */

    {"htatyp",2}, /* sweep trace taper type code:
		     1 = linear
		     2 = cos-squared
		     3 = other */

    {"hcorr",2},  /* correlated data traces code:
		     1 = no
		     2 = yes */

    {"bgrcv",2},  /* binary gain recovered code:
		     1 = yes
		     2 = no */

    {"rcvm",2},	  /* amplitude recovery method code:
		     1 = none
		     2 = spherical divergence
		     3 = AGC
		     4 = other */

    {"mfeet",2},  /* measurement system code:
		     1 = meters
		     2 = feet */

    {"polyt",2},  /* impulse signal polarity code:
		     1 = increase in pressure or upward
		     geophone case movement gives
		     negative number on tape
		     2 = increase in pressure or upward
		     geophone case movement gives
		     positive number on tape */

    {"vpol",2}	  /* vibratory polarity code:
		     code	seismic signal lags pilot by
		     1	337.5 to  22.5 degrees
		     2	 22.5 to  67.5 degrees
		     3	 67.5 to 112.5 degrees
		     4	112.5 to 157.5 degrees
		     5	157.5 to 202.5 degrees
		     6	202.5 to 247.5 degrees
		     7	247.5 to 292.5 degrees
		     8	293.5 to 337.5 degrees */
};

static segy segy_key[] = {
    {"tracl",  4},  /* trace sequence number within line 0 */

    {"tracr",  4},  /* trace sequence number within reel 4 */

    {"fldr",   4},  /* field record number 8 */

    {"tracf",  4},  /* trace number within field record 12 */

    {"ep",     4},  /* energy source point number 16 */

    {"cdp",    4},  /* CDP ensemble number 20 */

    {"cdpt",   4},  /* trace number within CDP ensemble 24 */

    {"trid",   2},  /* trace identification code:
                        1 = seismic data
                        2 = dead
                        3 = dummy
                        4 = time break
                        5 = uphole
                        6 = sweep
                        7 = timing
                        8 = water break
                        9---, N = optional use (N = 32,767) 28 */

    {"nvs",    2},  /* number of vertically summed traces (see
		    vscode in bhed structure) 30 */

    {"nhs",    2},  /* number of horizontally summed traces (see
		    vscode in bhed structure) 32 */

    {"duse",   2},  /* data use:
		    1 = production
		    2 = test 34 */

    {"offset", 4},  /* distance from source point to receiver
		    group (negative if opposite to direction
		    in which the line was shot) 36 */

    {"gelev",  4},  /* receiver group elevation from sea level
		    (above sea level is positive) 40 */

    {"selev",  4},  /* source elevation from sea level
		     (above sea level is positive) 44 */

    {"sdepth", 4},  /* source depth (positive) 48 */

    {"gdel",   4},  /* datum elevation at receiver group 52 */

    {"sdel",   4},  /* datum elevation at source 56 */

    {"swdep",  4},  /* water depth at source 60 */

    {"gwdep",  4},  /* water depth at receiver group 64 */

    {"scalel", 2},  /* scale factor for previous 7 entries
		     with value plus or minus 10 to the
		     power 0, 1, 2, 3, or 4 (if positive,
		     multiply, if negative divide) 68 */

    {"scalco", 2},  /* scale factor for next 4 entries
		     with value plus or minus 10 to the
		     power 0, 1, 2, 3, or 4 (if positive,
		     multiply, if negative divide) 70 */

    {"sx",     4},  /* X source coordinate 72 */

    {"sy",     4},  /* Y source coordinate 76 */

    {"gx",     4},  /* X group coordinate 80 */

    {"gy",     4},  /* Y source coordinate 84 */

    {"counit", 2},  /* coordinate units code:
		     for previoius four entries
		     1 = length (meters or feet)
		     2 = seconds of arc (in this case, the
		     X values are unsigned intitude and the Y values
		     are latitude, a positive value designates
		     the number of seconds east of Greenwich
		     or north of the equator 88 */

    {"wevel",   2},  /* weathering velocity 90 */

    {"swevel",  2},  /* subweathering velocity 92 */

    {"sut",     2},  /* uphole time at source 94 */

    {"gut",     2},  /* uphole time at receiver group 96 */

    {"sstat",   2},  /* source static correction 98 */

    {"gstat",   2},  /* group static correction 100 */

    {"tstat",   2},  /* total static applied 102 */

    {"laga",    2},  /* lag time A, time in ms between end of 240-
		      byte trace identification header and time
		      break, positive if time break occurs after
		      end of header, time break is defined as
		      the initiation pulse which maybe recorded
		      on an auxiliary trace or as otherwise
		      specified by the recording system 104 */

    {"lagb",    2},  /* lag time B, time in ms between the time
		      break and the initiation time of the energy source,
		      may be positive or negative 106 */

    {"delrt",   2},  /* delay recording time, time in ms between
		      initiation time of energy source and time
		      when recording of data samples begins
		      (for deep water work if recording does not
		      start at zero time) 108 */

    {"muts",    2},  /* mute time--start 110 */

    {"mute",    2},  /* mute time--end 112 */

    {"ns",      2},  /* number of samples in this trace 114 */

    {"dt",      2},  /* sample interval, in micro-seconds 116 */

    {"gain",    2},  /* gain type of field instruments code:
		      1 = fixed
		      2 = binary
		      3 = floating point
		      4 ---- N = optional use 118 */

    {"igc",    2},   /* instrument gain constant 120 */

    {"igi",    2},   /* instrument early or initial gain 122 */

    {"corr",   2},   /* correlated:
		      1 = no
		      2 = yes 124 */    

    {"sfs",    2},   /* sweep frequency at start 126 */

    {"sfe",    2},   /* sweep frequency at end 128 */

    {"slen",   2},   /* sweep length in ms 130 */

    {"styp",   2},   /* sweep type code:
		      1 = linear
		      2 = cos-squared
		      3 = other 132 */   

    {"stas",   2},   /* sweep trace length at start in ms 134 */

    {"stae",   2},   /* sweep trace length at end in ms 136 */

    {"tatyp",  2},   /* taper type: 1=linear, 2=cos^2, 3=other 138 */

    {"afilf",  2},   /* alias filter frequency if used 140 */

    {"afils",  2},   /* alias filter slope 142 */

    {"nofilf", 2},   /* notch filter frequency if used 144 */

    {"nofils", 2},   /* notch filter slope 146 */

    {"lcf",    2},   /* low cut frequency if used 148 */

    {"hcf",    2},   /* high cut frequncy if used 150 */

    {"lcs",    2},   /* low cut slope 152 */

    {"hcs",    2},   /* high cut slope 154 */

    {"year",   2},   /* year data recorded 156 */

    {"day",    2},   /* day of year 158 */

    {"hour",   2},   /* hour of day (24 hour clock) 160 */

    {"minute", 2},   /* minute of hour 162 */

    {"sec",    2},   /* second of minute 164 */

    {"timbas", 2},   /* time basis code:
		      1 = local
		      2 = GMT
		      3 = other 166 */   

    {"trwf",   2},   /* trace weighting factor, defined as 1/2^N
		      volts for the least sigificant bit 168 */

    {"grnors", 2},   /* geophone group number of roll switch
		      position one 170 */

    {"grnofr", 2},   /* geophone group number of trace one within
		      original field record 172 */

    {"grnlof", 2},   /* geophone group number of last trace within
		      original field record 174 */

    {"gaps",   2},   /* gap size (total number of groups dropped) 176 */

    {"otrav",  2},   /* overtravel taper code:
		      1 = down (or behind)
		      2 = up (or ahead) 71/178 */
    {"cdpx",   4},   /* X coordinate of CDP 180 */
    {"cdpy",   4},   /* Y coordinate of CDP 184 */
    {"iline",  4},   /* in-line number 188 */
    {"xline",  4},   /* cross-line number 192 */
    {"shnum",  4},   /* shotpoint number 196 */
    {"shsca",  2},   /* shotpoint scalar 200 */
    {"tval",   2},   /* trace value meas. 202 */
    {"tconst4",4},   /* transduction const 204 */
    {"tconst2",2},   /* transduction const 208 */
    {"tunits", 2},   /* transduction units 210 */
    {"device", 2},   /* device identifier 212 */
    {"tscalar",2},   /* time scalar 214 */
    {"stype",  2},   /* source type 216 */
    {"sendir", 4},   /* source energy dir. 218 */
    {"unknown",2},   /* unknown 222 */
    {"smeas4", 4},   /* source measurement 224 */
    {"smeas2", 2},   /* source measurement 228 */
    {"smeasu", 2},   /* source measurement unit 230 */
    {"unass1", 4},   /* unassigned 232 */
    {"unass2", 4}    /* unassigned 236 */
};

static void swapb(byte *x, byte *y) 
/* swap two bytes */
{
    byte tmp; 

    tmp = *x; 
    *x = *y; 
    *y = tmp;
}

static int convert2(const char* buf)
/* convert buf to 2-byte int */
{
    union {
	byte b[2];
	short s;
    } x;

    memcpy(x.b,buf,2);

	little_endian=sf_endian(); //add

    if (little_endian) swapb(x.b,x.b+1);

    return (int) x.s;
}

static void insert2(int y, char* buf)
/* convert 2-byte int to buf */
{
    union {
	byte b[2];
	short s;
    } x;

    x.s = (short) y;

	little_endian=sf_endian(); //add

    if (little_endian) swapb(x.b,x.b+1);
    
    memcpy(buf,x.b,2);
}

static int convert4(const char* buf)
/* convert buf to 4-byte int */
{
    union {
	byte b[4];
	int s;
    } x;

    memcpy(x.b,buf,4);

	little_endian=sf_endian(); //add

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }

    return x.s;
}

static float fconvert4(const char* buf)
/* convert buf to 4-byte float */
{
    union {
	byte b[4];
	float s;
    } x;

    memcpy(x.b,buf,4);

	little_endian=sf_endian(); //add

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }

    return x.s;
}

static void insert4(int y, char* buf)
/* convert 4-byte int to buf */
{
    union {
	byte b[4];
	int s;
    } x;

    x.s=y;

	little_endian=sf_endian(); //add

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }
    
    memcpy(buf,x.b,4);
}

static void finsert4(float y, char* buf)
/* convert 4-byte float to buf */
{
    union {
	byte b[4];
	float s;
    } x;

    x.s=y;

	little_endian=sf_endian(); //add

    if (little_endian) {
	swapb(x.b,x.b+3);
	swapb(x.b+1,x.b+2);
    }
    
    memcpy(buf,x.b,4);
}

void ebc2asc (int narr, char* arr)
/*< Convert char array arrr[narr]: EBC to ASCII >*/
{
    int i;
    unsigned char j;

    for (i=0; i < narr; i++) {
	j = (unsigned char) arr[i];
	arr[i] = (char) EBCtoASC[j];
    }
}

void asc2ebc (int narr, char* arr)
/*< Convert char array arrr[narr]: ASCII to EBC >*/
{
    int i;
    unsigned char j;

    for (i=0; i < narr; i++) {
	j = (unsigned char) arr[i];
	arr[i] = (char) ASCtoEBC[j];
    }
}

int segyformat (const char* bhead)
/*< extracts SEGY format from binary header >*/
{
    return convert2(bhead+SF_SEGY_FORMAT);
}

void set_segyformat (char* bhead, int format)
/*< set SEGY format in binary header >*/
{
    insert2(format,bhead+SF_SEGY_FORMAT);
}

int segyns (const char* bhead)
/*< extracts ns (number of samples) from binary header >*/
{
    return convert2(bhead+SF_SEGY_NS);
}

void set_segyns(char* bhead, int ns)
/*< set ns (number of samples) in binary header >*/
{
    insert2(ns,bhead+SF_SEGY_NS);
}

float segydt (const char* bhead)
/*< extracts dt (sampling) from binary header >*/
{
    return (float) (convert2(bhead+SF_SEGY_DT)/1000000.);
}

void set_segydt(char* bhead, float dt)
/*< set dt (sampling) in binary header >*/    
{
    float scale;

    /* input in seconds or miliseconds? */
    scale = (dt < 1.0)? 1000000.:1000.;

    insert2((int) (scale*dt),bhead+SF_SEGY_DT);
}
int segyntrace (const char* bhead)
/*< extracts ntrace (number of traces for one source) from binary header >*/
{
	return (float) convert2(bhead+SF_SEGY_NTRC);
}

void set_segyntrace(char* bhead, int ntrace)
/*< set ntrace (number of traces for one source) in binary header >*/
{
	insert2(ntrace,bhead+SF_SEGY_NTRC);
}

static void float2ibm (float y, char* num)
/* floating-point conversion to IBM format */
{
    unsigned int x, s, f;
    int e;
    const unsigned int fMAXIBM = 0x7FFFFFFF;

    memcpy (&x,&y,4);

    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) {
	insert4(x,num);
	return; 
    }

    /* fetch the sign, exponent (removing excess 127), and fraction */
    s =   x & 0x80000000;
    e = ((x & 0x7f800000) >> 23) - 127;
    f =   x & 0x007fffff;

    /* convert 23 bit fraction to 24 bit fraction */
    f <<= 1; 

    /* restore the '1' preceeded the IEEE binary point */
    f |= 0x01000000; 
    
    /* convert scale factor from base-2 to base-16 */
    if (e >= 0) {
      f <<= (e & 3); 
      e >>= 2;
    } else {
	f >>= ((-e) & 3); 
	e = -((-e) >> 2);
    }

    /* reduce fraction to 24 bits */
    if (f & 0x0f000000) {
	f >>= 4;
	e += 1;
    }

    /* convert exponent to excess 64 and store the number */
    if ((e += 64) > 127) {
	s |= fMAXIBM;
    } else if (e >= 0) {
	s |= (e << 24) | f;
    }

    insert4(s,num);
}

static float ibm2float (const char* num)
/* floating point conversion from IBM format */
{
    unsigned int x, s, f;
    const unsigned int fMAXIEEE = 0x7F7FFFFF;
    int e;         
    float y;
                                                                     
    x = convert4 (num);
    
    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) return 0.0; 

    /* fetch the sign, exponent (removing excess 64), and fraction */   
    s =   x & 0x80000000;                                               
    e = ((x & 0x7f000000) >> 24) - 64;                                   
    f =   x & 0x00ffffff;                                                
                                                                    
    /* convert scale factor from base-16 to base-2 */        
    if (e >= 0) {
	e <<= 2;  
    } else { 
	e = -((-e) << 2); 
    }
                                                                        
    /* convert exponent for 24 bit fraction to 23 bit fraction */           
    e -= 1;                                                               
                                                                            
    /* normalize the fraction */                                            
    if (0 != f) {
	while ((f & 0x00800000) == 0) {         
	    f <<= 1;
	    e -= 1;
	}       
    }                                                               
	
    /* drop the '1' preceeding the binary point */                       
    f &= 0x007fffff;                                                         
    
    /* convert exponent to excess 127 and store the number */
    if ((e += 127) >= 255) {
	s |= fMAXIEEE;
    } else if (e > 0) {
	s |= (e << 23) | f; 	    
    }    

    memcpy (&y,&s,4);
    return y;
}

void segy2trace(const char* buf, float* trace, int ns, int format)
/*< Extract a floating-point trace[nt] from buffer buf.
---
format: 1: IBM, 2: int4, 3: int2, 5: IEEE
>*/
{
    int i, nb;

    nb = (3==format)? 2:4;

    for (i=0; i < ns; i++, buf += nb) {
	switch (format) {
	    case 1: trace[i] = ibm2float (buf);       break; /* IBM float */
	    case 2: trace[i] = (float) convert4(buf); break; /* int4 */
	    case 3: trace[i] = (float) convert2(buf); break; /* int2 */
	    case 5: trace[i] = fconvert4(buf);        break; /* IEEE float */	  
		//default: sf_error("Unknown format %d",format); break;
		default: printf("Unknown format %d",format); break;
	}
    }
}

void trace2segy(char* buf, const float* trace, int ns, int format)
/*< Convert a floating-point trace[ns] to buffer buf.
---
format: 1: IBM, 2: int4, 3: int2, 5: IEEE
>*/
{
    int i, nb;

    nb = (3==format)? 2:4;

    for (i=0; i < ns; i++, buf += nb) {
      switch (format) {
      case 1: float2ibm(trace[i],buf);     break; /* IBM float */
      case 2: insert4((int) trace[i],buf); break; /* int4 */
      case 3: insert2((int) trace[i],buf); break; /* int2 */
      case 5: finsert4((float) trace[i],buf); break; /* IEEE float */
      //default: sf_error("Unknown format %d",format); break;
	  default: printf("Unknown format %d",format); break;
      }
    }
}

void segy2head(const char* buf, int* trace, int nk)
/*< Create an integer trace header trace[nk] from buffer buf >*/
{
    int i, byte;
    const char *buf0, *bufi;
	
    buf0 = buf;
    for (i=0; i < nk; i++) {
	/* allow to remap header keys */
	//bufi = sf_getint(segy_key[i].name,&byte)? buf0+byte:buf;
	bufi=buf;

	switch (segy_key[i].size) {
	    case 2:
		trace[i] = convert2(bufi);
		buf += 2;
		break;
	    case 4:
		trace[i] = convert4(bufi);
		buf += 4;
		break;
	    default:
		//sf_error("Unknown size %d",segy_key[i].size);
		printf("Unknown size %d",segy_key[i].size);
		break;
	}
    }
}

int segykey (const char* key) 
/*< Extract a SEGY key value >*/
{
    
    int i;

    for (i=0; i < SF_NKEYS; i++) {
	if (0==strcmp(key,segy_key[i].name)) return i;
    }
    //sf_error("no such key %s",key);
	printf("no such key %s",key);
    return 0;
}

const char* segykeyword (int k) 
/*< Find a SEGY key from its number >*/
{
    return segy_key[k].name;
}

void head2segy(char* buf, const int* trace, int nk)
/*< Convert an integer trace[nk] to buffer buf >*/
{
    int i;

    for (i=0; i < nk; i++) {
	if (i < SF_NKEYS && 2 == segy_key[i].size) {
	    insert2(trace[i],buf);
	    buf += 2;
	} else {
	    insert4(trace[i],buf);
	    buf += 4;
	}
    }
}

void binary_head(char* buf)
/*< Create a binary header for SEGY >*/
{
    int i, val, size;

    for (i=0; i < 3; i++) {
	if (!sf_getint(bheadkey[i].name,&val)) val=1;
	insert4(val,buf);
	buf += 4;
    }

    for (i=3; i < SF_BHKEYS; i++) {
	size = bheadkey[i].size;
	if (sf_getint(bheadkey[i].name,&val)) 
	    insert2(val,buf);
	buf += 2;
    }

}

bool endian (void)
/*< Set endianness >*/
{
    little_endian = sf_endian();    
    return little_endian;
}

bool sf_endian (void)
/*< Endianness test, returns true for little-endian machines >*/
{
    bool little_endian;
	
    union {
		unsigned char c[4];
		int i;
    } test;
	
    test.i=0;
    test.c[0] = (unsigned char) 1;
    
    assert (2 == sizeof(short) && 4 == sizeof(int)); /* fix this later */
    little_endian = (bool) (0 != (test.i << 8));
    
    return little_endian;
}
/* 	$Id: segy.c 7107 2011-04-10 02:04:14Z ivlad $	 */

static sf_simtab pars;

bool sf_getint (const char* key,/*@out@*/ int* par) 
/*< get an int parameter from the command line >*/
{
    return sf_simtab_getint (pars,key,par);
}


