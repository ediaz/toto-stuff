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

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int nt, ntopt;
static float dt,ot;

static int nw;
static float dw,ow;
static bool sym;


void fft1_init (int n, float d, float o, bool opt, bool symin)
/*< >*/
{
  nt = n;
  dt = d;
  ot = o;

  if (opt){
	  ntopt = 2*kiss_fft_next_fast_size((n+1)/2);
  } else {
    ntopt = nt;
  }
  int n1;
  n1 = ntopt;
	if (nt%2) n1=ntopt+1;
	nw = n1/2+1;
	dw = 1./(ntopt*d);
  sym = symin;
}



void fft1 (float *field, float *Field, sf_file in, bool inv, bool sym, bool opt)
/*< Fast Fourier Transform along the first axis >*/
{
    /*bool inv, sym, opt;*/
    int n1, nt, nw, i1, i2, n2;
    float dw, *p, d1, o1, wt, shift;
    kiss_fft_cpx *pp, ce;
    char *label;
    sf_file out=NULL;
    kiss_fftr_cfg cfg;
	bool verb;

	verb=1;

	if (verb) fprintf(stderr,"Beginning of fft1. inv:%d sym: %d opt:%d \n",inv,sym,opt);
    
	/*sf_init(argc, argv);
    in  = sf_input("in");
    out = sf_output("out");*/

    /*if (!sf_getbool("inv",&inv)) inv=false;*/
    /* if y, perform inverse transform */
    /*if (!sf_getbool("sym",&sym)) sym=false;*/
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */
    /*if (!sf_getbool("opt",&opt)) opt=true;*/
    /* if y, determine optimal size for efficiency */
	
    
	if (inv) {
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	/*sf_settype (out,SF_FLOAT);*/
    } else {
	if (SF_FLOAT   != sf_gettype(in)) sf_error("Need float input");
    	/*sf_settype (out,SF_COMPLEX);*/
    }

    n2 = sf_leftsize(in,1);

    if (!inv) {
	if (!sf_histint  (in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (verb) fprintf(stderr,"n1: %d n2: %d nt: %d nw: %d dw:%f o1:%f\n",n1,n2,nt,nw,dw,o1);
	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);
  sf_warning(" Fil NW=%d dw=%g d1=%g",nw,dw,d1);


	if (verb) fprintf(stderr,"n1: %d n2: %d nt: %d nw: %d d1: %f dw:%f o1:%f\n",n1,n2,nt,nw,d1,dw,o1);

	/*sf_putint  (out,"n1",nw);
	sf_putfloat(out,"o1",0.);
	sf_putfloat(out,"d1",dw);

	sf_putfloat(out,"fft_o1",o1);
	sf_putfloat(out,"fft_n1",n1);*/

	/* fix label */
	/*if (NULL != (label = sf_histstring(in,"label1"))) {
	    sf_putstring(out,"fft_label1",label);
	    if (!sf_fft_label(1,label,out))
		sf_putstring(out,"label1","Wavenumber");
	}*/
    } else {
	if (!sf_histint  (in,"n1",&nw)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"fft_o1",&o1)) o1=0.; 

	nt = 2*(nw-1);
	d1 = 1./(nt*dw);

	if (!opt || !sf_histint  (in,"fft_n1",&n1)) n1 = nt;
	
	if (verb) fprintf(stderr,"n1: %d n2: %d nt: %d nw: %d d1: %f dw:%f o1:%f\n",n1,n2,nt,nw,d1,dw,o1);
	
	/*sf_putint  (out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);*/

	/* fix label */
	/*if (NULL != (label = sf_histstring(in,"fft_label1"))) {
	    sf_putstring(out,"label1",label);
	} else if (NULL != (label = sf_histstring(in,"label1"))) {
	    (void) sf_fft_label(1,label,out);
	}*/

    }	
    
	/*sf_fft_unit(1,sf_histstring(in,"unit1"),out);*/

    p = sf_floatalloc(nt);
    pp = (kiss_fft_cpx*) sf_complexalloc(nw);

    cfg = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
    wt = sym? 1./sqrtf((float) nt): 1.0/nt;

    
  sf_warning("inv=%d wt=%g n1=%d nt=%d nw=%d ot=%g",1,wt,n1,nt,nw,ot);







	for (i2=0; i2 < n2; i2++) {
	if (!inv) {
	    /*sf_floatread (p,n1,in);*/
	    memcpy(p,&field[i2*n1],n1*sizeof(float));

	    if (sym) {
		for (i1=0; i1 < n1; i1++) {
		    p[i1] *= wt;
		}
	    }

	    for (i1=n1; i1 < nt; i1++) {
		p[i1]=0.0;
	    }
	    kiss_fftr (cfg,p,pp);

	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    shift = -2.0*SF_PI*i1*dw*o1;
		    ce.r = cosf(shift);
		    ce.i = sinf(shift);
		    pp[i1]=sf_cmul(pp[i1],ce);
		}
	    }
		
		memcpy(&Field[i2*2*nw],pp,2*nw*sizeof(float));

	    /*sf_floatwrite((float*) pp,2*nw,out);*/
	} else {
	    /*sf_floatread((float*) pp,2*nw,in);*/
		memcpy(pp,&field[i2*2*nw],2*nw*sizeof(float));

	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    shift = +2.0*SF_PI*i1*dw*o1;
		    ce.r = cosf(shift);
		    ce.i = sinf(shift);
		    pp[i1]=sf_cmul(pp[i1],ce);
		}
	    }

	    kiss_fftri(cfg,pp,p);

	    for (i1=0; i1 < n1; i1++) {
		p[i1] *= wt;
	    }

		memcpy(&Field[i2*n1],p,n1*sizeof(float));
	    /*sf_floatwrite (p,n1,out);*/
	}
    }

	free(p);
	free(pp);

    
	/*exit (0);*/
}

/* 	$Id: Mfft1.c 7107 2011-04-10 02:04:14Z ivlad $	 */




void fft1_2D (float *field, float *Field, int n2, bool inv)
/*< Fast Fourier Transform along the first axis >*/
{
    /*bool inv, sym, opt;*/
    int n1, i1, i2;
    float  *p,  wt, shift;
    kiss_fft_cpx *pp, ce;
    char *label;
    sf_file out=NULL;
    kiss_fftr_cfg cfg;
	bool verb;

	verb=1;

    
    
  p = sf_floatalloc(ntopt);
  pp = (kiss_fft_cpx*) sf_complexalloc(nw);

  cfg = kiss_fftr_alloc(ntopt,inv?1:0,NULL,NULL);
  wt = sym? 1./sqrtf((float) ntopt): 1.0/ntopt;

    
  if(!inv) sf_warning("gforward fft");
  if(inv) sf_warning("inverse fft");

	for (i2=0; i2 < n2; i2++) {

    	if(inv) {
    	    /*sf_floatread((float*) pp,2*nw,in);*/
    		  memcpy(pp,&field[i2*2*nw],2*nw*sizeof(float));
    
    	    if (0. != ot) {
    		    for (i1=0; i1 < nw; i1++) {
    		      shift = +2.0*SF_PI*i1*dw*ot;
    		      ce.r = cosf(shift);
    		      ce.i = sinf(shift);
    		      pp[i1]=sf_cmul(pp[i1],ce);
    		    }
    	    }
    
    	    kiss_fftri(cfg,pp,p);
    
    	    for (i1=0; i1 < n1; i1++) {
    		    p[i1] *= wt;
    	    }
    
    		  memcpy(&Field[i2*n1],p,n1*sizeof(float));
    	    /*sf_floatwrite (p,n1,out);*/
    	}else{
    	    /*sf_floatread (p,n1,in);*/
          n1 = ntopt;
          sf_warning("i2=%d n1=%d",i2,n1); 
    	    memcpy(p,&field[i2*n1],n1*sizeof(float));
    	    if (sym) {
        		for (i1=0; i1 < n1; i1++) {
      		    p[i1] *= wt;
    	    	}
    	    }
    
    	    for (i1=n1; i1 < ntopt; i1++) {
    		    p[i1]=0.0; // padding
    	    }
    
    	    kiss_fftr (cfg,p,pp);
    
    	    if (0. != ot) {
    		    for (i1=0; i1 < nw; i1++) {
    		      shift = -2.0*SF_PI*i1*dw*ot;
    		      ce.r = cosf(shift);
    		      ce.i = sinf(shift);
    		      pp[i1]=sf_cmul(pp[i1],ce);
    		    }
    	    }
    		
    		  memcpy(&Field[i2*2*nw],pp,2*nw*sizeof(float));
      }
    
    }

	free(p);
	free(pp);

    
	/*exit (0);*/
}




void fft1_2D_fwd (float *input, float *output, int n2)
/*< Fast Fourier Transform along the first axis forward>*/
{
    /*bool inv, sym, opt;*/
    int n1, i1, i2;
    float  *p,  wt, shift;
    kiss_fft_cpx *pp, ce;
    char *label;
    sf_file out=NULL;
    kiss_fftr_cfg cfg;
	bool verb;

	verb=1;

    
    
  wt = sym? 1./sqrtf((float) ntopt): 1.0/ntopt;

  p = sf_floatalloc(ntopt);
  pp = (kiss_fft_cpx*) sf_complexalloc(nw);

  cfg = kiss_fftr_alloc(ntopt,0,NULL,NULL);
    

	for (i2=0; i2 < n2; i2++) {

    /*sf_floatread (p,n1,in);*/
    memcpy(p,&input[i2*nt],nt*sizeof(float));
    if (sym) {
      for (i1=0; i1 < nt; i1++) {
    	  p[i1] *= wt;
      }
    }
    for (i1=nt; i1 < ntopt; i1++) {
      p[i1]=0.0; // padding
    }
    
      kiss_fftr (cfg,p,pp);
    
    if (0. != ot) {
      for (i1=0; i1 < nw; i1++) {
          shift = -2.0*SF_PI*i1*dw*ot;
          ce.r = cosf(shift);
          ce.i = sinf(shift);
          pp[i1]=sf_cmul(pp[i1],ce);
      }
    }
    
    memcpy(&output[i2*2*nw],pp,2*nw*sizeof(float));
  }
    free(p);
    free(pp);

	/*exit (0);*/
}


void fft1_2D_inv (float *input, float *output, int n2)
/*< Fast Fourier Transform along the first axis >*/
{
    /*bool inv, sym, opt;*/
    int n1, i1, i2;
    float  *p,  wt, shift;
    kiss_fft_cpx *pp, ce;
    char *label;
    sf_file out=NULL;
    kiss_fftr_cfg cfg;
	bool verb;

	verb=1;

    
    
  p = sf_floatalloc(ntopt);
  pp = (kiss_fft_cpx*) sf_complexalloc(nw);

  cfg = kiss_fftr_alloc(ntopt,1,NULL,NULL);
  wt = sym? 1./sqrtf((float) ntopt): 1.0/ntopt;


	for (i2=0; i2 < n2; i2++) {
	    /*sf_floatread((float*) pp,2*nw,in);*/
		  memcpy(pp,&input[i2*2*nw],2*nw*sizeof(float));


		  for (i1=0; i1 < nw; i1++) {
        shift = +2.0*SF_PI*i1*dw*ot;
        ce.r = cosf(shift);
        ce.i = sinf(shift);
		    pp[i1]=sf_cmul(pp[i1],ce);
		  }

	    kiss_fftri(cfg,pp,p);

	    for (i1=0; i1 < ntopt; i1++) {
		    p[i1] *= wt;
	    }

		  memcpy(&output[i2*nt],p,nt*sizeof(float));
  }

	free(p);
	free(pp);

    
	/*exit (0);*/
}
