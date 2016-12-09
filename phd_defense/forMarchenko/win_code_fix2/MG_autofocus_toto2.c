/* Marchenko-Wapenaar-Satyan iterative scheme

sfmarchenko < downgoing.rsf refl=REFL_000.rsf conj=y verb=n Gtot=y niter=21 nshots=401 scale=1 eps=1e-4 shift=5 Gm=Gm.rsf G=G.rsf> Gp.rsf

	======= INPUTS ============

	p0plus.rsf: initial downgoing wavefield

	REFL_000.rsf: Fourier transform of the reflection response

	======= PARAMETERS ========

	conj  = [y]/n	- complex-conjugation of the first input (corresponds to time-reversal in time)
	verb = y/[n]	- verbosity flag
	twin  = y/[n]	- returns the timewindow as one of the outputs (window=window.rsf)
	pandq  = y/[n]	- pandq=true: returns p and q, pandq=false returns Gp and Gm
	Gtot  = y/[n]	- Gtot=true returns G=Gp+Gm
	Htot  = y/[n]	- Htot=true returns H=Gp-Gm
	niter  = 1		- number of iterations
	nshots  = 1		- number of shots in the reflection response
	scale  = 1.0	- scale factor (often due to resampling)
	eps  = 1e-4		- threshold for the timewindow
	shift  = 5		- shift in samples for the timewindow
	*/

	/*
	  Copyright (C) 2012 Colorado School of Mines
	  
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
	#include "fft1.h"
  #include "autofocus.h"

	int main(int argc, char* argv[])
	{

    Marchenko mar;


		bool verb,conj,twin,Pf1,PG;
		
	    /* OMP parameters */
		#ifdef _OPENMP
	    int ompnth;
		#endif 	
		
		float	*F_arrival, *Refl, *G, *Gm, *Gp, *f1pS, *f1p, *F1pS, *F1p;
		float	*MS, *MS1, *F1m_0, *f1m_0, *F1m, *F1m1, *f1m, *ms, *gm, *gp;
    float *MS_0, *ms_0, *ms_2, *MS2;
    float *gp1, *gm1;
		float	*window, *taper, pi;
		int		*tw;
    float *convergence;

	    /* I/O files */
	    sf_file FF_arrival;
	    sf_file FRefl;
	    sf_file FGp;
	    sf_file FGm;
	    sf_file FG;
	    sf_file Ff1m;
	    sf_file Ff1p;
	    sf_file Ftwin;
      sf_file Fconvergence;

		char *filename1, filename2[256], filename3[256];
		
		/* Cube axes */
	    sf_axis at,af,ax,at1;

	    int     nt,nf,ntr,mode,nshots,niter,len;
	    int     i,it,ix,ishot,iter,i0;
	    int		twc, twa, shift, n[2], rect[2], s[2], tap;
	    float   scale,eps,dt,df,dx,ot,of,a,b,c,d,e,f,r;

		sf_triangle tr;

	    /*------------------------------------------------------------*/
	    /* Initialize RSF parameters 								  */
	    /*------------------------------------------------------------*/
	    sf_init(argc,argv);	
		
	    /*------------------------------------------------------------*/
	    /* Initialize OMP parameters */
	    /*------------------------------------------------------------*/
		#ifdef _OPENMP
	    ompnth=omp_init();
		#endif	

		/*------------------------------------------------------------*/
		/* Flags 													  */
		/*------------------------------------------------------------*/
	    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
	    if(! sf_getbool("conj",&conj)) conj=false; /* complex conjugation (time-reversal) flag */
	    if(! sf_getbool("twin",&twin)) twin=false; /* returns the timewindow as one of the outputs */
	    if(! sf_getbool("Pf1",&Pf1)) Pf1=false; /* Htot=true: returns H=Gp-Gm */
	    if(! sf_getbool("PG",&PG)) PG=false; /* Htot=true: returns H=Gp-Gm */
	    if(! sf_getint("niter",&niter)) niter=1; /* number of iterations */
	    if(! sf_getint("nshots",&nshots)) nshots=1; /* number of shots */
	    if(! sf_getfloat("r",&r)) r=-1; /* reflection coefficient if flux 
						normalised r=-1 */
	    if(! sf_getfloat("scale",&scale)) scale=1.0; /* scale factor */
  		if(! sf_getfloat("eps",&eps)) eps=1e-4; /* threshold for the timewindow */
	  	if(! sf_getint("shift",&shift)) shift=5; /* shift in samples for the timewindow */
	  	if(! sf_getint("tap",&tap)) tap=20; /* taper of R */

	  	if (verb) {
		  	fprintf(stderr,"This program was called with \"%s\".\n",argv[0]);
			  /*fprintf(stderr,"Nr: %d Nx: %d Nt:%d\n",nr,nx,nt);*/
		  }	

	    mar = init_Marchenko( verb, conj, twin, Pf1, PG, niter, nshots, r, scale,eps, shift, tap);

      set_files(mar);


      set_axes(mar);

      set_output_files(mar);


	/*	memcpy(FA,2*nf*ntr*sizeof(float));*/

		
		/* Time-reversal flag */
		if (mar->conj) {
			mar->mode = -1;
		}
		else {
			mar->mode = +1;
		}
	    
		/* Load the reflection response into the memory */
		if (verb) fprintf(stderr,"Before loading R %d\n",2*mar->nf*mar->ntr);
		Refl = (float *)calloc(2*mar->nf*mar->ntr*mar->nshots,sizeof(float));
		/* Read REFL_000.rsf */
    mar->FRefl = sf_input("refl");
    sf_warning("reading refl");
    sf_floatread(Refl,2*mar->nf*mar->nshots*mar->ntr,mar->FRefl); 
    sf_warning("read refl");

  



    /*------------------------------------------------------------*/
    /* Allocate arrays											  */
    /*------------------------------------------------------------*/
		/* First arrival - Time */
    allocate_Arrays(mar);

		/* The three flags of fft1 are: inv, sym, and opt */

    buildWindow(mar); // window belongs to virtual source
    
    taper = buildTaper(mar); // global taper

    initSolutions(mar,Refl,taper);
    
    iterations(mar,Refl,taper);

    buildF1p(mar);

    buildGpGmG(mar,Refl,taper);

    
		sf_floatwrite(mar->gp,mar->nt2*mar->ntr,mar->FGp);
		sf_floatwrite(mar->gm,mar->nt2*mar->ntr,mar->FGm);
		if (mar->PG) {
			sf_floatwrite(mar->G,mar->nt2*mar->ntr,mar->FG);
		}
		if (mar->Pf1) {
			sf_floatwrite(mar->f1p,mar->nt*mar->ntr,mar->Ff1p);
			sf_floatwrite(mar->f1m,mar->nt*mar->ntr,mar->Ff1m);
		}
		
		if (mar->twin) {
			sf_floatwrite(mar->window,mar->nt*mar->ntr,mar->Ftwin);
  	}
    marchenko_close(mar);

    exit (0);
}



