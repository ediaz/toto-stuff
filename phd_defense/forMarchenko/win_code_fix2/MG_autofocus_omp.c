/* Marchenko-Wapenaar-Satyan iterative scheme

sfmarchenko < downgoing.rsf refl=REFL_000.rsf conj=y verb=n Gtot=y
               niter=21 nshots=401 scale=1 eps=1e-4 shift=5 Gm=Gm.rsf G=G.rsf> Gp.rsf

	======= INPUTS ============

	p0plus.rsf: initial downgoing wavefield

	REFL_000.rsf: Fourier transform of the reflection response 

	======= PARAMETERS ========

	conj  = [y]/n	- complex-conjugation of the first input (corresponds to
                  time-reversal in time)
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

    Marchenko *mar;


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

		/* Cube axes */
	    sf_axis at,af,ax,at1,axf,azf;

	    int     nt,nf,ntr,mode,nshots,niter,len;
	    int     i,it,ix,ishot,iter,i0;
	    int		twc, twa, shift, n[2], rect[2], s[2], tap;
	    float   scale,eps,dt,df,dx,ot,of,a,b,c,d,e,f,r;


	    /*------------------------------------------------------------*/
	    /* Initialize RSF parameters 								  */
	    /*------------------------------------------------------------*/
	    sf_init(argc,argv);	
		
	    /*------------------------------------------------------------*/
	    /* Initialize OMP parameters */
	    /*------------------------------------------------------------*/
		#ifdef _OPENMP
	    ompnth=omp_init();
    #else
      ompnth=1;
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
		  }	

      mar = (Marchenko*) sf_alloc(ompnth,sizeof(Marchenko));

      for (int ithread=0; ithread<ompnth; ithread++ )
	      mar[ithread] = init_Marchenko( verb, conj, twin, Pf1, 
                                       PG, niter, nshots, r, 
                                       scale,eps, shift, tap);

      set_files(mar[0]);
      set_axes(mar[0]);
      for (int ithread=0; ithread<ompnth; ithread++ )
	      set_axes_copy(mar[ithread],mar[0]);




      set_output_files(mar[0]);
  
      azf = sf_iaxa(mar[0]->FF_arrival,3);
      axf = sf_iaxa(mar[0]->FF_arrival,4);

      sf_oaxa(mar[0]->FF_arrival,axf,4);
      sf_oaxa(mar[0]->FF_arrival,azf,3);


		/* Time-reversal flag */
		if (mar[0]->conj) {
			mar[0]->mode = -1;
		}
		else {
			mar[0]->mode = +1;
		}
	    
		/* Load the reflection response into the memory */
		if (verb) fprintf(stderr,"Before loading R %d\n",2*mar[0]->nf*mar[0]->ntr);
		Refl = (float *)calloc(2*mar[0]->nf*mar[0]->ntr*mar[0]->nshots,sizeof(float));
		/* Read REFL_000.rsf */
    mar[0]->FRefl = sf_input("refl");
    sf_warning("reading refl");
    sf_floatread(Refl,2*mar[0]->nf*mar[0]->nshots*mar[0]->ntr,mar[0]->FRefl); 
    sf_warning("read refl");

    int nsources = (sf_n(axf)*sf_n(azf));
    int ngroups = SF_MIN(nsources/(ompnth)+1,nsources);

    float *input = sf_floatalloc(mar[0]->nt*mar[0]->ntr*ompnth);
    
    for (int ithread=0; ithread<ompnth; ithread++){
      just_allocate_Arrays(mar[ithread]);
    }

    sf_warning("ngroups=%d nsources=%d ",ngroups,nsources);
    taper = buildTaper(mar[0]); // global taper
  
    for (int group=0; group < ngroups; ++group){
      int ntraces = SF_MIN(ompnth,nsources-group*ompnth);
      // read the first arrival for ntraces virtual sources:
      sf_floatread(&input[0],mar[0]->nt*mar[0]->ntr*ntraces,mar[0]->FF_arrival);


      int trace=0; 
      // 
			#ifdef _OPENMP
			#pragma omp parallel for private(trace) \
		  	shared(mar,Refl,taper)
			#endif
      for (trace=0; trace<ntraces; ++trace){
          sf_warning("procesing trace %d of %d, batch %d",trace, ntraces,group);

          setFirstArrival(mar[trace],&input[trace*mar[0]->nt*mar[0]->ntr]);

          buildWindow(mar[trace]); // window belongs to virtual source
          taper = buildTaper(mar[trace]); // global taper
          initSolutions(mar[trace],Refl,taper);
          iterations(mar[trace],Refl,taper);
          buildF1p(mar[trace]);
          buildGpGmG(mar[trace],Refl,taper);
      }


      // serial writing:
      for (int trace=0; trace<ntraces; ++trace){
          sf_warning("writing trace=%d ",trace);
    		  sf_floatwrite(mar[trace]->gp,mar[0]->nt2*mar[0]->ntr,mar[0]->FGp);
    		  sf_floatwrite(mar[trace]->gm,mar[0]->nt2*mar[0]->ntr,mar[0]->FGm);
    		  if (mar[0]->PG) {
    		  	sf_floatwrite(mar[trace]->G,mar[0]->nt2*mar[0]->ntr,mar[0]->FG);
    		  }
    		  if (mar[0]->Pf1) {
    		  	sf_floatwrite(mar[trace]->f1p,mar[0]->nt*mar[0]->ntr,mar[0]->Ff1p);
    		  	sf_floatwrite(mar[trace]->f1m,mar[0]->nt*mar[0]->ntr,mar[0]->Ff1m);
    		  }
    		  if (mar[0]->twin) {
    		  	sf_floatwrite(mar[trace]->window,mar[0]->nt*mar[0]->ntr,mar[0]->Ftwin);
      	  }
         // zeroArrays(mar[trace]);
          sf_warning("writing done trace=%d",trace);
      }
    }

    sf_warning("closing struct");
    for (int ithread=0; ithread<ompnth; ithread++)
      marchenko_close(mar[ithread]);

    close_files(mar[0]);
    exit (0);
}
