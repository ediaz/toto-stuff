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


	int main(int argc, char* argv[])
	{

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
		

      Fconvergence = sf_output("conv");
      convergence = sf_floatalloc(niter);
      sf_putint (Fconvergence, "n1", niter);
      sf_putint (Fconvergence, "n2", 1);
      sf_putfloat (Fconvergence, "d1", 1.0f);
      sf_putfloat (Fconvergence, "o1", 0.0f);


	    /*------------------------------------------------------------*/
	    /* I/O files 												  */
	    /*------------------------------------------------------------*/	
		/* "in" is the transposed version of p00plus_xxxx_xxxx.rsf
		   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=ntr,n2=nt
		   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=nt,n2=ntr */
		FF_arrival = sf_input("in");
		
		/* refl is REFL_000.rsf
		   It is used to read nf, df, of
		  Dimensions are: n1=nf,n2=ntr */
		/*FRefl = (sf_file)sf_alloc(1,sizeof(sf_file));*/	
		FRefl = sf_input("refl");

		FGp = sf_output("out");
		FGm = sf_output("Gm");
		
		if (PG) {
			FG  = sf_output("G");
		}
		if (Pf1) {
			Ff1p  = sf_output("f1p");
			Ff1m  = sf_output("f1m");
		}
		
		if (twin) {
			Ftwin = sf_output("window"); /* time window */
		}

   
		







		/*------------------------------------------------------------*/
		/* Axes */
		/*------------------------------------------------------------*/    
		at = sf_iaxa(FF_arrival,1); sf_setlabel(at,"Time"); if(verb) sf_raxa(at); /* time */
		at1 = sf_iaxa(FF_arrival,1); sf_setlabel(at,"Time"); if(verb) sf_raxa(at); /* time */
		af = sf_iaxa(FRefl,1); sf_setlabel(af,"Frequency"); if(verb) sf_raxa(af); /* frequency */
		ax = sf_iaxa(FF_arrival,2); sf_setlabel(ax,"r"); if(verb) sf_raxa(ax); /* space */
	    
		nt = sf_n(at); dt = sf_d(at); ot = sf_o(at);
		nf = sf_n(af); df = sf_d(af); of = sf_o(af);
		ntr = sf_n(ax); dx = sf_d(ax);

    int nt2=(nt/2);
    sf_setn(at1,nt2);


		if (verb) fprintf(stderr,"nt: %d nf: %d ntr:%d\n",nt,nf,ntr);

		sf_fileclose(FRefl);

	    /*------------------------------------------------------------*/
	    /* Setup output data and wavefield header					  */
	    /*------------------------------------------------------------*/
		sf_oaxa(FGp,at1,1);
		sf_oaxa(FGp,ax,2);
		sf_oaxa(FGm,at1,1);
		sf_oaxa(FGm,ax,2);
		
		if (PG) {
			sf_oaxa(FG,at1,1);
			sf_oaxa(FG,ax,2);
		}
		if (Pf1) {
			sf_oaxa(Ff1p,at,1);
			sf_oaxa(Ff1p,ax,2);
			sf_oaxa(Ff1m,at,1);
			sf_oaxa(Ff1m,ax,2);
		}
		if (twin) {
			sf_oaxa(Ftwin,at,1);
			sf_oaxa(Ftwin,ax,2);
		}

	    /*------------------------------------------------------------*/
	    /* Allocate arrays											  */
	    /*------------------------------------------------------------*/
		/* First arrival - Time */
		F_arrival = (float *)calloc(nt*ntr,sizeof(float));
		sf_floatread(F_arrival,nt*ntr,FF_arrival);
		ms   = (float *)calloc(nt*ntr,sizeof(float));
		ms_0 = (float *)calloc(nt*ntr,sizeof(float));
		ms_2 = (float *)calloc(nt*ntr,sizeof(float));
		f1m_0= (float *)calloc(nt*ntr,sizeof(float));
		f1m  = (float *)calloc(nt*ntr,sizeof(float));
		f1pS = (float *)calloc(nt*ntr,sizeof(float));
		f1p  = (float *)calloc(nt*ntr,sizeof(float));
		memcpy(ms,F_arrival,nt*ntr*sizeof(float));

		/* Allocate for coda M of f2 - Frequency */
		MS   = (float *)calloc(2*nf*ntr,sizeof(float));
		MS_0 = (float *)calloc(2*nf*ntr,sizeof(float));
		MS1  = (float *)calloc(2*nf*ntr,sizeof(float));
		MS2  = (float *)calloc(2*nf*ntr,sizeof(float));
		F1m_0= (float *)calloc(2*nf*ntr,sizeof(float));
		F1m  = (float *)calloc(2*nf*ntr,sizeof(float));
		F1m1 = (float *)calloc(2*nf*ntr,sizeof(float));
		F1pS = (float *)calloc(2*nf*ntr,sizeof(float));
		F1p  = (float *)calloc(2*nf*ntr,sizeof(float));
		/* The three flags of fft1 are: inv, sym, and opt */

    fft1_init (nt, dt, ot, true, false);
    fft1_2D_fwd(F_arrival,MS,ntr);


	/*	memcpy(FA,2*nf*ntr*sizeof(float));*/

      fprintf(stderr,"nt2 is %d\n",nt2);
		/* Output wavefields */
		G = (float *)calloc(nt2*ntr,sizeof(float));
		gp1= (float *)calloc(nt*ntr,sizeof(float));
		gp= (float *)calloc(nt2*ntr,sizeof(float));
		gm1= (float *)calloc(nt*ntr,sizeof(float));
		gm= (float *)calloc(nt2*ntr,sizeof(float));
		Gp= (float *)calloc(2*nf*ntr,sizeof(float));
		Gm= (float *)calloc(2*nf*ntr,sizeof(float));


		
		/* Time-reversal flag */
		if (conj) {
			mode = -1;
		}
		else {
			mode = +1;
		}
	    
		/* Load the reflection response into the memory */
		if (verb) fprintf(stderr,"Before loading R %d\n",2*nf*ntr);
		Refl = (float *)calloc(2*nf*ntr*nshots,sizeof(float));
		/* Read REFL_000.rsf */
    FRefl = sf_input("refl");
    sf_warning("reading refl");
    sf_floatread(Refl,2*nf*nshots*ntr,FRefl); 
    sf_warning("read refl");

		/* Build time-window */
		tw = (int *)calloc(ntr,sizeof(int));
		window = (float *)calloc(nt*ntr,sizeof(float));
		/*memset(window,0,nt*ntr*sizeof(float));*/
	    /* I am not sure why I set it to this value */
		/*for (ix=0; ix<ntr; ix++) {
			tw[ix] = nt*dt+ot+0.15; 
		}*/
		
		if (verb) fprintf(stderr,"---> Build time-window?\n");
    // checking time sample corresponding to muting time
		for (ix=0; ix<ntr; ix++) {
      for (it=0; it<nt; it++) {
        if ((F_arrival[it+ix*nt]*F_arrival[it+ix*nt])>eps*eps) {
          /*tw[ix] = it*dt+ot;*/
          tw[ix] = it;
          break;
        }
			}
		}
		if (verb) fprintf(stderr,"---> Build time-window1\n");
		for (ix=0; ix<ntr; ix++) {
			twc = (int)(tw[ix]-shift-10);
			twa = (int)(-twc+nt);
			/*if (verb) fprintf(stderr,"%d %d\n",twc,twa);*/
			for (it=0; it<nt; it++) {
			/*	if ((it>twa) || (it<twc)) {*/
				if ((it>twa) && (it<twc)) { 
					window[it+ix*nt] = 1.0; // building windowing function W from Filippo's paper
				}
			}
		}

		if (verb) fprintf(stderr,"---> Build time-window2\n");
		/* Smoothing of the window */
		/* Should I implement flags for rect and iter? */
		/* Look at Msmooth.c to understand below */
		n[0] = nt;
		n[1] = ntr;
		s[0] = 1;
		s[1] = nt;
		rect[0] = 5;
		rect[1] = 5;

		for (ix=0; ix <= 1; ix++) {
			if (rect[ix] <= 1) continue;
			tr = sf_triangle_init (rect[ix],n[ix],false);
			for (it=0; it < (nt*ntr/n[ix]); it++) {
				i0 = sf_first_index (ix,it,1+1,n,s);
				for (iter=0; iter < 2; iter++) {
					sf_smooth2 (tr,i0,s[ix],false,window );
				}
			}
			sf_triangle_close(tr);
		}
		if (verb) fprintf(stderr,"---> Here\n");
		
		/* Tapering */
		pi = 4.0*atan(1.0);
		
		taper = (float *)calloc(ntr,sizeof(float));
		memset(taper,0,ntr*sizeof(float));

		for (ix=0; ix<tap; ix++) {
			taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/(2*tap))));
			taper[ntr-ix-1] = taper[ix];
		}
		for (ix=tap; ix<(ntr-tap); ix++) {
			taper[ix] = 1.0;
		}
		if (verb) fprintf(stderr,"---> taper finish\n");
		
		FRefl = sf_input("refl");

		/*------------------------------------------------------------*/
		/* Loop over iterations */
		/*------------------------------------------------------------*/
		if (verb) fprintf(stderr,"---> Begin to iterative solve for f1p and f1m\n");
		/*starting iteration for f1m */
		memset(F1m_0,0,2*nf*ntr*sizeof(float));

			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d) \
					shared(MS,taper,Refl,F1m_0,MS2)
				#endif
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = MS[ishot*2*nf+it];
						d = MS[ishot*2*nf+it+1];

						F1m_0[ix*2*nf+it]   += (a*c - mode*b*d);
						F1m_0[ix*2*nf+it+1] += (mode*a*d + b*c);
						MS2  [ix*2*nf+it]   += r*(a*c - b*d); // rTd* R
						MS2  [ix*2*nf+it+1] += r*(a*d + b*c);

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      // Inverse fft
      fft1_2D_inv (F1m_0, f1m_0,ntr);
      fft1_2D_inv (MS2, ms_2,ntr);


			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it) \
				shared(f1m, f1m_0, window,ms_2)
			#endif
			/* window to get f1m_0 */
			for (ix=0; ix<ntr; ix++) {
				#pragma ivdep
				for (it=0; it<nt; it++) {
					f1m_0[it+ix*nt] = scale*window[it+ix*nt]*f1m_0[it+ix*nt];  
					ms_2[it+ix*nt] = scale*window[it+ix*nt]*ms_2[it+ix*nt];  
					//f1m_0[it+ix*nt] = scale*f1m_0[it+ix*nt];  
					f1m[it+ix*nt] = f1m_0[it+ix*nt];
				}	
			}
			
      fft1_2D_fwd (f1m, F1m,ntr);


	/* initialise MS the coda for f1+ */
		memset(MS_0,0,2*nf*ntr*sizeof(float));
		memset(MS,0,2*nf*ntr*sizeof(float));

			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d) \
					shared(MS_0,taper,Refl,F1m)
				#endif
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = F1m[ishot*2*nf+it];
						d = F1m[ishot*2*nf+it+1];

						MS_0[ix*2*nf+it]   += (a*c - mode*b*d);
						MS_0[ix*2*nf+it+1] += (mode*a*d + b*c);

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      fft1_2D_inv (MS_0,ms_0,ntr);


			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it) \
				shared(ms, ms_0, window,ms_2)
			#endif
			/* window to get f1m_0 */
			for (ix=0; ix<ntr; ix++) {
				#pragma ivdep
				for (it=0; it<nt; it++) {
					ms[it+ix*nt] =-ms_2[it+ix*nt] +scale*window[it+ix*nt]*ms_0[it+ix*nt];  
					//f1m_0[it+ix*nt] = scale*f1m_0[it+ix*nt];  
				  //ms[it+ix*nt] = ms_0[it+ix*nt];
				}	
			}
      fft1_2D_fwd (ms,MS,ntr);


	if (verb) fprintf(stderr,"---> Beginning Iteration\n");
		for (iter=0; iter<niter; iter++) {

	/* initialise MS1 and f1m1 the coda for f1+ */	
		memset(MS1,0,2*nf*ntr*sizeof(float));
		memset(F1m1,0,2*nf*ntr*sizeof(float));

			/*------------------------------------------------------------*/
			/* Loop over shot positions */
			/*------------------------------------------------------------*/
			for (ishot=0; ishot<nshots; ishot++) {
		
				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d,e,f) \
					shared(MS,taper,Refl,F1m,F1m1)
				#endif 	
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {
						
						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						/*(a + bi)(e + fi) = (ae - bf) + (af + be)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = MS[ishot*2*nf+it];
						d = MS[ishot*2*nf+it+1];
						e = F1m[ishot*2*nf+it];
						f = F1m[ishot*2*nf+it+1];
						
						F1m1[ix*2*nf+it]   += (a*c - mode*b*d) - r*(a*e - b*f);
						F1m1[ix*2*nf+it+1] += (mode*a*d + b*c) - r*(a*f + b*e);
					
					} /* End of loop over frequencies */	
				} /* End of loop over receivers (traces) */
				
			} /* End of loop over shot positions */

			/* Get time domain output of f1m and ms */
      fft1_2D_inv (F1m1,f1m,ntr);
			
			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it) \
				shared(f1m, f1m_0, window)
			#endif
			for (ix=0; ix<ntr; ix++) {
				#pragma ivdep
				for (it=0; it<nt; it++) {
					f1m[it+ix*nt] = f1m_0[it+ix*nt] + scale*window[it+ix*nt]*(f1m[it+ix*nt]);  
				}	
			}
			
      fft1_2D_fwd (f1m,F1m,ntr);


			for (ishot=0; ishot<nshots; ishot++) {
		
				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d,e,f) \
					shared(MS,MS1,taper,Refl,F1m)
				#endif 	
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {
						
						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						/*(a + bi)(e + fi) = (ae - bf) + (af + be)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = MS[ishot*2*nf+it];
						d = MS[ishot*2*nf+it+1];
						e = F1m[ishot*2*nf+it];
						f = F1m[ishot*2*nf+it+1];

						MS1[ix*2*nf+it]    += (a*e - mode*b*f) - r*(a*c - b*d);
						MS1[ix*2*nf+it+1]  += (mode*a*f + b*e) - r*(a*d + b*c);
					
					} /* End of loop over frequencies */	
				} /* End of loop over receivers (traces) */
				
			} /* End of loop over shot positions */

			/* Get time domain output of f1m and ms */
      fft1_2D_inv (MS1,ms,ntr);
			
			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it) \
				shared(ms, window)
			#endif
			for (ix=0; ix<ntr; ix++) {
				#pragma ivdep
				for (it=0; it<nt; it++) {
					ms[it+ix*nt] =-ms_2[it+ix*nt]+ scale*window[it+ix*nt]*(ms[it+ix*nt]);  
				}	
			}
      fft1_2D_fwd (ms,MS,ntr);

			if(iter%4==0) fprintf(stderr,"Iteration %d\n",iter);
      convergence[iter] = 0.0f;
			for (ix=0; ix<ntr; ix++) {
				for (it=0; it<nt; it++) {
					convergence[iter] += ms[it+ix*nt]*ms[it+ix*nt];
				}	
			}      


		} /* End of loop over iterations */ 
    

		/* Build f1p* by adding Tinv to coda M */
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(f1pS,F_arrival,ms)
		#endif
		for (ix=0; ix<ntr; ix++) {
			#pragma ivdep
			for (it=0; it<nt; it++) {
				f1pS[it+ix*nt] =  F_arrival[it+ix*nt] + ms[it+ix*nt];
				/* note  this is the time reverse version of f1p */
			}	
		}
    fft1_2D_fwd (f1pS,F1pS,ntr);

	/* to get G by looping over shots */
		memset(Gp,0,2*nf*ntr*sizeof(float));
		memset(Gm,0,2*nf*ntr*sizeof(float));
			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d,e,f) \
					shared(F1pS, F1m, taper, Refl, Gp, Gm)
				#endif
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = F1pS[ishot*2*nf+it];
						d = F1pS[ishot*2*nf+it+1];
						e = F1m[ishot*2*nf+it];
						f = F1m[ishot*2*nf+it+1];

						Gm[ix*2*nf+it]   += (a*c -mode* b*d) -r*(a*e - b*f);
						Gm[ix*2*nf+it+1] += (mode*a*d + b*c) -r*(a*f + b*e);

						Gp[ix*2*nf+it]   += -(a*e - mode*b*f) + r*(a*c - b*d);
						Gp[ix*2*nf+it+1] += -(mode*a*f + b*e) + r*(a*d + b*c);   
						

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      fft1_2D_inv (Gp,gp1,ntr);
      fft1_2D_inv (Gm,gm1,ntr);

		if (Pf1) { 
			if (verb) fprintf(stderr,"---> Build f1p\n");
			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,c,d) \
					shared(F1pS, F1p)
				#endif
					#pragma ivdep
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						c = F1pS[ishot*2*nf+it];
						d = F1pS[ishot*2*nf+it+1];

						F1p[ishot*2*nf+it]   =  c;
						F1p[ishot*2*nf+it+1] = -d;   
						

					} /* End of loop over frequencies */
			}
		}
      fft1_2D_inv (F1p,f1p,ntr);

	 


			if (verb) fprintf(stderr,"Build Gp, Gm and G\n");

		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(f1m, f1pS, gp, gm, gp1, gm1, G)
		#endif
		for (ix=0; ix<ntr; ix++) {
			#pragma ivdep
			for (it=0; it<nt2; it++) {
				gm[it+ix*nt2] =  ((scale*gm1[it+nt2+ix*nt]) - f1m[ it+nt2+ix*nt])*(1.0 - window[it+nt2+ix*nt]); 	
				gp[it+ix*nt2] =  ((scale*gp1[it+nt2+ix*nt]) + f1pS[it+nt2+ix*nt])*(1.0 - window[it+nt2+ix*nt]); 	
				G[ it+ix*nt2] = 0.5*(gp[it+ix*nt2] +    gm[it+ix*nt2]);
			}
		}

		fprintf(stderr,"Build Gp, Gm and G\n");



		/* Write the final result */
	    /*FRefl = sf_input(argv[1]);*/
		/*fft1(Gp,,FRefl,1,0,0);
		fft1(Gm,,FRefl,1,0,0);*/
		sf_floatwrite(convergence,niter,Fconvergence);
    free(convergence);

		sf_floatwrite(gp,nt2*ntr,FGp);
		sf_fileclose(FGp);
		sf_floatwrite(gm,nt2*ntr,FGm);
		sf_fileclose(FGm);
		if (PG) {
			sf_floatwrite(G,nt2*ntr,FG);
			sf_fileclose(FG);
		}
		if (Pf1) {
			sf_floatwrite(f1p,nt*ntr,Ff1p);
			sf_fileclose(Ff1p);
			sf_floatwrite(f1m,nt*ntr,Ff1m);
			sf_fileclose(Ff1m);
		}
		
		if (twin) {
			sf_floatwrite(window,nt*ntr,Ftwin);
		sf_fileclose(Ftwin);
	}
	sf_fileclose(FRefl);
	sf_fileclose(FF_arrival);
	
	free(G);
	free(Gm);
	free(Gp);
	free(Refl);
	free(f1pS);
	free(F1pS);
	free(f1p);
	free(F1p);
	free(tw);
	free(filename1);
	free(MS);
	free(MS_0);
	free(MS1);
	free(MS2);
	free(F1m);
	free(F1m1);
	free(f1m_0);
	free(F1m_0);
	free(f1m);
	free(gp);
	free(gp1);
	free(gm);
	free(gm1);
	free(ms);
	free(ms_0);
	free(ms_2);
	free(F_arrival);
	free(window);
	
    exit (0);
}

