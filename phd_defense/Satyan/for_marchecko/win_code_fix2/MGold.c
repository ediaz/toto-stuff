/* Marchenko-Wapenaar-Broggini iterative scheme

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

void fft1 (float *, float *, sf_file, bool, bool, bool);

int main(int argc, char* argv[])
{

	bool verb,conj,twin,Pf2;
	
    /* OMP parameters */
	#ifdef _OPENMP
    int ompnth;
	#endif 	
	
	float	*F_arrival, *Refl, *G, *G1, *Gf, *f2S, *F2S;
	float	*MS, *MS1, *ms, *ms1, *ms0;
	float	*window, *taper, pi;
	int		*tw;

    /* I/O files */
    sf_file FF_arrival;
    sf_file FRefl;
    sf_file FG;
    sf_file Ff2;
    sf_file Ftwin;

	char *filename1, filename2[256], filename3[256];
	
	/* Cube axes */
    sf_axis at,af,ax;

    int     nt,nf,ntr,mode,nshots,niter,len;
    int     i,it,ix,ishot,iter,i0;
    int		twc, twa, shift, n[2], rect[2], s[2];
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
    if(! sf_getbool("Pf2",&Pf2)) Pf2=false; /* Htot=true: returns H=Gp-Gm */
    if(! sf_getint("niter",&niter)) niter=1; /* number of iterations */
    if(! sf_getint("nshots",&nshots)) nshots=1; /* number of shots */
    if(! sf_getint("r",&r)) r=-1; /* reflection coefficient if flux 
					normalised r=-1 */
    if(! sf_getfloat("scale",&scale)) scale=1.0; /* scale factor */
	if(! sf_getfloat("eps",&eps)) eps=1e-4; /* threshold for the timewindow */
	if(! sf_getint("shift",&shift)) shift=5; /* shift in samples for the timewindow */
    
	if (verb) {
		fprintf(stderr,"This program was called with \"%s\".\n",argv[0]);
		/*fprintf(stderr,"Nr: %d Nx: %d Nt:%d\n",nr,nx,nt);*/
		
		if (argc > 1) {
			for (i = 1; i<argc; i++) {
				fprintf(stderr,"argv[%d] = %s\n", i, argv[i]);
			}
		}
		else {
			fprintf(stderr,"The command had no other arguments.\n");
    	}	
	}

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

	FG = sf_output("out");
	
	if (Pf2) {
		Ff2  = sf_output("f2S");
	}
	
	if (twin) {
		Ftwin = sf_output("window"); /* time window */
	}
	
	/*------------------------------------------------------------*/
	/* Axes */
	/*------------------------------------------------------------*/    
	at = sf_iaxa(FF_arrival,1); sf_setlabel(at,"Time"); if(verb) sf_raxa(at); /* time */
	af = sf_iaxa(FRefl,1); sf_setlabel(af,"Frequency"); if(verb) sf_raxa(af); /* frequency */
	ax = sf_iaxa(FF_arrival,2); sf_setlabel(ax,"r"); if(verb) sf_raxa(ax); /* space */
    
	nt = sf_n(at); dt = sf_d(at); ot = sf_o(at);
	nf = sf_n(af); df = sf_d(af); of = sf_o(af);
	ntr = sf_n(ax); dx = sf_d(ax);
	
	if (verb) fprintf(stderr,"nt: %d nf: %d ntr:%d\n",nt,nf,ntr);

	sf_fileclose(FRefl);

    /*------------------------------------------------------------*/
    /* Setup output data and wavefield header					  */
    /*------------------------------------------------------------*/
	sf_oaxa(FG,at,1);
	sf_oaxa(FG,ax,2);
	
	if (Pf2) {
		sf_oaxa(Ff2,at,1);
		sf_oaxa(Ff2,ax,2);
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
	ms = (float *)calloc(nt*ntr,sizeof(float));
	ms1= (float *)calloc(nt*ntr,sizeof(float));
	ms0= (float *)calloc(nt*ntr,sizeof(float));
	memcpy(ms,F_arrival,nt*ntr*sizeof(float));
	memcpy(ms0,F_arrival,nt*ntr*sizeof(float));


	
	/* Allocate for coda M of f2 - Frequency */
	MS = (float *)calloc(2*nf*ntr,sizeof(float));
	MS1= (float *)calloc(2*nf*ntr,sizeof(float));
	/* The three flags of fft1 are: inv, sym, and opt */
	fft1(F_arrival,MS,FF_arrival,0,0,1);
/*	memcpy(FA,2*nf*ntr*sizeof(float));*/

	/* Output wavefields */
	G = (float *)calloc(nt*ntr,sizeof(float));
	G1= (float *)calloc(nt*ntr,sizeof(float));
	Gf = (float *)calloc(2*nf*ntr,sizeof(float));
	
	f2S = (float *)calloc(nt*ntr,sizeof(float)); /*time f2*/
	F2S = (float *)calloc(2*nf*ntr,sizeof(float)); /*frequency f2*/
	
	
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
	filename1 = sf_getstring("refl");
	/* 000.rsf are 7 characters */
	len = strlen(filename1)-7;
	/* copy the filename without 000.rsf */
	strncpy(filename2,filename1,len);
	filename2[len] = '\0';
	if (verb) fprintf(stderr,"filename2 is: %s and len is: %d\n",filename2,len);
	/*if (NULL == filename1) {
		fprintf(stderr,"Cannot read header file %s",filename1);
	}*/
  
	for (ishot=0; ishot<nshots; ishot++) {
		
	  	/* write xxx.rsf in the string filename3 */
		sprintf(filename3,"%03d.rsf\0",ishot);
		for (i=0; i<7; i++)
			filename2[len+i] = filename3[i];
			filename2[len+7] = '\0';
		if (verb) fprintf(stderr,"Loading %s in memory\n",filename2);
	  	FRefl = sf_input(filename2);
    	sf_floatread(&Refl[ishot*2*nf*ntr],2*nf*ntr,FRefl);
		sf_fileclose(FRefl);
		/*if (verb) fprintf(stderr,"Iteration %d\n",ishot);*/
	}

	/* Build time-window */
	tw = (int *)calloc(ntr,sizeof(int));
	window = (float *)calloc(nt*ntr,sizeof(float));
	/*memset(window,0,nt*ntr*sizeof(float));*/
    /* I am not sure why I set it to this value */
	/*for (ix=0; ix<ntr; ix++) {
		tw[ix] = nt*dt+ot+0.15; 
	}*/
	
	if (verb) fprintf(stderr,"Build the time-window\n");
	for (ix=0; ix<ntr; ix++) {
		for (it=0; it<nt; it++) {
			if (F_arrival[it+ix*nt]>eps) {
				/*tw[ix] = it*dt+ot;*/
				tw[ix] = it;
				/*if (verb) fprintf(stderr,"%d %d\n",ix,it);*/
				break;
			}
		}
	}
	for (ix=0; ix<ntr; ix++) {
		twc = (int)(tw[ix]-shift);
		twa = (int)(-twc+shift+nt);
		/*if (verb) fprintf(stderr,"%d %d\n",twc,twa);*/
		for (it=0; it<nt; it++) {
			if ((it<twc) || (it>twa)) {
				window[it+ix*nt] = 1.0;
			}
		}
	}

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
		tr = sf_triangle_init (rect[ix],n[ix]);
		for (it=0; it < (nt*ntr/n[ix]); it++) {
			i0 = sf_first_index (ix,it,1+1,n,s);
			for (iter=0; iter < 2; iter++) {
				sf_smooth2 (tr,i0,s[ix],false,false,window);
			}
		}
		sf_triangle_close(tr);
	}
	
	/* Tapering */
	pi = 4.0*atan(1.0);
	/*fprintf(stderr,"pi: %f\n",pi);
	fprintf(stderr,"ntr: %d\n",ntr);*/
	
	taper = (float *)calloc(ntr,sizeof(float));
	memset(taper,0,ntr*sizeof(float));

		
	for (ix=0; ix<151; ix++) {
		taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/300)));
		taper[ntr-ix-1] = taper[ix];
	}
	/*for (ix=(ntr-1); ix>(701-151-1); ix--) {
		taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/300)));

	}*/
	for (ix=151; ix<(ntr-151); ix++) {
		taper[ix] = 1.0;
	}
	/*for (ix=0; ix<ntr; ix++) {
		fprintf(stderr,"taper[%d]: %f\n",ix,taper[ix]);
	}*/
	
	FRefl = sf_input("refl");

	/*------------------------------------------------------------*/
	/* Loop over iterations */
	/*------------------------------------------------------------*/
	if (verb) fprintf(stderr,"Beginning of loop over iterations\n");
/*starting iteration for M0* */
 	memset(MS1,0,2*nf*ntr*sizeof(float));

                for (ishot=0; ishot<nshots; ishot++) {

                        /* Loop over receivers (traces) */
                        #ifdef _OPENMP
                        #pragma omp parallel for private(ix,it,a,b,c,d) \
                                shared(MS,taper,Refl,MS1)
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

                                        MS1[ix*2*nf+it]   += (a*c - mode*b*d);
                                        MS1[ix*2*nf+it+1] += (mode*a*d + b*c);

                                } /* End of loop over frequencies */
                        } /* End of loop over receivers (traces) */
                }
		fft1(MS1,ms1,FRefl,1,0,1);
		if (verb) fprintf(stderr,"Build the next iteration of pplus and qplus\n");
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(ms,ms0,window,ms1)
		#endif
		for (ix=0; ix<ntr; ix++) {
			#pragma ivdep
			for (it=0; it<nt; it++) {
				ms0[it+ix*nt] = -scale*window[it+ix*nt]*ms1[it+ix*nt];  
				ms[it+ix*nt] = ms0[it+ix*nt];
			}	
		}
		
		
		if (verb) fprintf(stderr,"%d %d\n",ix,it);
		
		fft1(ms,MS,FF_arrival,0,0,1);

		if(iter%10==0) fprintf(stderr,"Iteration %d\n",iter);



    for (iter=0; iter<niter; iter++) {
	
	/* Set Pminus and Qminus to 0 */
 	memset(MS1,0,2*nf*ntr*sizeof(float));

		/*------------------------------------------------------------*/
		/* Loop over shot positions */
		/*------------------------------------------------------------*/
		if (verb) fprintf(stderr,"Beginning of loop over shot positions\n");
		for (ishot=0; ishot<nshots; ishot++) {
   	
			/* Loop over receivers (traces) */
			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it,a,b,c,d) \
				shared(MS,taper,Refl,MS1)
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

					MS1[ix*2*nf+it]   += (a*c - mode*b*d) + (-1*(a*c - b*d));
					MS1[ix*2*nf+it+1] += (mode*a*d + b*c) + (-1*(a*d + b*c));
				
				} /* End of loop over frequencies */	
			} /* End of loop over receivers (traces) */
			
			if (verb) if(ishot%50==0) fprintf(stderr,"Trace %d\n",ishot);

		} /* End of loop over shot positions */

		/* Save a copy of pplus and qplus before creating their next iteration */
		/* Build the next iteration of Pplus and Qplus */
		fft1(MS1,ms1,FRefl,1,0,1);
		
		
		if (verb) fprintf(stderr,"Build the next iteration of pplus and qplus\n");
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(ms,ms0,window,ms1)
		#endif
		for (ix=0; ix<ntr; ix++) {
			#pragma ivdep
			for (it=0; it<nt; it++) {
				ms[it+ix*nt] = ms0[it+ix*nt] - scale*window[it+ix*nt]*ms1[it+ix*nt];  
			}	
		}
		
		
		if (verb) fprintf(stderr,"%d %d\n",ix,it);
		
		fft1(ms,MS,FF_arrival,0,0,1);

		if(iter%10==0) fprintf(stderr,"Iteration %d\n",iter);

	} /* End of loop over iterations */ 


	/* Build Gp and Gm */
	if (verb) fprintf(stderr,"Build Gp and Gm\n");
	#ifdef _OPENMP
	#pragma omp parallel for private(ix,it) \
		shared(f2S,F_arrival,ms)
	#endif
	for (ix=0; ix<ntr; ix++) {
		#pragma ivdep
		for (it=0; it<nt; it++) {
			f2S[it+ix*nt] =  F_arrival[it+ix*nt] + ms[it+ix*nt];
			/* note  ms1 is the time reversal of the first arrival*/
		}	
	}
	memset(F2S,0,2*nf*ntr*sizeof(float));
		fft1(f2S,F2S,FF_arrival,0,0,1);

/* to get G by looping over shots */
	memset(Gf,0,2*nf*ntr*sizeof(float));
                if (verb) fprintf(stderr,"Beginning of loop over shot positions for G\n");
                for (ishot=0; ishot<nshots; ishot++) {

                        /* Loop over receivers (traces) */
                        #ifdef _OPENMP
                        #pragma omp parallel for private(ix,it,a,b,c,d) \
                                shared(F2S,taper,Refl,r,Gf)
                        #endif
                        for (ix=0; ix<ntr; ix++) {
                                /* Loop over frequencies */
                                #pragma ivdep
                                for (it=0; it<2*nf; it=it+2) {

                                        /*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
                                        a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
                                        b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
                                        c = F2S[ishot*2*nf+it];
                                        d = F2S[ishot*2*nf+it+1];

                                        Gf[ix*2*nf+it]  +=  (a*c -mode* b*d) + (-1*(a*c - b*d));
                                        Gf[ix*2*nf+it+1]+=  (mode*a*d + b*c) + (-1*(a*d + b*c));


                                } /* End of loop over frequencies */
                        } /* End of loop over receivers (traces) */
		}
                        if (verb) if(ishot%50==0) fprintf(stderr,"Trace %d\n",ishot);
		fft1(Gf,G1,FRefl,1,0,1);

        #ifdef _OPENMP
        #pragma omp parallel for private(ix,it) \
                shared(f2S,G,G1)
        #endif
        for (ix=0; ix<ntr; ix++) {
                #pragma ivdep
                for (it=0; it<nt; it++) {
                        G[it+ix*nt] =  f2S[it+ix*nt] + scale*G1[it+ix*nt];
                        /* note this is conjugate of f2 */
                }
        }




	/* Write the final result */
    /*FRefl = sf_input(argv[1]);*/
	/*fft1(Gp,,FRefl,1,0,0);
	fft1(Gm,,FRefl,1,0,0);*/
	
	sf_floatwrite(G,nt*ntr,FG);
	if (Pf2) {
		sf_floatwrite(f2S,nt*ntr,Ff2);
		sf_fileclose(Ff2);
	}
	
	if (twin) {
		sf_floatwrite(window,nt*ntr,Ftwin);
		sf_fileclose(Ftwin);
	}
	sf_fileclose(FRefl);
	sf_fileclose(FG);
	sf_fileclose(FF_arrival);
	
	free(G);
	free(Gf);
	free(Refl);
	free(f2S);
	free(tw);
	free(filename1);
	free(MS);
	free(MS1);
	free(G1);
	free(ms1);
	free(ms);
	free(ms0);
	free(F2S);
	free(F_arrival);
	free(window);
	
    exit (0);
}

