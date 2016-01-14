/*
t
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
#include "fft1.h"




typedef struct _Mar *Marchenko;
/*^*/

struct _Mar{
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
	    sf_axis at,af,ax,at1;
	    int     nt,nt2,nf,ntr,mode,nshots,niter,len;
	    int     i,it,ix,ishot,iter,i0;
	    int		twc, twa, shift, n[2], rect[2], s[2], tap;
	    float   scale,eps,dt,df,dx,ot,of,a,b,c,d,e,f,r;
    float *convergence;
		sf_triangle tr;
		bool verb,conj,twin,Pf1,PG;
    /* independent arrays */
    float  *F_arrival, *ms, *ms_0, *ms_2, *f1m_0, *f1m, *f1pS, *f1p; // time
    float  *MS, *MS_0, *MS1, *MS2, *F1m_0, *F1m, *F1m1, *F1pS, *F1p; // frequency
    /* output gf */
    float *G, *gp1, *gp, *gm1, *gm, *Gp, *Gm;
    float *window;
};
/*^*/


Marchenko  init_Marchenko(
     bool verb, bool conj,bool twin, bool Pf1, bool PG, int niter, int nshots, float r, float scale, float eps, int shift, int tap)
/*< >*/
{
    Marchenko mar = (Marchenko) sf_alloc(1,sizeof(*mar));
    mar->verb = verb;
    mar->conj = conj;
    mar->twin = twin;
    mar->Pf1  = Pf1;
    mar->PG   = PG;
    mar->niter = niter;
    mar->nshots = nshots;
    mar->r = r;
    mar->scale = scale;
    mar->eps = eps;
    mar->shift = shift;
    mar->tap = tap;
    sf_warning("exiting init");
    return mar;
}


void set_files(Marchenko mar)
/*< >*/
{
    sf_warning("mar->niter = %d",mar->niter);
      mar->Fconvergence = sf_output("conv");
      mar->convergence = sf_floatalloc(mar->niter);
    sf_warning("mar->niter = %d",mar->niter);
      sf_putint (mar->Fconvergence, "n1", mar->niter);
      sf_putint (mar->Fconvergence, "n2", 1);
      sf_putfloat (mar->Fconvergence, "d1", 1.0f);
      sf_putfloat (mar->Fconvergence, "o1", 0.0f);


	    /*------------------------------------------------------------*/
	    /* I/O files 												  */
	    /*------------------------------------------------------------*/	
		/* "in" is the transposed version of p00plus_xxxx_xxxx.rsf
		   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=ntr,n2=nt
		   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=nt,n2=ntr */
		mar->FF_arrival = sf_input("in");
		
		/* refl is REFL_000.rsf
		   It is used to read nf, df, of
		  Dimensions are: n1=nf,n2=ntr */
		/*FRefl = (sf_file)sf_alloc(1,sizeof(sf_file));*/	
		mar->FRefl = sf_input("refl");

		mar->FGp = sf_output("out");
		mar->FGm = sf_output("Gm");
		
		if (mar->PG) {
			mar->FG  = sf_output("G");
		}
		if (mar->Pf1) {
			mar->Ff1p  = sf_output("f1p");
			mar->Ff1m  = sf_output("f1m");
		}
		
		if (mar->twin){
      sf_warning("twin!");
			mar->Ftwin = sf_output("window"); /* time window */
		}
}




void set_axes_copy(Marchenko mar,Marchenko mar1)
/*< >*/
{

		/*------------------------------------------------------------*/
		/* Axes */
		/*------------------------------------------------------------*/    
		mar->at = mar1->at;
    sf_setlabel(mar->at,"Time"); if(mar->verb) sf_raxa(mar->at); /* time */
    
		mar->at1 = mar1->at1; 
    sf_setlabel(mar->at,"Time"); if(mar->verb) sf_raxa(mar->at); /* time */

		mar->af = mar1->af;
    sf_setlabel(mar->af,"Frequency"); if(mar->verb) sf_raxa(mar->af); /* frequency */

		mar->ax = mar1->ax;
    sf_setlabel(mar->ax,"r"); if(mar->verb) sf_raxa(mar->ax); /* space */
	    
		mar->nt = sf_n(mar->at); mar->dt = sf_d(mar->at); mar->ot = sf_o(mar->at);
		mar->nf = sf_n(mar->af); mar->df = sf_d(mar->af); mar->of = sf_o(mar->af);
		mar->ntr = sf_n(mar->ax); mar->dx = sf_d(mar->ax);

    mar->nt2=(mar->nt/2);
    sf_setn(mar->at1,mar->nt2);
}

void set_axes(Marchenko mar)
/*< >*/
{

		/*------------------------------------------------------------*/
		/* Axes */
		/*------------------------------------------------------------*/    
		mar->at = sf_iaxa(mar->FF_arrival,1); sf_setlabel(mar->at,"Time"); if(mar->verb) sf_raxa(mar->at); /* time */
		mar->at1 = sf_iaxa(mar->FF_arrival,1); sf_setlabel(mar->at,"Time"); if(mar->verb) sf_raxa(mar->at); /* time */
		mar->af = sf_iaxa(mar->FRefl,1); sf_setlabel(mar->af,"Frequency"); if(mar->verb) sf_raxa(mar->af); /* frequency */
		mar->ax = sf_iaxa(mar->FF_arrival,2); sf_setlabel(mar->ax,"r"); if(mar->verb) sf_raxa(mar->ax); /* space */
	    
		mar->nt = sf_n(mar->at); mar->dt = sf_d(mar->at); mar->ot = sf_o(mar->at);
		mar->nf = sf_n(mar->af); mar->df = sf_d(mar->af); mar->of = sf_o(mar->af);
		mar->ntr = sf_n(mar->ax); mar->dx = sf_d(mar->ax);

    mar->nt2=(mar->nt/2);
    sf_setn(mar->at1,mar->nt2);
}


void set_output_files(Marchenko mar)
/*< >*/
{

	    /*------------------------------------------------------------*/
	    /* Setup output data and wavefield header					  */
	    /*------------------------------------------------------------*/
		sf_oaxa(mar->FGp,mar->at1,1);
		sf_oaxa(mar->FGp,mar->ax,2);
		sf_oaxa(mar->FGm,mar->at1,1);
		sf_oaxa(mar->FGm,mar->ax,2);
		
		if (mar->PG) {
			sf_oaxa(mar->FG,mar->at1,1);
			sf_oaxa(mar->FG,mar->ax,2);
		}
		if (mar->Pf1) {
			sf_oaxa(mar->Ff1p,mar->at,1);
			sf_oaxa(mar->Ff1p,mar->ax,2);
			sf_oaxa(mar->Ff1m,mar->at,1);
			sf_oaxa(mar->Ff1m,mar->ax,2);
		}
		if (mar->twin) {
			sf_oaxa(mar->Ftwin,mar->at,1);
			sf_oaxa(mar->Ftwin,mar->ax,2);
		}
}


void allocate_Arrays(Marchenko mar)
/*< >*/
{
		mar->ms   = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_0 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_2 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m_0= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1pS = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1p  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));

    mar->F_arrival = (float *)calloc(mar->nt*mar->ntr,sizeof(float));

		sf_floatread(mar->F_arrival,mar->nt*mar->ntr,mar->FF_arrival);


		memcpy(mar->ms,mar->F_arrival,mar->nt*mar->ntr*sizeof(float));

		/* Allocate for coda M of f2 - Frequency */
		mar->MS   = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS_0 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS1  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS2  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m_0= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m1 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1pS = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1p  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));

		/* Output wavefields */
		mar->G = (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gp1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gp= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gm1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gm= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->Gp= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->Gm= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
    
    mar->window =  sf_floatalloc(mar->nt*mar->ntr);
}


void just_allocate_Arrays(Marchenko mar)
/*< >*/
{
    sf_warning("ntr=%d nt=%d",mar->ntr,mar->nt);
		mar->ms   = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_0 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_2 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m_0= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1pS = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1p  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));

    mar->F_arrival = (float *)calloc(mar->nt*mar->ntr,sizeof(float));

		/* Allocate for coda M of f2 - Frequency */
		mar->MS   = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS_0 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS1  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS2  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m_0= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m1 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1pS = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1p  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));

		/* Output wavefields */
		mar->G = (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gp1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gp= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gm1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gm= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->Gp= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->Gm= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
    
    mar->window =  sf_floatalloc(mar->nt*mar->ntr);
}


void setFirstArrival( Marchenko mar, float *input)
/*< >*/
{
    for (int it=0; it<mar->nt*mar->ntr ; ++it)
      mar->F_arrival[it] = input[it];
		memcpy(mar->ms,mar->F_arrival,mar->nt*mar->ntr*sizeof(float));
}


void allocate_Arrays_F(Marchenko mar, float *F_arrival)
/*< >*/
{
		mar->ms   = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_0 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->ms_2 = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m_0= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1m  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1pS = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->f1p  = (float *)calloc(mar->nt*mar->ntr,sizeof(float));

    mar->F_arrival = (float *)calloc(mar->nt*mar->ntr,sizeof(float));
    for (int it=0; it<mar->nt ; ++it)
      mar->F_arrival[it] = F_arrival[it];

		memcpy(mar->ms,mar->F_arrival,mar->nt*mar->ntr*sizeof(float));

		/* Allocate for coda M of f2 - Frequency */
		mar->MS   = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS_0 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS1  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->MS2  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m_0= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1m1 = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1pS = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->F1p  = (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));

		/* Output wavefields */
		mar->G = (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gp1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gp= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->gm1= (float *)calloc(mar->nt*mar->ntr,sizeof(float));
		mar->gm= (float *)calloc(mar->nt2*mar->ntr,sizeof(float));
		mar->Gp= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
		mar->Gm= (float *)calloc(2*mar->nf*mar->ntr,sizeof(float));
    
    mar->window =  sf_floatalloc(mar->nt*mar->ntr);
}

void zeroArrays(Marchenko mar)
/*< >*/
{
	  memset(mar->ms,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->ms_0,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->ms_2,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->f1m_0,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->f1m,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->f1pS,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->f1p,0,mar->nt*mar->ntr*sizeof(float));

	  memset(mar->MS,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->MS_0,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->MS1,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->MS2,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->F1m_0,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->F1m,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->F1m1,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->F1pS,0,mar->nf*mar->ntr*sizeof(float));
	  memset(mar->F1p,0,mar->nf*mar->ntr*sizeof(float));
  
	  memset(mar->G,0,mar->nt2*mar->ntr*sizeof(float));
	  memset(mar->gp1,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->gp,0,mar->nt2*mar->ntr*sizeof(float));
	  memset(mar->gm1,0,mar->nt*mar->ntr*sizeof(float));
	  memset(mar->gm,0,mar->nt2*mar->ntr*sizeof(float));
	  memset(mar->Gp,0,2*mar->nf*mar->ntr*sizeof(float));
	  memset(mar->Gm,0,2*mar->nf*mar->ntr*sizeof(float));

	  memset(mar->window,0,mar->nt*mar->ntr*sizeof(float));
}



void  buildWindow(Marchenko mar)
/*< >*/
{

    int ix,it,iter;
    int ntr = mar->ntr;
    int nt = mar->nt;
    int shift = mar->shift;
    int n[2],s[2],rect[2];
    int i0;
    float eps = mar->eps;
    sf_triangle tr; 


		/* Build time-window */
		int *tw = sf_intalloc(ntr);
		/*memset(window,0,nt*ntr*sizeof(float));*/
	    /* I am not sure why I set it to this value */
		/*for (ix=0; ix<ntr; ix++) {
			tw[ix] = nt*dt+ot+0.15; 
		}*/
		
		if (mar->verb) fprintf(stderr,"---> Build time-window?\n");
    // checking time sample corresponding to muting time
		for (ix=0; ix<ntr; ix++) {
      for (it=0; it<nt; it++) {
        if ((mar->F_arrival[it+ix*nt]*mar->F_arrival[it+ix*nt])>eps*eps) {
          /*tw[ix] = it*dt+ot;*/
          tw[ix] = it;
          break;
        }
			}
		}
		if (mar->verb) fprintf(stderr,"---> Build time-window1\n");
		for (ix=0; ix<ntr; ix++) {
			int twc = (int)(tw[ix]-shift-10);
			int twa = (int)(-twc+nt);
			/*if (verb) fprintf(stderr,"%d %d\n",twc,twa);*/
			for (it=0; it<nt; it++) {
			/*	if ((it>twa) || (it<twc)) {*/
				if ((it>twa) && (it<twc)) { 
					mar->window[it+ix*nt] = 1.0; // building windowing function W from Filippo's paper
				}
			}
		}

		if (mar->verb) fprintf(stderr,"---> Build time-window2\n");
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
					sf_smooth2 (tr,i0,s[ix],false,mar->window );
				}
			}
			sf_triangle_close(tr);
		}
    free(tw);
}




float* buildTaper(Marchenko mar)
/*< >*/
{


		/* Tapering */
		float pi = 4.0*atan(1.0);
		int ntr = mar->ntr;
    int ix;
    int tap = mar->tap;

		float *taper = (float *)calloc(ntr,sizeof(float));
		memset(taper,0,ntr*sizeof(float));

		for (ix=0; ix<mar->tap; ix++) {
			taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/(2*tap))));
			taper[ntr-ix-1] = taper[ix];
		}
		for (ix=tap; ix<(ntr-tap); ix++) {
			taper[ix] = 1.0;
		}
		if (mar->verb) fprintf(stderr,"---> taper finish\n");
    return taper;
}




void initSolutions(Marchenko mar,float *Refl,float *taper)
/*< >*/
{
   float a,b,c,d;
   int ix,it,ishot;

   int nshots = mar->nshots;
   int ntr = mar->ntr;
   int nf = mar->nf;
   int nt = mar->nt;
   float r = mar->r;

  int mode = mar->mode;
  float scale= mar->scale;

    fft1_init (mar->nt, mar->dt, mar->ot, true, false);
    fft1_2D_fwd(mar->F_arrival,mar->MS,mar->ntr);


		/*------------------------------------------------------------*/
		/* Loop over iterations */
		/*------------------------------------------------------------*/
		if (mar->verb) fprintf(stderr,"---> Begin to iterative solve for f1p and f1m\n");
		/*starting iteration for f1m */
		memset(mar->F1m_0,0,2*nf*ntr*sizeof(float));

			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = mar->MS[ishot*2*nf+it];
						d = mar->MS[ishot*2*nf+it+1];

						mar->F1m_0[ix*2*nf+it]   += (a*c - mode*b*d);
						mar->F1m_0[ix*2*nf+it+1] += (mode*a*d + b*c);
						mar->MS2  [ix*2*nf+it]   += r*(a*c - b*d); // rTd* R
						mar->MS2  [ix*2*nf+it+1] += r*(a*d + b*c);

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      // Inverse fft
      fft1_2D_inv (mar->F1m_0, mar->f1m_0,ntr);
      fft1_2D_inv (mar->MS2, mar->ms_2,ntr);


			/* window to get f1m_0 */
			for (ix=0; ix<ntr; ix++) {
				for (it=0; it<nt; it++) {
					mar->f1m_0[it+ix*nt] = mar->scale*mar->window[it+ix*nt]*mar->f1m_0[it+ix*nt];  
					mar->ms_2[it+ix*nt] = mar->scale*mar->window[it+ix*nt]*mar->ms_2[it+ix*nt];  
					//f1m_0[it+ix*nt] = scale*f1m_0[it+ix*nt];  
					mar->f1m[it+ix*nt] = mar->f1m_0[it+ix*nt];
				}	
			}
			
      fft1_2D_fwd (mar->f1m, mar->F1m,ntr);


	/* initialise MS the coda for f1+ */
		memset(mar->MS_0,0,2*nf*ntr*sizeof(float));
		memset(mar->MS,0,2*nf*ntr*sizeof(float));

			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = mar->F1m[ishot*2*nf+it];
						d = mar->F1m[ishot*2*nf+it+1];

						mar->MS_0[ix*2*nf+it]   += (a*c - mode*b*d);
						mar->MS_0[ix*2*nf+it+1] += (mode*a*d + b*c);

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      fft1_2D_inv (mar->MS_0,mar->ms_0,ntr);


			/* window to get f1m_0 */
			for (ix=0; ix<ntr; ix++) {
				for (it=0; it<nt; it++) {
					mar->ms[it+ix*nt] =-mar->ms_2[it+ix*nt]
                             +mar->scale*mar->window[it+ix*nt]*mar->ms_0[it+ix*nt];  
				}	
			}
      fft1_2D_fwd (mar->ms,mar->MS,ntr);
}





void iteration(Marchenko mar,float *Refl,float *taper)
/*< >*/
{
   float a,b,c,d,e,f;
   int ix,it,ishot;

   int nshots = mar->nshots;
   int ntr = mar->ntr;
   int nf = mar->nf;
   int nt = mar->nt;
   float r = mar->r;

  int mode = mar->mode;
  float scale= mar->scale;

    fft1_init (mar->nt, mar->dt, mar->ot, true, false);



	/* initialise MS1 and f1m1 the coda for f1+ */	
		memset(mar->MS1,0,2*nf*ntr*sizeof(float));
		memset(mar->F1m1,0,2*nf*ntr*sizeof(float));

			/*------------------------------------------------------------*/
			/* Loop over shot positions */
			/*------------------------------------------------------------*/
			for (ishot=0; ishot<nshots; ishot++) {
		
				/* Loop over receivers (traces) */
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					for (it=0; it<2*nf; it=it+2) {
						
						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						/*(a + bi)(e + fi) = (ae - bf) + (af + be)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = mar->MS[ishot*2*nf+it];
						d = mar->MS[ishot*2*nf+it+1];
						e = mar->F1m[ishot*2*nf+it];
						f = mar->F1m[ishot*2*nf+it+1];
						
						mar->F1m1[ix*2*nf+it]   += (a*c - mode*b*d) - r*(a*e - b*f);
						mar->F1m1[ix*2*nf+it+1] += (mode*a*d + b*c) - r*(a*f + b*e);
					
					} /* End of loop over frequencies */	
				} /* End of loop over receivers (traces) */
				
			} /* End of loop over shot positions */

			/* Get time domain output of f1m and ms */
      fft1_2D_inv (mar->F1m1,mar->f1m,ntr);
			
			for (ix=0; ix<ntr; ix++) {
				for (it=0; it<nt; it++) {
					mar->f1m[it+ix*nt] = mar->f1m_0[it+ix*nt] +
                               scale*mar->window[it+ix*nt]*(mar->f1m[it+ix*nt]);  
				}	
			}
			
      fft1_2D_fwd (mar->f1m,mar->F1m,ntr);


			for (ishot=0; ishot<nshots; ishot++) {
		
				/* Loop over receivers (traces) */
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					for (it=0; it<2*nf; it=it+2) {
						
						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						/*(a + bi)(e + fi) = (ae - bf) + (af + be)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = mar->MS[ishot*2*nf+it];
						d = mar->MS[ishot*2*nf+it+1];
						e = mar->F1m[ishot*2*nf+it];
						f = mar->F1m[ishot*2*nf+it+1];

						mar->MS1[ix*2*nf+it]    += (a*e - mode*b*f) - r*(a*c - b*d);
						mar->MS1[ix*2*nf+it+1]  += (mode*a*f + b*e) - r*(a*d + b*c);
					
					} /* End of loop over frequencies */	
				} /* End of loop over receivers (traces) */
				
			} /* End of loop over shot positions */

			/* Get time domain output of f1m and ms */
      fft1_2D_inv (mar->MS1,mar->ms,ntr);
			
			for (ix=0; ix<ntr; ix++) {
				for (it=0; it<nt; it++) {
					mar->ms[it+ix*nt] =-mar->ms_2[it+ix*nt]+ scale*mar->window[it+ix*nt]*(mar->ms[it+ix*nt]);  
				}	
			}
      fft1_2D_fwd (mar->ms,mar->MS,ntr);
}


void iterations(Marchenko mar,float *Refl,float *taper)
/*< >*/
{
int iter=0;

if (mar-> verb) fprintf(stderr,"---> Beginning Iteration\n");
  for (iter=0; iter<mar->niter; iter++) {
    sf_warning("iter=%d",iter);
    iteration(mar,Refl,taper);
  }
}





void buildF1p(Marchenko mar)
/*< >*/
{
  int ix, it;
   int nshots = mar->nshots;
   int ntr = mar->ntr;
   int nf = mar->nf;
   int nt = mar->nt;

  fft1_init (mar->nt, mar->dt, mar->ot, true, false);
	/* Build f1p* by adding Tinv to coda M */


		/* Build f1p* by adding Tinv to coda M */
		for (ix=0; ix<ntr; ix++) {
			for (it=0; it<nt; it++) {
				mar->f1pS[it+ix*nt] =  mar->F_arrival[it+ix*nt] + mar->ms[it+ix*nt];
				/* note  this is the time reverse version of f1p */
			}	
		}
    fft1_2D_fwd (mar->f1pS,mar->F1pS,ntr);
}



void buildGpGmG(Marchenko mar, float *Refl, float *taper)
/*< >*/
{
   float a,b,c,d,e,f;
   int ix,it,ishot;

   int nshots = mar->nshots;
   int ntr = mar->ntr;
   int nf = mar->nf;
   int nt = mar->nt;
   float r = mar->r;
   int nt2 = mar->nt2;
  int mode = mar->mode;
  float scale= mar->scale;
  bool verb = mar->verb;

  fft1_init (mar->nt, mar->dt, mar->ot, true, false);

	/* to get G by looping over shots */
		memset(mar->Gp,0,2*nf*ntr*sizeof(float));
		memset(mar->Gm,0,2*nf*ntr*sizeof(float));
			for (ishot=0; ishot<nshots; ishot++) {

				/* Loop over receivers (traces) */
				for (ix=0; ix<ntr; ix++) {
					/* Loop over frequencies */
					for (it=0; it<2*nf; it=it+2) {

						/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
						a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
						b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
						c = mar->F1pS[ishot*2*nf+it];
						d = mar->F1pS[ishot*2*nf+it+1];
						e = mar->F1m[ishot*2*nf+it];
						f = mar->F1m[ishot*2*nf+it+1];

						mar->Gm[ix*2*nf+it]   += (a*c -mode* b*d) -r*(a*e - b*f);
						mar->Gm[ix*2*nf+it+1] += (mode*a*d + b*c) -r*(a*f + b*e);

						mar->Gp[ix*2*nf+it]   += -(a*e - mode*b*f) + r*(a*c - b*d);
						mar->Gp[ix*2*nf+it+1] += -(mode*a*f + b*e) + r*(a*d + b*c);   
						

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */
			}
      fft1_2D_inv (mar->Gp,mar->gp1,ntr);
      fft1_2D_inv (mar->Gm,mar->gm1,ntr);

		if (mar->Pf1) { 
			if (verb) fprintf(stderr,"---> Build f1p\n");
			for (ishot=0; ishot<nshots; ishot++) {
				/* Loop over receivers (traces) */
				for (it=0; it<2*nf; it=it+2) {
					/*(a + bi)(c + di) = (ac - bd) + (ad + bc)i*/
					c = mar->F1pS[ishot*2*nf+it];
					d = mar->F1pS[ishot*2*nf+it+1];

					mar->F1p[ishot*2*nf+it]   =  c;
					mar->F1p[ishot*2*nf+it+1] = -d;   

				} /* End of loop over frequencies */
			}
		}
    fft1_2D_inv (mar->F1p,mar->f1p,ntr);

		if (verb) fprintf(stderr,"Build Gp, Gm and G\n");

		for (ix=0; ix<ntr; ix++) {
			for (it=0; it<nt2; it++) {
				mar->gm[it+ix*nt2] =  ((scale*mar->gm1[it+nt2+ix*nt]) - mar->f1m[ it+nt2+ix*nt])*(1.0 - mar->window[it+nt2+ix*nt]); 	
				mar->gp[it+ix*nt2] =  ((scale*mar->gp1[it+nt2+ix*nt]) + mar->f1pS[it+nt2+ix*nt])*(1.0 - mar->window[it+nt2+ix*nt]); 	
				mar->G[ it+ix*nt2] = 0.5*(mar->gp[it+ix*nt2] +    mar->gm[it+ix*nt2]);
			}
		}

}

void marchenko_close(Marchenko mar)
/*< >*/
{
		free(mar->ms);
		free(mar->ms_0);
		free(mar->ms_2);
		free(mar->f1m_0);
		free(mar->f1m);
		free(mar->f1pS);
		free(mar->f1p);

    free(mar->F_arrival);


		/* Allocate for coda M of f2 - Frequency */
		free(mar->MS);
		free(mar->MS_0);
		free(mar->MS1);
		free(mar->MS2);
		free(mar->F1m_0);
		free(mar->F1m);
		free(mar->F1m1);
		free(mar->F1pS);
		free(mar->F1p);

		/* Output wavefields */
		free(mar->G);
		free(mar->gp1);
		free(mar->gp);
		free(mar->gm1);
		free(mar->gm);
		free(mar->Gp);
		free(mar->Gm);
    free(mar->window);
}


void close_files(Marchenko mar)
/*< >*/
{
		sf_fileclose(mar->FF_arrival);
		
		sf_fileclose(mar->FRefl);

		sf_fileclose(mar->FGp);
		sf_fileclose(mar->FGm);
		
		if (mar->PG) {
			sf_fileclose(mar->FG);
		}
		if (mar->Pf1) {
			sf_fileclose(mar->Ff1p);
			sf_fileclose(mar->Ff1m);
		}
		
		if (mar->twin){
			sf_fileclose(mar->Ftwin);
		}
}


