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



typedef struct Mar *Marchenko;
/*^*/

struct Mar{
    int nt,nf,ntr,mode,nshots,niter;
    int   nz,nzpad;
    int   nx,nxpad;
    float eps,dt,scale,df,dx,ot,of;


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

	    int     nt,nf,ntr,mode,nshots,niter,len;
	    int     i,it,ix,ishot,iter,i0;
	    int		twc, twa, shift, n[2], rect[2], s[2], tap;
	    float   scale,eps,dt,df,dx,ot,of,a,b,c,d,e,f,r;


		/* Cube axes */
	    sf_axis at,af,ax,at1;

	    int     nt,nf,ntr,mode,nshots,niter,len;
	    int     i,it,ix,ishot,iter,i0;
	    int		twc, twa, shift, n[2], rect[2], s[2], tap;
	    float   scale,eps,dt,df,dx,ot,of,a,b,c,d,e,f,r;

		sf_triangle tr;

};
/*^*/





