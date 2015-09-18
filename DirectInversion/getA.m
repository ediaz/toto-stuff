function [A,S,M] = getA(f,m,h,n,freeSurface)
% 5-point discretization of the 2D Helmholtz operator with absorbing
% boundary conditions.
%
% use:
%   A = getA(f,m,h,n);
%
% input:
%   f - frequency [Hz]
%   m - model [m^2/s^2]
%   h - [dz,dx] gridspacing in z and x direction [m]
%   n - [nz,nx] number of gridpoints in z and x direction
%
% ouput:
%   A - sparse matrix
%
% This program is part of the paper
% "Mitigating local minima in full-waveform inversion by expanding the search space",
% T. van Leeuwen and F.J. Herrmann, 2013 (submitted to GJI).
%
% Copyright (C) 2013 Tristan van Leeuwen (tleeuwen@eos.ubc.ca)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
switch nargin
    case 4
        freeSurface = 0;
end

% angular frequency
omega = 2*pi*f;

% number of gridpoints
N     = prod(n);

% Stiffness matrices
S = LaplacianOperator(h,n);

% absorving 1 grid point layers (w =1 inside domain, 0 at the boundary layer) 

M = Mass(omega,m,n,h,freeSurface );

% 
A = M + S;




