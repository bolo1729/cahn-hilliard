% Version: $Id: initch1d.m,v 1.6 2006/07/11 20:37:56 bolo Exp $
% INITCH1D  Initial state for Cahn-Hilliard solver.
%
% STATE = INITCH1D generates a function and sets default
% parameters for solver.
%
% STATE = INITCH1D(U) sets the initial function to U.
%
% STATE = INITCH1D(U, EPS) sets the initial function to U,
% and epsilon to EPS.
%
% STATE = INITCH1D(U, EPS, DX, DT, A, T0) sets the initial
% function to U, epsilon to EPS, and space and time distances
% to DX and DT, stabilization parameter to A, and initial
% time to T0.
%
% See also CH1D.

function state = initch1d(u, eps, dt, dx, a, t0)
	if (nargin < 1)
		N = 128; x = (0:N-1)/(N-1);
		u = 0.1*sin(2*pi*x) + 0.01*cos(4*pi*x) + 0.06*sin(4*pi*x) + 0.02*cos(10*pi*x);
	end

	N = size(u); N = N(2);

	if (nargin < 2) eps = 0.0005; end
	if (nargin < 3) dx = 1.0/(N-1); end
	if (nargin < 4) dt = 0.01; end
	if (nargin < 5) a = 3; end
	if (nargin < 6) t0 = 0.0; end

	x = 0:dx:(N-1)*dx;

	state = {x u [eps dx dt a] t0};

