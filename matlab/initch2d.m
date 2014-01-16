% $Id: initch2d.m,v 1.5 2006/07/11 20:37:56 bolo Exp $

function ret = initch2d(u, eps, dx, dt, a, t)
	if (nargin < 1)
		M = 512; N = 4;
		M = 128; N = 128;
		x = (0:M-1)/(M-1); y = (0:N-1)/(N-1);
		u = 0.01*sin(3*pi*y'*x) + 0.04*cos(4*pi*ones(N,1)*x) + 0.06*sin(5*pi*ones(N,1)*y) + 0.01*cos(10*pi*y'*x);
	end

	MN = size(u); M = MN(2); N = MN(1);

	if (nargin < 2) eps = 0.0001; end
	if (nargin < 3) dx = 1.0/(M-1); end
	if (nargin < 4) dt = 0.0001; end
	if (nargin < 5) a = 3; end
	if (nargin < 6) t = 0.0; end

	x = 0:dx:(M-1)*dx;
	y = 0:dx:(N-1)*dx;

	ret = {{x y} u [eps dx dt a] t};

