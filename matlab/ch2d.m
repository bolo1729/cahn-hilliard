% $Id: ch2d.m,v 1.3 2006/03/22 21:33:56 bolo Exp $

function ret = ch2d(s, T, vf)
	if (nargin < 2) T = 1.0; end
	if (nargin < 3) vf = inf; end

	x = s{1}{1};
	y = s{1}{2};
	u = s{2};
	eps = s{3}(1);
	dx = s{3}(2);
	dt = s{3}(3);
	a = s{3}(4);
	t0 = s{4};

	MN = size(u); M = MN(2); N = MN(1);

	lambda1 = dt/dx^2;
	lambda2 = eps*lambda1/dx^2;

	L = 2*cos(pi*(0:(N-1))/(N-1))'*ones(1,M)  +  ones(N,1)*2*cos(pi*(0:(M-1))/(M-1))  -  4;
	C = 1 + (1-a)*lambda1*L + lambda2*L.^2;
	P = lambda1*L;

	nextv = t0 + vf;
	uhat = dct2(u);
	for t = t0:dt:t0+T
		uhat = (uhat + P.*dct2(u.^3 - a*u)) ./ C;
		u = idct2(uhat);

		if (t >= nextv)
			l = idct2(L.*uhat/dx^2);

			subplot(3,3,[1 2 4 5])
			pcolor(x, y, u), shading interp, axis('off'); %, axis('equal');
			
			subplot(3,3,[7 8])
			plot(x, u(N/2, :));

			subplot(3,3,3)
			pcolor(uhat), shading interp, axis('off'), axis('equal');

			subplot(3,3,6)
			pcolor(x, y, l), shading interp, axis('off'), axis('equal');

			nextv = t + vf;
			pause(0.1);
		end
	end

	ret = {s{1} u s{3} t};

