% Version: $Id: ch1d.m,v 1.6 2006/03/22 21:33:56 bolo Exp $
% CH1D  Calculate Cahn-Hilliard equation in 1D.
%
% STATE = CH1D(S) calculates C-H equation for 1 unit of time.
% The initial state S should be initialized by INITCH1D.
% 
% STATE = CH1D(S, T) calculates C-H equation for T units
% of time.
%
% STATE = CH1D(S, T, VF) calculates C-H equation for T units
% of time and plots the results eqery VF units of time.
%
% See also INITCH1D.

function state = ch1d(s, T, vf)
	if (nargin < 2) T = 1.0; end
	if (nargin < 3) vf = inf; end

	x = s{1}; N = size(x); N = N(2);
	u = s{2};
	eps = s{3}(1);
	dx = s{3}(2);
	dt = s{3}(3);
	a = s{3}(4);
	t0 = s{4};

	lambda1 = dt/dx^2;
	lambda2 = eps*lambda1/dx^2;

	L = 2*cos(pi*(0:(N-1))/(N-1)) - 2;
	C = 1 + (1-a)*lambda1*L + lambda2*L.^2;
	P = lambda1*L;

	nextv = t0 + vf;
	uhat = dct(u);
	for t = t0:dt:t0+T
		uold = u; uhatold = uhat;
		uhat = (uhat + P.*dct2(u.^3 - a*u)) ./ C;
		u = idct(uhat);

		if (t >= nextv)
			e = eps*0.5*(u(2:N)-u(1:N-1)).^2/dx^2 + 0.25*(u(1:N-1).^2-1).^2;
			g = (u(2:N)-u(1:N-1))/dx;
			l = idct(L.*uhat)/dx^2;
			lold = idct(L.*uhatold/dx^2);
%			mf = -eps*(l(2:N) - l(1:N-1))/dx ...
%				+ (3*(uold(2:N)+uold(1:N-1)).^2/4 - a).*(uold(2:N)-uold(1:N-1))/dx ...
%				- (1-a)*(u(2:N)-u(1:N-1))/dx;
%			mf = -eps*(l(2:N) - l(1:N-1))/dx ...
%				+ (3*(u(2:N)+u(1:N-1)).^2/4 - 1).*(u(2:N)-u(1:N-1))/dx;
%			mf = -eps*(l(2:N) - l(1:N-1))/dx ...
%				+ (3/4*(uold(2:N)+uold(2:N)).^2 - a).*(uold(2:N)-uold(1:N-1))/dx ...
%				- (1-a)*(u(2:N)-u(1:N-1))/dx;

% 			mf = -eps*(l(3:N) - l(1:N-2))/2/dx ...
% 				+ (3*uold(2:N-1).^2 - a).*(uold(3:N)-uold(1:N-2))/2/dx ...
% 				- (1-a)*(u(3:N)-u(1:N-2))/2/dx;

 			mf = -eps*(l(2:N) - l(1:N-1))/dx ...
 				+ (uold(2:N).^3-uold(1:N-1).^3)/dx ...
 				- a*(uold(2:N)-uold(1:N-1))/dx ...
 				- (1-a)*(u(2:N)-u(1:N-1))/dx;

			stabW = idct(-eps/dx^2*L.*uhat  + dct(uold.^3 - a*uold - (1-a)*u));
			unstabW = idct(-eps/dx^2*L.*uhat  + dct(u.^3 - u));

			stabRHS = idct(-lambda2*L.^2.*uhat  + lambda1*L.*dct(uold.^3 - a*uold) - (1-a)*lambda1*L.*uhat);
			% stabRHS = (u - uold);
			unstabRHS = idct(-lambda2*L.^2.*uhat  + lambda1*L.*dct(u.^3 - u));

%			mf = cumsum(stabRHS);

			subplot(3,3,1)
			plot(x, u, x, sqrt(1/3), x, -sqrt(1/3)), axis([x(1) x(N) -1 1]), title(['u  (t=' num2str(t) ')']);
			subplot(3,3,4)
			% plot(x, stabRHS), title('u_t - u_{t-1}');
			plot(x, stabRHS), title('\Delta(-\epsilon \Delta u + u^3 - u)  (stable)');
			subplot(3,3,7)
			plot(x, stabW), title('-\epsilon \Delta u + u^3 - u  (stable)');

			subplot(3,3,5)
			plot(x, unstabRHS), title('\Delta(-\epsilon \Delta u + u^3 - u)  (unstable)');
			subplot(3,3,8)
			plot(x, unstabW), title('\Delta(-\epsilon \Delta u + u^3 - u)  (unstable)');

			subplot(3,3,9)
%			plot(e), title('\epsilon/2 |\nabla u|^2 + 1/4 (u^2 - 1)^2');

			subplot(3,3,9)
			plot(x(1:N-1), mf), title('-\epsilon (\Delta u)'' + (3 u^2 - 1) u'' (wrong?)');

			subplot(3,3,6)
			bar(1:N, uhat(1:N)-uhatold(1:N)), title('uhat_t - uhat_{t-1}');
			subplot(3,3,2)
			plot(x(1:N-1), g), title('\nabla u');
			subplot(3,3,3)
			bar(uhat), title('uhat');

			nextv = t + vf;
			pause(0.1);
		end
	end

	state = {s{1} u s{3} t {L C P}};

