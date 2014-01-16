% Version: $Id: net.m,v 1.1 2006/07/11 20:37:56 bolo Exp $
% NET  Calculate heat equation on 1-D net.
%
% STATE = NET(S) calculates heat equation for 1 unit of time.
% The initial state S should be initialized by INITNET.
% 
% STATE = NET(S, T) calculates heat equation for T units
% of time.
%
% STATE = NET(S, T, VF) calculates heat equation for T units
% of time and plots the results eqery VF units of time.
%
% See also INITNET

function state = heat(s, T, vf)

if (nargin < 2) T = 1.0; end
if (nargin < 3) vf = inf; end

M = s{1}{1};
S = s{1}{2};

u = s{2};

f = s{3};

lambda = s{4}(1);
dx = s{4}(2);
dt = s{4}(3);
N = s{4}(4);

c1 = s{5}(1);
c2 = s{5}(2);
c3 = s{5}(3);
c4 = s{5}(4);

x = 0:dx:(N-1)*dx;

u = u';
f = f';

nextV = vf;

clf;

for t = 0:dt:T
	uold = u;
	A = (M + dt*lambda*S);
	b = M * (uold + f);
	u = inv(A) * b;
	if (t >= nextV)
		minH = min(u);
		maxH = max(u);
		subplot(2,4,[1 2])
		plot(x(1:c1), u(1:c1)), axis([x(1) x(c1) minH maxH]), title(['u  (t=' num2str(t) ')']);

		subplot(2,4,3)
		plot(x(c1+1:c2), u(c1+1:c2)), axis([x(c1+1) x(c2) minH maxH]), title(['u  (t=' num2str(t) ')']);

		subplot(2,4,4)
		plot(x(c2+1:c3), u(c2+1:c3)), axis([x(c2+1) x(c3) minH maxH]), title(['u  (t=' num2str(t) ')']);

		subplot(2,4,7)
		plot(x(c3+1:c4), u(c3+1:c4)), axis([x(c3+1) x(c4) minH maxH]), title(['u  (t=' num2str(t) ')']);

		subplot(2,4,8)
		plot(x(c4+1:N), u(c4+1:N)), axis([x(c4+1) x(N) minH maxH]), title(['u  (t=' num2str(t) ')']);

		nextV = (floor(t/vf)+1)*vf;
		pause(0.1);
	end

end

u = u';

state = {s{1} u s{3} s{4} s{5}};

