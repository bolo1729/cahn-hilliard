% Version: $Id: netwave.m,v 1.1 2006/07/14 16:29:14 bolo Exp $
% NETWAVE  Calculate heat equation on 1-D net.
%
% STATE = NETWAVE(S) calculates heat equation for 1 unit of time.
% The initial state S should be initialized by INITGRAPH.
% 
% STATE = NETWAVE(S, T) calculates heat equation for T units
% of time.
%
% STATE = NETWAVE(S, T, VF) calculates heat equation for T units
% of time and plots the results eqery VF units of time.
%
% See also INITNET

function state = netheat(s, T, vf)

if (nargin < 2) T = 1.0; end
if (nargin < 3) vf = inf; end

M = s{1}{1};
S = s{1}{2};

N = s{2}{1};
horig = s{2}{2};
dt = s{2}{3};
nameout = s{2}{4};

t0 = s{3};
u = s{4}{1};
w = s{4}{2};

f = s{5};

lambda = s{6}(1);



nextV = 0;

iV = round(T/vf);

dlmwrite(nameout, [iV+1], '-append', 'delimiter', ' ');

AA = speye(2*N);
AA(1:N, 1:N) = M;
AA(1:N, N+1:2*N) = -M*dt;
AA(N+1:2*N, 1:N) = S*dt;
AA(N+1:2*N, N+1:2*N) = M;

uu = [u' w']';

BB = speye(2*N);
BB(1:N,1:N) = M;
BB(N+1:2*N, N+1:2*N) = M;

ff = zeros(2*N, 1); ff(1:N) = f;

iT = round(T/dt);
for i = 0:iT
	t = t0 + 1.0*i*dt;
	uuold = uu;
	bb = BB * (uuold + dt*dt*ff);
	[uu, flag] = bicg(AA, bb);
	if (t >= t0 + nextV)
		t
		u = uu(1:N);
		dlmwrite(nameout, u, '-append', 'delimiter', ' ');
		nextV = (floor(t/vf)+1)*vf;
	end

end

state = {s{1} s{2} t0+T {u, w} s{5} s{6}};

