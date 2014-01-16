% Version: $Id: netheat.m,v 1.1 2006/07/14 16:29:14 bolo Exp $
% NETHEAT  Calculate heat equation on 1-D net.
%
% STATE = NETHEAT(S) calculates heat equation for 1 unit of time.
% The initial state S should be initialized by INITGRAPH.
% 
% STATE = NETHEAT(S, T) calculates heat equation for T units
% of time.
%
% STATE = NETHEAT(S, T, VF) calculates heat equation for T units
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

f = s{5};

lambda = s{6}(1);



nextV = 0;

iV = round(T/vf);

dlmwrite(nameout, [iV+1], '-append', 'delimiter', ' ');

iT = round(T/dt);
for i = 0:iT
	t = t0 + 1.0*i*dt;
	uold = u;
	A = (M + dt*lambda*S);
	b = M * (uold + dt*f);
	[u, flag] = bicg(A, b);
	if (t >= t0 + nextV)
		t
		dlmwrite(nameout, u, '-append', 'delimiter', ' ');
		nextV = (floor(t/vf)+1)*vf;
	end

end

state = {s{1} s{2} t0+T {u} s{5} s{6}};

