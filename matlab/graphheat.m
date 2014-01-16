% Version: $Id: graphheat.m,v 1.1 2006/07/13 12:13:35 bolo Exp $
% GRAPHHEAT  Calculate heat equation on 1-D net.
%
% STATE = GRAPHHEAT(S) calculates heat equation for 1 unit of time.
% The initial state S should be initialized by INITGRAPH.
% 
% STATE = GRAPHHEAT(S, T) calculates heat equation for T units
% of time.
%
% STATE = GRAPHHEAT(S, T, VF) calculates heat equation for T units
% of time and plots the results eqery VF units of time.
%
% See also INITGRAPH

function state = graphheat(s, T, vf)

if (nargin < 2) T = 1.0; end
if (nargin < 3) vf = inf; end

name = s{6};

name = [name '.out'];

M = s{1}{1};
S = s{1}{2};

u = s{2};

f = s{3};


dx = s{4}(1);
dt = s{4}(2);
N = s{4}(3);

lambda = s{5}(1);

x = 0:dx:(N-1)*dx;

nextV = vf;

iV = round(T/vf);

dlmwrite(name, [N], ' ');
dlmwrite(name, [iV], '-append', 'delimiter', ' ');

iT = round(T/dt);
for i = 0:iT
	t = i*dt;
	uold = u;
	A = (M + dt*lambda*S);
	b = M * (uold + f);
	u = inv(A) * b;
	if (t >= nextV)
		t
		dlmwrite(name, u, '-append', 'delimiter', ' ');
		nextV = (floor(t/vf)+1)*vf;
	end

end

state = {s{1} u s{3} s{4} s{5} s{6}};

