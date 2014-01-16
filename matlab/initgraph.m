% Equation u' = lambda * laplace(u)
function state = initgraph(name, consts, h, nameout)

if (nargin < 3) h = 0.01; end
if (nargin < 4) nameout = name; end

dt = 1/1000;

g = rdgraph([name '.in'], h, [nameout '.graph']);

N = g{1};

M = g{2};
S = g{3};

% initial function
u = random('unif', 0, 0.1, N, 1);

% source term
f = g{4};

state = {{M, S}, u, f, [h, dt, N], consts, nameout};

