function state = initnet(nametopo, consts, nameresults, rdev, h, dt)

if (nargin < 2) consts = []; end
if (nargin < 3) nameresults = nametopo; end
if (nargin < 4) rdev = 0; end
if (nargin < 5) h = 0.01; end
if (nargin < 6) dt = 0.001; end

namein = [nametopo '.in'];
namegraph = [nameresults '.graph'];
nameout = [nameresults '.out'];

file = dlmread(namein, ' ');

nnodes = file(1, 1);
nodes = file(2:nnodes+1, 1:3);
nedges = file(nnodes+2, 1);
edges = file(nnodes+3:nnodes+nedges+2, 1:2) + 1;

% resize file to 6 columns if necessary
s = size(file); if (s(2) < 6) file(1, 6) = 0; end;

% forces applied to nodes
nf = file(2:nnodes+1, 4);
% initial values at nodes
nu = file(2:nnodes+1, 5);
% first derivatives of values at nodes
nuu = file(2:nnodes+1, 6);

% a parameter defined on edges
eparam = file(nnodes+3:nnodes+nedges+2, 3) + 1.0;
nparam = zeros(nnodes,1);

% store node degrees in 4th column
nodes(1:nnodes, 4) = 0;

% first sweep through edges in order to calculate the total amount of points (N) in the output graph
N = nnodes;
for e = 1:nedges
	i1 = edges(e,1);
	i2 = edges(e,2);
	nodes(i1,4) = nodes(i1,4) + 1;
	nodes(i2,4) = nodes(i2,4) + 1;
	x1 = nodes(i1,1:3);
	x2 = nodes(i2,1:3);
	d = sqrt(sum((x1 - x2).^2));
	nint = ceil(d/h);
	hint = d/nint;
	edges(e,3) = N + 1;
	N = N + nint - 1;
	edges(e,4) = N;
	edges(e,5) = hint;

	nparam(i1) = nparam(i1) + eparam(e);
	nparam(i2) = nparam(i2) + eparam(e);
end

horig = h;

f = zeros(N,1);
f(1:nnodes) = nf(1:nnodes);

u = zeros(N,1);
u(1:nnodes) = nu(1:nnodes);

uu = zeros(N,1);
uu(1:nnodes) = nuu(1:nnodes);

M = spdiags(zeros(N,1), [0], speye(N));
S = spdiags(zeros(N,1), [0], speye(N));

points = nodes(1:nnodes, 1:3);

for e = 1:nedges
	i1 = edges(e,1);
	i2 = edges(e,2);
	x1 = nodes(i1,1:3);
	x2 = nodes(i2,1:3);
	d = sqrt(sum((x1 - x2).^2));
	ps = edges(e,3);
	pe = edges(e,4);
	h = edges(e,5);

	for p = ps:pe
		points(p,1:3) = x1 + (x2-x1)*(p-ps+1)/(pe-ps+2);
		u(p) = u(i1) + (u(i2)-u(i1))*(p-ps+1)/(pe-ps+2);
		uu(p) = uu(i1) + (uu(i2)-uu(i1))*(p-ps+1)/(pe-ps+2);
		f(p) = f(i1) + (f(i2)-f(i1))*(p-ps+1)/(pe-ps+2);

		M(p,p) = 2/3 * h;
		S(p,p) = 2 / h * eparam(e);
		if (p ~= ps)
			M(p,p-1) = 1/6 * h;
			S(p,p-1) = -1/h * eparam(e);
		end
		if (p ~= pe)
			M(p,p+1) = 1/6 * h;
			S(p,p+1) = -1/h * eparam(e);
		end
	end

	if (ps <= pe)
		M(i1, ps) = 1/6 * h;
		M(ps, i1) = 1/6 * h;
		S(i1, ps) = -1/h * eparam(e);
		S(ps, i1) = -1/h * eparam(e);

		M(i2, pe) = 1/6 * h;
		M(pe, i2) = 1/6 * h;
		S(i2, pe) = -1/h * eparam(e);
		S(pe, i2) = -1/h * eparam(e);
	else
		M(i1, i2) = 1/6 * h;
		M(i2, i1) = 1/6 * h;
		S(i1, i2) = -1/h * eparam(e);
		S(i2, i1) = -1/h * eparam(e);
	end

	M(i1,i1) = M(i1,i1) + h/3;
	S(i1,i1) = S(i1,i1) + 1/h * eparam(e);

	M(i2,i2) = M(i2,i2) + h/3;
	S(i2,i2) = S(i2,i2) + 1/h * eparam(e);
end



% find connections between points
[r,c] = find(M);
nconn = 0;
conn = [];
for i = [r,c]'
	if (i(1) < i(2))
		nconn = nconn+1;
		conn(nconn, 1:2) = i-1;
	end
end

% store points and connections in graph file
dlmwrite(namegraph, [N], ' ');
dlmwrite(namegraph, points, '-append', 'delimiter', ' ');
dlmwrite(namegraph, [nconn], '-append', 'delimiter', ' ');
dlmwrite(namegraph, conn, '-append', 'delimiter', ' ');

% initialize output file
dlmwrite(nameout, [N], ' ');

if (rdev > 0)
	u = u + random('unif', 0, rdev, N, 1);
	uu = uu + random('unif', 0, rdev, N, 1);
end;

% unew = u + uu*dt;

T = 0;

state = {{M, S}, {N, horig, dt, nameout}, T, {u, uu}, f, consts};

