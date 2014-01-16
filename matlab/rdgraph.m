function graph = initgraph(filenamein, h, filenameout)

if (nargin < 2) h = 0.01; end
if (nargin < 3) filenameout = [filenamein '.graph']; end

file = dlmread(filenamein, ' ');

nnodes = file(1, 1);
nodes = file(2:nnodes+1,1:3);
nedges = file(nnodes+2, 1);
edges = file(nnodes+3:nnodes+nedges+2,1:2)+1;

file(nnodes+2,4) = 0;


nodes(1:nnodes,4) = 0;

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
end

f = zeros(N,1);
f(1:nnodes) = file(2:nnodes+1,4);

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
		M(p,p) = 2/3 * h;
		S(p,p) = 2 / h;
		if (p ~= ps)
			M(p,p-1) = 1/6 * h;
			S(p,p-1) = -1/h;
		end
		if (p ~= pe)
			M(p,p+1) = 1/6 * h;
			S(p,p+1) = -1/h;
		end
	end

	if (ps <= pe)
		M(i1, ps) = 1/6 * h;
		M(ps, i1) = 1/6 * h;
		S(i1, ps) = -1/h;
		S(ps, i1) = -1/h;

		M(i2, pe) = 1/6 * h;
		M(pe, i2) = 1/6 * h;
		S(i2, pe) = -1/h;
		S(pe, i2) = -1/h;
	else
		M(i1, i2) = 1/6 * h;
		M(i2, i1) = 1/6 * h;
		S(i1, i2) = -1/h;
		S(i2, i1) = -1/h;
	end

	M(i1,i1) = M(i1,i1) + h/3;
	S(i1,i1) = S(i1,i1) + 1/h;

	M(i2,i2) = M(i2,i2) + h/3;
	S(i2,i2) = S(i2,i2) + 1/h;
end

[r,c] = find(M);
nconn = 0;
conn = [];
for i = [r,c]'
	if (i(1) < i(2))
		nconn = nconn+1;
		conn(nconn, 1:2) = i-1;
	end
end

graphname = filenameout;
dlmwrite(graphname, [N], ' ');
dlmwrite(graphname, points, '-append', 'delimiter', ' ');
dlmwrite(graphname, [nconn], '-append', 'delimiter', ' ');
dlmwrite(graphname, conn, '-append', 'delimiter', ' ');

graph = {N, M, S, f, nodes, edges};

