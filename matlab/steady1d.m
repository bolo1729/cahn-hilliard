function f = steady1d(a, l, N, eps, dx)

% 	x(3) = a;
% 	x(4) = a + dx^2*l/2;
% 	x(2) = x(4);
% 
% 	x(1) = 0;
% 	n = 3;
% 	x(n+2) = -(x(n-2) - 4*x(n-1) + 6*x(n) -4*x(n+1)) + dx^2/eps*(x(n-1)^3 -2*x(n)^3 + x(n+1)^3 - x(n-1) + 2*x(n) - x(n+1));
% 	x(5) = x(5)/2;
% 	x(1) = x(5);
% 
% 	for n = 4:N
% 	x(n+2) = -(x(n-2) - 4*x(n-1) + 6*x(n) -4*x(n+1)) + dx^2/eps*(x(n-1)^3 -2*x(n)^3 + x(n+1)^3 - x(n-1) + 2*x(n) - x(n+1));
% 	end
% 
% 	f = x;


	u = 0.3*ones(1, N);

	for it=1:10

	l1 = 1/dx^2;
	l2 = eps^2/dx^4;

	A = -l2*ones(1,N);
	B = 4*l2 + l1*u.^2 - l1;
	C = -6*l2 - 2*l1*u.^2 + 2*l1;
	D = 4*l2 + l1*u.^2 - l1;
	E = -l2*ones(1,N);

	E(3) = 2*E(3);

	M = spdiags([A(3:N) ; B(3:N) ; C(3:N) ; D(3:N) ; E(3:N)]', [-4 -3 -2 -1 0], speye(N-2,N-2));

	b = a + l/2/l1;

	m = zeros(1,N-2);
	m(1) = l2*(6*a - 8*b) + l1*(-2*b^3 + 2*a^3) + l1*(2*b - 2*a);
	m(2) = l2*(7*b - 4*a) + l1*(-a^3 + 2*b^3) + l1*(a - 2*b);
	m(3) = l2*(a - 4*b) - l1*b^3 + l1*b;
	m(4) = l2*b;

	u = [a b (inv(M)*m')'];

	end

	f = u;

