
% Given an analytic function f, initial root z0, initial box radius
% r0 = width + height * i
% Return a root near z0

function [z, numRoots] = ArgumentPrincipleRootFinder(f, z, r0, tol)

	% Cardinal directions in complex plane
	w = [1; 0; -1; 0] * abs(real(r0)) + [0; 1i; 0; -1i] * abs(imag(r0));

	% Useful for computation
	K = 1:4;
	Kp1 = mod(K, 4) + 1;

	% Compute max # of iterations
	maxIt = log2( abs(r0) / tol);
	numRoots = 1;

	% Bisection iteration
	for i = 1:maxIt

		% Check if tolerance is good enough
		if ( abs(f(z)) < tol)
			break
		end

		% STEP: 1 Compute all of the line integrals

		% Line integrals that go through our guess z
		Zint = zeros(4, 1);
		for j = 1:4
			Zint(j) = LineIntegrate(z, z + w(j));
		end

		% Line integrals to go around the "outer perimeter"
		Wint = zeros(4, 2);
		for j = 1:4
			a = w(j);
			b = w(mod(j, 4) + 1);
			Wint(j, 1) = LineIntegrate(z + a, z + a + b);
			Wint(j, 2) = LineIntegrate(z + a, z + a - b);
		end

		% STEP: 2 Use line integrals to compute argument of each box
		BoxInts = real( -1i * ( Zint(K) - Zint(Kp1) + Wint(K,1) - Wint(Kp1,2) ) / (2 * pi) );

		% Get box with max number of roots
		[~, kMax] = max(BoxInts);
		numRoots = round(sum(BoxInts));
		if (i == 1 && numRoots == 0)
			return
		end

		kMaxp1 = mod(kMax,4) + 1;

		% Update z to be at center of new box
		z = z + (w(kMax) + w(kMaxp1)) / 2;

		% Update radiuses to be smaller
		w = w / 2;
	end

	return;

	% Function to call for line integration of f'/f from a to b
	% Using trapezoidal integration method
	function I = LineIntegrate(a, b)
		N = max(8, ceil(10 * abs(b - a)));
		dz = (b - a) / N;

		zVals = linspace(a - dz, b + dz, N + 3);
		fVals = f(zVals);
		Kint = 3:(N + 1);
		Bint = [2, N + 2];

		I = sum( (fVals(Kint + 1) - fVals(Kint - 1)) ./ (2 * fVals(Kint)) ) + ...
			sum( (fVals(Bint + 1) - fVals(Bint - 1)) ./ (2 * fVals(Bint)) )/2;

	end

end
