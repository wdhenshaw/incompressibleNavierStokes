
% Secant method for complex function
function zn = SecantMethodRootFinder(f, znm1, tol, maxIt)


	fnm1 = f(znm1);

	zn = znm1 * (1 + 1e-3) + eps;
	fn = f(zn);

	converged = 0;
	for k = 1:maxIt

		dz = -fn * (zn - znm1) / (fn - fnm1);
		
		znm1 = zn;
		fnm1 = fn;
		
		zn = zn + dz;
		fn = f(zn);

		if( abs(dz) < tol*abs(zn) )
			return
		end
	end

	warning('Secant Iteration did not converge past tolerance!');

end
