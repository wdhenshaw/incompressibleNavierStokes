
function z = findComplexRoot(f, z, APchecks, tol, maxIt)

	if (APchecks)
		r0 = 0.95 * z;

		tolAP = 1e-2*abs(r0);

		for n = 1:APchecks
			
			[z, numRoots] = ArgumentPrincipleRootFinder(f, z, r0, tolAP);

			if (numRoots) 
				break; 
			end
			z = z + r0;
			r0 = 2 * r0;

		end
		if (~numRoots)
			error('Did not find root with Argument Princple!');
		end
	end

	z = SecantMethodRootFinder(f, z, tol, maxIt);
	return;

end
