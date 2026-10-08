% ----------------------
% Assign default parameters before assignCommandLineOptions to par
% To account for different cases

function par = assignDefaults( temp_par, par )

	% Assign default options based on type of Known Solution
	if (isfield(temp_par, 'knownSolution'))
		switch temp_par.knownSolution		
			case 'GravityCapillaryWave'
				par.xa=0;
				par.xb=1;
				par.ya=-1;
				par.yb=0;

				par.bcs = 'ppnt';
				par.map = 'freeSurface';
				par.motion = 'freeSurfaceMotion';

				par.icfs = 'cos';
				par.ampfs= 1e-8;
				par.gravity = -1;

			case 'Poiseuille'
				par.bcs = 'Ionn';
				par.outflowPressureCoeff=1;
				par.outflowPressureCoeffpn=0;

		end
	end

	% Finally "merge" temp_par with par
	names = fieldnames(temp_par);
	for n = 1:length(names)
		name = names{n};
		par = setfield(par, name, getfield(temp_par, name));
	end

end
