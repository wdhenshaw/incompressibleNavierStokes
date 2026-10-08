
function par = defineKnownSolution( par )

knownSolution = par.knownSolution;
nu = par.nu;
kx = par.kx;
ky = par.ky;

if( strcmp(knownSolution,'TaylorGreen') )

  %       printF("--- Taylor Green vortex is an exact solution ---\n"
  %             " u = sin(k x) cos(k y) F(t)\n"
  %             " v =-cos(k x) sin(k y) F(t)\n"
  %             " p = (rho/4)*( cos(2 k x) + cos(2 k y) ) F(t)*F(t)\n"
  %             "   where F(t) = exp( -2 nu k^2 t )\n");

  par.computeErrors=1;

  ampu=1.; ampp=ampu^2/4; 
  ftg = @(t) exp( -2.*nu*kx^2*t);
  par.ue = @(x,y,t)  ampu*sin(kx*x).*cos(kx*y)*ftg(t); 
  par.ve = @(x,y,t) -ampu*cos(kx*x).*sin(kx*y)*ftg(t);
  par.pe = @(x,y,t)  ampp*( cos(2*kx*x) + cos(2.*kx*y) )*ftg(t)*ftg(t);

  ftgt = @(t) (-2.*nu*kx^2)*exp( -2.*nu*kx^2*t);
  par.uet = @(x,y,t)  ampu*sin(kx*x).*cos(kx*y)*ftgt(t); 
  par.vet = @(x,y,t) -ampu*cos(kx*x).*sin(kx*y)*ftgt(t);

  par.pex = @(x,y,t)  ampp*( -2*kx*sin(2.*kx*x) )*ftg(t)*ftg(t);
  par.pey = @(x,y,t)  ampp*( -2*ky*sin(2.*kx*y) )*ftg(t)*ftg(t);

  par.uexy = @(x,y,t) (-ampu*kx^2)*cos(kx*x).*sin(kx*y)*ftg(t); 
  par.vexy = @(x,y,t) ( ampu*kx^2)*sin(kx*x).*cos(kx*y)*ftg(t);
  

  par.ufe = @(x,y,t) 0.;
  par.vfe = @(x,y,t) 0.;
  par.pfe = @(x,y,t) 0.;

  % par.u0 = @(x,y) par.ue(x,y,0);
  % par.v0 = @(x,y) par.ve(x,y,0);

  % pe2  = @(x,y,t) pe(x,y,t); % to avoid matlab function pe

elseif( strcmp(knownSolution,'Poiseuille') )
  % p = ampp*x; 
  % nu*u.yy = p.x 
  par.computeErrors=1;

  ampu =  (par.pressureInflowValue - par.pOutflow) / ( 2 * par.mu );
  ampp = -(par.pressureInflowValue - par.pOutflow);
  par.ue = @(x,y,t)  ampu*(1.-y).*y; 
  par.ve = @(x,y,t)  0.*x;
  par.pe = @(x,y,t)  ampp*x + par.pressureInflowValue; 

  par.uet = @(x,y,t)  0.*x;
  par.vet = @(x,y,t)  0.*x;

  par.pex = @(x,y,t)  ampp;
  par.pey = @(x,y,t)  0*x;

  par.uexy = @(x,y,t) 0.*x;
  par.vexy = @(x,y,t) 0.*x;
  

 par.ufe = @(x,y,t) 0.;
 par.vfe = @(x,y,t) 0.;
 par.pfe = @(x,y,t) 0.;

 %  pe2  = @(x,y,t) pe(x,y,t); % to avoid matlab function pe

 % Capillary wave where 2 * pi * n * (H / L) = eps, 0 < eps << 1
 % WARN: Currently giving linear convergence in the pressure error
elseif( strcmp(knownSolution,'GravityCapillaryWave') )
	
	% Set the amplitude of the wave so that quadratic terms of PDE are at machine percision
	ampu = par.ampfs;
	
	% Define length and height scale of the problem
	H = abs(par.yb - par.ya);
	L = abs(par.xb - par.xa);

	% Store necessary prameters of the problem for ease of readability
	rho = par.rho;
	% NOTE: Main code uses gravity as negative value, our TW solution uses positive g. Oops!
	g = -par.gravity;

	% Store important dimensionless parameters

	% Inverse of the bond number, ratio between surface tension & bouyant forces
	% BoInv = par.gamma / ( ( rho  * g * L^2 ) ) ;

	% Wave number times different spatial scales
	kH = kx * H;

	% Store commonly occuring dimensionless quantity
	Pi = 1 + (kx^2 * par.gamma) / ( ( rho  * g  ) ) ;

	% The dispersion relation for the specified wave number
	omega = sqrt( g * kx * Pi * tanh( kH ) );

	if (~~imag(omega))
		warning(sprintf('Dispersion Relation Has Imaginary Component %g!!', imag(omega)));
	end

	% Store Dispersion Relation in parameter struct
	par.DispersionRelation = omega;

	% Complex wave of the free surface w/ specified wave#
	wave = @(x, t) ampu * exp( 1i * ( kx * x - omega * t) );

	% Partial derivatives of wave function
	wavex  = @(x, t)   1i *    kx * wave(x, t);
	wavet  = @(x, t)  -1i * omega * wave(x, t);

	wavext = @(x, t) kx * omega * wave(x, t);
	wavexx = @(x, t)    -(kx^2) * wave(x, t);

	% Coefficient function in y for fluid potential
	phiHat   = @(y)      -1i * (g / omega) * Pi * (cosh(kx * y + kH) / cosh(kH));
	phiHaty  = @(y) kx * -1i * (g / omega) * Pi * (sinh(kx * y + kH) / cosh(kH));

	phiHatyy = @(y) kx^2 * phiHat(y);

	% Define exact solution & partial derivatives
	par.ue = @(x, y, t)          real( phiHat  (y) .* wavex (x, t) );
	par.ve = @(x, y, t)          real( phiHaty (y) .* wave  (x, t) );

	par.uex = @(x, y, t)         real( phiHat  (y) .* wavexx(x, t) );
	par.uey = @(x, y, t)         real( phiHaty (y) .* wavex (x, t) );

	par.vex = @(x, y, t)         real( phiHaty (y) .* wavex (x, t) );
	par.vey = @(x, y, t)         real( phiHatyy(y) .* wave  (x, t) );

	par.uet = @(x, y, t)         real( phiHat  (y) .* wavext(x, t) );
	par.vet = @(x, y, t)         real( phiHaty (y) .* wavet (x, t) );

	par.pe  = @(x, y, t) - rho * real( phiHat  (y) .* wavet (x, t) + g * y );
	par.pex = @(x, y, t) - rho * real( phiHat  (y) .* wavext(x, t) );
	par.pey = @(x, y, t) - rho * real( phiHaty (y) .* wavet (x, t) + g );

	% No forcing required as solution is exact
	par.ufe = @(x,y,t) 0.;
	par.vfe = @(x,y,t) 0.;
	par.pfe = @(x,y,t) 0.;


elseif( strcmp(knownSolution,'none') )

  if( strcmp(par.ms,'none') )
    par.uet = @(x,y,t)  0.*x;
    par.vet = @(x,y,t)  0.*x;
    par.ufe = @(x,y,t) 0.;
    par.vfe = @(x,y,t) 0.;
    par.pfe = @(x,y,t) 0.;
  end 
else
 fprintf('ERROR: unknown known-solution: knownSolution=[%s]\n',knownSolution);
 pause; pause; 
end 

return
end
