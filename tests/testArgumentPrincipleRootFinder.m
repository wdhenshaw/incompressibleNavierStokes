%
% Tests of ArgumentPrincipleRootFinder as a refiner of the initial guess omega0 for the secant
% iteration in getViscousCapillaryWaveCoeffs.
%
%   addpath(genpath(pwd)); testArgumentPrincipleRootFinder
%

% Sanity check: after n iterations the root must still be inside the current box, |z-z*| <= sqrt(2) r0/2^n.
% tol=0 so that exactly n iterations are taken.
fprintf('\nSanity check: f = z - z*, z0=0, r0=1, n=1..8 iterations\n');
zStarValues = [ 0.3137+0.2171i; -0.3137+0.2171i; -0.3137-0.2171i; 0.3137-0.2171i ];
for q=1:4
  zStar = zStarValues(q);
  f = @(z) z - zStar;
  inside = 1;
  for n=1:8
    z = ArgumentPrincipleRootFinder( f, 0, 1, 0, n );
    inside = inside && abs(z - zStar) <= sqrt(2)/2^n;
  end
  fprintf('  quadrant %d, z*=%7.4f%+7.4fi: root stayed in the box: %s\n',q,real(zStar),imag(zStar),passFail(inside));
end

% Several roots: the wanted root is the least damped one (largest imaginary part), but here two
% more damped roots share a quadrant
zRoots = [ 0.4+0.6i; -0.3-0.4i; -0.6-0.7i ];
f = @(z) (z - zRoots(1)).*(z - zRoots(2)).*(z - zRoots(3));
z = ArgumentPrincipleRootFinder( f, 0, 1, 0, 8 );
[~, nearest] = min( abs(z - zRoots) );
fprintf('\nSeveral roots: z*=%g%+gi (least damped), %g%+gi, %g%+gi\n',real(zRoots(1)),imag(zRoots(1)),...
        real(zRoots(2)),imag(zRoots(2)),real(zRoots(3)),imag(zRoots(3)));
fprintf('  after 8 iterations z=%7.4f%+7.4fi is nearest root %d, least damped root found: %s\n',real(z),imag(z),nearest,passFail(nearest==1));

% Root on an integration path: f=0 at a sample point gives f'/f = Inf
fprintf('\nRoot on an integration path (diagnostic), z0=0, r0=1, 4 iterations\n');
for zStar = [ 0.5, 0.5+0.3i ]
  f = @(z) z - zStar;
  z = ArgumentPrincipleRootFinder( f, 0, 1, 0, 4 );
  fprintf('  z*=%g%+gi: z=%g%+gi, |z-z*|=%9.2e\n',real(zStar),imag(zStar),real(z),imag(z),abs(z-zStar));
end

% Main test: does refining omega0 let the secant converge where the raw guess fails?
par = struct('kx',2*pi,'xa',0,'xb',1,'ya',-1,'yb',0,'nu',0,'rho',1,'gravity',-1,'gamma',0.05,'ampfs',1e-4,'omega0',0,'idebug',0);
k     = par.kx;
H     = par.yb - par.ya;
rho   = par.rho;
g     = -par.gravity;
gamma = par.gamma;
nuValues = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
numItValues = [2, 4, 6];

fprintf('\nSecant from raw and refined omega0 (gamma=%g, k=2*pi, H=1, rho=1, g=1), r0 = 0.5|omega0|\n',gamma);
fprintf('ok = secant converges to the reference omega (relative difference < 1e-8)\n');
fprintf('%5s %-24s %-24s %4s','nu','reference omega','omega0','raw');
for numIt = numItValues
  fprintf('   n=%d: |z-ref|  ',numIt);
end
fprintf('\n');

for nu = nuValues
  par.nu = nu;

  % Reference from the findRootBox search
  par.omega0 = 0;
  try
    [~, vcw] = getViscousCapillaryWaveCoeffs( par );
    omegaRef = vcw.omega;
  catch
    omegaRef = NaN;
  end

  omega0 = sqrt( (g*k + gamma*k^3/rho)*tanh(k*H) ) - 1i*2*nu*k^2;
  r0 = 0.5*abs(omega0);
  F = @(omega) arrayfun( @(w) viscousCapillaryCountFunction( w, k, H, nu, rho, g, gamma ), omega );

  fprintf('%5.2f %10.4f%+10.4fi   %10.4f%+10.4fi   %4s',nu,real(omegaRef),imag(omegaRef),...
          real(omega0),imag(omega0),okFail(secantFrom(par,omega0,omegaRef)));
  for numIt = numItValues
    z = ArgumentPrincipleRootFinder( F, omega0, r0, 0, numIt );
    fprintf('   %9.2e %4s',abs(z - omegaRef),okFail(secantFrom(par,z,omegaRef)));
  end
  fprintf('\n');
end
fprintf('Cost: 12 line integrals of N+3=7 points, about 84 evaluations of F per iteration\n\n');


% Run the secant in getViscousCapillaryWaveCoeffs from omega0 and check it reaches omegaRef
function ok = secantFrom( par, omega0, omegaRef )
  par.omega0 = omega0;
  try
    [~, vcw] = getViscousCapillaryWaveCoeffs( par );
    ok = abs(vcw.omega - omegaRef) < 1e-8*abs(omegaRef);
  catch
    ok = 0;
  end
end


function s = passFail( pass )
  if( pass )
    s = 'PASS';
  else
    s = 'FAIL';
  end
end


function s = okFail( ok )
  if( ok )
    s = 'ok';
  else
    s = 'fail';
  end
end


% Copy of getCountFunction in getViscousCapillaryWaveCoeffs: F = omega^2 det(M unscaled), scaled by exp(-Re(alpha) H)
function F = viscousCapillaryCountFunction( omega, k, H, nu, rho, g, gamma )
  mu = nu*rho;
  G  = rho*g + gamma*k^2;

  alpha   = sqrt( k^2 - 1i*omega/nu );
  scaleAH = exp( -real(alpha)*H );
  phaseAH = exp( 1i*imag(alpha)*H );
  decayAH = exp( -2*real(alpha)*H )/phaseAH;
  coshAHs = 0.5*( phaseAH + decayAH );
  sinhAHs = 0.5*( phaseAH - decayAH );
  coshKHs = cosh(k*H)*scaleAH;
  sinhKHs = sinh(k*H)*scaleAH;

  vHatA   = (k/alpha)*sinhAHs - sinhKHs;
  vHatB   = coshAHs - coshKHs;
  vHatyA  = k*( coshAHs - coshKHs );
  vHatyB  = alpha*sinhAHs - k*sinhKHs;
  vHatyyA = k*( alpha*sinhAHs - k*sinhKHs );
  vHatyyB = alpha^2*coshAHs - k^2*coshKHs;

  F11 = vHatyyA + k^2*vHatA;
  F12 = vHatyyB + k^2*vHatB;
  F21 = -omega^2*coshKHs + (1i*k*omega/rho)*2*mu*vHatyA - (k*G/rho)*vHatA;
  F22 = -omega^2*sinhKHs + (1i*k*omega/rho)*2*mu*vHatyB - (k*G/rho)*vHatB;

  F = F11*F22 - F12*F21;
end
