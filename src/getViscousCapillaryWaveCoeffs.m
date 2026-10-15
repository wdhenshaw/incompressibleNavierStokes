%
% Solve the dispersion relation and find the coefficients of the linearized viscous
% gravity-capillary traveling wave (Dan's thesis thesis, Appendix C.2.1).
%
function [par, vcw] = getViscousCapillaryWaveCoeffs( par )

  k   = par.kx;
  H   = abs(par.yb - par.ya);
  nu  = par.nu;
  rho = par.rho;
  mu  = nu*rho;
  g   = -par.gravity;   % main code uses gravity<0, the TW solution uses g>0
  gamma = par.gamma;
  G   = rho*g + gamma*k^2;

  % These do not depend on omega and are computed once
  coshKH = cosh(k*H);
  sinhKH = sinh(k*H);

  % Starting point for the secant iteration: use -omega0 if given, otherwise locate the root with
  % the argument principle (see findRootBox) inside a box placed using the inviscid frequency and
  % the weak viscous decay rate 2 nu k^2
  rootBox = [];
  if( par.omega0~=0 )
    omega0 = par.omega0;
  else
    omegaInviscid = sqrt( (g*k + gamma*k^3/rho)*tanh(k*H) );
    if( ~(omegaInviscid > 0) )
      error('getViscousCapillaryWaveCoeffs: no inviscid frequency to place the search box (g=%g, gamma=%g). Set -omega0.',g,gamma);
    end
    imTop  = -1e-6*omegaInviscid;    % keep the box off the real axis
    boxTol = 1e-3*omegaInviscid;

    searchBox = [0.25*omegaInviscid, 2*omegaInviscid, -4*nu*k^2, imTop];
    [omega0, rootBox, searchInfo] = findRootBox( @getCountFunction, searchBox, boxTol, par.idebug );
    if( isempty(omega0) )
      % Enlarge the box once
      searchBox = [0.125*omegaInviscid, 4*omegaInviscid, -8*nu*k^2, imTop];
      [omega0, rootBox, searchInfo] = findRootBox( @getCountFunction, searchBox, boxTol, par.idebug );
    end
    if( isempty(omega0) )
      error('getViscousCapillaryWaveCoeffs: no root found in Re(omega) in [%g,%g], Im(omega) in [%g,%g]. Set -omega0.',searchBox(1),searchBox(2),searchBox(3),searchBox(4));
    end
  end

  % Secant iteration for det M(omega) = 0 in the complex plane
  tol = 1e-13;
  maxIt = 100;
  w0 = omega0;
  w1 = omega0*(1 + 1e-3);
  f0 = det( getMatrix(w0) );
  f1 = det( getMatrix(w1) );
  converged = 0;
  for it=1:maxIt
    dw = -f1*(w1-w0)/(f1-f0);
    w0 = w1;
    f0 = f1;
    w1 = w1 + dw;
    f1 = det( getMatrix(w1) );
    if( abs(dw) < tol*abs(w1) )
      converged = 1;
      break;
    end
  end
  if( ~converged )
    error('getViscousCapillaryWaveCoeffs: secant iteration did not converge from omega0=%g%+gi. Set -omega0 to a better guess.',real(omega0),imag(omega0));
  end
  omega = w1;

  % Reject roots that left the root box, are not finite, grow in time, or where det M is not
  % small compared to its terms
  if( ~isempty(rootBox) )
    boxWidth = max( rootBox(2)-rootBox(1), rootBox(4)-rootBox(3) );
    if( abs(omega - omega0) > 2*boxWidth )
      error('getViscousCapillaryWaveCoeffs: secant iteration left the root box, omega0=%g%+gi, omega=%g%+gi.',real(omega0),imag(omega0),real(omega),imag(omega));
    end
  end
  M = getMatrix(omega);
  relDet = abs(det(M))/( abs(M(1,1)*M(2,2)) + abs(M(1,2)*M(2,1)) );
  if( ~isfinite(omega) || imag(omega) > 0 || ~(relDet < 1e-8) )
    error('getViscousCapillaryWaveCoeffs: omega=%g%+gi is not an acceptable root (relative det=%8.2e). Set -omega0 to a better guess.',real(omega),imag(omega),relDet);
  end

  % Null vector [a;b] of M, taken from the row with the larger norm
  if( norm(M(1,:)) >= norm(M(2,:)) )
    ab = [ -M(1,2); M(1,1) ];
  else
    ab = [ -M(2,2); M(2,1) ];
  end
  a = ab(1);
  b = ab(2);

  % Normalize so the surface displacement etaHat = (i/omega) vHat(H) = ampfs is real
  alpha = sqrt( k^2 - 1i*omega/nu );
  q = 1i*k/(rho*omega);
  coshAH = cosh(alpha*H);
  sinhAH = sinh(alpha*H);
  vHatH = q*( a*( (k/alpha)*sinhAH - sinhKH ) + b*( coshAH - coshKH ) );
  etaHat = (1i/omega)*vHatH;

  scale = par.ampfs/etaHat;
  a = a*scale;
  b = b*scale;
  vHatH = vHatH*scale;
  c1 = q*b;
  c2 = q*(k/alpha)*a;

  % Thesis constants (C.26)-(C.29), scaled so the surface displacement |(uHat,vHat)|/|omega| = 1 at y=0 as in (C.36), for comparison with Table 4.1.
  % The thesis phase is not known so only ratios such as B/A, C/A, D/A are directly comparable.
  uHatH = (1i/k)*q*( a*k*( coshAH - coshKH ) + b*( alpha*sinhAH - k*sinhKH ) );
  dispMag = sqrt( abs(uHatH)^2 + abs(vHatH)^2 )/abs(omega);
  A = 0.5*(a + b)*exp( k*H)/dispMag;
  B = 0.5*(a - b)*exp(-k*H)/dispMag;
  C = 0.5*(c1 + c2)*exp( alpha*H)/dispMag;
  D = 0.5*(c1 - c2)*exp(-alpha*H)/dispMag;

  par.DispersionRelation = omega;

  % Data needed by defineKnownSolution to build the exact solution
  vcw.omega = omega;
  vcw.alpha = alpha;
  vcw.k     = k;
  vcw.H     = H;
  vcw.q     = q;
  vcw.a     = a;
  vcw.b     = b;
  vcw.c1    = c1;
  vcw.c2    = c2;
  vcw.ABCD  = [A; B; C; D];

  if( par.idebug>0 )
    fprintf('getViscousCapillaryWaveCoeffs: k=%g, H=%g, nu=%g, rho=%g, g=%g, gamma=%g\n',k,H,nu,rho,g,gamma);
    if( ~isempty(rootBox) )
      fprintf('  search box: Re(omega) in [%g,%g], Im(omega) in [%g,%g], roots in box=%d\n',searchBox(1),searchBox(2),searchBox(3),searchBox(4),searchInfo.numRoots);
      fprintf('  subdivision levels=%d, evaluations of F=%d\n',searchInfo.numLevels,searchInfo.numFEvals);
    end
    fprintf('  omega0 = %10.4e %+10.4ei, omega = %10.4e %+10.4ei (secant its=%d)\n',real(omega0),imag(omega0),real(omega),imag(omega),it);
    fprintf('  A/u0=%10.4e%+10.4ei, B/u0=%10.4e%+10.4ei, C/u0=%10.4e%+10.4ei, D/u0=%10.4e%+10.4ei\n',...
            real(A),imag(A),real(B),imag(B),real(C),imag(C),real(D),imag(D));
    fprintf('  B/A=%10.4e%+10.4ei, C/A=%10.4e%+10.4ei, D/A=%10.4e%+10.4ei\n',...
            real(B/A),imag(B/A),real(C/A),imag(C/A),real(D/A),imag(D/A));
  end

  return


  % F(omega) = omega^2 * det( M unscaled ), an entire function of omega whose zeros are the roots of
  % the dispersion relation (the 1/cosh(alpha H) and 1/omega factors in M are removed so that F has
  % no poles, as needed by the argument principle in findRootBox). Every term is even in alpha, so
  % the branch of the square root does not matter. All terms are multiplied by the positive factor
  % exp(-Re(alpha) H), which avoids overflow and does not change arg(F).
  function F = getCountFunction( omega )

    alphaw = sqrt( k^2 - 1i*omega/nu );
    scaleAH = exp( -real(alphaw)*H );
    phaseAH = exp( 1i*imag(alphaw)*H );
    decayAH = exp( -2*real(alphaw)*H )/phaseAH;
    coshAHs = 0.5*( phaseAH + decayAH );   % cosh(alpha H) exp(-Re(alpha) H)
    sinhAHs = 0.5*( phaseAH - decayAH );   % sinh(alpha H) exp(-Re(alpha) H)
    coshKHs = coshKH*scaleAH;
    sinhKHs = sinhKH*scaleAH;

    vHatA   = (k/alphaw)*sinhAHs - sinhKHs;
    vHatB   = coshAHs - coshKHs;
    vHatyA  = k*( coshAHs - coshKHs );
    vHatyB  = alphaw*sinhAHs - k*sinhKHs;
    vHatyyA = k*( alphaw*sinhAHs - k*sinhKHs );
    vHatyyB = alphaw^2*coshAHs - k^2*coshKHs;

    % shear stress row (factor q dropped)
    F11 = vHatyyA + k^2*vHatA;
    F12 = vHatyyB + k^2*vHatB;
    % normal stress row times omega^2: omega^2 q = i k omega/rho and omega^2 q (i G/omega) = -k G/rho
    F21 = -omega^2*coshKHs + (1i*k*omega/rho)*2*mu*vHatyA - (k*G/rho)*vHatA;
    F22 = -omega^2*sinhKHs + (1i*k*omega/rho)*2*mu*vHatyB - (k*G/rho)*vHatB;

    F = F11*F22 - F12*F21;

  end


  % Build the 2x2 matrix M(omega), both rows scaled by 1/cosh(alpha H).
  % Only one complex exponential, exp(-alpha H), is needed per evaluation.
  function Mw = getMatrix( omega )

    alphaw = sqrt( k^2 - 1i*omega/nu );    % Re(alphaw) >= 0 so eAH is bounded
    qw     = 1i*k/(rho*omega);

    eAH    = exp(-alphaw*H);
    e2AH   = eAH*eAH;
    tanhAH = (1 - e2AH)/(1 + e2AH);
    sechAH = 2*eAH/(1 + e2AH);
    coshKHsechAH = coshKH*sechAH;      % cosh(k H)/cosh(alpha H)
    sinhKHsechAH = sinhKH*sechAH;      % sinh(k H)/cosh(alpha H)

    % At s=H, scaled by 1/cosh(alpha H):
    %   vHat   = qw*( a*vHatA   + b*vHatB   )
    %   vHat'  = qw*( a*vHatyA  + b*vHatyB  )
    %   vHat'' = qw*( a*vHatyyA + b*vHatyyB )
    vHatA   = (k/alphaw)*tanhAH - sinhKHsechAH;
    vHatB   = 1 - coshKHsechAH;
    vHatyA  = k*(1 - coshKHsechAH);
    vHatyB  = alphaw*tanhAH - k*sinhKHsechAH;
    vHatyyA = k*(alphaw*tanhAH - k*sinhKHsechAH);
    vHatyyB = alphaw^2 - k^2*coshKHsechAH;

    Mw = zeros(2,2);
    % shear stress: vHat'' + k^2 vHat = 0  (common factor qw dropped)
    Mw(1,1) = vHatyyA + k^2*vHatA;
    Mw(1,2) = vHatyyB + k^2*vHatB;
    % normal stress: -pHat + 2 mu vHat' + (i G/omega) vHat = 0
    Mw(2,1) = -coshKHsechAH + qw*( 2*mu*vHatyA + (1i*G/omega)*vHatA );
    Mw(2,2) = -sinhKHsechAH + qw*( 2*mu*vHatyB + (1i*G/omega)*vHatB );

  end

end
