%
% Compute the dispersion relation and coefficients for the linearized viscous
% gravity-capillary traveling wave (Serino thesis, Appendix C.2.1).
%
%   Domain: [xa,xb] x [-H,0], no-slip at y=-H, free surface linearized about y=0.
%
%   v(x,y,t) = vHat(y) exp( i(k x - omega t) ),   p(x,y,t) = -rho g y + pHat(y) exp( i(k x - omega t) )
%
% The thesis writes the coefficient functions with exponentials e^{+-k y}, e^{+-alpha y}. Here we
% use hyperbolic functions shifted to the bottom wall, s = y + H, so that the no-slip conditions
% at s=0 eliminate two of the four constants in closed form:
%
%   alpha^2 = k^2 - i omega/nu          (C.27; note (4.75) has the sign reversed)
%   q       = i k/(rho omega)
%   pHat    = a cosh(k s) + b sinh(k s)
%   v2Hat   = q*[ a( (k/alpha) sinh(alpha s) - sinh(k s) ) + b( cosh(alpha s) - cosh(k s) ) ]
%   v1Hat   = (i/k) v2Hat'
%
% (equivalently c1 = q b, c2 = q (k/alpha) a are the coefficients of cosh(alpha s), sinh(alpha s)).
% The free surface conditions at s=H give a 2x2 system M(omega)[a;b] = 0:
%
%   shear :  v2Hat'' + k^2 v2Hat = 0
%   normal: -pHat + 2 mu v2Hat' + (i G/omega) v2Hat = 0,     G = rho g + gamma k^2
%
% Both rows are divided by cosh(alpha H) so the entries stay O(1) for small nu. All entries are
% even in alpha so the branch of the square root does not matter.
%
% The constants are normalized so the surface displacement at t=0 is real:
%    eta(x,0) = ampfs*cos(k x),   eta = Re( etaHat e^{i(kx-omega t)} ),  etaHat = (i/omega) v2Hat(H) = ampfs
%
% Input (par): kx (already scaled by 2*pi), xa,xb,ya,yb, nu, rho, gravity, gamma, ampfs,
%              omega0 (optional initial guess for the complex frequency).
% Output: par.DispersionRelation = omega, and par.vcw with fields
%    omega, alpha, k, H, a, b, c1, c2, etaHat, ABCD (thesis exponential constants, see below)
%
function par = getViscousCapillaryWaveCoeffs( par )

  k   = par.kx;
  H   = abs(par.yb - par.ya);
  nu  = par.nu;
  rho = par.rho;
  mu  = nu*rho;
  g   = -par.gravity;   % main code uses gravity<0, the TW solution uses g>0
  gamma = par.gamma;
  G   = rho*g + gamma*k^2;

  % --- omega independent quantities (computed once) ---
  Ck = cosh(k*H);
  Sk = sinh(k*H);

  % --- initial guess ---
  if( isfield(par,'omega0') && par.omega0~=0 )
    omega0 = par.omega0;
  else
    % inviscid frequency (same as GravityCapillaryWave) plus viscous decay 2 nu k^2
    omega0 = sqrt( (g*k + gamma*k^3/rho)*tanh(k*H) ) - 1i*2*nu*k^2;
  end

  detM = @(omega) det( getMatrix(omega) );

  % --- secant iteration for det M(omega) = 0 in the complex plane ---
  tol = 1e-13;  maxIt = 100;
  w0 = omega0;  w1 = omega0*(1 + 1e-3);
  f0 = detM(w0); f1 = detM(w1);
  converged = 0;
  for it=1:maxIt
    dw = -f1*(w1-w0)/(f1-f0);
    w0 = w1;  f0 = f1;
    w1 = w1 + dw;
    f1 = detM(w1);
    if( abs(dw) < tol*abs(w1) )
      converged = 1;
      break;
    end
  end
  if( ~converged )
    error('getViscousCapillaryWaveCoeffs: secant iteration did not converge, omega=%g%+gi, |dw|=%8.2e\n',real(w1),imag(w1),abs(dw));
  end
  omega = w1;

  % --- null vector of M ---
  M = getMatrix(omega);
  if( norm(M(1,:)) >= norm(M(2,:)) )
    ab = [ -M(1,2); M(1,1) ];
  else
    ab = [ -M(2,2); M(2,1) ];
  end
  a = ab(1); b = ab(2);

  % --- normalize so etaHat = (i/omega) v2Hat(H) = ampfs ---
  alpha = sqrt( k^2 - 1i*omega/nu );
  q  = 1i*k/(rho*omega);
  Ca = cosh(alpha*H);
  Sa = sinh(alpha*H);
  v2H = q*( a*( (k/alpha)*Sa - Sk ) + b*( Ca - Ck ) );
  etaHat = (1i/omega)*v2H;

  scale = par.ampfs/etaHat;
  a = a*scale;  b = b*scale;
  c1 = q*b;
  c2 = q*(k/alpha)*a;

  % --- thesis constants (C.26)-(C.29) for comparison with Table 4.1 ---
  %   pHat = A e^{k y} + B e^{-k y},  homogeneous v2Hat = C e^{alpha y} + D e^{-alpha y}
  %   NOTE: the thesis normalizes |uHat(0)| = u0 with an unknown phase, so only ratios such as
  %   B/A, C/A, D/A are directly comparable. C and D swap if alpha -> -alpha.
  A = 0.5*(a + b)*exp( k*H);
  B = 0.5*(a - b)*exp(-k*H);
  C = 0.5*(c1 + c2)*exp( alpha*H);
  D = 0.5*(c1 - c2)*exp(-alpha*H);
  % scale to |uHat(0)| = 1 as in (C.36):  uHat = vHat/(-i omega)
  v1H = (1i/k)*q*( a*k*( Ca - Ck ) + b*( alpha*Sa - k*Sk ) );
  u0 = sqrt( abs(v1H)^2 + abs(v2H*scale)^2 )/abs(omega);
  ABCD = [A; B; C; D]/u0;

  par.DispersionRelation = omega;

  par.vcw.omega  = omega;
  par.vcw.alpha  = alpha;
  par.vcw.k      = k;
  par.vcw.H      = H;
  par.vcw.a      = a;
  par.vcw.b      = b;
  par.vcw.c1     = c1;
  par.vcw.c2     = c2;
  par.vcw.etaHat = par.ampfs;
  par.vcw.ABCD   = ABCD;

  if( par.idebug>0 )
    fprintf('--- getViscousCapillaryWaveCoeffs: k=%g, H=%g, nu=%g, rho=%g, g=%g, gamma=%g ---\n',k,H,nu,rho,g,gamma);
    fprintf('    omega0 = %10.4e %+10.4ei, omega = %10.4e %+10.4ei (its=%d)\n',real(omega0),imag(omega0),real(omega),imag(omega),it);
    fprintf('    A/u0=%10.4e%+10.4ei, B/u0=%10.4e%+10.4ei, C/u0=%10.4e%+10.4ei, D/u0=%10.4e%+10.4ei\n',...
            real(ABCD(1)),imag(ABCD(1)),real(ABCD(2)),imag(ABCD(2)),real(ABCD(3)),imag(ABCD(3)),real(ABCD(4)),imag(ABCD(4)));
    fprintf('    B/A=%10.4e%+10.4ei, C/A=%10.4e%+10.4ei, D/A=%10.4e%+10.4ei\n',...
            real(B/A),imag(B/A),real(C/A),imag(C/A),real(D/A),imag(D/A));
  end

  return


  % ---------------------------------------------------------------------
  % Build the 2x2 matrix M(omega), both rows scaled by 1/cosh(alpha H).
  % Only one complex exponential, exp(-alpha H), is needed per evaluation.
  % ---------------------------------------------------------------------
  function M = getMatrix( omega )

    alpha = sqrt( k^2 - 1i*omega/nu );    % Re(alpha) >= 0 so em is bounded
    q     = 1i*k/(rho*omega);

    em  = exp(-alpha*H);
    em2 = em*em;
    Ta  = (1 - em2)/(1 + em2);            % tanh(alpha H)
    sech= 2*em/(1 + em2);                 % 1/cosh(alpha H)
    ck  = Ck*sech;                        % cosh(k H)/cosh(alpha H)
    sk  = Sk*sech;                        % sinh(k H)/cosh(alpha H)

    % v2Hat   = q*[ a*P2 + b*Q2 ]   (all scaled by 1/cosh(alpha H))
    % v2Hat'  = q*[ a*P1 + b*Q1 ]
    % v2Hat'' = q*[ a*P0 + b*Q0 ]
    P2 = (k/alpha)*Ta - sk;      Q2 = 1 - ck;
    P1 = k*(1 - ck);             Q1 = alpha*Ta - k*sk;
    P0 = k*(alpha*Ta - k*sk);    Q0 = alpha^2 - k^2*ck;

    M = zeros(2,2);
    % shear stress: v2Hat'' + k^2 v2Hat = 0  (common factor q dropped)
    M(1,1) = P0 + k^2*P2;
    M(1,2) = Q0 + k^2*Q2;
    % normal stress: -pHat + 2 mu v2Hat' + (i G/omega) v2Hat = 0
    M(2,1) = -ck + q*( 2*mu*P1 + (1i*G/omega)*P2 );
    M(2,2) = -sk + q*( 2*mu*Q1 + (1i*G/omega)*Q2 );

  end

end
