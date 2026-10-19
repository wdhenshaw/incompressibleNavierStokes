%
% Stage 1 test: compare getViscousCapillaryWaveCoeffs with Table 4.1 of Dan's thesis.
%
%   addpath(genpath(pwd)); testViscousCapillaryWaveCoeffs
%

gammaValues = [0.05; 0.1; 0.2];

% Table 4.1: omega, A/u0, B/u0, C/u0, D/u0 (B and D are listed in the table scaled by 1e5)
omegaTable = [ 3.045-1.773i;  4.309-2.043i;  6.277-2.337i ];
ATable     = [-1.879-0.4893i; -3.247+0.1056i; -6.099+0.7879i ];
BTable     = [ 1.296+1.793i;  -2.521-0.3937i; -3.843-3.253i  ]*1e-5;
CTable     = [ 3.606+0.5887i;  3.587-0.02137i; 3.514-0.08220i ];
DTable     = [ 3.835-4.201i;  -1.968+0.2356i;  0.2423-0.8369i ]*1e-5;

% Parameters of Table 4.1: L=H=1, nu=0.05, kx=2*pi, rho=1, g=1
par = struct();
par.kx      = 2*pi;
par.xa      = 0;
par.xb      = 1;
par.ya      = -1;
par.yb      = 0;
par.nu      = 0.05;
par.rho     = 1;
par.gravity = -1;
par.ampfs   = 1e-4;
par.idebug  = 0;

k = par.kx;
H = par.yb - par.ya;
g = -par.gravity;

fprintf('\nStage 1: getViscousCapillaryWaveCoeffs vs Table 4.1 (L=H=1, nu=%g, k=2*pi, rho=%g, g=%g)\n',par.nu,par.rho,g);
fprintf('The thesis phase is unknown, so only |A|,|B|,|C|,|D| and the ratios B/A, C/A, D/A are compared.\n');

for m=1:length(gammaValues)

  par.gamma = gammaValues(m);
  tableABCD = [ATable(m); BTable(m); CTable(m); DTable(m)];

  % Default guess (same formula as getViscousCapillaryWaveCoeffs)
  defaultGuess = sqrt( (g*k + par.gamma*k^3/par.rho)*tanh(k*H) ) - 1i*2*par.nu*k^2;

  guessNames  = {'default guess', 'table omega as guess'};
  guessValues = [0, omegaTable(m)];

  for n=1:2

    par.omega0 = guessValues(n);

    lastwarn('');
    [~, vcw] = getViscousCapillaryWaveCoeffs( par );
    warned = ~isempty(lastwarn);

    ABCD = vcw.ABCD;

    fprintf('\ngamma=%g, %s',par.gamma,guessNames{n});
    if( n==1 )
      fprintf(' (omega0=%.4f%+.4fi), warning issued: %s',real(defaultGuess),imag(defaultGuess),yesNo(warned));
    end
    fprintf('\n');
    fprintf('  %-10s %23s %23s %10s\n','quantity','computed','table','rel. err');
    printRow('omega', vcw.omega, omegaTable(m));
    printRow('|A|',   abs(ABCD(1)), abs(tableABCD(1)));
    printRow('|B|',   abs(ABCD(2)), abs(tableABCD(2)));
    printRow('|C|',   abs(ABCD(3)), abs(tableABCD(3)));
    printRow('|D|',   abs(ABCD(4)), abs(tableABCD(4)));
    printRow('B/A',   ABCD(2)/ABCD(1), tableABCD(2)/tableABCD(1));
    printRow('C/A',   ABCD(3)/ABCD(1), tableABCD(3)/tableABCD(1));
    printRow('D/A',   ABCD(4)/ABCD(1), tableABCD(4)/tableABCD(1));
    % In case the thesis used the other branch of alpha, which swaps C and D
    printRow('D/A (swap)', ABCD(4)/ABCD(1), tableABCD(3)/tableABCD(1));
    printRow('C/A (swap)', ABCD(3)/ABCD(1), tableABCD(4)/tableABCD(1));

  end
end
fprintf('\n');


function printRow( name, computed, tableValue )
  relErr = abs(computed - tableValue)/abs(tableValue);
  fprintf('  %-10s %11.4e%+11.4ei %11.4e%+11.4ei %10.2e\n',name,real(computed),imag(computed),real(tableValue),imag(tableValue),relErr);
end

function s = yesNo( flag )
  if( flag )
    s = 'yes';
  else
    s = 'no';
  end
end
