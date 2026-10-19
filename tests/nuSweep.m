%
% nu sweep through getViscousCapillaryWaveCoeffs with omega0=0 (argument principle + secant).
%
%   addpath(genpath(pwd)); nuSweep
%

par = struct('kx',2*pi,'xa',0,'xb',1,'ya',-1,'yb',0,'nu',0,'rho',1,'gravity',-1,'gamma',0.05,'ampfs',1e-4,'omega0',0,'idebug',0);
nuValues = [0.01, 0.02, 0.05, 0.1, 0.12, 0.14, 0.15, 0.2, 0.5, 1];

% Expected omega for nu<=0.1 from the earlier findRootBox search (4 digits); NaN where none is known.
% For nu>=0.15 the wave is overdamped and no root is expected.
omegaExpected = [4.1620-0.6131i, 3.8952-1.0588i, 3.0446-1.7730i, 2.0163-2.4054i, NaN, NaN, NaN, NaN, NaN, NaN];

fprintf('\nnu sweep through getViscousCapillaryWaveCoeffs (gamma=0.05, k=2*pi, H=1)\n');
fprintf('%5s %-24s %-24s %10s %s\n','nu','omega','expected','|diff|','note');
for m=1:length(nuValues)
  par.nu = nuValues(m);
  note = '';
  try
    tic;
    [~, vcw] = getViscousCapillaryWaveCoeffs( par );
    omega = vcw.omega;
    note = sprintf('%.3fs',toc);
  catch e
    omega = NaN;
    note = e.message;
  end
  fprintf('%5.2f %10.4f%+10.4fi   %10.4f%+10.4fi   %10.2e %s\n',par.nu,real(omega),imag(omega),...
          real(omegaExpected(m)),imag(omegaExpected(m)),abs(omega-omegaExpected(m)),note);
end
