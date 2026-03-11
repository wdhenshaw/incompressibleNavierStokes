
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

  ampu=1.; ampp=-nu*ampu*8;
  par.ue = @(x,y,t)  (4.*ampu)*(1.-y).*y; 
  par.ve = @(x,y,t)  0.*x;
  par.pe = @(x,y,t)  ampp*x; 

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
