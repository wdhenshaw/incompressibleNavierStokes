%
% IMEX scheme -- Grid Function version 
%
function [gf,par,ut,vt] = advanceIM( t,dt, gf,cur,next, ut,vt,par )

  % declare operators 
  % --- Difference Operators ---
  % Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  % Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  dx = par.dx;
  dy = par.dy;

  % DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  % DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy


  tnp1 = t + dt;
  nu = par.nu;
 
  [I1,I2] = getIndex( par.gid );

  utm=ut; vtm=vt; % save old du/dt 

  % ====== Predictor =====
  nuScaleFactor=0.;  % leave off viscous terms
  [ ut,vt,par,uLap,vLap ] = getUt( t,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, nuScaleFactor,par );
  
  gf{next}.u(I1,I2) = gf{cur}.u(I1,I2) + par.ab1*ut(I1,I2) +par.ab2*utm(I1,I2) + (.5*dt*nu)*uLap(I1,I2);
  gf{next}.v(I1,I2) = gf{cur}.v(I1,I2) + par.ab1*vt(I1,I2) +par.ab2*vtm(I1,I2) + (.5*dt*nu)*vLap(I1,I2);

  if( par.gridMotion~=par.noMotion) par.factorImplicitMatrix=1; end
  [gf{next}.u,gf{next}.v,par] = solveImplicitTimeStep( gf{next}.u,gf{next}.v,tnp1, gf,next, par );

  [gf{next}.u,gf{next}.v,par] = applyBoundaryConditions( gf{next}.u,gf{next}.v,tnp1, gf,next, par );   

  % --- solve the pressure equation ---
  if( par.gridMotion~=par.noMotion) par.factorPressureMatrix=1; end   
  [gf{next}.p,par] = pressureEquation( tnp1, gf{next}.u,gf{next}.v,dt, gf,next,par ); 

  % ====== CORRECTOR : Crank-Nicolson =====
  nuScaleFactor=.0;  % leave off viscous terms
  [ utp,vtp,par ] = getUt( tnp1,gf{next}.u,gf{next}.v,gf{next}.p, gf,next, nuScaleFactor,par );
 
  gf{next}.u(I1,I2) = gf{cur}.u(I1,I2) + (.5*dt)*(utp(I1,I2)+ut(I1,I2)) + (.5*dt*nu)*uLap(I1,I2);
  gf{next}.v(I1,I2) = gf{cur}.v(I1,I2) + (.5*dt)*(vtp(I1,I2)+vt(I1,I2)) + (.5*dt*nu)*vLap(I1,I2);     


  [gf{next}.u,gf{next}.v,par] = solveImplicitTimeStep( gf{next}.u,gf{next}.v,tnp1, gf,next, par );
  [gf{next}.u,gf{next}.v,par] = applyBoundaryConditions( gf{next}.u,gf{next}.v,tnp1, gf,next, par );

  % --- solve the pressure equation ---
  [gf{next}.p,par] = pressureEquation( tnp1, gf{next}.u,gf{next}.v,dt, gf,next, par ); 



  return
end