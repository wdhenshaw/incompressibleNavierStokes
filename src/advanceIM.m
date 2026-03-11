%
% IMEX scheme
%
function [unp1,vnp1,pnp1,par,ut,vt] = advanceIM( t,dt, un,vn,pn,unp1,vnp1,pnp1, ut,vt,par )

  % declare operators 
  % --- Difference Operators ---
  % Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  % Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  dx = par.dx;
  dy = par.dy;

  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy


  tnp1 = t + dt;
  nu = par.nu;
 
  [I1,I2] = getIndex( par.gid );

  utm=ut; vtm=vt; % save old du/dt 

  % ====== Predictor =====
  nuScaleFactor=0.;  % leave off viscous terms
  [ ut,vt,par ] = getUt( t,un,vn,pn,nuScaleFactor,par );
  
  unp1(I1,I2) = un(I1,I2) + par.ab1*ut(I1,I2) +par.ab2*utm(I1,I2) + (.5*dt*nu)*( DpxDmx(un,I1,I2) + DpyDmy(un,I1,I2) );
  vnp1(I1,I2) = vn(I1,I2) + par.ab1*vt(I1,I2) +par.ab2*vtm(I1,I2) + (.5*dt*nu)*( DpxDmx(vn,I1,I2) + DpyDmy(vn,I1,I2) ); 

  [unp1,vnp1,par] = solveImplicitTimeStep( unp1,vnp1,tnp1, par );
  [unp1,vnp1,par] = applyBoundaryConditions( unp1,vnp1,tnp1,par );   

  % --- solve the pressure equation ---
  factorMatrix=0;
  [pnp1,par] = pressureEquation( tnp1, unp1,vnp1,dt,factorMatrix,par ); 

  % ====== CORRECTOR : Crank-Nicolson =====
  nuScaleFactor=.0;  % leave off viscous terms
  [ utp,vtp,par ] = getUt( tnp1,unp1,vnp1,pnp1,nuScaleFactor,par );
 
  unp1(I1,I2) = un(I1,I2) + (.5*dt)*(utp(I1,I2)+ut(I1,I2)) + (.5*dt*nu)*( DpxDmx(un,I1,I2) + DpyDmy(un,I1,I2) );
  vnp1(I1,I2) = vn(I1,I2) + (.5*dt)*(vtp(I1,I2)+vt(I1,I2)) + (.5*dt*nu)*( DpxDmx(vn,I1,I2) + DpyDmy(vn,I1,I2) ); 

  [unp1,vnp1,par] = solveImplicitTimeStep( unp1,vnp1,tnp1, par );
  [unp1,vnp1,par] = applyBoundaryConditions( unp1,vnp1,tnp1,par );

  % --- solve the pressure equation ---
  [pnp1,par] = pressureEquation( tnp1, unp1,vnp1,dt,factorMatrix,par ); 



  return
end