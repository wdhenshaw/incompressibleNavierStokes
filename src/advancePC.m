%
% AB2+AM2 Predictor corrector -- Grid Function version 
%
function [gf,par,ut,vt] = advancePC( t,dt, gf,cur,next, ut,vt,par )

  % fprintf('advancePC: start un=[%d,%d] vn=[%d,%d]\n',...
  %        size(un,1),size(un,2), size(vn,1),size(vn,2) );

  tnp1 = t + dt;
 
  [I1,I2] = getIndex( par.gid );

  utm=ut; vtm=vt; % save old du/dt 

  % ====== AB2 Predictor =====
  [ut,vt,par] = getUt( t,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, par.nuScaleFactor,par );

  gf{next}.u(I1,I2) = gf{cur}.u(I1,I2) + par.ab1*ut(I1,I2) + par.ab2*utm(I1,I2); 
  gf{next}.v(I1,I2) = gf{cur}.v(I1,I2) + par.ab1*vt(I1,I2) + par.ab2*vtm(I1,I2); 

  [gf{next}.u,gf{next}.v,par] = applyBoundaryConditions( gf{next}.u,gf{next}.v,tnp1, gf,next, par );

  % --- solve the pressure equation ---
  if( par.gridMotion~=par.noMotion) par.factorPressureMatrix=1; end    
  [gf{next}.p,par] = pressureEquation( tnp1, gf{next}.u,gf{next}.v,dt, gf,next, par );

  % ====== AM2 CORRECTOR =====
  [ utp,vtp,par ] = getUt( tnp1,gf{next}.u,gf{next}.v,gf{next}.p, gf,next, par.nuScaleFactor,par );
   
  gf{next}.u(I1,I2) = gf{cur}.u(I1,I2) + (.5*dt)*(utp(I1,I2)+ut(I1,I2));
  gf{next}.v(I1,I2) = gf{cur}.v(I1,I2) + (.5*dt)*(vtp(I1,I2)+vt(I1,I2)); 

  [gf{next}.u,gf{next}.v,par] = applyBoundaryConditions( gf{next}.u,gf{next}.v,tnp1, gf,next, par );

  % --- solve the pressure equation ---
  [gf{next}.p,par] = pressureEquation( tnp1, gf{next}.u,gf{next}.v,dt, gf,next,par ); 


  return
end