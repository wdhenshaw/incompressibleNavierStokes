%
% AB2+AM2 Predictor corrector
%
function [unp1,vnp1,pnp1,par,ut,vt] = advancePC( t,dt, un,vn,pn,unp1,vnp1,pnp1, ut,vt,par )

  % fprintf('advancePC: start un=[%d,%d] vn=[%d,%d]\n',...
  %        size(un,1),size(un,2), size(vn,1),size(vn,2) );

  tnp1 = t + dt;
 
  [I1,I2] = getIndex( par.gid );

  utm=ut; vtm=vt; % save old du/dt 

  % ====== AB2 Predictor =====
  [ut,vt,par] = getUt( t,un,vn,pn,par.nuScaleFactor,par );
  unp1(I1,I2) = un(I1,I2) + par.ab1*ut(I1,I2) + par.ab2*utm(I1,I2); 
  vnp1(I1,I2) = vn(I1,I2) + par.ab1*vt(I1,I2) + par.ab2*vtm(I1,I2); 

  [unp1,vnp1,par] = applyBoundaryConditions( unp1,vnp1,tnp1,par );

  % --- solve the pressure equation ---
  factorMatrix=0;
  [pnp1,par] = pressureEquation( tnp1, unp1,vnp1,dt,factorMatrix,par ); 


  % ====== AM2 CORRECTOR =====
  [ utp,vtp,par ] = getUt( tnp1,unp1,vnp1,pnp1,par.nuScaleFactor,par );
   
  unp1(I1,I2) = un(I1,I2) + (.5*dt)*(utp(I1,I2)+ut(I1,I2));
  vnp1(I1,I2) = vn(I1,I2) + (.5*dt)*(vtp(I1,I2)+vt(I1,I2)); 

  [unp1,vnp1,par] = applyBoundaryConditions( unp1,vnp1,tnp1,par );

  % --- solve the pressure equation ---
  [pnp1,par] = pressureEquation( tnp1, unp1,vnp1,dt,factorMatrix,par ); 


  return
end