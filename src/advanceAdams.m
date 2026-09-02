%
% Adams Bashforth -- Grid Function version 
%
function [gf,par,ut,vt] = advanceAdams( t,dt, gf,cur,next, ut,vt,par )

  tnp1 = t + dt;

  [I1,I2] = getIndex( par.gid );

  utm=ut; vtm=vt; % save old du/dt 
  [ut,vt,par] = getUt( t,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, par.nuScaleFactor,par );

  % if( 1==1 )
  %   utErr = max(max(abs( ut(I1,I2) - par.uet(gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),t) ) ));
  %   vtErr = max(max(abs( vt(I1,I2) - par.vet(gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),t) ) ));
  %   fprintf('advanceAdams: t=%9.3e cur=%d [utErr,vtErr]=[%9.3e,%9.3e]\n',t,cur, utErr,vtErr);
  % end

  gf{next}.u(I1,I2) = gf{cur}.u(I1,I2) + par.ab1*ut(I1,I2) + par.ab2*utm(I1,I2); 
  gf{next}.v(I1,I2) = gf{cur}.v(I1,I2) + par.ab1*vt(I1,I2) + par.ab2*vtm(I1,I2); 

  if( 1==0 ) %  && t<=1*dt )
    fprintf('advanceAdams: Set (u,v) to exact at tnp1=%9.3e \n',tnp1);
    gf{next}.u(I1,I2) = par.ue(gf{next}.x(I1,I2,1),gf{next}.x(I1,I2,2),tnp1);
    gf{next}.v(I1,I2) = par.ve(gf{next}.x(I1,I2,1),gf{next}.x(I1,I2,2),tnp1);
  end

  [gf{next}.u,gf{next}.v,par] = applyBoundaryConditions( gf{next}.u,gf{next}.v,tnp1, gf,next, par );

  % if( 1==1 ) %  && t<=1*dt )
  %   fprintf('advanceAdams: Set (u,v) to exact at tnp1=%9.3e \n',tnp1);
  %   [I1,I2] = getIndex( par.dim );
  %   gf{next}.u(I1,I2) = par.ue(gf{next}.x(I1,I2,1),gf{next}.x(I1,I2,2),tnp1);
  %   gf{next}.v(I1,I2) = par.ve(gf{next}.x(I1,I2,1),gf{next}.x(I1,I2,2),tnp1);
  % end

  % --- solve the pressure equation ---
  if( par.gridMotion~=par.noMotion) par.factorPressureMatrix=1; end
  [gf{next}.p,par] = pressureEquation( tnp1, gf{next}.u,gf{next}.v,dt, gf,next, par );   




  return
end