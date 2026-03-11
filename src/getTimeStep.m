%
%  Compute the time-step dt 
%
function [dt,par] = getTimeStep( step,un,vn, par)

  dx = par.dx;
  dy = par.dy;
  nu = par.nu;

  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y  

  % Choose dt: assume stability region is an ellipse:
  %      (reLambda*dt/aStab)^2 + (imLambda*dt/bStab)^2 =1 
  if( strcmp(par.ts,'fe') ) 
    aStab=-2; bStab=1;      % approximate FE - not true if kappa=0
  elseif( strcmp(par.ts,'rk2') )
    aStab=-2; bStab=1.5;    % approximate RK2
  elseif( strcmp(par.ts,'rk4') )
    aStab=-2.7; bStab=2.7;  % approximate RK4
  elseif( strcmp(par.ts,'ab2')' )
    aStab=-1; bStab=.8;     % approximate AB2
  elseif( strcmp(par.ts,'pc2')' )
    aStab=-1.75; bStab=1.3; % approximate PC2
  elseif( strcmp(par.ts,'im2')' )
    aStab=-1e20; bStab=1.7; % approximate IM2 *check me*
  else
    fprintf('ERROR: unknown par.ts=[%s]\n',par.ts);
    pause;
  end;


  reLambda=4*nu*(1/dx^2+1/dy^2);    % real part of time-stepping eigenvalue
  if( strcmp(par.ts,'im2') )
    reLambda=0.; 
  end 
  if( par.ad )
    % See cg/ins/src/indts.bf

    ad21 = par.ad21;
    cd22 = par.ad22/(par.nd^2);

    [I1,I2] = getIndex( par.gid );
    maxGradU = max( abs(Dzx(un,I1,I2)) + abs(Dzy(un,I1,I2)) + abs(Dzx(vn,I1,I2)) + abs(Dzy(vn,I1,I2)) , [], "all");
    reLambda=reLambda + 8.*( ad21 + cd22*maxGradU );

  end 
  uMax = max(max(abs(un)));
  vMax = max(max(abs(vn)));
  imLambda= ( uMax/dx + vMax/dy );


  dt = par.cfl/sqrt( (reLambda/aStab)^2 + (imLambda/bStab)^2 ); % adjusted a bit below to reach tFinal

  dt = min( dt,par.dtMax );

  if( par.idebug >0 && step==0 )
    fprintf('getTimeStep: dt = %9.3e (dtMax=%9.2e)\n',dt,par.dtMax);
  end

  par.dtOld = dt;

end
