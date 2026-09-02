%
%  Compute the time-step dt 
%
function [dt,par] = getTimeStep( step,un,vn, gf,cur, par)

  dx = par.dx;
  dy = par.dy;
  dr = par.dr(1);
  ds = par.dr(2);  
  nu = par.nu;

  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y  

 % derivatives of entreis in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s

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

 [I1,I2] = getIndex( par.gid );

  if( par.isCartesian )

    reLambda=4*nu*(1/dx^2+1/dy^2);    % real part of time-stepping eigenvalue

    % uMax = max(max(abs(un)));
    % vMax = max(max(abs(vn)));

    imLambda= max( abs(un(I1,I2))/dx + abs(vn(I1,I2))/dy , [], 'all' );

  else
    % --- curvilinear ---
    imLambda = 0;
    reLambda = 0.;

    ds1Min = 1e50; ds2Min = ds1Min;
    ds1Max = 0;    ds2Max = ds1Max;

    for( i1=I1 )
    for( i2=I2 ) 
      rx = gf{cur}.rx(i1,i2,1,1);
      ry = gf{cur}.rx(i1,i2,1,2);
      sx = gf{cur}.rx(i1,i2,2,1);
      sy = gf{cur}.rx(i1,i2,2,2);

      % grid velocity 
      gv(1) = gf{cur}.gv(i1,i2,1); 
      gv(2) = gf{cur}.gv(i1,i2,2); 

      % u*u.x + v*u.y = u (rx*ur + sx*us) + v ( ry*ur + sy*us ) = (u*rx + v*ry) u.r + (u.rx + v*ry)* u.s
      im = abs( (un(i1,i2)-gv(1))*rx + (vn(i1,i2)-gv(2))*ry )/dr + abs( (un(i1,i2)-gv(1))*sx + (vn(i1,i2)-gv(2))*sy )/ds;

      rxr = DJzr(gf{cur}.rx,i1,i2,1,1);
      rxs = DJzs(gf{cur}.rx,i1,i2,1,1);
      ryr = DJzr(gf{cur}.rx,i1,i2,1,2);
      rys = DJzs(gf{cur}.rx,i1,i2,1,2);
      sxr = DJzr(gf{cur}.rx,i1,i2,2,1);
      sxs = DJzs(gf{cur}.rx,i1,i2,2,1);
      syr = DJzr(gf{cur}.rx,i1,i2,2,2);
      sys = DJzs(gf{cur}.rx,i1,i2,2,2);        

      rxx = rx*rxr + sx*rxs;
      ryy = ry*ryr + sy*rys;
      sxx = rx*sxr + sx*sxs;
      syy = ry*syr + sy*sys;

      % Lap(u) = c20*u.rr + c11*u.rs + c02*u.ss + c10*u.rr + c01*u.ss 
      c20 = (rx^2 + ry^2)/dr^2;
      c02 = (sx^2 + sy^2)/ds^2;
      c11 = 2*( rx*sx + ry*sy )/(4*dr*ds); 
      c10 = (rxx+ryy)/(2.*dr);
      c01 = (sxx+syy)/(2.*ds);

      re = nu*( c20*4 + abs(c11)*4 + c02*4 );

      im = im + nu*( abs(c10)*2 + abs(c01)*2 );

      imLambda = max( imLambda, im );
      reLambda = max( reLambda, re );

      % fprintf(' (i1,i2)=(%3d,%3d) rx/dr = %9.3e, sy/ds = %9.3e, im=%9.3e c10=%9.2e c01=%9.3e ry=%9.3e sx=%9.3e\n',i1,i2,rx/dr,sy/ds,im,c10,c01,ry,sx);
      % pause

      % estimate grid spacings in r and s directions
      ds1 = sqrt( (gf{cur}.x(i1+1,i2,1)-gf{cur}.x(i1,i2,1))^2 + (gf{cur}.x(i1+1,i2,2)-gf{cur}.x(i1,i2,2))^2 );
      ds2 = sqrt( (gf{cur}.x(i1,i2+1,1)-gf{cur}.x(i1,i2,1))^2 + (gf{cur}.x(i1,i2+1,2)-gf{cur}.x(i1,i2,2))^2 );

      ds1Min = min(ds1Min,ds1); ds2Min = min(ds2Min,ds2);
      ds1Max = max(ds1Max,ds1); ds2Max = max(ds2Max,ds2);

    end
    end

  end



  if( par.ad )
    % See cg/ins/src/indts.bf

    ad21 = par.ad21;
    cd22 = par.ad22/(par.nd^2);

    % [I1,I2] = getIndex( par.gid );
    % maxGradU = max( abs(Dzx(un,I1,I2)) + abs(Dzy(un,I1,I2)) + abs(Dzx(vn,I1,I2)) + abs(Dzy(vn,I1,I2)) , [], "all");

    [maxDivU,maxGradU,par] = getMaxDivergence( un,vn,par );

    reLambda=reLambda + 8.*( ad21 + cd22*maxGradU );

  end 

  % dt including viscous term: 
  dte = par.cfl/sqrt( (reLambda/aStab)^2 + (imLambda/bStab)^2 ); % adjusted a bit below to reach tFinal

  % dt with no viscous term for IMEX
  dti = par.cfl/sqrt( (imLambda/bStab)^2 );

  if( strcmp(par.ts,'im2') )
    dt = dti; 
  else
    dt = dte;
  end 


  dt = min( dt,par.dtMax );

  if( par.idebug >0 && step==0 )
    fprintf('getTimeStep: dt = %9.3e (dte=%9.3e, dti=%9.3e, dti/dte=%8.2e, dtMax=%9.2e)\n',dt,dte,dti,dti/dte,par.dtMax);
    fprintf('  imLambda=%9.3e, reLambda=%9.3e\n',imLambda,reLambda);
    if( ~ par.isCartesian )
      fprintf('  grid spacings: [ds1Min,ds2Min]=[%8.2e,%8.2e] [ds1Max,ds2Max]=[%8.2e,%8.2e]\n',ds1Min,ds2Min,ds1Max,ds2Max);
    end
    % pause
  end

  par.dtOld = dt;
 

end
