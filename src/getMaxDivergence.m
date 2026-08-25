function [maxDivU,maxGradU,par] = getMaxDivergence( un,vn,par )


  dx = par.dx;
  dy = par.dy;
  dr = par.dr(1);
  ds = par.dr(2);   

  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r to second order 
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s    

  [I1,I2] = getIndex( par.gid );

  if( par.isCartesian )

    ux = zeros(par.Ngx,par.Ngy);
    uy = zeros(par.Ngx,par.Ngy);


    ux(I1,I2) = Dzx(un,I1,I2);
    vy(I1,I2) = Dzy(vn,I1,I2);

    maxDivU = max( abs( ux(I1,I2) + vy(I1,I2) ), [], "all");

    maxGradU = max( abs(ux(I1,I2)) + abs(Dzy(un,I1,I2)) + abs(Dzx(vn,I1,I2)) + abs(vy(I1,I2)) , [], "all");

  else

    % --- curvilinear ---
    maxDivU  = 0.;
    maxGradU = 0.;
    for( i1=I1 )
    for( i2=I2 ) 
      rx = par.rx(i1,i2,1,1);
      ry = par.rx(i1,i2,1,2);
      sx = par.rx(i1,i2,2,1);
      sy = par.rx(i1,i2,2,2);

      ur = Dr2(un,i1,i2); us = Ds2(un,i1,i2);
      vr = Dr2(vn,i1,i2); vs = Ds2(vn,i1,i2);

      ux = rx*ur + sx*us; 
      vx = rx*vr + sx*vs;

      vy = ry*vr + sy*vs;
      uy = ry*ur + sy*us;

      div = ux + vy;

      maxDivU = max( maxDivU,div );
      maxGradU = max( maxGradU, abs(ux) + abs(uy) + abs(vx) + abs(vy) );

    end
    end

  end

  return
end