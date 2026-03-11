function [maxDivU,maxGradU,par] = getMaxDivergence( un,vn,par )


  dx = par.dx;
  dy = par.dy;

  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  [I1,I2] = getIndex( par.gid );
  ux = zeros(par.Ngx,par.Ngy);
  uy = zeros(par.Ngx,par.Ngy);

  ux(I1,I2) = Dzx(un,I1,I2);
  vy(I1,I2) = Dzy(vn,I1,I2);

  maxDivU = max( abs( ux(I1,I2) + vy(I1,I2) ), [], "all");

  maxGradU = max( abs(ux(I1,I2)) + abs(Dzy(un,I1,I2)) + abs(Dzx(vn,I1,I2)) + abs(vy(I1,I2)) , [], "all");


  return
end