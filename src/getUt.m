%
% Compute the RHS to the momentum equations
%
function [ ut,vt,par ] = getUt( t,un,vn,pn,nuScaleFactor,par )

  cpu0 = cputime;

  dx = par.dx;
  dy = par.dy;
  nu = par.nu;

  % declare operators 
  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y
  
  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy
  

  [I1,I2] = getIndex( par.gid );

  ux = zeros(par.Ngx,par.Ngy);
  uy = zeros(par.Ngx,par.Ngy);
  vx = zeros(par.Ngx,par.Ngy);
  vy = zeros(par.Ngx,par.Ngy);

  ux(I1,I2) = Dzx(un,I1,I2);
  uy(I1,I2) = Dzy(un,I1,I2);
  vx(I1,I2) = Dzx(vn,I1,I2);
  vy(I1,I2) = Dzy(vn,I1,I2);

  if( nuScaleFactor~=0 )
    ut(I1,I2) = -( un(I1,I2).*ux(I1,I2) + vn(I1,I2).*uy(I1,I2) + Dzx(pn,I1,I2) ) ...
  	            + nu*( DpxDmx(un,I1,I2) + DpyDmy(un,I1,I2) ) ...
                + par.ufe(par.x(I1,I2,1),par.x(I1,I2,2),t);

    vt(I1,I2) = -( un(I1,I2).*vx(I1,I2) + vn(I1,I2).*vy(I1,I2) + Dzy(pn,I1,I2) ) ...
  	            + nu*( DpxDmx(vn,I1,I2) + DpyDmy(vn,I1,I2) ) ...
                + par.vfe(par.x(I1,I2,1),par.x(I1,I2,2),t);  

  else
    ut(I1,I2) = -( un(I1,I2).*ux(I1,I2) + vn(I1,I2).*uy(I1,I2) + Dzx(pn,I1,I2) ) ...
                +  par.ufe(par.x(I1,I2,1),par.x(I1,I2,2),t);

    vt(I1,I2) = -( un(I1,I2).*vx(I1,I2) + vn(I1,I2).*vy(I1,I2) + Dzy(pn,I1,I2) ) ...
                + par.vfe(par.x(I1,I2,1),par.x(I1,I2,2),t);  

  end
  % ut(I1,I2) =  -( un(I1,I2).*Dzx(un,I1,I2) + vn(I1,I2).*Dzy(un,I1,I2) + Dzx(pn,I1,I2) ) ...
  %       + (nu*nuScaleFactor)*( DpxDmx(un,I1,I2) + DpyDmy(un,I1,I2) ) + par.ufe(par.x(I1,I2,1),par.x(I1,I2,2),t);

  % vt(I1,I2) =  -( un(I1,I2).*Dzx(vn,I1,I2) + vn(I1,I2).*Dzy(vn,I1,I2) + Dzy(pn,I1,I2) ) ...
  %       + (nu*nuScaleFactor)*( DpxDmx(vn,I1,I2) + DpyDmy(vn,I1,I2) ) + par.vfe(par.x(I1,I2,1),par.x(I1,I2,2),t);  

  if( par.ad )
    %    --- 2nd order 2D artificial diffusion ---
    %    see cg/ins/src/insdt.bf
    ad21 = par.ad21;
    cd22 = par.ad22/(par.nd^2);
    adc = zeros(par.Ngx,par.Ngy);
    adc(I1,I2)=ad21 + cd22*( abs(ux(I1,I2))+abs(uy(I1,I2))  ...
                            +abs(vx(I1,I2))+abs(vy(I1,I2)) );

    ut(I1,I2) = ut(I1,I2) + adc(I1,I2)*(un(I1+1,I2  )-4.*un(I1,I2)+un(I1-1,I2) ... 
                                       +un(I1  ,I2+1)             +un(I1,I2-1));

    vt(I1,I2) = vt(I1,I2) + adc(I1,I2)*(vn(I1+1,I2  )-4.*vn(I1,I2)+vn(I1-1,I2) ... 
                                       +vn(I1  ,I2+1)             +vn(I1,I2-1));
  end 

 par.cpuGetUt = par.cpuGetUt + cputime - cpu0;
return
end


