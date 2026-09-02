%
% Compute the RHS to the momentum equations **NEW VERSION**
%
function [ ut,vt,par,uLap,vLap ] = getUt( t,un,vn,pn, gf,cur, nuScaleFactor,par )

  cpu0 = cputime;

  dx = par.dx;
  dy = par.dy;
  dr = par.dr(1);
  ds = par.dr(2);  
  nu = par.nu;

  if( ~strcmp(par.ms,'none') ) manufacturedSolution=1; else manufacturedSolution=0; end

  % declare operators 
  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y
  
  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy
  
  % derivatives of entries in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r to second order 
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s

  Drr2 = @(u,I1,I2) ( u(I1+1,I2) -2*u(I1,I2) +u(I1-1,I2) )/(dr^2);                             % u.rr to second order
  Dss2 = @(u,I1,I2) ( u(I1,I2+1) -2*u(I1,I2) +u(I1,I2-1) )/(ds^2);                             % u.ss
  Drs2 = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) + u(I1-1,I2-1) )/(4*dr*ds);   % u.rs

  % Dx2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,1)*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,1)*Ds2(u,I1,I2);  % u.x to order 2 
  % Dy2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,2)*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,2)*Ds2(u,I1,I2);  % u.y to order 2 


  [I1,I2] = getIndex( par.gid );

  ux = zeros(par.Ngx,par.Ngy);
  uy = zeros(par.Ngx,par.Ngy);
  vx = zeros(par.Ngx,par.Ngy);
  vy = zeros(par.Ngx,par.Ngy);

  uLap = zeros(par.Ngx,par.Ngy);
  uLap = zeros(par.Ngx,par.Ngy);

  % precompute forcing at all points (This is much faster)
  uf = zeros(par.Ngx,par.Ngy);
  vf = zeros(par.Ngx,par.Ngy);
  if( manufacturedSolution )
    uf(I1,I2) = par.ufe(gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),t);
    vf(I1,I2) = par.vfe(gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),t);

    % FIX ME FOR MOVING GRIDS 
  end


  % Note: moving grids are always treated as non-Cartesian
  if( par.isCartesian && par.gridMotion==par.noMotion )

    % --- Cartesian ---

    ux(I1,I2) = Dzx(un,I1,I2);
    uy(I1,I2) = Dzy(un,I1,I2);
    vx(I1,I2) = Dzx(vn,I1,I2);
    vy(I1,I2) = Dzy(vn,I1,I2);

    uLap(I1,I2) = DpxDmx(un,I1,I2) + DpyDmy(un,I1,I2);
    vLap(I1,I2) = DpxDmx(vn,I1,I2) + DpyDmy(vn,I1,I2);

    if( nuScaleFactor~=0 )
      ut(I1,I2) = -( un(I1,I2).*ux(I1,I2) + vn(I1,I2).*uy(I1,I2) + Dzx(pn,I1,I2) ) + nu*uLap(I1,I2) + uf(I1,I2);
      vt(I1,I2) = -( un(I1,I2).*vx(I1,I2) + vn(I1,I2).*vy(I1,I2) + Dzy(pn,I1,I2) ) + nu*vLap(I1,I2) + vf(I1,I2);

    else
      % leave off viscous terms 
      ut(I1,I2) = -( un(I1,I2).*ux(I1,I2) + vn(I1,I2).*uy(I1,I2) + Dzx(pn,I1,I2) ) + uf(I1,I2);
      vt(I1,I2) = -( un(I1,I2).*vx(I1,I2) + vn(I1,I2).*vy(I1,I2) + Dzy(pn,I1,I2) ) + vf(I1,I2);

    end

  else

    % --- curvilinear ---

    for( i1=I1 )
    for( i2=I2 ) 
      rx = gf{cur}.rx(i1,i2,1,1);
      ry = gf{cur}.rx(i1,i2,1,2);
      sx = gf{cur}.rx(i1,i2,2,1);
      sy = gf{cur}.rx(i1,i2,2,2);

      rxr = DJzr(gf{cur}.rx,i1,i2,1,1);
      ryr = DJzr(gf{cur}.rx,i1,i2,1,2);
      sxr = DJzr(gf{cur}.rx,i1,i2,2,1);
      syr = DJzr(gf{cur}.rx,i1,i2,2,2);

      rxs = DJzs(gf{cur}.rx,i1,i2,1,1);
      rys = DJzs(gf{cur}.rx,i1,i2,1,2);
      sxs = DJzs(gf{cur}.rx,i1,i2,2,1);
      sys = DJzs(gf{cur}.rx,i1,i2,2,2);        

      rxx = rx*rxr + sx*rxs;
      ryy = ry*ryr + sy*rys;
      sxx = rx*sxr + sx*sxs;
      syy = ry*syr + sy*sys;

      ur = Dr2(un,i1,i2); us = Ds2(un,i1,i2);
      vr = Dr2(vn,i1,i2); vs = Ds2(vn,i1,i2);
      pr = Dr2(pn,i1,i2); ps = Ds2(pn,i1,i2);

      urr = Drr2(un,i1,i2); urs = Drs2(un,i1,i2); uss = Dss2(un,i1,i2);
      vrr = Drr2(vn,i1,i2); vrs = Drs2(vn,i1,i2); vss = Dss2(vn,i1,i2);

      ux(i1,i2) = rx*ur + sx*us;  % save for use below in AD
      vx(i1,i2) = rx*vr + sx*vs;
      px        = rx*pr + sx*ps;

      uy(i1,i2) = ry*ur + sy*us;
      vy(i1,i2) = ry*vr + sy*vs;
      py        = ry*pr + sy*ps;

      uLap(i1,i2) = (rx^2+ry^2)*urr + 2*(rx*sx+ry*sy)*urs + (sx^2+sy^2)*uss + (rxx+ryy)*ur + (sxx+syy)*us;
      vLap(i1,i2) = (rx^2+ry^2)*vrr + 2*(rx*sx+ry*sy)*vrs + (sx^2+sy^2)*vss + (rxx+ryy)*vr + (sxx+syy)*vs;

      % grid velocity 
      gv(1) = gf{cur}.gv(i1,i2,1); 
      gv(2) = gf{cur}.gv(i1,i2,2); 

      % fprintf('getUt: gv=[%g,%g]\n',gv(1),gv(2));
      % pause

      % ut(i1,i2) = -( un(i1,i2)*ux(i1,i2) + vn(i1,i2)*uy(i1,i2) + px ) + nu*uLap + par.ufe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
      if( nuScaleFactor~=0 )
        ut(i1,i2) = -( (un(i1,i2)-gv(1))*ux(i1,i2) + (vn(i1,i2)-gv(2))*uy(i1,i2) + px ) + nu*uLap(i1,i2) + uf(i1,i2);
        vt(i1,i2) = -( (un(i1,i2)-gv(1))*vx(i1,i2) + (vn(i1,i2)-gv(2))*vy(i1,i2) + py ) + nu*vLap(i1,i2) + vf(i1,i2);
      else
        % leave off viscous terms 
        ut(i1,i2) = -( (un(i1,i2)-gv(1))*ux(i1,i2) + (vn(i1,i2)-gv(2))*uy(i1,i2) + px ) + uf(i1,i2);
        vt(i1,i2) = -( (un(i1,i2)-gv(1))*vx(i1,i2) + (vn(i1,i2)-gv(2))*vy(i1,i2) + py ) + vf(i1,i2);
      end


    end
    end   


  end

  if( par.ad )
    %    --- 2nd order 2D artificial diffusion ---
    %    see cg/ins/src/insdt.bf
    ad21 = par.ad21;
    cd22 = par.ad22/(par.nd^2);
    adc = zeros(par.Ngx,par.Ngy);
    adc(I1,I2)=ad21 + cd22*( abs(ux(I1,I2))+abs(uy(I1,I2))  ...
                            +abs(vx(I1,I2))+abs(vy(I1,I2)) );

    ut(I1,I2) = ut(I1,I2) + adc(I1,I2).*(un(I1+1,I2  )-4.*un(I1,I2)+un(I1-1,I2) ... 
                                        +un(I1  ,I2+1)             +un(I1,I2-1));

    vt(I1,I2) = vt(I1,I2) + adc(I1,I2).*(vn(I1+1,I2  )-4.*vn(I1,I2)+vn(I1-1,I2) ... 
                                        +vn(I1  ,I2+1)             +vn(I1,I2-1));
  end 

 par.cpuGetUt = par.cpuGetUt + cputime - cpu0;
return
end


