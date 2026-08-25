%
%  Apply boundary conditions
% 
function [u,v,par] = applyBoundaryConditions( u,v,t,par )

  cpu0 = cputime;

  % fprintf('Entering applyBC u=[%d,%d] v=[%d,%d]\n',...
  %        size(u,1),size(u,2), size(v,1),size(v,2) );

  Ngx = par.Ngx;
  Ngy = par.Ngy;

  dx = par.dx;
  dy = par.dy;
  dr = par.dr(1);
  ds = par.dr(2);  


  % declare operators 
  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy

  DzxDzy = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) +u(I1-1,I2-1) )*(1./(4.*dx*dy)); % u.xy 

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s

  Drr2 = @(u,I1,I2) ( u(I1+1,I2) -2*u(I1,I2) +u(I1-1,I2) )/(dr^2);                             % u.rr to second order
  Dss2 = @(u,I1,I2) ( u(I1,I2+1) -2*u(I1,I2) +u(I1,I2-1) )/(ds^2);                             % u.ss
  Drs2 = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) + u(I1-1,I2-1) )/(4*dr*ds);   % u.rs
  
  Dx2 = @(u,I1,I2) par.rx(I1,I2,1,1).*Dr2(u,I1,I2) + par.rx(I1,I2,2,1).*Ds2(u,I1,I2);  % u.x to order 2 
  Dy2 = @(u,I1,I2) par.rx(I1,I2,1,2).*Dr2(u,I1,I2) + par.rx(I1,I2,2,2).*Ds2(u,I1,I2);  % u.y to order 2   

  % Define extrapolations: (is1=+1/-1 and is2=+1/-1 defines the direction ("shift") of extrapolation)
  extrap2 = @(u,I1,I2,is1,is2) (2.*u(I1+is1,I2+is2) -    u(I1+2*is1,I2+2*is2)                       );   % 2nd-order extrapolation
  extrap3 = @(u,I1,I2,is1,is2) (3.*u(I1+is1,I2+is2) - 3.*u(I1+2*is1,I2+2*is2) + u(I1+3*is1,I2+3*is2));   % 3rd-order extrapolation

  if( ~strcmp(par.ms,'none') ) manufacturedSolution=1; else manufacturedSolution=0; end

  % --- STAGE 1: Dirichlet Type boundary conditions ---
  for side=1:2
    for axis=1:2
      
      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);
      mbc = side+2*(axis-1);  % pointer into gu and gv arrays
      [I1b,I2b]=getBoundaryIndex(side,axis,par);

      if( par.bc(side,axis)==par.dirichlet || par.bc(side,axis)==par.noSlipWall || par.bc(side,axis)==par.inflow )

        u(I1b,I2b)=par.gu{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
        v(I1b,I2b)=par.gv{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);

      elseif( par.bc(side,axis)==par.slipWall )
        % set n.uv = nv.gv 
        %  Set the normal component using: 
        %    uv = uv - nv.uv + nv.gv 

        % Start by extrapolating ghost IS THIS NEEDED ??
        if( 1==0 )
          I1g=I1b-is1; I2g=I2b-is2; % ghost points 
          u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
          v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2);  
        end

        n1 = -is1; n2 = -is2; % outward normal
        nDotU = n1*u(I1b,I2b) + ...
                n2*v(I1b,I2b);
        if( manufacturedSolution )
          nDotU = nDotU - n1*par.ue(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ...
                        - n2*par.ve(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
          % nDotU = nDotU - ( n1*( par.gu{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) + ...
          %                   n2*( par.gv{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) );          
        end
        u(I1b,I2b)= u(I1b,I2b) - nDotU*n1;
        v(I1b,I2b)= v(I1b,I2b) - nDotU*n2;


        % nDotU = zeros(par.Ngx,par.Ngy);
        % nDotU(I1b,I2b) = n1*u(I1b,I2b) + ...
        %                  n2*v(I1b,I2b);
        % if( ~strcmp(par.ms,'none') )
        %   nDotU(I1b,I2b) = nDotU(I1b,I2b) - ( n1*( par.gu{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) + ...
        %                                       n2*( par.gv{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) );          
        % end
        % u(I1b,I2b)= u(I1b,I2b) - nDotU(I1b,I2b)*n1;
        % v(I1b,I2b)= v(I1b,I2b) - nDotU(I1b,I2b)*n2;

      elseif( par.bc(side,axis)==par.pressureInflow )
        %
        %     p = given 
        %     tv.uv = tv.gv 
        %     set div(u)=0 
        %     Extrap( tv.uv )
        % 
        %  Set the tangential component to zero: 
        %    uv = (nv.uv) nv  

        n1 = -is1; n2 = -is2; % outward normal
        nDotU = n1*u(I1b,I2b) + ...
                n2*v(I1b,I2b);
        u(I1b,I2b)= nDotU*n1;
        v(I1b,I2b)= nDotU*n2;
        if( manufacturedSolution )

          if( axis ==1 )
            v(I1b,I2b)=par.ve(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
          else
            u(I1b,I2b)=par.ue(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
          end
          % nDotU = nDotU - ( n1*( par.gu{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) + ...
          %                   n2*( par.gv{mbc}(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t) ) );          
        end        


      elseif( par.bc(side,axis)==par.periodic || ...
              par.bc(side,axis)==par.outflow )
        % done below
      else
        fprintf('applyBoundaryConditions:ERROR: unknown BC=%d\n',par.bc(side,axis)); pause;
      end
    end
  end 

  % --- STAGE 2: Neumann type boundry conditons ----

  % -- First extrapolate all ghost points ---
  if( 1==1 )
    for side=1:2
      for axis=1:2
        % isv(1:2), is1, is2 : index shifts
        isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);

        [I1b,I2b]=getBoundaryIndex(side,axis,par);
        I1g=I1b-is1; I2g=I2b-is2; % ghost points 
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2); 
      end
    end
  end

  for side=1:2
    for axis=1:2
      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);

      [I1b,I2b]=getBoundaryIndex(side,axis,par);
      I1g=I1b-is1; I2g=I2b-is2; % ghost points 
      % -- set ghost ---
      if( par.bc(side,axis)==par.outflow )
        % -- extrapolate ghost to order 3 ? --
        % Also possible: apply a Neumann BC 
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2);   

      elseif( par.bc(side,axis)==par.slipWall ) 
        % ( tv.uv).n = 0  
        n1 = -is1;
        n2 = -is2;
     
        if( axis==1 )
          % v.x = given
          if( manufacturedSolution )
            gvx = n1*par.vex(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
          else
            gvx = 0.; 
          end 
          v(I1g,I2g) = v(I1b+is1,I2b+is2) + (2.*dx)*( gvx );
        else
          % u.y = given 
          if( manufacturedSolution )
            guy = n2*par.uey(par.x(I1b,I2b,1),par.x(I1b,I2b,2),t);
          else
            guy = 0.; % fix me for MS
          end 
          u(I1g,I2g) = u(I1b+is1,I2b+is2) + (2.*dy)*( guy );
        end 

      elseif( par.bc(side,axis)==par.pressureInflow )
        % Extrap tv.uv 
        % Extrap Both: (div(u) will be set below)
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2);        

      elseif( par.bc(side,axis)==par.periodic )
        % done below 

      elseif( 1==1 )
        % -- extrapolate ghost --
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2);

      else
        mbc = side+2*(axis-1);  % pointer into gu and gv arrays
        % Fot testing set the exact values at the ghost
        u(I1g,I2g)=par.gu{mbc}(par.x(I1g,I2g,1),par.x(I1g,I2g,2),t);
        v(I1g,I2g)=par.gv{mbc}(par.x(I1g,I2g,1),par.x(I1g,I2g,2),t);
      end
    end
  end

  % --- Assign div(u) =0 along walls ---
  %   Note: at a corner with two noSlipWalls, the same divergence condition
  %   is applied twice to get the two ghost -- fix me --
  for side=1:2
    for axis=1:2
      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);
      [I1b,I2b]=getBoundaryIndex(side,axis,par);
      I1g=I1b-is1; I2g=I2b-is2; % ghost points 

      % -- Use divergence to set ghost ---
      %   **** WATCH OUT FOR CORNERS WHERE div(u)=0 is applied twice -- FIX ME 
      if( par.bc(side,axis)==par.noSlipWall     || ...
          par.bc(side,axis)==par.slipWall       || ...
          par.bc(side,axis)==par.inflow         || ...
          par.bc(side,axis)==par.pressureInflow || ...
          par.bc(side,axis)==par.outflow )
        if( par.isCartesian )

          % --- Cartesian grid ---
          if( axis==1 )
            % u.x = - v-y 
            u(I1g,I2g) = u(I1b+is1,I2b+is2) + (2.*dx*is1)*Dzy(v,I1b,I2b); 
          else
            % v.y = - u.x 
            v(I1g,I2g) = v(I1b+is1,I2b+is2) + (2.*dy*is2)*Dzx(u,I1b,I2b); 
          end
        else
          % ---- Curvilinear ---
          %  div(uv) = 0 SETS THE NORMAL COMPONENT OF uv
          % 
          %   div(uv) = ar*u.r + br*v.r + as*u.s + bs*v.s 
          %    ar = rx, as=sx
          %    br = ry, bs=sy
          % is = 1-2*(side-1);

          % n1 = -is*par.rx(I1b,I2b,axis,1); % outward normal is (n1,n2)
          % n2 = -is*par.rx(I1b,I2b,axis,2); 
          % rxNorm = sqrt( n1.^2 + n2.^2 ); 
          % n1 = n1./rxNorm;
          % n2 = n2./rxNorm;
          [n1,n2,rxNorm] = getBoundaryNormal( side,axis,I1b,I2b,par );

          nDotUg = n1.*u(I1g    ,I2g    ) + n2.*v(I1g    ,I2g    ); % nv.uv on ghost to start
          nDotUp = n1.*u(I1b+is1,I2b+is2) + n2.*v(I1b+is1,I2b+is2); % nv.uv on first line in
          if( axis==1)
            nDotUg = nDotUp - (2*dr)*(par.rx(I1b,I2b,2,1).*Ds2(u,I1b,I2b) + par.rx(I1b,I2b,2,2).*Ds2(v,I1b,I2b))./rxNorm - nDotUg;
          else
            nDotUg = nDotUp - (2*ds)*(par.rx(I1b,I2b,1,1).*Dr2(u,I1b,I2b) + par.rx(I1b,I2b,1,2).*Dr2(v,I1b,I2b))./rxNorm - nDotUg;
          end 
          % Set normal component of the ghost value
          u(I1g,I2g) = u(I1g,I2g) + nDotUg.*n1; 
          v(I1g,I2g) = v(I1g,I2g) + nDotUg.*n2; 

          if( par.idebug>3 )
            res = Dx2(u,I1b,I2b) + Dy2(v,I1b,I2b);
            fprintf('applyBC: (side,axis)=(%d,%d) after set div(uv)=0 : max(abs(div))=%9.2e\n',side,axis,max(abs(res)));
            pause
          end
   
        end

      end
    end
  end

  % --- extrapolate corners along the diagonal ---
  for( side1=1:2 )
    for( side2=1:2 )
      is1 = 1-2*(side1-1); % is1 = +1 or -1 
      is2 = 1-2*(side2-1); % is2 = +1 or -1 
      ix=par.gid(side1,1)-is1; iy=par.gid(side2,2)-is2; % corner ghost point
      u(ix,iy) = extrap3(u,ix,iy,is1,is2);
      v(ix,iy) = extrap3(v,ix,iy,is1,is2);
    end
  end 

  % --- periodic BC's --- 
  % fprintf('applyBC: par.bc=[%d,%d,%d,%d]\n',par.bc(1,1),par.bc(2,1),par.bc(1,2),par.bc(2,2));
  for axis=1:2
    if( par.bc(1,axis)==par.periodic || par.bc(2,axis)==par.periodic )
      if( ~(par.bc(1,axis)==par.periodic && par.bc(2,axis)==par.periodic ) )
        fprintf('applyBC:ERROR: bc is periodic on one side but not the other for axis=%d\n',axis); pause;
      end
      % fprintf('applyBC: periodic...\n');

      isv(1)=0; isv(2)=0; isv(axis)=1;  is1=isv(1); is2=isv(2);

      [I1a,I2a]=getBoundaryIndex(1,axis,par); if( axis==1 ) I2a=1:Ngy; else I1a=1:Ngx; end 
      [I1b,I2b]=getBoundaryIndex(2,axis,par); if( axis==1 ) I2b=1:Ngy; else I1b=1:Ngx; end 

      % fprintf('apply periodic BC: axis=%d\n',axis); 
      u(I1b    ,I2b    )=u(I1a    ,I2a    ); v(I1b    ,I2b    )=v(I1a    ,I2a    );  % right/top side = left/bot
      u(I1b+is1,I2b+is2)=u(I1a+is1,I2a+is2); v(I1b+is1,I2b+is2)=v(I1a+is1,I2a+is2);  % right/top ghost = 
      u(I1a-is1,I2a-is2)=u(I1b-is1,I2b-is2); v(I1a-is1,I2a-is2)=v(I1b-is1,I2b-is2);  % left/bottom ghost = 
     
    end
  end 

  par.cpuBC = par.cpuBC + cputime - cpu0;

end



