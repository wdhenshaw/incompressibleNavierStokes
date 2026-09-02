%
%  Apply boundary conditions
% 
function [u,v,par] = applyBoundaryConditions( u,v,t, gf,cur, par )

  cpu0 = cputime;

  % fprintf('Entering applyBC u=[%d,%d] v=[%d,%d]\n',...
  %        size(u,1),size(u,2), size(v,1),size(v,2) );

  nu = par.nu;
  mu = par.mu;

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
  
  Dx2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,1).*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,1).*Ds2(u,I1,I2);  % u.x to order 2 
  Dy2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,2).*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,2).*Ds2(u,I1,I2);  % u.y to order 2   

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

        u(I1b,I2b)=par.gu{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
        v(I1b,I2b)=par.gv{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);

      elseif( par.bc(side,axis)==par.slipWall )
        % done below now 
        
        %   % set n.uv = nv.gv 
        %   %  Set the normal component using: 
        %   %    uv = uv - nv.uv + nv.gv 

        %   if( par.isCartesian )
        %     n1 = -is1; n2 = -is2; % outward normal
        %   else
        %     [n1,n2,rxNorm] = getBoundaryNormal( side,axis,I1b,I2b, gf,cur, par );
        %   end
        %   nDotU = n1*u(I1b,I2b) + ...
        %           n2*v(I1b,I2b);
        %   if( manufacturedSolution )
        %     nDotU = nDotU - n1*par.ue(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) ...
        %                   - n2*par.ve(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
        %     % nDotU = nDotU - ( n1*( par.gu{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) ) + ...
        %     %                   n2*( par.gv{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) ) );          
        %   end
        %   u(I1b,I2b)= u(I1b,I2b) - nDotU*n1;
        %   v(I1b,I2b)= v(I1b,I2b) - nDotU*n2;


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
            v(I1b,I2b)=par.ve(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
          else
            u(I1b,I2b)=par.ue(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
          end
          % nDotU = nDotU - ( n1*( par.gu{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) ) + ...
          %                   n2*( par.gv{mbc}(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) ) );          
        end        


      elseif( par.bc(side,axis)==par.periodic || ...
              par.bc(side,axis)==par.outflow  || ...
              par.bc(side,axis)==par.traction )
        % done below
      else
        fprintf('applyBoundaryConditions:ERROR: unknown BC=%d\n',par.bc(side,axis)); pause;
      end
    end
  end 

  % --- STAGE 2: Neumann type boundary conditons ----

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

  % --- APPLY periodic BC's I --- 
  %  (these are done again below )
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

      elseif( par.bc(side,axis)==par.slipWall && par.isCartesian && par.combinedImplicitSolverNeeded==0 ) 
        % --- Slip wall ---
        % ( tv.uv).n = 0  
        % NOTE: curvlinear case is doen below with traction BC

        n1 = -is1;
        n2 = -is2;
     
        if( axis==1 )
          % v.x = given
          if( manufacturedSolution )
            gvx = n1*par.vex(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
          else
            gvx = 0.; 
          end 
          v(I1g,I2g) = v(I1b+is1,I2b+is2) + (2.*dx)*( gvx );
        else
          % u.y = given 
          if( manufacturedSolution )
            guy = n2*par.uey(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
          else
            guy = 0.;
          end 
          u(I1g,I2g) = u(I1b+is1,I2b+is2) + (2.*dy)*( guy );
        end 

      elseif( par.bc(side,axis)==par.pressureInflow )
        % Extrap tv.uv 
        % Extrap Both: (div(u) will be set below)
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2); 

      elseif( par.bc(side,axis)==par.traction || par.bc(side,axis)==par.slipWall ) 
        %  ---- Traction BC or Slip Wall ----
        %   
        %   tv^T sigmav nv = 0 
        % 
        if( par.isCartesian && par.combinedImplicitSolverNeeded==0 )
          % mu*( u.y + v.x ) = 0 
          if( manufacturedSolution )
            g = par.uey(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t) + ...
                par.vex(gf{cur}.x(I1b,I2b,1),gf{cur}.x(I1b,I2b,2),t);
          else
            g = 0;
          end
          if( axis==1 )
            % v.x = -u.y + g 
            v(I1g,I2g) = v(I1b+is1,I2b+is2) - is1*(2.*dx)*( -Dzy(u,I1b,I2b) + g  );
          else
            % u.y = -v.x + g 
            u(I1g,I2g) = u(I1b+is1,I2b+is2) - is2*(2.*dy)*( -Dzx(v,I1b,I2b) + g  );
          end

        else
          % TRACTION CURVLINEAR:
          % IMPOSE
          %    div(uv)=0 and tv tauv nv = 0
          % Write these two equations as
          %  (div(uv)) nv = (I-v nv^T) tauv nv=0 
 
          % fprintf('applyBC: ERROR: finish traction BC for curvilinear\n');
          % error();   

          
          f = zeros(2,1); b=zeros(2,1); A=zeros(2,2); gv=zeros(2,1);
          for i2=I2b
            for i1=I1b

              % Use Newton's method to solve the equations f(1:2)=0 for the unknown ghost values u(-1), v(-1)
              if( par.isCartesian )
                n1=-is1;
                n2=-is2;
                rx=1; ry=0; sx=0; sy=1; 
              else
                [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );

                rx = gf{cur}.rx(i1,i2,1,1); sx = gf{cur}.rx(i1,i2,2,1);
                ry = gf{cur}.rx(i1,i2,1,2); sy = gf{cur}.rx(i1,i2,2,2);
              end 
              px = -is1*rx/(2*dr) -is2*sx/(2*ds); % coeff of ghost in D_x = rx*D0r + sx*D0s 
              py = -is1*ry/(2*dr) -is2*sy/(2*ds); % coeff of ghost in D_y  

              % A(1,1) = coeff of u(ghost) in f(1) 
              % A(1,2) = coeff of v(ghost) in f(1) 
              % A(2,1) = coeff of u(ghost) in f(2) 
              % A(2,2) = coeff of v(ghost) in f(2) 
              A(1,1) = px*n1 + (1-n1^2)*mu*( 2*px*n1 + py*n2 )    -n1*n2*mu*(           py*n1 ); 
              A(1,2) = py*n1 + (1-n1^2)*mu*(           px*n2 )    -n1*n2*mu*( 2*py*n2 + px*n1 ); 
              A(2,1) = px*n2  -   n1*n2*mu*( 2*px*n1 + py*n2 ) +(1-n2^2)*mu*(           py*n1 ); 
              A(2,2) = py*n2  -   n1*n2*mu*(           px*n2 ) +(1-n2^2)*mu*( 2*py*n2 + px*n1 ); 

              if( manufacturedSolution )
                ux = par.uex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                uy = par.uey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                vx = par.vex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                vy = par.vey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);    

                tauvnv1 =  mu*( 2*ux*n1 + (uy+vx)*n2 );
                tauvnv2 =  mu*( 2*vy*n2 + (uy+vx)*n1 );   
                gv(1) = ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )    -n1*n2*( tauvnv2 );
                gv(2) = ( ux+vy )*n2 +   -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 );                                         
              end                                                   

              ux =  rx*Dr2(u,i1,i2) + sx*Ds2(u,i1,i2);
              vx =  rx*Dr2(v,i1,i2) + sx*Ds2(v,i1,i2);
              uy =  ry*Dr2(u,i1,i2) + sy*Ds2(u,i1,i2);
              vy =  ry*Dr2(v,i1,i2) + sy*Ds2(v,i1,i2);

              % tauv*nv : 
              %   tau nv = mu [ 2*ux uy+vx ][ n1 ]
              %               [ uy+vx 2*vy ]  n2 ]           
              tauvnv1 =  mu*( 2*ux*n1 + (uy+vx)*n2 );
              tauvnv2 =  mu*( 2*vy*n2 + (uy+vx)*n1 );
              %   I - nv nv^T = [ 1-n1^2 -n1*n2 ]
              %                 [ -n1*n2 1-n2^2 ]
              f(1) = ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )    -n1*n2*( tauvnv2 ) - gv(1);
              f(2) = ( ux+vy )*n2 +   -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 ) - gv(2);
    


              b(1) = u(i1-is1,i2-is2); % solve for these ghost point 
              b(2) = v(i1-is1,i2-is2);

              b = b - A\f;

              u(i1-is1,i2-is2) = b(1);
              v(i1-is1,i2-is2) = b(2);

              if( 1==0 )
                % check residual and errors
                ux =  rx*Dr2(u,i1,i2) + sx*Ds2(u,i1,i2);
                vx =  rx*Dr2(v,i1,i2) + sx*Ds2(v,i1,i2);
                uy =  ry*Dr2(u,i1,i2) + sy*Ds2(u,i1,i2);
                vy =  ry*Dr2(v,i1,i2) + sy*Ds2(v,i1,i2);

                % tauv*nv : 
                tauvnv1 =  mu*( 2*ux*n1 + (uy+vx)*n2 );
                tauvnv2 =  mu*( 2*vy*n2 + (uy+vx)*n1 );
                f(1) = ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )   -n1*n2*( tauvnv2 ) - gv(1);
                f(2) = ( ux+vy )*n2 +  -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 ) - gv(2);    
                
                uErr = abs( u(i1-is1,i2-is2) - par.ue(gf{cur}.x(i1-is1,i2-is2,1),gf{cur}.x(i1-is1,i2-is2,2),t));
                vErr = abs( v(i1-is1,i2-is2) - par.ve(gf{cur}.x(i1-is1,i2-is2,1),gf{cur}.x(i1-is1,i2-is2,2),t));

                fprintf('applyBC:tractionBC: [side,axis]=[%d,%d] (i1,i2)=%3d,%3d)  f=(%9.2e,%9.2e) err=(%9.2e,%9.2e)\n',side,axis,i1,i2,f(1),f(2),uErr,vErr);  
                % pause;

              end


            end % end for i1
          end % end for i2 


        end

      elseif( par.bc(side,axis)==par.periodic )
        % done below 

      elseif( 1==1 )
        % -- extrapolate ghost --
        u(I1g,I2g) = extrap3(u,I1g,I2g,is1,is2);
        v(I1g,I2g) = extrap3(v,I1g,I2g,is1,is2);

      else
        mbc = side+2*(axis-1);  % pointer into gu and gv arrays
        % Fot testing set the exact values at the ghost
        u(I1g,I2g)=par.gu{mbc}(gf{cur}.x(I1g,I2g,1),gf{cur}.x(I1g,I2g,2),t);
        v(I1g,I2g)=par.gv{mbc}(gf{cur}.x(I1g,I2g,1),gf{cur}.x(I1g,I2g,2),t);
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
          par.bc(side,axis)==par.traction       || ...
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

          % n1 = -is*gf{cur}.rx(I1b,I2b,axis,1); % outward normal is (n1,n2)
          % n2 = -is*gf{cur}.rx(I1b,I2b,axis,2); 
          % rxNorm = sqrt( n1.^2 + n2.^2 ); 
          % n1 = n1./rxNorm;
          % n2 = n2./rxNorm;
          [n1,n2,rxNorm] = getBoundaryNormal( side,axis,I1b,I2b, gf,cur, par );

          nDotUg = n1.*u(I1g    ,I2g    ) + n2.*v(I1g    ,I2g    ); % nv.uv on ghost to start
          nDotUp = n1.*u(I1b+is1,I2b+is2) + n2.*v(I1b+is1,I2b+is2); % nv.uv on first line in
          if( axis==1)
            nDotUg = nDotUp - (2*dr)*(gf{cur}.rx(I1b,I2b,2,1).*Ds2(u,I1b,I2b) + gf{cur}.rx(I1b,I2b,2,2).*Ds2(v,I1b,I2b))./rxNorm - nDotUg;
          else
            nDotUg = nDotUp - (2*ds)*(gf{cur}.rx(I1b,I2b,1,1).*Dr2(u,I1b,I2b) + gf{cur}.rx(I1b,I2b,1,2).*Dr2(v,I1b,I2b))./rxNorm - nDotUg;
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

  % --- CHECK FOR CORNER CASES WE DO NOT TREAT YET ---
  for( side1=1:2 )
    for( side2=1:2 )
      if( par.bc(side1,1)==par.traction && ...
          par.bc(side2,2)==par.traction ) 
        fprintf('\n ***** applyBoundaryConditions:ERROR -- a traction-traction corner is not implemented yet ***\n\n');
        error('ERROR')
      end 
    end
  end


  par.cpuBC = par.cpuBC + cputime - cpu0;

end



