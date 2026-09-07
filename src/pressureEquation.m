%
%  Solve the Poisson equation for p 
%
% Input:
%  factorMatrix =1 : factor the matrix
%
function [p,par] = pressureEquation( t,u,v,dt, gf,cur, par )

 % fprintf('Entering pressureEquation u=[%d,%d] v=[%d,%d]\n',...
 %         size(u,1),size(u,2), size(v,1),size(v,2) );

  nu  = par.nu;
  mu  = par.mu;
  cdv = par.cdv;
  ms  = par.ms;
  rho = par.rho;

  % iax = par.iax;  ibx = par.ibx; 
  % iay = par.iay;  iby = par.iby; 

  % i1x = par.i1x;  i2x = par.i2x; 
  % i1y = par.i1y;  i2y = par.i2y; 

  iax = par.gid(1,1); ibx=par.gid(2,1);
  iay = par.gid(1,2); iby=par.gid(2,2);

  Ngx = par.Ngx;
  Ngy = par.Ngy;

  dx = par.dx;
  dy = par.dy;  
  dr = par.dr(1);
  ds = par.dr(2);  

  numGhost = par.numGhost;

  idebug = par.idebug;

  % --- declare difference operators ---
  % declare operators 
  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy

  DzxDzy = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) +u(I1-1,I2-1) )*(1./(4.*dx*dy)); % u.xy 

  % derivatives of entreis in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s

  Dx2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,1).*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,1).*Ds2(u,I1,I2);  % u.x to order 2 
  Dy2 = @(u,I1,I2) gf{cur}.rx(I1,I2,1,2).*Dr2(u,I1,I2) + gf{cur}.rx(I1,I2,2,2).*Ds2(u,I1,I2);  % u.y to order 2 



  % Define extrapolations: (is1=+1/-1 and is2=+1/-1 defines the direction ("shift") of extrapolation)
  extrap3 = @(u,I1,I2,is1,is2) (3.*u(I1+is1,I2+is2) - 3.*u(I1+2*is1,I2+2*is2) + u(I1+3*is1,I2+3*is2));   % 3rd-order extrapolation

  % convert (i1,i2) to equation number in the matrix:  
  eqn = @(i1,i2)  1 + i1-1 + Ngx*( i2-1 );

    
  % persistent NgSave L U ; % holds LU factors of the pressure matrix

  Ng = Ngx*Ngy; % total number of grid points

  % check if the Pressure matrix is singular
  isSingular = 1;
  for side=1:2
    for axis=1:2
      if( par.bc(side,axis)==par.dirichlet      || ...
          par.bc(side,axis)==par.traction       || ...
          par.bc(side,axis)==par.pressureInflow || ...
          par.bc(side,axis)==par.outflow )
        isSingular=false;
      end
    end
  end
  Ngs =Ng;
  if( isSingular==1 ) Ngs=Ng+1; end % add an extra equation for the singular case
  


  if( par.factorPressureMatrix==1 )  
    % ----- Form the pressure matrix ------

    par.factorPressureMatrix=0;

    % Allocate the sparse matrix, vector-solution and RHS

    if( par.plotOption>=0 && ( mod(floor(idebug/2),2)==1 || (par.gridMotion~=par.noMotion && t<=2*dt) )  )
      fprintf('*** FORM THE PRESSURE MATRIX : Ng=%d t=%12.4e cur=%d **\n',Ng,t,cur);
    end
    
    % NgSave=Ng; % save 

    rightNullValue=0;
    if( isSingular==1 )
      if( par.idebug>0 && t<=2*dt ) fprintf('pressureEquation: **MATRIX IS SINGULAR** add an extra equation.\n'); end
      rightNullValue=1; 
    end



    % -------------------------------
    % ---- OPTIMIZED FILL MATRIX ----
    % -------------------------------

    nzzEst = Ngs*5; % estimated number of non-zeros
    if( isSingular==1 )
      nzzEst=Ngs*6; 
    end
    ia=zeros(nzzEst,1); ja=zeros(nzzEst,1); aa=zeros(nzzEst,1);

    nzz=0; % counts non-zeros


    % ---- NEW OPTIMIZED WAY TO FILL THE MATRIX ----
    [I1,I2] = getIndexInterior( par.gid,par.pc,par );

    % fprintf('form pressure matrixImplicitTimeSteppingMatrix:  I1=[%d,%d] I2=[%d,%d] (interior eqn)\n',I1(1),I1(end), I2(1),I2(end));

    if( par.isCartesian )

      % --- Cartesian grid ---
      for( i1=I1 )
      for( i2=I2 )
        ie = eqn(i1,i2) ; % eqn number for pt (i1,i2) 
        setValue( ie,eqn(i1  ,i2-1),(               1/dy^2   ) );
        setValue( ie,eqn(i1-1,i2  ),(      1/dx^2            ) );
        setValue( ie,eqn(i1  ,i2  ),( -2*( 1/dx^2 + 1/dy^2 ) ) );
        setValue( ie,eqn(i1+1,i2  ),(      1/dx^2            ) );
        setValue( ie,eqn(i1  ,i2+1),(               1/dy^2   ) );
        if( isSingular==1 ) 
          setValue(ie,Ngs,1.);
        end 
      end
      end
      
    else
      % --- curvilinear ---
 

      for( i1=I1 )
      for( i2=I2 )
        ie = eqn(i1,i2) ; % eqn number for pt (i1,i2) 
        rx = gf{cur}.rx(i1,i2,1,1);
        ry = gf{cur}.rx(i1,i2,1,2);
        sx = gf{cur}.rx(i1,i2,2,1);
        sy = gf{cur}.rx(i1,i2,2,2);

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

        setValue( ie,eqn(i1-1,i2-1),(                            +c11  ) );
        setValue( ie,eqn(i1  ,i2-1),(            c02        -c01       ) );
        setValue( ie,eqn(i1+1,i2-1),(                            -c11  ) );
        setValue( ie,eqn(i1-1,i2  ),(      c20        -c10             ) );
        setValue( ie,eqn(i1  ,i2  ),( -2*( c20 + c02  )                ) );
        setValue( ie,eqn(i1+1,i2  ),(      c20        +c10             ) );
        setValue( ie,eqn(i1-1,i2+1),(                            -c11  ) );
        setValue( ie,eqn(i1  ,i2+1),(            c02       +c01        ) );
        setValue( ie,eqn(i1+1,i2+1),(                            +c11  ) );        
        if( isSingular==1 ) 
          setValue(ie,Ngs,1.);
        end 
      end
      end      
    end
    
    % add extra equation in the last row -- all ones except that last column
    if( isSingular==1 )
      ie=Ngs; % final equation
      for( i1=I1 )
      for( i2=I2 )
        setValue(ie,eqn(i1,i2),1.);
      end
      end
    end
     
    % ---- Boundary conditions for the matrix ----

    dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

    % cornerIsSet = zeros(2,2,par.nd); 

    for side=1:2
      for axis=1:2
        % isv(1:2), is1, is2 : index shifts
        isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
  
        [I1b,I2b]=getBoundaryIndex(side,axis,par);

        if( par.bc(side,axis)==par.dirichlet ||  par.bc(side,axis)==par.traction )
          % -- special case for Dirichlet or Traction BC (FIX ME for some other cases)
          J1b = I1b;
          J2b = I2b;
          if( par.bc(side,axis)==par.dirichlet && axis==2 ) % top or bottom   ** FIX ME**
            % skip Dirichlet-Dirichlet, D-T, T-T corners  
            i1a = par.gid(1,1);
            i1b = par.gid(2,1);
            if( par.bc(1,1)==par.dirichlet || par.bc(1,1)==par.pressureInflow || par.bc(1,1)==par.traction ) i1a=par.gid(1,1)+1; end
            if( par.bc(2,1)==par.dirichlet || par.bc(2,1)==par.pressureInflow || par.bc(2,1)==par.traction ) i1b=par.gid(2,1)-1; end
            J1b = i1a:i1b;
          end
          % fprintf(' side=%d axis=%d, I1b=[%d,%d] I2b=[%d,%d] J1b=[%d,%d] J2b=[%d,%d]\n',...
          %      side,axis,I1b(1),I1b(end), I2b(1), I2b(end), J1b(1),J1b(end), J2b(1), J2b(end) );
          % WARNING : Do not set the corner value more than once since Matlab ACCUMULATES repeated entries
          for i2=J2b
          for i1=J1b
            ie=eqn(i1,i2);   % boundary point
            setValue(ie,ie,1.); 
          end 
          end       
        end 


        for i2=I2b
          for i1=I1b

            if( par.bc(side,axis)==par.dirichlet || par.bc(side,axis)==par.traction )

              % Extrapolate ghost for now 
              % >>> We could use a CBC here 

              % Boundary points are now done above
              % % WARNING : Do not set the corner value more than once since Matlab ACCUMULATES repeated entries
              % ie=eqn(i1,i2);   % boundary point
              % setValue(ie,ie,1.);

              % Extrapolate ghost
              ie=eqn(i1-is1,i2-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(i1      ,i2      ),-3.); 
              setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
              setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 

            elseif( par.bc(side,axis)==par.pressureInflow )

              % % WARNING : Do not set the corner value more than once since Matlab ACCUMULATES repeated entries
              ie=eqn(i1,i2);   % boundary point
              setValue(ie,ie,1.);

              % Extrapolate ghost
              ie=eqn(i1-is1,i2-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(i1      ,i2      ),-3.); 
              setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
              setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 


            elseif( par.bc(side,axis)==par.noSlipWall || ...
                    par.bc(side,axis)==par.slipWall   || ...
                    par.bc(side,axis)==par.inflow )
              % Neumann BC :
                %    (+-)*Dz( p ) = RHS  : use outward normal 
              % nSign = 2*(side-1)-1; % sign of normal, -1 on left and +1 on right
              % fprintf('pressureEqn: fill-in a Neumann BC (side,axis)=(%d,%d) (i1,i2)=(%d,%d) (is1,is2)=(%d,%d) nSign=%g\n',side,axis,i1,i2,is1,is2,nSign); 

              ie=eqn(i1-is1,i2-is2); % ghost point

              if( par.isCartesian )
                setValue(ie,eqn(i1-is1,i2-is2), 1/(2.*dxv(axis))); % ghost point 
                setValue(ie,eqn(i1+is1,i2+is2),-1/(2.*dxv(axis))); % first interior point
              else
                % p.n = n1*p.x + n2*p.y 
                %     = (n1*rx + n2*ry) p.r + (n1*sx + n2*sy)*p.s = ...
                rx = gf{cur}.rx(i1,i2,1,1);
                ry = gf{cur}.rx(i1,i2,1,2);
                sx = gf{cur}.rx(i1,i2,2,1);
                sy = gf{cur}.rx(i1,i2,2,2);

                % ---- get outward normal (n1,n2) ----
                is = 1-2*(side-1);
                n1 = -is*gf{cur}.rx(i1,i2,axis,1); 
                n2 = -is*gf{cur}.rx(i1,i2,axis,2); 
                rxNorm = sqrt( n1.^2 + n2.^2 ); 
                n1 = n1./rxNorm;
                n2 = n2./rxNorm;                  
                % --- done get normal ---

                ar = n1*rx + n2*ry;
                as = n1*sx + n2*sy;
                setValue(ie,eqn(i1-1,i2  ),-ar/(2.*dr)); 
                setValue(ie,eqn(i1+1,i2  ), ar/(2.*dr)); 
                setValue(ie,eqn(i1  ,i2-1),-as/(2.*ds)); 
                setValue(ie,eqn(i1  ,i2+1), as/(2.*ds)); 

              end

            elseif( par.bc(side,axis)==par.outflow )
              % a0*p + a1*p.n = RHS 

              % This equation is put into the matrix at the ghost point: 
              ie=eqn(i1-is1,i2-is2);     % ghost point 
              a0 = par.outflowPressureCoeffp;
              setValue(ie,eqn(i1,i2),a0); % 

              ie=eqn(i1-is1,i2-is2); % ghost point
              a1 = par.outflowPressureCoeffpn;
              setValue(ie,eqn(i1-is1,i2-is2), a1/(2.*dxv(axis))); % ghost point 
              setValue(ie,eqn(i1+is1,i2+is2),-a1/(2.*dxv(axis))); % first interior point                

              
            elseif( par.bc(side,axis)==par.periodic )
              % done below now 
              % % Periodic: p(iax-1,.) = p(ibx-1,.) ...
              % ie=eqn(i1-is1,i2-is2); % ghost point
              % setValue(ie,ie                                            , 1.); 
              % setValue(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2),-1.); 


            else
              fprintf('pressure matrx: fill-BC: finish me...\n'); pause; 
            end
          end
        end % end for i2

        % ---- periodic boundary conditions ---
        %  include ghost points 
        % [I1b,I2b]=getBoundaryIndex(side,axis,par);
        if( par.bc(side,axis)==par.periodic )
          if( axis==1 )
            I1b=par.gid(side,axis);
            I2b=par.dim(1,2):par.dim(2,2);
          elseif( axis==2 )
            I1b=par.dim(1,1):par.dim(2,1);  
            I2b=par.gid(side,axis);
          end        

          for i2=I2b
            for i1=I1b
                % Periodic: p(iax-1,.) = p(ibx-1,.) ...
                ie=eqn(i1-is1,i2-is2); % ghost point
                setValue(ie,ie                                            , 1.); 
                setValue(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2),-1.);  
            end 
          end   
        end      
      end % end for axis
    end % end for side


    % extrapolate corners along the diagonal 
    % par.gid(side,axis) : grid index range
    % par.gid(1,1)=iax; par.gid(2,1)=ibx; par.gid(1,2)=iay; par.gid(2,2)=iby; 
    for( side1=0:1 )
      for( side2=0:1 )
        if( par.bc(side1+1,1)~=par.periodic &&  par.bc(side2+1,2)~=par.periodic )
          is1 = 1-2*side1;
          is2 = 1-2*side2;
          i1=par.gid(side1+1,1)-is1; i2=par.gid(side2+1,2)-is2; % corner ghost point
          ie=eqn(i1,i2); 
          setValue(ie,ie                    , 1.);
          setValue(ie,eqn(i1+  is1,i2+  is2),-3.); 
          setValue(ie,eqn(i1+2*is1,i2+2*is2), 3.);  
          setValue(ie,eqn(i1+3*is1,i2+3*is2),-1.); 
        end
      end
    end 

    % --- CHECK FOR CORNER CASES WE DO NOT TREAT YET ---
    for( side1=1:2 )
      for( side2=1:2 )
        if( par.bc(side1,1)==par.traction && ...
            par.bc(side2,2)==par.traction ) 
          fprintf('\n ***** pressureEquation:ERROR -- a traction-traction corner is not implemented yet ***\n\n');
          error('ERROR')
        end 
      end
    end    

    if( par.idebug>0 && t<=dt )
      fprintf('Optimized fill of pressure matrix: nzzEst=%d, nzz=%d\n',nzzEst,nzz);
    end
    A = sparse(ia(1:nzz), ja(1:nzz), aa(1:nzz), Ngs, Ngs); % Creates the sparse matrix 

    

    % [L,U]=lu(A);   % factor the matrix: A=LU
    par.dA = decomposition(A);   % factor the matrix

    if( mod(floor(idebug/2),2)==1 )
      fprintf('+++ cond(A)=%10.2e\n',condest(A));
      figure(3); 
      spy(A);
      pause
    end

  else
    if( mod(floor(idebug/4),2)==1 )
      fprintf('*** DO NOT FORM THE PRESSURE MATRIX : Ng=%d **\n',Ng);
    end
  end % end build matrix

  cpu0 = cputime;

  % -------------------------------------
  % ---------- Assign the RHS  ----------
  % -------------------------------------
  rhs = zeros(Ngs,1);     % right-hand-side

  % I = iax:ibx; J=iay:iby;  % interior and boundary points
  [I,J] = getIndex( par.gid );

  pf = zeros(Ngx,Ngy);
  % fprintf('Ngx=%d, Ngy=%d, [iax,ibx]=[%d,%d] [iay,iby]=[%d,%d] factorMatrix=%d x=[%d,%d]\n',Ngx,Ngy,iax,ibx,iay,iby,factorMatrix,size(gf{cur}.x,1),size(gf{cur}.x,2));
  pf(I,J) =  par.pfe(gf{cur}.x(I,J,1),gf{cur}.x(I,J,2),t);

  
  ux = zeros(Ngx,Ngy);
  uy = zeros(Ngx,Ngy);
  vx = zeros(Ngx,Ngy);
  vy = zeros(Ngx,Ngy);
  if( par.isCartesian )
    ux(I,J) = Dzx(u,I,J); 
    uy(I,J) = Dzy(u,I,J); 
    vx(I,J) = Dzx(v,I,J); 
    vy(I,J) = Dzy(v,I,J);  
  else
    % these next could be made faster by re-using u.r and u.s etc.
    ux(I,J) = Dx2(u,I,J); 
    uy(I,J) = Dy2(u,I,J); 
    vx(I,J) = Dx2(v,I,J); 
    vy(I,J) = Dy2(v,I,J);      
  end 

  % I = i1x:i2x; J=i1y:i2y; 
  [I1,I2] = getIndexInterior( par.gid,par.pc,par );
  for( i2=I2)
  for( i1=I1 )
    ie = eqn(i1,i2); % eqn number for pt (i1,i2) 

    % ux = Dzx(u,i1,i2); uy=Dzy(u,i1,i2); 
    % vx = Dzx(v,i1,i2); vy=Dzy(v,i1,i2); 
    divDamping = (cdv/dt)*(ux(i1,i2)+vy(i1,i2)); % divergence damping 

    rhs(ie)= -( ux(i1,i2)^2 + 2.*uy(i1,i2)*vx(i1,i2) + vy(i1,i2)^2) + divDamping + pf(i1,i2);
    % rhs(ie)= -( ux^2 + 2.*uy*vx + vy^2) + divDamping + pfe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);

    if( 1==0 ) % For TESTING 
      rhs(ie) = par.pexx(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) +  par.peyy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
    end

  end
  end

  % --- Boundary condition RIGHT HAND SIDES  ---

  %    NOTE: rhs for extrapolation and periodic equations are zero 
  if( ~strcmp(par.ms,'none') ) manufacturedSolution=1; else manufacturedSolution=0; end    
  if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'none') ) addBoundaryForcing=1; else addBoundaryForcing=0; end

  for side=1:2
    for axis=1:2

      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
      mbc = side+2*(axis-1);  % pointer into gp array
      [I1b,I2b]=getBoundaryIndex(side,axis,par);

      if( par.bc(side,axis)==par.traction && par.gamma~=0 )
        kappa = getCurvature( t,gf,cur, par );
      end 
      
      % Loop over boundary points on face=(side,axis)
      for i2=I2b
        for i1=I1b
          if( par.bc(side,axis)==par.dirichlet )
            ie=eqn(i1,i2); % boundary pt
            % rhs(ie)=par.gp{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t); 
            rhs(ie)=par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);

          elseif( par.bc(side,axis)==par.pressureInflow  )

            ie=eqn(i1,i2); % boundary pt
            
            if( addBoundaryForcing==1 )
              rhs(ie)=par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
            else
              rhs(ie)=par.pressureInflowValue;
            end

          elseif( par.bc(side,axis)==par.noSlipWall || ...
                  par.bc(side,axis)==par.slipWall   || ...
                  par.bc(side,axis)==par.inflow )

            % Neumann BC :
            ie=eqn(i1-is1,i2-is2); % ghost point
      
            if( par.isCartesian )

              % -- cartesian ---
              n1 = -is1; % outward normal [n1,n2]
              n2 = -is2; 
              if( axis==1 )
                % p.n = (+/-) nu*( u.xx ) = (+/-) nu*( -v.xy )
                rhs(ie)= -is1*nu*( -DzxDzy(v,i1,i2) ) + rho*n1*par.gravityVector(1);

                if( addBoundaryForcing==1 )
                  % vxy = DzxDzy(v,i1,i2); vxyTrue=vexy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t); 
                  % fprintf('--Pressure: par.noSlipWall : add boundary forcing, vxy=%8.2e vxye=%8.2e diff=%8.2e\n',vxy,vxyTrue,vxy-vxyTrue);
                  
                  rhs(ie)= rhs(ie) -is1*( par.pex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) + nu*par.vexy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ); 
                  %% rhs(ie)= -is1*( pex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ); % TEST 
                end
              else
                % p.n = (+/-) nu*( v.yy ) = (+/-) nu*( -u.xy )
                rhs(ie)= -is2*nu*( -DzxDzy(u,i1,i2) ) + rho*n2*par.gravityVector(2);
                if( addBoundaryForcing==1 )
                  rhs(ie)= rhs(ie) -is2*( par.pey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) + nu*par.uexy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ); 
                end
              end

            else 
              % -- curvilinear --
              %  p.n = rho * nu* [ n1 * ( nu*Delta(u) ) + n2 * ( nu Delta v ) ] + rho*( n1*gv(1) + n2*gv(2) )
              % CURL-CURL BC
              %   u.xx + u.yy -> -v.xy + u.yy
              %   v.xx + v.yy ->  v.xx - u.xy 
              % 

              % ---- get outward normal (n1,n2) ----
              is = 1-2*(side-1);
              n1 = -is*gf{cur}.rx(i1,i2,axis,1); 
              n2 = -is*gf{cur}.rx(i1,i2,axis,2); 
              rxNorm = sqrt( n1.^2 + n2.^2 ); 
              n1 = n1./rxNorm;
              n2 = n2./rxNorm;                  
              % --- done get normal ---

              if( 1==0 )
                % ****** TESTING set p.n = exact p.n ***
                rhs(ie)= n1*par.pex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) + n2*par.pey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
              else

                [UX,UXX] = getDerivatives( i1,i2,u,v, gf,cur, par );

                vxy = UXX(2,1,2);
                uyy = UXX(1,2,2);
                vxx = UXX(2,1,1);
                uxy = UXX(1,1,2);

                rhs(ie)= rho*nu*( n1*(-vxy+uyy) + n2*(vxx-uxy) ) + rho*(n1*par.gravityVector(1) + n2*par.gravityVector(2)); % curl-curl BC 

                if( addBoundaryForcing==1 )
                  rhs(ie)= rhs(ie) ...
                      + n1*par.pex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) + n2*par.pey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ...
                      - rho*nu*( n1*( -par.vexy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) + par.ueyy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ) + ...
                                 n2*(  par.vexx(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) - par.uexy(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) ) );
                end
              end
              
            end

          elseif( par.bc(side,axis)==par.outflow )

            % Outflow: a0*p + a1*p.n = a0*pOutflow
            ie=eqn(i1-is1,i2-is2); % ghost point
            a0 = par.outflowPressureCoeffp;
            a1 = par.outflowPressureCoeffpn;
            if( manufacturedSolution )
              rhs(ie) =            a0* par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t)  ...
                        -isv(axis)*a1*par.pex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
            else
              rhs(ie)= a0*par.pOutflow;  
            end  


          elseif( par.bc(side,axis)==par.traction )

            % ----- TRACTION BC ----
            %   nv^T sigmav nv = 0 

            ie=eqn(i1,i2); % boundary pt

            if( par.isCartesian )
              if( axis==1 )
                % p = 2*mu*u.x 
                rhs(ie) = 2*mu*Dzx(u,i1,i2);
                if( addBoundaryForcing==1 )
                   rhs(ie)= rhs(ie) + ...
                     par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) - 2*mu*par.uex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                end 
              else
                % p = 2*mu*v.y 
                rhs(ie) = 2*mu*Dzy(v,i1,i2);
                if( addBoundaryForcing==1 )
                   rhs(ie)= rhs(ie) + ...
                     par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t) - 2*mu*par.vey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                end                 
              end

            else
              % ------------------------------
              % --- traction : curvilinear ---
              % ------------------------------

              % nv^T sigmav nv = -p + 2*mu*( ux*n1^2 + vy*n2^2 + (uy+vx)*n1*n2 )   + gamma*kappa

              % ---- get outward normal (n1,n2) ----
              [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );
              % is = 1-2*(side-1);
              % n1 = -is*gf{cur}.rx(i1,i2,axis,1); 
              % n2 = -is*gf{cur}.rx(i1,i2,axis,2); 
              % rxNorm = sqrt( n1.^2 + n2.^2 ); 
              % n1 = n1./rxNorm;
              % n2 = n2./rxNorm;                  
              % --- done get normal ---

              rx = gf{cur}.rx(i1,i2,1,1); sx = gf{cur}.rx(i1,i2,2,1);
              ry = gf{cur}.rx(i1,i2,1,2); sy = gf{cur}.rx(i1,i2,2,2);  
              ux =  rx*Dr2(u,i1,i2) + sx*Ds2(u,i1,i2);
              vx =  rx*Dr2(v,i1,i2) + sx*Ds2(v,i1,i2);
              uy =  ry*Dr2(u,i1,i2) + sy*Ds2(u,i1,i2);
              vy =  ry*Dr2(v,i1,i2) + sy*Ds2(v,i1,i2);

              rhs(ie) = 2*mu*( ux*n1^2 + vy*n2^2 + (uy+vx)*n1*n2 );  

              if( par.gamma~=0 )
                % Note the sign of the curvatutre term here since the equation is
                %    p = nv^T \tauv nv - gamma kappa
                rhs(ie) = rhs(ie) - par.gamma*kappa(i1); 
              end   
                                  
              if( addBoundaryForcing==1 )
                uex = par.uex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                uey = par.uey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                vex = par.vex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);
                vey = par.vey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t);

                rhs(ie)= rhs(ie) + ...
                     par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t)  ...
                     - 2*mu*( uex*n1^2 + vey*n2^2 + (uey+vex)*n1*n2 );            
              end
               % fprintf('pressure matrx: fill-BC: finish me for traction BC and curvlinear\n');
               % error('error');
            end 


          elseif( par.bc(side,axis)==par.periodic )
            % rhs should be zero already
            ie=eqn(i1-is1,i2-is2); % ghost point
            rhs(ie)=0.; 
          else
            fprintf('pressure matrx: fill-BC: finish me...\n'); pause; 
          end

        end
      end
    end
  end


   if( isSingular && addBoundaryForcing==1 )
      % set the "mean" value in the extra equation to the "mean" value of the exact solution
      extraVal=0; 
      % I1=iax:ibx; I2=iay:iby; 
      [I1,I2]=getIndex( par.gid );
      for i2=I2
        for i1=I1
          extraVal = extraVal + par.pe(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),t); 
        end
      end
      rhs(Ngs)=extraVal;
      if( mod(floor(idebug/2),2)==1 ) fprintf('RHS to extra equation=%12.4e\n',rhs(Ngs)); end 
   end

   % rhs = U\( L\rhs );  % solve LU (rhs) = rhs 
   rhs = par.dA\rhs;  % solve LU (rhs) = rhs 

   p   = zeros(Ngx,Ngy);   % holds p

   % copy vector solution to the grid function solution
   % I1g = iax-numGhost:ibx+numGhost; I2g=iay-numGhost:iby+numGhost; % include ghost points
   [I1g,I2g] = getIndex( par.gid,numGhost ); % include ghost points
   for( i1=I1g )
   for( i2=I2g )
     ie = eqn(i1,i2); % eqn number for pt (i1,i2) 
     p(i1,i2)=rhs(ie);
   end
   end
   if( isSingular==1 && mod(floor(idebug/2),2)==1 )
     fprintf('Solution to extra equation=%12.4e\n',rhs(Ngs));
   end
    
   % check errors
   if( mod(floor(idebug/2),2)==1 )
     pTrue = par.pe(gf{cur}.x(I1g,I2g,1),gf{cur}.x(I1g,I2g,2),t);

     [I1,I2] = getIndex(par.gid);
     err = p(I1,I2) - pTrue(I1,I2);
     % maxErr = norm(err,inf); % TROUBLE %
     maxErr=max(max(abs(err)));
     fprintf('>>> Poisson: t=%10.3e, dx=%g dy=%g max error=%9.2e\n',t,dx,dy,maxErr); 
  
     if( numGhost>0 )
       % check errors on ghost too
       err = p(I1g,I2g) - pTrue(I1g,I2g);
       maxErr=max(max(abs(err)));
       fprintf('>>> Poisson: t=%10.3e, dx=%g dy=%g max error=%9.2e (with ghost)\n',t,dx,dy,maxErr); 
     end
     
     if( 1==1 || mod(floor(par.plotOption/2),2)== 1 )
       figure(4)
       surf(gf{cur}.x(:,:,1),gf{cur}.x(:,:,2),p); hold on;
       contour3( gf{cur}.x(:,:,1),gf{cur}.x(:,:,2),p,'k-' ); 
       colormap(par.rainbowMap); colorbar; shading interp; 
       xlabel('x'); ylabel('y'); view(0,90);  % top view 
       title(sprintf('p : t=%8.2e Nx=%d Ny=%d',t,Ngx,Ngy)); xlabel('x'); ylabel('y'); set(gca,'FontSize',14);
       hold off;
       setAspectRatio();       
    
       % figure(2)
       % surf(x,y,pTrue);
       % title('p-exact'); xlabel('x'); ylabel('y');
       % drawnow; commandwindow; 
  
       figure(5)
       surf(gf{cur}.x(:,:,1),gf{cur}.x(:,:,2),err); hold on;
       contour3( gf{cur}.x(:,:,1),gf{cur}.x(:,:,2),err,'k-' ); 
       colormap(par.rainbowMap); colorbar; shading interp; 
       xlabel('x'); ylabel('y'); view(0,90);  % top view 

       title(sprintf('p-err : t=%8.2e Nx=%d Ny=%d',t,Ngx,Ngy)); xlabel('x'); ylabel('y'); set(gca,'FontSize',14);
       hold off;
       setAspectRatio();   

       drawnow; commandwindow; 

       pause; 
     end 
   end

  par.cpuPressure = par.cpuPressure + cputime - cpu0; % time for pressure SOLVES (not factor)


  % fill matrix entries
  function setValue( ii,jj,val )
   nzz=nzz+1;
   ia(nzz)=ii; ja(nzz)=jj; aa(nzz)=val; 
  end  


 % return
end
