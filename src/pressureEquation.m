%
%  Solve the Poisson equation for p 
%
% Input:
%  factorMatrix =1 : factor the matrix
%
function [p,par] = pressureEquation( t,u,v,dt,factorMatrix,par )

 % fprintf('Entering pressureEquation u=[%d,%d] v=[%d,%d]\n',...
 %         size(u,1),size(u,2), size(v,1),size(v,2) );

  nu  = par.nu;
  cdv = par.cdv;
  ms  = par.ms;

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

  % Define extrapolations: (is1=+1/-1 and is2=+1/-1 defines the direction ("shift") of extrapolation)
  extrap3 = @(u,I1,I2,is1,is2) (3.*u(I1+is1,I2+is2) - 3.*u(I1+2*is1,I2+2*is2) + u(I1+3*is1,I2+3*is2));   % 3rd-order extrapolation

  % convert (ix,iy) to equation number in the matrix:  
  eqn = @(ix,iy)  1 + ix-1 + Ngx*( iy-1 );

    
  % persistent NgSave L U ; % holds LU factors of the pressure matrix

  Ng = Ngx*Ngy; % total number of grid points

  % check if the Pressure matrix is singular
  isSingular = 1;
  for side=1:2
    for axis=1:2
      if( par.bc(side,axis)==par.dirichlet || ...
          par.bc(side,axis)==par.pressureInflow || ...
          par.bc(side,axis)==par.outflow )
        isSingular=false;
      end
    end
  end
  Ngs =Ng;
  if( isSingular==1 ) Ngs=Ng+1; end % add an extra equation for the singular case
  


  if( factorMatrix==1 )  % isempty(L) || Ng ~= NgSave 

    % ----- Form the pressure matrix ------
    

    % Allocate the sparse matrix, vector-solution and RHS

    if( mod(floor(idebug/2),2)==1 )
      fprintf('*** FORM THE PRESSURE MATRIX : Ng=%d t=%12.4e **\n',Ng,t);
    end
    
    % NgSave=Ng; % save 

    rightNullValue=0;
    if( isSingular==1 )
      if( par.idebug>0 ) fprintf('pressureEquation: **MATRIX IS SINGULAR** add an extra equation.\n'); end
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

    for( ix=I1 )
    for( iy=I2 )
      ie = eqn(ix,iy) ; % eqn number for pt (ix,iy) 
      setValue( ie,eqn(ix  ,iy-1),(               1/dy^2   ) );
      setValue( ie,eqn(ix-1,iy  ),(      1/dx^2            ) );
      setValue( ie,eqn(ix  ,iy  ),( -2*( 1/dx^2 + 1/dy^2 ) ) );
      setValue( ie,eqn(ix+1,iy  ),(      1/dx^2            ) );
      setValue( ie,eqn(ix  ,iy+1),(               1/dy^2   ) );
      if( isSingular==1 ) 
        setValue(ie,Ngs,1.);
      end 
    end
    end
    
    % add extra equation in the last row -- all ones except that last column
    if( isSingular==1 )
      ie=Ngs; % final equation
      for( ix=I1 )
      for( iy=I2 )
        setValue(ie,eqn(ix,iy),1.);
      end
      end
    end
     
    % Boundary conditions for the matrix 
    dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

    for side=1:2
      for axis=1:2
        % isv(1:2), is1, is2 : index shifts
        isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
  
        [I1b,I2b]=getBoundaryIndex(side,axis,par);
        for iy=I2b
          for ix=I1b
            if( par.bc(side,axis)==par.dirichlet || ...
                par.bc(side,axis)==par.pressureInflow )
              ie=eqn(ix,iy); 
              setValue(ie,ie,1.);
              % Extrapolate ghost
              ie=eqn(ix-is1,iy-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(ix      ,iy      ),-3.); 
              setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
              setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.); 

            elseif( par.bc(side,axis)==par.noSlipWall || ...
                    par.bc(side,axis)==par.slipWall   || ...
                    par.bc(side,axis)==par.inflow )
              % Neumann BC :
                %    (+-)*Dz( p ) = RHS  : use outward normal 
              % nSign = 2*(side-1)-1; % sign of normal, -1 on left and +1 on right
              % fprintf('pressureEqn: fill-in a Neumann BC (side,axis)=(%d,%d) (ix,iy)=(%d,%d) (is1,is2)=(%d,%d) nSign=%g\n',side,axis,ix,iy,is1,is2,nSign); 

              ie=eqn(ix-is1,iy-is2); % ghost point
              setValue(ie,eqn(ix-is1,iy-is2), 1/(2.*dxv(axis))); % ghost point 
              setValue(ie,eqn(ix+is1,iy+is2),-1/(2.*dxv(axis))); % first interior point

            elseif( par.bc(side,axis)==par.outflow )
              % a0*p + a1*p.n = RHS 

              % This equation is put into the matrix at the ghost point: 
              ie=eqn(ix-is1,iy-is2);     % ghost point 
              a0 = par.outflowPressureCoeffp;
              setValue(ie,eqn(ix,iy),a0); % 

              ie=eqn(ix-is1,iy-is2); % ghost point
              a1 = par.outflowPressureCoeffpn;
              setValue(ie,eqn(ix-is1,iy-is2), a1/(2.*dxv(axis))); % ghost point 
              setValue(ie,eqn(ix+is1,iy+is2),-a1/(2.*dxv(axis))); % first interior point                

              
            elseif( par.bc(side,axis)==par.periodic )
              % Periodic: p(iax-1,.) = p(ibx-1,.) ...
              ie=eqn(ix-is1,iy-is2); % ghost point
              setValue(ie,ie                                            , 1.); 
              setValue(ie,eqn(ix+(ibx-iax)*is1-is1,iy+(iby-iay)*is2-is2),-1.); 
            else
              fprintf('pressure matrx: fill-BC: finish me...\n'); pause; 
            end
          end
        end
      end
    end


    % extrapolate corners along the diagonal 
    % par.gid(side,axis) : grid index range
    % par.gid(1,1)=iax; par.gid(2,1)=ibx; par.gid(1,2)=iay; par.gid(2,2)=iby; 
    for( side1=0:1 )
      for( side2=0:1 )
        is1 = 1-2*side1;
        is2 = 1-2*side2;
        ix=par.gid(side1+1,1)-is1; iy=par.gid(side2+1,2)-is2; % corner ghost point
        ie=eqn(ix,iy); 
        setValue(ie,ie                    , 1.);
        setValue(ie,eqn(ix+  is1,iy+  is2),-3.); 
        setValue(ie,eqn(ix+2*is1,iy+2*is2), 3.);  
        setValue(ie,eqn(ix+3*is1,iy+3*is2),-1.); 
      end
    end 

    if( par.idebug>0 )
      fprintf('Optimized fill of pressure matrix: nzzEst=%d, nzz=%d\n',nzzEst,nzz);
    end
    A = sparse(ia(1:nzz), ja(1:nzz), aa(1:nzz), Ngs, Ngs); % Creates the sparse matrix 

    
    % else 
    %   % ---- OLD WAY TO FILL MATRIX ---
    %   A = sparse(Ngs,Ngs);   % implicit matrix 
    %   % L = sparse(Ngs,Ngs);   % implicit matrix 
    %   % U = sparse(Ngs,Ngs);   % implicit matrix 
      


    %   % Fill in the interior equations (and also the boundary pts for Neumann BC's )

    %   [I1,I2] = getIndexInterior( par.gid,par.pc,par );

    %   for( ix=I1 )
    %   for( iy=I2 )
    %     ie = eqn(ix,iy) ; % eqn number for pt (ix,iy) 
    %     A(ie,eqn(ix  ,iy-1)) = (               1/dy^2   );
    %     A(ie,eqn(ix-1,iy  )) = (      1/dx^2            );
    %     A(ie,eqn(ix  ,iy  )) = ( -2*( 1/dx^2 + 1/dy^2 ) );
    %     A(ie,eqn(ix+1,iy  )) = (      1/dx^2            );
    %     A(ie,eqn(ix  ,iy+1)) = (               1/dy^2   );
    %     if( isSingular==1 ) A(ie,Ngs)=1.; end 
    %   end
    %   end
      
    %   % add extra equation in the last row -- all ones except that last column
    %   if( isSingular==1 )
    %     ie=Ngs; % final equation
    %     for( ix=I1 )
    %     for( iy=I2 )
    %       A(ie,eqn(ix,iy))=1;       
    %     end
    %     end
    %   end
       
    %   % Boundary conditions for the matrix 
    %   dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

    %   for side=1:2
    %     for axis=1:2
    %       % isv(1:2), is1, is2 : index shifts
    %       isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
    
    %       [I1b,I2b]=getBoundaryIndex(side,axis,par);
    %       for iy=I2b
    %         for ix=I1b
    %           if( par.bc(side,axis)==par.dirichlet )
    %             ie=eqn(ix,iy); A(ie,ie)=1.;
    %             % Extrapolate ghost
    %             ie=eqn(ix-is1,iy-is2); A(ie,ie)=1.; A(ie,eqn(ix,iy))=-3.; A(ie,eqn(ix+is1,iy+is2))=3.;  A(ie,eqn(ix+2*is1,iy+2*is2))=-1.; 

    %           elseif( par.bc(side,axis)==par.noSlipWall || ...
    %                   par.bc(side,axis)==par.slipWall   || ...
    %                   par.bc(side,axis)==par.inflow )
    %             % Neumann BC :
    %               %    (+-)*Dz( p ) = RHS  : use outward normal 
    %             nSign = 2*(side-1)-1; % sign of normal, -1 on left and +1 on right
    %             % fprintf('pressureEqn: fill-in a Neumann BC (side,axis)=(%d,%d) (ix,iy)=(%d,%d) (is1,is2)=(%d,%d) nSign=%g\n',side,axis,ix,iy,is1,is2,nSign); 

    %             ie=eqn(ix-is1,iy-is2); % ghost point
    %             A(ie,eqn(ix-is1,iy-is2))=  1/(2.*dxv(axis)); % ghost point 
    %             A(ie,eqn(ix+is1,iy+is2))= -1/(2.*dxv(axis)); % first interior point

    %           elseif( par.bc(side,axis)==par.periodic )
    %             % Periodic: p(iax-1,.) = p(ibx-1,.) ...
    %             ie=eqn(ix-is1,iy-is2); % ghost point
    %             A(ie,ie)=1.; A(ie,eqn(ix+(ibx-iax)*is1-is1,iy+(iby-iay)*is2-is2))=-1; 
    %           else
    %             fprintf('pressure matrx: fill-BC: finish me...\n'); pause; 
    %           end
    %         end
    %       end
    %     end
    %   end


    %   % extrapolate corners along the diagonal 
    %   % par.gid(side,axis) : grid index range
    %   % par.gid(1,1)=iax; par.gid(2,1)=ibx; par.gid(1,2)=iay; par.gid(2,2)=iby; 
    %   for( side1=0:1 )
    %     for( side2=0:1 )
    %       is1 = 1-2*side1;
    %       is2 = 1-2*side2;
    %       ix=par.gid(side1+1,1)-is1; iy=par.gid(side2+1,2)-is2; % corner ghost point
    %       ie=eqn(ix,iy); 
    %       A(ie,ie)=1.; A(ie,eqn(ix+is1,iy+is2))=-3.; A(ie,eqn(ix+2*is1,iy+2*is2))=3.;  A(ie,eqn(ix+3*is1,iy+3*is2))=-1.; 
    %     end
    %   end 
  
    % end % end old way 


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

  % -- Assign the RHS  --
  rhs = zeros(Ngs,1);     % right-hand-side

  % I = iax:ibx; J=iay:iby;  % interior and boundary points
  [I,J] = getIndex( par.gid );

  pf = zeros(Ngx,Ngy);
  % fprintf('Ngx=%d, Ngy=%d, [iax,ibx]=[%d,%d] [iay,iby]=[%d,%d] factorMatrix=%d x=[%d,%d]\n',Ngx,Ngy,iax,ibx,iay,iby,factorMatrix,size(par.x,1),size(par.x,2));
  pf(I,J) =  par.pfe(par.x(I,J,1),par.x(I,J,2),t);

  useOpt=0;
  if(  useOpt )
    % THIS WAS ABOUT THE SAME TIME 
    [I1,I2] = getIndexInterior( par.gid,par.pc,par );
    if( 1==0 )
      r=zeros(Ngx,Ngy);
      for( i2=I2)
      for( i1=I1 )
        ux(i1,i2) = Dzx(u,i1,i2); 
        uy(i1,i2) = Dzy(u,i1,i2); 
        vx(i1,i2) = Dzx(v,i1,i2); 
        vy(i1,i2) = Dzy(v,i1,i2);               
        r(i1,i2)= -( ux(i1,i2).^2 + 2.*uy(i1,i2).*vx(i1,i2) + vy(i1,i2).^2) + (cdv/dt)*(ux(i1,i2)+vy(i1,i2)) + pf(i1,i2); 
      end
      end
      r = reshape( r, [Ng,1]) ;
      rhs(1:Ng) = r(1:Ng);      
    else
      r=zeros(Ngx,Ngy);
      ux(I,J) = Dzx(u,I,J); 
      uy(I,J) = Dzy(u,I,J); 
      vx(I,J) = Dzx(v,I,J); 
      vy(I,J) = Dzy(v,I,J);       
      r(I1,I2)= -( ux(I1,I2).^2 + 2.*uy(I1,I2).*vx(I1,I2) + vy(I1,I2).^2) + (cdv/dt)*(ux(I1,I2)+vy(I1,I2)) + pf(I1,I2); 
      r = reshape( r, [Ng,1]) ;
      rhs(1:Ng) = r(1:Ng);
    end

  else
    ux = zeros(Ngx,Ngy);
    uy = zeros(Ngx,Ngy);
    vx = zeros(Ngx,Ngy);
    vy = zeros(Ngx,Ngy);
    ux(I,J) = Dzx(u,I,J); 
    uy(I,J) = Dzy(u,I,J); 
    vx(I,J) = Dzx(v,I,J); 
    vy(I,J) = Dzy(v,I,J);   

    % I = i1x:i2x; J=i1y:i2y; 
    [I1,I2] = getIndexInterior( par.gid,par.pc,par );
    for( iy=I2)
    for( ix=I1 )
      ie = eqn(ix,iy); % eqn number for pt (ix,iy) 

      % ux = Dzx(u,ix,iy); uy=Dzy(u,ix,iy); 
      % vx = Dzx(v,ix,iy); vy=Dzy(v,ix,iy); 
      divDamping = (cdv/dt)*(ux(ix,iy)+vy(ix,iy)); % divergence damping 

      rhs(ie)= -( ux(ix,iy)^2 + 2.*uy(ix,iy)*vx(ix,iy) + vy(ix,iy)^2) + divDamping + pf(ix,iy);
      % rhs(ie)= -( ux^2 + 2.*uy*vx + vy^2) + divDamping + pfe(par.x(ix,iy,1),par.x(ix,iy,2),t);

    end
    end
  end

  % --- Boundary conditions ---
  %    NOTE: rhs for extrapolation and periodic equations are zero 
  if( ~strcmp(par.ms,'none') ) manufacturedSolution=1; else manufacturedSolution=0; end    
  if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'none') ) addBoundaryForcing=1; else addBoundaryForcing=0; end

  for side=1:2
    for axis=1:2

      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
      mbc = side+2*(axis-1);  % pointer into gp array
      [I1b,I2b]=getBoundaryIndex(side,axis,par);

      
      % Loop over boundary points on face=(side,axis)
      for iy=I2b
        for ix=I1b
          if( par.bc(side,axis)==par.dirichlet )
            ie=eqn(ix,iy); % boundary pt
            % rhs(ie)=par.gp{mbc}(par.x(ix,iy,1),par.x(ix,iy,2),t); 
            rhs(ie)=par.pe(par.x(ix,iy,1),par.x(ix,iy,2),t);

          elseif( par.bc(side,axis)==par.pressureInflow  )

            ie=eqn(ix,iy); % boundary pt
            
            if( addBoundaryForcing==1 )
              rhs(ie)=par.pe(par.x(ix,iy,1),par.x(ix,iy,2),t);
            else
              rhs(ie)=par.pressureInflowValue;
            end

          elseif( par.bc(side,axis)==par.noSlipWall || ...
                  par.bc(side,axis)==par.slipWall   || ...
                  par.bc(side,axis)==par.inflow )

            % Neumann BC :
            ie=eqn(ix-is1,iy-is2); % ghost point
            if( 1==0 ) % TEMP ***
              rhs(ie)=par.gp{mbc}(par.x(ix,iy,1),par.x(ix,iy,2),t);  % test: give exact p.n 
            else
      
              if( axis==1 )
                % p.n = (+/-) nu*( u.xx ) = (+/-) nu*( -v.xy )
                rhs(ie)= -is1*nu*( -DzxDzy(v,ix,iy) );
                if( addBoundaryForcing==1 )
                  % vxy = DzxDzy(v,ix,iy); vxyTrue=vexy(par.x(ix,iy,1),par.x(ix,iy,2),t); 
                  % fprintf('--Pressure: par.noSlipWall : add boundary forcing, vxy=%8.2e vxye=%8.2e diff=%8.2e\n',vxy,vxyTrue,vxy-vxyTrue);
                  
                  rhs(ie)= rhs(ie) -is1*( par.pex(par.x(ix,iy,1),par.x(ix,iy,2),t) + nu*par.vexy(par.x(ix,iy,1),par.x(ix,iy,2),t) ); 
                  %% rhs(ie)= -is1*( pex(par.x(ix,iy,1),par.x(ix,iy,2),t) ); % TEST 
                end
              else
                % p.n = (+/-) nu*( v.yy ) = (+/-) nu*( -u.xy )
                rhs(ie)= -is2*nu*( -DzxDzy(u,ix,iy) );
                if( addBoundaryForcing==1 )
                  rhs(ie)= rhs(ie) -is2*( par.pey(par.x(ix,iy,1),par.x(ix,iy,2),t) + nu*par.uexy(par.x(ix,iy,1),par.x(ix,iy,2),t) ); 
                end
              end
            end
          elseif( par.bc(side,axis)==par.outflow )

            % Outflow: a0*p + a1*p.n = a0*pOutflow
            ie=eqn(ix-is1,iy-is2); % ghost point
            a0 = par.outflowPressureCoeffp;
            a1 = par.outflowPressureCoeffpn;
            if( manufacturedSolution )
              rhs(ie) =            a0* par.pe(par.x(ix,iy,1),par.x(ix,iy,2),t)  ...
                        -isv(axis)*a1*par.pex(par.x(ix,iy,1),par.x(ix,iy,2),t);
            else
              rhs(ie)= a0*par.pOutflow;  
            end  

          elseif( par.bc(side,axis)==par.periodic )
            % rhs should be zero already
            ie=eqn(ix-is1,iy-is2); % ghost point
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
      for iy=I2
        for ix=I1
          extraVal = extraVal + par.pe(par.x(ix,iy,1),par.x(ix,iy,2),t); 
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
   for( ix=I1g )
   for( iy=I2g )
     ie = eqn(ix,iy); % eqn number for pt (ix,iy) 
     p(ix,iy)=rhs(ie);
   end
   end
   if( isSingular==1 && mod(floor(idebug/2),2)==1 )
     fprintf('Solution to extra equation=%12.4e\n',rhs(Ngs));
   end
    
   % check errors
   if( mod(floor(idebug/2),2)==1 )
     pTrue = par.pe(par.x(I1g,I2g,1),par.x(I1g,I2g,2),t);
     err = p(I,J) - pTrue(I,J);
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
       surf(par.x(:,:,1),par.x(:,:,2),p);
       title(sprintf('p : t=%8.2e Nx=%d Ny=%d',t,Ngx,Ngy)); xlabel('x'); ylabel('y'); set(gca,'FontSize',14);
    
       % figure(2)
       % surf(x,y,pTrue);
       % title('p-exact'); xlabel('x'); ylabel('y');
       % drawnow; commandwindow; 
  
       figure(5)
       surf(par.x(:,:,1),par.x(:,:,2),err);
       title(sprintf('p-err : t=%8.2e Nx=%d Ny=%d',t,Ngx,Ngy)); xlabel('x'); ylabel('y'); set(gca,'FontSize',14);
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
