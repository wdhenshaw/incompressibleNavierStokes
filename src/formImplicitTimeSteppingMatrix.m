function par = formImplicitTimeSteppingMatrix( dt, par )


  cpu0=cputime;

  nu = par.nu;

  Ngx = par.Ngx;
  Ngy = par.Ngy;

  % iax = par.iax;  ibx = par.ibx; 
  % iay = par.iay;  iby = par.iby;  

  iax = par.gid(1,1); ibx=par.gid(2,1);
  iay = par.gid(1,2); iby=par.gid(2,2);   

  numGhost = par.numGhost;  

  % i1x = par.i1x;  i2x = par.i2x; 
  % i1y = par.i1y;  i2y = par.i2y; 

  dx = par.dx;
  dy = par.dy;   

  % convert (ix,iy) to equation number in the matrix:  
  eqn = @(ix,iy)  1 + ix-1 + Ngx*( iy-1 );
  
  Ng = Ngx*Ngy; % total number of grid points


  if( mod(floor(par.idebug/2),2)==1 )
    fprintf('...FORM the IMPLICIT TIME-STEPPING MATRIX: Ng=%d\n',Ng); 
  end 

  if( par.multipleImplicitSolversNeeded )
    numImplicitSolvers=2;
    if( par.idebug>1 )
      fprintf('*** formImplicitTimeSteppingMatrix: INFO -- multiple implicit solvers are needed sincs BCs do not match.\n');
    end
  else
    numImplicitSolvers=1;
  end
  

  for( iuv=1:numImplicitSolvers )

    % -------------------------------
    % ---- OPTIMIZED FILL MATRIX ----
    % -------------------------------

    nzzEst = Ng*5; % estimated number of non-zeros
    ia=zeros(nzzEst,1); ja=zeros(nzzEst,1); aa=zeros(nzzEst,1);

    nzz=0; % counts non-zeros

    if( iuv==1 ) component = par.uc; else component=par.vc; end
    [I1a,I2a] = getIndexInterior( par.gid,component,par );

    if( par.multipleImplicitSolversNeeded )
      if( par.idebug>1 )
        fprintf('formImplicitTimeSteppingMatrix: solver %d: I1a=[%d,%d] I2a=[%d,%d] (interior eqn)\n',iuv,I1a(1),I1a(end), I2a(1),I2a(end));
      end
    end 

    for( iy=I2a )
    for( ix=I1a )
      ie = eqn(ix,iy) ; % eqn number for pt (ix,iy) 
      setValue(ie,eqn(ix  ,iy-1),    - (.5*nu*dt)*(               1/dy^2   ));
      setValue(ie,eqn(ix-1,iy  ),    - (.5*nu*dt)*(      1/dx^2            ));
      setValue(ie,eqn(ix  ,iy  ), 1. - (.5*nu*dt)*( -2*( 1/dx^2 + 1/dy^2 ) ));
      setValue(ie,eqn(ix+1,iy  ),    - (.5*nu*dt)*(      1/dx^2            ));
      setValue(ie,eqn(ix  ,iy+1),    - (.5*nu*dt)*(               1/dy^2   ));
    end
    end
    
    % Boundary conditions for the matrix 
    if( 1==1 )
      % **NEW WAY**
      % Boundary conditions for the matrix 
      dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

      for side=1:2
      for axis=1:par.nd
          % isv(1:2), is1, is2 : index shifts
          isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
    
          [I1b,I2b]=getBoundaryIndex(side,axis,par);  

          if( par.bc(side,axis)==par.dirichlet   || ...
              par.bc(side,axis)==par.inflow      || ...
              par.bc(side,axis)==par.noSlipWall )
            % Give u = g
            for iy=I2b
            for ix=I1b
              ie=eqn(ix,iy); 
              setValue(ie,ie,1.);
              % Extrapolate ghost
              ie=eqn(ix-is1,iy-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(ix      ,iy      ),-3.); 
              setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
              setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.); 
            end
          end
        elseif( par.bc(side,axis)==par.slipWall )
          %   nv.uv   = given 
          %  (tv.uv).n = given          
          for iy=I2b
          for ix=I1b   
            if( (axis==1 && iuv==1) || (axis==2 && iuv==2) )
              % Dirichlet condition on this component
              ie=eqn(ix,iy);         % boundary pt
              setValue(ie,ie,1.); 

              % Extrapolate ghost
              ie=eqn(ix-is1,iy-is2); % ghost pt
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(ix      ,iy      ),-3.); 
              setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
              setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.); 

            else
              % Neumann condition on this component
              ie=eqn(ix-is1,iy-is2); % ghost point
              setValue(ie,eqn(ix-is1,iy-is2), 1/(2.*dxv(axis))); % ghost point 
              setValue(ie,eqn(ix+is1,iy+is2),-1/(2.*dxv(axis))); % first interior point              
            end            

          end
          end       

        elseif( par.bc(side,axis)==par.pressureInflow )
          %   tv.uv = given 
          %   nv.uv = extrapolated           
          for iy=I2b
          for ix=I1b   
            if( (axis==1 && iuv==2) || (axis==2 && iuv==1) )
              % Dirichlet condition on this component
              ie=eqn(ix,iy); 
              setValue(ie,ie,1.);  
              % Extrapolate ghost
              ie=eqn(ix-is1,iy-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(ix      ,iy      ),-3.); 
              setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
              setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.); 

            else
              % Extrapolate ghost
              ie=eqn(ix-is1,iy-is2); 
              setValue(ie,ie,                     1.); 
              setValue(ie,eqn(ix      ,iy      ),-3.); 
              setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
              setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.);               
            end                        
          end
          end            

        elseif( par.bc(side,axis)==par.outflow )
          % extrapolation
          % **CHECK ME** Maybe should use [1,-2,1] extrpolation
          for iy=I2b
          for ix=I1b   
            % Extrapolate ghost
            ie=eqn(ix-is1,iy-is2); 
            setValue(ie,ie,                     1.); 
            setValue(ie,eqn(ix      ,iy      ),-3.); 
            setValue(ie,eqn(ix+  is1,iy+  is2), 3.);  
            setValue(ie,eqn(ix+2*is1,iy+2*is2),-1.); 
          end
          end            

        elseif( par.bc(side,axis)==par.periodic )
          for iy=I2b
          for ix=I1b   
            % Periodic: u(iax-1,.) = u(ibx-1,.) ...
            ie=eqn(ix-is1,iy-is2); % ghost point
            setValue(ie,ie                                            , 1.); 
            setValue(ie,eqn(ix+(ibx-iax)*is1-is1,iy+(iby-iay)*is2-is2),-1.);             
          end
          end

        else 
          fprintf('formImplicitTimeSteppingMatrix: Error: unknown bc=%d\n',par.bc(side,axis));
          pause
        end         

      end % end for axis
      end % end for side

    else
      % *** OLD WAY ***
      if( par.bc(1,1)~=par.periodic ) % not par.periodic in x 
        for( iy=iay:iby )
          ix=iax; ie=eqn(ix,iy);  
          setValue(ie,ie,1.);  % BC at x=ax
          ix=ibx; ie=eqn(ix,iy);  
          setValue(ie,ie,1.);  % BC at x=bx 
          % extrapolate ghost 
          ix=iax-1; ie=eqn(ix,iy); 
          setValue(ie,ie          , 1.); 
          setValue(ie,eqn(ix+1,iy),-3.); 
          setValue(ie,eqn(ix+2,iy), 3.);  
          setValue(ie,eqn(ix+3,iy),-1.); 
          ix=ibx+1; ie=eqn(ix,iy); 
          setValue(ie,ie          , 1.); 
          setValue(ie,eqn(ix-1,iy),-3.); 
          setValue(ie,eqn(ix-2,iy), 3.);  
          setValue(ie,eqn(ix-3,iy),-1.); 
        end
      else
        % par.periodic in x 
        for( iy=iay:iby )
          ix=iax-1; ie=eqn(ix,iy); 
          setValue(ie,ie           , 1.); 
          setValue(ie,eqn(ibx-1,iy),-1.);   % u(iax-1,.) = u(ibx-1,.) 
          ix=ibx+1; ie=eqn(ix,iy); 
          setValue(ie,ie           , 1.); 
          setValue(ie,eqn(iax+1,iy),-1.);   % u(ibx+1,.) = u(iax+1,.)
        end    
      end

      if( par.bc(1,2)~=par.periodic )  % not par.periodic in y 
        for( ix=iax:ibx )
          iy=iay; ie=eqn(ix,iy);  
          setValue(ie,ie,1.);  % BC at y=ay
          iy=iby; ie=eqn(ix,iy);  
          setValue(ie,ie,1.);  % BC at y=by
          % extrapolate ghost
          iy=iay-1; ie=eqn(ix,iy); 
          setValue(ie,ie          , 1.); 
          setValue(ie,eqn(ix,iy+1),-3.); 
          setValue(ie,eqn(ix,iy+2), 3.);  
          setValue(ie,eqn(ix,iy+3),-1.);
          iy=iby+1; ie=eqn(ix,iy); 
          setValue(ie,ie          , 1.); 
          setValue(ie,eqn(ix,iy-1),-3.); 
          setValue(ie,eqn(ix,iy-2), 3.);  
          setValue(ie,eqn(ix,iy-3),-1.);
        end
      else
        % par.periodic in y 
        for( ix=iax:ibx )
          iy=iay-1; ie=eqn(ix,iy); 
          setValue(ie,ie           , 1.); 
          setValue(ie,eqn(ix,iby-1),-1.);  % u(.,iay-1)=u(.,iby-1) 
          iy=iby+1; ie=eqn(ix,iy); 
          setValue(ie,ie           , 1.); 
          setValue(ie,eqn(ix,iay+1),-1.);  % u(.,iby+1)=u(.,iay+1) 
        end
      end


    end % end old way

    % ----- CORNERS ---
    if( par.bc(1,1)~=par.periodic ||  par.bc(1,2)~=par.periodic )
      % extrapolate corners along the diagonal 
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
    else
      % par.periodic in both x and y 
      for( side1=0:1 )
        for( side2=0:1 )
          is1 = 1-2*side1;
          is2 = 1-2*side2;
          ix=par.gid(side1+1,1)-is1; iy=par.gid(side2+1,2)-is2; % corner ghost point
          ie=eqn(ix,iy);  
          setValue(ie,ie                                    , 1.); 
          setValue(ie,eqn(ix+(ibx-iax)*is1,iy+(iby-iay)*is2),-1.); 
        end
      end 
    end     

    if( par.idebug>0 )
      fprintf('Optimized fill of implicit time-stepping matrix: nzzEst=%d, nzz=%d\n',nzzEst,nzz);
    end

    AA = sparse(ia(1:nzz), ja(1:nzz), aa(1:nzz), Ng, Ng); % Creates the sparse matrix     




    if( mod(floor(par.idebug/2),2)==1 )
      fprintf('+++ cond(AA)=%10.2e\n',condest(AA));
      figure(3); 
      spy(AA);
      pause
    end 

    cpuFillImpMatrix = cputime - cpu0;

    cpu0 = cputime;
    par.dAimp{iuv} =decomposition(AA);   % factor the matrix: A=LU  

    cpuFactor = cputime-cpu0;
    par.cpuFactorImpMatrix = par.cpuFactorImpMatrix + cpuFactor;

    if( par.idebug>0 )
      fprintf('formImplicitMatrix: time to fill matrix = %8.2e(s), factor matrix = %8.2e(s)\n',cpuFillImpMatrix,cpuFactor);
    end

  end % end for iuv 

  % fill matrix entries
  function setValue( ii,jj,val )
   nzz=nzz+1;
   ia(nzz)=ii; ja(nzz)=jj; aa(nzz)=val; 
  end  

return
end