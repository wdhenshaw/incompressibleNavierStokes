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
  dr = par.dr(1);
  ds = par.dr(2);   

  % derivatives of entreis in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  % Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r
  % Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s  

  % convert (i1,i2) to equation number in the matrix:  
  eqn = @(i1,i2)  1 + i1-1 + Ngx*( i2-1 );
  
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

    if( par.isCartesian )

      % --- Cartesian grid ---

      for( i2=I2a )
      for( i1=I1a )
        ie = eqn(i1,i2) ; % eqn number for pt (i1,i2) 
        setValue(ie,eqn(i1  ,i2-1),    - (.5*nu*dt)*(               1/dy^2   ));
        setValue(ie,eqn(i1-1,i2  ),    - (.5*nu*dt)*(      1/dx^2            ));
        setValue(ie,eqn(i1  ,i2  ), 1. - (.5*nu*dt)*( -2*( 1/dx^2 + 1/dy^2 ) ));
        setValue(ie,eqn(i1+1,i2  ),    - (.5*nu*dt)*(      1/dx^2            ));
        setValue(ie,eqn(i1  ,i2+1),    - (.5*nu*dt)*(               1/dy^2   ));
      end
      end
    
    else

      for( i1=I1 )
      for( i2=I2 )
        ie = eqn(i1,i2) ; % eqn number for pt (i1,i2) 
        rx = par.rx(i1,i2,1,1);
        ry = par.rx(i1,i2,1,2);
        sx = par.rx(i1,i2,2,1);
        sy = par.rx(i1,i2,2,2);

        rxr = DJzr(par.rx,i1,i2,1,1);
        rxs = DJzs(par.rx,i1,i2,1,1);
        ryr = DJzr(par.rx,i1,i2,1,2);
        rys = DJzs(par.rx,i1,i2,1,2);
        sxr = DJzr(par.rx,i1,i2,2,1);
        sxs = DJzs(par.rx,i1,i2,2,1);
        syr = DJzr(par.rx,i1,i2,2,2);
        sys = DJzs(par.rx,i1,i2,2,2);        

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

        setValue( ie,eqn(i1-1,i2-1),    - (.5*nu*dt)*(                            +c11  ) );
        setValue( ie,eqn(i1  ,i2-1),    - (.5*nu*dt)*(            c02        -c01       ) );
        setValue( ie,eqn(i1+1,i2-1),    - (.5*nu*dt)*(                            -c11  ) );
        setValue( ie,eqn(i1-1,i2  ),    - (.5*nu*dt)*(      c20        -c10             ) );
        setValue( ie,eqn(i1  ,i2  ), 1. - (.5*nu*dt)*( -2*( c20 + c02  )                ) );
        setValue( ie,eqn(i1+1,i2  ),    - (.5*nu*dt)*(      c20        +c10             ) );
        setValue( ie,eqn(i1-1,i2+1),    - (.5*nu*dt)*(                            -c11  ) );
        setValue( ie,eqn(i1  ,i2+1),    - (.5*nu*dt)*(            c02       +c01        ) );
        setValue( ie,eqn(i1+1,i2+1),    - (.5*nu*dt)*(                            +c11  ) );        
        if( isSingular==1 ) 
          setValue(ie,Ngs,1.);
        end 
      end
      end      


    end

    % ---- Boundary conditions for the matrix ----

    dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

    for axis=1:par.nd
    for side=1:2
      % isv(1:2), is1, is2 : index shifts
      isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  
  
      [I1b,I2b]=getBoundaryIndex(side,axis,par);  

      if( par.bc(side,axis)==par.dirichlet     || ...
          par.bc(side,axis)==par.inflow      || ...
          par.bc(side,axis)==par.noSlipWall )

        % Set  u = g
        [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 

        % fprintf('impMatrix: (side,axis)=(%d,%d) : NOSLIP getAdjustedBoundaryIndex: [I1b,I2b]=[%d,%d][%d,%d]\n', side,axis, I1b(1),I1b(end), I2b(1),I2b(end)); 


        for i2=I2b
        for i1=I1b
          ie=eqn(i1,i2); 
          % fprintf('Set (i1,i2)=(%3d,%3d) ie=%4d to %10.2e\n',i1,i2,ie,1.);
          setValue(ie,ie,1.);
        end 
        end 

        [I1b,I2b]=getBoundaryIndex(side,axis,par); 
        %  --- extrapolate ghost ---
        for i2=I2b
        for i1=I1b         
          ie=eqn(i1-is1,i2-is2);  % ghost point
          setValue(ie,ie,                     1.); 
          setValue(ie,eqn(i1      ,i2      ),-3.); 
          setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
          setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 
        end
        end

      elseif( par.bc(side,axis)==par.slipWall )
        %   nv.uv   = given 
        %  (tv.uv).n = given 

        [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
        % fprintf('impMatrix: (side,axis)=(%d,%d) : SLIP   getAdjustedBoundaryIndex: [I1b,I2b]=[%d,%d][%d,%d]\n', side,axis, I1b(1),I1b(end), I2b(1),I2b(end)); 

        for i2=I2b
        for i1=I1b   
          if( (axis==1 && iuv==1) || (axis==2 && iuv==2) )
            % Dirichlet condition on this component
            ie=eqn(i1,i2);         % boundary pt
            setValue(ie,ie,1.); 
          end            
        end
        end  

        [I1b,I2b]=getBoundaryIndex(side,axis,par); 
        for i2=I2b
        for i1=I1b   
          if( (axis==1 && iuv==1) || (axis==2 && iuv==2) )
            % Extrapolate ghost
            ie=eqn(i1-is1,i2-is2); % ghost pt
            setValue(ie,ie,                     1.); 
            setValue(ie,eqn(i1      ,i2      ),-3.); 
            setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
            setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 

          else
            % Neumann condition on this component
            ie=eqn(i1-is1,i2-is2); % ghost point
            setValue(ie,eqn(i1-is1,i2-is2), 1/(2.*dxv(axis))); % ghost point 
            setValue(ie,eqn(i1+is1,i2+is2),-1/(2.*dxv(axis))); % first interior point              
          end            

        end
        end       

      elseif( par.bc(side,axis)==par.pressureInflow )
        %   tv.uv = given 
        %   nv.uv = extrapolated 

        [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index           
        for i2=I2b
        for i1=I1b   
          if( (axis==1 && iuv==2) || (axis==2 && iuv==1) )
            % Dirichlet condition on this component
            ie=eqn(i1,i2); 
            setValue(ie,ie,1.);  
          end                        
        end
        end   

        [I1b,I2b]=getBoundaryIndex(side,axis,par); 
        for i2=I2b
        for i1=I1b   
          if( (axis==1 && iuv==2) || (axis==2 && iuv==1) )
            % Extrapolate ghost
            ie=eqn(i1-is1,i2-is2); 
            setValue(ie,ie,                     1.); 
            setValue(ie,eqn(i1      ,i2      ),-3.); 
            setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
            setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 

          else
            % Extrapolate ghost
            ie=eqn(i1-is1,i2-is2); 
            setValue(ie,ie,                     1.); 
            setValue(ie,eqn(i1      ,i2      ),-3.); 
            setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
            setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.);               
          end                        
        end
        end  

      elseif( par.bc(side,axis)==par.outflow )
        % extrapolation
        % **CHECK ME** Maybe should use [1,-2,1] extrpolation
        for i2=I2b
        for i1=I1b   
          % Extrapolate ghost
          ie=eqn(i1-is1,i2-is2); 
          setValue(ie,ie,                     1.); 
          setValue(ie,eqn(i1      ,i2      ),-3.); 
          setValue(ie,eqn(i1+  is1,i2+  is2), 3.);  
          setValue(ie,eqn(i1+2*is1,i2+2*is2),-1.); 
        end
        end            

      elseif( par.bc(side,axis)==par.periodic )
        for i2=I2b
        for i1=I1b   
          % Periodic: u(iax-1,.) = u(ibx-1,.) ...
          ie=eqn(i1-is1,i2-is2); % ghost point
          setValue(ie,ie                                            , 1.); 
          setValue(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2),-1.);             
        end
        end

      else 
        fprintf('formImplicitTimeSteppingMatrix: Error: unknown bc=%d\n',par.bc(side,axis));
        pause
      end         

    end % end for axis
    end % end for side



    % ----- CORNERS ---
    if( par.bc(1,1)~=par.periodic ||  par.bc(1,2)~=par.periodic )
      % extrapolate corners along the diagonal 
      for( side1=0:1 )
        for( side2=0:1 )
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
    else
      % par.periodic in both x and y 
      for( side1=0:1 )
        for( side2=0:1 )
          is1 = 1-2*side1;
          is2 = 1-2*side2;
          i1=par.gid(side1+1,1)-is1; i2=par.gid(side2+1,2)-is2; % corner ghost point
          ie=eqn(i1,i2);  
          setValue(ie,ie                                    , 1.); 
          setValue(ie,eqn(i1+(ibx-iax)*is1,i2+(iby-iay)*is2),-1.); 
        end
      end 
    end 

    %
    % ----- FIX UP CORNERS ---
    % 
    for( side1=1:2 )
    for( side2=1:2 )

      if( par.bc(side1,1)==par.slipWall && ...
          par.bc(side2,2)==par.slipWall ) 
        % slipwall - slipwall corner: give both components
        i1=par.gid(side1,1); i2=par.gid(side2,2); % corner point
        ie=eqn(i1,i2); 
        setValue(ie,ie, 1.); 
      end

    end
    end 

    if( par.idebug>0 )
      fprintf('Optimized fill of implicit time-stepping matrix: nzzEst=%d, nzz=%d\n',nzzEst,nzz);
    end

    AA = sparse(ia(1:nzz), ja(1:nzz), aa(1:nzz), Ng, Ng); % Creates the sparse matrix     


    if( 1==0 )
      % printArray(ia,'ia','%3d');
      % printArray(ja,'ja','%3d');
      % printArray(aa,'aa','%3.0f');


      % print corner values:
      Af = full(AA);
      printArray(Af','Af','%3.0f');
      % [I1,I2] = getIndex(par.gid)
      % for i2=I2
      % for i1=I1
      %   fprintf('(i1,i2)=(%d,%d)');
      %   ie = eqn(i1,i2);
      %   for( je=1:Ng )
      %     fprint
      % end 
      % end 
      for( side1=1:2 )
      for( side2=1:2 )
        i1=par.gid(side1,1); i2=par.gid(side2,2); % corner point
        ie = eqn(i1,i2);
        fprintf('impMatrix: iuv=%d: corner (i1,i2)=(%3d,%3d) ie=%3d  AA(%3d,%3d)=%e\n',iuv,i1,i2,ie,ie,ie,Af(ie,ie));
      end
      end       
    end

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
   % fprintf('setValue: (ii,jj,aa)=(%3d,%3d,%10.3e)\n',ii,jj,val);
  end  

return
end