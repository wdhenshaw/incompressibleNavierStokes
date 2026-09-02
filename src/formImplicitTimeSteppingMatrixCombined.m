%
% Form the implicit time-stepping matrix for the combined case when all velocity components are coupled
%
function par = formImplicitTimeSteppingMatrixCombined( t,dt, gf,cur, par )


  cpu0=cputime;

  nu = par.nu;
  mu = par.mu;

  uc = par.uc;
  vc = par.vc;

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

  % derivatives of entries in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  % Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r
  % Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s  

 
  Ng = Ngx*Ngy;    % total number of grid points 
  Ngc = Ng*par.nd; % dimension of the combined matrix 


  if(  par.plotOption>=0 &&  ( 1==1 || mod(floor(par.idebug/2),2)==1 || (par.gridMotion~=par.noMotion && t<=2*dt) ) )
    fprintf('...FORM the COMBINED IMPLICIT TIME-STEPPING MATRIX: Ng=%d Ngc=%d t=%9.3e cur=%d\n',Ng,Ngc,t,cur); 
  end 

  % convert (i1,i2,ic) to equation number in the matrix:  (ic=component number, 1,2)
  eqn = @(i1,i2,ic)  1 + i1-1 + Ngx*( i2-1 + Ngy*(ic-1) );

  nzzEst = Ngc*5; % estimated number of non-zeros
  ia=zeros(nzzEst,1); ja=zeros(nzzEst,1); aa=zeros(nzzEst,1);

  nzz=0; % counts non-zeros
  [I1a,I2a] = getIndexInterior( par.gid,-1,par );  % FIX ME 

  if( par.isCartesian )

    % --- Cartesian grid ---

    for( ic=1:par.nd ) % velocity components
      for( i2=I2a )
      for( i1=I1a )
        ie = eqn(i1,i2,ic) ; % eqn number for pt (i1,i2,ic) 
        setValue(ie,eqn(i1  ,i2-1,ic),    - (.5*nu*dt)*(               1/dy^2   ));
        setValue(ie,eqn(i1-1,i2  ,ic),    - (.5*nu*dt)*(      1/dx^2            ));
        setValue(ie,eqn(i1  ,i2  ,ic), 1. - (.5*nu*dt)*( -2*( 1/dx^2 + 1/dy^2 ) ));
        setValue(ie,eqn(i1+1,i2  ,ic),    - (.5*nu*dt)*(      1/dx^2            ));
        setValue(ie,eqn(i1  ,i2+1,ic),    - (.5*nu*dt)*(               1/dy^2   ));
      end
      end
    end % end for ic 
  
  else
    % ---- CURVLINEAR -----

    for( ic=1:par.nd ) % velocity components -- would be more efficient to move this loop down 
      for( i1=I1a )
      for( i2=I2a )
        ie = eqn(i1,i2,ic) ; % eqn number for pt (i1,i2) 
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

        setValue( ie,eqn(i1-1,i2-1,ic),    - (.5*nu*dt)*(                            +c11  ) );
        setValue( ie,eqn(i1  ,i2-1,ic),    - (.5*nu*dt)*(            c02        -c01       ) );
        setValue( ie,eqn(i1+1,i2-1,ic),    - (.5*nu*dt)*(                            -c11  ) );
        setValue( ie,eqn(i1-1,i2  ,ic),    - (.5*nu*dt)*(      c20        -c10             ) );
        setValue( ie,eqn(i1  ,i2  ,ic), 1. - (.5*nu*dt)*( -2*( c20 + c02  )                ) );
        setValue( ie,eqn(i1+1,i2  ,ic),    - (.5*nu*dt)*(      c20        +c10             ) );
        setValue( ie,eqn(i1-1,i2+1,ic),    - (.5*nu*dt)*(                            -c11  ) );
        setValue( ie,eqn(i1  ,i2+1,ic),    - (.5*nu*dt)*(            c02       +c01        ) );
        setValue( ie,eqn(i1+1,i2+1,ic),    - (.5*nu*dt)*(                            +c11  ) );        
      end
      end     
    end % end for ic  


  end    

  
  % ---- Boundary conditions for the matrix ----

  dxv(1)=dx; dxv(2)=dy; % save grid spacing in an array

  for axis=1:par.nd
  for side=1:2
    % isv(1:2), is1, is2 : index shifts
    isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);  

    [I1b,I2b]=getBoundaryIndex(side,axis,par);  

 
    if( par.bc(side,axis)==par.slipWall )
      % ---- SPECIAL CASE : slip wall -----
      % The slip wall condition on the boundary
      %   nv.uv   = given 
      % is combined with the tangential component of the interior PDE
      % 
      %   (nv.uv) nv +  tv tv^T  [ I - nu*dt/2 Delta_h ]
      % or 
      %    nv nv^T  +  (I-nv nv^T) [ I - nu*dt/2 Delta_h ]
      % 
      %  nv nv^T = [ n1^2  n1*n2 ],    I - nv nv^T = [ 1-n1^2  -n1*n2]
      %            [ n1*n2 n2^2  ]                   [ -n1*n2  1-n2^2 ]

      % NOTE: equations for ghost points are set below with traction BC 

      if( par.isCartesian )
        % There is a bug somwhere, maybe in another routine, that prevents slip wall working with a Cartesian grid and the combined matrix 

        fprintf('\n formImpTimeSteppingCombined:ERROR: slip wall does not work for a Cartesian grid\n');
        fprintf('There is a bug somwhere, maybe in another routine, that prevents slip wall working with a Cartesian grid and the combined matrix\n\n');
        error('error');
      end

      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 

      for i2=I2b
      for i1=I1b   
        
        ieu=eqn(i1,i2,uc);         % u boundary pt
        iev=eqn(i1,i2,vc);         % v boundary pt

        if( par.isCartesian )
          n1 = -is1;
          n2 = -is2;
          c20 = 1/(dx^2); c02=1/(dy^2); c11=0; c10=0; c01=0; 
        else 
          [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );
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
        end 

        % Equation 1: 
        %   n1^2 * u + n1*n2*v + (1-n^2)*[ I - nu*dt/2 Delta_h ] u  - n1*n2 * [ I - nu*dt/2 Delta_h ] v
        setValue( ieu,eqn(i1-1,i2-1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(                            +c11  )) );
        setValue( ieu,eqn(i1  ,i2-1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(            c02        -c01       )) );
        setValue( ieu,eqn(i1+1,i2-1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( ieu,eqn(i1-1,i2  ,uc),         (1-n1^2)*(    - (.5*nu*dt)*(      c20        -c10             )) );
        setValue( ieu,eqn(i1  ,i2  ,uc), n1^2 +  (1-n1^2)*( 1. - (.5*nu*dt)*( -2*( c20 + c02  )                )) );
        setValue( ieu,eqn(i1+1,i2  ,uc),         (1-n1^2)*(    - (.5*nu*dt)*(      c20        +c10             )) );
        setValue( ieu,eqn(i1-1,i2+1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( ieu,eqn(i1  ,i2+1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(            c02       +c01        )) );
        setValue( ieu,eqn(i1+1,i2+1,uc),         (1-n1^2)*(    - (.5*nu*dt)*(                            +c11  )) );             

        setValue( ieu,eqn(i1-1,i2-1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(                            +c11  )) );
        setValue( ieu,eqn(i1  ,i2-1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(            c02        -c01       )) );
        setValue( ieu,eqn(i1+1,i2-1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( ieu,eqn(i1-1,i2  ,vc),        - (n1*n2)*(    - (.5*nu*dt)*(      c20        -c10             )) );
        setValue( ieu,eqn(i1  ,i2  ,vc), n1*n2  - (n1*n2)*( 1. - (.5*nu*dt)*( -2*( c20 + c02  )                )) );
        setValue( ieu,eqn(i1+1,i2  ,vc),        - (n1*n2)*(    - (.5*nu*dt)*(      c20        +c10             )) );
        setValue( ieu,eqn(i1-1,i2+1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( ieu,eqn(i1  ,i2+1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(            c02       +c01        )) );
        setValue( ieu,eqn(i1+1,i2+1,vc),        - (n1*n2)*(    - (.5*nu*dt)*(                            +c11  )) );             

        % Equation 2:
        % n1*n2 * u + n2^2*v - n1*n2*[ I - nu*dt/2 Delta_h ] u  +(1-n2^2) * [ I - nu*dt/2 Delta_h ] v 
        setValue( iev,eqn(i1-1,i2-1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(                            +c11  )) );
        setValue( iev,eqn(i1  ,i2-1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(            c02        -c01       )) );
        setValue( iev,eqn(i1+1,i2-1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( iev,eqn(i1-1,i2  ,vc),         (1-n2^2)*(    - (.5*nu*dt)*(      c20        -c10             )) );
        setValue( iev,eqn(i1  ,i2  ,vc), n2^2 +  (1-n2^2)*( 1. - (.5*nu*dt)*( -2*( c20 + c02  )                )) );
        setValue( iev,eqn(i1+1,i2  ,vc),         (1-n2^2)*(    - (.5*nu*dt)*(      c20        +c10             )) );
        setValue( iev,eqn(i1-1,i2+1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( iev,eqn(i1  ,i2+1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(            c02       +c01        )) );
        setValue( iev,eqn(i1+1,i2+1,vc),         (1-n2^2)*(    - (.5*nu*dt)*(                            +c11  )) );             

        setValue( iev,eqn(i1-1,i2-1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(                            +c11  )) );
        setValue( iev,eqn(i1  ,i2-1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(            c02        -c01       )) );
        setValue( iev,eqn(i1+1,i2-1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( iev,eqn(i1-1,i2  ,uc),        - (n1*n2)*(    - (.5*nu*dt)*(      c20        -c10             )) );
        setValue( iev,eqn(i1  ,i2  ,uc), n1*n2  - (n1*n2)*( 1. - (.5*nu*dt)*( -2*( c20 + c02  )                )) );
        setValue( iev,eqn(i1+1,i2  ,uc),        - (n1*n2)*(    - (.5*nu*dt)*(      c20        +c10             )) );
        setValue( iev,eqn(i1-1,i2+1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(                            -c11  )) );
        setValue( iev,eqn(i1  ,i2+1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(            c02       +c01        )) );
        setValue( iev,eqn(i1+1,i2+1,uc),        - (n1*n2)*(    - (.5*nu*dt)*(                            +c11  )) );      

      end
      end  

 
    end % end if slip wall     


    if( par.bc(side,axis)==par.dirichlet     || ...
        par.bc(side,axis)==par.inflow        || ...
        par.bc(side,axis)==par.noSlipWall )

      % Set  u = g
      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 

      % fprintf('impMatrix: (side,axis)=(%d,%d) : NOSLIP getAdjustedBoundaryIndex: [I1b,I2b]=[%d,%d][%d,%d]\n', side,axis, I1b(1),I1b(end), I2b(1),I2b(end)); 

      for( ic=1:par.nd ) 
        for i2=I2b
        for i1=I1b
          ie=eqn(i1,i2,ic); 
          setValue(ie,ie,1.);
        end 
        end 

        [I1b,I2b]=getBoundaryIndex(side,axis,par); 
        %  --- extrapolate ghost ---
        for i2=I2b
        for i1=I1b         
          ie=eqn(i1-is1,i2-is2,ic);  % ghost point
          setValue(ie,ie,                        1.); 
          setValue(ie,eqn(i1      ,i2      ,ic),-3.); 
          setValue(ie,eqn(i1+  is1,i2+  is2,ic), 3.);  
          setValue(ie,eqn(i1+2*is1,i2+2*is2,ic),-1.); 
        end
        end
      end % end for ic 


    elseif( par.bc(side,axis)==par.pressureInflow )
      %   tv.uv = given 
      %   nv.uv = extrapolated 

      fprintf('fillCombined: finish me - pressureInflow\n'); pause

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
      for( ic=1:par.nd ) 
        for i2=I2b
        for i1=I1b   
          % Extrapolate ghost
          ie=eqn(i1-is1,i2-is2,ic); 
          setValue(ie,ie,                        1.); 
          setValue(ie,eqn(i1      ,i2      ,ic),-3.); 
          setValue(ie,eqn(i1+  is1,i2+  is2,ic), 3.);  
          setValue(ie,eqn(i1+2*is1,i2+2*is2,ic),-1.); 
        end
        end  
      end % end for ic 

    elseif( par.bc(side,axis)==par.traction || par.bc(side,axis)==par.slipWall )

      % tauv*nv : 
      %   tau nv = mu [ 2*ux uy+vx ][ n1 ]
      %               [ uy+vx 2*vy ]  n2 ]           
      % tauvnv1 =  mu*( 2*ux*n1 + (uy+vx)*n2 );
      % tauvnv2 =  mu*( 2*vy*n2 + (uy+vx)*n1 );
      %   I - nv nv^T = [ 1-n1^2 -n1*n2 ]
      %                 [ -n1*n2 1-n2^2 ]
      % Traction BCs 
      %  A: ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )    -n1*n2*( tauvnv2 )
      %  B: ( ux+vy )*n2 +   -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 )

      % fprintf('fillImpMat: finish me for traction bc\n');
      % error('error');
      [I1b,I2b]=getBoundaryIndex(side,axis,par);  % use adjusted index ??

      for i2=I2b
        for i1=I1b   
          % fprintf('formImpCombined: traction BC (i1,i2)=(%3d,%3d)\n',i1,i2);
          if( par.isCartesian )
            n1 = -is1;
            n2 = -is2;    
            rx=1; ry=0; sx=0; sy=1; 
          else      
            [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );

            rx = gf{cur}.rx(i1,i2,1,1);
            ry = gf{cur}.rx(i1,i2,1,2);
            sx = gf{cur}.rx(i1,i2,2,1);
            sy = gf{cur}.rx(i1,i2,2,2);
          end

          ieu=eqn(i1-is1,i2-is2,uc);         % u ghost pt
          iev=eqn(i1-is1,i2-is2,vc);         % v ghost pt

          % --- Equation A: ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )    -n1*n2*( tauvnv2 ) = 0  ---
          cru = n1*( rx/(2*dr) ) + (1-n1^2)*mu*( 2*n1*rx/(2*dr) + n2*ry/(2*dr) ) -n1*n2*mu*(                   n1*ry/(2*dr) ); % coeff of u(i1+1,i2)
          crv = n1*( ry/(2*dr) ) + (1-n1^2)*mu*(                  n2*rx/(2*dr) ) -n1*n2*mu*(  2*n2*ry/(2*dr) + n1*rx/(2*dr) ); % coeff of v(i1+1,i2)

          csu = n1*( sx/(2*ds) ) + (1-n1^2)*mu*( 2*n1*sx/(2*ds) + n2*sy/(2*ds) ) -n1*n2*mu*(                   n1*sy/(2*ds) ); % coeff of u(i1,i2+1)
          csv = n1*( sy/(2*ds) ) + (1-n1^2)*mu*(                  n2*sx/(2*ds) ) -n1*n2*mu*(  2*n2*sy/(2*ds) + n1*sx/(2*ds) ); % coeff of v(i1,i2+1)

          setValue( ieu,eqn(i1  ,i2-1,uc),-csu );
          setValue( ieu,eqn(i1  ,i2-1,vc),-csv );          
          setValue( ieu,eqn(i1-1,i2  ,vc),-crv );
          setValue( ieu,eqn(i1-1,i2  ,uc),-cru );
          setValue( ieu,eqn(i1+1,i2  ,uc), cru );
          setValue( ieu,eqn(i1+1,i2  ,vc), crv );
          setValue( ieu,eqn(i1  ,i2+1,uc), csu );
          setValue( ieu,eqn(i1  ,i2+1,vc), csv );

          % --- Equation B: ( ux+vy )*n2 +   -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 ) ---
          cru = n2*( rx/(2*dr) ) -n1*n2*mu*( 2*n1*rx/(2*dr) + n2*ry/(2*dr) ) +(1-n2^2)*mu*(                   n1*ry/(2*dr) ); % coeff of u(i1+1,i2)
          crv = n2*( ry/(2*dr) ) -n1*n2*mu*(                  n2*rx/(2*dr) ) +(1-n2^2)*mu*(  2*n2*ry/(2*dr) + n1*rx/(2*dr) ); % coeff of v(i1+1,i2)
          csu = n2*( sx/(2*ds) ) -n1*n2*mu*( 2*n1*sx/(2*ds) + n2*sy/(2*ds) ) +(1-n2^2)*mu*(                   n1*sy/(2*ds) ); % coeff of u(i1,i2+1)
          csv = n2*( sy/(2*ds) ) -n1*n2*mu*(                  n2*sx/(2*ds) ) +(1-n2^2)*mu*(  2*n2*sy/(2*ds) + n1*sx/(2*ds) ); % coeff of v(i1,i2+1)

          setValue( iev,eqn(i1  ,i2-1,uc),-csu );
          setValue( iev,eqn(i1  ,i2-1,vc),-csv );          
          setValue( iev,eqn(i1-1,i2  ,vc),-crv );
          setValue( iev,eqn(i1-1,i2  ,uc),-cru );
          setValue( iev,eqn(i1+1,i2  ,uc), cru );
          setValue( iev,eqn(i1+1,i2  ,vc), crv );
          setValue( iev,eqn(i1  ,i2+1,uc), csu );
          setValue( iev,eqn(i1  ,i2+1,vc), csv );


        end
      end        


    elseif( par.bc(side,axis)==par.periodic )
      % done below now

      % for( ic=1:par.nd ) 
      %   for i2=I2b
      %   for i1=I1b   
      %     % Periodic: u(iax-1,.) = u(ibx-1,.) ...
      %     ie=eqn(i1-is1,i2-is2,ic); % ghost point
      %     setValue(ie,ie                                              , 1.); 
      %     setValue(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2,ic),-1.);             
      %   end
      %   end
      % end 

    else 
      fprintf('formImplicitTimeSteppingMatrixCombined: Error: unknown bc=%d\n',par.bc(side,axis));
      pause
    end         


    if( par.bc(side,axis)==par.periodic )
      % ---- periodic boundary conditions ---
      %  include ghost points 
      if( par.bc(side,axis)==par.periodic )
        if( axis==1 )
          I1b=par.gid(side,axis);
          I2b=par.dim(1,2):par.dim(2,2);
        elseif( axis==2 )
          I1b=par.dim(1,1):par.dim(2,1);  
          I2b=par.gid(side,axis);
        end        
        for( ic=1:par.nd ) 
          for i2=I2b
          for i1=I1b   
            % Periodic: u(iax-1,.) = u(ibx-1,.) ...
            ie=eqn(i1-is1,i2-is2,ic); % ghost point
            setValue(ie,ie                                              , 1.); 
            setValue(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2,ic),-1.);             
          end
          end
        end 
      end
    end 

  end % end for axis
  end % end for side



  % ----- CORNERS ---
  if( par.bc(1,1)~=par.periodic &&  par.bc(1,2)~=par.periodic )
    % extrapolate corners along the diagonal
    for( ic=1:par.nd )  
      for( side1=0:1 )
        for( side2=0:1 )
          is1 = 1-2*side1;
          is2 = 1-2*side2;
          i1=par.gid(side1+1,1)-is1; i2=par.gid(side2+1,2)-is2; % corner ghost point
          ie=eqn(i1,i2,ic); 
          setValue(ie,ie                       , 1.); 
          setValue(ie,eqn(i1+  is1,i2+  is2,ic),-3.); 
          setValue(ie,eqn(i1+2*is1,i2+2*is2,ic), 3.);  
          setValue(ie,eqn(i1+3*is1,i2+3*is2,ic),-1.); 
        end
      end
    end 
    % else
    % % par.periodic in both x and y 
    % for( ic=1:par.nd )  
    %   for( side1=0:1 )
    %     for( side2=0:1 )
    %       is1 = 1-2*side1;
    %       is2 = 1-2*side2;
    %       i1=par.gid(side1+1,1)-is1; i2=par.gid(side2+1,2)-is2; % corner ghost point
    %       ie=eqn(i1,i2,ic);  
    %       setValue(ie,ie                                       , 1.); 
    %       setValue(ie,eqn(i1+(ibx-iax)*is1,i2+(iby-iay)*is2,ic),-1.); 
    %     end
    %   end 
    % end
  end 

  %
  % ----- FIX UP CORNERS ---
  % 
  for( ic=1:par.nd ) 
    for( side1=1:2 )
    for( side2=1:2 )

      if( par.bc(side1,1)==par.slipWall && ...
          par.bc(side2,2)==par.slipWall ) 
        % slipwall - slipwall corner: give both components
        i1=par.gid(side1,1); i2=par.gid(side2,2); % corner point
        ie=eqn(i1,i2,ic); 
        setValue(ie,ie, 1.); 
      end

    end
    end 
  end 

  % --- CHECK FOR CORNER CASES WE DO NOT TREAT YET ---
  for( side1=1:2 )
    for( side2=1:2 )
      if( par.bc(side1,1)==par.traction && ...
          par.bc(side2,2)==par.traction ) 
        fprintf('\n ***** formImplicitTimeSteppingMatrixCombined:ERROR -- a traction-traction corner is not implemented yet ***\n\n');
        error('ERROR')
      elseif( par.bc(side1,1)==par.slipWall && ...
              par.bc(side2,2)==par.slipWall ) 
        fprintf('\n ***** formImplicitTimeSteppingMatrixCombined:ERROR -- a slipWall-slipWall corner is not implemented yet ***\n\n');
        error('ERROR')
      end 
    end
  end

  if(  par.plotOption>=0 && ( mod(floor(par.idebug/2),2)==1 || (par.gridMotion~=par.noMotion && t<=2*dt) ) )
    fprintf('Optimized fill of implicit time-stepping matrix: nzzEst=%d, nzz=%d\n',nzzEst,nzz);
  end



  AA = sparse(ia(1:nzz), ja(1:nzz), aa(1:nzz), Ngc, Ngc); % Creates the sparse matrix     


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
    % for( side1=1:2 )
    % for( side2=1:2 )
    %   i1=par.gid(side1,1); i2=par.gid(side2,2); % corner point
    %   ie = eqn(i1,i2);
    %   fprintf('impMatrix: iuv=%d: corner (i1,i2)=(%3d,%3d) ie=%3d  AA(%3d,%3d)=%e\n',iuv,i1,i2,ie,ie,ie,Af(ie,ie));
    % end
    % end       
  end

  if( mod(floor(par.idebug/2),2)==1 )
    fprintf('+++ cond(AA)=%10.2e\n',condest(AA));
    figure(3); 
    spy(AA);
    pause
  end 

  cpuFillImpMatrix = cputime - cpu0;

  cpu0 = cputime;
  par.dAimp =decomposition(AA);   % factor the matrix: A=LU  

  cpuFactor = cputime-cpu0;
  par.cpuFactorImpMatrix = par.cpuFactorImpMatrix + cpuFactor;

  if( par.plotOption>=0 && (mod(floor(par.idebug/2),2)==1 || (par.gridMotion~=par.noMotion && t<=2*dt) ) )
    fprintf('formImplicitMatrixCombined: time to fill matrix = %8.2e(s), factor matrix = %8.2e(s)\n',cpuFillImpMatrix,cpuFactor);
  end

  % fprintf('formImplicitMatrixCombined: stop here for now ...\n'); 
  % pause;
  % pause;

  % function to fill matrix entries
  function setValue( ii,jj,val )
   nzz=nzz+1;
   ia(nzz)=ii; ja(nzz)=jj; aa(nzz)=val; 
   % fprintf('setValue: (ii,jj,aa)=(%3d,%3d,%10.3e)\n',ii,jj,val);
  end  

return
end