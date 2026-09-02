%
%  Solve the Implicit-explicit time-stepping system when the velocity components are coupled
%
% NOTE:
%  (unp1,vnp1) : holds the right-hand-side to the implicit equations at interior equations
%
function [unp1,vnp1,par] = solveImplicitTimeStepCombined( unp1,vnp1,tnp1, gf,cur, par )


  if( par.factorImplicitMatrix )
     par = formImplicitTimeSteppingMatrixCombined( tnp1, par.dt, gf,cur, par );
     par.factorImplicitMatrix=0; 
  end 

  cpu0 = cputime;

  nu = par.nu;
  mu = par.mu;

  Ngx = par.Ngx;
  Ngy = par.Ngy;

  uc = par.uc; % u component index (=1)
  vc = par.vc; % v component index (=2)

  % iax = par.iax;  ibx = par.ibx; 
  % iay = par.iay;  iby = par.iby; 

  iax = par.gid(1,1); ibx=par.gid(2,1);
  iay = par.gid(1,2); iby=par.gid(2,2);    

  numGhost = par.numGhost;  

  % convert (i1,i2,ic) to equation number in the matrix:  (ic=component number, 1,2)
  eqn = @(i1,i2,ic)  1 + i1-1 + Ngx*( i2-1 + Ngy*(ic-1) );

    
  Ng = Ngx*Ngy; % total number of grid points
  Ngc = Ng*par.nd; % dimension of the combined matrix 
  
  rhs = zeros(Ngc,1);     % right-hand-side

  % copy grid-function to the combined vector solution
  for( ix=iax:ibx )
  for( iy=iay:iby )
    ieu = eqn(ix,iy,uc); 
    iev = eqn(ix,iy,vc); 
    rhs(ieu) = unp1(ix,iy);
    rhs(iev) = vnp1(ix,iy);
  end
  end

  % Boundary conditions for implicit equations
  if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'known')  ) manufacturedSolution=1; else manufacturedSolution=0; end

  for axis=1:2
  for side=1:2
    isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);
    mbc = side+2*(axis-1);  % pointer into gu and gv arrays
    [I1b,I2b]=getBoundaryIndex(side,axis,par);


    if( par.bc(side,axis)==par.slipWall )
      % -- slip wall  : BOUNDARY POINTS --
      % NOTE: RHS values for GHOST POINTS are set below with traction BC 

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


      % fprintf('solveImpCombined: finish me - slip wall\n'); pause

      % [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      gu=0; gv=0; 
      for i2=I2b
      for i1=I1b 
        ieu = eqn(i1,i2,uc);   % boundary point 
        iev = eqn(i1,i2,vc);        
        if( par.isCartesian )
          n1 = -is1;
          n2 = -is2;
        else
          [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );
        end  
        if( manufacturedSolution )
          gu = par.ue(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          gv = par.ve(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
        end         
        %   n1^2 * u + n1*n2*v + (1-n^2)*[ I - nu*dt/2 Delta_h ] u  - n1*n2 * [ I - nu*dt/2 Delta_h ] v = RHS
        rhs(ieu) =  n1^2*gu + n1*n2*gv + (1-n1^2)*unp1(i1,i2)   - n1*n2*vnp1(i1,i2);
        %   n1*n2 * u + n2^2*v - n1*n2*[ I - nu*dt/2 Delta_h ] u  +(1-n2^2) * [ I - nu*dt/2 Delta_h ] v = RHS
        rhs(iev) = n1*n2*gu +  n2^2*gv     -n1*n2*unp1(i1,i2) +(1-n2^2)*vnp1(i1,i2);

      end % end for i1
      end % end for i2     
    end % end if slipWall

    if( par.bc(side,axis)==par.dirichlet  || ...
        par.bc(side,axis)==par.noSlipWall || ...
        par.bc(side,axis)==par.inflow )

      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      % u = given
      for i2=I2b
      for i1=I1b 
        ieu = eqn(i1,i2,uc);   % boundary point 
        iev = eqn(i1,i2,vc);         
        if( manufacturedSolution )
          rhs(ieu) = par.ue(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          rhs(iev) = par.ve(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);

          % fprintf(' imp rhs: (i1,i2)=(%3d,%3d) ie=%3d ue=%10.3e\n',i1,i2,ie,rhsu(ie));
        else
          rhs(ieu) = 0;
          rhs(iev) = 0;
        end
      end
      end



    elseif( par.bc(side,axis)==par.pressureInflow )

      %   tv.uv = given 
      %   nv.uv = extrapolated 

      fprintf('solveImpCombined: finish me - pressureInflow\n'); pause

      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      for i2=I2b
      for i1=I1b 
        ie=eqn(i1,i2); % boundary point
        if( axis==1 )
          % left/right : set v  
          if( manufacturedSolution )
            rhsv(ie) = par.gv{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsv(ie)=0.; 
          end
        else
          % bottom/top: set u 
          if( manufacturedSolution )
            rhsu(ie) = par.gu{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsu(ie)=0.; 
          end
        end
      end
      end

    elseif( par.bc(side,axis)==par.traction || par.bc(side,axis)==par.slipWall )

      % fprintf('solveImpCombined: finish me for traction bc\n');
      % error('error');
      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      for i2=I2b
      for i1=I1b 
        ieu = eqn(i1-is1,i2-is2,uc);   % u ghost 
        iev = eqn(i1-is1,i2-is2,vc);   % v ghost       
        if( manufacturedSolution )
          if( par.isCartesian )
            n1 = -is1;
            n2 = -is2;    
          else 
            [n1,n2,rxNorm] = getBoundaryNormal( side,axis,i1,i2, gf,cur, par );
          end

          % rx = gf{cur}.rx(i1,i2,1,1); sx = gf{cur}.rx(i1,i2,2,1);
          % ry = gf{cur}.rx(i1,i2,1,2); sy = gf{cur}.rx(i1,i2,2,2);

          ux = par.uex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          uy = par.uey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          vx = par.vex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          vy = par.vey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);    

          tauvnv1 =  mu*( 2*ux*n1 + (uy+vx)*n2 );
          tauvnv2 =  mu*( 2*vy*n2 + (uy+vx)*n1 );   
          rhs(ieu) = ( ux+vy )*n1 + (1-n1^2)*( tauvnv1 )    -n1*n2*( tauvnv2 );
          rhs(iev) = ( ux+vy )*n2 +   -n1*n2*( tauvnv1 ) +(1-n2^2)*( tauvnv2 );             

          % fprintf(' imp rhs: (i1,i2)=(%3d,%3d) ie=%3d ue=%10.3e\n',i1,i2,ie,rhsu(ie));
        else
          rhs(ieu) = 0;
          rhs(iev) = 0;
        end
      end
      end      


    elseif( par.bc(side,axis)==par.outflow )
      % velocity is extrapolated -- do nothing here

    elseif( par.bc(side,axis)==par.periodic )
      % do nothing here 

    else
      fprintf('solveImplicitTimeStep: ERROR: unexpected bc=%d\n',par.bc(side,axis));
      pause;

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
      ieu = eqn(i1,i2,uc);
      iev = eqn(i1,i2,vc);
      if( manufacturedSolution )
        rhsu(ieu) = par.gu{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
      else
        rhsu(ieu)=0.;
      end 
      if( manufacturedSolution )
        rhsv(iev) = par.gv{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
      else
        rhsv(iev)=0.; 
      end
    end
  end
  end 

  % --- solve ---
  rhs = par.dAimp\rhs; 

  % copy vector solution to the grid function solution
  I1g = iax-numGhost:ibx+numGhost; I2g=iay-numGhost:iby+numGhost; % include ghost points 
  for( ix=I1g )
  for( iy=I2g )
     ieu = eqn(ix,iy,uc); 
     iev = eqn(ix,iy,vc); 
     unp1(ix,iy)=rhs(ieu);
     vnp1(ix,iy)=rhs(iev);
  end
  end

  if( 1==0 || par.idebug>3 )
    % -- check the error --
    pnp1 = zeros(par.Ngx,par.Ngy);
    fprintf('ERRORS AFTER IMPLICIT SOLVE\n');
    par = plotSolution( tnp1,unp1,vnp1,pnp1, par);
    pause
  end

  par.cpuImplicit = par.cpuImplicit + cputime - cpu0;
end
