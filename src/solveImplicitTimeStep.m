%
%  Solve the Implicit-explicit time-stepping system 
%
% NOTE:
%  (unp1,vnp1) : holds the right-hand-side to the implicit equations at interior equations
%
function [unp1,vnp1,par] = solveImplicitTimeStep( unp1,vnp1,tnp1, gf,cur, par )


  if( par.combinedImplicitSolverNeeded )
    % The next function handles the case when the velcoity components are coupled
    [unp1,vnp1,par] = solveImplicitTimeStepCombined( unp1,vnp1,tnp1, gf,cur, par );
    return
  end

  if( par.factorImplicitMatrix )
     par = formImplicitTimeSteppingMatrix( tnp1, par.dt, gf,cur, par );
     par.factorImplicitMatrix=0; 
  end 

  cpu0 = cputime;

  Ngx = par.Ngx;
  Ngy = par.Ngy;

  % iax = par.iax;  ibx = par.ibx; 
  % iay = par.iay;  iby = par.iby; 

  iax = par.gid(1,1); ibx=par.gid(2,1);
  iay = par.gid(1,2); iby=par.gid(2,2);    

  numGhost = par.numGhost;  

  % --- declare difference operators ---
  % declareOperators;

  % convert (ix,iy) to equation number in the matrix:  
  eqn = @(ix,iy)  1 + ix-1 + Ngx*( iy-1 );
  
  Ng = Ngx*Ngy; % total number of grid points
  
  rhsu = zeros(Ng,1);     % right-hand-side
  rhsv = zeros(Ng,1);     % right-hand-side

  % copy grid-function to vector solution
  for( ix=iax:ibx )
  for( iy=iay:iby )
    ie = eqn(ix,iy) ; % eqn number for pt (ix,iy) 
    rhsu(ie) = unp1(ix,iy);
    rhsv(ie) = vnp1(ix,iy);
  end
  end

  % Boundary conditions for implicit equations
  if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'known')  ) manufacturedSolution=1; else manufacturedSolution=0; end

  for axis=1:2
  for side=1:2
    isv(1)=0; isv(2)=0; isv(axis)=1-2*(side-1);  is1=isv(1); is2=isv(2);
    mbc = side+2*(axis-1);  % pointer into gu and gv arrays
    [I1b,I2b]=getBoundaryIndex(side,axis,par);

    if( par.bc(side,axis)==par.dirichlet  || ...
        par.bc(side,axis)==par.noSlipWall || ...
        par.bc(side,axis)==par.inflow )

      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      % u = given
      for i2=I2b
      for i1=I1b 
        ie=eqn(i1,i2); % boundary point 
        if( manufacturedSolution )
          rhsu(ie) = par.ue(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          rhsv(ie) = par.ve(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);

          % fprintf(' imp rhs: (i1,i2)=(%3d,%3d) ie=%3d ue=%10.3e\n',i1,i2,ie,rhsu(ie));
        else
          rhsu(ie) = 0;
          rhsv(ie) = 0;
        end
      end
      end

    elseif( par.bc(side,axis)==par.slipWall )
      % -- slip wall  --
      %   nv.uv   = given 
      %  (tv.uv).n = given
      % (n1,n2) = outward normal 
      n1 = -is1;
      n2 = -is2; 

      [I1b,I2b]=getAdjustedBoundaryIndex(side,axis,par);  % adjusted index 
      for i2=I2b
      for i1=I1b 

        ie=eqn(i1,i2); % boundary point
        if( axis==1 )
          % left/right : set u 
          if( manufacturedSolution )
            rhsu(ie) = par.gu{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsu(ie)=0.;
          end 
        else
          % bottom/top: set v 
          if( manufacturedSolution )
            rhsv(ie) = par.gv{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsv(ie)=0.; 
          end
        end

      end % end for i1
      end % end for i2 

      [I1b,I2b]=getBoundaryIndex(side,axis,par);
      for i2=I2b
      for i1=I1b       
        % (tv.uv).n = given
        ie=eqn(i1-is1,i2-is2); % ghost point 
        if( axis==1 )
           % left/right : v.x = 
          if( manufacturedSolution )
            rhsv(ie) = n1*par.vex(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsv(ie)=0.;
          end
        else
          % bottom top: u.y = 
          if( manufacturedSolution )
            rhsu(ie) = n2*par.uey(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
          else
            rhsu(ie)=0.;
          end               
        end 

      end % end for i1
      end % end for i2 

    elseif( par.bc(side,axis)==par.pressureInflow )

      %   tv.uv = given 
      %   nv.uv = extrapolated 
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

    elseif( par.bc(side,axis)==par.traction )

      fprintf('solveImp: finish me for traction bc\n');
      error('error');


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
      ie = eqn(i1,i2);
      if( manufacturedSolution )
        rhsu(ie) = par.gu{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
      else
        rhsu(ie)=0.;
      end 
      if( manufacturedSolution )
        rhsv(ie) = par.gv{mbc}(gf{cur}.x(i1,i2,1),gf{cur}.x(i1,i2,2),tnp1);
      else
        rhsv(ie)=0.; 
      end
    end
  end
  end 

  % --- solve ---
  rhsu = par.dAimp{1}\rhsu; 
  if( par.multipleImplicitSolversNeeded ) 
    rhsv = par.dAimp{2}\rhsv;    
  else 
    rhsv = par.dAimp{1}\rhsv;   % we can use the same implicit solver for u and v
  end

  % copy vector solution to the grid function solution
  I1g = iax-numGhost:ibx+numGhost; I2g=iay-numGhost:iby+numGhost; % include ghost points 
  for( ix=I1g )
  for( iy=I2g )
    ie = eqn(ix,iy); % eqn number for pt (ix,iy) 
     unp1(ix,iy)=rhsu(ie);
     vnp1(ix,iy)=rhsv(ie);
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
