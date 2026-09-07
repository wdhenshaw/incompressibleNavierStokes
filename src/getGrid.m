%
%  Evaluate the grid, metrics and gridVelocity 
% 
function [gf,par] = getGrid( t, gf,cur, par )


  if( ~isfield(gf{cur},'x') )
  	if( par.idebug>0 )
  	  fprintf('*** Create the grid for cur=%d t=%9.3e ***\n',cur,t);
  	end

    % first time: create the x and rx arrays as needed
	  gf{cur}.x = par.x(:,:,:);
	  if( ~strcmp(par.map,'Cartesian') || par.gridMotion~=par.noMotion )
	    gf{cur}.rx = par.rx(:,:,:,:);
	    gf{cur}.gv = zeros(par.Ngx,par.Ngy,2);
	  end
	end

  if( par.gridMotionOption==par.correctGrid && ~strcmp(par.map,'freeSurface') )
    % --- we only corrext the grid for a free surface problem ---
    return
  end

  if( par.gridMotion~=par.noMotion )


    if( par.gridMotion==par.translate )

      % NOTE: par.x holds the grid points at t=0 
      gf{cur}.x(:,:,1) = par.x(:,:,1) + par.transVect(1)*t;
      gf{cur}.x(:,:,2) = par.x(:,:,2) + par.transVect(2)*t;

      % grid velocity
      gf{cur}.gv(:,:,1) = par.transVect(1);  
      gf{cur}.gv(:,:,2) = par.transVect(2);

      % metrics are unchanged   

    elseif( par.gridMotion==par.rotate )

      % counter-clockwise rotation about the origin
      %     x(t) = R(t) x(0)
      %     R(t) = [ cos(omega*t) -sin(omega*t) ]
      %            [ sin(omega*t)  cos(omega*t) ]

      % NOTE: par.x holds the grid points at t=0 
      omega=1; % 2*pi 
      cw = cos(omega*t);
      sw = sin(omega*t);
      gf{cur}.x(:,:,1) = cw*par.x(:,:,1) - sw*par.x(:,:,2);
      gf{cur}.x(:,:,2) = sw*par.x(:,:,1) + cw*par.x(:,:,2);

      % grid velocity
      gf{cur}.gv(:,:,1) = -omega*sw*par.x(:,:,1) - omega*cw*par.x(:,:,2);
      gf{cur}.gv(:,:,2) =  omega*cw*par.x(:,:,1) - omega*sw*par.x(:,:,2);

      % metrics are rotated
      %    x.r(t) = R(t) x.r(t=0)
      % =>
      %    r.x(t) = r.x(t=0) R^{-1}(t) 
      %    [ rx ry ] = [ rx0 ry0 ] [ c s ]
      %    [ sx sy ]   [ sx0 sy0 ] [-s c ]
      gf{cur}.rx(:,:,1,1) =  cw*par.rx(:,:,1,1) - sw*par.rx(:,:,1,2);
      gf{cur}.rx(:,:,1,2) =  sw*par.rx(:,:,1,1) + cw*par.rx(:,:,1,2);

      gf{cur}.rx(:,:,2,1) = cw*par.rx(:,:,2,1) - sw*par.rx(:,:,2,2);
      gf{cur}.rx(:,:,2,2) = sw*par.rx(:,:,2,1) + cw*par.rx(:,:,2,2);

    elseif( par.gridMotion==par.deform )

      if( strcmp(par.map,'TFI' ) )

        % DEFORM A transfinite interpolation mapping 

        % --- Hard code bottom and top curves for now ---

        % --- bottom curve ---
        ampb=0.05;
        curveb =  @(r) ampb*sin(2*pi*r);
        curvebr = @(r) (2*pi*ampb)*cos(2*pi*r);
        curvebt = @(r) 0.;                       % time derivative of curveb

        % --- top curve ---
        omegat=pi;
        ampt=         0.1*cos(omegat*t);
        amptt=-omegat*0.1*sin(omegat*t);
        cty=0.5; 
        curvet =  @(r)  cty + ampt*sin(2*pi*r);
        curvetr = @(r) (2*pi*ampt)*cos(2*pi*r);
        curvett = @(r)       amptt*sin(2*pi*r); % time derivative of curvet 


        i1a = par.gid(1,1);
        i2a = par.gid(1,2);
        [I1,I2]=getIndex(par.dim);
        for i2=I2
        for i1=I1
          r1 = (i1-i1a)*par.dr(1);
          r2 = (i2-i2a)*par.dr(2);          

          gf{cur}.x(i1,i2,1)=r1; 
          gf{cur}.x(i1,i2,2)=(1-r2)*curveb(r1) + r2*curvet(r1);

          % grid velocity
          gf{cur}.gv(i1,i2,1) = 0;
          gf{cur}.gv(i1,i2,2) = (1-r2)*curvebt(r1) + r2*curvett(r1);


          xr(1,1)= 1;
          xr(2,1)= (1-r2)*curvebr(r1) + r2*curvetr(r1);
          xr(1,2)= 0; 
          xr(2,2)= -curveb(r1) + curvet(r1);
       
          rx = inv(xr);

          for m2=1:par.nd
          for m1=1:par.nd
            gf{cur}.rx(i1,i2,m1,m2) = rx(m1,m2);
          end
          end

        end
        end  


      else  
        fprintf('getGrid"ERROR: finish motion=deform for map=%s\n',par.map);
        error();  
      end

    elseif( par.gridMotion==par.freeSurfaceMotion  )

        % ===== FREE SURFACE =======
        % freeSurface : TFI with top curve defined by data points 

        i1a = par.gid(1,1);   i1b=par.gid(2,1);
        i2a = par.gid(1,2);   i2b=par.gid(2,2);
        [I1,I2]=getIndex(par.dim);
        [J1,J2]=getIndex(par.gid);
        for i1=I1
          rv(i1) = (i1-i1a)*par.dr(1);
        end 


        % --- bottom curve ---
        % FLAT FOR NOW 
        cby=-0.5;
        ampb=0.0;
        curveb =  @(r) cby + ampb*sin(2*pi*r);
        curvebr = @(r) (2*pi*ampb)*cos(2*pi*r);
        curvebt = @(r) 0;

        % --- top curve ---
        if( t<=0  )

          % --- initialize the free surface -----

          amp=par.ampfs; 
          xv = par.xa + (par.xb-par.xa)*rv;
          if( strcmp(par.icfs,'sine') )
            yv = amp*sin(par.kx*rv);
          elseif( strcmp(par.icfs,'gaussian') )
            yv = amp*exp( - (par.betag*(xv-par.x0g)).^2 );
          else
            fprintf('\n getGrid:ERROR: unknown free surface initial condition, icsf=[%s]\n',par.icsf);
          end

          % Initial time derivative: (zero for now)
          xvt = 0*xv;
          yvt = 0*yv(I1);


        else
          % ADVANCE THE FREE SURFACE 

          prev =mod(cur-2+par.numberOfGridFunctions,par.numberOfGridFunctions)+1;
          prev2=mod(cur-3+par.numberOfGridFunctions,par.numberOfGridFunctions)+1;

          if( par.gridMotionOption==par.predictGrid )

            % Solve: d(eta)/dt = (u,v)

            % AB2: 
            % -- FIX THESE FOR VARIABLE dt 
            xv = gf{prev}.eta(I1,1) + par.dt*( 1.5*gf{prev}.u(I1,i2b) -.5*gf{prev2}.u(I1,i2b) );
            yv = gf{prev}.eta(I1,2) + par.dt*( 1.5*gf{prev}.v(I1,i2b) -.5*gf{prev2}.v(I1,i2b) );

            xvt = 2*gf{prev}.u(I1,i2b) - gf{prev2}.u(I1,i2b);
            yvt = 2*gf{prev}.v(I1,i2b) - gf{prev2}.v(I1,i2b);





            if( par.plotEveryStep )
              figure(7)
              plot( xv,yv,'-x', gf{prev}.eta(I1,1),gf{prev}.eta(I1,2),'-' );
              grid on; xlabel('x');
              title(sprintf('Predict grid t=%10.2e cur=%d prev=%d prev2=%d\n',t,cur,prev,prev2));
              legend('\eta (predicted)', '\eta(t-dt)')
              pause
            end 


          elseif( par.gridMotionOption==par.correctGrid )

            % IM2 : (trap)
            xv = gf{prev}.eta(I1,1) + par.dt*( .5*gf{cur}.u(I1,i2b) + .5*gf{prev}.u(I1,i2b) );
            yv = gf{prev}.eta(I1,2) + par.dt*( .5*gf{cur}.v(I1,i2b) + .5*gf{prev}.v(I1,i2b) );

            xvt = gf{cur}.u(I1,i2b);
            yvt = gf{cur}.v(I1,i2b);

            if( par.plotEveryStep )
              figure(8)
              plot( xv,yv,'-x', gf{prev}.eta(I1,1),gf{prev}.eta(I1,2),'-' );
              grid on; xlabel('x');
              title(sprintf('Correct grid t=%10.2e cur=%d prev=%d prev2=%d\n',t,cur,prev,prev2));
              legend('\eta (corrected)', '\eta(t-dt)')
              pause
            end             

          else
            fprintf('Unknown gridMotionOption=%d\n',par.gridMotionOption);
          end

        end  % end correctGrid 

        % ----- BOUNDARY CONDITIONS ON THE FREE SURFACE POSITIONS ----

        % **** FINISH ME *****


        % PERIODICITY is done below with the periodic spline 

        % if( 1==0 && par.bc(1,1)==par.periodic )
        %   % periodic boundary conditions **FIX ME** x-periodicity needs to shift by xb-xa 
        %   xv(i1a-1) = xv(i1b-1);  yv(i1a-1) = yv(i1b-1);
        %   xv(i1b  ) = xv(i1a  );  yv(i1b  ) = yv(i1a  );
        %   xv(i1b+1) = xv(i1a+1);  yv(i1b+1) = yv(i1a+1);

        %   xvt(i1a-1) = xvt(i1b-1);  yvt(i1a-1) = yvt(i1b-1);
        %   xvt(i1b  ) = xvt(i1a  );  yvt(i1b  ) = yvt(i1a  );
        %   xvt(i1b+1) = xvt(i1a+1);  yvt(i1b+1) = yvt(i1a+1);              
        % end


        % ---- Create splines and evaluate tangential derivatives ---

        conds='second'; % second derivative zero ('natural')
        if( par.bc(1,1)==par.periodic )
          conds='periodic';
        end   
        if( par.idebug>0 )
          fprintf('*** getGrid: create splines for the free surface for cur=%d t=%9.3e ***\n',cur,t);
        end

        gf{cur}.ycs = splineInit( rv(J1),yv(J1),conds ); % spline is periodic using J1
        [yv,yvr] = splineEval( gf{cur}.ycs, rv );        % This should set the periodic images in yv 

        gf{cur}.xcs = splineInit( rv(J1),xv(J1),conds );
        [xv,xvr] = splineEval( gf{cur}.xcs, rv );                    

        gf{cur}.rvs =rv; % save the parameterization

        % Do we need to save the free surface coordinates ??
        gf{cur}.eta(I1,1) = xv(I1);
        gf{cur}.eta(I1,2) = yv(I1);


        for i2=I2
        for i1=I1
          r1 = (i1-i1a)*par.dr(1);
          r2 = (i2-i2a)*par.dr(2);          


           gf{cur}.x(i1,i2,1)=xv(i1); 
           gf{cur}.x(i1,i2,2)=(1-r2)*curveb(r1) + r2*yv(i1);

           % grid velocity
           gf{cur}.gv(i1,i2,1) = xvt(i1);
           gf{cur}.gv(i1,i2,2) = (1-r2)*curvebt(r1) + r2*yvt(i1);           

           xr(1,1)= xvr(i1);
           xr(2,1)= (1-r2)*curvebr(r1) + r2*yvr(i1);
           xr(1,2)= 0; 
           xr(2,2)= -curveb(r1) + yv(i1);
       
           rx = inv(xr);

           for m2=1:par.nd
           for m1=1:par.nd
             gf{cur}.rx(i1,i2,m1,m2) = rx(m1,m2);
           end
           end

        end
        end





    else

    	fprintf('getGrid:ERROR:finish me for gridMotion=%d\n',par.gridMotion);
    	error();
    end




	  
	end  % end if gridMotion


  return
end