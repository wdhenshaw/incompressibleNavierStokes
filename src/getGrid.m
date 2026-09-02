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


      elseif( strcmp(par.map,'freeSurface') )

        % freeSurface : TFI with top curve defined by data points 

        i1a = par.gid(1,1);
        i2a = par.gid(1,2);
        [I1,I2]=getIndex(par.dim);
        for i1=I1
          rv(i1) = (i1-i1a)*par.dr(1);
        end 


        % --- bottom curve ---
        cby=-0.5;
        ampb=0.0;
        curveb =  @(r) cby + ampb*sin(2*pi*r);
        curvebr = @(r) (2*pi*ampb)*cos(2*pi*r);
        curvebt = @(r) 0;

        % --- top curve ---
        omegat=pi;
        ampt=         0.1*cos(omegat*t);
        amptt=-omegat*0.1*sin(omegat*t);

        yv  =  ampt*sin(3*pi*rv);
        yvt = amptt*sin(3*pi*rv); % time derivative of top curve 

        bcLeft=0; gLeft=0; bcRight=0; gRight=0; % natural BCs
        coeff = splineCoeff( rv,yv, bcLeft,gLeft, bcRight,gRight );

        % eval spline and derivative at points rv
        [yv0,yvr] = splineEval( rv,yv,coeff, rv );    

        for i2=I2
        for i1=I1
          r1 = (i1-i1a)*par.dr(1);
          r2 = (i2-i2a)*par.dr(2);          


           gf{cur}.x(i1,i2,1)=r1; 
           gf{cur}.x(i1,i2,2)=(1-r2)*curveb(r1) + r2*yv(i1);

           % grid velocity
           gf{cur}.gv(i1,i2,1) = 0;
           gf{cur}.gv(i1,i2,2) = (1-r2)*curvebt(r1) + r2*yvt(i1);           

           xr(1,1)= 1;
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
        fprintf('getGrid"ERROR: finish motion=deform for map=%s\n',par.map);
        error();  
      end

    else
    	fprintf('getGrid:ERROR:finish me for gridMotion=%d\n',par.gridMotion);
    	error();
    end


	  
	end 


  return
end