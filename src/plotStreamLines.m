function par = plotStreamLines( u,v,par )

  cpu0 = cputime;

  plotStreamLineArrows = 1;

  % numArrowsPerSegment = 3;

  % streamline stopping tol: 
  par.streamLineStoppingTolerance=1e-3;   

  xa = par.xa;
  xb = par.xb;
  ya = par.ya;
  yb = par.yb;

  xba=xb-xa;
  yba=yb-ya;  


  % streamline grid is 40x40 by default
  par.nxg=40; % ceil(par.Nx/2);
  par.nyg=40; % ceil(par.Ny/2);
  if( xba>yba )
    par.nxg=ceil( par.nxg*xba/yba);
  else
    par.nyg=ceil( par.nyg*yba/xba);
  end

  par.maskForStreamLines = zeros(par.nxg,par.nyg);

  par.uMin = min( min(min(u)), min(min(v)) );
  par.uMax = max( max(max(u)), max(max(v)) );

  % uScale = par.uMax - par.uMin;

  nxg = par.nxg;
  nyg = par.nyg;  

  % plot a fake line off screen
  h = plot( [xa-1,xa-1], [ya,yb], 'k-' );

  size = max(xb-xa,yb-ya); % over-all domain size : determines size of arrows 
  scale=0.5; % scale size of arrows

  par.maxSpeed=0;
  par.numArrows=0; 

  hold on;
  for( i=1:par.nxg )
    for( j=1:par.nyg )
      if( par.maskForStreamLines(i,j)==0 )  % no streamline has passed through this point
        xtp=xa+xba*(i-.5)/nxg;  % starting point for streamline
        ytp=ya+yba*(j-.5)/nyg;
        % first integrate backwards in time from this spot
        par.cflStreamLine = -.5; % plot in backward direction
        [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  

        % yp(end)=NaN;
        % fprintf('xp=%d, yp=%d, cp=%d\n',length(xp),length(yp),length(cp));
        % pause
        if( length(xp)>3 )
          patch([xp,xp(end)],[yp,NaN],[cp,cp(end)],'EdgeColor','interp','LineWidth',2,'LineJoin','round');
   
          % if( plotStreamLineArrows && length(xp)>10 )
          %   % *** THIS IS EXPENSIVE*** ARROWS DONE BELOW NOW
          %   % Note: flip (reverse) order
          %   plotArrows(h,flip(xp), flip(yp), 'number', numArrowsPerSegment, 'color', 'k', 'LineWidth', 1,'scale', scale, 'ratio', 'equal', 'size',size ); 
          % end
        end
        % plot( xp,yp,'k-'); 
        hold on

        % ---  Now plot the streamline in the forward direction ---
        par.cflStreamLine = .5; % plot in forward direction
        [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par ); 
        % fprintf('xp=%d, yp=%d, cp=%d\n',length(xp),length(yp),length(cp));
        if( length(xp)>3 )
          patch([xp,xp(end)],[yp,NaN],[cp,cp(end)],'EdgeColor','interp','LineWidth',2,'LineJoin','round');

          % if( plotStreamLineArrows && length(xp)>10 )
          %   % *** THIS IS EXPENSIVE*** ARROWS DONE BELOW NOW
          %   % h = plot( [xa-1,xb-1], [ya,yb], 'k-' ); % plot a fake line off screen
          %   plotArrows(h, xp, yp, 'number', numArrowsPerSegment, 'color', 'k', 'LineWidth', 1,'scale', scale, 'ratio', 'equal', 'size',size );
          % end         
        end       
        % plot( xp,yp,'k-');
      end
    end
  end

  if( par.numArrows>0 )
    % plot(par.arrows(:,1),par.arrows(:,2),'kx');

    if( 1==1 )
      %-- the template triangle (it points "east", centered in (0,0)):
      clr = [.4,.4,.4]; 'k';
      % xyScale = .01*max(par.xb-par.xa,par.yb-par.ya);
      xScale=.01*(par.xb-par.xa);
      yScale=.01*(par.yb-par.ya);
      ArWidth = .5;  % aspect ratio, smaller -> thinner triangle
      xt = [1 -sin(pi/6)  -sin(pi/6)];
      yt = ArWidth*[0  cos(pi/6)  -cos(pi/6)];

      for( i=1:par.numArrows )
        x     = par.arrows(i,1);
        y     = par.arrows(i,2);
        theta = par.arrows(i,3);
        % for i=1:3
        %   xd(i) = x + cos(theta)*xyScale*xt(i) - sin(theta)*xyScale*yt(i);
        %   yd(i) = y + sin(theta)*xyScale*xt(i) + cos(theta)*xyScale*yt(i);
        % end  
        for i=1:3
          xd(i) = x + xScale*(cos(theta)*xt(i) - sin(theta)*yt(i));
          yd(i) = y + yScale*(sin(theta)*xt(i) + cos(theta)*yt(i));
        end                
        patch(xd,yd,clr,'EdgeColor',clr,'facecolor','none', 'LineWidth',2);
      end 

    % %-- the arrow is made of a transformed "template triangle".
    % %-- it will be created, rotated, moved, resized and shifted.
    % %-- the template triangle (it points "east", centered in (0,0)):
    % xt = [1 -sin(pi/6)  -sin(pi/6)];
    % yt = ArWidth*[0  cos(pi/6)  -cos(pi/6)];
    % clr='k';

    % for( i=1:par.numArrows )
    %   theta = par.arrows(i,3);
    %       %-- rotate it by the angle determined above:
    %   xd = []; yd = [];
    %   for i=1:3
    %     xd(i) = cos(p)*xt(i) - sin(p)*yt(i);
    %     yd(i) = sin(p)*xt(i) + cos(p)*yt(i);
    %   end


    %   %-- move the triangle so that its "head" lays in (0,0):
    %   xd = xd - cos(theta);
    %   yd = yd - sin(theta);
    %   %-- stretch/deform the triangle to look good on the current axes:
    %   xd = xd*Xextend*ArSize/10000;
    %   yd = yd*Yextend*ArSize/10000;
    %   %-- move the triangle to the location where it's needed
    %   xd = xd + x2;
    %   yd = yd + y2;
    %   %-- draw the actual triangle
    %   patch(xd,yd,clr,'EdgeColor',clr);

    %   %wdh handle(Loop) = patch(xd,yd,clr,'EdgeColor',clr);
    %   %wdh if nonsolid, set(handle(Loop),'facecolor','none'); end
    end

    % u = cos(par.arrows(:,3));
    % v = sin(par.arrows(:,3));
    % quiver(par.arrows(:,1),par.arrows(:,2), u, v, 'r', 'LineWidth', 2);
  end

  hold off;
  colormap(par.rainbowMap);  colorbar; % shading interp; 
  xlabel('x'); ylabel('y');

  % clim([0,uScale]);
  if( par.maxSpeed>0 )
    clim([0,par.maxSpeed]);
  end
  xlim([xa,xb]);
  ylim([ya,yb]);
  title('streamlines');

  cpuStreamLines = cputime-cpu0;

  if( 1==0 || par.idebug> 1) fprintf('Plot streamlines: cpu=%9.2e(s) plotStreamLineArrows=%d\n',cpuStreamLines,plotStreamLineArrows); end

  return
end