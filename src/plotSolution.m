%
% plot the INS solution
%
function par = plotSolution( tn,un,vn,pn, gf,cur, par)

  cpu0 = cputime;

  % fprintf('plotSolution: tn=%g, par.savePlotThisStep=%d, par.savePlots=%d\n',tn,par.savePlotThisStep,par.savePlots); 

  Nx  = par.Ngx;
  cfl = par.cfl;
  ts  = par.ts;
  dx  = par.dx;
  dy  = par.dy;

  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y  

  if( par.plotSolutionOnGhost ) extra=par.numGhost; else extra=0; end
  [I1,I2] = getIndex( par.gid,extra );

  if 1==1 
  	figure(1)
    tl =tiledlayout('flow','TileSpacing','Compact');

    extra='';
    if( par.gravity~=0 ) extra=strcat(sprintf('g=%g, ',par.gravity )); end 
    if( par.gamma~=0 )   extra=strcat(extra,sprintf(' \\gamma=%g, ',par.gamma )); end 

    title(tl,sprintf('INS: ts=%s t=%8.2e step=%3d nu=%g, %sN=[%d,%d] cfl=%5.2f bc=%s ad=%d',par.ts,tn,par.step,par.nu,extra,par.Ngx,par.Ngy,par.cfl,par.bcLabel,par.ad),'FontSize',16);
    % title(tl,sprintf('%s: t=%5.2f, m=%d, N=[%d,%d], G=%s, cfl=%5.2f, bc=%s',gp.schemeName,t,gp.m,gp.N1,gp.N2,gp.mapName,gp.cfl,gp.bcLabel),'FontSize',16);

    nexttile;
    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),un(I1,I2) ); hold on
    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),un(I1,I2),'k-' ); 
    colormap(par.rainbowMap); colorbar; shading interp; 
    xlabel('x'); ylabel('y'); view(0,90);  % top view 
    title('u'); hold off;
    setAspectRatio();

    nexttile;
    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),vn(I1,I2) ); hold on
    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),vn(I1,I2),'k-' ); 
    colormap(par.rainbowMap); colorbar; shading interp; 
    xlabel('x'); ylabel('y'); view(0,90);  % top view 
    title('v'); hold off;
    setAspectRatio();

    nexttile;
    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),pn(I1,I2) ); hold on
    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),pn(I1,I2),'k-' ); 
    colormap(par.rainbowMap); colorbar; shading interp; 
    xlabel('x'); ylabel('y'); view(0,90);  % top view 
    title('p'); hold off;
    setAspectRatio();

    % ratio = (par.xb-par.xa)/(par.yb-par.ya); % **************************** FIX ME *********************
    
    if( par.plotAspectRatio<=0 )
      xWidth = max(max(gf{cur}.x(I1,I2,1))) - min(min(gf{cur}.x(I1,I2,1)));
      yWidth = max(max(gf{cur}.x(I1,I2,2))) - min(min(gf{cur}.x(I1,I2,2)));
      par.plotAspectRatio=xWidth/yWidth;
    end
    if( par.plotAspectRatio==1 )
      xwidth = 1200;
      ywidth = 350; 
    elseif( par.plotAspectRatio>=1 )
      xwidth = 600*par.plotAspectRatio; 
      ywidth = 300;  
    else
      xwidth = 900; % *(gf{cur}.xb-gf{cur}.xa);
      ywidth = 250*par.plotAspectRatio;  % *(par.yb-par.ya);
    end
    if( 1==0 )
      pos = get(gcf,'position');  pos(3) = xwidth; pos(4) = ywidth; set(gcf,'position',pos);  
    end

    % Hardcopy 
    timeLabel=sprintf('t%.1f',tn);
    timeLabel = replace(timeLabel, ".", "p");
    % fprintf('timeLabel=[%s]\n',timeLabel);
    myPlotName = sprintf('%s%s',par.plotName,timeLabel);    
    if( par.savePlotThisStep && par.savePlots )
      savePlotFile( sprintf('%s/%sVelocityAndPressure',par.figDir,myPlotName),'pdf' );
    end 


    % --- plot vorticity and streamlines  ---
    if( par.plotVorticity && strcmp(par.ms,'none') )
   	  figure(2)
      tl =tiledlayout('flow','TileSpacing','Compact');

      title(tl,sprintf('ts=%s t=%8.2e step=%3d nu=%g N=[%d,%d] cfl=%5.2f bc=%s ad=%d',par.ts,tn,par.step,par.nu,par.Ngx,par.Ngy,par.cfl,par.bcLabel,par.ad),'FontSize',16);    

      vor = zeros(par.Ngx,par.Ngy);
      [J1,J2] = getIndex( par.gid );
      vor(J1,J2) = Dzy(un,J1,J2) - Dzx(vn,J1,J2);
      nexttile;
      surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),vor(I1,I2) ); hold on
      contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),vor(I1,I2),'k-' ); 
      colormap(par.rainbowMap); colorbar; shading interp; 
      xlabel('x'); ylabel('y'); view(0,90);  % top view 
      title('vorticity'); hold off;
      setAspectRatio();    

      % -- plot streamlines ---
      nexttile;
      par = plotStreamLines( un,vn,par );
      setAspectRatio();

      if( ratio==1 )
        xwidth = 800*ratio; 
        ywidth = 400; 
      elseif( ratio>=1 )
        xwidth = 600*ratio; 
        ywidth = 300; 
      else
        xwidth = 600; 
        ywidth = 300*ratio; 
      end

      pos = get(gcf,'position');  pos(3) = xwidth; pos(4) = ywidth; set(gcf,'position',pos);    

      % Horizontal placement
      % xwidth =  350; % *(par.bx-par.ax);
      % ywidth = 1100; % *(par.by-par.ay);
      % pos = get(gcf,'position');  pos(3) = xwidth; pos(4) = ywidth; set(gcf,'position',pos);

      % fprintf('plotSolution: par.savePlotThisStep=%d, par.savePlots=%d\n',par.savePlotThisStep,par.savePlots);
      % Hardcopy 
      if( par.savePlotThisStep && par.savePlots )
         savePlotFile( sprintf('%s/%sVorticityAndStreamLines',par.figDir,myPlotName),'pdf' );
      end 
    end

    % ---- ERRORS ----
    if( par.computeErrors && par.plotErrors )
	    [maxErr,perr,uerr,verr,div] = getErrors( tn,un,vn,pn, gf,cur, par );

      par.maxErr(1) = max(abs(perr(I1,I2)),[],'all');
      par.maxErr(2) = max(abs(uerr(I1,I2)),[],'all');
      par.maxErr(3) = max(abs(verr(I1,I2)),[],'all');
      par.maxErr(4) = max(abs(div (I1,I2)),[],'all');

	 	  figure(3)
	    tl =tiledlayout('flow','TileSpacing','Compact');

	    title(tl,sprintf('INS: ts=%s t=%5.2f step=%3d nu=%g N=[%d,%d] cfl=%5.2f bc=%s ad=%d',par.ts,tn,par.step,par.nu,par.Ngx,par.Ngy,par.cfl,par.bcLabel,par.ad),'FontSize',16);

	    nexttile;
	    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),uerr(I1,I2) ); hold on
	    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),uerr(I1,I2),'k-' ); 
	    colormap(par.rainbowMap); colorbar; shading interp; 
	    xlabel('x'); ylabel('y'); view(0,90);  % top view 
	    title('u-err'); hold off;
	    setAspectRatio();

	    nexttile;
	    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),verr(I1,I2) ); hold on
	    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),verr(I1,I2),'k-' ); 
	    colormap(par.rainbowMap); colorbar; shading interp; 
	    xlabel('x'); ylabel('y'); view(0,90);  % top view 
	    title('v-err'); hold off;
	    setAspectRatio();

	    nexttile;
	    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),perr(I1,I2) ); hold on
	    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),perr(I1,I2),'k-' ); 
	    colormap(par.rainbowMap); colorbar; shading interp; 
	    xlabel('x'); ylabel('y'); view(0,90);  % top view 
	    title('p-err'); hold off;
	    setAspectRatio();

	    nexttile;
	    surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),div(I1,I2) ); hold on
	    contour3( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),div(I1,I2),'k-' ); 
	    colormap(par.rainbowMap); colorbar; shading interp; 
	    xlabel('x'); ylabel('y'); view(0,90);  % top view 
	    title('div'); hold off;
	    setAspectRatio();
    end
  else 
  	% OLD WAY 

		colormap(par.rainbowMap);
		figure(1); 
		subplot(2,2,1); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),un(I1,I2) ); setContourShading();
		title(sprintf('u: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 
		subplot(2,2,2); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),vn(I1,I2) ); setContourShading();
		title(sprintf('v: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 
		subplot(2,2,3); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),pn(I1,I2) ); setContourShading();
		title(sprintf('p: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 

		% divergence
		div(I1,I2) = Dzx(un,I1,I2)+Dzy(vn,I1,I2);
		divMax=max(max(abs(div(I1,I2)))); 
		subplot(2,2,4); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),div(I1,I2) ); setContourShading();
		title(sprintf('div: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 

		[maxErr,perr,uerr,verr,div] = getErrors( tn,un,vn,pn, gf,cur, par );

		figure(2); 
		colormap(par.rainbowMap);
		subplot(2,2,1); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),uerr(I1,I2) ); setContourShading();
		title(sprintf('u-err: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 
		subplot(2,2,2); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),verr(I1,I2) ); setContourShading();
		title(sprintf('v-err: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 
		subplot(2,2,3); surf( gf{cur}.x(I1,I2,1),gf{cur}.x(I1,I2,2),perr(I1,I2) ); setContourShading();
		title(sprintf('p-err: t=%9.3e, ts=%s (Nx=%d) cfl=%4.2f',tn,ts,Nx,cfl)); xlabel('x'); ylabel('y'); drawnow; commandwindow; 

  end

  if( par.plotGrid )
    plotGrid(  tn, gf,cur, par );
  end


  par.cpuPlot = par.cpuPlot + cputime - cpu0;
return 
end
