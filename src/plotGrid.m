 % --- plot the grid points --
function plotGrid( t, gf,cur, par )


  [I1,I2]=getIndex(par.dim);

  figure(5)
  for i2=I2
    if( i2<par.gid(1,2) || i2>par.gid(2,2) )
      plot( gf{cur}.x(I1,i2,1),gf{cur}.x(I1,i2,2),'b--','LineWidth',1); % ghost line
    else
      plot( gf{cur}.x(I1,i2,1),gf{cur}.x(I1,i2,2),'b-'); 
    end
    hold on;
  end
  for i1=I1
    if( i1<par.gid(1,1) || i1>par.gid(2,1) )
      plot( gf{cur}.x(i1,I2,1),gf{cur}.x(i1,I2,2),'k--','LineWidth',1);
    else
      plot( gf{cur}.x(i1,I2,1),gf{cur}.x(i1,I2,2),'k-');
    end
  end

  title(sprintf('Grid: t=%9.2e',t)); xlabel('x'); ylabel('y');

  % set plot bounds
  xMin = min(gf{cur}.x(:,:,1), [], "all");
  xMax = max(gf{cur}.x(:,:,1), [], "all");

  yMin = min(gf{cur}.x(:,:,2), [], "all");
  yMax = max(gf{cur}.x(:,:,2), [], "all");

  delta=.1; 
  xPlotMin = xMin- delta*(xMax-xMin);
  xPlotMax = xMax+ delta*(xMax-xMin);

  yPlotMin = yMin- delta*(yMax-yMin);
  yPlotMax = yMax+ delta*(yMax-yMin);

  xlim([xPlotMin,xPlotMax]); ylim([yPlotMin,yPlotMax]);
  grid on;
  axis equal;

  hold off;
  % pause

  return
end
