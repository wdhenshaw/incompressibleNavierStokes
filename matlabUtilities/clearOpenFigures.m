  %
  % Clear a list of open figures:
  %
  function clearOpenFigures( figList ) 
    if( nargin<1 )
    	figList=1:100; % by default clear figs 1..100
    end
    % Clear all open figures in figList 
    for( ifig = figList )
      % fprintf('Clear fig=%d ?\n',ifig);  
      try
        clf(ifig);
      catch
            % Nothing to do
      end
    end
  % --- OLD: This also clears GUI's !!
  % FigList = findall(groot, 'Type', 'figure');
  % for iFig = 1:numel(FigList)
  %     try
  %         clf(FigList(iFig));
  %     catch
  %         % Nothing to do
  %     end
  % end
  return
end
