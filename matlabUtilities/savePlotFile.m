%
%  Save a matlab plot as a pdf 
%
% plotName (input) :  plot-name (without suffix '.pdf')
% type (input) : 'pdf' (default)
% res = resolution: 
%    '-r0' = screen resolution
%    '-r200' = 200 dpots per inch
%
function savePlotFile( plotName, type, res )

  if( nargin<2 ) type ='pdf'; end 
  if( nargin<3 ) res ='default'; end 

  if( strcmp(type,'pdf') )
    epsPlotName=sprintf('%s.eps',plotName); 

    %% print('-depsc2',epsPlotName); % save as a .eps file (temporary)
    % new way: Aug, 2020
    if( strcmp(res,'default') )
      saveas(gcf,epsPlotName,'epsc2'); 
    else
      print('-depsc2',epsPlotName,res);  % try this -- June 15, 2024 NOT SURE THIS WORKS 
    end
    system(sprintf('/usr/local/bin/ps2pdf -dEPSCrop %s.eps %s.pdf',plotName,plotName),'-echo');
    system(sprintf('rm %s',epsPlotName)); % now remove the .eps 
    fprintf('savePlotFile:saved plot [%s.pdf]. (resolution=%s)\n',plotName,res);

  elseif( strcmp(type,'png') )
    
    % for surf plots use png
    if( strcmp(res,'default') )
      print('-dpng',plotName);
    else
      print('-dpng',plotName,res); % try this -- June 15, 2024 : SEEMS TO WORK
    end
      fprintf('savePlotFile:saved plot [%s.png]. (resolution=%s)\n',plotName,res);

  else
    fprintf('savePlotFile:ERROR: unknown file type=[%s]\n',type);
  end 

end 
