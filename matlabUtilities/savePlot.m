%
%  Save a matlab plot as a pdf 
%
% plotName (input) :  plot-name (without suffix '.pdf')
% type (input) : 'pdf' (default)
%
function  savePlot( plotName, type )

  if( nargin<2 ) type ='pdf'; end 

  if( strcmp(type,'pdf') )
    epsPlotName=sprintf('%s.eps',plotName); 
    print('-depsc2',epsPlotName); % save as a .eps file (temporary)
    system(sprintf('/usr/local/bin/ps2pdf -dEPSCrop %s.eps %s.pdf',plotName,plotName),'-echo');
    system(sprintf('rm %s',epsPlotName)); % now remove the .eps 
    fprintf('Saved plot [%s.pdf].\n',plotName);
  else
    fprintf('savePlot:ERROR: unknown file type=[%s]\n',type);
  end 

end 
