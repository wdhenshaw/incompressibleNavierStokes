%
%  Save a plot an an eps or png file (based on suffix being .eps or .png)
%
function savePlotFile( plotName )
   drawnow; 
   if( endsWith(plotName,'.png') )
     print('-dpng',plotName); fprintf('Wrote file %s\n',plotName);
   else 
     print('-depsc2',plotName); fprintf('Wrote file %s\n',plotName);
   end 

end
