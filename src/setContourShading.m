% -----------------------------------------------------------------------
% set shading for contours based on the 'shade' global variable
% -----------------------------------------------------------------------
function setContourShading(par)
  % globalDeclarations;
  if( strcmp(par.shade,'faceted') ) shading faceted; elseif( strcmp(par.shade,'flat') ) shading flat; else  shading interp; end; 
end
