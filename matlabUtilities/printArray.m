%
%  Output an array to a file in a nice format 
%
% title : title string 
% format : format string
% file (optional input) : by default, write to command window 
% 
function printArray( u,title,format,file )

  if( nargin<4 ) file = 1; end 
  
  nd = size(u);
  
  fprintf(file,'%s\n',title);
  for( i2=nd(2):-1:1 )
    fprintf('i2=%3d:  ',i2); 
    for( i1=1:nd(1) )
      fprintf(file,format,u(i1,i2)); 
    end
    fprintf(file,'\n');
  end   


return
