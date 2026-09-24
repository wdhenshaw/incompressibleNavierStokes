% Function getString: read a command line argument for a string variable
function [ val ] = getString( line,name,val)

 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%s',name)); 
 end
end
