% Function getInt: read a command line argument for an integer variable
function [ val ] = getInt( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%d',name)); 
 end
end
