% Function getReal: read a command line argument for a real variable
function [ val ] = getReal( line,name,val)
 % fprintf('getReal: val=%g line=[%s] name=[%s]\n',val,line,name);
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%e',name)); 
   % fprintf('getReal: scan for val=%g\n',val);
 end
end
