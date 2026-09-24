%
% Initialize the polynomial TZ function 
%
% Input:
%  par.degreex : degree in space of the poly
%  par.degreet : degree in time of the poly
%  par.nc      : number of components
%
function par = initPolyTZ( par )

  degreex    = par.degreex;
  degreet    = par.degreet;
  nc         = par.nc; 
    
  % define the coefficients in the poly

  if( isfield(par,'xpoly') )
    % User has supplied the polynomial coefficients in x 
  else
    % generate a polynomial of the desired degree
    par.xpoly = zeros(degreex+1,degreex+1,nc);
    % general polynomial 
    for( ix=0:degreex )
      for( iy=0:degreex )
        if( ix+iy <= degreex )
          for( ic=1:nc )
            par.xpoly(ix+1,iy+1,ic) = 1./( 1+ ix*(1.5-ic) + 2*iy*(ic+1) + 1.5*ic  ); 
            % par.xpoly(ix+1,iy+1,ic,d) = 1./( 1+ ix*(1.5) + 2*iy*(1) + d ); 
            % par.xpoly(ix+1,iy+1,ic,d) = 1./( 1+ ix*(1.5) + 2*iy*(1) );   
          end 
        end
      end
    end     
  end

  if( isfield(par,'tpoly') )
    % User has supplied the polynomial coefficients in t 
  else
    par.tpoly = zeros(degreet+1,nc);
    for( it=0:degreet )
     for( ic=1:nc )
       par.tpoly(it+1,ic) = 1./( 2 + it*3 + 1.5*ic );
     end
    end   
  end 

  fprintf('+++ initPolyTZ : DEFINE POLY COEFF degreex=%d degreet=%d nc=%d ++++\n',...
        degreex,degreet,nc );

 if( 0==1 )
   printArray(par.xpoly,'xpoly(degreex+1,degreex+1,nc)','%5.2f ');
   pause
 end
 


return
end
