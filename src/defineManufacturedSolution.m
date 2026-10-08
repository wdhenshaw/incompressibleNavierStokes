%
%  Define the manufactured solution for INS
%
function par = defineManufacturedSolution( par )

  ms = par.ms;
  kx      = par.kx;
  ky      = par.ky;
  kt      = par.kt;
  nu      = par.nu;
  tzScale = par.tzScale;

  if( ~strcmp(ms,'none') )
    if( strcmp(ms,'trig') )
      % --- trigonometric exact solution ---
      %     This function is divergence free.
      par.computeErrors=1;
      if( tzScale==0 )
        ampu = ky/sqrt(kx^2+ky^2); 
        ampv = kx/sqrt(kx^2+ky^2); 
      else
        ampu = ky/(kx^2+ky^2); 
        ampv = kx/(kx^2+ky^2); 
      end;
      ampp=1.; 
      par.ue = @(x,y,t) ampu*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      par.ve = @(x,y,t) ampv*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      par.pe = @(x,y,t) ampp*cos(kx*x).*sin(ky*y)*cos(kt*t); 

      % derivatives
      par.uex  = @(x,y,t) ampu*(-kx)*sin(kx*x).*cos(ky*y)*cos(kt*t); 
      par.vex  = @(x,y,t) ampv*( kx)*cos(kx*x).*sin(ky*y)*cos(kt*t); 
      par.pex  = @(x,y,t) ampp*(-kx)*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      
      par.uey  = @(x,y,t) ampu*(-ky)*cos(kx*x).*sin(ky*y)*cos(kt*t); 
      par.vey  = @(x,y,t) ampv*( ky)*sin(kx*x).*cos(ky*y)*cos(kt*t); 
      par.pey  = @(x,y,t) ampp*( ky)*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      
      par.uet = @(x,y,t) ampu*(-kt)*cos(kx*x).*cos(ky*y)*sin(kt*t); 
      par.vet = @(x,y,t) ampv*(-kt)*sin(kx*x).*sin(ky*y)*sin(kt*t); 
      par.pet = @(x,y,t) ampp*(-kt)*cos(kx*x).*sin(ky*y)*sin(kt*t); 

      par.uexx = @(x,y,t) ampu*(-kx^2)*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      par.vexx = @(x,y,t) ampv*(-kx^2)*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      par.pexx = @(x,y,t) ampp*(-kx^2)*cos(kx*x).*sin(ky*y)*cos(kt*t); 
      
      par.uexy = @(x,y,t) ampu*( kx*ky)*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      par.vexy = @(x,y,t) ampv*( kx*ky)*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      par.pexy = @(x,y,t) ampp*(-kx*ky)*sin(kx*x).*cos(ky*y)*cos(kt*t); 

      par.ueyy = @(x,y,t) ampu*(-ky^2)*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      par.veyy = @(x,y,t) ampv*(-ky^2)*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      par.peyy = @(x,y,t) ampp*(-ky^2)*cos(kx*x).*sin(ky*y)*cos(kt*t); 

      par.uett = @(x,y,t) ampu*(-kt^2)*cos(kx*x).*cos(ky*y)*cos(kt*t); 
      par.vett = @(x,y,t) ampv*(-kt^2)*sin(kx*x).*sin(ky*y)*cos(kt*t); 
      par.pett = @(x,y,t) ampp*(-kt^2)*cos(kx*x).*sin(ky*y)*cos(kt*t); 

    elseif( strcmp(ms,'poly') )

		if (~isfield(par,'degreex'))
			par.degreex=2; 
		end
		if (~isfield(par,'degreet'))
			par.degreet=2;    % degree of poly MS % WARN: Maybe case specific?
		end
  		
		degreex = par.degreex;
  		degreet = par.degreet;

      % --- polynomial manufactured solution ---
      %    Note: choose the solution to be divergence free  
      par.computeErrors=1;        
      cu0=1.; cu1=0; cu2=0; 
      if( degreet>0 ) cu1=.5; end; 
      if( degreet>1 ) cu2=.25; end; 
      timePoly = @(t) (cu0 + t*(cu1 + t*cu2)); 
      timePolyt = @(t) (cu1 + 2.*cu2*t);
      timePolytt = @(t) 2.*cu2;

      cu0=1; cuxx=0; cuyy=0; cux=0; cuy=0; cuxy=0; 
      if( degreex>0 ) cux=1.; cuy=.5; end;
      if( degreex>1 ) cuxx=1; cuyy=1.; cuxy=2.; end; 
      par.ue = @(x,y,t) (cu0 + cux*x + cuy*y + cuxx*x.*x + cuxy*x.*y + cuyy*y.*y)*timePoly(t);

      par.uex  = @(x,y,t) ( cux + cuxy*y + 2.*cuxx*x )*timePoly(t);
      par.vey  = @(x,y,t) ( cvy + cvxy*x + 2.*cvyy*y )*timePoly(t);

      par.uey  = @(x,y,t) ( cuy + cuxy*x + 2.*cuyy*y )*timePoly(t);
      % *wdh* mistake found Aug 5, 2021
      % uexx = @(x,y,t) (       cuxy               )*timePoly(t);
      % uexy = @(x,y,t) ( 2.*cuxx                  )*timePoly(t);    
      par.uexx = @(x,y,t) ( 2.*cuxx                  )*timePoly(t);
      par.uexy = @(x,y,t) (       cuxy               )*timePoly(t);
      par.ueyy = @(x,y,t) ( 2.*cuyy                  )*timePoly(t);
      par.uet  = @(x,y,t) (cu0 + cux*x + cuy*y + cuxx*x.*x + cuxy*x.*y + cuyy*y.*y)*timePolyt(t);
      par.uett = @(x,y,t) (cu0 + cux*x + cuy*y + cuxx*x.*x + cuxy*x.*y + cuyy*y.*y)*timePolytt(t);

      cv0=.5; cvxx=0; cvyy=0; cvx=0; cvy=0; cvxy=0.; 
      if( degreex>0 ) cvx=.75; cvy=-cux; end;
      if( degreex>1 ) cvxx=1; cvyy=-.5*cuxy; cvxy=-2.*cuxx; end; 
      par.ve = @(x,y,t) (cv0 + cvx*x + cvy*y + cvxx*x.*x + cvxy*x.*y + cvyy*y.*y)*timePoly(t);
      par.vex  = @(x,y,t) ( cvx + cvxy*y + 2.*cvxx*x )*timePoly(t);
      par.vey  = @(x,y,t) ( cvy + cvxy*x + 2.*cvyy*y )*timePoly(t);
      par.vexx = @(x,y,t) ( 2.*cvxx                  )*timePoly(t);
      par.vexy = @(x,y,t) (    cvxy                  )*timePoly(t);
      par.veyy = @(x,y,t) ( 2.*cvyy                  )*timePoly(t);
      par.vet  = @(x,y,t) (cv0 + cvx*x + cvy*y + cvxx*x.*x + cvxy*x.*y + cvyy*y.*y)*timePolyt(t);
      par.vett = @(x,y,t) (cv0 + cvx*x + cvy*y + cvxx*x.*x + cvxy*x.*y + cvyy*y.*y)*timePolytt(t);


      cp0=1; cpx=0.; cpy=.0; cpxx=0; cpxy=0.; cpyy=0; 
      if( degreex>0 ) cpx=.50; cpy=.25; end;
      if( degreex>1 ) cpxx=1; cpyy=1.5; cpxy=-.5; end; 
      par.pe   = @(x,y,t) cp0 + cpx*x + cpy*y + cpxx*x.*x + cpxy*x.*y + cpyy*y.*y;
      par.pex  = @(x,y,t) cpx + 2.*cpxx*x + cpxy*y;
      par.pey  = @(x,y,t) cpy + 2.*cpyy*y + cpxy*x;
      par.pexx = @(x,y,t) 2.*cpxx;
      par.peyy = @(x,y,t) 2.*cpyy;

    else
      fprintf('defineManufacturedSolution: ERROR: unknown manufactured solution: ms=[%s]\n',ms);
      pause; pause; 
    end;

    % NOTE: DO NOT ADJUST THE FORCING FOR THE GRID VELOCITY 
    %       SINCE THIS IS THE FORCING IN THE NON-MOVING FRAME 
    % par.gvu = @(x,y,t) 0.;
    % par.gvv = @(x,y,t) 0.;    
    % par.ufe = @(x,y,t) par.uet(x,y,t) + (par.ue(x,y,t)-par.gvu(x,y,t)).*par.uex(x,y,t) + (par.ve(x,y,t)-par.gvv(x,y,t)).*par.uey(x,y,t) + par.pex(x,y,t) - nu*( par.uexx(x,y,t) + par.ueyy(x,y,t) );
    % par.vfe = @(x,y,t) par.vet(x,y,t) + (par.ue(x,y,t)-par.gvu(x,y,t)).*par.vex(x,y,t) + (par.ve(x,y,t)-par.gvv(x,y,t)).*par.vey(x,y,t) + par.pey(x,y,t) - nu*( par.vexx(x,y,t) + par.veyy(x,y,t) );
    % par.pfe = @(x,y,t) par.pexx(x,y,t) + par.peyy(x,y,t) + ( par.uex(x,y,t).*par.uex(x,y,t) + 2.*par.uey(x,y,t).*par.vex(x,y,t) + par.vey(x,y,t).*par.vey(x,y,t) ) ; 

    % if( par.gridMotion==par.noMotion )
    %   par.gvu = @(x,y,t) 0.;
    %   par.gvv = @(x,y,t) 0.;
    % elseif( par.gridMotion==par.translate )
    %   par.gvu = @(x,y,t) 0.; % par.transVect(1);
    %   par.gvv = @(x,y,t) 0.; % par.transVect(2);
    % else
    %   fprintf('defineManufacturedSolution: ERROR: finish grid velocity\n');
    %   error();
    % end

    % --- MS FORCING ---
    par.ufe = @(x,y,t) par.uet(x,y,t) + par.ue(x,y,t).*par.uex(x,y,t) + par.ve(x,y,t).*par.uey(x,y,t) + par.pex(x,y,t) - nu*( par.uexx(x,y,t) + par.ueyy(x,y,t) );
    par.vfe = @(x,y,t) par.vet(x,y,t) + par.ue(x,y,t).*par.vex(x,y,t) + par.ve(x,y,t).*par.vey(x,y,t) + par.pey(x,y,t) - nu*( par.vexx(x,y,t) + par.veyy(x,y,t) );
    par.pfe = @(x,y,t) par.pexx(x,y,t) + par.peyy(x,y,t) + ( par.uex(x,y,t).*par.uex(x,y,t) + 2.*par.uey(x,y,t).*par.vex(x,y,t) + par.vey(x,y,t).*par.vey(x,y,t) ) ; 

    % par.pe2  = @(x,y,t) pe(x,y,t); % to avoid matlab function pe

  end

end
