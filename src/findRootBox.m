%
% Locate the root with the largest imaginary part of an analytic function F inside a box in the
% complex plane using the argument principle (a 2D analogue of bisection).
%
%   searchBox = [ReLo, ReHi, ImLo, ImHi], F must have no poles in the box.
%   Boxes containing roots are split into four until they are smaller than boxTol.
%   If there is no root in searchBox then omegaStart=[] and info.numRoots=0.
%
function [omegaStart, rootBox, info] = findRootBox( F, searchBox, boxTol, idebug )

  omegaStart = [];
  rootBox = [];
  info.numLevels = 0;

  [info.numRoots, info.numFEvals] = countRoots( F, searchBox );
  if( info.numRoots==0 )
    return
  end

  % Each row of boxes is [ReLo, ReHi, ImLo, ImHi, numRoots]
  boxes = [searchBox(1:4), info.numRoots];
  maxLevels = 60;
  while( max( max(boxes(:,2)-boxes(:,1)), max(boxes(:,4)-boxes(:,3)) ) > boxTol )
    info.numLevels = info.numLevels + 1;
    if( info.numLevels > maxLevels )
      error('findRootBox: root isolation did not finish in %d levels.',maxLevels);
    end

    newBoxes = zeros(0,5);
    for j=1:size(boxes,1)
      reMid = 0.5*( boxes(j,1) + boxes(j,2) );
      imMid = 0.5*( boxes(j,3) + boxes(j,4) );
      children = [ boxes(j,1), reMid,      boxes(j,3), imMid,      0;
                   reMid,      boxes(j,2), boxes(j,3), imMid,      0;
                   boxes(j,1), reMid,      imMid,      boxes(j,4), 0;
                   reMid,      boxes(j,2), imMid,      boxes(j,4), 0 ];
      for c=1:4
        [children(c,5), numEvals] = countRoots( F, children(c,1:4) );
        info.numFEvals = info.numFEvals + numEvals;
      end
      if( sum(children(:,5)) ~= boxes(j,5) )
        error('findRootBox: inconsistent root count while splitting a box (root on a box edge?).');
      end
      newBoxes = [ newBoxes; children(children(:,5)>0,:) ];
    end
    boxes = newBoxes;

    if( idebug>1 )
      fprintf('  findRootBox level %2d: boxes with roots=%d, roots=%d, box size=%9.2e x %9.2e, evaluations of F=%d\n',...
              info.numLevels,size(boxes,1),sum(boxes(:,5)),boxes(1,2)-boxes(1,1),boxes(1,4)-boxes(1,3),info.numFEvals);
    end
  end

  % Keep the root with the largest imaginary part
  [~, j] = max( boxes(:,4) );
  rootBox = boxes(j,1:4);
  omegaStart = 0.5*( rootBox(1) + rootBox(2) ) + 1i*0.5*( rootBox(3) + rootBox(4) );

end


% Number of roots of F inside the box: the winding number of F around the boundary. The boundary
% is traversed counter-clockwise with parameter t in [0,4), one unit per edge. Points are added
% where arg(F) changes by more than pi/4 so that no winding around zero is missed.
function [numRoots, numEvals] = countRoots( F, box )

  corners = [ box(1)+1i*box(3), box(2)+1i*box(3), box(2)+1i*box(4), box(1)+1i*box(4), box(1)+1i*box(3) ];
  boundaryPoint = @(t) corners(floor(t)+1) + ( t - floor(t) ).*( corners(floor(t)+2) - corners(floor(t)+1) );

  numPerEdge = 16;
  t = (0:4*numPerEdge-1)/numPerEdge;
  Fv = arrayfun( F, boundaryPoint(t) );
  numEvals = length(t);

  maxRefinements = 20;
  for r=1:maxRefinements
    dAngle = angle( Fv([2:end,1])./Fv );
    big = find( abs(dAngle) > pi/4 );
    if( isempty(big) )
      break;
    end
    tNext = [ t(2:end), 4 ];
    tMid = 0.5*( t(big) + tNext(big) );
    FMid = arrayfun( F, boundaryPoint(tMid) );
    numEvals = numEvals + length(tMid);
    [t, order] = sort( [t, tMid] );
    Fv = [Fv, FMid];
    Fv = Fv(order);
  end

  dAngle = angle( Fv([2:end,1])./Fv );
  if( any( abs(dAngle) > pi/2 ) )
    error('findRootBox: a root lies on or very near the edge of the box [%g,%g]x[%g,%g].',box(1),box(2),box(3),box(4));
  end

  windingNumber = sum(dAngle)/(2*pi);
  numRoots = round( windingNumber );
  if( abs(windingNumber - numRoots) > 0.1 )
    error('findRootBox: inaccurate winding number %g.',windingNumber);
  end

end
