function visualizeMesh2D(m)
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
  [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
  [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
  plot(X, Y, 'or', ...
       Xf, Yf, '-b', Xf', Yf', '-b');
end
