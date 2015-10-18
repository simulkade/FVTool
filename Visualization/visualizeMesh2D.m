function visualizeMesh2D(m)
  [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
  [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
  plot(X, Y, 'or', ...
       Xf, Yf, '-b', Xf', Yf', '-b');
end