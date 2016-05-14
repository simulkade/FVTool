% How to use cellLocations and faceLocations functions to define cell and
% face values, variable in space
% Cell Variable:
m=createMesh2D(30, 20, 5.0, 3.0);
[X, Y, Z]=cellLocations(m);
c=createCellVariable(m, sin(X.value).*cos(Y.value));
figure(1);
visualizeCells(c);

% face variable:
[Xf, Yf, Zf]=faceLocations(m);
v=createFaceVariable(m, 0);
v.xvalue=sin(Xf.xvalue).*Xf.yvalue;
v.yvalue=cos(Yf.xvalue).*sin(Yf.yvalue);
figure(2); 
subplot(2,1,1);
pcolor(v.xvalue')
subplot(2,1,2);
pcolor(v.yvalue');