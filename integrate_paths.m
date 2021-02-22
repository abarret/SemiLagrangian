function [xx_new,yy_new] = integrate_paths(xx, yy, gcw, u, v, dt)
uu = u(xx(gcw+1:end-gcw, gcw+1:end-gcw),yy(gcw+1:end-gcw, gcw+1:end-gcw));
vv = v(xx(gcw+1:end-gcw, gcw+1:end-gcw),yy(gcw+1:end-gcw, gcw+1:end-gcw));
xxhat = xx(gcw+1:end-gcw, gcw+1:end-gcw) - 0.5*dt * uu;
yyhat = yy(gcw+1:end-gcw, gcw+1:end-gcw) - 0.5*dt * vv;
uu = u(xxhat, yyhat);
vv = v(xxhat, yyhat);
xx_new = xx(gcw+1:end-gcw, gcw+1:end-gcw) - dt * uu;
yy_new = yy(gcw+1:end-gcw, gcw+1:end-gcw) - dt * vv;
end