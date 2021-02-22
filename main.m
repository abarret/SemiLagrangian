clear all;
save_figs = 0;
draw_num = 12;
base_path = 'Graphs/semi_lag_adv_';
CFL_MAX = 0.5;
Nx = 64;
Ny = 64;

Lx = 1.0;
Ly = 1.0;

poly_order = 2;
gcw = get_width_of_spline(poly_order) + 2;

x_side  = linspace(0, Lx, Nx);
dx = x_side(2) - x_side(1);
x_side = linspace(0 - gcw*dx, Lx + gcw*dx, Nx + 2*gcw);
y_side = linspace(0, Ly, Ny);
dy = y_side(2) - y_side(1);
y_side = linspace(0 - gcw*dy, Ly + gcw*dy, Ny + 2*gcw);

x_cent = 0.5 * (x_side(2:end) + x_side(1:end-1));
y_cent = 0.5 * (y_side(2:end) + y_side(1:end-1));

[xx, yy] = meshgrid(x_cent, y_cent);

u = 1.0;
v = 1.0;
dt = CFL_MAX*dx;

xx_evals = xx(gcw+1:end-gcw, gcw+1:end-gcw);% + 0.5*dx*rand(Nx-1);
yy_evals = yy(gcw+1:end-gcw, gcw+1:end-gcw);% + 0.5*dy*rand(Nx-1);


f = @(x,y) sin(2*pi*x).^4 .* sin(2*pi*y).^4;
%f = @(x,y) exp(-100*((x-0.35).^2+(y-0.35).^2));

% Project solution onto interpolation subspace
h = f(xx,yy);
fvals = 0*xx_evals;

t = 0;
T = 2;
tt = t:dt:(T-dt);
figure(1); clf;
timestep_num = 0;
fig = pcolor(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw, gcw+1:end-gcw), h(gcw+1:end-gcw, gcw+1:end-gcw));
colorbar();
set(fig, 'EdgeColor', 'none');
set(gca, 'fontsize', 16);
title(['t = ', num2str(0.0)], 'fontsize', 16);
pause(0.05);
fig_num = 0;
if (save_figs)
    saveas(gcf, [base_path, num2str(i), '.png'], 'png');
    fig_num = fig_num+1;
end
while (t < T)
    dt = min(dt, T - t);
    fprintf('At beginning of timestep  %i, the time is t = %f\n', timestep_num, t);
    % Fill boundaries for h
    h = fillBoundaries(h(gcw+1:end-gcw, gcw+1:end-gcw), gcw);
    % Compute pre image of particle positions 
    [xx_evals, yy_evals] = integrate_paths(xx, yy, gcw, @(x,y) -(y - 0.5), @(x,y) x-0.5, dt);
    % Now evaluate mapping using values on pre-image
    fvals = evaluateMapping(xx_evals, yy_evals, h, xx, yy, gcw, [0 0], [Lx Ly], [dx dx]);
%    for i = 1:length(xx_evals(gcw+1:end-gcw, gcw+1:end-gcw))
%        for j = 1:length(yy_evals(gcw+1:end-gcw, gcw+1:end-gcw))
%            fvals(i,j) = compute_Z_spline_sum(xx_evals(gcw+i, gcw+j), yy_evals(gcw+i, gcw+j), xx, yy, dx, dy, h, poly_order);
%        end
%    end

    timestep_num = timestep_num + 1;
    if (mod(timestep_num, draw_num) == 0)
        clf;
        fig = pcolor(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw, gcw+1:end-gcw), fvals);
        hold on;
        contour(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw, gcw+1:end-gcw), fvals, 'k');
        colorbar();
        set(fig, 'EdgeColor', 'none');
        title(['t = ', num2str(t+dt)], 'fontsize', 16);
        set(gca, 'fontsize', 16);
        pause(0.05);
        if (save_figs)
            fprintf('saving figure as %s\n', [base_path, num2str(fig_num), '.png']);
            saveas(gcf, [base_path, num2str(fig_num), '.png'], 'png');
            fig_num = fig_num+1;
        end
    end
    % Now update for next iteration
    % Reset lag grid
    %xx_evals = xx(gcw+1:end-gcw, gcw+1:end-gcw);
    %yy_evals = yy(gcw+1:end-gcw, gcw+1:end-gcw);
    % Reset function values
    h(gcw+1:end-gcw, gcw+1:end-gcw) = fvals;
    t = t+dt;
    fprintf('At end of timestep  %i, the time is t = %f\n', timestep_num, t);
end

fig = pcolor(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw, gcw+1:end-gcw), fvals);
colorbar();
set(fig, 'EdgeColor', 'none');
title(['t = ', num2str(t)], 'fontsize', 16);
set(gca, 'fontsize', 16);
pause(0.05);
if (save_figs)
    fprintf('saving figure as %s\n', [base_path, num2str(fig_num), '.png']);
    saveas(gcf, [base_path, num2str(fig_num), '.png'], 'png');
    fig_num = fig_num+1;
end

% Compute error
err = abs(fvals - f(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw,gcw+1:end-gcw)));
figure(2); clf;
fig = pcolor(xx(gcw+1:end-gcw, gcw+1:end-gcw), yy(gcw+1:end-gcw, gcw+1:end-gcw), err);
colorbar();
set(fig, 'EdgeColor', 'none');
L1_err = sum(sum(abs(err)*dx*dy));
L2_err = sqrt(sum(sum(err.^2*dx*dy)));
max_err = max(max(abs(err)));
fprintf("Error norms: \n L1-norm:  %1.12f \n L2-norm:  %1.12f \n max-norm: %1.12f \n", L1_err, L2_err, max_err);