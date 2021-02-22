function val = compute_Z_spline_sum(x, y, xx, yy, dx, dy, h, order)
% Determine width of box
width = get_width_of_spline(order)+1;

x_min = x - width*dx;
x_max = x + width*dx;
y_min = y - width*dy;
y_max = y + width*dy;

% Get indices that we need to sum over
indices = find(xx > x_min & xx < x_max & yy > y_min & yy < y_max);


val = 0.0;
for index = 1 : length(indices)
    i = mod(indices(index), length(xx))+1;
    j = fix(indices(index) / length(xx)) + 1;
    val = val + h(i,j) * evaluate_Z_spline(x-xx(i,j), y-yy(i,j), dx, dy, order);
end
end