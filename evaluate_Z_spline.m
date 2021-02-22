function val = evaluate_Z_spline(x, y, dx, dy, order)
val = Z_spline(x / dx, order) * Z_spline(y / dy, order);
end