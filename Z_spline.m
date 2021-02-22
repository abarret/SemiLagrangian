function val = Z_spline(x, order)
val = x;
switch (order)
    case 1
        for i = 1:length(val)
            xi = abs(x(i));
            if(xi < 1.0)
                val(i) = 1.0 - 2.5 * xi*xi + 1.5 * xi^3.0;
            elseif (xi < 2.0)
                val(i) = 0.5*(2.0 - xi)^2.0*(1.0 - xi);
            else
                val(i) = 0.0;
            end
        end
    case 2
        for i = 1:length(val)
            xi = abs(x(i));
            if(xi < 1.0)
                val(i) = 1.0 - 15.0/12.0*xi*xi - 35.0/12.0*xi^3.0 + 63.0/12.0*xi^4.0 - 25.0/12.0*xi^5.0;
            elseif (xi < 2.0)
                val(i) = -4.0 + 75.0/4.0*xi - 245.0/8.0*xi*xi + 545.0/24.0*xi^3.0 - 63.0/8.0*xi^4.0 + 25.0/24.0*xi^5.0;
            elseif (xi < 3.0)
                val(i) = 18.0 - 153.0/4.0*xi + 255/8.0*xi*xi- 313.0/24.0*xi^3.0 + 21.0/8.0*xi^4.0 - 5.0/24.0*xi^5.0;
            else
                val(i) = 0.0;
            end
        end
    otherwise
        error("Unknown interpolation order");
end
end