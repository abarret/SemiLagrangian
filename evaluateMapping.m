function [q] = evaluateMapping(xx_new, yy_new, Q, xx, yy, gcw, x_low, x_up, dx)
q = 0*Q(gcw+1:end-gcw, gcw+1:end-gcw);

% First wrap all xx_new and yy_new points to account for periodic
% boundaries
Lx = x_up(1) - x_low(1);
Ly = x_up(2) - x_low(2);
for i = 1:length(xx_new)
    for j = 1:length(xx_new(i,:))
        while (xx_new(i,j) > x_up(1))
            xx_new(i,j) = xx_new(i,j) - Lx;
        end
        while (xx_new(i,j) < x_low(1))
            xx_new(i,j) = xx_new(i,j) + Lx;
        end
        while (yy_new(i,j) > x_up(2))
            yy_new(i,j) = yy_new(i,j) - Ly;
        end
        while(yy_new(i,j) < x_low(2))
            yy_new(i,j) = yy_new(i,j) + Ly;
        end
    end
end

% Reconstruct function at new indices
for i = 1:length(xx_new)
    for j = 1:length(xx_new(i,:))
        % Determine index in the reference solution
        j_ref = round((xx_new(i,j) - x_low(1)) / dx(1)) + 1;
        i_ref = round((yy_new(i,j) - x_low(2)) / dx(2)) + 1;
        
        % Reconstruct function
        % Collect x and Q information
        index_tot = 9;
        X_vals = zeros(1,index_tot);
        Y_vals = X_vals;
        Q_vals = X_vals;
        index = 1;
        for ii = (i_ref - 1) : (i_ref + 1)
            for jj = (j_ref - 1) : (j_ref + 1)
                X_vals(index) = xx(ii + gcw, jj + gcw);
                Y_vals(index) = yy(ii + gcw, jj + gcw);
                Q_vals(index) = Q(ii + gcw, jj + gcw);
                index = index + 1;
            end
        end
        
        poly_size = 6;
        A = zeros(index_tot, index_tot);
        B = ones(index_tot, poly_size);
        U = zeros(index_tot + poly_size, 1);
        for ii = 1:index_tot
            for jj = 1:index_tot
                X_dist = [X_vals(ii) - X_vals(jj), Y_vals(ii) - Y_vals(jj)];
                A(ii,jj) = rbf(norm(X_dist));
            end
            B(ii, 2) = X_vals(ii);
            B(ii, 3) = Y_vals(ii);
            B(ii, 4) = X_vals(ii)*X_vals(ii);
            B(ii, 5) = Y_vals(ii)*Y_vals(ii);
            B(ii, 6) = X_vals(ii)*Y_vals(ii);
            U(ii, 1) = Q_vals(ii);
        end
        final_mat = [A, B; transpose(B), zeros(poly_size, poly_size)];
        x = linsolve(final_mat, U);
        
        % Now sum solution
        val = 0.0;
        for ii = 1:index_tot
            X_dist = [X_vals(ii) - xx_new(i,j), Y_vals(ii) - yy_new(i,j)];
            val = val + x(ii) * rbf(norm(X_dist));
        end
        val = val + x((index_tot + 1):end)' * [1; xx_new(i,j); yy_new(i,j); xx_new(i,j)^2.0; yy_new(i,j)^2.0; xx_new(i,j) * yy_new(i,j)];
        q(i,j) = val;
    end
end
end

function val = rbf(r)
    val = r * r * r * r * r;
end