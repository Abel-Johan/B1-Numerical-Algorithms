% QN 1
function [fit, condition] = regress1(x, y, M, condition)
    % Fits and plots Mth order polynomial to data points (x, y)
    % Plot fitted data at points described by x1

    % There is a MATLAB function vander(v) which automatically returns a
    % Vandermonde matrix based on x-points 'v' but that is only for square
    % Vandermondes (thus order of polynomial is fixed to length(x) - 1).
    % Here we allow the user to choose the order of polynomial themselves.

    V = zeros(length(x), M + 1);
    % Populate Vandermonde matrix by taking integer powers of x along the
    % columns
    for i = 0:M
        V(:, i+1) = x.^i;
    end

    V_inv = inv(V'*V)*V';   % Pseudoinverse of Vandermonde

    condition(M) = cond(V'*V);

    % Compute polynomial coefficients
    % y is by default a row vector. Turn it into column vector
    alpha = V_inv*y;       

    fit = V*alpha;

end

% Generate synthetic data y = sin(omega*t + phi)
t = [0:0.001:2*pi]';
omega = 2.0;
y = sin(omega*t + phi);
noise = normrnd(0, sigma, [1, length(t)]);
y_noise = y + noise';

% Initialise array to store polynomial fit data and conditioning
M = 10;     % max polynomial order
fits = zeros(M, length(t));
condition = zeros(1, M);

figure
tiledlayout(5, 2)
for m = 1:M
    [fits(m, :), condition] = regress1(t, y_noise, m, condition);

    nexttile
    plot(t, fits(m,:), 'r', t , y_noise, 'b')
    legend('fitted poly', 'original noisy y')
    xlabel("t")
    ylabel("y")
end

figure
semilogy(1:M, condition, '-x')
xlabel("polynomial order m")
ylabel("condition number (log-axis)")
% QN 1

% QN 2



% QN 2