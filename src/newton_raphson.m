function [ X y ] = newton_raphson( f, X0, rate, tol, varargin )
%NEWTON_RAPHSON Newton-Raphson optimization for scalar f
%   [X y] = newton_raphson(f, X0, tol, varargin) returns a stationary point
%   of f at X and its value y = f(X) for a scalar function f. X0 is a
%   vector of initial guesses. Extra arguments to varargin are passed to f.
%   The function returns when norm(Xnew - Xold) < tol.

while true
    [y g H] = f(X0, varargin{:});
    X = X0 - rate * H \ g
    if norm(X - X0) < tol
        return
    else
        X0 = X;
    end
end


end

