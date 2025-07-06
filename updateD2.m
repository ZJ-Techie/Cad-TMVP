function [d, struct_out] = updateD2(beta, struct_in, lnorm)
% --------------------------------------------------------------------
group = 0;
if nargin == 1
    d = 0.5 ./ sqrt(beta.^2 + eps);
    group = sum(abs(beta));
elseif strcmpi(lnorm, 'GGL') %  Graph group lasso
    p = length(beta);
    w = beta .^ 2;
    Gp = struct_in*w+eps;
    Gp = sqrt(Gp);
    Gp = 0.5 ./ Gp;
    d  = sum(reshape(Gp,p-1,[]));
elseif strcmpi(lnorm, 'FGL') % Fused group lasso
    p = length(beta);
    w = beta .^ 2;
    Gp = struct_in * w + eps;
    Gp = sqrt(Gp);
    Gp = 0.5 ./ Gp;
    tmp1 = Gp(1);
    tmpp = Gp(2 * (p - 1));
    Gp = Gp(2 : 2 * (p - 1) - 1);
    d = sum(reshape(Gp, 2, []));
    d = [tmp1 d tmpp];
end
struct_out = group;