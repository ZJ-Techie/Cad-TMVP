function [E] = updateGraph2(p, flag)
% --------------------------------------------------------------------
% Update the weight matrix 
% -----------------------------------------
ncol = p;
iedge = 0;
switch flag
    case 'GGL'
        E = zeros(p*(p-1),p);
        for i = 1:p
            for j = 1:p
                if i ~= j
                    iedge = iedge+1;
                    E(iedge,i) = 1;
                    E(iedge,j) = 1;
                end
            end
        end
    case 'FGL' % Fused group lasso
        E = zeros(2*(p-1),p);
        for i = 1:p-1
            j = i+1;
            E(2*i-1,i) = 1;
            E(2*i-1,j) = 1;
            E(2*i,i) = 1;
            E(2*i,j) = 1;
        end
end