function [W, u, v, w, obj] = FastCadTMVP(Data, opts)
% --------------------------------------------------------------------
% example code for FAST Causality-driven Trustworthy Multi-View maPping approach (FAST Cad-TMVP)
%------------------------------------------
% Author: jinzhang@mail.nwpu.edu.cn
% Date created: 06-07-2025
% @Northwestern Ploytechnical University.
% -----------------------------------------% --------------------------------------------------------------------
% data
X = Data.X{1}; p = size(X, 2);
Y = Data.X{2}; q = size(Y, 2);
Z = Data.X{3}; r = size(Z, 2);
n = size(X, 1);
DX = Data.DX;

X = getNormalization(X, 'normalize');
Y = getNormalization(Y, 'normalize');
Z = getNormalization(Z, 'normalize');

XX = X'* X;
YY = Y'* Y;
ZZ = Z'* Z;

u0 = ones(p, 1);
v0 = ones(q, 1);
w0 = ones(r, 1);
u = u0;
v = v0;
w = w0;

% scale u v and w
scale = sqrt(u' * XX * u);
u = u ./ scale;
scale = sqrt(v' * YY * v);
v = v ./ scale;
scale = sqrt(w' * ZZ * w);
w = w ./ scale;
iu = 1; iv = 2; iw = 3;
W{iu} = u; W{iv} = v; W{iw} = w;

lambda1 = opts.lambda1; 
lambda2 = opts.lambda2; 

max_Iter = 100;
i = 0;
tol = 1e-5;
obj = [];
tv = inf;
tu = inf;
tw = inf;
res1 = [];
res2 = [];
res3 = [];

Q = ones(n,1)./n;
lambda0 = 0.1;
% lambda1 = 0.1;
% lambda2 = 0.1;
lambda3 = 0.1;
lambda4 = 0.1;
lambda5 = 0.1;
lambda_Q = 0.001;
alpha = 0.1;

% get the structure
Gu = updateGraph2(p,'FGL');
Gv = updateGraph2(q,'GGL');
Gw = updateGraph2(r,'GGL');

block = 20;
Xall = X;
Yall = Y;
Zall = Z;

% initialize loss weights
Xu = X * u; Yv = Y * v; Zw = Z * w;
gamma = 2;    
p_c12 =  abs(corr(Xu, Yv));
p_c13 =  abs(corr(Xu, Zw));
p_c23 =  abs(corr(Yv, Zw));
omega12 = - ((1 - p_c12)^gamma) * log(p_c12);
omega13 = - ((1 - p_c13)^gamma) * log(p_c13);
omega23 = - ((1 - p_c23)^gamma) * log(p_c23);

while (i < max_Iter && tu > tol && tv > tol && tw > tol) % default 100 times of iteration
    i = i + 1;
    % update u
    % -------------------------------------
    nblock = ceil(p/block);
    ut = [];
    s1 = 0;
    u_old = u;
    for iu = 1 : nblock
        if iu * block <= p
            X = Xall(:, 1 + ( iu - 1) * block : iu * block);
            sub_u = u(1 + (iu - 1) * block : iu * block);
        else
            X = Xall(:, 1 + (iu - 1) * block : end);
            sub_u = u(1 + (iu - 1) * block : end);
        end
        XX = X' * X;
        % D1 = updateD(sub_u);
        % D1 = diag(D1);
        D1 = updateD2(sub_u);
        D1 = diag(D1);
        pp = size(X, 2);
        % dfglu = updateD(sub_u, 'FGL');
        % DFGLu = diag(dfglu);
        Gu = updateGraph2(pp,'FGL'); % FGL
        dfglu = updateD2(sub_u, Gu, 'FGL'); % calculate FGL gradient
        DFGLu = diag(dfglu);
        
        F1 = (omega12 + omega13) *XX + lambda1 * D1 + lambda2 * DFGLu;
        b1 = omega12 * X' * diag(Q) * Yall  * v + omega13 * X' * diag(Q) * Zall * w + alpha * X' * (DX - Yall * v - Zall * w);
        sub_u = F1 \ b1;
        ut = [ut; sub_u];
        s1 = s1 + sub_u' * XX * sub_u;
    end
    % scale u
    u = ut./s1;
%   Xu = Xall * u;
    
    % -------------------------------------
    nblock = ceil(q/block);
    vt = [];
    s1 = 0;
    v_old = v;
    for iv = 1 : nblock
        if iv * block <= q
            Y = Yall(:, 1 + (iv - 1) * block : iv * block);
            sub_v = v(1 + (iv - 1) * block : iv * block);
        else
            Y = Yall(:, 1 + (iv - 1) * block : end);
            sub_v = v(1 + (iv - 1) * block : end);
        end
        YY = Y' * Y;
        
        D1 = updateD2(sub_v);
        D1 = diag(D1);
        qq = size(Y, 2);
        Gv = updateGraph2(qq,'GGL'); % GGL
        dgglv = updateD2(sub_v, Gv, 'GGL');% calculate GGL gradient
        DGGLv = diag(dgglv);

        F2 = (omega12 + omega23) * YY + lambda1 * D1 + lambda2 * DGGLv;
        b2 =  omega12 * Y' * diag(Q) * Xall  * u +  omega23 * Y' * diag(Q) * Zall * w + alpha * Y' * (DX - Xall * u - Zall * w);
        sub_v = F2\b2;
        vt = [vt; sub_v];
        s1 = s1 + sub_v' * YY * sub_v;
    end
    
      v = vt./s1;
%     Yv = Yall * v;

    % -------------------------------------
    nblock = ceil(r/block);
    wt = [];
    s1 = 0;
    w_old = w;
    for iw = 1 : nblock
        if iw * block <= q
            Z = Zall(:, 1 + (iw - 1) * block : iw * block);
            sub_w = w(1 + (iw - 1) * block : iw * block);
        else
            Z = Zall(:, 1 + (iw - 1) * block : end);
            sub_w = w(1 + (iw - 1) * block : end);
        end
        ZZ = Z' * Z;
      
        D1 = updateD2(sub_w);
        D1 = diag(D1);
        rr = size(Z, 2);
        Gw = updateGraph2(rr, 'GGL'); % GGL
        dgglw = updateD2(sub_w, Gw, 'GGL');% calculate GGL gradient
        DGGLw = diag(dgglw);
        
        F3 = (omega13 + omega23) * ZZ + lambda1 * D1 + lambda2 * DGGLw;
        b3 = omega13 * Z' * diag(Q) * Xall  * u +  omega23 * Z' * diag(Q) * Yall * v + alpha * Z' * (DX - Xall * u - Yall * v);
        sub_w = F3\b3;
        wt = [wt; sub_w];
        s1 = s1 + sub_w' * ZZ * sub_w;
    end
    % scale w
    w = wt./s1;
%   Zw = Zall * w;

    Xu = Xall * u; Yv = Yall * v; Zw = Zall * w;
    gamma = 2;    
    p_c12 =  abs(corr(Xu, Yv));
    p_c13 =  abs(corr(Xu, Zw));
    p_c23 =  abs(corr(Yv, Zw));
    omega12 = - ((1 - p_c12)^gamma) * log(p_c12);
    omega13 = - ((1 - p_c13)^gamma) * log(p_c13);
    omega23 = - ((1 - p_c23)^gamma) * log(p_c23);

    obju = omega12 * norm(Xall * u - Yall * v)^ 2 + omega13 * norm(Xall * u - Zall * w)^ 2;
    objv = omega12 * norm(Xall * u - Yall * v)^ 2 + omega23 * norm(Yall * v - Zall * w)^ 2;
    objw = omega13 * norm(Xall * u - Zall * w)^ 2 + omega23 * norm(Yall * v - Zall * w)^ 2;

%     obju = -omega12 * u' * Xall'*Yall * v - omega13 * u' * Xall'*Zall * w;
%     objv = -omega12 * u' * Xall'*Yall * v - omega23 * v' * Yall'*Zall * w;
%     objw = -omega13 * u' * Xall'*Zall * w - omega23 * v' * Yall'*Zall * w;

    grad_Q = 2 * lambda0 * obju.*Q...
            + lambda3 * balance_grad(Q, Xall) * ones(p,1)...
            + 2 * lambda0 * objv.*Q...
            + lambda3 * balance_grad(Q, Yall) * ones(q,1)...
            + 2 * lambda0 * objw.*Q...
            + lambda3 * balance_grad(Q, Zall) * ones(r,1)...
            + 4 * lambda4 * Q.*Q.*Q...           
            + 4 * lambda5 * (sum(Q.*Q) - 1) * Q;
        Q = Q - lambda_Q * grad_Q;
        Q = Q./sqrt(sum(Q.*Q));
        Q = Q.*Q;  
    % ------------------------------
    % stopping condition
    if i > 1
        tu = max(abs(u - u_old));
        tv = max(abs(v - v_old));
        tw = max(abs(w - w_old));
    else
        tu = tol * 10;
        tv = tol * 10;
        tw = tol * 10;
    end
end

W{iu} = u; W{iv} = v; W{iw} = w;

end

function g_w = balance_grad(Q, X)
    n = size(X, 1); 
    m = size(X, 2); 
    g_w = zeros(n, m); 
            for i = 1:m
                for j =1:m
                    J1 = ((Q.*Q)'*X(:,i).*X(:,j))./n...
                        -X(:,i)'*(Q.*Q)*X(:,j)'*(Q.*Q)./(n^2);
                    dJ1Q = (X(:,i).*X(:,j))./n-(X(:,i)*X(:,j)'*(Q.*Q)+ X(:,j)*X(:,i)'*(Q.*Q))./(n^2);
                    dJ1Q =  dJ1Q.*Q;
                    g_w(:,i) = 2*dJ1Q'*J1;
                end
            end           
end

% function g_w = balance_grad(Q, X)
%     n = size(X, 1); 
%     m = size(X, 2); 
%     g_w = zeros(n, m); 
%     block = 10;
%     nblock = ceil(m/block);
%     for iu = 1:nblock
%         if iu*block <= m
%             for i = (1 + (iu - 1) * block : iu * block)
%                 for j =( 1 + (iu - 1) * block : iu * block)
%                     J1 = ((Q.*Q)'*X(:,i).*X(:,j))./n...
%                         -X(:,i)'*(Q.*Q)*X(:,j)'*(Q.*Q)./(n^2);
%                     dJ1Q = (X(:,i).*X(:,j))./n-(X(:,i)*X(:,j)'*(Q.*Q) + X(:,j)*X(:,i)'*(Q.*Q))./(n^2);
%                     dJ1Q =  dJ1Q.*Q;
%                     g_w(:,i) = 2 * dJ1Q' * J1;
%                 end
%             end        
%         else
%              for i = (1 + (iu - 1) * block : m)
%                 for j = (1 + (iu - 1) * block : m)               
%                     J1 = ((Q.*Q)'*X(:,i).*X(:,j))./n...
%                         -X(:,i)'*(Q.*Q)*X(:,j)'*(Q.*Q)./(n^2);                    
%                     dJ1Q = (X(:,i).*X(:,j))./n-(X(:,i)*X(:,j)'*(Q.*Q)+ X(:,j)*X(:,i)'*(Q.*Q))./(n^2);                   
%                     dJ1Q =  dJ1Q.*Q;
%                     g_w(:,i) = 2 * dJ1Q' * J1;
%                 end
%             end
%                  
%         end
%         
%     end
% end
