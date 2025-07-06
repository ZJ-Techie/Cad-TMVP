function x = soft_Threshold(x,lambda)
p = length(x);
for i = 1:p
    if x(i) > lambda
        x(i) = x(i)-lambda;
    elseif x(i) < -lambda
        x(i) = x(i)+lambda;
    else
        x(i) = 0;
    end
end