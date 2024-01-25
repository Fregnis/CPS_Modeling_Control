clear all
close all

q = 10; %sensori
k = 2; %attacchi
p = 20; % celle

Tmax = 100;
 e = 1e-8;
 n = 20;
diff = zeros(5,1);
correct_est = zeros(5,1);

C = randn(q,p); 
 sigma = 1e-2;
noise = sigma * randn(q,1);
percentage = zeros(1,5);

for j=1:5

    C=[C;randn(1,p)];
    
    tau = norm(C)^-2 - e;
    
    lambda = 1/(100*tau);
    
    Lambda = lambda * ones(p,1);
    
    S_x =randperm(q,2);
    
    x_reale = zeros(p,1);
    
    for i=1:length(S_x)
    
        x_reale(S_x(i)) = unifrnd(-1,1);
    
        if x_reale(S_x(i)) > 0
            x_reale(S_x(i)) = x_reale(S_x(i)) + 1;
        else
            x_reale(S_x(i)) = x_reale(S_x(i)) -1;
        end
    
    
    end
    
    
   
    noise = [noise;sigma * randn * ones(1,1)];
    
    x_hat = zeros(p,1);
    
   
    gamma = tau *Lambda;
    y = C*x_reale + noise;
    
    for i = 1:Tmax
    
    
        x_prev = x_hat;
        x_hat = shrinkage_thresholding(x_prev + tau * C' * (y - C * x_prev),gamma,p); 
    
        if norm(x_hat - x_prev) < 1e-12
            break;
        end
    
    end

    diff(j)= norm(x_reale - x_hat)^2;

    if diff(j) <=1 
        correct_est(j) = 1;
    end
percentage(j) = sum(correct_est) / 5 * 100;
end

%conclusione 1
% Avendo stesse misurazioni e aumentando q non è detto che la stima migliori poichè le misurazioni aggiunte possono essere affette 
% da errori