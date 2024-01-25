clear all
close all

q = 10; %sensori
k = 2; %attacchi
p = 20; % celle

Tmax = 100;
 e = 1e-8;
 n = 20;

percentage = zeros(6,1);


for k = 0:5
    diff = zeros(n,1);
    correct_est = zeros(n,1);


    for j=1:n
    
        C = randn(q+k,p); 
        
       
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
        
        
        sigma = 1e-2;
        noise = sigma * randn(q+k,1);
        
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
    
    end

percentage(k+1) = sum(correct_est) / n * 100;
end

%Con situazioni diverse aumentando q si ottiene una migliore stima 

