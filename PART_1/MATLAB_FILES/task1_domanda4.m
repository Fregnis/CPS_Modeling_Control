clear all
close all
clc

q = 10; %sensori
k = 2; %attacchi
p = 20; % celle

Tmax = 10000;
 e = 1e-8;
 n = 20;



num_iters = zeros(n,1);
idx = 1;
taus = zeros(0.016/0.0005,1);
perc = zeros(0.016/0.0005,1);
for tau=0.016:-0.0005:0

    correct_est = zeros(n,1);
    diff = zeros(n,1);
    for j=1:n
    
        C = randn(q,p); 
        
       
        %tau = norm(C)^-2 -e-0.02;
        
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
        noise = sigma * randn * ones(q,1);
        
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
    
        num_iters(j) = i;
        diff(j)= norm(x_reale - x_hat)^2;
    
        if diff(j) <=1 
            correct_est(j) = 1;
        end
    
    
        
    end
    
    result = double([mean(num_iters) max(num_iters) min(num_iters)]);
    percentage = sum(correct_est) / n * 100;
    
   taus(idx) = tau;
   perc(idx) = percentage;
   idx = idx +1;
end

figure(1)
plot(taus,perc,'-x')
grid on


%se tau cambia segno è impossibile avere una stima corretta perchè nelle
%condizioni dell'IST il vincolo sul valore assoluto non si verifica mai
%quindi non avremo mai un vettore sparso

%se tau è positivo è molto piccolo influisce sulla stima perchè ad ogni
%passo togliamo una piccola quantità ad x e quindi ci vuole più tempo per
%la convergenza al valore reale

