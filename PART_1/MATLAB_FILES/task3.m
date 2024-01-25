clear all
close all
clc

D = -[46 58 76 69 67 80 61;
      54 49 72 63 56 65 59;
      52 50 56 58 58 62 42;
      61 62 49 61 60 65 44;
      73 78 65 69 55 57 61;
      72 65 69 47 53 44 62];

y= -[62 58 41 46 64 63]';

q = 6;
p = 7;

Tmax = 5000; % lo aumentiamo per avere dei risultati più accurati in quanto in quetso modo la stima converge al risultato reale

counter = 0;
norms = zeros(20,1); 

 lambda= 1;
  G = [D eye(q)];
  G = normalize(G);
  tau = norm(G)^-2 - 1e-8;
  

 tau_Lambda = tau*lambda*ones(q+p,1);


  z_prev = zeros(q+p,1);
  z_hat = zeros(q+p,1);

  h = 1;

  S_a=randperm(q,1);

  y(S_a) =y(S_a) - y(S_a)*1/5;
  

for i = 1:Tmax
    
    
        z_prev = z_hat;
        z_hat = shrinkage_thresholding(z_prev + tau * G' * (y - G * z_prev),tau_Lambda,q+p); 
    
        if norm(z_hat(1:p) - z_prev(1:p)) < 1e-12
            break;
        end
   
end

  x_hat = z_hat(1:p);

  targets = find(x_hat);

  w = find(lasso(D,y,'Lambda',lambda));

% modificando la y con 1/5*y è tantissima la differenza quindi non
% riconosce perfettamente i target

% a seconda di quale sensore è preso sotto attacco la stima dei target può
% essere più o meno sbagliata

