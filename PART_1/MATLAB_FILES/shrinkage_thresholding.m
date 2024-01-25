function [y]=shrinkage_thresholding(x,gamma,p)
   
    y=zeros(p,1);

    for i=1:size(x)

        if x(i) > gamma(i)
            y(i) = x(i) - gamma(i);
        
        elseif x(i) < -gamma

            y(i) = x(i)+gamma(i);
        
        else

            y(i)=0;

        end
    
    end

    
end