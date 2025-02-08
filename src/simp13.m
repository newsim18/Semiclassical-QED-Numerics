function simpson13_integral = simp13(elements,d)
    % 1/3 Simpson's method
    
    f1 = 0;
    f2 = 0;
        
    % Calculating the odd k terms
    for j = 2:2:length(elements(1,:))-1
        f1 = f1 + elements(1,j);
    end
    
    % Calculating the even k terms    
    for j = 3:2:length(elements(1,:))-2
        f2 = f2 + elements(1,j);
    end
    
    simpson13_integral = (d/3)*(elements(1,1) + 4*f1 + 2*f2 + elements(1,length(elements(1,:))));