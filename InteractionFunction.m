function[y] = InteractionFunction(ftype,a,b,c,d,x)

if ftype == 1 %Linear
    z=a*x+b;
    
elseif ftype == 2 %Logistic
    z = a / (1 + exp(-b * (c - x))) + d;
    
elseif ftype==3 %Polynomial
    z=a+b*x+(c*x^2)+(d*x^3);
    
elseif ftype==4 %Gaussian
    z = a * exp(-((x - b)^2) / (2 * c^2)) + d;
    
elseif ftype == 5 %Slope-only, with y=1 when x=1
    b = 1 - a;
    z = a * x + b;
    
elseif ftype == 6 %Slope-only, with y=0 when x=1
    b = -a;
    z = a * x + b;
    
elseif ftype == 7 %Slope-only, with y=0 when x=0
    z = a * x;
    
elseif ftype == 8 %Slope-only, with y=0 when x=1; with max y=0
    b = -a;
    z0 = a * x + b;
    z = min(0, z0);
    
elseif ftype == 9 %Gaussian with death threshold
    test_z = a * exp(-((x - b)^2) / (2 * c^2)) + d;
    z = test_z;
    if test_z <= -1
        z = -100;
    end
    
elseif ftype == 10 %constant
    z = a;
    
elseif ftype == 11 %Slope-only, with y=0 when x=1; with max y=0
    b = -a;
    z0 = a * x + b;
    z = min(0, z0);
    
elseif ftype == 12 %Slope-only (slope = a), with y=0 when F=b
    z = (x - b) * a;
    
end

y = z;

end

