maxMin = zeros(2000,3)
numMaxMin = 1;
for k = 1:size(v)
    m = mod(k, 10000);
    if(m == 0)
        numMaxMin = numMaxMin + 1;
    end
    
    if(maxMin(numMaxMin, 1) < v(k))    %Finds the max for the time period
       maxMin(numMaxMin, 1) = v(k); 
    end
    
    if(maxMin(numMaxMin, 2) > v(k))    %Finds the min for the time period
       maxMin(numMaxMin, 2) = v(k); 
    end
end
maxMin(numMaxMin + 1:end, :) = [];


for k = 1:size(maxMin)
    maxMin(k, 3) = maxMin(k, 1) - maxMin(k, 2);
end
averageDifference = mean(maxMin); %averageDifference(3) is the one that holds the real value
display(averageDifference(3))
