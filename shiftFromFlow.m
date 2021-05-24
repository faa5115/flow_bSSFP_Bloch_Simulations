function Mb = shiftFromFlow(Ma, NShift)

%Ma:  length along slice x 3.  the isochromat arrangement along the initial
%slice. 
%NShift:  simply the number of elements that each isochromats shift. 

%Mb:  the new isochromat arrangement.  
Mb = zeros(size(Ma));
if(NShift >= 0)
    for n = size(Ma,1) : -1 : 1

        if (n <= NShift)
            Mb(n, :) = [0 0 1];
        else
            Mb(n, :) = Ma(n-NShift, :);
        end

    end
    
else
    for n = 1 : 1 : size(Ma, 1)
        if (n >= (size(Ma,1) + 1 - -NShift) )
            Mb(n, :) = [0 0 1];
        else
            Mb(n, :) = Ma(n - NShift, :);
        end
    end
end