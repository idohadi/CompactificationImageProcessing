function [bispInds, symComp] = utilityFunc(L)
% Calculate the bispetral invariant indices that are non-zero
bispInds = [];
c = 1;
for l1=0:L
    for l2=0:l1
        for l=abs(l1-l2):min(l1+l2, L)
            if mod(l1+l2-l, 2) == 0
                bispInds(c) = 1;
                c = c + 1;
                bispInds(c) = 0;
                c = c + 1;
            else
                bispInds(c) = 0;
                c = c + 1;
                bispInds(c) = 1;
                c = c + 1;
            end
        end
    end
end
bispInds = logical(bispInds);

% Calcualte the derivative of the mapping from {f_{l,m}} (modulo 
% symmetries) to {f_{l,m}}
symComp = sparse(2*(L+1)^2, (L+1)^2);

bandlimit = L;
for order=0:bandlimit
    for degree=0:order
        if degree==0
            symComp(2*(order^2+order)+1, order^2+1) = 1;
        else
            % Handling +degree
            % Real part of coefficient
            symComp(2*(order^2+order+degree)+1, order^2+1+2*(degree-1)+1) = 1;
            % Imaginary part of coefficeint
            symComp(2*(order^2+order+degree)+2, order^2+1+2*(degree-1)+2) = 1;
            
            % Handling -degree
            % Real part of coefficient
            symComp(2*(order^2+order-degree)+1, order^2+1+2*(degree-1)+1) = (-1)^degree;
            % Imaginary part of coefficeint
            symComp(2*(order^2+order-degree)+2, order^2+1+2*(degree-1)+2) = (-1)^(degree+1);
        end
    end
end
end
