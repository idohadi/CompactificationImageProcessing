function complexCoeffs = r2c(realCoeffs, bandlimit)
complexCoeffs = zeros(2*(bandlimit+1)^2, 1);

for order=0:bandlimit
    for degree=0:order
        if degree==0
            complexCoeffs(2*(order^2+order)+1) ...
                = realCoeffs(order^2+1);
        else
            % Handling +degree
            % Real part of coefficient
            complexCoeffs(2*(order^2+order+degree)+1) ...
                = realCoeffs(order^2+1+2*(degree-1)+1);
            % Imaginary part of coefficeint
            complexCoeffs(2*(order^2+order+degree)+2) ...
                = realCoeffs(order^2+1+2*(degree-1)+2);
            
            % Handling -degree
            % Real part of coefficient
            complexCoeffs(2*(order^2+order-degree)+1) ...
                = (-1)^degree * realCoeffs(order^2+1+2*(degree-1)+1);
            % Imaginary part of coefficeint
            complexCoeffs(2*(order^2+order-degree)+2) ...
                = (-1)^(degree+1) * realCoeffs(order^2+1+2*(degree-1)+2);
        end
    end
end

complexCoeffs = reshape(complexCoeffs, [2, (bandlimit+1)^2]).';
complexCoeffs = complexCoeffs(:, 1) + 1i*complexCoeffs(:, 2);

end
