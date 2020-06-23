function realCoeffs = c2r(complexCoeffs, bandlimit)
complexCoeffs = reshape([real(complexCoeffs), imag(complexCoeffs)].', [2*length(complexCoeffs), 1]);

realCoeffs = zeros((bandlimit+1)^2, 1);
for order=0:bandlimit
    for degree=0:order
        if degree==0
            realCoeffs(order^2+1) ...
                = complexCoeffs(2*(order^2+order)+1);
        else
            % Handling +degree
            % Real part of coefficient
            realCoeffs(order^2+1+2*(degree-1)+1) ...
                = complexCoeffs(2*(order^2+order+degree)+1);
            
            % Imaginary part of coefficeint
             realCoeffs(order^2+1+2*(degree-1)+2) ...
                 = complexCoeffs(2*(order^2+order+degree)+2);
            
        end
    end
end
end
