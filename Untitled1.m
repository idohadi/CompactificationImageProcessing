L = 100;

c = 0;
for l1=0:L
    for l2=0:l1
        for l=abs(l1-l2):min(l1+l2, L)
            for m=-l:l
                c = c + 1;
            end
        end
    end
end

tab = zeros(c, 6);
c = 0;
c1 = 0;
for l1=0:L
    for l2=0:l1
        for l=abs(l1-l2):min(l1+l2, L)
            c1 = c1 + 1;
            for m=-l:l
                c = c + 1;
                nmin = max(-l1, m-l2);
                nmax = min(l1, m+l2);
                tab(c, :) = [nmax-nmin+1, l1, l2, l, m, c1];
            end
        end
    end
end

tab = sortrows(tab, [1, 2, 3, 4, 5, 6]);

blockSize = 64;
mem = zeros(ceil(size(tab, 1)/blockSize), 1);
for J=1:size(mem, 1)
    mem(J) = 4*6*blockSize;
    for M=2:3
        U = unique(tab(blockSize*J:min(blockSize*(J+1), size(tab, 1)), M));
        mem(J) = mem(J) + 4*length(U)*2*(2*L+1);
    end
end

% evDer = zeros(2*L+1, 1);
% for J=1:(2*L+1)
%     inds = tab(:, 1)==J;
%     tmp = sort(unique(tab(inds, 2)));
%     evDer(J) = all(tmp==(tmp(1):tmp(end))');
% end
% 
% evDer2 = cell(2*L+1, 1);
% for J=1:(2*L+1)
%     inds = tab(:, 1)==J;
%     tmp = sort(unique(tab(inds, 2)));
%     evDer2{J} = zeros(tmp(end)-tmp(1)+1, 1);
%     for M=tmp(1):tmp(end)
%         inds2 = tab(:, 2)==M;
%         inds2 = inds2 & inds;
%         tmp2 = sort(unique(tab(inds2, 3)));
%         evDer2{J}(M-tmp(1)+1) = all(tmp2==(tmp2(1):tmp2(end))');
%     end
% end
