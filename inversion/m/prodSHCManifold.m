function M = prodSHCManifold(L)
mans = struct();
for l=1:L
    mans.(['M', num2str(l)]) = spherefactory(2*l+1, 1);
end
M = productmanifold(mans);
