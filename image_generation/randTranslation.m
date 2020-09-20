function [translation, translationAngle, translationSize] ...
    = randTranslation(maxPixels, lambda)
% TODO: Docs

if nargin==0
    maxPixels = 5;
    lambda = 1.5;
end

translationAngle = 2*pi*rand();
translationSize = Inf;
while translationSize>maxPixels
    translationSize = exprnd(lambda);
end
translationSize = (2*randi(2)-3)*translationSize;
translation ...
    = translationSize*[cos(translationAngle), sin(translationAngle)];
