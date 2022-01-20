%% File setup
fnNOEXT = 'ClassificationCryoFollow-2022-01-20-1422';
fn = [fnNOEXT, '.mat'];
k = 50;
diary([fnNOEXT, '.log']); % Log file

save(fn, 'k', '-v7.3');

%% Follow up computation
printBegEndMsg('Running test', true);

printBegEndMsg('Finding k nearest neighrbors', true);
[B, I] = mink(distanceMatrix, k+1);
save(fn, 'B', 'I', '-append');
printBegEndMsg('Finding k nearest neighrbors', false);
printBegEndMsg('Running test', false);

%% Calcualting node specificity

printBegEndMsg('Calculating node specificity', true);
classSpecificity = zeros(sampleSize, 1);
for J=1:sampleSize
    sameClassNo = classMembership(I(:, J))==classMembership(J);
    classSpecificity(J) = (sum(sameClassNo)-1)/k;
end
printBegEndMsg('Calculating node specificity', false);
save(fn, 'classSpecificity', '-append');

%% Produce figure
fig = figure;

bins = 0:0.025:1;
fa = 0.5;
ea = 0.6;

histogram(classSpecificity, bins, 'FaceAlpha', fa, 'EdgeAlpha', ea);
set(gca, 'yscale', 'log');

ylim([10^0, 10^4]);
xlim([0, 1]);
xticks(0:0.2:1);
xlabel('Node score');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;