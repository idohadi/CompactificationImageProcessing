%% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2022
% ***********************************************************

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['ClassificationRepresentatives-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
repsToTake = [1, 2, 3];

%% Produce figure
fig = figure;
t = tiledlayout(ceil(length(repsToTake)/3), 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

caxisvals = classRepresentatives(:, :, repsToTake);
caxismin = min(caxisvals, [], 'all');
caxismax = max(caxisvals, [], 'all');

for I=repsToTake
    nexttile;
    imagesc(classRepresentatives(:, :, I));
    colormap('hot');
    colorbar;
    caxis([caxismin, caxismax]);
    title(num2str(I, 'Representative %d'))
end

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
