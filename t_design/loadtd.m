function td = loadtd(fn)
% TODO: update this function so that its only input will be the bandlimit
% (as a number) and add properly formated docs
% 
% td = obtd(t, N)   obtains t design from a prexisting file. 
%                   Assumes the file name is 
%                   TD{t=#1,N=#2}.mat and that it 
%                   is saved in the folder TD_files.
% 

path = fullfile(fileparts(which('loadtd.m')), 'tDesigns', fn);
td = importdata(path);

if any(any(isnan(td)))==true
    error('Errors in loading t-designs.');
end
