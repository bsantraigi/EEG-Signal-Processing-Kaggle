function [ fselect ] = GetTestFiles( user )
%GETTESTFILES Summary of this function goes here
%   Detailed explanation goes here
fileLimit = 3000;
fselect = {};
    
h = waitbar(0, '', 'Name', 'Looking for files...');
steps = fileLimit;

for i = 1:fileLimit
    fname = sprintf('data/test_%d_new/new_%d_%d.mat', user, user, i);
    waitbar(i / steps, h, sprintf('data/test-%d-new/new-%d-%d.mat', user, user, i));
    if(exist(fname) == 2)
        fselect = [fselect {fname}];
    end
end
newData = zeros(80, 29*length(fselect));

close(h);

end

