
function X = XUpdate_MissingData(IMin,PatchSize,idex,MuX,Yflag,ImageType)
%Update the training data taken into consideration
%Version 1: 09/12/2009
%Version 2: 10/26/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
if strcmp(ImageType,'rgb')==1
    X =  [im2col(IMin(:,:,1)/255,[PatchSize,PatchSize],'sliding');
        im2col(IMin(:,:,2)/255,[PatchSize,PatchSize],'sliding');
        im2col(IMin(:,:,3)/255,[PatchSize,PatchSize],'sliding')];
    X = X(:,idex).*Yflag;
    X = X - repmat(MuX,1,length(idex)).*Yflag;
else
    X = im2col(IMin/255,[PatchSize,PatchSize],'sliding');
    %X = X(:,idex)-repmat(MuX,1,length(idex));
    X = X(:,idex).*Yflag-MuX.*Yflag;    
end