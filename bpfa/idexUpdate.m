
function idexNew = idexUpdate(sizeIMin, PatchSize,colj,rowi)
%Update the index of the patches in the training data
%Version 1: 09/12/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
idexMat=zeros(sizeIMin-PatchSize+1);
if colj==1 && rowi==1
    idexMat([rowi:PatchSize:end-1,end],[colj:PatchSize:end-1,end])=1;
elseif colj==1 && rowi~=1
    idexMat(rowi:PatchSize:end,[colj:PatchSize:end-1,end])=1;
elseif colj~=1 && rowi==1
    idexMat([rowi:PatchSize:end-1,end],colj:PatchSize:end)=1;
else
    idexMat(rowi:PatchSize:end,colj:PatchSize:end)=1;
end
idexNew = find(idexMat);
end