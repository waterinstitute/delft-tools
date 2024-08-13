function [T, Ti] = getTri(faceNodes)


%FACENODES is an mxn matrix. where n is the number of faces in a mesh with
%mixed triangles and other shapes. m is the maximum number of edges on any
%face. Vertices for each face are listed, clockwise, in the columns.

%T is a triangulation, with larger shapes split into triangles
%Ti is the corresponding face indices for each triangle

numFaces=size(faceNodes,2); %number of faces
sides=sum(~isnan(faceNodes),1); %number of edges on each face

ii3=find(sides==3); %indexes of triangle faces
ii4=find(sides==4); %indexes of quad faces
newTris=numel(ii3)+2*numel(ii4);    %number of triangles in the new triangulation

T=NaN(3,newTris);   %vertexes of new triangulation
Ti=NaN(1,newTris);  %corresponding face index in mesh for each new triangle

ct=0;   %counter for triangles in new triangulation
for i=1:numFaces
    if sides(i)==3 %triangle face
        ct=ct+1;
        T(1:3,ct)=faceNodes([1 2 3],i);
        Ti(ct)=i;

    elseif sides(i)==4  %quad face
        ct=ct+1;
        T(1:3,ct)=faceNodes([1 2 3],i);
        Ti(ct)=i;
        
        ct=ct+1;
        T(1:3,ct)=faceNodes([3 4 1],i);
        Ti(ct)=i;
    else
        error('does not yet handle shapes larger than quads')
        %make sure to allocate the size of T for the larger number of
        %triangles that would be needed in this case.  
    end
end
