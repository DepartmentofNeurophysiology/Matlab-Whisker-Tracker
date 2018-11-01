function Result = Costumbackground(Objects, Edges)
%%
%Insert costum processing script to extract an object mask (ie further
%processing or load an external file. Output must be a binary matrix with
%similar dimensions as standard frames, non-zero valued entries are
%considered an object
%%

% 
% 
% % The gap crossing task is characterized by its gap, which is measured by
% % taking intensity profiles along the y-axis at the borders. Any pixel
% % within the gap is not an object.
% ncols = 5;
% y1 = sum(Objects(:,1:ncols),2);
% y1 = y1./ncols;
% y2 = sum(Objects(:,end-ncols+1:end),2);
% y2 = y2./ncols;
% y1(1:10) = 1;
% y1(end-10:end) = 1;
% y2(1:10) = 1;
% y2(end-10:end) = 1;
% 
% 
% idx = find(y1 < 0.5 |  y2 < 0.5);
% Objects(idx,:) = 0;
% 
% 
% y1 = sum(Objects, 1)./max(sum(Objects,1));
% y1 = conv(y1, ones(1,10)./10,'same');
% Objects(:, y1 < 0.3) = 0;
% 
% y1 = sum(Objects, 1)./max(sum(Objects,1));
% i1 = find(y1==0,1,'first');
% i2 = find(y1==0,1,'last');
% Objects(:,i1:i2) = 0;
% 
% 


y1 = sum(Edges, 2);
y1 = y1./max(y1);

n = length(y1);
I1 = find(y1(1:round(n/2)) > 0, 1,'last');
I2 = round(n/2) + find(y1(round(n/2):end) > 0, 1,'first');
if isempty(I2)
    I2 = length(y1);
end

if isempty(I1)
    I1 = 1;
end


Objects(I1:I2,:) = 0;
Result.gapinfo.main_ax = 2;
Result.gapinfo.edge_1 = I1;
Result.gapinfo.edge_2 = I2;
Result.Objects = Objects;
Result.Edges = Edges;

