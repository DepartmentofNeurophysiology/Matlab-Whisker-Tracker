function Objects = Costumbackground(Settings, Objects)
%%
%Insert costum processing script to extract an object mask (ie further
%processing or load an external file. Output must be a binary matrix with
%similar dimensions as standard frames, non-zero valued entries are
%consdiered an object
%%



% The gap crossing task is characterized by its gap, which is measured by
% taking intensity profiles along the y-axis at the borders. Any pixel
% within the gap is not an object.
ncols = 5;
y1 = sum(Objects(:,1:ncols),2);
y1 = y1./ncols;
y2 = sum(Objects(:,end-ncols+1:end),2);
y2 = y2./ncols;
y1(1:10) = 1;
y1(end-10:end) = 1;
y2(1:10) = 1;
y2(end-10:end) = 1;



idx = find(y1 < 0.5 |  y2 < 0.5);
Objects(idx,:) = 0;