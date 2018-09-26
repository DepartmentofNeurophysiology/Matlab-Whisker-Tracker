function gapinfo = detectGap(Objects)
main_dim = 2;
s = sum(Objects, main_dim);

idx = find(s ~= 0);
didx = diff(idx);
[~, id] = max(didx);

gapinfo.main_ax = 2;
gapinfo.edge_1 = idx(id);
gapinfo.edge_2 = idx(id+1);