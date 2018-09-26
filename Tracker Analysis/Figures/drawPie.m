function drawPie(ax, center, theta_1, theta_2, radius, clr_in)

%%
T = Annotations.Janelia.Traces_clean{fnr};



base_width = 200;
tip_width = 30;
d = zeros(size(f)*10);

for i = 1:size(T,2)
    t = round(T{i});
    width_ax = round(linspace(base_width, tip_width, size(t,1)));
    for ii = 1:size(t,1)
        
        x1 = t(ii,1)*10 - width_ax(ii);
        if x1 < 1 
            x1 = 1;
        end
        
        x2 = t(ii,1)*10 +width_ax(ii);
        if x2 > size(d,2)
            x2 = size(d,2);
        end
        
        y1 = t(ii,2)*10 - width_ax(ii);
        if y1 < 1
            y1 = 1;
        end
        
        y2 = t(ii,2)*10 + width_ax(ii);
        if y2 > size(d,1)
            y2 = size(d,1);
        end
        
        d( x1:x2, y1:y2) =  d( x1:x2, y1:y2) + 1;
        
        
    end
   
    
end
figure(2)
clf
imagesc(d)
