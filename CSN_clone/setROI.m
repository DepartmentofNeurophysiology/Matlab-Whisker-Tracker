function ROI = setROI(Frame)
% Cut out a square ROI set by the user

flag = 0;

while flag == 0
%%
    figure(1)
    clf
    imshow(Frame)
    title('select 4 points, starting left top corner, clockwise')
    [x,y] = getpts();
    if length(x) == 4
        hold on
        scatter(x,y,'b','filled')



        for i = 1:3
            plot(x(i:i+1),y(i:i+1),'r')
        end
        plot([x(4) x(1)], [y(4) y(1)], 'r')
        
        resp = inputdlg('correctROI (y/n)?');
    
        if resp{1} == 'y'
            flag = 1;
        end
    else
        disp('select 4 points')


    end      
    
    
end

%%

%{
% assume its a rectangle
ymean1 = round(mean(y(1:2)));
Objects(1:ymean1,:) = 1;
ymean2 = round(mean(y(3:4)));
Objects(ymean2:end,:) = 1;

xmean1 = round(mean([x(1), x(4)]));
xmean2 = round(mean(x(2:3)));

Objects(:,1:xmean1) = 1;
Objects(:,xmean2:end) = 1;
%}

ROI = [x,y];
close(figure(1))
