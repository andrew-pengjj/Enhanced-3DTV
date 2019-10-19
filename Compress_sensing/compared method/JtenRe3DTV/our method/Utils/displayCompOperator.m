function displayCompOperator(D1,D2)

szD1 = size(D1);
szD2 = size(D2);
if ~isequal(szD1,szD2)
    error('The size is not equal!\n');
end

figure(1);
for i = 1 : szD1(3)
    s = sprintf('%s-th frm',num2str(i));
    subplot(121); surf(D1(:,:,i));colormap(hot);
    title(['D1:',s, ' frontal slice']);
    
    subplot(122); surf(D2(:,:,i));colormap(hot);
    title(['D2:',s, ' frontal slice']);
    pause(0.3);
end

figure(5);
for i = 1 : szD1(2)
    s = sprintf('%s-th frm',num2str(i));
    subplot(121); surf(squeeze(D1(:,i,:)));colormap(hot);
    title(['D1:',s, ' horizontal slice']);
     
    subplot(122); surf(squeeze(D2(:,i,:)));colormap(hot);
    title(['D2:',s, ' horizontal slice']);
    pause(0.3);
end

figure(3)
for i = 1 : szD1(1)
    s = sprintf('%s-th frm',num2str(i));
    subplot(121); surf(squeeze(D1(i,:,:))); colormap(hot);
    title(['D1:',s, ' lateral slice']);
    
    subplot(122); surf(squeeze(D2(i,:,:))); colormap(hot);
    title(['D2:',s, ' lateral slice']);
    pause(0.3);
end

close all;
end