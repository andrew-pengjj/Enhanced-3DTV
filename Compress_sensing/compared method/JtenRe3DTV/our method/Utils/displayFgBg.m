function displayFgBg(Orig,Mask,Bg,Fg,Fscore)
figure,
subplot(221),imshow(Orig,[]),title('Original image');
subplot(222),imshow(Mask,[]),title('True mask');
subplot(223),imshow(Bg,[]),title('Rec. Bg');
subplot(224),imshow(Fg,[]),title(['Rec. Fg. FMeasure:',num2str(Fscore)]);
end