subplot(241);  imagesc(cconvfft2(s1{1}.A{1}, s1{1}.X{1}));
subplot(242);  imagesc(cconvfft2(s2{1}.A{1}, s2{1}.X{1}));
subplot(243);  imagesc(s1{1}.A{1});
subplot(244);  imagesc(s2{1}.A{1});
drawnow;