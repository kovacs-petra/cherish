%%
clear all
k7=2.5; 
set(0,'DefaultAxesFontSize',10);
d=['015';'019';'025';'038';'050';'075';'100'];
FS=44100;
for di=1:7
    load(['imp_nkce_090_' d(di,:)]);
    implv(:,di)=impl*10^((-LEVEL-(di==7)*k7)/20);
    imprv(:,di)=impr*10^((-LEVEL-(di==7)*k7)/20);
    meandiffv(:,di)=meandiff; 
    LEVELv(:,di)=LEVEL;
end

implv = implv.*10.^((repmat(11.82,32767,7))/20);

for dummy = 1:0 % play the sounds
%%
soundsc([zeros(40000,2);reshape(implv,32767*7,1) reshape(imprv,32767*7,1);zeros(40000,2)],FS)

%%
soundsc([zeros(40000,2);reshape(imprv,32767*7,1) reshape(implv,32767*7,1);zeros(40000,2)],FS)
%v=[zeros(40000,2);reshape(imprv,32767*7,1) reshape(implv,32767*7,1);zeros(40000,2)];
%figure; plot(v)
end

%% find onsets -> tmpiv
figure
for di=1:7
    subplot(2,5,di); 
    %plot([implv(:,di) imprv(:,di)]);
    [tmp,tmpi]=min(imprv(:,di));
    [tmp,tmpli]=min(implv(1:400,di));
    R=tmpi+[-5:40];%[-50:200];
    tmpiv(di)=tmpi;
    plot([implv(R,di) imprv(R,di)]);
    axis([1 max(R)-min(R)+1  [-1 1]*5e10]);
    %axis([min(R) max(R)  [-1 1]*2e9]);
    legend(num2str(20*log10(rms([implv(R,di) imprv(R,di)]))))
    num2str(20*log10(rms([implv(R,di) imprv(R,di)])))
    title(num2str([tmpi tmpli tmpli-tmpi round((tmpli-tmpi)/.0441)]));
end

% create on-off ramps
% convolve stimuli

%%
dists = [15, 19, 25, 38, 50, 75, 100];
% old EdBr = [203.65 200.3 196.48 191.41 188.52 184.14 183.81-k]; % EdBl = [169.61 171.21 169.43 168.23 167.1 164.4 163.27-k];
EdBr = [203.2718 199.9065 196.0995 191.0197 188.1348 183.7738 183.4541-k7]; % R=tmpi+[-5:40];
EdBl = [157.7934 159.3903 157.6123 156.4116 155.2812 152.5913 151.4548-k7]; %R=tmpli+[-5:40];
figure; semilogx(dists-8,EdBr,'r'); hold on; semilogx(dists+8,EdBl,'b')

%%
D = [111 116 125 142 156 191 222 267 305];
D2 = (D-D(1))/44100*341;
(D2'+.15)*100
d
% the last two values should be very similar; and they are at 44100
%%
figure
for di=1:7
    subplot(2,4,di); 
    R=115:450;
    plot(20*log10(abs(fft([implv(R,di) imprv(R,di)],44100)))); 
    axis([0 22050 160 220]);    %axis([0 22050 120 200])
    legend(num2str(20*log10(rms([implv(R,di) imprv(R,di)]))),'location','south')
end



for dummy = 1:0 % analysis for 9 distances
%%
clear all
k7=2.5; k89=16;%K7 = 0; %2.5;
set(0,'DefaultAxesFontSize',10);
d=['015';'019';'025';'038';'050';'075';'100';'140';'170'];
FS=44100;
for di=1:9
    load(['imp_nkce_090_' d(di,:)]);
    implv(:,di)=impl*10^((-LEVEL-(di==7)*k7-(di>7)*k89)/20);
    imprv(:,di)=impr*10^((-LEVEL-(di==7)*k7)/20);
    meandiffv(:,di)=meandiff; 
    LEVELv(:,di)=LEVEL;
%     dir(:,di)=20*log10(abs(fft(impl(1:450),44100)));
%     rev(:,di)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));
%     dirr(:,di)=20*log10(abs(fft(impr(1:450),44100)));
%     revr(:,di)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));
end

implv = implv.*10.^((repmat(11.82,32767,9))/20);

implv(:,8:9)=-implv(:,8:9); imprv(:,8:9)=-imprv(:,8:9);
%%
soundsc([zeros(40000,2);reshape(implv,32767*9,1) reshape(imprv,32767*9,1);zeros(40000,2)],FS)

%%
soundsc([zeros(40000,2);reshape(imprv,32767*7,1) reshape(implv,32767*7,1);zeros(40000,2)],FS)
%v=[zeros(40000,2);reshape(imprv,32767*7,1) reshape(implv,32767*7,1);zeros(40000,2)];
%figure; plot(v)
%%
R=100:250;
figure
for di=1:9
    subplot(2,5,di); 
    %plot([implv(:,di) imprv(:,di)]);
    [tmp,tmpi]=min(imprv(:,di));
    [tmp,tmpli]=min(implv(1:400,di));
    R=tmpli+[-5:40];%[-50:200];
    plot([implv(R,di) imprv(R,di)]);
    axis([1 max(R)-min(R)+1  [-1 1]*5e8]);
    %axis([min(R) max(R)  [-1 1]*2e9]);
    legend(num2str(20*log10(rms([implv(R,di) imprv(R,di)]))))
    num2str(20*log10(rms([implv(R,di) imprv(R,di)])))
    title(num2str([tmpi tmpli tmpli-tmpi round((tmpli-tmpi)/.0441)]));
end

%%
dists = [15, 19, 25, 38, 50, 75, 100, 140, 170];
% old EdBr = [203.65 200.3 196.48 191.41 188.52 184.14 183.81-k]; % EdBl = [169.61 171.21 169.43 168.23 167.1 164.4 163.27-k];
EdBr = [203.2718 199.9065 196.0995 191.0197 188.1348 183.7738 183.4541-k7 176.0622 174.323]; % R=tmpi+[-5:40];
EdBl = [157.7934 159.3903 157.6123 156.4116 155.2812 152.5913 151.4548-k7 163.8182-k89 162.4626-k89]; %R=tmpli+[-5:40];
figure; semilogx(dists-8,EdBr,'r'); hold on; semilogx(dists+8,EdBl,'b')

%%
D = [111 116 125 142 156 191 222 267 305];
D2 = (D-D(1))/44100*341;
(D2'+.15)*100
d
% the last two values should be very similar; and they are at 44100
%%
figure
for di=1:7
    subplot(2,4,di); 
    R=115:450;
    plot(20*log10(abs(fft([implv(R,di) imprv(R,di)],44100)))); 
    axis([0 22050 160 220]);    %axis([0 22050 120 200])
    legend(num2str(20*log10(rms([implv(R,di) imprv(R,di)]))),'location','south')
end
end

%figure; 
% plot(20*log10(rms(imprv))-LEVELv); %axis([80 300 [-1 1]*1.5e8]);
% plot(20*log10(max(abs(imprv)))-LEVELv); %axis([80 300 [-1 1]*1.5e8]);
% plot(20*log10(max(abs(imprv)))-20*log10(max(abs(implv)))); %axis([80 300 [-1 1]*1.5e8]);
% plot([20*log10(rms(implv(80:300,:))); 20*log10(rms(imprv(80:300,:)))]'); %axis([80 300 [-1 1]*1.5e8]);
% %plot(0*20*log10(rms(imprv))-LEVELv); %axis([80 300 [-1 1]*1.5e8]);
% set(gca,'xticklabel',d);
%figure; plot(meandiffv); %plot(LEVELv);

% this line from apper_figs.m is probably relevant rvdb = 20*log10(rv);
% lvdb = 20*log10(lv); ilddb = lvdb-rvdb;ilddb(:,1:7,:,:,6)=ilddb(:,1:7,:,:,6)-11.82;
%lv(azi,dii,si,i,namei)
%names = ['jb';'ai';'mm';'ts';'gd';'nk'];azv=[0, 90];dsts=[15 19 25 38 50 75 100 140 170]';   


%%
if 0
load imp_nkce_000_015
dir(:,1)=20*log10(abs(fft(impl(1:450),44100)));
rev(:,1)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));

load imp_nkce_000_150
dir(:,2)=20*log10(abs(fft(impl(1:450),44100)));
rev(:,2)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));

load imp_nkce_090_015
dir(:,3)=20*log10(abs(fft(impl(1:450),44100)));
rev(:,3)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));

load imp_nkce_090_150
dir(:,4)=20*log10(abs(fft(impl(1:450),44100)));
rev(:,4)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));

load imp_nkce_000_015
dirr(:,1)=20*log10(abs(fft(impr(1:450),44100)));
revr(:,1)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));

load imp_nkce_000_150
dirr(:,2)=20*log10(abs(fft(impr(1:450),44100)));
revr(:,2)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));

load imp_nkce_090_015
dirr(:,3)=20*log10(abs(fft(impr(1:450),44100)));
revr(:,3)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));

load imp_nkce_090_150
dirr(:,4)=20*log10(abs(fft(impr(1:450),44100)));
revr(:,4)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));

load vals

% subplot(211);
% semilogx(dir-rev);
% hold on;
% semilogx(mean(dir-rev,2),'k','linewidth',2);
% axis([300 5000 -20 40])
% title('direct-reverberant energy in the left ear');
% legend('000 015','000 150','090 015','090 150','mean',2);
% 
% subplot(212);
% semilogx(dirr-revr);
% hold on;
% semilogx(mean(dirr-revr,2),'k','linewidth',2);
% axis([300 5000 -20 40])
% title('direct-reverberant energy in the right ear');

end % if 0

if 0
d=['015';'019';'025';'038';'050';'075';'100';'150'];
for di=1:8
    load(['imp_nkce_000_' d(di,:)]);
    dir(:,di)=20*log10(abs(fft(impl(1:450),44100)));
    rev(:,di)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));
    dirr(:,di)=20*log10(abs(fft(impr(1:450),44100)));
    revr(:,di)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));
end
save vals000 dir rev dirr revr

d=['015';'019';'025';'038';'050';'075';'100';'150'];
for di=1:8
    load(['imp_nkce_090_' d(di,:)]);
    dir(:,di)=20*log10(abs(fft(impl(1:450),44100)));
    rev(:,di)=smooth_3rd_oct(20*log10(abs(fft(impl(451:end),44100))));
    dirr(:,di)=20*log10(abs(fft(impr(1:450),44100)));
    revr(:,di)=smooth_3rd_oct(20*log10(abs(fft(impr(451:end),44100))));
end
save vals090 dir rev dirr revr
end % if 0

figure
load vals000
subplot(221);
semilogx(dir-rev);
hold on;
semilogx(mean(dir-rev,2),'k','linewidth',2);
axis([300 5000 -20 40])
title('D-R energy in the left ear 0°');

subplot(222);
semilogx(dirr-revr);
hold on;
semilogx(mean(dirr-revr,2),'k','linewidth',2);
axis([300 5000 -20 40])
title('D-R energy in the right ear 0°');

load vals090
subplot(223);
semilogx(dir-rev);
hold on;
semilogx(mean(dir-rev,2),'k','linewidth',2);
axis([300 5000 -20 40])
title('D-R energy in the left ear 90°');

subplot(224);
semilogx(dirr-revr);
hold on;
semilogx(mean(dirr-revr,2),'k','linewidth',2);
axis([300 5000 -20 40])
title('D-R energy in the right ear 90°');

legend(d(1,:),d(2,:),d(3,:),d(4,:),d(5,:),d(6,:),d(7,:),d(8,:),'mean');
% legend('000 015','000 150','090 015','090 150','mean',2);

figure;
load vals000
semilogx((dir(:,1)-rev(:,1))-(dir(:,8)-rev(:,8)),'b');
hold on;
semilogx((dirr(:,1)-revr(:,1))-(dirr(:,8)-revr(:,8)),'r');
axis([300 5000 0 30])

load vals090
semilogx((dir(:,1)-rev(:,1))-(dir(:,8)-rev(:,8)),'b:');
semilogx((dirr(:,1)-revr(:,1))-(dirr(:,8)-revr(:,8)),'r:');

title('D-R energy @ 15 cm - D-R energy @ 150 cm');
legend('left 0°','right 0°','left 90°','right 90°');