%% part 1 : for randomly scattered peaks, see what the tadmap would look like.
HiCFolder = '/Users/AndrewTon/Documents/YalePhD/projects/chromatin/HiCData/mouseMaps/';
readsFilename = append(HiCFolder,'MouseChr12.mat');
colormap(redblue)
im = load(readsFilename).d;
ii=im(800:1000,800:1000);
jj=0*ii;
%jj : make random points equal to 1, rest are 0
jj(sub2ind(size(ii),fix(201*rand(1,40))+1,fix(201*rand(1,40))+1))=1;
jj=tril(jj+jj');
%kk : turn jj into a TADmap? 
kk = tril(cumsum(tril(cumsum(jj,1,'reverse')),2));
figure(1)
imagesc(kk)
mm=flipud(diff(diff(flipud(ii),[],1),[],2));
dmm = diff(flipud(ii),[],1);
ddmm = diff(dmm,[],2);
flipii = flipud(ii);
%recData is an actual recreation of the data
recData = flipud(cumsum([flipii(1,:) ; cumsum([dmm(:,1) ddmm],2)],1));

%tril of diffsVsData is processed and unprocessed, triu is raw data
diffsVsData = tril(recData) + triu(ii) - diag(diag(ii));
%% part 2 : display hi-C rate map (right) vs log of original data (left)
%todo: identify the peaks somehow. we haven't handled that yet. - above
%code generates some random peaks and demonstrates that we can recapitulate
%from there.
minInd = 50;
highInd = 100;
diffii = diff(flipud(ii),[],1);
diff2ii =  diff(diff(flipud(ii),[],1),[],2);
flipd2ii = flipud(diff2ii);

figure(2)
ax(1) = axes('Units', 'normalized', 'Position', [ .1 .1 .35 .8]);
ax(2) = axes('Units', 'normalized', 'Position', [ .5 .1 .35 .8]);
imagesc(log(ii(minInd:highInd, minInd:highInd)), 'Parent', ax(1))
imagesc(flipd2ii(minInd:highInd, minInd:highInd), 'Parent', ax(2))
set(ax, 'dataAspectRatio', [1 1 1])
linkaxes(ax)
caxis(ax(2),[-0.005,0.005])
%% part 3 : extract peaks from rate map, turn into tad map, compare to original data
%data = diff2ii(minInd:highInd, minInd:highInd);
%flipdata = flipud(data);
%orig_data = flipud(cumsum([ii(minInd,:) ; cumsum([diffii ,flipud(data)],1),2]));

valuedTads = diff2ii;
valuedTads(diff2ii<-0.005)=0;
valuedTads(diff2ii>0.005)=0;
valuedTads(abs(diff2ii)<0.003)=0;
valuedTads = tril(flipud(valuedTads));
tadmap = tril(cumsum(tril(cumsum(valuedTads,1,'reverse')),2));
figure(3)
ax(1) = axes('Units', 'normalized', 'Position', [ .1 .1 .35 .8]);
ax(2) = axes('Units', 'normalized', 'Position', [ .5 .1 .35 .8]);
imagesc(valuedTads,'Parent',ax(1))
imagesc(tril(tadmap)+triu(ii(2:end,2:end))-diag(diag(ii(2:end,2:end))),'Parent',ax(2))
linkaxes(ax)
%%
figure(4)
imagesc(flipud(diff2ii))
colorbar
caxis([-0.005,0.005])

figure(5) 

