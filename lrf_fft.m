clc;
clear;
close all;

Y1=importdata('E:\聚类\聚类代码\LGOD-master\dataset\data_06.csv');%%%算整个数据集的引力，密度差别大的影响太大，归一化后又不太对；明天要不再把距离归一化一下
Y1=Y1.data;
% label=importdata('wineGroup.mat');
% function [AUC,F1] = lrf_fft(Y1)
label=Y1(:,3);
A=Y1(:,1:2);

% ncluster= length(unique(label));
% ncluster=5;
d=size(A,2);%维度
ND=size(A,1);

% A=A.';
% A=mapminmax(A,0,1);
% A=A.';
% plot(A(:,1),A(:,2),'.','MarkerSize',5);
% hold on

dij=pdist(A);
dist= squareform(dij);

K=121;
lrf=zeros(ND,d);%lrf is the local resultant force
lrfnorm=zeros(ND,K);

NS=createns(A,'NSMethod','kdtree');
for i=1:ND
    LRF=zeros(ND,d);
    [idx,Dist]=knnsearch(NS,A(i,:),'k',K);
    for k=2:K
        deltaDist=A(idx(k),:)-A(i,:);
        if norm(deltaDist)~=0
            LRF(i,:)=LRF(i,:)+(deltaDist/norm(deltaDist));
        end
        m=1/sum(Dist(1:k));
        lrf(i,:)=LRF(i,:)*(1/m);
        lrfnorm(i,k)=norm(lrf(i,:));
    end
end
% quiver(A(:,1),A(:,2),lrf(:,1),lrf(:,2));
% for neighborNums=2:K
% %     neighborNums=round(ND*0.015);%总数的10%
%     % neighborNums=;
%     lrf=zeros(ND,d);%lrf is the local resultant force
%     for i=1:ND
%         [idx,Dist]=knnsearch(NS,A(i,:),'k',neighborNums);%对每个点i距离k近的点
%         for k=2:neighborNums   %
%             deltaDist=A(idx(k),:)-A(i,:);
%             if norm(deltaDist)~=0
%                 lrf(i,:)=lrf(i,:)+(deltaDist/norm(deltaDist));
%             end
%         end
%         m(i)=1/sum(Dist);
%         lrf(i,:)=lrf(i,:)*(1/m(i));
%         lrfnorm(i,neighborNums)=norm(lrf(i,:));
%     end
% end

%% LGOD
for i=1:ND
    for k=1:K-1
        delta(i,k)=abs(lrfnorm(i,k)-lrfnorm(i,k+1));
%         theta(i,k)=sum(delta(i,1:k))
    end
%     theta(i)=sum(delta(i,:));
%     noise(i)=sum(abs(fft(delta(i,:))));
end
num=numel(find(label==1));
for j=1:K-1
    for i=1:ND
        noise(i,j)=sum(abs(fft(delta(i,1:j))));
    end
    [noise_sorted,~]=sort(noise(:,j),'descend');
    index=zeros(ND,1);
%     threld(j)=noise_sorted(num);
    threld(j)=2.5*std(noise(:,j))+mean(noise(:,j));
    for i=1:ND
        if noise(i,j)>=threld(j)
           index(i)=1; 
        end
    end
%    y=index;t=label;
%    [tp,fp] = roc(t,y);
%    AUC(j)=auroc(tp,fp);
%     F1(j)=Rand_Jaccard_FM_AdjustedRand(t,y);
% end
% 
% plot(AUC(2:80))
% axis([0 80 0 1])
% end