% clc;
% clear;
% close all;
% % Y=importdata('E:\����\���ݼ�\dataset\data_06.csv');
% % Y=Y.data;
% Y=importdata('E:\����\���ݼ�\��ά���ݼ�\sizes5.mat');
% % label=Y(:,3);
% A=Y(:,1:2);
% 
% % Y=importdata('E:\����\���ݼ�\libras.mat');
% % label=importdata('E:\����\���ݼ�\librasgroup.mat');
% % A=Y;
% 
% % ncluster=length(unique(label));
% ncluster=4;

function [ cl, Center ] = better( A,ncluster )
Y1=A.';
Y1=mapminmax(Y1,0,1);
A=Y1.';
% subplot(1,2,2)
% plot(A(:,1),A(:,2),'.k','MarkerSize',5);
[ND,d]=size(A);%ά��
rowNums=ND;
colNums=d;

neighborNums=ceil(rowNums*0.1);
if rowNums>100 && rowNums<=1000
    neighborNums=10+ceil(rowNums*0.02);
end
if rowNums>1000
    neighborNums=ceil(rowNums*0.04);
end
if neighborNums>60
    neighborNums=40;
end
K=neighborNums

dij=pdist(A);
dist= squareform(dij);
distMatrix=dist;
 
sortDistMatrix=zeros(size(distMatrix));
sortDistMatrixIdx=zeros(size(distMatrix));
Y=zeros(rowNums,colNums);
mass=zeros(rowNums,1);
Ynorm=zeros(rowNums,1);
kkk=neighborNums*ones(rowNums,1);
X=A;

for m=1:rowNums
    [sortDistMatrix(m,:),sortDistMatrixIdx(m,:)]=sort(distMatrix(m,:));
    a=sortDistMatrix(m,:);
    idx=sortDistMatrixIdx(m,:);
    NormY=zeros(neighborNums,1);
    possibleY=zeros(neighborNums,colNums);
    for k=2:neighborNums
        deltaDist=X(idx(k),:)-X(m,:);
        if norm(deltaDist)~=0
            Y(m,:)=Y(m,:)+(deltaDist/norm(deltaDist));%����λʸ����
        end
    end
    deltaDist=X(idx(neighborNums+1),:)-X(m,:);
    if norm(deltaDist)~=0
        possibleY(1,:)=Y(m,:)+(deltaDist/norm(deltaDist));%�����һ������Ϊ��һ�������Լ�
    end
    NormY(1)=norm(possibleY(1,:));
    for k=2:neighborNums
        deltaDist=X(idx(neighborNums+k),:)-X(m,:);
        if norm(deltaDist)~=0
            possibleY(k,:)=possibleY(k-1,:)+(deltaDist/norm(deltaDist));
        end
        NormY(k)=norm(possibleY(k,:));
    end
    [minNorm,idx2]=min(NormY);%k��[neighbornNums,2*neighborNums]��ȡ�����˵�k��ʹ�õõ��ĺ�����С
    if minNorm<norm(Y(m,:));%���k��neighborNumsʱ��������С�������2-k�ľ���
        mass(m)=sum(a(2:neighborNums+idx2(1)));
        kkk(m)=neighborNums+idx2(1);
        Y(m,:)=possibleY(idx2(1),:);
    else%���򣬼��㵽neighborNums�ľ���
        mass(m)=sum(a(2:neighborNums));
    end
    Y(m,:)=Y(m,:)*mass(m);
    Ynorm(m)=norm(Y(m,:));
end
lrf=Y;
% hold on
% quiver(A(:,1),A(:,2),Y(:,1),Y(:,2));
K=round(mean(kkk));
% K=mode(kkk);
maxd=max(max(dist));
cosld=zeros(ND,K);
C=zeros(ND,1);
for i=1:ND
    idx=sortDistMatrixIdx(i,1:K+1);
    for j=1:K+1
        if(i==idx(j))
            continue;
        end
        a=lrf(idx(j),:);b=A(i,:)-A(idx(j),:);
        if(norm(a)==0)||(norm(b)==0)
            continue;
        end
        coss=(dot(a,b)/norm(a))/norm(b);
        cosld(i,j)=coss;
    end
    C(i)=sum(cosld(i,:));
end
[~,ordC]=sort(C,'descend');

ddelta(ordC(1))=-1;
nneigh(ordC(1))=0;
for i=2:ND
    ddelta(ordC(i))=maxd;
    for j=1:i-1
        if(dist(ordC(i),ordC(j))<ddelta(ordC(i)))
            ddelta(ordC(i))=dist(ordC(i),ordC(j));
            nneigh(ordC(i))=ordC(j);
        end
    end
end
ddelta(ordC(1))=max(ddelta(:));
gama=zeros(ND,1);
% ddelta=mapminmax(ddelta,0.1,1);


% scatter(A(find(C>0),1),A(find(C>0),2))
% subplot(1,2,1)
% tt=plot(C(:),ddelta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
% xlabel ('���϶�CE','FontName','����','FontSize',20)
% ylabel ('����\delta','FontName','����','FontSize',20)

for i=1:ND
    gama(i)=ddelta(i)*C(i);
end
[gama_sorted,ordgama]=sort(gama,'descend');
for i=1:ncluster
    cl(ordgama(i))=i;
    icl(i)=ordgama(i);
end

% cmap=colormap;
% for i=1:ncluster
%    ic=int8((i*64.)/(ncluster*1.));
%    hold on
%    plot(C(icl(i)),ddelta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
% end

% C=C.';
% C=mapminmax(C,0.1,1);
% C=C';
for i=1:ND
    gama(i)=ddelta(i)*C(i);
end
[~,ordgama]=sort(gama,'descend');
% ce=mapminmax(C,0.1,1);

cl=-1*ones(ND,1);
q=javaObject('java.util.LinkedList');
r=zeros(ncluster,1);
for i=1:ncluster
    cl(ordgama(i))=i;
    icl(i)=ordgama(i);
%     Dist=sort(dist(ordgama(i),:));%��ÿ�����������������ľ����������
    Dist=sortDistMatrix(ordgama(i),:);
    r(i)=Dist(K);%ѡ���K����ľ��� ��Ϊ�뾶
    rho(i)=numel(find(dist(ordgama(i),:)<=r(i)));
    CE(i)=sum(C(sortDistMatrixIdx(ordgama(i),1:K)));
end

Center=A(ordgama(1:ncluster),:);
[~,ord]=sort(r);
for i=1:ncluster-1%ÿ���ص��ܶȲ�һ�����뾶rҲ��һ�����Ӹ��ܶȽ���ɨ��
    per(i)=r(i)/r(i+1);
end
% per(1)=0.5;
per(ncluster)=1;
% per=min(0.5,per);
for i=1:ncluster
    q.add(ordgama(ord(i)));%�ȰѾ������ķ������
    while(size(q)>0)
        now=q.remove();
        D=sortDistMatrix(now,1:K);idx=sortDistMatrixIdx(now,1:K);
        rho1=numel(find(D<=r(ord(i))));%�ܶ�Ϊ �뾶�ڵ����е� 
        CEsum=sum(C(idx(1:rho1)));
        if(CEsum>0.1*CE(ord(i))) %�����ǰ�����  �ȴ��������ܶȵ�50%
            for k=2:rho1
                if(cl(idx(k))==-1)%ֻҪ��ǰ��뾶��Χ�ڵĵ�û�й���
                q.add(idx(k));%�ͷ������
                cl(idx(k))=ord(i);
                end

            end
        end
        for ii=1:K %�������е㣬���վ��뵱ǰ���˳�����û������С��ɨ��뾶
            if(cl(idx(ii))==-1)&&(D(ii)<=r(ord(i)))
                cl(idx(ii))=ord(i); %�͹���
            end
        end
    end
end
for i=1:ND
    if(cl(ordC(i))==-1)
        cl(ordC(i))=cl(nneigh(ordC(i)));
    end
end
% halo=cl;
% camp=colormap;
% Y1=zeros(ND,2);
% for i=1:ncluster
%     nn=0;
%     ic=int8((i*64.)/(ncluster*1.)); 
%     for j=1:ND
%         if(halo(j)==i)
%             nn=nn+1;
%             Y1(nn,1)=A(j,1);
%             Y1(nn,2)=A(j,2);
%         end
%     end
%     hold on
%     plot(Y1(1:nn,1),Y1(1:nn,2),'o','MarkerSize',3,'MarkerFaceColor',camp(ic,:),'MarkerEdgeColor',camp(ic,:));
% end
% subplot(1,2,2)
% hold on
% xlabel ('ά��1','FontName','����','FontSize',20);
% ylabel ('ά��2','FontName','����','FontSize',20);
% for i=1:ncluster
%    ic=int8((i*64.)/(ncluster*1.));
%    hold on
% %    plot(C(icl(i)),ddelta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%    plot(A(icl(i),1),A(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
% end

end