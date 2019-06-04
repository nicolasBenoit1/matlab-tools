function [xyz,gril] = geologicCrossSection(data,p1,p2,rgb,discretisation,window,base)
 
% Function creates geological cross-section from known interfaces of geological units
% Author: Nicolas Benoit (2019), Geological Survey of Canada, nicolas.benoit@canada.ca

% input: 
% data (table format) = One interface for each unit on a regular grid with following attributes:
%       XPT: x coordinates
%       YPT: y coordinates
%       ZPT: z coordinates
%       SEQNUM: stratigraphic sequence number, 1 (top) to n (bottom) unit
%       STRATUM: unit name, e.g., till
%       THICK: thickness of units
% p1: 1x2 xy coordinate for starting point of cross-section 
% p2: 1x2 xy coordinate for ending point of cross-section 
% code: rgb colour code (n units x 3) 
% discretisation: resolution along x-coordinate
% window: interface smoothing when > 1
% base: minimum thickness of basal unit


%%
x = data.XPT; 
y = data.YPT;
ux=unique(x);
dx=ux(2)-ux(1); 
uy=unique(y);
dy=uy(2)-uy(1);
if dx~=dy
    disp('non-regular grid input')
    return
end

x1=p1(1);y1=p1(2);x2=p2(1);y2=p2(2); 
[m,b] = FindLineEquation(x1,y1,x2,y2);
xl=(x1:discretisation:x2); 
if isempty(xl)
    xl=(x1:-1*discretisation:x2);
end

yl = m * xl + b; 
if length(xl)==1
    yl=(y1:discretisation:y2);
    xl=yl*0+xl;
end

data1=[];
dmax=sqrt((dx/2)^2+(dy/2)^2); 
for i=1:length(yl) 
    d = sqrt((x-xl(i)).^2+(y-yl(i)).^2);
    mind=min(d);
    if mind<=dmax 
        I = d==mind; 
        t=data(I,:); 
        t.XPT=repmat(xl(i),size(t,1),1); 
        t.YPT=repmat(yl(i),size(t,1),1);
        data1=[data1;t]; 
    end
end

x1 = data1.XPT;
y1 = data1.YPT;
z1 = data1.ZPT;
xy = [x1 y1];

d = [0;sqrt([diff(xy(:,1)).^2 + diff(xy(:,2)).^2])];

d=cumsum(d);
xyz = [x1 y1 z1 d d*0];
x1=d;

fn1 = data1.SEQNUM;
f1 = data1.STRATUM;

figure
hold on
u=sort(unique(fn1));
k=1;
pls=[];
dxx=200;dzz=1;
gril=grille2(min(x1),max(x1),dxx,ceil(min(z1)-base),max(z1),dzz);
gril=[gril gril(:,1)*0];

for i=1:length(u)
    id1 = fn1==u(i);
    if sum(data1.THICK(id1))>0
        p1=sortrows([x1(id1) z1(id1)],1);
        if i<length(u)
            id2 = fn1==u(i+1);
            p2=sortrows([x1(id2) z1(id2)],1);
        else
            p2=p1;
            p2(:,2)=min(p2(:,2))-base;
        end
        n=window;
        p1= [repmat(p1(1,:),n,1); p1]; p1(1:n,1)=(-1*n:1:-1)';
        p1= [p1;repmat(p1(end,:),n,1)]; p1(end-(n-1):end,1)=(max(p1(:,1))+1:1:max(p1(:,1))+1*n)';
        p2= [repmat(p2(1,:),n,1); p2]; p2(1:n,1)=(-1*n:1:-1)';
        p2= [p2;repmat(p2(end,:),n,1)]; p2(end-(n-1):end,1)=(max(p2(:,1))+1:1:max(p2(:,1))+1*n)';
        p1(:,2) = conv(p1(:,2),ones(n,1)/n,'same');
        p2(:,2) = conv(p2(:,2),ones(n,1)/n,'same');
        p1=p1(n+1:end-n,:);
        p2=p2(n+1:end-n,:);
        
        xyz(id1,3)=p1(:,2); 
        xyz(id1,5)=i;
        unique(p1-p2);
        p=[p1;flip(p2)];
        in = inpolygon(gril(:,1),gril(:,2),p(:,1),p(:,2));
        gril(in,3)=i;
        j=0;
        if j==1
            plot(gril(in,1),gril(in,2),'s','MarkerEdgeColor','none','MarkerFaceColor',rgb(i,:));hold on
        end
        temp=f1(id1);
        
        j=1;
        if j==1
            pl=patch(p(:,1),p(:,2),rgb(i,:),'EdgeColor','none'); hold on;
            legendInfo{k}=strcat(char(temp(1)));
            k=k+1;
            pls=[pls,pl];
        end
    end
end
xlim([min(p1(:,1)) max(p1(:,1))]);
ylim([min(z1)-base max(z1)]);
ylabel('Depth (m)')
xlabel('Distance (m)')

if j==1
    legend(pls,legendInfo,'location','northeastoutside'); 
end


function [m,b] = FindLineEquation(x1,y1,x2,y2)

m = (y2 - y1) / (x2 - x1); 
b = y1 - m*x1;
