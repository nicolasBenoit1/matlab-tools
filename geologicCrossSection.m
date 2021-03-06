function [xyz,gril] = geologicCrossSection(data,p1,p2,rgb,discretisation,window,base)

% Function creates geological cross-section from known interfaces of geological units
% Author: Nicolas Benoit (2019), Geological Survey of Canada, nicolas.benoit@canada.ca
% https://github.com/nicolasBenoit1/matlab-tools
% version 0.1

% input:
% data (table format) = One interface for each unit on a regular grid with following attributes:
%       XPT: x coordinates
%       YPT: y coordinates
%       ZPT: z coordinates
%       SEQNUM: stratigraphic sequence number, 1 (top) to n (bottom) unit
%       STRATUM: unit name, e.g., till
%       THICK: thickness of units
% p1: nx2 is the n xy coordinate for starting points of cross-section 
% p2: nx2 is the n xy coordinate for ending point of cross-section
%     ex.: cross-section containing 1 segment, n=1
%          cross-section containing 2 segments of different direction, n=2
% code: rgb colour code (n units x 3)
% discretisation: resolution along x-coordinate
% window: interface smoothing when > 1
% base: minimum thickness of basal unit below the lower elevation

nseg=size(p1,1);



%%
x = data.XPT;
y = data.YPT;
ux=unique(x);
dx=ux(2)-ux(1);
uy=unique(y);
dy=uy(2)-uy(1);
if dx~=dy
    xyz=[];gril=[];
    disp('non-regular grid input')
    return
end

p1a=p1;p2a=p2;
x0=0;
x0cum=0;
z1m1=inf;
z1m2=-inf;

figure

for n=1:nseg
    
    x1=p1a(n,1);y1=p1a(n,2);x2=p2a(n,1);y2=p2a(n,2);
    
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
                pl=patch(p(:,1)+x0,p(:,2),rgb(i,:),'EdgeColor','none'); hold on;
                legendInfo{k}=strcat(char(temp(1)));
                k=k+1;
                pls=[pls,pl];
            end
        end
    end
    z1m1=min([min(z1) z1m1]);
    z1m2=max([max(z1) z1m2]);
    x0=max(p1(:,1));
    x0cum=x0cum+x0;
end
xlim([0 x0cum]);
ylim([min(z1m1)-base max(z1m2)]);
ylabel('Depth (m)')
xlabel('Distance (m)')

if j==1
    warning off
    legend(pls,legendInfo,'location','northeastoutside');
end


function [m,b] = FindLineEquation(x1,y1,x2,y2)

m = (y2 - y1) / (x2 - x1);
b = y1 - m*x1;

function [gril,X,Y,nx,ny]=grille2(xmin,xmax,dx,ymin,ymax,dy);

x=[xmin:dx:xmax]';
y=[ymin:dy:ymax]';
nx=length(x);
ny=length(y);
gril=[kron(ones(ny,1),x), kron(y,ones(nx,1))];
X=reshape(gril(:,1),nx,ny)';
yy=sortrows([(length(gril):-1:1)' gril(:,2)],1);
Y=reshape(yy(:,2),nx,ny)';
