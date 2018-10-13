function bestori = orifit_lsq(BraggPk, a, b, c, alpha, beta, gamma, angle_adj, idpk, pxs, res, trueori, part_corr)

%% Inputs:
% BraggPk: A table from REDp indicating where the max intensity of a
% reflection is (frame number)
% a, b, c, alpha, beta, gamma: unit cell dimensions, in angstrom and deg
% angle_adj: experimentally set oscillation angle
% idpk: a list of reflections observed on each frame (from REDp)
% pxs: pixel size
% res: resolution cutoff, in angstrom

%% Calculation of crystal fixed coordinate system
alpha = alpha/180*pi;
beta = beta/180*pi;
gamma = gamma/180*pi;

unitcell = [a;b;c;alpha;beta;gamma];

cosastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma));
cosbstar = (cos(alpha)*cos(gamma)-cos(beta))/(sin(alpha)*sin(gamma));
coscstar = (cos(beta)*cos(alpha)-cos(gamma))/(sin(beta)*sin(alpha));
sinastar = sqrt(1-cosastar^2);
sinbstar = sqrt(1-cosbstar^2);
%ori_mat = [1/a 1/b*cos(gamma) 1/c*cos(beta); 0 1/b*sin(gamma) -1/c*sin(beta)*cosastar; 0 0 1/c*sin(beta)*sinastar];

V = a*b*c*sqrt(1-(cos(alpha))^2-(cos(beta))^2-(cos(gamma))^2+2*cos(alpha)*cos(beta)*cos(gamma));
astar = b*c*sin(alpha)/V;
bstar = a*c*sin(beta)/V;
cstar = a*b*sin(gamma)/V;

ori_mat = [astar bstar*cosastar cstar*cosbstar; 0 bstar*sinastar -cstar*sinbstar*cos(alpha); 0 0 cstar*sinbstar*sin(alpha)];

resolution = @(h, k, l) 1/sqrt(h^2*astar^2 + k^2*bstar^2 + l^2*cstar^2 +2*k*l*bstar*cstar*cosastar + 2*l*h*cstar*astar*cosbstar + 2*h*k*astar*bstar*coscstar);

radii_ewald=1/0.02508;

pxs = pxs;

BraggPk = BraggPk;
bpk=fscanf(BraggPk,'%f');
bpk=reshape(bpk,[18 length(bpk)/18]);
bpk=bpk';
fclose(BraggPk);

idpk = idpk;
pk=fscanf(idpk,'%f');
pk=reshape(pk,[9 length(pk)/9]);
pk=pk';
fn = pk(end,1);
reflections_3dposition=zeros(size(pk,1),7);
reflections_3dposition(:,1:4) = pk(:,1:4);
for i = 1:size(reflections_3dposition,1)
    reflections_3dposition(i,5:7)=(ori_mat'*pk(i,2:4)')';
end
fclose(idpk);

reflections_3dposition=horzcat(reflections_3dposition,zeros(size(reflections_3dposition,1),2));

% ----------------append alpha-RA and distance to bragg condition to the
% reflection_3dpos matrix
disp('Calculating deviation from Bragg peak for every reflection...')
for i=1:size(reflections_3dposition,1)
    for j=1:size(bpk,1)
        if reflections_3dposition(i,2:4)==bpk(j,2:4)
            reflections_3dposition(i,8)=bpk(j,10);
            reflections_3dposition(i,9)=bpk(j,8)-reflections_3dposition(i,1);
            break;
        end
    end
end

n=[1; 1; 1]; %using primary rotation axis from previous cycle
n2=[0; 0; 0];

bestori=zeros(fn,4);
thetaset=zeros(fn,1);
phiset=zeros(fn,1);

ag=ones(fn-1,2);
ag(:,1)=2:fn;
ag(:,2)=ag(:,2)*angle_adj;

indrefdiff=zeros(size(reflections_3dposition,1),5);
numb_run = 0;
failed_run = 0;

warning('off', 'all');
warning

if part_corr == 1
    disp('Partiality correction: on');
else
    disp('Partiality correction: off');
end

f = @(x,xdata)sqrt((x(1)-xdata(1,:)).^2+(x(2)-xdata(2,:)).^2+(x(3)-xdata(3,:)).^2);
std_hist = [];
startpoint = 0;
tic
while norm(n-n2)>0.00018 || numb_run <= 10
    numb_run = numb_run + 1;
    disp(['Number of run: ', num2str(numb_run)]);
    
    indrefdiff=zeros(size(reflections_3dposition,1),5);
    n=n2;
    ind = 1;
    reflectionsonframe=[0 0 0 0 0];
    expori=bestori;
    expori(:,1)=1:fn;
    ff=0;
    fn=fn;
    
    fname=strcat('square_distance_sum_12.07_bpkModified_leastsq.txt');
    fom=fopen(fname,'w');
    
    for i = 1:fn
        while reflections_3dposition(ind,1)==i
            if resolution(reflections_3dposition(ind, 2),reflections_3dposition(ind, 3), reflections_3dposition(ind, 4)) > res
                reflectionsonframe=vertcat(reflectionsonframe,reflections_3dposition(ind,5:9));
                indrefdiff(ind,1:4)=reflections_3dposition(ind,1:4);
            end
            ind=ind+1;
            if ind == size(reflections_3dposition,1)+1
                break;
            end
        end
        reflectionsonframe=reflectionsonframe(2:end,:);
        
        %---------------------------------------------end
        if part_corr == 1
            if n ~= [0; 0; 0]
                for k=1:size(reflectionsonframe,1)
                    if mod(failed_run,2) == 0
                        rotmat=rotationmat3D(angle_adj/180*pi*reflectionsonframe(k,5),n); %not sure if it is rotation clockwise or counterclockwise
                    else
                        rotmat=rotationmat3D(angle_adj/180*pi*reflectionsonframe(k,5),-n);
                    end
                    newposition=rotmat*reflectionsonframe(k,1:3)';
                    reflectionsonframe(k,1:3)=newposition';
                end
            end
        end
        
        ydata = ones(1, size(reflectionsonframe,1))*radii_ewald;
        xdata = reflectionsonframe(:,1:3);
        xdata = xdata';
        
        fomb=1000;
        prevfomb=10000;
        k=0;
        
        if size(reflectionsonframe,1) < 4
            bestori(i,:)=[i 0 0 0];
            disp(['Frame ', num2str(i), ' too few reflections.'])
            fprintf(fom,'%4d %4d\r\n',i,size(reflectionsonframe,1));
        else
            if ff==0 % I know it's the first calculated frame
                ff=ff+1;
                if startpoint
                    theta = theta;
                    phi = phi;
                else
                    %-----------------Simulated Annealing
                    theta=pi;
                    phi=pi;
                end
                T=pi/2;
                sphere_center=[radii_ewald*sin(theta)*sin(phi) radii_ewald*sin(theta)*cos(phi) radii_ewald*cos(theta)];
                figuremerit=0;
                weight=0;
                for j=1:size(reflectionsonframe,1)
                    figuremerit=figuremerit+(radii_ewald-norm(reflectionsonframe(j,1:3)-sphere_center))^2*norm(reflectionsonframe(j,1:3));
                    weight=weight+norm(reflectionsonframe(j,1:3));
                end
                figuremerit=figuremerit/weight;
                fomb_sa=figuremerit;
                while T>0
                    theta_cd=normrnd(theta,T);
                    phi_cd=normrnd(phi,T);
                    sphere_center=[radii_ewald*sin(theta_cd)*sin(phi_cd) radii_ewald*sin(theta_cd)*cos(phi_cd) radii_ewald*cos(theta_cd)];
                    figuremerit=0;
                    weight=0;
                    for j=1:size(reflectionsonframe,1)
                        figuremerit=figuremerit+(radii_ewald-norm(reflectionsonframe(j,1:3)-sphere_center))^2*norm(reflectionsonframe(j,1:3));
                        weight=weight+norm(reflectionsonframe(j,1:3));
                    end
                    figuremerit=figuremerit/weight;
                    if figuremerit<fomb_sa
                        theta=theta_cd;
                        phi=phi_cd;
                        fomb_sa=figuremerit;
                    end
                    T=T-0.01;
                end
                prevtheta=theta;
                prevphi=phi;
                
                while prevfomb-fomb > 0.000001 || k<1000
                    k=k+1;
                    theta=prevtheta+(rand()-0.5)*2*pi/180*angle_adj;
                    phi=prevphi+(rand()-0.5)*2*pi/180*angle_adj;
                    sphere_center=[radii_ewald*sin(theta)*sin(phi) radii_ewald*sin(theta)*cos(phi) radii_ewald*cos(theta)];
                    figuremerit=0;
                    weight=0;
                    avgdistpx=0;
                    for j=1:size(reflectionsonframe,1)
                        avgdistpx=avgdistpx+(radii_ewald-norm(reflectionsonframe(j,1:3)-sphere_center))^2;
                        figuremerit=figuremerit+(radii_ewald-norm(reflectionsonframe(j,1:3)-sphere_center))^2*norm(reflectionsonframe(j,1:3));
                        weight=weight+norm(reflectionsonframe(j,1:3));
                    end
                    figuremerit=figuremerit/weight;
                    avgdistpx=avgdistpx/size(reflectionsonframe,1);
                    avgdistpx=avgdistpx^0.5/pxs;
                    if figuremerit<fomb
                        prevfomb=fomb;
                        fomb=figuremerit;
                        bestori(i,:)=[i sphere_center(1) sphere_center(2) sphere_center(3)];
                        prevtheta=theta;
                        prevphi=phi;
                        bpx=avgdistpx;
                    end
                end
                sphere_center=[radii_ewald*sin(prevtheta)*sin(prevphi) radii_ewald*sin(prevtheta)*cos(prevphi) radii_ewald*cos(prevtheta)];
                
                for j=1:size(reflectionsonframe,1)
                    indrefdiff(ind-size(reflectionsonframe,1)+j-1,5)=abs(radii_ewald-norm(reflectionsonframe(j,1:3)-sphere_center))/pxs;
                end
                fprintf(fom,'%4d %4d %9.7f %4.2f\r\n',i,size(reflectionsonframe,1),fomb,bpx);
            else
                options = optimoptions('lsqcurvefit', 'TolFun', 1e-9, 'display', 'off');
                bestori(i,2:4) = lsqcurvefit(f, bestori(i-1,2:4), xdata, ydata,[],[], options);
                bestori(i,1) = i;
            end
        end
        if n ~= [0; 0; 0]
            if ff==1
                expori(i,:)=bestori(i,:);
            elseif ff>1 && (i-1)/21~=floor((i-1)/21)
                rotmat=rotationmat3D(0.1/180*pi,n);
                expori(i,2:4)=(rotmat*expori(i-1,2:4)')';
            elseif ff>1 && (i-1)/21==floor((i-1)/21)
                expori(i,2:4)=expori(i-1,2:4);
            end
        end
        %thetaset(i)=prevtheta;
        %phiset(i)=prevphi;
        reflectionsonframe=[0 0 0 0 0];
    end
    %figure,scatter3(bestori(:,2),bestori(:,3),bestori(:,4),'.');
    
    %---------------finding the normal vector n of the fitted plane bestori
    [n2 v p]=affine_fit(bestori(2:end,2:4));
    
    ag=zeros((fn-1),2);
    
    for i=1:(fn-1) %calculating angles between adjacent frames
        u=bestori(i,2:4);
        v=bestori(i+1,2:4);
        if nnz(u) == 0 || nnz(v) == 0
            disp(['Frame ', num2str(i), ' no orientation computed.'])
            ag(i,:)=[i+1 angle_adj];
        else
            ag(i,1)=i+1;
            ag(i,2)=atan2d(norm(cross(u,v)),dot(u,v));
        end
    end
    %figure,plot(ag(:,1),ag(:,2),'.')
    
    a=[0];
    for i = 1:size(ag,1)
        if ag(i,2) > angle_adj-0.1 && ag(i,2) < angle_adj+0.1
            a=horzcat(a,ag(i,2));
        end
    end
    a=a(2:end);
    mean1_ = mean(a);
    std1_ = std(a);
    disp(['Running Result: Mean of angle:',num2str(mean1_),'  STD of angle:',num2str(std1_)]);
    
    %-----------------------------------Record figure of merit
    merits=fscanf(fom,'%f');
    merits=reshape(merits,[3,length(merits)/3]);
    merits=merits';
    
    %% saving the result for the first fitting result w/o partiality correction
    if n == [0; 0; 0]
        a=[0];
        for i = 1:size(ag,1)
            if ag(i,2) > angle_adj-0.1 && ag(i,2) < angle_adj+0.1
                a=horzcat(a,ag(i,2));
            end
        end
        a=a(2:end);
        mean1 = mean(a);
        std1 = std(a);
        n_first = n2;
        disp(['FirstRun: Mean of angle:',num2str(mean1),'  STD of angle:',num2str(std1)]);
        std_hist(end+1) = std1;
    end
    
    if norm(n-n2) > 0.00018
        a=[0];
        for i = 1:size(ag,1)
            if ag(i,2) > angle_adj*0.5 && ag(i,2) < angle_adj*2
                a=horzcat(a,ag(i,2));
            end
        end
        a=a(2:end);
        mean2 = mean(a);
        std2 = std(a);
        std_hist(end+1) = std2;
        
        if std2 <= min(std_hist)
            nbest = n2;
            besttheta = prevtheta;
            bestphi = prevphi;
        end
    end
    
    if norm(n-n2) <= 0.00018 && numb_run > 10
        a=[0];
        for i = 1:size(ag,1)
            if ag(i,2) > angle_adj*0.5 && ag(i,2) < angle_adj*2
                a=horzcat(a,ag(i,2));
            end
        end
        a=a(2:end);
        mean2 = mean(a);
        std2 = std(a);
        std_hist(end+1) = std2;
        
        if std2 > min(std_hist)+0.001
            try
                n = nbest;
                theta = besttheta;
                phi = bestphi;
                startpoint = 1;
            catch
                n = n_first;
            end
            disp('Angle did not converge, rerunning the partiality correction process...');
            disp(['Mean of angle:',num2str(mean2),'  STD of angle:',num2str(std2)])
            failed_run = failed_run +1;
        end
    end
    
    if failed_run >= 10
        disp(['Failed too many times. Program terminated.'])
        disp(['Check parameters and rerun the program.'])
        break
    end
    
end
toc

if failed_run < 10 && trueori ~= -1
    %scatter3(bestori(:,2),bestori(:,3),bestori(:,4),'.');
    figure,histogram(ag(:,2));
    %% calculate deviations from true orientation
    delta=acos((cos(gamma)-cos(alpha)*cos(beta))/(sin(alpha)*sin(beta)));
    a = unitcell(1);
    b = unitcell(2);
    c = unitcell(3);
    M=[a*sin(beta) b*sin(alpha)*cos(delta) 0; 0 b*sin(alpha)*sin(delta) 0; a*cos(beta) b*cos(alpha) c];
    
    a=M*[1;0;0];
    b=M*[0;1;0];
    c=M*[0;0;1];
    Mv=dot(a,cross(b,c));
    a_=cross(b,c)/Mv;
    b_=cross(c,a)/Mv;
    c_=cross(a,b)/Mv;
    
    uvw=zeros(size(bestori,1),4);
    
    for i=1:size(bestori,1)
        f=i;
        hkl=bestori(f,2:4);
        uvw(i,1)=i;
        uvw(i,2:4)=hkl(1)*a_+hkl(2)*b_+hkl(3)*c_;
    end
    uvw = uvw(:, 2:4);
    s = uvw(1, :);
    s = abs(s);
    [A, I] = sort(s);
    uvw = uvw(:, I);
    s0 = uvw(1,:);
    
    absori=fscanf(trueori,'%f');
    absori=reshape(absori,[5 length(absori)/5]);
    absori=absori';
    
    absori=absori(:,3:5);
    
    s = absori(1, :);
    s = abs(s);
    [A, I] = sort(s);
    absori = absori(:, I);
    s1 = absori(1,:);
    
    for i = 1:3
        if s0(i)*s1(i) < 0
            uvw(:, i) = -uvw(:,i);
        end
    end
    
    errorori=zeros(size(bestori,1),3);
    errorang=zeros(size(bestori,1),1);
    for i=1:size(bestori,1)
        u=uvw(i,:);
        v=absori(i,:);
        errorori(i,:)=u-v;
        errorang(i)=atan2d(norm(cross(u,v)),dot(u,v));
    end
    
    devi=mean(errorang);
    disp(['Tilt angle: ', num2str(angle_adj), '; Resolution: ', num2str(res)])
    disp(['Mean deviation to true orientation in degrees: ', num2str(devi)])
    standd = std(errorang);
    disp(['standard deviation: ', num2str(standd)])
end
fclose('all');