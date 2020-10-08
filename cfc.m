a = 4.49;
C = 200;
m = 0.344*10^-24;
l=pi/a;
% CFC:
an = a/2 * [0,0,0; 1,1,0; 1,0,1; 0,1,1; -1,-1,0; -1,0,-1; 0,-1,-1; 1,-1,0; 1,0,-1; 0,1,-1;-1,1,0; -1,0,1; 0,-1,1];
nk = 1000;
k_GAMMA_X=[linspace(0,0,nk);linspace(0,2*pi/a,nk);linspace(0,0,nk)];
k_X_W=[linspace(0,pi/a,nk);linspace(2*pi/a,2*pi/a,nk);linspace(0,0,nk)];
k_W_K=[linspace(pi/a,3*pi/(2*a),nk);linspace(2*pi/a,3*pi/(2*a),nk);linspace(0,0,nk)];
k_K_GAMMA=[linspace(3*pi/(2*a),0,nk);linspace(3*pi/(2*a),0,nk);linspace(0,0,nk)];
k_GAMMA_L=[linspace(0,pi/a,nk);linspace(0,pi/a,nk);linspace(0,pi/a,nk)];

k=[k_GAMMA_X k_X_W k_W_K k_K_GAMMA k_GAMMA_L];

kk=[linspace(0,2*pi/a,nk), 2*pi/a+linspace(0,pi/a,nk),2*pi/a+pi/a+linspace(0,pi/(sqrt(2)*a),nk),2*pi/a+pi/a+pi/(sqrt(2)*a)+linspace(0,3*pi/(sqrt(2)*a),nk),2*pi/a+pi/a+pi/(sqrt(2)*a)+3*pi/(sqrt(2)*a)+linspace(0,sqrt(3)*pi/a,nk)];

w = NaN(3,size(k,2));
for i = 1:size(k,2)
k_aktuell = [k(1,i) k(2,i) k(3,i)]; 
w(:,i) = disp_solver(an, m, C, k_aktuell)/sqrt(C/m); 

end
figure(1)
plot(kk, w(1,:),'b-',kk, w(2,:),'b-',kk, w(3,:),'b-')
hold on

ymax = round(10*1.2*max(w(:)))/10;
plot([0,0],[0,ymax],'r-')
plot([2*pi/a,2*pi/a],[0,ymax],'r-')
plot([2*pi/a+pi/a,2*pi/a+pi/a],[0,ymax],'r-')
plot([2*pi/a+pi/a+pi/(sqrt(2)*a),2*pi/a+pi/a+pi/(sqrt(2)*a)],[0,ymax],'r-')
plot([2*pi/a+pi/a+pi/(sqrt(2)*a)+3*pi/(sqrt(2)*a),2*pi/a+pi/a+pi/(sqrt(2)*a)+3*pi/(sqrt(2)*a)],[0,ymax],'r-')
xlim([kk(1) kk(end)])
ylim([0 ymax])
hold off

function w = disp_solver(an, m, C, k)

ar =an(2:end,:);
arnorm = zeros(size(ar,1),1);
for n = 1:size(ar,1)
arnorm(n) = norm(ar(n,[1,2,3]));
ar(n,[1,2,3]) = ar(n,[1,2,3])./arnorm(n);
end

for n = 1:size(ar,1)

    delta_uxx = (ar(:,1).*ar(:,1)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uxy = (ar(:,1).*ar(:,2)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uxz = (ar(:,1).*ar(:,3)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);

    delta_uyx = (ar(:,2).*ar(:,1)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uyy = (ar(:,2).*ar(:,2)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uyz = (ar(:,2).*ar(:,3)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);

    delta_uzx = (ar(:,3).*ar(:,1)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uzy = (ar(:,3).*ar(:,2)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
    delta_uzz = (ar(:,3).*ar(:,3)).*(exp(1i*sum(an(2:end,:).*repmat(k,size(ar,1),1),2))-1);
end
M = -real([sum(delta_uxx,1),sum(delta_uxy,1),sum(delta_uxz,1);sum(delta_uyx,1),sum(delta_uyy,1),sum(delta_uyz,1);sum(delta_uzx,1),sum(delta_uzy,1),sum(delta_uzz,1)]);
w = real(sqrt(eigs(M)*C/m));
end