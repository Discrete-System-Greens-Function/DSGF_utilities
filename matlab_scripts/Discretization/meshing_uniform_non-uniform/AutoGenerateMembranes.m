%% Helping Sweetie Pie

clear; clc;


x=30;
y=20;
z=4;
nm=1e-9;
foil=[4,5,6,7,7,8,9,10];
gap=15; %You will need to explore this gap function a bit

n=0;
for i = 1:x
    for j = 1:y
        for k = 1:z
            n=n+1;
            Left_main(n,:)=nm*[i,j,k];
            
            
            
        end
    end
end

nk=0;
n=0;
for nii=1:length(foil)
    k=length(foil)-nii;
    nk=nk+1;
    for j=1:y
        ni=0;
        for i=1:foil(nii)
            ni=ni+1;
            n=n+1;
            Left_foil(n,:)=nm*[ni/2+max(Left_main(:,1))/nm,j,k/2+1/2];
            ni
        end
    end
end
Left=[Left_main; Left_foil];

Right=[-Left(:,1),Left(:,2),Left(:,3)]+nm*[x*2+gap,0,0];


figure()
plot3(Left(:,1),Left(:,2),Left(:,3),'.b')

figure()
hold on
plot3(Left(:,1),Left(:,2),Left(:,3),'.b')
plot3(Right(:,1),Right(:,2),Right(:,3),'.r')