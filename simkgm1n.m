function [ ss mpes   mpep lambda a1  phi m   ] = simkgm1n( x,y,sigma,c,n )

% This is the source code of the paper "X. Ma, Z.B. Liu. The kernel-based
% nonlinear multivariate grey model. Applied Mathematical modellig, xxx."

% To use this code please see the comments below:
% Input:
% x --  the input series, column vectors 
% y --  the output series, column vectors
% sigma -- the kernel parameter sigma of the Gaussian kernel
% c  --  the regularized parameter 
% n  --  number of samples used for modelling

% Output:
% ss -- all the simulated values, including the fitting and prediction
% values
% mpes -- MAPE for the simulation (fitting)
% mpep -- MAPE for the prediction
% a1 -- the optimized parameter a
% lambda -- the Lagrangian multipliers lambda in the paper
% phi -- the nonlinear function phi 
% m -- the block matrix in the linear system (the right-down side matrix)

% Example:
% Firstly set the input and output series values 
% Here is the first case in the case study 1 in the paper 

%-------------------You can just copy the following code 

% %step 1: set the sample data, notice that all should be column vectors
% >> x=[0.2772	0.2808	0.2708	0.2862	0.2809	0.2536	0.2737	0.2656	0.2809	0.271	0.2801	0.2831	0.2684	0.2807	0.2703	0.2792	0.2828	0.2529	0.28	0.2717 ]';
% >> y=[0.0425	0.048	0.039	0.0226	0.0575	0.0489	0.0483	0.0493	0.0421	0.0512	0.0478	0.0475	0.0458	0.0473	0.0468	0.0443	0.0472	0.0416	0.0459	0.045 ]';
% % step 2:set the parameters sigma , c, and remember to set the number for modelling
% % 
%>> sigma=0.3;
%>> c=5;
%>> n=15; % in my paper, I firstly choose the first 15 ones for modelling
% % Step 3: run it with only one line 
%[ ss mpes   mpep lambda a1  phi m   ] = simkdgm1n( x,y,sigma,c,n );

%------------------You can just copy the above code and directly run it 
% mpes =
% 
%     6.0602
% 
% 
% mpep =
% 
%     3.5672

% %  If you want to see the plot, just do it
%>> plot([ss y]);
%>> legend('Prediction','Raw data');
% Good luck!

x1 = cumsum(x);
y1 = cumsum(y);
zy = 0.5* ( y1(2:end) + y1(1:end-1) );

matrix = getgreykernel(x(1:n,:),y(1:n),sigma,c);

m = [0 ones(1,n-1); ones(n-1,1)  matrix];
Y = [0;y(2:n)];

% solve the linear system 
params = [0, ones(1,n-1); ones(n-1,1), matrix] \ [0;y(2:n)];

% compute all the system parameters
u=params(1);

lambda = params(2:end);

a = lambda'* zy(1:n-1);

l=length(y);
ss(1)=y(1);
temp=0;

% the first parameter a in the model 
a1 = (1-0.5*a)/(1+0.5*a);

% compute the predicted 1-AGO values 

for i =2:l
    for j=2:n       
         temp= temp + lambda(j-1)*Gaussian(x1(j,:),x1(i,:),sigma) ;
    end   
         ss(i) = temp / (1+0.5*a) + a1 * ss(i-1) + u / (1+0.5*a);
    phi(i)=temp;
         temp=0;
end

% get the restored values 
ss = [ss(1) ss(2:end)-ss(1:end-1)];

sim=ss(1:n);
pre=ss(n+1:end);

% Compute the MAPEs 

mpes = mean(abs(sim'-y(1:n))./y(1:n))*100
mpep = mean(abs(pre'-y(n+1:end))./y(n+1:end))*100

testmape = mean(abs(pre(1:end-5)'-y(n+1:end-5))./y(n+1:end-5))*100 ;

ss=ss';

end

%  construct the block matrix in the linear system

function [ matrix ] = getgreykernel(x,y,sigma,c)

n = length(y);
x1 = cumsum(x);
y1 = cumsum(y);
zy = 0.5* ( y1(2:end) + y1(1:end-1) );
matrix = [];
for i=1:n-1
     for j=1:n-1
         if i==j
           matrix(i,j)= - zy(i)*zy(j)+ Gaussian(x1(i+1,:),x1(j+1,:),sigma) + 1/c;
         else
           matrix(i,j)= - zy(i)*zy(j)+ Gaussian(x1(i+1,:),x1(j+1,:),sigma);
         end         
     end 
end

end


