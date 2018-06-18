% def CTR_f = solver(Dp)
function [cost, grad] = objective_F(Dp)
L = 1;
Nx = 11;
Nt = 1001;
T_end = 100;

x = linspace(0,L,Nx);

left_cond = 1;
right_cond = -1;
CTR_0 = -10;

CTR = zeros(Nx,Nt);
D = zeros(Nx,Nt);

D = reshape(Dp,[Nx,Nt]);

CTR(:,1) = CTR_0;
CTR(1,:) = left_cond;
CTR(Nx,:) = right_cond;

A = zeros(Nx,Nx,Nt);
b = zeros(Nx,Nt);
b1 = zeros(Nx,Nx);

dd = (L / (Nx - 1)^2) / (T_end / (Nt - 1));

for j = 1:(Nt-1)
    A(1,1,j) = 1;
    A(Nx,Nx,j) = 1;
    b(1,j) = CTR(1,j);
    b(Nx,j) = CTR(Nx,j);
    b1(1,1) = 1;
    b1(Nx,Nx) = 1;
    
    for i = 2:(Nx-1)
        A(i,i,j) = - (D(i,j) + D(i+1,j) + dd);
        A(i,i+1,j) = D(i+1,j);
        b(i,j) = - dd * CTR(i,j);
        b1(i,i) = - dd;
        
    end
    for i = 1:(Nx-2)
        A(i+1,i,j) = D(i,j);
    end

    CTR(:,j+1) = A(:,:,j) \ b(:,j);
end

CTR_f = reshape(CTR,[Nx*Nt,1]);
%CTR_f = CTR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_fix = zeros(Nx,Nt);
ssin = sin(x);
for j = 1:Nt
    C_fix(:,j) = ssin;
end

C_fix_f = reshape(C_fix,[Nx*Nt,1]);
lam = zeros(Nx,Nt);

lam(:,Nt) = - (CTR(:,Nt) - C_fix(:,Nt))' * inv(A(:,:,Nt-1));

for j = Nt-1:-1:2
    lam(:,j) = reshape(inv(A(:,:,j-1)), [Nx,Nx]) * reshape(- (CTR(:,j) - C_fix(:,j) + (-b1) * lam(:,j+1)), [Nx,1]);
end

gu = zeros(Nx,Nx,Nt);

for j = 1:Nt-1
    gu(1,1,j) = - CTR(1,j+1);
    gu(Nx,Nx,j) = -CTR(Nx,j+1);
    for i = 2:Nx-1
        gu(i,i,j) = CTR(i-1,j+1) - CTR(i,j+1);
        gu(i,i+1,j) = CTR(i+1,j+1) - CTR(i,j+1);
    end
end
grad = zeros(Nx,Nt);

for j = 1:Nt-1
    grad(:,j) = reshape(gu(:,:,j), [Nx,Nx]) * reshape(lam(:,j+1), [Nx,1]);
end

grad = reshape(grad,[Nx*Nt,1]);
cost = (CTR_f - C_fix_f)' * (CTR_f - C_fix_f) / 2;



end