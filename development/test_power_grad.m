d = 10;

X = rand(d)
rho = X*X'
V = rand(d)
V = V+V';

grad1 = power_grad(1.1,rho,V)
grad2 = (mpower(rho+1e-5*V,1.1)-mpower(rho,1.1))*1e5
grad1./grad2