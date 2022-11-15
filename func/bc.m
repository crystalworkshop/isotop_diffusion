
function Eq=bc(X,par)
ct=par.alpha.*X+par.beta;
grad=-par.D.*(ct-X)';
Eq(1)=sum(X(1:5)) - par.csat;
Eq(2:5)=grad(2:5)*X(1)-X(2:5)'*grad(1);
KD_Hf=kdHf(par.T,par);
CHfs=X(6)*KD_Hf;
Cz=par.Cz*X(1)/par.csat;
Eq(6)=Cz*grad(6)-grad(1)*(CHfs-X(6));
end
