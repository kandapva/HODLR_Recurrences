clear variables;
clc;
choice = 1;
% Choice can be -1,0,1
% -1    -> Hmatrix
% 0     -> HODLR2D
% 1     -> HODLR in 2D
syms z x(n) y(n) v(n) w(n) a b L
assume(n>=0 & in(n, 'integer') & in(a, 'integer') & in(b, 'integer') & L>0)

% x(n) -> Edge sharing
% y(n) -> Complete recurrence
% a -> 2*\tau*Nmax
% b -> Nmax^2
% Explanation on how the recurrence relation is solved symbolically through
% z-transform:
% https://in.mathworks.com/help/symbolic/compute-z-transforms-and-inverse-z-transforms.html

% Well-separated interactions
syms Wpzt
W       = w(n) - a*4^n
Wzt     = ztrans(W,n,z);
Wzt     = subs(Wzt,ztrans(w(n),n,z),Wpzt);
Wpzt    = solve(Wzt,Wpzt);
Wn      = iztrans(Wpzt,z,n);
Wn      = simplify(subs(Wn,w(0), a))

% Vertex-sharing interactions
syms Vpzt
if choice <0
    V       = v(n+1) - v(n) - 15 * Wn
    Vzt     = ztrans(V,n,z);
    Vzt     = subs(Vzt,ztrans(v(n),n,z),Vpzt);
    Vpzt    = solve(Vzt,Vpzt);
    Vn      = iztrans(Vpzt,z,n);
    Vn      = simplify(subs(Vn,v(0), b))
else
    Vn = Wn
end

% Edge-sharing interactions
%E       = x(n+1) - 2 * x(n) - 12 * Wn - 2 * Vn
syms Epzt
if choice <1
    E       = x(n+1) - 2 * x(n) - 12 * Wn - 2 * Vn
    Ezt     = ztrans(E,n,z);
    Ezt     = subs(Ezt,ztrans(x(n),n,z),Epzt);
    Epzt    = solve(Ezt,Epzt);
    En      = iztrans(Epzt,z,n);
    En      = simplify(subs(En,x(0), b))
else
    En       = Wn
end



% Complete hierarchical strucutre
syms Tpzt
T       = y(n+1) - 4 * y(n) - 4 * Vn - 8 * En
Tzt     = ztrans(T,n,z);
Tzt     =  subs(Tzt,ztrans(y(n),n,z),Tpzt);
Tpzt    = solve(Tzt,Tpzt);
Tn      = iztrans(Tpzt,z,n);
Tn      = simplify(subs(Tn,y(0), b))
latex(Tn)

% a = 3;
% b = 1000;
% 
% Tans = subs(Tn)
% fplot(log10(Tans), [1,8])
% hold on
