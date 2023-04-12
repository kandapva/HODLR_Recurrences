clear variables;
clc;
choice = 0;
% Choice can be -1,0,1,2
% -1    -> Hmatrix
% 0     -> HODLR2D
% 1     -> HODLR with d'=1
% 2     -> HODLR with d'=2 (HODLR in 3D)

syms z f(n) e(n) y(n) v(n) w(n) a b
assume(n>=0 & in(n, 'integer') & in(a, 'integer') & in(b, 'integer'))

% x(n) -> Edge sharing
% y(n) -> Complete recurrence
% a -> 2*\tau*Nmax
% b -> Nmax^2
% Explanation on how the recurrence relation is solved symbolically through
% z-transform:
% https://in.mathworks.com/help/symbolic/compute-z-transforms-and-inverse-z-transforms.html

% Well-separated interactions
syms Wpzt
W       = w(n) - a*8^n
Wzt     = ztrans(W,n,z);
Wzt     = subs(Wzt,ztrans(w(n),n,z),Wpzt);
Wpzt    = solve(Wzt,Wpzt);
Wn      = iztrans(Wpzt,z,n);
Wn      = simplify(subs(Wn,w(0), a))

% Vertex-sharing interactions
syms Vpzt
if choice <0
    V       = v(n+1) - v(n) - 63 * Wn
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
    E       = e(n+1) - 2 * e(n) - 60 * Wn - 2 * Vn
    Ezt     = ztrans(E,n,z);
    Ezt     = subs(Ezt,ztrans(e(n),n,z),Epzt);
    Epzt    = solve(Ezt,Epzt);
    En      = iztrans(Epzt,z,n);
    En      = simplify(subs(En,e(0), b))
else
    En       = Wn
end

% FACE-sharing interactions
%E       = x(n+1) - 2 * x(n) - 12 * Wn - 2 * Vn
syms Fpzt
if choice <2
    F       = f(n+1) - 4 * f(n) - 8 * En - 4 * Vn - 48 * Wn
    Fzt     = ztrans(F,n,z);
    Fzt     = subs(Fzt,ztrans(f(n),n,z),Fpzt);
    Fpzt    = solve(Fzt,Fpzt);
    Fn      = iztrans(Fpzt,z,n);
    Fn      = simplify(subs(Fn,f(0), b))
else
    Fn       = Wn
end

% Complete hierarchical strucutre
syms Tpzt
T       = y(n+1) - 8 * y(n) - 8 * Vn - 16 * En - 32 * Fn
Tzt     = ztrans(T,n,z);
Tzt     =  subs(Tzt,ztrans(y(n),n,z),Tpzt);
Tpzt    = solve(Tzt,Tpzt);
Tn      = iztrans(Tpzt,z,n);
Tn      = simplify(subs(Tn,y(0), b))
latex(Tn)