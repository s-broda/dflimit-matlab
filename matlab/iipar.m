function P = iipar(s,c,cas)
%IIPAR  The Girsanov pieces (A2, A1, A0, det Sig, ...) of cases DU and DT under the tilt, shared by TR_LIMDIST and
%   DFX_CORE (formerly duplicated in both; DET3, ADJ3, ADJ3N are their own files); see tr_limdist and the closed forms free of cancellation in the Supplementary Material.
% Cases DU, DT, zero start, by Girsanov: under the tilt exp(-th^2/2 int W^2)
% of Brownian motion (th = sqrt(c^2 - 2s)), X = (W(1), int W[, int r W]) is
% Gaussian with covariance Sig(th); e^(sB) gains exp(X'QX/2), Q = diag(c, -2s Gm),
% and R = X'KX/2.  Returns, elementwise in s:
%   lcd = log(cosh th det Sig), D = det(I - Q Sig) (cancellation-free),
%   A2 = det P' = D/det Sig, A1, A0: det(K - nu P') = 0  <=>  A2 nu^2 - A1 nu + A0 = 0
%   (P' = Sig^-1 - Q; nu = 2 kap the weights of R under the tilt), logMB.
T=[1.0 -0.33333333333333333 0.13333333333333333 -0.053968253968253968 0.021869488536155203 -0.0088632355299021966 0.003592128036572481 -0.0014558343870513183 0.00059002744094558598 -0.00023912911424355248 9.6915379569294503e-5 -3.9278323883316834e-5 1.5918905069328965e-5 -6.4516892156554308e-6 2.6147711512907546e-6 -1.0597268320104654e-6 4.2949110782738059e-7 -1.7406618963571648e-7 7.0546369464009683e-8 -2.8591366623052539e-8 1.1587644432798852e-8 -4.6962953982309016e-9 1.9033368339312759e-9 -7.7139336353590623e-10 3.126339545892087e-10 -1.2670576930305401e-10 5.1351914080393677e-11 -2.0812146867700474e-11 8.4348454190943383e-12 -3.4185140868111558e-12 1.3854715742948469e-12 -5.6151047924146801e-13 2.2757162553728748e-13 -9.2231305851395317e-14 3.7379940310967389e-14 -1.5149519187148605e-14 6.1398688626168138e-15 -2.4883951222762789e-15 1.0085085566354096e-15 -4.0873312268690138e-16 1.6565329513786284e-16 -6.7136751750487117e-17 2.7209500613044588e-17 -1.1027595233722371e-17 4.4693160072374125e-18];
E=[1.0 -0.5 0.20833333333333333 -0.084722222222222222 0.034350198412698413 -0.013922233245149912 0.0056424968100315323 -0.0022868190951648293 0.00092681292737742184 -0.0003756231338525945 0.00015223432221797662 -6.1698246877700525e-5 2.5005357609459247e-5 -1.0134289721572026e-5 4.1072729198567004e-6 -1.664615015128028e-6 6.7464305456632013e-7 -2.7342253129896734e-7 1.1081397802278207e-7 -4.4911213669537136e-8 1.820182931124577e-8 -7.3769235610848765e-9 2.9897545073926764e-9 -1.2117018619561619e-9 4.9108426750009155e-10 -1.9902895700495881e-10 8.066339801136952e-11 -3.2691643852499818e-11 1.3249424201396147e-11 -5.3697893706595739e-12 2.1762936597810882e-12 -8.8201859824934e-13 3.5746867347671487e-13 -1.4487659644686842e-13 5.871627293628006e-14 -2.3796809091881837e-14 9.6444834564008509e-15 -3.9087619176859166e-15 1.5841615363041244e-15 -6.4203648775599251e-16 2.6020758752402585e-16 -1.0545816304260601e-16 4.2740583616893926e-17 -1.7322106086512011e-17 7.0203851674541611e-18 -2.8452549391662691e-18];
q=c^2-2*s; th=sqrt(q);
if c>0, del=-2*s./(th+c); else, del=th-c; end
sm=abs(th)<1; e2=exp(-2*th);
% g1 = tanh/th, g2 = (1-sech)/th^2, g3 = (th-tanh)/th^3, g4 = (th^2-2+2sech)/(2th^4),
% g5 = (th^3/3-th+tanh)/th^5, t1 = 1 - tanh, lch = log cosh
tn=(1-e2)./(1+e2); sc=2*exp(-th)./(1+e2);
g1=tn./th; g2=(1-sc)./th.^2; g3=(th-tn)./th.^3; g4=(th.^2-2+2*sc)./(2*th.^4); g5=(th.^3/3-th+tn)./th.^5;
t1=2*e2./(1+e2); lch=th+log((1+e2)/2);
if any(sm(:))
  Q=q(sm); h=@(cf) polyval(fliplr(cf),Q);
  g1(sm)=h(T); g2(sm)=-h(E(2:end)); g3(sm)=-h(T(2:end)); g4(sm)=h(E(3:end)); g5(sm)=h(T(3:end));
  t1(sm)=1-th(sm).*g1(sm); lch(sm)=log(cosh(th(sm)));
end
% N11 = 1 - c g1 = (th - c tanh)/th = (del + c t1)/th (cancellation-free); small th: direct
N11=(del+c*t1)./th; N11(sm)=1-c*g1(sm);
if cas(2)=='U'
  s11=g1; s12=g2; s22=g3;
  dS=s11.*s22-s12.^2;
  D=N11.*(1+2*s.*s22)+2*c*s.*s12.^2;           % det(I - Q Sig), Q = diag(c,-2s)
  A2=D./dS;
  A1=(s11-2*s12)./dS+2*s;                       % tr(adj K P'), K = [1 -1; -1 0]
  A0=-ones(size(s));
else
  % X = (W1, m0, m1); Gm = [4 -6; -6 12]
  S=cell(3); S{1,1}=g1; S{1,2}=g2; S{1,3}=g3; S{2,2}=g3; S{2,3}=g4; S{3,3}=g5;
  S{2,1}=S{1,2}; S{3,1}=S{1,3}; S{3,2}=S{2,3};
  Qd={c,0,0;0,-8*s,12*s;0,12*s,-24*s};          % Q = blkdiag(c, -2s Gm)
  N=cell(3); for i=1:3, for j=1:3, N{i,j}=(i==j)-(Qd{i,1}.*S{1,j}+Qd{i,2}.*S{2,j}+Qd{i,3}.*S{3,j}); end, end
  N{1,1}=N11;                                   % (row 1 of Q is (c,0,0))
  D=det3(N); dS=det3(S);
  A2=D./dS;
  Si=adj3(S); Pp=cell(3); for i=1:3, for j=1:3, Pp{i,j}=Si{i,j}./dS-Qd{i,j}; end, end
  % P'11 = (Sig^-1)11 - c: from P' = N Sig^-1 (row 1 of N is cancellation-free)
  Pp{1,1}=(N{1,1}.*Si{1,1}+N{1,2}.*Si{2,1}+N{1,3}.*Si{3,1})./dS;
  K=[1 2 -6; 2 -12 12; -6 12 0]; aK=adj3n(K);
  aP=adj3(Pp);
  A1=0; A0=0; for i=1:3, for j=1:3, A1=A1+K(i,j)*aP{j,i}; A0=A0+aK(i,j)*Pp{j,i}; end, end
end
P=struct('th',th,'del',del,'lcd',lch+log(dS),'D',D,'A2',A2,'A1',A1,'A0',A0,'dS',dS,'lch',lch);
P.logMB=-c/2-lch/2-log(D)/2;
end
