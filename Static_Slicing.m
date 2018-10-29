function [rates,fractions,btd]=Static_Slicing(NetSettings, OpSettings,c_u,bs)
for u=1:NetSettings.users
    o_in=OpSettings.ops_belongs(u);
    qd=([OpSettings.ops_belongs==o_in].*[bs==bs(u)]);
    fractions(u)=OpSettings.s_o(o_in)/sum(qd);
end
rates=fractions.*c_u;
btd=1./rates;
end
