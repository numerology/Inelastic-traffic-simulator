function [rates,fractions,btd]=Static_Slicing(NetSettings, OpSettings,c_u,bs)
shareVec = OpSettings.s_o / sum(OpSettings.s_o);
for u=1:NetSettings.users
    o_in=OpSettings.ops_belongs(u);
    qd=([OpSettings.ops_belongs==o_in].*[bs==bs(u)]);
    fractions(u) = shareVec(o_in) / sum(qd);
end
rates=fractions.*c_u;
btd=1./rates;
end
