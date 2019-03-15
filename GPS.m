function [rates,fractions,btd]= GPS (NetSettings, OpSettings, c_u, bs)
nUsers = length(bs);
for u = 1 : nUsers
    ops=unique(OpSettings.ops_belongs([bs==bs(u)]));
    o_in=OpSettings.ops_belongs(u);
    qd=([OpSettings.ops_belongs==o_in].*[bs==bs(u)]);
    fractions(u)=[OpSettings.s_o(OpSettings.ops_belongs(u))/...
                                sum(OpSettings.s_o(ops))]/...
                                sum(qd);
end
rates=fractions.*c_u;
btd=1./rates;
end

