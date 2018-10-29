function [rates,fractions,btd]=Social_optimum(NetSettings, OpSettings,c_u,bs)
for u=1:NetSettings.users
    fractions(u)=OpSettings.w_i(u)./sum(OpSettings.w_i(bs(:)==bs(u)));
end
rates=fractions.*c_u;
btd=1./rates;
end

