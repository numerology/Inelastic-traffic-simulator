function [GenSpots] = makeSlawMap(size_max, n_wp, v_Hurst)
%
% makeSlawMap.m
%
% SLAW Waypoint Map Generator
% Written by Kyunghan Lee, KAIST, Korea
%
% function [GenSpots] = makeSlawMap(xmax, ymax, n_wp, v_Hurst)
%
% Input arguments
%	size_max: a side of a right square simulation area
%	n_wp: number of the waypoints
%	v_Hurst: Hurst parameter (0.5 < v_Hurst < 1)
%
% Based on the method of Kyunghan Lee (KAIST), Seongik Hong (NCSU),
%	Seong Joon Kim (NCSU), Injong Rhee (NCSU) and Song Chong (KAIST),
%	SLAW: A Mobility Model for Human Walks, The 28th IEEE
%	Conference on Computer Communications (INFOCOM), Rio de Janeiro,
%	Brazil, Apr. 2009.

%%
% slot levels, use 9 for higher resolution
Levels = 8;
% convert Hurst parameter
Valpha = 2 - 2*v_Hurst;
% initial variance  
v_initial_slot = 0.8;

Display = 1;        % Display Generated Spots

xmax = size_max;
ymax = size_max;

XY = [0 xmax 0 ymax];

%%
% ------------------------------------------------------------------------
% Map Generation
% ------------------------------------------------------------------------
for i=1:Levels
	xval(i) = 4^i;
	%Vslots_emul(i) = Vslots(1)*(xval(i)/xval(1))^Valpha;
	Vslots(i) = v_initial_slot*(xval(i)/xval(1))^Valpha;
end

Xmin = XY(1);
Xside = XY(2) - XY(1);
Ymin = XY(3);
Yside = XY(4) - XY(3);

%Gslots(1,1) = length(pausept);
Gslots(1,1) = n_wp;
Xoffset(1,1) = 0;
Yoffset(1,1) = 0;

for ia = 1:Levels
	disp_str = ['Level ' num2str(ia) ' of ' num2str(Levels) ' started'];
	disp(disp_str);

	for ib = 1:4^(ia-1)
		if mod(ib, 2000) == 0
			disp('.');
		end
		Xwind = Xside/(2^(ia-1));
		Ywind = Yside/(2^(ia-1));
		Xoffset(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Xoffset(ia, ib)+ [0 Xwind/2 0 Xwind/2];
		Yoffset(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Yoffset(ia, ib)+ [0 0 Ywind/2 Ywind/2];

		if Gslots(ia,ib) == 0
			Gslots(ia+1,4*(ib-1)+1:4*(ib-1)+4) = [0 0 0 0];
		else
			[ia ib];
			%Gslots(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Get4numsRand2(Gslots(ia,ib), Vslots(ia));

			if ia == 1
				Gslots(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Get4numsRand2(Gslots(ia,ib), Vslots(ia));
			else
				%Gslots(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Get4numsRand2(Gslots(ia,ib), ((Vslots(ia)+1)/(Vslots(ia-1)+1))-1 );
				Gtemp = Get4numsRand2(Gslots(ia,ib), ( (Vslots(ia)+1)/(var(Gslots(ia,1:4^(ia-1))/mean(Gslots(ia,1:4^(ia-1))),1)+1) )-1 );
				Gslots(ia+1,4*(ib-1)+1:4*(ib-1)+4) = Gtemp;
				if Gslots(ia,ib) ~= sum(Gtemp)
					Gtemp
					Gslots(ia,ib)
					( (Vslots(ia)+1)/(var(Gslots(ia,1:4^(ia-1))/mean(Gslots(ia,1:4^(ia-1))),1)+1) )-1
					disp('Error:Gtemp, please ask the authors.');
					return;
				end
			end
		end
	end
end

% ------------------------------------------------------------------------
% Map Plot
% ------------------------------------------------------------------------
if Display == 1

	for ie=1:Levels
		yval(ie) = var(Gslots(ie+1,1:4^(ie)) / mean(Gslots(ie+1,1:4^(ie))),1 );
		yval2(ie) = var(Gslots(ie+1,1:4^(ie)),1);
		yval3rd(ie) = mean( (Gslots(ie+1,1:4^(ie))/mean(Gslots(ie+1,1:4^(ie)))-1).^3 );
		xval(ie) = 4^ie; %Xside / (2^ie);
		xval_area(ie) = Xside*Yside / (4^ie);
	end

	Xmin = XY(1);
	Xside = XY(2) - XY(1);
	Ymin = XY(3);
	Yside = XY(4) - XY(3);

	Xwind = Xside / sqrt(4^(size(Gslots,1)-1));
	Ywind = Yside / sqrt(4^(size(Gslots,1)-1));

	ij= 0;
	for id = 1:4^(size(Gslots,1)-1)
		if Gslots(size(Gslots,1),id) == 0
		else
			for ie = 1:Gslots(size(Gslots,1),id)
				ij = ij + 1;
				theta = 2*pi*rand;
				LocX(ij) = Xoffset(size(Gslots,1),id) + Xwind/2 + (rand*Xwind/2)*cos(theta);
				LocY(ij) = Yoffset(size(Gslots,1),id) + Ywind/2 + (rand*Xwind/2)*sin(theta);
			end
		end
	end

% 	figure;
% 	plot(LocX+Xmin, LocY+Ymin, 'b.');
% 	xlabel('X (m)');
% 	ylabel('Y (m)');
% 	title('Distribution of Generated Spots');

	GenSpots(:,1) = LocX;
	GenSpots(:,2) = LocY;

end

%%
function [numbers] = Get4numsRand2(tot, vars)

% =======================
gran = 0.01;
Error = 0.03;
Thresh = 30;
% =======================

numbers = [0 0 0 0];
ii=0;
if vars >= 4
	aa = ceil(rand*4);
	numbers(aa) = tot;
elseif vars <= 0
	numbers = floor(tot/4)*ones(1,4);
	for ia = 1:4
		if sum(numbers) == tot
			return;
		end
		numbers(ia) = numbers(ia) + 1;
	end
else

	%=== Flatter Method ===================
	for ia = 0:gran:1
		for ib = (1-ia)/3:gran:1-ia
			for ic=(1-ia-ib)/2:gran:1-ia-ib
				%======================================
				id=(1-ia-ib-ic);
				zz=[ia ib ic id];

				if abs(vars - var(zz/mean(zz),1)) < Error*vars
					ii=ii+1;
					yy(ii,1:4) = zz;
					if ii > Thresh
						break;
					end
				end
			end
			if ii > Thresh
				break;
			end
		end
		if ii > Thresh
			break;
		end
	end

end

% ================================
% Randomize Solution
% ================================
yylen = size(yy,1);
xx = yy(ceil(yylen*rand),1:4);

% ================================
% Randomize Probability Slot (xx)
% ================================
temp = rand(1,4);
[tmp ord] = sort(temp, 'descend');
xx = xx(ord);

% ==========================================
% Distribute Total spots in Prability Slots
% ==========================================
xcum = [xx(1) sum(xx(1:2)) sum(xx(1:3)) sum(xx)];
for iz=1:tot
	dice = rand;
	for iy = 1:4
		if dice <= xcum(iy)
			numbers(iy) = numbers(iy) + 1;
			break;
		end
	end
end
%var(xx/mean(xx)) % NV verification
return;
