function mel_dydt = modifiedModel_bclDefault(mel_time,mel_diffvars,varargin)


   mel_dydt = zeros(5,1);
   % define time
   t = mel_time;
   % make copies of the differential vars
   V = mel_diffvars(1);
   h = mel_diffvars(2);
   j = mel_diffvars(3);
   m = mel_diffvars(4);
   potassium_channel_n = mel_diffvars(5);
   narginchk(2,3);
   if (nargin >= 3)
       mel_params = varargin{1};
   else
      mel_params = struct();
   end


   % define the differential update
   if (isfield(mel_params, 'Cm'))
      Cm = mel_params.Cm;
   else
      Cm = 1;
      mel_params.Cm = Cm;
   end
   E_R = -75;
   iStim = 0;
   if (isfield(mel_params, 'g_Na'))
      g_Na = mel_params.g_Na;
   else
      g_Na = 14.8380000000000;
      mel_params.g_Na = g_Na;
   end
   if (isfield(mel_params, 'duration'))
      duration = mel_params.duration;
   else
      duration = 1;
      mel_params.duration = duration;
   end
   if (isfield(mel_params, 'strength'))
      strength = mel_params.strength;
   else
      strength = 60;
      mel_params.strength = strength;
   end
   if (isfield(mel_params, 'n'))
      n = mel_params.n;
   else
      n = 1;
      mel_params.n = n;
   end
   g_K = 36;
   if (isfield(mel_params, 'g_L'))
      g_L = mel_params.g_L;
   else
      g_L = 0.300000000000000;
      mel_params.g_L = g_L;
   end
   if (isfield(mel_params, 'offset'))
      offset = mel_params.offset;
   else
      offset = 0;
      mel_params.offset = offset;
   end
   if (isfield(mel_params, 'bcl'))
      bcl = mel_params.bcl;
   else
      bcl = 1000;
      mel_params.bcl = bcl;
   end
   E_L = E_R + 10.613;
   mel_temp_009 = t < offset;
   E_Na = E_R + 115;
   E_K = E_R - 12;
   if (mel_temp_009)
      bcl_time = 1000*bcl - offset + t;
   else
      beat = floor((-offset + t)./bcl);
      mel_temp_008 = beat >= n;
      if (mel_temp_008)
         bcl_time_beat = n - 1;
      else
         bcl_time_beat = 0;
      end
      bcl_time = -bcl.*bcl_time_beat - offset + t;
   end
   bcl_time;
   mel_temp_006 = V < -40;
   beta_m = 0.1./(exp(V/5 + 7) + 1) + 0.1./(exp(V/200 - 1/4) + 1);
   mel_temp_004 = V < -40;
   beta_n = 0.125*exp(V/80 + 15/16);
   i_Na = m.^3.*g_Na.*h.*j.*(-E_Na + V);
   m_inf = (1 + 0.00184221158116513*exp(-0.110741971207087*V)).^(-2);
   mel_temp_002 = V < -40;
   j_inf = (15212.5932856544*exp(0.134589502018843*V) + 1).^(-2);
   alpha_n = (-0.01*V - 0.65)./(exp(-V/10 - 13/2) - 1);
   i_L = g_L.*(V - E_L);
   alpha_m = 1./(exp(-V/5 - 12) + 1);
   h_inf = (15212.5932856544*exp(0.134589502018843*V) + 1).^(-2);
   i_K = potassium_channel_n.^4.*g_K.*(V - E_K);
   mel_temp_010 = bcl_time_beat < duration;
   mel_temp_000 = V < -40;
   if (mel_temp_010)
      squareStim_iStim = strength;
   else
   end
   i_Natot = i_Na;
   tau_m = alpha_m.*beta_m;
   Iion_area = i_K;
   if (mel_temp_004)
      mel_temp_005 = (V + 37.78).*(-25428*exp(0.2444*V) - 6.948e-6*exp(-0.04391*V))./(50262745825.954*exp(0.311*V) + 1);
   else
      mel_temp_005 = 0;
   end
   if (mel_temp_006)
      mel_temp_007 = 0.02424*exp(-0.01052*V)./(1 + 0.00396086833990426*exp(-0.1378*V));
   else
      mel_temp_007 = 0.6*exp(0.057*V)./(1 + 0.0407622039783662*exp(-0.1*V));
   end
   n_diff = -potassium_channel_n.*beta_n + alpha_n.*(-potassium_channel_n + 1);
   if (mel_temp_002)
      mel_temp_003 = 2.7*exp(0.079*V) + 310000*exp(0.3485*V);
   else
      mel_temp_003 = 5.92307692307692 + 2.26708690145756*exp(-0.0900900900900901*V);
   end
   leakage_current_Iion_area = i_L;
   if (mel_temp_000)
      mel_temp_001 = 4.43126792958051e-7*exp(-0.147058823529412*V);
   else
      mel_temp_001 = 0;
   end
   alpha_h = mel_temp_001;
   m_diff = (-m + m_inf)./tau_m;
   beta_j = mel_temp_007;
   HH_Iion_area = Iion_area + leakage_current_Iion_area;
   beta_h = mel_temp_003;
   alpha_j = mel_temp_005;
   tau_j = 1./(alpha_j + beta_j);
   tau_h = 1./(alpha_h + beta_h);
   Iion = HH_Iion_area./Cm;
   h_diff = (-h + h_inf)./tau_h;
   j_diff = (-j + j_inf)./tau_j;
   Iion_1 = Iion + i_Natot;
   V_diff = -Iion_1 + squareStim_iStim;
   % stuff the differential update into an array
   mel_dydt(1) = V_diff;
   mel_dydt(2) = h_diff;
   mel_dydt(3) = j_diff;
   mel_dydt(4) = m_diff;
   mel_dydt(5) = n_diff;
end
