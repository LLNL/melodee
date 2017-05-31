function [mel_y_init, mel_ordering, mel_params] = modifiedModel_bclDefault_init(varargin)
   narginchk(0,1);
   if (nargin >= 1)
      mel_params = varargin{1};
   else
      mel_params = struct();
   end


   % define time
   t = 0;


   %define the initial conditions
   if (isfield(mel_params, 'Cm'))
      Cm = mel_params.Cm;
   else
      Cm = 1;
      mel_params.Cm = Cm;
   end
   m_init = 0.00155000000000000;
   if (isfield(mel_params, 'g_Na'))
      g_Na = mel_params.g_Na;
   else
      g_Na = 14.8380000000000;
      mel_params.g_Na = g_Na;
   end
   j_init = 0.722500000000000;
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
   h_init = 0.757300000000000;
   if (isfield(mel_params, 'n'))
      n = mel_params.n;
   else
      n = 1;
      mel_params.n = n;
   end
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
   V_init = -75;
   n_init = 0.325000000000000;
   if (isfield(mel_params, 'bcl'))
      bcl = mel_params.bcl;
   else
      bcl = 1000;
      mel_params.bcl = bcl;
   end
   h = h_init;
   potassium_channel_n = n_init;
   j = j_init;
   m = m_init;
   Vm_V_init = V_init;
   V = Vm_V_init;
   mel_y_init = zeros(5, 1);
   mel_y_init(1) = V;
   mel_y_init(2) = h;
   mel_y_init(3) = j;
   mel_y_init(4) = m;
   mel_y_init(5) = potassium_channel_n;
   mel_ordering = struct();
   mel_ordering.V = 1;
   mel_ordering.h = 2;
   mel_ordering.j = 3;
   mel_ordering.m = 4;
   mel_ordering.potassium_channel_n = 5;
end


