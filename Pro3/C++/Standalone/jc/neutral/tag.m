function [t,n] = tag(i0, E)
  e = size(E,1);
  r = [ones(1,E(i0,5)), zeros(1, 1000)];
  t = [];
  n = [];
  for i=i0:e
      sp = E(i, 5) + 1;
      if E(i,6) == -1
          % extinction
          r(sp) = [];
      else
          r = [r(1:sp-1),0,r(sp:end)];
      end
      n(length(t)+1) = sum(r);
      t(length(t)+1) = E(i,1);
      if n(end) == 1
          break;
      end
  end
end
