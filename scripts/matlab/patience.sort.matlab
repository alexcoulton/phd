function subsequence = PatienceLIS(sequence)
   % Use patience sort to generate a (not necessarily unique)
   %   longest increasing subsequence from a given sequence

   subInds = [];   % indices of subsequence elements
   for s = 1:numel(sequence)
      % put new element on first stack with top element > sequence(s)
      newStack = find(sequence(s) <= [sequence(subInds) inf], 1);
      subInds(newStack) = s;   % put current index on top of newStack
      if (newStack > 1)
         % point to the index currently to the left
         pred(s) = subInds(newStack - 1);
      end
   end
   % recover the subsequence indices
   % last element in subsequence is found from last subInds
   pathInds = subInds(end);
   while (pred(pathInds(1)) > 0)
      pathInds = [pred(pathInds(1)) pathInds];   % add predecessor index to list
   end
   subsequence = sequence(pathInds);   % recover subsequence from indices
end
