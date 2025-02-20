"""C-COMPASS


Acronyms / terms:

* profile: here, a set of protein amounts across different fractions or
  compartments.
* CA: class abundance (median of protein amounts of a given class)
* CC: class contribution (contribution of a compartment to a profile,
  CC ∈ [0, 1])
* fCC: filtered class contribution. False positive (according to some
  percentile value) CC values are set to 0 and renormalized.
  fCC ∈ [0, 1]
* DS: distance score
* RL: relocalization (difference between two class contributions, RL ∈ [-1, 1])
* RLS: relocalization score (sum of RL values across all compartments)
  RLS ∈ [0, 2] (no relocalization .. full relocalization)
* nCC: normalized class contribution (= CC * CA)
* TPA: total protein amount
* CPA: compartment protein amount (= CC * TPA)
* nCPA: normalized CPA (= nCC * TPA)
"""
