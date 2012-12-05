(ns snippets-test
  (:require [clojure.test :as test]))

(defn Z-old
  "Calculates the partition function for a structure (from sequences
  S) from i (starts at 0) to j. Currently assume homopolymer (any base
  can bind an y other base). Energy function = 1 to count structures
  but needs to be changed to actual energies.

  The probability of a base pair i,j occuring is given by:
  Pr(i,j) = (* Z(1,i-1) E(s1 s2) Z(i+1,j-1) Z(j+1,n))/Z(1,n)"

  [i j S]
  (let [S (.toUpperCase S)
        n (- j i -1)
        bp? (fn [b1 b2] ;;b1=base1 b2=base2
              (let [bp #{"AU" "UA" "GC"  "CG" "GU" "UG"}]
                (contains? bp (str b1 b2))))
        e (fn [i j S] (let [s1 (subs S i (inc i))
                           s2 (subs S j (inc j))
                           R 2
                           T 310
                           E (fn [b1 b2] (if (bp? b1 b2) 1 0))] ;;Energy of basepair, E(basepair)
                       (E s1 s2) #_(Math/exp (/ (E s1 s2) R T -1))))
        u 0 ;;min loop size
        ]
    (if (<= (- j i) u) 1
        (+ (Z i (dec j) S) ;;j unpaired
           (* (e i j S) (Z (inc i) (dec j) S)) ;;i,j pair
           (reduce (fn [x k]  ;;k,j paired for an intermediate a<k<b
                     (+ x (* (e k j S) (Z i (dec k) S) (Z (inc k) (dec j) S))))
                   0 (range (inc i) (- j u)))))))

(test/deftest simple-partition (test/is (= 6 (Z-old 0 4 "CGAGC"))))
