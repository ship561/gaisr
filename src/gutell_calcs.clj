(ns gutell_calcs
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.stats :as stats]
            [incanter.core :as math]
            [clojure.set :as set])
  (:use [consensus_seq :only [read-sto profile]]))

(defn fractpairs_contain_nogaps
  "Checks pairs of columns x and y to see how many contain gaps. returns the percent of gaps in
   each pair of columns."
  
  [profile]
  (let [inseqs (profile :seqs)
        len (count (first inseqs))
        freqs (partition 2 (interleave (range (count (first inseqs)))
                                       (apply map vector (map #(rest (str/split #"" %)) inseqs))))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])]
    (reduce (fn [m [i j]]
              (let [pairs (map (fn [b1 b2]
                               (str b1 b2))
                             (second (nth freqs i)) (second (nth freqs j)))
                    fr (map #(re-find #"\.|-" %) pairs)] ;;fr = # of pairs containing gaps
                (assoc m [i j]
                       (double (/ (count (filter nil? fr)) (count fr)))
                       )))
            {}  all-loc)
    ))

(defn remove_gap_col
  "removes the columns from a caculation where there are greater than 'per' gaps.
   m = map, gap_col = column gap percents, per = threshold below which columns are removed"

  [m gap_col per]
  (select-keys m (keys (remove #(< (val %) per) gap_col)))
  )

(defn entropy [profile]
  (let [fract-map (profile :fract)
        len (count (first (profile :seqs)))
        ;;sum over all bases using -p*log2(p) at position i
        entropy_i (fn [p-baseb i] 
                    (math/sum
                     (map #(* -1 % (math/log %)) (vals p-baseb))))]
    (reduce (fn [m i]
              (assoc m i (entropy_i (second (nth fract-map i)) i)))
            {} (range len))
    ))

(defn joint_entropy [profile]
  (let [inseqs (profile :seqs)
        len (count (first inseqs))
        freqs (partition 2 (interleave (range (count (first inseqs)))
                                       (apply map vector (map #(rest (str/split #"" %)) inseqs))))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])]
    (reduce (fn [m [i j]]
              (let [fij (frequencies
                         (map (fn [b1 b2]
                                (str b1 b2))
                              (second (nth freqs i)) (second (nth freqs j))))
                    fr (map #(/ (second %) (math/sum (vals fij))) fij)] ;;fr = percentages
                (assoc m [i j]
                       (math/sum
                        (map #(* -1 % (math/log %)) fr))
                       )))
            {}  all-loc)
    ))

(defn mutual_info_ij [Hi Hij i j]
  (let [Hx (get Hi i)
        Hy (get Hi j)
        Hxy (get Hij [i j])]
    (+ Hx Hy (* -1 Hxy) ))
  )

(defn mutual_info [profile]
  (let [len (count (first (profile :seqs)))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])
        fract-freqs (profile :fract)
        Hx (entropy profile)
        Hxy (joint_entropy profile)]
    (reduce (fn [m [i j]]
              (assoc m [i j] (mutual_info_ij Hx Hxy i j)))
            {}  all-loc)))

(defn R [profile]
  (let [len (count (first (profile :seqs)))
        fract-map (profile :fract)
        Mxy (mutual_info profile)
        H (entropy profile)]
    (reduce (fn [m [i j]]
              (let [Hx (get H i)
                    Hy (get H j)
                    M (get Mxy [i j])]
                (assoc m [i j]
                       [Hx Hy M
                        (if (zero? Hx) 0 (/ M Hx))
                        (if (zero? Hy) 0 (/ M Hy))])))
            {} (for [i (range len)
                     j (range (+ 4 i) len)]
                 [i j]))))
