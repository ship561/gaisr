(ns robustness.utils
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.core.reducers :as r]
            [edu.bc.utils.snippets-utils :as snip-utils]
            [edu.bc.bio.sequtils.alignment-info :as info])
  (:use edu.bc.utils.fold-ops
        refold))

(defn equiv= [x y delta]
  (<= (- x delta) y (+ x delta)))

(defn mutant-neighbor
  "Takes a string s and finds the 3L 1-mer mutants. String can only
   contain letters A, C, G, U. with-names returns a vector with the
   mutant name as well"
  
  [s & {:keys [with-names]
        :or {with-names false}}]
  (let [s (.toUpperCase s)]
    (apply concat
           (for [i (range (count s))]
             (map (fn [r]
                    ;;makes new sequence with the base substitution r
                    (if with-names
                      ;;returns the mutation name and seq otherwise just the seq
                      [(str (subs s i (inc i)) i r)
                       (str (subs s 0 i) r (subs s (inc i) (count s)))]
                      (str (subs s 0 i) r (subs s (inc i) (count s)))))
                  (keys (dissoc {"A" 1 "G" 1 "U" 1 "C" 1} ;;3 other bases to sub
                                (subs s i (inc i)))))))))

(defn degap-conskeys
  "Produces a vector of a degapped seq and structure and corresponding
   cons-keys. Takes a seq (s) and structure (st)"

  [s st]
  (let [s (.toUpperCase s)
        [s st] (remove-gaps s st)
        cons-keys (set (keys (struct->matrix st)))]
    [s st cons-keys]))

(defn  build-mut-neighbors
  "Takes a state [wt st conskeys] and returns a new state containing
  the mutant neighbors"

  [[wt st conskeys]]
  [wt (mutant-neighbor wt) st conskeys])

(defn calculate-dist

  [[wt muts st conskeys] distfn]
  (let [state [wt st conskeys]]
    (r/fold 15
            (fn ([] [])
              ([l r] (concat l r)))
            (fn ([] [])
              ([V i]
                 (conj V (distfn state i))))
            (vec muts))))

(defn neut-data->csv
  "takes the output data from the main-subopt-overlap function which
  should be a clj data structure. The function appends the data in csv
  format to a specific csv"

  [data dist-metric trainset outcsv]
  (let [;outcsv
        ;"/home/peis/bin/gaisr/robustness/compare-neutrality/neutrality-distribution-trainset3.csv"
        ]
    (snip-utils/with-out-appender  outcsv
      (doseq [[filename vs] data
              v vs ;neutrality value
              :let [sto (fs/basename filename)
                    {nm :name type :type}
                    (-> (str (re-find #"RF\d+\-seed" sto) ".sto")
                        info/parse-sto-function )]]
        (->> [sto v trainset nm (apply str type) dist-metric]
             (str/join "," )
             println)))))
