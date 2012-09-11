(ns robustness
 (:require [clojure.contrib.string :as str]
            [clojure.java.shell :as shell]
            [incanter.stats :as stats]
            [incanter.charts :as charts])
  (:use edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.snippets-math
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.core :only (view)]
        refold
        ))

(defn fold
  "Folds a sequence of RNA and returns only the target structure."
  
  [s]
  (->> ((shell/sh "RNAfold"
                  "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                  "--noPS"
                  :in s)
        :out)
       (str/split-lines)
       second
       (str/split #" ")
       first))

(defn mutant-neighbor
  "Takes a string s and finds the 3L 1-mer mutants. String can only contain
   letters A, C, G, U."
  
  [s]
  (let [s (.toUpperCase s)]
    (flatten
     (for [i (range (count s))]
       (map (fn [r]
              ;;makes new sequence with the substitution r
              (str (subs s 0 i) r (subs s (inc i) (count s))))
            (keys (dissoc {"A" 1 "G" 1 "U" 1 "C" 1} ;;3 other bases to sub
                          (subs s i (inc i)))))))))

(defn inverse-fold
  "Given a target structure, it will use RNAinverse to find n sequences which fold into an identical structure"
  
  [target n]
  (loop [c 0
         cand []]
    (if (< c n)
      (let [x (remove nil?
                      (flatten 
                       (map (fn [[s ensemble]]
                              (when-not (re-find #"d=" s) (re-find #"\w+" s)))
                            (->> ((shell/sh "RNAinverse"
                                            "-Fmp"
                                            (str "-R" (- n c))
                                            "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                            :in target)
                                  :out)
                                 (str/split-lines)
                                 (partition-all 2))
                            )))]
        (recur (count (distinct cand))
               (concat (distinct cand) x)))
      (take n (distinct cand)))))

(defn neutrality
  "takes a string s and  returns neutrality of each of the seq compared
   to each of the 3L 1 neighbors mutants"

  [s]
  (let [st (fold s)
        neighbors (mutant-neighbor s) ]
    (map (fn [neighbor]
           (/ (- (count s) (levenshtein st (fold neighbor)))
              (count s)))
         neighbors)))

(defn neutrality_rand
  "neutrality from finding a random seq with the same structure. if neutrality
   for native seq > rand seq then the native seq is more robust."
  
  ([target]
     (neutrality_rand target 100))
  
  ([target n]
     (let [cand (inverse-fold target n)]
       #_(prn distinct? cand)
       ;;(filter #(= target (fold %)) cand)
       (map #(stats/mean (neutrality %)) cand))))


(defn suboptimals
  "Finds the centroid structure of suboptimal structures and a vector
   representation using 0's and 1's using a RNAmutants or RNAsubopt. s
   is the RNA sequence (case doesn't matter as it will be all
   upper-cased) and n is the number of suboptimal structures to
   consider."
  
  [s n]
  (let [;;s "AACGAUCCCGGG"
        ;;n 10
        s (.toUpperCase s)
        RNAmutants 0 #_((shell/sh "./RNAmutants"
                              "-l" "./lib/"
                              "--mutation" "1"
                              "-n" (str n)
                              "--input-string" s
                              :dir "/home/kitia/Desktop/RNAmutants/")
                    :out)
        RNAsubopt ((shell/sh "RNAsubopt"
                             "-p" (str n)
                             :in s)
                   :out)
        out (->> RNAsubopt str/split-lines)
        structures (->> (if (some #(re-find #"\w+" %) out)
                          (drop-until #(re-find #"\> sampling \d+" out))
                          out)
                        (remove #(re-find #"[^\(\)\.]" %))
                        )
        struct->matrix (fn [st] (reduce (fn [m kv]
                                         (assoc m kv 1))
                                       {} (make_pair_table st)))
        Z->centroid (fn [matrix]
                      (sort-by key
                               (reduce (fn [m [[i j] p]]
                                         (assoc m i "(" j ")"))
                                       (into {}
                                             (map #(vector %1 %2) (range (count s)) (repeat ".")))
                                       (filter (fn [[[i j] p]]
                                                 (and (< i j)
                                                      (>= p 0.5)))
                                               matrix))))
        struct->vector (fn [st] (map #(if (= \. %) 0 1) (seq st)))
        partition-function (reduce  (fn [m [k v]]
                                      (assoc m k (/ v n)))
                                    {}  (apply merge-with + (map struct->matrix structures)))
        centroid (apply str (vals (Z->centroid partition-function)))]
    ;;(doseq [i structures] (prn i))
    ;;(prn (apply str (vals (Z->centroid partition-function))))
    [centroid (struct->vector centroid)]
    ))




