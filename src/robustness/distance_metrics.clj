(ns robustness.distance-metrics
  (:require [clojure.contrib.string :as str]
            [clojure.set :as sets]
            )
  (:use edu.bc.bio.seq-utils
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.fold-ops                
        robustness.utils
        ))

(def ^:dynamic *globals*
  {:nsamples 1000
   :ncore 6
   :bpdist true})

(defn bpsomething
  "uses the <1-d/L> method mentioned in Borenstein miRNA
  paper. Distance measured from given structure"
  
  [s st & {:keys [bp]
           :or {bp true}}]
  (let [[s st cons-keys] (degap-conskeys s st)
        neighbors (mutant-neighbor s)
        len (count st)
        dist (fn [neighbor]
               (- 1 (/ (bpdist st (fold neighbor) :bpdist bp)
                       len)))]
    (mean (map dist neighbors))))

(defn bpsomething2
  "uses the <1-d/L> method mentioned in Borenstein miRNA
  paper. Distance measured from given structure"
  
  [df neighbor]
  (let [st (second df)
        bp (*globals* :bpdist)
        len (count st)
        dist (fn [neighbor]
               (- 1
                  (/ (bpdist st (fold neighbor) :bpdist bp)
                     len)))]
    (dist neighbor)))

(defn pccsomething
  "uses the pearsons correlation coefficient comparing the wt
  structure to the mutant structure to find the disruption to the
  structure. Returns a value bounded by [0 1]."

  [s st & {:keys [ncore]
           :or {ncore 1}}]
  (let [[s st cons-keys] (degap-conskeys s st)
        n 10000
        stvec (map #(if (= \. %) 0 1) (seq st))
        neighbors (mutant-neighbor s)
        dist (fn [neighbor]
               (let [mutst (second (suboptimals neighbor n))]
                 (if (not-every? zero? mutst)
                   (- 1 (* 0.5
                           (- 1 (pearson-correlation stvec mutst))))
                   0)))]
    (mean (pxmap dist ncore neighbors))))

(defn pccsomething2
  "uses the pearsons correlation coefficient comparing the wt
  structure to the mutant structure to find the disruption to the
  structure. Returns a value bounded by [0 1]."

  [df neighbor]
  (let [n (*globals* :nsamples)
        [s st] df
        stvec (map #(if (= \. %) 0 1) (seq st))
        dist (fn [neighbor]
               (let [mutst (second (suboptimals neighbor n))]
                 (if (not-every? zero? mutst)
                   (- 1 (* 0.5
                           (- 1 (pearson-correlation stvec mutst))))
                   0)))]
    (dist neighbor)))

(defn expected-subopt-overlapsomething
    "finds the neutrality of a wt sequence s by using the expected-subopt-overlap"

    [s _ & args]
    (let [expected-subopt-overlap
          (fn [s1 s2]
            (let [P (fn [s]
                      (fold s {:foldmethod :RNAfoldp})
                      (parse-dotps "dot.ps"))
                  Ps (P s1)
                  Pt (P s2)
                  pairs (set (sets/union (keys Ps) (keys Pt)))]
              (/ (reduce #(+ %1
                             (* (get Ps %2 0)
                                (get Pt %2 0)))
                         0 pairs)
                 (sum Ps))))
          s (str/replace-re #"[\.\-]" "" s)
          neighbors (mutant-neighbor s)]
      (mean (map #(expected-subopt-overlap s %) neighbors))))

(defn subopt-overlap-seq
  "Determine the percent overlap of each n suboptimal structure of a
  sequence s to the consensus structure (cons-keys). compares against
  n suboptimal structures.  returns a frequency map where k=dist and
  v=frequency.

  average is subopt overlap"

  [s cons-keys n]
  (let [[_ substruct] (suboptimals s n :centroid-only false)
        dist (fn [st1 st2]
               (/ (count (sets/intersection st1
                                            (set (keys st2))))
                  (count st1)))]
    ;;takes percent overlap and
    ;;reduces it to a freqmap to
    ;;save memeory
    (frequencies (map (fn [ks]
                        ;;percent overlap 
                        (dist cons-keys ks));dist(s,t)
                      substruct))))

(defn subopt-seq
  "Determine the percent overlap of each n suboptimal structure of a
  sequence s to the consensus structure (cons-keys). compares against
  n suboptimal structures.  returns a frequency map where k=dist and
  v=frequency.

  average is subopt overlap"

  [df mut]
  (let [n (*globals* :nsamples)
        cons-keys (third df)
        [_ substruct] (suboptimals mut n :centroid-only false)
        dist (fn [st1 st2]
               (/ (count (sets/intersection st1
                                            (set (keys st2))))
                  (count st1)))]
    ;;takes percent overlap and
    ;;reduces it to a freqmap to
    ;;save memeory
    (->> substruct
         (map #(dist cons-keys %))
         frequencies
         mean)))
