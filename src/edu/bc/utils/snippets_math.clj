(ns edu.bc.utils.snippets-math
  (:require [incanter.stats :as stats])
  (:use edu.bc.utils
        edu.bc.utils.probs-stats))


(defn cov

  [x y]
  (- (stats/mean (map * x y))
     (* (stats/mean x)
        (stats/mean y))))

(defn pearsonsCC
  "Finds the pearson correlation coefficient between 2 cols. If either
   sdev is 0 then the correlation is returned as 0"

  [x y]
  (let [;;x [18 25 57 45 26 64 37 40 24 33]
        ;;y [15000 29000 68000 52000 32000 80000 41000 45000 26000 33000]
        cov (cov x y)
        ;;standard deviation of population
        sd (fn [x] (Math/sqrt (- (stats/mean (map #(* % %) x)) 
                                (* (stats/mean x) 
                                   (stats/mean x)))))]
    (if (or (= 0 (sd x))
            (= 0 (sd y)))
      0
      (/ (cov x y) (sd x) (sd y)))))

(defn rand-gauss
  "Generate a normally distributed random number centered on mu with a standard
   deviation of sigma."
  
  [mu sigma]
  (let [r (fn [] (- (rand 2) 1))
        norm (fn [] (loop [x (r)
                          y (r)]
                     (let [s (+ (* x x) (* y y))]
                       (if (>= s 1)
                         (recur (r)
                                (r))
                         (* x (Math/sqrt (/ (* -2 (log s)) s)))))))]
    (+ mu (* sigma (norm)))))

(defn mean
  "Takes a frequency map where k=value and v=frequency. Returns a mean for the frequency map"
  
  [m]
  (let [[n d] (reduce (fn [v [val n]]
                        [(+ (first v) (* val n))
                         (+ (second v) n)])
                      [0 0] m)]
    (/ n d)))

(defn variance
  "Takes a frequency map where k=value and v=frequency. Returns variance for the frequency map."
  [m]
  (let [sumsq (mean (reduce (fn [m [val n]]
                              (assoc m (* val val) n))
                            {} m))
        mu (mean m)]
    (- sumsq (* mu mu))))

(defn sd
  "Takes a frequency map where k=value and v=frequency. Returns standard deviation for the frequency map."
  
  [m]
  (Math/sqrt (variance m)))

(defn freqn->list
  "Takes a frequency map and makes a list of it so that the values (x)
  are represented n-times in the list. List can then be used to find
  summary statistics."

  [m]
   (flatten (map (fn[[x n]] (repeat n x)) m)))

(defn median
  "Takes a frequency map where k=value and v=frequency. Returns median for the frequency map."
  
  [m]
  (let [;m {4 1 2 1 3 1 1 1}
        m (->> m freqn->list sort vec)
        c (count m)]
    (if (odd? c)
      (m (int (/ c 2)))
      (mean [[(m (int (dec (/ c 2)))) 1]
             [(m (int (/ c 2))) 1]])
    )))

(defn median-est
  "Takes a frequency map where k=value and v=frequency. Returns and estimate of the median for the frequency map."
  
  [m]
  (let [m (probs m) ;frequency to probability distribution
        m (into {} (map (fn [k cdf] ;change pdf into cdf
                          [k cdf])
                        (keys m) (reductions + (vals m))))
        ;;sample 1000 values from the probability distribution
        m (frequencies (for [n (repeatedly 1000 rand)] 
                         (ffirst (filter #(< n (second %)) m))))
        ]
    ;;find the median of 1000 values sampled from the pdf
    (stats/median (freqn->list m)))) 
