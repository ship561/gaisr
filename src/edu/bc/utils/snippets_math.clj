(ns edu.bc.utils.snippets-math
  (:use edu.bc.utils
        edu.bc.utils.probs-stats))

(defn rand-gauss
  "Generate a normally distributed random number centered on mu with a
   standard deviation of sigma."
  
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

(defn sum-map
  ""
  [m]
  (reduce (fn [cur-sum [next-val freq]]
            (+ cur-sum (* next-val freq)))
          0 m))

(defn median-est
  "Takes a frequency map where k=value and v=frequency. Returns and
   estimate of the median for the frequency map."
  
  [m]
  (let [m (probs m) ;frequency to probability distribution
        m (into {} (map (fn [k cdf] ;change pdf into cdf
                          [k cdf])
                        (keys m) (reductions + (vals m))))
        ;;sample 1000 values from the probability distribution
        m (frequencies (for [n (repeatedly 1000 rand)]
                         ;;takes number when rand falls into the
                         ;;cdf interval
                         (ffirst (filter #(< n (second %)) m))))
        ]
    ;;find the median of 1000 values sampled from the pdf
    (median m)))

(defn cov
  "Finds the covariation between a list x and a list y."
  
  [x y]
  (- (mean (frequencies (map * x y)))
     (* (mean (frequencies x))
        (mean (frequencies y)))))

(defn pearsonsCC
  "Finds the pearson correlation coefficient between 2 cols. If either
   sdev is 0 then the correlation is returned as 0"

  [x y]
  (let [;;x [18 25 57 45 26 64 37 40 24 33]
        ;;y [15000 29000 68000 52000 32000 80000 41000 45000 26000 33000]
        sdx (sd (frequencies x))
        sdy (sd (frequencies y))]
    (if (or (= 0 sdx)
            (= 0 sdy))
      0 ;return 0 when no standard deviation
      (/ (cov x y) sdx sdy))))
