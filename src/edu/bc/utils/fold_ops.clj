(ns edu.bc.utils.fold-ops
  (:require [clojure.contrib.string :as str]
           [clojure.java.shell :as shell]
           )
  (:use refold))

(defn inverse-fold
  "Given a target structure, it will use RNAinverse to find n
   sequences which fold into a similar structure. If :perfect? is
   true, only returns sequences which fold into identical structures
   else returns the first n sequences. Returns a list of sequences."
  
  [target n & {:keys [perfect? ncore]
               :or {perfect? false ncore 2}}]
  (let [inv-fold (fn [target n perfect?]
                   (->> (map (fn [[s ensemble]]
                               (if perfect?
                                 (when-not (re-find #"d=" s) (re-find #"\w+" s)) ;perfect match
                                 (re-find #"\w+" s))) ;take all output
                             ;;calls the RNAinverse to generate inverse-fold seqs
                             (->> ((shell/sh "RNAinverse"
                                             "-Fmp"
                                             (str "-R" n)
                                             "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                             :in target)
                                   :out)
                                  str/split-lines
                                  (partition-all 2)))
                        flatten
                        (remove nil? ))) ;imperfect matches removed if
                                        ;they were nil
        ;;generate the proper number of distinct inverse-fold sequences
        inv-seq (loop [c 0
                       cand []]
                  (if (< c n)    
                    (recur (count cand) ;distinct candidate seqs
                           ;;add current list to newly generated ones
                           (->> (pmap (fn [_] (inv-fold target (min 10 (quot n ncore)) perfect?)) (range ncore))
                                (apply concat cand )
                                distinct))
                    (take n cand)))] 
    inv-seq))

(defn struct->matrix
  "creates array of bp locations. Array resembles a hash-map where the
  keys are base-pair locations and the value is 1 if
  present. Locations not present are not represented in the hash-map."
  
  [st]
  (reduce (fn [m kv] ;creates array of bp locations
            (assoc m kv 1))
          {} (make_pair_table st)))

(defn suboptimals
  "Finds the centroid structure of suboptimal structures and a vector
   representation using 0's and 1's using a RNAmutants or RNAsubopt. s
   is the RNA sequence (case doesn't matter as it will be all
   upper-cased) and n is the number of suboptimal structures to
   consider."
  
  [s n & {:keys [centroid-only]
          :or {centroid-only true}}]
  (let [;;s "AACGAUCCCGGG"
        ;;n 10
        s (.toUpperCase s)
        structures (do (declare fold)
                       (fold s :foldtype "RNAsubopt" :n n))
        Z->centroid (fn [matrix] ;converts the list of subopt
                                ;structures into a centroid
                      (sort-by key
                               (reduce (fn [m [[i j] p]]
                                         (assoc m i "(" j ")"))
                                       (into {}
                                             (map #(vector %1 %2) (range (count s)) (repeat ".")))
                                       ;;keeps bases that have over 50%
                                       ;;representation in suboptimal structures
                                       (filter (fn [[[i j] p]] 
                                                 (and (< i j)
                                                      (>= p 0.5)))
                                               matrix))))
        struct->vector (fn [st] (map #(if (= \. %) 0 1) (seq st))) ;change structure representation
                                                                  ;to a vector of 0's and 1's
        map-structures (map struct->matrix structures) ;sparse matrix
                                                       ;of each structure
        ]
    ;;(doseq [i structures] (prn i))
    ;;(prn (apply str (vals (Z->centroid partition-function))))
    (if centroid-only
      (let [partition-function (reduce  (fn [m [k v]]
                                          (assoc m k (/ v n)))
                                        {}  (apply merge-with + map-structures))
            centroid (apply str (vals (Z->centroid partition-function)))]
        [centroid (struct->vector centroid)]) ;returns centroid and vector representation
      [0 map-structures]) ;returns all suboptimal structures
    ))

(defn fold
  "Folds a sequence of RNA and returns only the target
   structure. Target structure can either be centroid or MFE."
  
  [s & {:keys [foldtype n]
        :or {foldtype "mfe" n 10000}}]
  (case foldtype
    "mfe"
    (->> ((shell/sh "RNAfold"
                    "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                    "--noPS"
                    :in s)
          :out)
         (str/split-lines)
         second
         (str/split #" ")
         first)
   
    "centroid"
    (first (suboptimals s n))
    
    "RNAmutants"
    0 #_(->> ((shell/sh "./RNAmutants"
                        "-l" "./lib/"
                        "--mutation" "1"
                        "-n" (str n)
                        "--input-string" s
                        :dir "/home/kitia/Desktop/RNAmutants/")
              :out)
             (drop-until #(re-find #"\> sampling \d+" ))
             (remove #(re-find #"[^\(\)\.]" %)))
    
    "RNAsubopt"
    (->> ((shell/sh "RNAsubopt"
                    "-p" (str n) ;samples according to
                                        ;Boltzmann distribution
                    :in s)
          :out)
         str/split-lines
         (remove #(re-find #"[^\(\)\.]" %)))
    ))
