(ns robustness
 (:require [clojure.contrib.string :as str]
           [clojure.java.shell :as shell]
           [clojure.contrib.io :as io]
           [incanter.stats :as stats]
           [incanter.charts :as charts]
           )
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
                                            (str "-R" n)
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
  "takes a string s and  returns neutrality of the seq  when compared
   to each of the 3L 1-mutant neighbors"

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

(defn psdc [dataset]
  (map (fn [wt mut]
         (map #(* (Math/sqrt (count wt))
                  (- 1 (pearsonsCC wt %))) mut))
       (first dataset) (second dataset)))

(defn generate-vectors
  "Generates suboptimal structures to find the centroid
   structure (structure and vector representation) of sequence s. n
   inverse folded sequences will be generated based on the centroid
   structure. Then the centroid structure for all 3L 1-mutant
   neighbors will be found for s and its n inverse folded sequences.
   Returns the structure of s then the mutuant in vector form."

  [s n]
  (let [s (.toUpperCase s)
        [stc stv] (suboptimals s 10000)
        wts  (inverse-fold stc n)
        muts (pxmap
              (fn [i] (doall (map #(second (suboptimals % 10000)) i)))
              10
              (map #(mutant-neighbor %) (cons s wts)))] 
    [(cons stv (map #(second (suboptimals % 10000)) wts)) muts]))


(defn makechart
  "creates a JSON datastructure for a chart from a map. The JSON
   output file in .js format can then be used to create a highcharts
   chart

   f is typically \"gaisr/robustness/highchart-test.js\""
  
  [f dataset title subtitle]
  (let [;;abc (map #(map (fn [x] (stats/mean x)) (partition-all 3 %))
   ;;barr)
        abc (map #(map (fn [x] (stats/mean x)) (partition-all 3 %)) dataset)
        wt (first abc)
        n (count wt)
        xy {:chart {
                    :renderTo "container",
                    :type "line",
                    :marginRight 130,
                    :marginBottom 50
                    :height 500}
            :title {:text title,
                    :x -20
                    }
            :subtitle {:text subtitle 
                       :x -20}
            :xAxis {:title {:text "position"} 
                    :categories (vec (range 1 (inc n)))}
            :yAxis {:title {:text "pSDC"} 
                    :min 0 :maxPadding 0.001
                    :plotLines [{:value 0
                                 :width 1
                                 :color "#0808080"}]}
            :legend {:layout "vertical"
                     :align "right"
                     :verticalAlign "top"
                     :x -10
                     :y 100
                     :borderWidth 0}
            :series (vec
                     (conj
                      (map (fn [i y]
                             {:name (str "neg" i)
                              :data (vec y)
                              :visible true})
                           (iterate inc 1) (rest abc))
                      {:name "wt"
                       :data (vec wt)}
                      ))}
        ]
    (clojure.contrib.io/with-out-writer f 
      (print "var foo=")
      (prn (clojure.contrib.json/json-str xy)))))


(defn combinedset []
  (let [avgpt (fn [dataset] (map #(map (fn [x] (stats/mean x)) (partition-all 3 %)) dataset))
        makemap (fn [names data] 
                  (into {}
                        (map (fn [nm x] 
                               [nm {:wt (first x) :con (rest x)}])
                             names data)))
        exp1 (map avgpt 
                  (map psdc 
                       (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con.clj"))))
        exp1 (makemap [:L10seq2 :L13 :L20 :L21] exp1)
        exp2 (map avgpt
                  (map psdc 
                       (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con2.clj"))))
        exp2 (makemap [:L10seq1 :L10seq2 :L13 :L20 :L21] exp2)                             
        exp3 (map avgpt (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con3.clj")))
        exp3 (makemap [:L10seq1 :L10seq2 :L13 :L20 :L21 :FMN] exp3)
        comboexp (merge-with (fn [a b]
                               (let [{curwt :wt curcon :con} a
                                     {newwt :wt newcon :con} b]
                                 (assoc {} :wt (map + curwt newwt) :con (concat curcon newcon))))
                             exp1 exp2 exp3)]
    (reduce (fn [m [k {wt :wt con :con}]]
              (let [n {:L10seq1 2 :L10seq2 3 :L13 3 :L20 3 :L21 3 :FMN 1}
                    w (map (fn [x] (/ x (n k))) wt)
                    c (map stats/mean con)]
                (assoc m k 
                       (assoc {} :wt w :wtavg (stats/mean w) :conavg c
                              :con con))))
            {} comboexp)))
