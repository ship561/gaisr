(ns consensus_seq
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.stats :as stats]
            [incanter.core :as math]
            [clojure.set :as set]
            [simplesvm :as ssvm]
            [clojure.java.shell :as shell])
  (:use [clojure.contrib.condition
         :only [raise handler-case *condition*
                print-stack-trace stack-trace-info]]

        [clojure.contrib.pprint
         :only (cl-format compile-format)]

        net.n01se.clojure-jna
        refold
        libsvm2weka))

(def possible_pairs {"AU" 1 "UA" 1 "GC" 1 "CG" 1 "GU" 1 "UG" 1} )

(defn jna-malloc
  "Create a 'C' level USB buffer of size SIZE.  Returns a pair (as a
   vector): [ptr-to-the-buffer the-buffer] Where ptr-to-the-buffer is
   a JNA/C ptr object and the-buffer is a java.nio.DirectByteBuffer
   object." [size]
  (let [buf (make-cbuf size)
        ptr (pointer buf)]
    [ptr buf]))

;;read stockholm file creates a map where the key=name val=sequence
(defn read-sto [f]
  (reduce (fn [m x] 
            (let [[name s] (str/split #"\s{2,}+" x)]
              (if-not (nil? s)
                (if (= name "#=GC SS_cons")
                  (assoc m :cons (conj (get m :cons) s))
                  (assoc m :seqs (conj (get m :seqs) s)))
                m)))
          {:seqs [] :cons []} (io/read-lines f)))

(defn write-svm [out-file m]
  (doseq [i m]
    (println i)))

;;change the stucture line to something that can be read by RNAfold
(defn change-parens [struct]
  (str/replace-re #"\<" "("
                  (str/replace-re #"\>" ")" struct)))

(defn energy-of-seq [profile]
  (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
  (jna-invoke Integer RNA/set_ribo_switch 1)
  (jna-invoke Void RNA/update_fold_params)
   (map (fn [inseq]
          (let [struct (profile :structure)
                [i st] (remove-gaps inseq struct)
                [ptr buf] (jna-malloc (inc (count i)))
                e (jna-invoke Float RNA/energy_of_structure i st)
                ]
            ;;(prn "i =" i)
            ;;(prn "st=" st)
            e
            ;; (println "E =" e)
            ;; (println "sequence  =" i)
            ;; (println "Structure =" (.getString ptr 0 false))
            ))
        (profile :seqs)))

(defn energy-of-aliseq [profile]
  (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
  (jna-invoke Integer RNA/set_ribo_switch 1)
  (jna-invoke Void RNA/update_alifold_params)
  (let [struct (profile :structure)
        inseqs (profile :seqs)
        [btr bbuf] (jna-malloc (count (first inseqs)))
        [etr ebuf] (jna-malloc (* 4 2))
        e (jna-invoke Float RNA/energy_of_alistruct (into-array inseqs) 
                               struct 
                               (count inseqs) 
                               etr)]
    ;; (println "e =" e)
    (println "e = " (+ (.getFloat etr 0) (.getFloat etr 4)) "=" (.getFloat etr 0) "+" (.getFloat etr 4))
    ;; (doseq [i inseqs] 
    ;;   (println i))
    ;; (println struct)
    (+ (.getFloat etr 0) (.getFloat etr 4))))

(defn sci [Ealign Eavg]
  (prn "Ea" Ealign "E" Eavg)
  (/ Ealign (stats/mean Eavg)))

(defn sum [m]
  (apply + (vals m)))

;;expected number of base pairs given frequencies of a base a position
;;i and j
(defn expected_qij [prob1 prob2]
  (for [k1 (keys prob1)
        k2 (keys prob2)]
    (* (get possible_pairs (str k1 k2) 0) (get prob1 k1) (get prob2 k2))))

(defn fraction-base-b [freq-map]
  (reduce (fn [m [k v]]
            ;;(prn k v (sum freq-map))
            (assoc m k (/ v (sum freq-map))))
          {} freq-map))

;;calculates the information for each base in a column and returns a
;;map where key=base value=information
(defn information_i [p fract-map i]
  (reduce (fn [m [b q]]
            (if-not (= b ".")
              (assoc m b (* 1 q (math/log2 (/ q (get p b 0.25)))))
              (assoc m b (* 1 q (math/log2 q)))))
          {} (second (nth fract-map i))))

;;returns only the information in columns that are base paired
;;according to the alignment
(defn information_only_bp [info bp-loc]
  (loop [y bp-loc
         v []]
    (if (seq y)
      (recur (rest y)
             (reduce (fn [v i] 
                       (if (= (first i) (ffirst y))
                         (conj v (second i))
                         v))
                     v info))
      v)))

;;calculates the information of all columns in an alignment of sequences
(defn information [profile]
  (let [fract-map (profile :fract)
        p (profile :background)
        len (count (first (profile :seqs)))]
    (for [i (range len)]
      [i (sum (information_i p fract-map i))])
    ))

;;calculates the fraction of base pairs at i and j where there is a
;;base pair and returns a map where key=[i j] value=faction
;;complementary bp
(defn fract_comp_ij [inseqs bp-loc]
  (let [freqs (partition 2 (interleave (range (count (first inseqs)))
                                       (apply map vector (map #(rest (str/split #"" %)) inseqs))))]
    (reduce (fn [m [i j]]
              (assoc m [i j]
                     (/ (apply + (map (fn [b1 b2]
                                        (get possible_pairs (str b1 b2) 0))
                                      (second (nth freqs i)) (second (nth freqs j))))
                        (count inseqs))))
            {}  bp-loc)))

(defn mutual_info_ij [q fract-freqs i j]
  (let [qij (get q [i j])
        Eqij (apply + (expected_qij (second (nth fract-freqs i)) (second (nth fract-freqs j))))
        ]
    (if (= qij 1)
      (* -1 (math/log2 Eqij))
      (+ (* qij (math/log2 (/ qij Eqij))) (* (- 1 qij) (math/log2 (/ (- 1 qij)(- 1 Eqij))))))))

(defn mutual_info [profile]
  (let [s (profile :seqs)
        bp-loc (profile :pairs)
        len (count (first s))
        all-loc (for [i (range len)
                      j (range len)]
                  [i j])
        q (fract_comp_ij s bp-loc) ;;fraction of bp ij where
        ;;there is complementarity
        fract-freqs (profile :fract)]
    (map (fn [[i j]]
           (mutual_info_ij q fract-freqs i j))
         bp-loc)))

(defn pairwise_identity [inseqs]
  (stats/mean (for [i (range (count inseqs))
                    j (range (inc i) (count inseqs))]
                (let [s1 (nth inseqs i)
                      s2 (nth inseqs j)
                      c (loop [k 0
                               c 0]
                          (if (< k (count s1))
                            (if (and (not= (.substring s1 k (inc k)) ".")
                                     (not= (.substring s2 k (inc k)) ".")
                                     (= (.substring s1 k (inc k)) (.substring s2 k (inc k))))
                              (recur (inc k)
                                     (inc c))
                              (recur (inc k)
                                     c))
                            c))]
                  (/ c (max (count (str/replace-re #"\." "" s1)) (count (str/replace-re #"\." "" s2))))))))

(defn ratio [b1 b2]
  (let [sum (+ b1 b2)]
    (if (zero? sum) 0 (double (/ b1 sum)))))

(defn base_comp_features [profile]
  (let [s (map #(str/replace-re #"\." "" %) (profile :seqs))
        base-comp (map (fn [m]
                         (fraction-base-b m))
                       (map frequencies
                            (map #(rest (str/split #"" %)) s)))]
    (map (fn [j i b]
          (str j " 1:"(count i) " 2:"(ratio (+ (get b "G" 0) (get b "C" 0))
                                                          (+ (get b "A" 0) (get b "U" 0)))
                   " 3:"(ratio (get b "A" 0) (get b "U" 0))
                   " 4:"(ratio (get b "G" 0) (get b "C" 0))))
         (iterate inc 1) s base-comp)))

(defn zscore [profile]
  (io/with-out-writer "/home/kitia/bin/libsvm-3.1/zscore-out1.txt"
    (doseq [i (base_comp_features profile)]
      (println i)))
  (libsvm2weka/txt2csv "/home/kitia/bin/libsvm-3.1/zscore-out1.txt" "/home/kitia/bin/libsvm-3.1/zscore-out1.csv")
  (let [mu (getpredicted (libsvm2weka/wekaNN "/home/kitia/mean1.csv" "/home/kitia/bin/libsvm-3.1/zscore-out1.csv"))
        sdev (getpredicted (libsvm2weka/wekaNN "/home/kitia/sd1.csv" "/home/kitia/bin/libsvm-3.1/zscore-out1.csv"))]
    (map (fn [x mean sdev]
           ;;(prn "x" x "mean" mean "sigma" sdev)
           (/ (- x mean) sdev))
         (energy-of-seq profile) mu sdev)))

(defn profile [m]
  (let [struct (change-parens (first (m :cons)))
        s (m :seqs)
        freqs (partition 2 (interleave (range (count (first s)))
                                       (map frequencies
                                            (apply map vector (map #(rest (str/split #"" %)) s)))))
        fract-freqs (sort-by key (reduce (fn [l [n freq-map]]
                                           (assoc l n (fraction-base-b freq-map)))
                                         {} freqs))
        q (fraction-base-b (frequencies (flatten (map #(rest (str/split #"" %)) s))))
        pairs (refold/make_pair_table struct)]
  {:seqs s
   :structure struct
   :fract fract-freqs
   :background q
   :pairs pairs}))

(defn main [f]
  (let [m (profile (read-sto f))]
    (prn "zscore" (stats/mean (zscore m)))
    (prn "sci" (sci (energy-of-aliseq m) (energy-of-seq m)))
    (prn "information" (stats/mean (information_only_bp (information m) (m :pairs))))
    (prn "mutual info" (stats/mean (mutual_info m)))
    (prn "pairwise identity" (pairwise_identity (m :seqs)))
    (prn "number of seqs" (count (m :seqs)))))


;; (let [s (find-gaps "GGUAUG.UAUUUC...AA..CCCCA..C..GAUA.AGCCCCGGAA..CU.UA...UU..........G..UGU......U....G.U........GAA.AUAG....AAC")
;;       [ptr buf] (jna-malloc (inc (count s)))]
;;   (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
;;   (println "e = " (jna-invoke Float RNA/fold s ptr))
;;   (println "struct = " (.getString ptr 0 false)))
