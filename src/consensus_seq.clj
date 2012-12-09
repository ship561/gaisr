(ns consensus_seq
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.stats :as stats]
            [incanter.core :as math]
            [clojure.set :as set]
            [clojure.java.shell :as shell])
  (:use [clojure.pprint :only (cl-format)]
        refold
        libsvm2weka
        edu.bc.bio.sequtils.files
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens profile)]
        edu.bc.utils
        [robustness :only (subopt-overlap-sto avg-overlap)]
        ))

(def possible_pairs {"AU" 1 "UA" 1 "GC" 1 "CG" 1 "GU" 1 "UG" 1} )
(def base #{"A" "C" "G" "U" "." "-"} )
(def all_pairs {"AA" 1 "AC" 1 "AG" 1 "AU" 1
                "CA" 1 "CC" 1 "CG" 1 "CU" 1
                "GA" 1 "GC" 1 "GG" 1 "GU" 1
                "UA" 1 "UC" 1 "UG" 1 "UU" 1} )

#_(defn jna-malloc
  "Create a 'C' level USB buffer of size SIZE.  Returns a pair (as a
   vector): [ptr-to-the-buffer the-buffer] Where ptr-to-the-buffer is
   a JNA/C ptr object and the-buffer is a java.nio.DirectByteBuffer
   object." [size]
  (let [buf (make-cbuf size)
        ptr (pointer buf)]
    [ptr buf]))


(defn write-svm [out-file m]
  (doseq [i m]
    (println i)))


#_(defn energy-of-seq [profile]
  (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
  (jna-invoke Integer RNA/set_ribo_switch 1)
  (jna-invoke Void RNA/update_fold_params)
   (map (fn [inseq]
          (let [struct (first (profile :structure))
                [i st] (remove-gaps inseq struct)
                [ptr buf] (jna-malloc (inc (count i)))
                e (jna-invoke Float RNA/fold i st)
                ;;e (jna-invoke Float RNA/energy_of_structure i st)
                ]
            ;;(prn "i =" i)
            ;;(prn "st=" st)
            e
            ;; (println "E =" e)
            ;; (println "sequence  =" i)
            ;; (println "Structure =" (.getString ptr 0 false))
            ))
        (profile :seqs)))

(defn energy-of-seq2 [profile]
  (map (fn [inseq]
         (let [struct (first (profile :structure))
               [i st] (remove-gaps inseq struct)
               foldout (-> (str i "\n" (str/replace-re #"\(|\)" "|" st))
                           ((fn [x] (shell/sh "echo" x)))
                           ((fn [x] (shell/sh "RNAfold" "-C" "--noPS" "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par" :in (get x :out))))
                           ((fn [x] (get x :out))))]
           (Double/parseDouble
                 (re-find #"\-*\d*.\d+" foldout))))
       (profile :seqs)))

#_(defn energy-of-aliseq [profile]
  (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
  (jna-invoke Integer RNA/set_ribo_switch 1)
  (jna-invoke Void RNA/update_alifold_params)
  (let [struct (first (profile :structure))
        inseqs (profile :seqs)
        [btr bbuf] (jna-malloc (count (first inseqs)))
        [etr ebuf] (jna-malloc (* 4 2))
        ali (jna-invoke Float RNA/alifold (into-array inseqs) btr)
        e (jna-invoke Float RNA/energy_of_alistruct (into-array inseqs) 
                               btr ;(.getString btr 0 false) ;struct 
                               (count inseqs) 
                               etr)]
    (println "e =" e ali)
    ;;(println "e = " (+ (.getFloat etr 0) (.getFloat etr 4)) "=" (.getFloat etr 0) "+" (.getFloat etr 4))
    ;; (doseq [i inseqs] 
    ;;   (println i))
    (println struct)
    (println (.getString btr 0 false))
    (prn (.getFloat etr 0) (.getFloat etr 4))))

(defn energy-of-aliseq2 [profile]
  (let [st (first (profile :structure))
        f (profile :filename)
        foldout (-> st
                    ((fn [x] (shell/sh "echo" x)))
                    ((fn [x] (shell/sh "RNAalifold"  "-P" "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par" "-r" "--noPS" "-C" f :in (get x :out))))
                    ((fn [x] (get x :out))))]
    (Double/parseDouble
     (re-find #"\-*\d*.\d+" foldout))
    ))

(defn sci [Ealign Eavg]
  ;;(prn "Ea" Ealign "E" Eavg)
  (/ Ealign (stats/mean Eavg)))




(defn expected_qij
  "expected number of base pairs given frequencies of a base a
  position i and j"

  [prob1 prob2]
  (for [k1 (keys prob1)
        k2 (keys prob2)]
    (* (get possible_pairs (str k1 k2) 0) (get prob1 k1) (get prob2 k2))))

;; (defn fraction-base-b [freq-map]
;;   (reduce (fn [m [k v]]
;;             ;;(prn k v (sum freq-map))
;;             (if (contains? base k)
;;               (assoc m k (/ v (sum freq-map)))
;;               m))
;;           {} freq-map))
  
(defn col->prob [col & {gaps :gaps :or {gaps false}}]
  (let [P (frequencies col)
        b (cond
           (char? (first (keys P)))
           (if gaps #{\A \C \G \U \.} #{\A \C \G \U})
           (= (count (first (keys P))) 1)
           (if gaps #{"A" "U" "G" "C" "."} #{"A" "U" "G" "C"})
           :else
           (if gaps
             (set (for [i [\A \C \G \U \.]
                        j [\A \C \G \U \.]]
                    (str i j)))
             (set (for [i [\A \C \G \U]
                        j [\A \C \G \U]]
                    (str i j)))))
        remove-gaps (fn [m]
                      (into {} (filter #(contains? b (key %)) m )))
        freq->prob (fn [f]
                     (reduce (fn [m [k v]]
                               (assoc m (str k) (/ (+ 0.0 v)
                                             (+ 0.0 (sum f)))))
                             {} f))]
    (freq->prob (remove-gaps (merge (into {} (map #(vector % 0) b))
                                    P)))))


;;calculates the information for each base in a column and returns a
;;map where key=base value=information
(defn information_i [p fract-map i]
  (reduce (fn [m [b q]]
            (if (contains? base b )
              (assoc m b (* 1 q (math/log2 (/ q (get p b 0.25)))))
              (assoc m b (* 1 q (math/log2 q)))))
          {} (second (nth fract-map i))))


(defn information_only_bp
  "returns only the information in column, i, that is base paired
   according to the alignment"
  
  [info bp-loc]
  (let [bp-loc (if (seq? bp-loc) (first bp-loc) bp-loc)]
    (map second (filter (fn [[i v]]
                          (contains? (set (map first bp-loc)) i))
                        info))))


(defn information
  "calculates the information in each column in an alignment of
   sequences"
  
  [profile]
  (let [fract-map (profile :fract)
        p (profile :background)
        len (count (first (profile :seqs)))]
    (reduce (fn [m i]
              (assoc m i (sum (information_i p fract-map i))))
            {} (range len))))


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


(defn fraction-base-b [freq-map]
  (reduce (fn [m [k v]]
            ;;(prn k v (sum freq-map))
            (if (contains? base k)
              (assoc m k (double (/ v (sum freq-map))))
              m))
          {} freq-map))
 
(defn mutual_info_only_bp
  "returns a vector of mutual information only at positions where
   there is base pairing provided by alignment."
  
  [info bp-loc]
  (let [bp-loc (if (seq? bp-loc) (first bp-loc) bp-loc)]
    (map second (filter (fn [[[i j] v]]
                          (contains? (set bp-loc) [i j]))
                        info))))

;;calculates the mutual information between 2 positions i and j
(defn mutual_info_ij [q fract-freqs i j]
  (let [qij (get q [i j])
        Eqij (apply + (expected_qij (second (nth fract-freqs i)) (second (nth fract-freqs j))))
        ]
    ;;(prn "i=" i "j=" j "qij" qij "Eqij=" Eqij)
    (cond
     (= qij 1)
     (* -1 (math/log2 Eqij))
     (= qij 0)
     (math/log2 (/ 1 (- 1 Eqij)))
     (= Eqij 0)
     0
     :else
     (+ (* qij (math/log2 (/ qij Eqij))) (* (- 1 qij) (math/log2 (/ (- 1 qij)(- 1 Eqij))))))))

;;mutual information for all possible combinations of i and j columns
;;in an alignment is calculated. A vector is returned consisting of
;;[[i j] mutual information at ij]. This measures the log likelihood
;;ratio of the observed base pairing occuring at random. higher = better
(defn mutual_info [profile]
  (let [s (profile :seqs)
        len (count (first s))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])
        q (fract_comp_ij s all-loc) ;;fraction of bp ij where
        ;;there is complementarity
        fract-freqs (profile :fract)]
    
    (map (fn [[i j]]
           [[i j] (mutual_info_ij q fract-freqs i j)])
         all-loc)))

(defn pairwise_identity [inseqs]
  (if (= 1 (count inseqs)) [[0 1]]
      (for [i (range (count inseqs))
            j (range (inc i) (count inseqs))]
        (let [s1 (nth inseqs i)
              s2 (nth inseqs j)
              c (loop [k 0
                       match 0
                       pairs 0]
                  (if (< k (count s1))
                    (if (or (and (not= (.substring s1 k (inc k)) "-")
                                 (not= (.substring s1 k (inc k)) "."))
                            (and (not= (.substring s2 k (inc k)) "-")
                                 (not= (.substring s2 k (inc k)) ".")))
                      (if (= (.substring s1 k (inc k)) (.substring s2 k (inc k)))
                        (recur (inc k)
                               (inc match)
                               (inc pairs))
                        (recur (inc k)
                               match
                               (inc pairs)))
                      (recur (inc k)
                             match
                             pairs))
                    [match pairs]))]
         c))))

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
  ;; (io/with-out-writer "/home/kitia/bin/gaisr/testset.csv"
  ;;   (doseq [i (base_comp_features profile)]
  ;;     (println i)))
  (txt2csv (base_comp_features profile) "/home/peis/bin/gaisr/testset.csv")
  (let [mu (getpredicted (wekalibsvm "/home/peis/bin/gaisr/mean.csv" "/home/peis/bin/gaisr/testset.csv" "mean"))
        sdev (getpredicted (wekalibsvm "/home/peis/bin/gaisr/sd.csv" "/home/peis/bin/gaisr/testset.csv" "sd"))]
    (map (fn [x mean sdev]
           ;;(prn "x" x "mean" mean "sigma" sdev)
           (/ (- x mean) sdev))
         (energy-of-seq2 profile) mu sdev)))



(defn main-file [f]
  (let [m (profile (read-sto f))]
    (prn "zscore" (stats/mean (zscore m)))
    (prn "sci" (sci (energy-of-aliseq2 m) (energy-of-seq2 m)))
    (prn "information" (stats/mean (information_only_bp (information m) (m :pairs))))
    ;; (prn "entropy" (stats/mean (information_only_bp (entropy m) (m :pairs))))
    ;; (prn "mutual info" (stats/mean (mutual_info_only_bp (gutell_mutual_info m) (m :pairs))))
    (prn "mutual info" (stats/mean (mutual_info_only_bp (mutual_info m) (m :pairs))))
    (prn "pairwise identity" (let [id (pairwise_identity (m :seqs))]
                               (double (/ (apply + (map first id))
                                          (apply + (map second id))))))
    (prn "number of seqs" (count (m :seqs)))))

#_(defn main-sto-old
    "DEPRECATED"

    [sto class]
    (let [m (profile sto)]
      (println (stats/mean (zscore m)) ","
               (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ","
               (stats/mean (information_only_bp (information m) (m :pairs))) ","
               ;; (stats/mean (information_only_bp (entropy m) (m :pairs))) ","
               (stats/mean (mutual_info_only_bp (mutual_info m) (m :pairs))) ","
               ;; (stats/mean (mutual_info_only_bp (gutell_mutual_info m) (m :pairs))) ","
               (let [id (pairwise_identity (m :seqs))]
                 (double (/ (apply + (map first id))
                            (apply + (map second id))))) ","
                            (count (m :seqs)) ","
                            class)))

#_(defn main-sto
  "Takes a sto and the class 0 or 1. Performs the various calculations
   for various features. Prints the resulting feat ure vector to the
   console.

   This version does not yet calculate the neutrality. Also need to
   add the name of file in to integrate with future data. "

  [sto class]
  (let [m (profile (read-sto sto))
        ;;gets only the base paired keys
        mi (fn [x]
             ((group-by (fn [[[i j] _]]
                          (contains? (set
                                      (apply concat (map #(vec %) (m :pairs))))
                                     [i j]))
                        x) true))]
    (print (-> (->> (subopt-overlap-sto sto)
                vector
                avg-overlap)
               first
               second
               (get :mean)))
    (print (stats/mean (zscore m)) ",")
    (print (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ",")
    (print (stats/mean (information_only_bp (gutell_calcs/entropy-sto m) (m :pairs))) ",")
    (print (stats/mean (->> (mi (gutell_calcs/mutual_info-sto m)) (into {}) vals)) ",")
    (print (stats/mean (information_only_bp
                        (reduce (fn [x i]
                                  (assoc x i (gutell_calcs/JS
                                              (col->prob (->> (transpose (m :seqs))
                                                              (drop i)
                                                              first)
                                                         :gaps true) :Q (m :background))))
                                {} (range (m :length))) (m :pairs))) ",")
    (print (let [id (pairwise_identity (m :seqs))]
             (double (/ (apply + (map first id))
                        (apply + (map second id))))) ",")
    (print (count (m :seqs)) ",")
    (print class "\n")))

#_(defn cur-make-svm-features
  "Current method to make the train2, train3 csvs. Takes no inputs as
   the program is setup by default. Reads the list of files in
   [pos|neg]/list.txt. Does the calculations for each of the features
   and then prints out the resulting feature vectors and class to the
   file train2.csv. "

  [outfile]
  (io/with-out-writer outfile ;"/home/kitia/bin/gaisr/trainset2/train2.csv"
    (println "zscore, sci, information, MI, JS, pairwise identity, number of seqs, class")
    (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/pos/list.txt")] 
      (let [fdir "/home/kitia/bin/gaisr/trainset2/pos/"]
        (main-sto (str fdir f) 1)))
    (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/neg/list.txt")] 
      (let [fdir "/home/kitia/bin/gaisr/trainset2/neg/"]
        (main-sto (str fdir f) 0)))))

;; (let [negfiles (io/read-lines "/home/kitia/bin/gaisr/trainset/negtrainset.txt")
;;       posfiles (io/read-lines "/home/kitia/bin/gaisr/trainset/postrainset.txt")
;;       remove-empty (fn [files] 
;;                      (remove (fn [x]
;;                                (empty? (get x :cons)))
;;                              (map #(read-sto (str "/home/kitia/bin/gaisr/trainset/" %)) files)))
;;       negf (remove-empty negfiles)
;;       posf (remove-empty posfiles)]
;;   (io/with-out-writer "/home/kitia/bin/gaisr/trainset/trainfeatures2.csv"
;;     (println "zscore, sci, information, mutual information, pairwise ident, number, class")
;;     (doseq [i negf] 
;;       (main-sto i -1))
;;     (doseq [i posf]
;;       (main-sto i 1))))

;; (let [s (find-gaps "GGUAUG.UAUUUC...AA..CCCCA..C..GAUA.AGCCCCGGAA..CU.UA...UU..........G..UGU......U....G.U........GAA.AUAG....AAC")
;;       [ptr buf] (jna-malloc (inc (count s)))]
;;   (jna-invoke Void RNA/read_parameter_file "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par")
;;   (println "e = " (jna-invoke Float RNA/fold s ptr))
;;   (println "struct = " (.getString ptr 0 false)))

#_(defn old-make-svm-features
    "DEPRECATED"
    
    [f]
    (let [m (profile (read-sto f))
          mi (fn [x]
               ((group-by (fn [[[i j] _]]
                            (contains? (set
                                        (apply concat (map #(vec %) (m :pairs))))
                                       [i j]))
                          x) true))]
      (prn "zscore" (stats/mean (zscore m)))
      (prn "sci" (sci (energy-of-aliseq2 m) (energy-of-seq2 m)))
      (prn "information" (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))))
      (prn "mutual info" (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)))
      (prn "pairwise identity" (let [id (pairwise_identity (m :seqs))]
                                 (double (/ (apply + (map first id))
                                            (apply + (map second id))))))
      (prn "number of seqs" (count (m :seqs)))))


#_(defn make-svm-features-old
    "DEPRECATED

     Current method to make the train2, train3 csvs. Takes no inputs
     as the program is setup by default. Reads the list of files in
     [pos|neg]/list.txt. Does the calculations for each of the
     features and then prints out the resulting feature vectors and
     class to the file train2.csv. "

  []
  (do (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/train2.csv"
      (println "zscore, sci, information, MI, JS, pairwise identity, number of seqs, class")
      (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/pos/list.txt")] 
        (let [fdir "/home/kitia/bin/gaisr/trainset2/pos/"
              m (profile (read-sto (str fdir f)))
              mi (fn [x]
                   ((group-by (fn [[[i j] _]]
                                (contains? (set
                                            (apply concat (map #(vec %) (m :pairs))))
                                           [i j]))
                              x) true))]
          (print (stats/mean (zscore m)) ",")
          (print (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ",")
          (print (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))) ",")
          (print (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)) ",")
          (print (stats/mean (information_only_bp (reduce (fn [x i]
                                                          (assoc x i (gutell_calcs/JS (col->prob (->> (transpose (m :seqs))
                                                                                         (drop i)
                                                                                         first)
                                                                                    :gaps true) :Q (m :background))))
                                                        {} (range (m :length))) (m :pairs))) ",")
          (print (let [id (pairwise_identity (m :seqs))]
                   (double (/ (apply + (map first id))
                              (apply + (map second id))))) ",")
          (print (count (m :seqs)) ",1\n")))
    (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/neg/list.txt")] 
      (let [dir "/home/kitia/bin/gaisr/trainset2/neg/"
            sto (str (subs f 0 (- (count f) 3)) "sto")
            m (profile (read-sto (snippet/aln->sto (str dir f) (str dir sto))))
            mi (fn [x]
                 ((group-by (fn [[[i j] _]]
                              (contains? (set
                                          (apply concat (map #(vec %) (m :pairs))))
                                         [i j]))
                            x) true))]
        (print (stats/mean (zscore m)) ",")
        (print (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ",")
        (print (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))) ",")
        (print (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)) ",")
        (print (stats/mean (information_only_bp (reduce (fn [x i]
                                                          (assoc x i (gutell_calcs/JS (col->prob (->> (transpose (m :seqs))
                                                                                         (drop i)
                                                                                         first)
                                                                                    :gaps true) :Q (m :background))))
                                                        {} (range (m :length))) (m :pairs))) ",")
        (print (let [id (pairwise_identity (m :seqs))]
                 (double (/ (apply + (map first id))
                            (apply + (map second id))))) ",")
        (print (count (m :seqs)) ",0\n"))))))
