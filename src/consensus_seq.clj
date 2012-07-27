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
        libsvm2weka
        edu.bc.bio.sequtils.files
        edu.bc.utils))

(def possible_pairs {"AU" 1 "UA" 1 "GC" 1 "CG" 1 "GU" 1 "UG" 1} )
(def base #{"A" "C" "G" "U" "." "-"} )
(def all_pairs {"AA" 1 "AC" 1 "AG" 1 "AU" 1
                "CA" 1 "CC" 1 "CG" 1 "CU" 1
                "GA" 1 "GC" 1 "GG" 1 "GU" 1
                "UA" 1 "UC" 1 "UG" 1 "UU" 1} )

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
  (let [[gc-lines seq-lines cons-lines] (join-sto-fasta-lines f "")
        cov (first (map #(last (str/split #"\s+" %))
                (filter #(.startsWith % "#=GC cov_SS_cons") gc-lines)))
        cl (map #(last (second %))
                (filter #(.startsWith (first %) "#=GC SS_cons") cons-lines))
        sl (reduce (fn [v [_ [_ sq]]]
                     (conj v (str/replace-re #"T" "U" (.toUpperCase sq))))
                [] seq-lines)]
    (assoc {} :seqs sl :cons cl :file f :cov cov)))

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

(defn energy-of-aliseq [profile]
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
              (assoc m k (/ v (sum freq-map)))
              m))
          {} freq-map))
 
(defn mutual_info_only_bp [info bp-loc]
  "returns a vector of mutual information only at positions where
   there is base pairing provided by alignment."

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

(defn profile [m]
  (let [struct (map #(change-parens %) (m :cons))
        s (m :seqs)
        freqs (partition 2 (interleave (range (count (first s)))
                                       (map frequencies
                                            (apply map vector (map #(rest (str/split #"" %)) s)))))
        fract-freqs (sort-by key (reduce (fn [l [n freq-map]]
                                           (assoc l n (fraction-base-b freq-map)))
                                         {} freqs))
        q (col->prob (flatten (map #(rest (str/split #"" %)) s)) :gaps true)
        pairs (map #(refold/make_pair_table %) struct)]
  {:seqs s
   :structure struct
   :fract fract-freqs
   :background q
   :pairs pairs
   :cov (m :cov)
   :filename (m :file)
   :length  (count (first s))}))

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

(defn main-sto [sto class]
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
