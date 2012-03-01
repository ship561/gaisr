(ns gutell_calcs
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.stats :as stats]
            [incanter.core :as math]
            [clojure.set :as set]
            [clojure-csv.core :as csv]
            [clj-shell.shell :as shell]
            [clojure.contrib.seq :as seq]
            [edu.bc.fs :as fs])
  (:use [consensus_seq
         :only [read-sto profile]]
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [edu.bc.utils]
        [edu.bc.bio.seq-utils]))

(defn fractpairs_contain_nogaps
  "Checks pairs of columns x and y to see how many contain
   gaps. returns a map of the percent of gaps in each pair of
   columns."
  
  [profile]
  (let [inseqs (profile :seqs)
        len (count (first inseqs))
        freqs (partition 2 (interleave (range (count (first inseqs)))
                                       (apply map vector (map #(rest (str/split #"" %)) inseqs))))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])]
    (reduce (fn [m [i j]]
              (let [pairs (map (fn [b1 b2]
                               (str b1 b2))
                             (second (nth freqs i)) (second (nth freqs j)))
                    fr (map #(re-find #"\.|-" %) pairs)] ;;fr = # of pairs containing gaps
                (assoc m [i j]
                       (double (/ (count (filter nil? fr)) (count fr)))
                       )))
            {}  all-loc)
    ))

(defn remove_gap_col
  "removes the columns from a caculation where there are greater than
   'per' gaps. m = map, gap_col = column gap percents, per = threshold
   below which columns are removed.  Returns a map"

  [m gap_col per]
  (select-keys m (keys (remove #(< (val %) per) gap_col)))
  )

(defn entropy [profile]
  (let [fract-map (profile :fract)
        len (count (first (profile :seqs)))
        ;;sum over all bases using -p*log2(p) at position i
        entropy_i (fn [p-baseb i] 
                    (sum
                     (map #(* -1 % (Math/log %)) (vals p-baseb))))]
    (reduce (fn [m i]
              (assoc m i (entropy_i (second (nth fract-map i)) i)))
            {} (range len))
    ))

(defn joint_entropy [profile]
  (let [inseqs (profile :seqs)
        len (count (first inseqs))
        freqs (vec (partition 2 (interleave (range (count (first inseqs)))
                                       (apply map vector (map #(rest (str/split #"" %)) inseqs)))))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])]
    (reduce (fn [m [i j]]
              (let [fij (frequencies
                         (map (fn [b1 b2]
                                (str b1 b2))
                              (second (nth freqs i)) (second (nth freqs j))))
                    tot (sum (vals fij))
                    fr (map #(/ (second %) tot ) fij)] ;;fr = percentages
                (assoc m [i j]
                       (sum
                        (map #(* -1 % (Math/log %)) fr)))))
            {}  all-loc)
    ))

(defn mutual_info_ij [Hi Hij i j]
  (let [Hx (get Hi i)
        Hy (get Hi j)
        Hxy (get Hij [i j])]
    (+ Hx Hy (* -1 Hxy) ))
  )

(defn mutual_info [profile]
  (let [len (count (first (profile :seqs)))
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])
        fract-freqs (profile :fract)
        Hx (entropy profile)
        Hxy (joint_entropy profile)]
    (reduce (fn [m [i j]]
              (assoc m [i j] (mutual_info_ij Hx Hxy i j)))
            {}  all-loc)))

(defn mutual_info_only_bp
  "returns a vector of mutual information only at positions where
   there is base pairing provided by alignment."
  
  [info bp-loc]
  (filter (fn [[[i j] v]]
            (contains? (set
                        (apply concat (map #(vec %) bp-loc)))
                       [i j]))
          info))

(defn mutual_info_not_bp
  "returns a vector of mutual information only at positions where
   there is base pairing provided by alignment."
  
  [info bp-loc]
  (remove (fn [[[i j] v]]
            (contains? (set
                        (apply concat (map #(vec %) bp-loc)))
                       [i j]))
          info))

(defn R
  "Function returns a map where k=[i j] and v=[Hx Hy M R1 R2]"
  
  [profile]
  (let [len (count (first (profile :seqs)))
        fract-map (profile :fract)
        Mxy (mutual_info profile)
        H (entropy profile)]
    (reduce (fn [m [i j]]
              (let [Hx (get H i)
                    Hy (get H j)
                    M (get Mxy [i j])]
                (assoc m [i j]
                       [Hx Hy M
                        (if (zero? Hx) 0 (/ M Hx))
                        (if (zero? Hy) 0 (/ M Hy))])))
            {} (for [i (range len)
                     j (range (+ 4 i) len)]
                 [i j]))))

(defn print_out
  "prints out the information, mutual information, R, outputs in a csv
   type format. if no out file is provided, then it prints to the
   repl, otherwise, it will print to specified file"
  
  ([info]
  (doseq [[[i j] x] info]
    (print (csv/write-csv (vector (map str (conj (apply list x) j i)))))
      ))

  ([info outfile]
     (io/with-out-writer outfile (print_out info))))

(defn main
  "launches the default way to calcuate Hx Hy M R1 R2. It will print
  out a default file at the current time."
  
  [f]
  (let [p (profile (read-sto f))
        mi (R p)
        loc (p :pairs)
        mi (remove_gap_col mi (fractpairs_contain_nogaps p) 0.5) ;;remove gapped columns
        x (map #(mutual_info_only_bp mi %) loc) ;;intersection between
        ;;mi and loc
        y (map #(mutual_info_not_bp mi %) loc)] ;;difference between
    ;;mi and loc
    (print_out (apply set/union (map set x)) (str (subs f 0 (- (count f) 3)) "onlybp.csv"))
    (print_out (apply set/intersection (map set y)) (str (subs f 0 (- (count f) 3)) "notbp.csv"))))

;; (let [p (profile (read-sto "/home/kitia/Downloads/S15_101711UBedit.sto"))
;;                    mi (R p)
;;                    loc (p :pairs)
;;                    mi (remove_gap_col mi (fractpairs_contain_nogaps p) 0.5) ;;removed gapped columns
;;                    x (map #(mutual_info_only_bp mi %) loc)]
;;                    (for [i x]
;;                         (let [mis (map (fn [[_ [_ _ a _]]] a) i)]
;;                              [(stats/mean mis) (stats/sd mis) (apply min mis) (apply max mis) (count mis)])))

(defn rand_aln
  "Generates n random alignments based on the input alignment, f, in Clustal format.
   This function produces a map similar to the read-sto function. The
   map contains {:seqs :cons :file}"
  
  [f n]
  (let [s (second
           (seq/separate (fn [x] (= x '("CLUSTAL W (SISSIz 0.1 simulation)"))) 
                         (partition-by #(= % "CLUSTAL W (SISSIz 0.1 simulation)")
                                       (str/split-lines
                                        ((shell/sh "SISSIz" "--rna" "-s" "-n" (str n) f) :out)))))
        recombine-lines (fn [x]
                          (reduce (fn [m l]
                                    (let [[nm sq]
                                          (str/split #"\s+" l)
                                          prev (get m nm [(gen-uid) ""])]
                                      (if-not (empty? nm)
                                        (assoc m  nm [(first prev)
                                                      (str (second prev) sq)])
                                        m)))
                                  {} x)) ]
    (map (fn [x]
           (let [m (recombine-lines x)]
             (assoc {} :seqs (map (fn [[n [uid sq]]]
                                    sq)
                                  m)
                    :cons (list (apply str
                                 (repeat (count (second (second (first m)))) ".")))
                    :file "")))
         s)))

(defn rand_aln2
  "Generates n random alignments based on the input alignment, f, in Clustal format.
   This function produces a map similar to the read-sto function. The
   map contains {:seqs :cons :file}"
  
  [f n]
  (let [temp (fs/tempfile)
        sissiz (while (fs/empty? t)
                 (runx "SISSIz" "--rna" "-s" "-n" (str 2) "/home/peis/bin/gaisr/l10-reduced.aln" :> t :?> t))
        sissiz (rest (str/split #"CLUSTAL W \(SISSIz 0.1 simulation\)" (sissiz :out)))
        recombine-lines (fn [x]
                          (reduce (fn [m l]
                                    (let [[nm sq]
                                          (str/split #"\s+" l)
                                          prev (get m nm [(gen-uid) ""])]
                                      (if-not (empty? nm)
                                        (assoc m  nm [(first prev)
                                                      (str (second prev) sq)])
                                        m)))
                                  {} (drop-while empty? (str/split-lines x))))]
    (io/with-out-writer temp
      (doseq [i sissiz]
        (doseq [[nm [uid sq]] (sort-by #(first (second %)) (recombine-lines i))]
          (println nm uid sq))
        (println ">")))
    temp))

(defn pval
  "Reads in a sto file and n = an integer. The sto file is used to
   generate a Clustal W type file with the same file name with an .aln
   extension instead of .sto.  Produces .aln and generates n random
   shuffled alignments and then calcuates p value by taking the
   fraction sampled alignment MI > true MI. returns a position-pair
   and the pvalue"

  [stoin n]
  (let [aln (sto->aln stoin (str (subs stoin 0 (- (count stoin) 3)) "aln"))
        p (profile (read-sto stoin))
        mi (mutual_info p)
        rand_mi (map #(mutual_info (profile %)) (rand_aln aln n))]
    (for [k (keys mi)]
      [k (double (/ (count
                     (filter (fn [x]
                               (let [sample-mi (get x k)]
                                 (> sample-mi (get mi k))))
                             rand_mi))
                    n))])))

(defn pval2
  "Reads in a sto file and n = an integer. The sto file is used to
   generate a Clustal W type file with the same file name with an .aln
   extension instead of .sto.  Produces .aln and generates n random
   shuffled alignments and then calcuates p value by taking the
   fraction sampled alignment MI > true MI. returns a position-pair
   and the pvalue"

  [stoin n]
  (let [aln (sto->aln stoin (str (subs stoin 0 (- (count stoin) 3)) "aln"))
        p (profile (read-sto stoin))
        mi (mutual_info p)
        rand_mi (apply concat (doall
                               (pmap (fn [randa]
                                       (doall (map #(mutual_info (profile %)) randa)))
                                     (partition-all (/ n 10) (rand_aln aln n)))))]
    (for [k (keys mi)]
      [k (double (/ (count
                     (filter (fn [x]
                               (let [sample-mi (get x k)]
                                 (> sample-mi (get mi k))))
                             rand_mi))
                    n))])))
