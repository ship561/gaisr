(ns gutell_calcs
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.stats :as stats]
            [clojure.set :as set]
            [clojure-csv.core :as csv]
            [clojure.java.shell :as shell])
  (:use [edu.bc.bio.sequtils.snippets-files
         :only (read-sto profile)]
        edu.bc.utils
        edu.bc.bio.sequtils.files
        edu.bc.utils.probs-stats))


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

#_(defn entropy-old
    "*DEPRECATED"
    [profile]
    (let [fract-map (profile :fract)
          len (profile :length)
          ;;sum over all bases using -p*log2(p) at position i
          H (fn [p-baseb] 
              (sum
               (map #(* -1 % (Math/log %)) (vals p-baseb))))]
      (reduce (fn [m i]
                (assoc m i (H (second (nth fract-map i)))))
              {} (range len))
      ))

(defn entropy-sto
  "Finds the entropy of each col of the sto file. Returns a map k=col
   number and v=entropy"

  [profile]
  (into {}
        (map (fn [i col]
               (vector i (entropy col :logfn log)))
             (iterate inc 0)
             (transpose (profile :seqs)))))

#_(defn joint_entropy-old
    "DEPRECATED"
    
    [profile]
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
                          (map #(* -1 % (log %)) fr)))))
              {}  all-loc)
      ))

(defn joint_entropy-sto
  "Finds the joint entropy of 2 cols in a sto file. "
  
  [profile]
  (let [inseqs (profile :seqs)
        freqs (->> (transpose inseqs) ;groups cols together
                   (interleave (iterate inc 0)) ;adds col number
                   (partition 2)
                   vec)
        ]
    (reduce (fn [m [[i coll1] [j coll2]]]
              (assoc m [i j]
                     (entropy (joint-probability (fn [a b]
                                                   (map #(str %1 %2) a b)) ;col pairs
                                                 false coll1 coll2)
                              :logfn log)))
            {}  (combins 2 freqs)) ;all combinations of 2 cols
    ))

#_(defn mutual_info_ij-old [Hi Hij i j]
    (let [Hx (get Hi i)
          Hy (get Hi j)
          Hxy (get Hij [i j])]
      (+ Hx Hy (- Hxy) ))
    )

#_(defn mutual_info-old [profile]
    (let [len (count (first (profile :seqs)))
          all-loc (for [i (range len)
                        j (range (+ 4 i) len)]
                    [i j])
          fract-freqs (profile :fract)
          Hx (into {}
                   (map (fn [i col]
                          (vector i (entropy col :logfn log)))
                        (iterate inc 0)
                        (transpose (profile :seqs))))
          Hxy (joint_entropy-old profile)]
      (reduce (fn [m [i j]]
                (assoc m [i j] (mutual_info_ij-old Hx Hxy i j)))
              {}  all-loc)))

(defn mutual_info-sto
  "Finds the mutual information of all pairwise cols in a sto. Returns a map of mutual information"

  [profile]
  (let [len (count (first (profile :seqs)))
        ;;Finds entropy of each col
        Hx (entropy-sto profile)
        Hxy (joint_entropy-sto profile)] ;joint entropy
    (reduce (fn [m [i j]]
              (assoc m [i j] (+ (Hx i) (Hx j) (- (Hxy [i j])))))
            {}  (combins 2 (range len)))))

(defn relative_mutual_info
  "Relative mutual information calculated based on only base pairing bases"

  [profile]
  (let [p profile
        bp {"AU" 1 "UA" 1
            "GU" 1 "UG" 1
            "GC" 1 "CG" 1}
        len (count (first (p :seqs)))
        inseqs (p :seqs)
        all-loc (for [i (range len)
                      j (range (+ 4 i) len)]
                  [i j])
        freqs (vec (partition 2 (interleave (range (count (first inseqs)))
                                            (apply map vector (map #(rest (str/split #"" %)) inseqs)))))
        pij (reduce (fn [m [i j]]
                      (let [fij (frequencies
                                 (map (fn [b1 b2]
                                        (str b1 b2))
                                      (second (nth freqs i)) (second (nth freqs j))))
                            tot (sum (vals fij))
                            fr (map #(/ (second %) tot ) (select-keys fij (keys bp)))] ;;fr = percentages
                        (assoc m [i j]
                               (sum fr))))
                    {}  all-loc)
        qij (reduce (fn [m [i j]]
                      (let [b1 (second (nth freqs i))
                            b2 (second (nth freqs j))
                            fb1 (frequencies b1)
                            fb2 (frequencies b2)]
                        (assoc m [i j]
                               (+ 0.00001 (sum
                                           (map (fn [x y]
                                                  (* (/ (get fb1 x 0) (sum fb1))
                                                     (/ (get fb2 y 0) (sum fb2))))
                                                ["A" "C" "G" "U" "U" "G"] ["U" "G" "C" "A" "G" "U"]))))))
                    {} all-loc)]
    (reduce (fn [m [i j]]
              (assoc m [i j]
                     (* (pij [i j]) (log2 (/ (pij [i j]) (qij [i j]))))))
            {} all-loc)))

(defn R
  "Function returns a map where k=[i j] and v=[covi covj bi bj Hx Hy
   Hxy M Mxyr R1 R2] cov is the covariance at a location"
  
  [profile]
  (let [len (count (first (profile :seqs)))
        fract-map (profile :fract)
        cov (profile :cov)
        Mxy (mutual_info-sto profile)
        Mxyr (relative_mutual_info profile)
        H (entropy profile)
        Hij (joint_entropy-sto profile)
        cov->int (fn [i]
                   (cond
                    (nil? cov)
                    -2
                    (= "." (str (.charAt cov i)))
                    -2
                    (= "?" (str (.charAt cov i)))
                    -1
                    :else
                    (Integer/parseInt (str (.charAt cov i)))))
        commonbase (fn [x] (first (last (sort-by val (second (nth fract-map x))))))
        ]
    (reduce (fn [m [i j]]
              (let [Hx (get H i)
                    Hy (get H j)
                    Hxy (get Hij [i j])
                    M (get Mxy [i j])
                    Mr (get Mxyr [i j])
                    covi (cov->int i)
                    covj (cov->int j)
                    bi (commonbase i)
                    bj (commonbase j)]
                (assoc m [i j]
                       [covi covj bi bj
                        Hx Hy Hxy M Mr
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
  "launches the default way to calcuate covi covj bi bj Hx Hy Hxy M
   Mxyr R1 R2 class pval. It will print out a default file at the
   current time."
  
  [sto fpval]
  (let [p (profile (read-sto sto))
        mi (R p)
        loc (p :pairs)
        mi (remove_gap_col mi (fractpairs_contain_nogaps p) 0.5) ;;remove gapped columns
        x (group-by (fn [[[i j] _]]
                      (contains? (set
                                  (apply concat (map #(vec %) loc)))
                                 [i j]))
                    mi)
        y (merge (reduce (fn [m a]
                           (assoc m (key a) (conj (val a) 0)))
                         {} (into {} (x false)))
                 (reduce (fn [m a]
                           (assoc m (key a) (conj (val a) 1)))
                         {} (into {} (x true))))
        pvals (reduce (fn [m [[x y z]]]
                        (assoc m [(Integer/parseInt x) (Integer/parseInt y)] z))
                      {} (map #(csv/parse-csv (str/replace-re #" " "" %)) (io/read-lines fpval)))]
    (print_out (merge-with conj y (select-keys pvals (keys y))) (str (subs sto 0 (- (count sto) 3)) "all.csv"))
    ))

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
  (let [s (-> (group-by (fn [x] (= x '("CLUSTAL W (SISSIz 0.1 simulation)"))) 
                        (partition-by #(= % "CLUSTAL W (SISSIz 0.1 simulation)")
                                      (str/split-lines
                                       ((shell/sh "SISSIz" "--rna" "-s" "-n" (str n) f) :out))))
              second second)
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


(defn pval
  "Reads in a sto file and n = an integer. The sto file is used to
   generate a Clustal W type file with the same file name with an .aln
   extension instead of .sto.  Produces .aln and generates n random
   shuffled alignments and then calcuates p value by taking the
   fraction sampled alignment MI > true MI. returns a position-pair
   and the pvalue"

  [stoin n]
  (let [aln (sto->aln stoin (str (subs stoin 0 (- (count stoin) 3)) "aln"))
        p (profile stoin)
        mi (mutual_info-sto p)
        rand_mi (apply concat (doall
                               (pmap (fn [randa]
                                       (doall (map #(mutual_info-sto (profile %)) randa)))
                                     (partition-all (/ n 10) (rand_aln aln n)))))]
    (for [k (keys mi)]
      [k (double (/ (count
                     (filter (fn [x]
                               (let [sample-mi (get x k)]
                                 (> sample-mi (get mi k))))
                             rand_mi))
                    n))])))


               
                                             
(defn KL
  "kullback leibler divergence takes in a column and a model
   distribution Q in the form of a key value map right answer for
   sample.sto = (0.8239592165010823 0.0 1.3862943611198906
   1.3862943611198906 1.3862943611198906 1.3862943611198906
   0.8239592165010823 0.34657359027997264)"

  [P & {Q :Q :or {Q (into {}
                          (map #(vector % 0.25) #{\A \C \G \U "A" "U" "G" "C"})) }}]                                              
  (sum (for [i (keys P)]
         (let [p (P i)
               q (if-not (nil? (Q i)) ;;determine background probabilites
                               ;;for 1 base or multiple bases
                   (Q i)
                   (reduce (fn [x y]
                             (* x (get Q y)))
                          1 (seq i)))]
           (if (or (zero? p) (zero? q))
             0
             (* p (Math/log (/ p q))))))))

(defn JS [P & {Q :Q :or {Q (into {} (map #(vector % 0.25) #{\A \C \G \U "A" "U" "G" "C"})) }}]
  (let [Q (if (char? (first (keys P))) ;;determine background probabilites for 1
                         ;;base or multiple bases
            Q
            (reduce (fn [m k]
                      (assoc m k
                             (reduce (fn [x y]
                                       (* x (get Q (str y))))
                                     1 (seq k))))
                    {} (keys P)))
        M (reduce (fn [m k]
                    (assoc m k (/ (+ (P k) (Q k)) 2)))
                  {} (keys P))]
    ;(prn M P)
    (+ (* 0.5 (KL P :Q M)) (* 0.5 (KL Q :Q M)))))

