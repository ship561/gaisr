(ns edu.bc.utils.fold-ops
  (:require [clojure.contrib.io :as io]
            [clojure.string :as str]
            [clojure.java.shell :as shell]
            [edu.bc.fs :as fs])
  (:use refold
        [edu.bc.bio.sequtils.files :only [join-sto-fasta-file]]
        [slingshot.slingshot :only [throw+]]))

(def param-file (let [viennadir (if (fs/directory? "/usr/local/ViennaRNA/")
                                  "/usr/local/ViennaRNA/" 
                                  (fs/join (fs/homedir) "/bin/ViennaRNA/"))
                      pfile (str viennadir "rna_andronescu2007.par")
                      pfile (if (fs/exists? pfile)
                              pfile
                              (str viennadir "/misc/rna_andronescu2007.par"))]
                  (if (fs/exists? pfile)
                    pfile
                    (throw+ {:file pfile} "parameter file not exist" ))))

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
                                             "-P" param-file
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
                                       (->> (map #(vector %1 %2)
                                                 (range (count s))
                                                 (repeat "."))
                                            (into {}))
                                       ;;keeps bases that have over 50%
                                       ;;representation in suboptimal structures
                                       (filter (fn [[[i j] p]] 
                                                 (and (< (+ i 3) j)
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

;(ns-unmap 'edu.bc.utils.fold-ops 'fold2)
(defmulti fold2 (fn [s & args]
                  ((or (first args) {}) :foldmethod)))

(defmethod fold2 :RNAfold [s args]
  (-> ((shell/sh "RNAfold"
                 "-P" param-file
                 "--noPS"
                 :in s )
       :out)
      (str/split-lines)
      second
      (str/split #" ")
      first))

(defmethod fold2 :RNAsubopt [s args] 
  (->> ((shell/sh "RNAsubopt"
                  "-p" (str (args :n))        ;samples according to
                                        ;Boltzmann distribution
                  :in s)
        :out)
       str/split-lines
       (remove #(re-find #"[^\(\)\.]" %))))

(defmethod fold2 :centroid [s args]
  (first (suboptimals s (args :n))))

(defmethod fold2 :default [s]
  (fold2 s {:foldmethod ::RNAfold}))

(defn fold
  "Folds a sequence of RNA and returns only the target
   structure. Target structure can either be centroid or MFE."
  
  [s & {:keys [foldtype n]
        :or {foldtype "mfe" n 10000}}]
  (case foldtype
    "mfe"
    (-> ((shell/sh "RNAfold"
                    "-P" param-file
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

(defn align-fold
  "Takes a fasta file and outputs an alignment with structure in
  stockholm file format. The outfile is the input file renamed to have
  a .sto extension. Uses mxscarna to do the fold/align"

  ([fna]
     (align-fold fna (fs/replace-type fna ".sto")))
  
  ([fna outfile]
     (let [tmp (fs/tempfile)
           call (shell/sh "mxscarna" "-stockholm" fna)]
       (io/with-out-writer tmp
         (-> call :out println))
       (join-sto-fasta-file tmp outfile :origin "#=GF AU mxscarna")
       outfile)))

;;;(ns-unmap 'edu.bc.utils.fold-ops 'fold-aln)
(defmulti

  ^{:doc "multimethod for folding an alignment. can use
  either :RNAalifold or :centroid_alifold as the first argument or
  none in args. Defaults to :RNAalifold. Always requires a file
  name (aln)."
    :arglists '([aln] [foldprogram aln])}
  
  fold-aln (fn [& args]
             (first args)))

(defmethod fold-aln :RNAalifold [_ aln]
  (-> ((shell/sh "RNAalifold"
                 "-P" param-file
                 "-r" "--noPS" aln)
       :out)
      str/split-lines
      second
      (str/split #"\s")
      first))

(defmethod fold-aln :centroid_alifold [_ aln]
  (-> ((shell/sh "centroid_alifold" aln)
       :out)
      str/split-lines
      last
      (str/split #"\s")
      first))

(defmethod fold-aln :default [aln]
  {:pre [(fs/exists? aln)]}
  (fold-aln :RNAalifold aln))

(comment 
  (defn fold-aln [aln]
    (-> ((shell/sh "RNAalifold"
                   "-P" param-file
                   "-r" "--noPS" aln)
         :out)
        (str/split-lines)
        second
        (str/split #"\s")
        first)))

(defn bpdist
  "finds the distance between 2 structures. uses tree edit distance by
  default. when bpdist = true, uses base pair distance"

  [st1 st2 & {:keys [bpdist]}]
  (->> ((shell/sh "RNAdistance" (if bpdist "-DP" "")
                  :in (str st1 "\n" st2))
        :out)
       (re-find #"\d+" )
       (Integer/parseInt)))
