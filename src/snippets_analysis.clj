(ns snippets-analysis
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.set :as sets]
            [edu.bc.fs :as fs]
            [clojure.core.reducers :as r]
           )
  (:use robustness
        refold
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.snippets-math
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto)]
        edu.bc.utils.fold-ops))

(def ^{:private true} homedir (fs/homedir))

(defmacro with-out-appender [f & body]
  `(with-open [w# (clojure.java.io/writer ~f :append true)]
     (binding [*out* w#] ; I forgot this bit before.
       ~@body)))

(defn lazy-file-lines [file]
  (letfn [(helper [rdr]
            (lazy-seq
             (if-let [line (.readLine rdr)]
               (cons line (helper rdr))
               (do (.close rdr) nil))))]
    (helper (clojure.java.io/reader file))))



;;;-----------------------------------------------------------------------------
;;;
;;;analytics for robustness.clj. trying to establish the average
;;;overlap of the suboptimal structures is a good representation of
;;;the distance between the WT and mutant structures.
;;;
;;;-----------------------------------------------------------------------------

(defn jsd-wt-neighbor
  "Compares the distribution of base-pairs and gap chars in each
   column between the wt and 1-mutant neighbor. Also computes the
   basepair distance between the wt and 1-mutant neighbor. Takes a seq
   and neighbors, consensus structure keys and
   n=number_suboptimal_structures. Returns a list of vectors where
   each vector is [mutant-name [ith-col jsd(wt,mut) %overlap
   bpdistance]]."
  
  [wt neighbors cons cons-keys n]
  (let [wt-probs (map #(probs 1 %) (transpose (fold wt :foldtype "RNAsubopt" :n n)))]
    (pxmap (fn [[nm neighbor]]
             (let [mut-probs (map #(probs 1 %) (transpose (fold neighbor :foldtype "RNAsubopt" :n n)))
                   overlap (double (mean (subopt-overlap-seq neighbor cons-keys n)))
                   bpdist (double (mean (subopt-bpdist-seq neighbor cons n)))]
               ;;find the jsd between the wt col_i and mut col_i
               (map (fn [i c1 c2]
                      [nm [i (jensen-shannon c1 c2) overlap bpdist]])
                    (iterate inc 0) wt-probs mut-probs)))
           100
           neighbors)))

(defn driver-jsd-wt-neighbor
  "for purposes of graphing."
  
  [sto]
  (let [;;sto (str homedir "/bin/gaisr/trainset2/pos/RF00555-seed.1.sto")
        {sqs :seqs cons :cons} (read-sto sto :with-names true)
        cons (change-parens (first cons))
        sto-nm (fs/basename sto)] 
    (io/with-out-writer (str homedir  "/bin/gaisr/robustness/temp.txt")
      (println "sto-name,seq-name,mutname, pos, jsd, overlap, bpdist")
      (doseq [[seq-nm s] sqs]
        (let [[wt st] (remove-gaps s cons)
              cons-keys (set (keys (struct->matrix st)))
              n 1000
              neighbors (into {} (mutant-neighbor wt :with-names true))]
          (doseq [i (jsd-wt-neighbor wt neighbors cons cons-keys n)
                  j i]
            (println (str/join "," (flatten (concat [sto-nm seq-nm] j))))))))))

(defn combine-jsd
  "read output from driver-jsd-wt-neighbor and finds the mean jsd of
   all cols for a given mut. the return is a vector of vectors where
   each vector contains [sto-name, seq-name, mut-name, %overlap,
   mean(jsd)]."
  
  []
  (->> (rest (io/read-lines (str homedir "/bin/gaisr/robustness/temp.txt")))
       (map #(str/split #"," %) )
       ;;merges the jsds for each col with the same mut-name and overlap
       (reduce (fn [m [sto-name seq-name mut-name _ jsd overlap bpdist]]
                 (let [jsd (Double/parseDouble jsd)
                       overlap (Double/parseDouble overlap)
                       bpdist (Double/parseDouble bpdist)
                       k [sto-name seq-name mut-name overlap bpdist jsd]
                       v (get-in m k 0)]
                   (assoc-in m k (inc v))))
               {} )
       ;;mean jsd
       (reduce (fn [m x]
                 (let [[sto-name x] x]
                   (->> (for [
                              [seq-name x] x
                              [mut-name x] x
                              [overlap x] x
                              [bpdist jsd] x]
                          [sto-name seq-name mut-name overlap bpdist (mean jsd)])
                        vec
                        (conj m ))))
               [] )
       first))

(defn compare-col-jsd
  "visualize data from an input sto. Compares a seq to its 1-mutant
   neighbor by comparing the jsd of a given column i in the suboptimal
   alignment. High jsd would indicate the column in the mutant is
   significantly different than the expected based on wt. Prints the
   data out to repl for examination."
  
  ([sto]
     (let [sto (str homedir "/bin/gaisr/trainset2/pos/RF00555-seed.1.sto")
           {sqs :seqs cons :cons} (read-sto sto :with-names true)
           cons (change-parens (first cons))]
       (prn sto)
       (doseq [[seqnm s] sqs]
         (let [[wt st] (remove-gaps s cons)
               cons-keys (set (keys (struct->matrix st)))
               n 1000
               neighbors (into {} (mutant-neighbor wt :with-names true))
               loi (->> (jsd-wt-neighbor wt neighbors cons-keys n)
                        ;;process the jsd to find col of interest
                        #_(map (fn [jsd]
                               (remove #(or (< (-> % second second) 0.8) ;jsd
                                            (> (-> % second third) 0.4)) ;%overlap
                                       jsd))
                             )
                        (remove #(empty? %) )
                        (map vec ))
               ]
           (prn seqnm)
           (doseq [x loi
                   [mutnm [pos h overlap]] x]
             (prn "wt" pos)
             (prn s)
             (prn (apply str (repeat 10 "0123456789")))
             (doseq [i (fold s :foldtype "RNAsubopt" :n 3)]
               (prn i))
             (prn seqnm mutnm "pos" pos "jsd" h "overlap" overlap)
             (prn (neighbors mutnm))
             (prn (apply str (repeat 10 "0123456789")))
             (doseq [i (fold (neighbors mutnm) :foldtype "RNAsubopt" :n 3)]
               (prn i)))))))
  
  ([sto outfile]
     (let [{l :seqs cons :cons} (read-sto sto :with-names true)
           cons (change-parens (first cons))]
       (io/with-out-writer outfile
         (compare-col-jsd sto)))))


(comment
  
  (let [[remaining-file & remaining-files]
        (filter #(re-find #".1.sto" %) 
                (fs/listdir (str homedir "/bin/gaisr/trainset2/pos/")))
        fdir (str homedir "/bin/gaisr/trainset2/pos/")
        odir (str homedir "/bin/gaisr/robustness/")
        outfn (fn [x] (doseq [out x] (println (str/join "," out))))]
    (driver-jsd-wt-neighbor (str fdir remaining-file))
    (io/with-out-writer (str odir "temp3.txt") 
      (println "sto-name,seq-name,mut-name,overlap,bpdist,jsd")
      (outfn (combine-jsd)))
    (doseq [sto remaining-files]
      (driver-jsd-wt-neighbor (str fdir sto))
      (with-out-appender (str odir "temp3.txt")
        (outfn (combine-jsd)))))


  )
