(ns edu.bc.bio.sequtils.alignment-info
  (:require [edu.bc.fs :as fs]
            [clojure.string :as str])
  (:use edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.sequtils.files
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto read-clj)]))

(def ^{:doc "parse out the function and designation of an RNA from the
    sto file. Typically only used in to parse the Rfam seed alignment
    file. Returns a map for quick lookup"}

  parse-sto-function
    (reduce (fn [m sto]
              (let [gc-lines (first (join-sto-fasta-lines sto ""))
                    get-comment (fn [re] (->> gc-lines
                                             (map #(-> (re-find re %) last) )
                                             (remove empty?)
                                             first))
                    tp (-> (get-comment #"GF TP\s*(.*)") (str/split #" ") set)
                    de (get-comment #"GF DE\s*(.*)")]
                (assoc m (fs/basename sto) {:name de :type tp})))
            {}
            (fs/directory-files (str (fs/homedir) "/bin/gaisr/trainset3/pos") "seed.sto")))

(defn lookup-sto-function
  "takes a sto name and returns the ID and function in a map with
  keywords :name and :type"
  
  [sto]
  (let [nm (str (re-find #"RF\d+.seed" sto)
                ".sto")]
    (parse-sto-function nm)))

(defn alignment-quality
    "find the alignment quality for the negative training
    examples. returns the fraction of base-pairs to the entire
    sequence and the mean number of contiguous bases in the alignment"

    [f]
    (let [{st :cons seqs :seqs} (doall (read-sto f))
          st (first st)
          ;_ (prn f (count (re-seq #"[<|\(]" st)))
          bp-len (double (/ (count (re-seq #"[<|\(]" st))
                            (count st)))
          mean-base-length (mean 
                            (mapcat (fn [s] 
                                      (->> (re-seq #"[ACGU]*" s)
                                           (remove empty? )
                                           (map count )))
                                    seqs))]
      [bp-len mean-base-length]))
