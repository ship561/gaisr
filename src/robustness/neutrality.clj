(ns robustness.neutrality
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure.core.reducers :as r]
            )
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.seq-utils
        ;edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto read-clj
                         parse-dotps)]
        edu.bc.utils.fold-ops
        robustness.utils
        robustness.distance-metrics
        ))

(defn subopt-bpdist-seq
  "Uses RNAdistance to find the difference between 2 structures"
  
  [mut-seq cons n]
  (let [substructs (fold mut-seq :foldtype "RNAsubopt" :n n)
        struct-similarity (fn [cons mut-subopt]
                            (- 1 (/ (bpdist cons mut-subopt) ;1-bpdist/L
                                    (count cons))))] ;length s
    ;;takes percent overlap and
    ;;reduces it to a freqmap to
    ;;save memeory
    (-> (map (fn [mut-subopt]
               (struct-similarity cons mut-subopt))
             substructs)
        frequencies)))

(defn subopt-bpdist-neighbors
  [s cons & {:keys [ncore nsubopt]
             :or {ncore 2 nsubopt 1000}}]
  (let [neighbors (mutant-neighbor s)] ;1-mut neighbors
    (pxmap (fn [neighbor]
             ;;a freqmap of % overlap for each neighbor
             (subopt-bpdist-seq neighbor cons nsubopt))
           ncore
           (concat (list s) neighbors))))

(defn subopt-overlap-neighbors
  "Finds nsubopt suboptimal structures and then finds the percent
  overlap of the suboptimal structures to the consensus
  structure. Returns a list-of-maps where each map is the freqmap of
  the subopt overlap for a 1-mutant neighbor.

  The average is the neutrality"

  [s st & {:keys [ncore nsubopt]
           :or {ncore 2 nsubopt 1000}}]
  (let [[s st cons-keys] (degap-conskeys s st)
        neighbors (mutant-neighbor s)] ;1-mut neighbors
    (->> (pxmap (fn [neighbor]
                  ;;a freqmap of % overlap for each neighbor
                  (subopt-overlap-seq neighbor cons-keys nsubopt))
                ncore
                (concat (list s)
                        neighbors)) ;first element is WT rest are mut neighbors
         (apply merge-with +)
         mean))) 

(defn subopt-overlap-sto
  "This is the main function in the robustness namespace.

  Takes a sto file and finds the suboptimal structure for the WT and
  each of its 1-mutant neighbors. Returns a vector [sto
  list-of-lists-of-maps] where each list-of-maps is the percent
  overlap for a particular sequence and its 1-mutant neighbors. The
  maps are the particular percent overlap for a 1-mutant neighbor."
  
  [sto & {:keys [ncore nsubopt altname]
          :or {ncore 6 nsubopt 1000 altname sto}}]
  (let [{l :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))
        cores (quot ncore
                    (count l))
        cores (if (< cores 1) 1 cores)]
    [altname ;return [filename data]
     (doall
      ;;go over each seq in the alignment
      (map (fn [[nm s]]
             ;;finds 1000 suboptimal structures and
             ;;finds the percent overlap of
             ;;suboptimal structures to the cons struct
             (mean (subopt-overlap-neighbors s cons
                                              :nsubopt nsubopt
                                              :ncore (inc ncore))));ncore
           l))] ;l=list of seqs in the sto
    ))

(defn generic-dist-sto
  "Takes a sto, distance function and optional args. Uses the distance
  function to find the mean distance of the sequences in the sto to
  their 1-mutant neighbors, respectively. "

  [sto distfun & {:keys [ncore nsubopt altname bp]
                  :or {ncore 6, nsubopt 1000,
                       altname sto, bp true}}]
  (let [{l :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))]
    [altname                          ;return [filename data]
     (doall
      ;;go over each seq in the alignment
      (map (fn [[nm s]]
             (apply distfun s cons
                    [:ncore ncore :nsubopt nsubopt :bp bp]))
           l))]))

(defn generic-dist-sto2
  "Takes a sto, distance function and optional args. Uses the distance
  function to find the mean distance of the sequences in the sto to
  their 1-mutant neighbors, respectively. "

  [sto distfn]
  (let [{l :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))
        ]
    (as-> l data
          (r/map #(degap-conskeys (second %) cons) data) 
          (r/map build-mut-neighbors data)
          (r/map #(calculate-dist % distfn) data)
          (r/map mean data)
          (into [] data)
          [sto data])))

(defn main-subopt-overlap
  "Main function for determining neutrality of all sequences and their
  1-mutant neighbors by finding the percent overlap between the
  consensus structure and its suboptimal structures.
  args are provided in a vector of [\"-flag\" value].

  ex
  (main-subopt-overlap \"-f\" \"/home/kitia/bin/gaisr/trainset2/pos/RF00167-seed.3.sto\")"

  [& args]
  (let [[opts _ usage] (cli args
                            ["-f" "--file" "file(s) to check neutrality for"
                             :parse-fn #(str/split #" " %) ;create list of files
                             :default nil]
                            ["-o" "--outfile" "file to write to" :default nil]
                            ["-di" "--dir" "dir in which files are located"
                             :default nil]
                            ["-dfn" "--distfn" "distance function" :parse-fn #(->> % symbol (ns-resolve 'robustness))
                             :default subopt-overlap-neighbors]
                            ["-nc" "--ncore" "number cores to use" :parse-fn #(Integer/parseInt %) :default 6]
                            ["-d" "--debug" "debug test"
                             :default nil :flag true]
                            ["-h" "--help" "usage" :default nil :flag true])
        _ (prn opts)
        fdir (opts :dir)
        fsto (if fdir
               (map #(fs/join fdir  %) (opts :file))
               (opts :file))
        stos (if (opts :debug)
               (take 2 (fs/re-directory-files (opts :dir) #"RF00555-seed*sto"))
               fsto)
        neutrality (for [sto stos] ;loop over stos
                                        ;(subopt-overlap-sto sto :ncore (opts :ncore))
                     (generic-dist-sto sto (opts :distfn) :ncore (opts :ncore))
                     )]
    (cond
     (or (nil? args) (opts :help)) (print usage) ;usage help
     (not (nil? (opts :outfile))) (io/with-out-writer (opts :outfile)
                                    (prn (vec neutrality))) ;data to file
     :else
     (doall neutrality))))
