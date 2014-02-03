(ns robustness.analysis
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [incanter.charts :as charts]
            )
  (:use edu.bc.utils.probs-stats
        [edu.bc.utils.snippets-math :only (median-est)]
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto read-clj
                         parse-dotps)]
        [incanter.core :only (view)]
        ))

;;;-----------------------------------
;;;Section for visualizing data
(defn seq-neutrality
  "Takes a list of sets of dists. The dist is defined as the fraction
  of base-pairs in the given structure retained in the structure of a
  mutant neighbor of the given seq. The subopt overlap is (mean dist1
  .. disti). Neutrality = (mean subopt-overlap1 .. subopt-overlap3n).
  The neutrality indicates the ability of the sequence to maintain
  structure despite mutations.

  Robustness = neutrality(wt) > mean neutrality(1-muts) "

  [list-dists]
  (let [subopt-overlap (fn [dist] (mean dist))]
    (-> (map subopt-overlap list-dists) ;calc subopt overlap of each 1-mut 
        mean) ;neutrality
    ))

(defn sto-neutrality
  "takes output from main-subopt-overlap. Finds the neutrality for
  each seq in the sto. returns summary stats for a list of sto as a
  whole."

  [x]
  (reduce (fn [m [nm vals]]
            (let [neut (map seq-neutrality vals)] ;neutrality of
                                        ;each seq in a sto
              ;;summary stats of neutrality for sto
              (assoc m nm {:median (median neut)
                           :mean (mean neut) ;sto mean neutrality
                           :sd (sd neut)
                           :raw neut})))
          {} x))

(defn avg-overlap
  "Takes a map of percent overlaps where it is organized in [k v]
   pairs. k=file name and v=list of lists of frequency maps of percent
   overlap of 1000 suboptimal structures for the WT and each of its
   1-mutant neighbors. Returns a map of maps of the summary
   stats. Equivalent of the weighted mean."
  
  [map-of-per-overlaps]
  (reduce (fn [m [k list-lists-maps]]
            (let [list-maps (->> list-lists-maps 
                                (apply concat) ;combines data/lists from all sequences
                                (apply merge-with +)) ;combines the
                                        ;freqmaps into 1 map
                  avg (mean list-maps)
                  sd (sd list-maps)
                  med (double (median-est list-maps))]
              (assoc m k {:med med :mean avg :sd sd})))
          {} map-of-per-overlaps))

(defn overlap-per-seq
  "Takes an entry from teh map-of-lists-of-lists-of-maps data
   structure. This entry has a format of [nm data]. It averages the
   %overlap for all mutations at a position. Returns a vector of
   average %overlap by position."

  [map-of-per-overlaps]
  (let [[nm data] map-of-per-overlaps]
    (map (fn [ea-seq]
           (let [wt (mean (first ea-seq))
                 ;;merges the 3 mutations at a position and finds the mean
                 muts (for [mut (partition-all 3 (rest ea-seq))]
                        (->> mut
                             (apply merge-with +) ;merge freqmaps together
                             mean))]
             (vec (cons wt muts)))) ;returns vector starting with wt
                                    ;then muts
         data)))

(defn chart-overlap-sto
  "Takes an entry from the map-of-lists-of-lists-of-maps data
   structure. It calls the overlap-per-seq to find the avarge
   %overlaps at each position. Returns a graph of %overlap for each
   sto. Each line on graph represents 1 sequence from the sto."

  [map-of-per-overlaps & {:keys [title]
                          :or {title "neg RF0555.4"}}]
  (let [lines (overlap-per-seq map-of-per-overlaps)
        l (charts/xy-plot (range 200) (first lines) :title title :series 1 :legend true
                          :x-label "position" :y-label "mut % overlap with cons")]
    (view l)
    (map (fn [i y]
           (charts/add-lines l (range 200) y :series-label i))
         (iterate inc 2) (rest lines))))

(defn chart-overlap-sto2
  "Similar to the chart-overlap-sto except it makes a data structure
   which contains the percent overlap at each point. The lengths are
   the same by adding back the gaps back into the seq. The gaps are
   the same value (overlap) as the point directly before it. It will
   either print out in csv format to a file specified (outcsv) or just
   return the vector-of-lists where each list contains the percent
   overlap of each point.

  Once the csv is created, it can be graphed using
  robustness/chart-overlap-sto.R In R, first source the file then
  'chartoverlap(file)' to execute to produce a graph."

  [map-of-per-overlaps & {:keys [outcsv]}]
  (let [fsdrop (fn [n file] (->> (fs/split file)
                                rest
                                (drop n)
                                (apply fs/join)))
        sto (first map-of-per-overlaps)
        sto-file (str (fs/homedir) "/" (fsdrop 2 sto)) 
        seq-names (map first (-> (read-sto sto-file :info :both) :seqs))
        lines (map (fn [original-seq pts-avg]
                     (loop [os original-seq
                            pa pts-avg
                            new-pts-avg []]
                       (if (seq os)
                         (if (= (first os) \.) 
                           (recur (rest os)
                                  pa
                                  (conj new-pts-avg (first pa)))
                           (recur (rest os)
                                  (rest pa)
                                  (conj new-pts-avg (first pa))))
                         new-pts-avg)))
                   (-> (read-sto sto-file) :seqs)
                   (overlap-per-seq map-of-per-overlaps))
        linecsv (map (fn [nm v]
                       (str/join "," (cons nm v)))
                     seq-names lines)]
    (if outcsv
      (io/with-out-writer outcsv
        (println (str/join "," (->> lines first count range (cons "name"))))
        (doseq [x linecsv] (println x)))
      lines)))



;;;---------------------------------------------------
