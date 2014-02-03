(ns robustness.robustness
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            ;[clojure.core.reducers :as r]
            )
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.seq-utils
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto read-clj
                         parse-dotps)]
        robustness.utils
        robustness.neutrality
        robustness.distance-metrics
        ))

(def ^{:private true} homedir (fs/homedir))

;;;-----------------------------------------------------------------------------
;;;suboptimal significance
;;;-----------------------------------------------------------------------------
(defn valid-seq-struct
  "Checks the sto to make sure that all sequences in the file can form
   part of the consensus structure. This is useful when ensuring that
   the shuffled stos will form valid structures. Occassionally,
   sequences in the sto will not fold into the consensus structure and
   cause other functions to fail. Returns true if all sequences can
   fold into part of the cons structure."
  
  [sto]
  (let [{sqs :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))
        valid? (fn [st] (pos? (count st)))]
    (every? true? (map (fn [[_ s]]
                         (let [[_ st cons-keys] (degap-conskeys s cons)]
                           (valid? cons-keys)))
                       sqs))))

(defn subopt-significance
  "Takes an input sto and estimates the significance of the suboptimal
   structure overlap when compared to n random simulated
   alignments. Returns the map-of-lists-of-lists-of-maps."
  
  [insto & {:keys [ncores nsamples]
            :or {ncores 3 nsamples 3}}]
  (concat [(subopt-overlap-sto insto :altname (fs/basename insto))] ;wt subopt overlap
          ;;shuffled version of sto
          (pxmap (fn [randsto]                   
                   (let [result (subopt-overlap-sto randsto :altname (fs/basename randsto))] ;find subopt overlap
                     (fs/rm randsto)
                     result))
                 ncores
                 (take nsamples
                       ;;only keep valid stos where the sequences can
                       ;;form a part of the consensus structure
                       (filter #(true? (valid-seq-struct %)) 
                               (repeatedly #(sto->randsto insto (fs/tempfile))))) ;create random stos
                 )))

;;;-----------------------------------------------------------------------------
;;;suboptimal robustness
;;;-----------------------------------------------------------------------------
(defn subopt-robustness-summary
  "avg-subopt = (subopt-robustness (str fdir insto) n). It is the
   list-of-lists average subopt overlap of 1-mut
   structures (neutrality)"

  [avg-subopt]
  (let [data (remove nil? (second avg-subopt))
        rank (map (fn [[wt & muts]] ;rank each individual sequence
                    {:rank (-> (remove #(< % wt) muts)
                               count
                               inc)
                     :num (count muts)})
                  data)
        neutrality (map (fn [[wt & muts]] ;;neut of each seq
                          {:wt wt :mut (mean muts)})
                        data)
        robustness (map (fn [m] (> (m :wt) (m :mut))) neutrality)]
    {:wt (let [wt (first (transpose data))]
           {:mean (-> wt mean) :sd (sd wt)})
     :muts (let [muts (map rest data)]
             {:mean (-> muts flatten mean)
              :sd (-> (map variance muts) mean Math/sqrt)})
     :rank rank
     :neutrality neutrality
     :robust? robustness}))

(defn subopt-neutrality-seq
  "Find the subopt overlap of a seq and all its 1-mutant
   neighbors. the seq (s) and structure (st) and number of suboptimal
   structures considered (n) must be given. Returns the neutrality
   <1-d/L>."
  
  [s st n]
  (->> (subopt-overlap-neighbors s st :nsubopt n)
       (map mean ) ;subopt-overlap
       mean)) ;neutrality

(defn subopt-robustness
  "Takes an input sto and estimates the significance of the robustness
   of the sequence. Takes a sequence and the consensus structure and
   generates n sequences with similar structure. Then finds the
   neutrality of each inverse-folded sequence. Returns the average
   %overlap-between-cons-and-suboptimal-structure for each
   sequence (cons wt muts)"

  [sto n & {:keys [ncore nsubopt invfile-ext distfun bp]
            :or {ncore 5
                 nsubopt 1000
                 invfile-ext ".inv.clj"
                 distfun subopt-overlap-neighbors
                 bp true}}]
  (let [ ;sto "/home/kitia/bin/gaisr/trainset2/pos/RF00555-seed.1.sto"
        inv-sto (read-clj (fs/replace-type sto invfile-ext))
        {l :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))
        ]
    [sto
     (map (fn [[nm s]]
            (when (get inv-sto nm)
              (let [[s st cons-keys] (degap-conskeys s cons)
                    inv-seq (remove (fn [iseq]
                                      (let [[iseq target target-keys] (degap-conskeys iseq st)]
                                        (empty? target-keys)))
                                    (take n (get inv-sto nm))) ;vector of n
                                        ;inverse-folded seqs
                    
                    neut (map (fn [x]
                                (subopt-overlap-neighbors x st
                                                          :ncore ncore
                                                          :nsubopt nsubopt))
                              (conj inv-seq s))]
                ;;average %overlap for each wt and inv-fold seq
                neut)))
          l)]))

(defn generic-robustness-sto
  "Takes an input sto and estimates the significance of the robustness
   of the sequence. Takes a sequence and the consensus structure and
   generates n sequences with similar structure. Then finds the
   neutrality of each inverse-folded sequence. Returns the average
   %overlap-between-cons-and-suboptimal-structure for each
   sequence (cons wt muts)"

  [sto n & {:keys [invfile-ext distfn]
            :or {invfile-ext ".inv.clj"
                 distfn subopt-seq}}]
  (let [ ;sto "/home/kitia/bin/gaisr/trainset2/pos/RF00555-seed.1.sto"
        inv-sto (read-clj (fs/replace-type sto invfile-ext))
        {l :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))
        ]
    [sto
     (map (fn [[nm s]]
            (when (get inv-sto nm)
              (let [[s st conskeys] (degap-conskeys s cons)
                    inv-seq (->> (get inv-sto nm)
                                 (remove (fn [iseq]
                                           (let [[iseq target target-keys]
                                                 (degap-conskeys iseq st)]
                                             (empty? target-keys))))
                                 (take n))]
                (map (fn [s]
                       (as-> [s st conskeys] state
                             (build-mut-neighbors state)
                             (calculate-dist state distfn)
                             (mean state)))
                     (conj inv-seq s))
                )))
          l)]))

;;;---------------------------------------------------
;;;driver functions which target specific datasets and use functions



(def ^{:private true} banner
  (let [parse (fn [s] (-> (str/split #" " s) vec))]
    [["-f" "--file" "REQUIRED. file(s) to check neutrality for"
      :parse-fn parse ;create list of files
      :default nil]
     ["-inv" "--invfile-ext" "extension of inverse folded seqs file"
      :default ".inv.clj"]
     ["-o" "--outfile" "REQUIRED. file to write to" :default nil]
     ["-di" "--dir" "dir in which files are located"
      :default nil]
     ["-dfn" "--distfn" "distance function"
      :parse-fn #(->> % symbol (ns-resolve 'robustness))
      :default subopt-overlap-neighbors]
     ["-d" "--debug" :default nil :flag true]
     ["-n" "--nseqs" "number of inverse seqs to create"
      :parse-fn #(Integer/parseInt %) :default 100]
     ["-nc" "--ncore" "number of cores to use"
      :parse-fn #(Integer/parseInt %) :default 2]
     ["-h" "--help" "usage" :default nil :flag true]]))
  
(defn main-subopt-robustness
  "Driver function for subopt-robustness. Takes input sto file and n
   inverse-folded structures. Compares the average
   value (subopt-overlap) of all seqs in the file for the wt seqs to
   the average value of all inverse-folded seqs. A seq in the sto is
   defined as robust when the wt average suboptimal overlap is greater
   than the average average-suboptimal-overlap of all inverse-folded
   seqs. The wt ranking defines the significance."
  
  [& args]
  (let [[opts _ usage] (apply cli args banner)
        ofile (opts :outfile) ;storage location
        infiles (filter #(fs/exists? (fs/replace-type % (opts :invfile-ext)))
                        (if (opts :dir)
                          (map #(fs/join (opts :dir) %) (opts :file))
                          (opts :file)))]
    (cond
     (opts :help) (print usage)
     (nil? args) (print usage)
     (nil? ofile) (println "requires outfile")
     (nil? (opts :file)) (println "requires in files")
     :else
     (doseq [instos (partition-all 2 infiles) 
             :let [cur (doall
                        (map (fn [insto]
                               [(keyword insto)
                                ;;list-of-lists average subopt overlap of 1-mut structures (neutrality)
                                (-> insto
                                    (generic-robustness-sto (opts :nseqs)
                                                            :invfile-ext (opts :invfile-ext)
                                                            :distfn (opts :distfn))
                                    )
                                #_(-> (subopt-robustness insto
                                                       (opts :nseqs)
                                                       :ncore (opts :ncore)
                                                       :invfile-ext (opts :invfile-ext)
                                                       :distfun (opts :distfn)) 
                                    subopt-robustness-summary)])
                             instos))
                   data (if (and (fs/exists? ofile)
                                 (not (fs/empty? ofile)))
                          (doall (concat (read-clj ofile) cur))
                          cur)]]
       (io/with-out-writer ofile
         (println ";;;generated using main-subopt-robustness. Estimate of the significance of the wild-type sto compared to the inverse folded version.")
         (prn data)))
     )))

#_(defn main-subopt-significance
  "Estimates the significance of the suboptimal overlaps seen when
   comparing the wild-type sto to the dinucleotide shuffed
   version. compares against 100 other stos. Since the process takes a
   long time to complete, the method constantly writes out to a file
   ofile.

   Memory usage could be high when it reads in all data and concats.

   Outputs data in the vector-map-of-list-of-lists-of maps

   analyze using: (reduce #(assoc %1 (first %2) (avg-overlap (second
   %2))) {} (read-string (slurp ofile)))"

  [& args]
  (let [[opts _ usage] (apply cli args banner)
        ofile (opts :file) ;storage location
        fdir (str homedir "/bin/gaisr/trainset2/pos/")
        ]
    (cond
     (opts :help) (print usage)
     (nil? args) (print usage)
     (not (fs/exists? ofile)) (println "requires outfile") 
     :else
     (doseq [instos (->> (remaining-files #(not (re-find #"\.3\.sto" %)) (fs/listdir fdir) (opts :ignore))
                         (partition-all 2 ) ;group into manageable chuncks
                         )]
       (let [cur (doall
                  (map (fn [insto]
                         [insto
                        ;;compare wild type sto against nsample shuffled versions
                          (avg-overlap
                           (subopt-significance (str fdir insto) (opts :ncores)
                                                (opts :nseqs) :ncores (opts :ncore)))])
                       instos))
             data (if (fs/exists? ofile) (concat (read-clj ofile) cur) cur)] ;add new data to existing
         (io/with-out-writer ofile
           (println ";;;generated using main-subopt-significance. Estimate of the significance of the wild-type sto compared to the dinucloetide shuffled version.")
           (prn data)) ;write to file
        )))))
;;;----------------------------------------------------








