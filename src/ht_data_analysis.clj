(ns ht-data-analysis
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [clojure-csv.core :as csv]
            [clojure.set :as sets]
            [clojure.core.reducers :as r]
            [clojure.pprint :as pp]
            [iota])
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        [edu.bc.bio.sequtils.snippets-files :only (read-aln)]
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.fold-ops
        edu.bc.bio.sequtils.tools))

(def t7 "TAATACGACTCACTATAgg")
(def prime5-const "TGCGTAACGTACACT")
(def variable-region "GGGATCGCTGAATTAGAGATCGGCGTCCTT");stem
(def prime3-const "TCATTCTATATACTTTGGAGTTTTAAA");after stem
(def cds "ATGTCTCTAAGTACT") ;pseduonoted partner
  
(def ref-seq
  (str prime5-const variable-region prime3-const cds))

(def parasite (str t7 prime5-const "AA" cds))

(def bob
  (str "TAATACGACTCACTATAGG"
       "TGCGTAACGTACACT"
       "TCCTTCGCTTATTCGGAGTAGATCACGTGA"
       "TCATTCTGTATGCTTTGGAGTTTTAAAATGTCTCTAAGTACT"))

(def S15-r1 "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.fastq")

(def S15-r2 "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R2_001.fastq")

;;;the files contain the pairs already in the same order for each of
;;;the fasta files
#_(every? true? 
          (apply map #(= (re-find #"\w+" (first %1))
                         (re-find #"\w+" (first %2))) 
                 (map read-fastq [S15-r1 S15-r2])))

(defn mutant-neighbor
  "Takes a string s and finds the 3L 1-mer mutants. String can only
   contain letters A, C, G, T. with-names returns a vector with the
   mutant name as well"
  
  [s & {:keys [with-names]
        :or {with-names false}}]
  (let [s (.toUpperCase s)]
    (for [i (range (count s))
          r [\A \C \G \T]
          :let [prefix (subs s 0 i)
                b (.charAt s i)
                suffix (subs s (inc i))
                name (str b i r)
                mut (str prefix r suffix)]
          :when  (not= r b)]
      mut)))

(defn n-mutant-neighbor [s n]
  (loop [n (range n)
         nmut [s]] 
    (if (seq n)
      (recur (rest n)
             (distinct (mapcat mutant-neighbor nmut)))
      (remove #(= s %) (distinct nmut)))))

(defn read-fastq [fastq]
  (partition-all 4 (io/read-lines fastq)))

(defn <=_< [a b c]
  (and (<= a b) (< b c)))

(defn pos->pseudo-num [pos rname]
  (+ (* 250375000
        (->> rname
             (str/tail 1 )
             (Integer/parseInt)
             dec))
     (Integer/parseInt pos)))

(defn pseudo-num->bob-num [pseudo-num]
  (quot pseudo-num 50))

(defn hts-name-parts
  "Breaks the name line from a fastq file into the components and
  stores it in a map"

  [nm]
  (let [[x y] (str/split #" " nm)
        xx (interleave [:inst-name :run-id
                        :flowcell-id :flowcell-lane
                        :flowcell-tile-num
                        :xcoord :ycoord]
                       (str/split #":" x))
        yy (when-not (nil? y)
             (interleave [:pair :filter?
                          :control-bits :index-seq]
                         (str/split #":" y)))]
    (reduce (fn [m [k v]]
              (assoc m k v))
            {} (partition-all 2 (concat xx yy)))))

(defn read-aln2 [f & {:keys [info]
                      :or {info :data}}]
  (let [p (read-seqs f)
        x (split-at (quot (count p) 2) p)
        [l r] [(first x) (rest (second x))]]
    (map (fn [l r]
           (let [[nm l] (str/split #"\s+" l)
                 s (str l (second (str/split #" +" r)))]
             (case info
               :data s
               :names nm
               :both [nm s])))
         l r)))
;;;===========================================================================


(defn forward?
  "Determines if the strand is forward or reverse by finding
  the 5'-constant region in the strand"

  [s]
  (re-find (re-pattern prime5-const) s))

(defn var-seeker
  "Finds the variable region in the strand by forcing the string to
  fit the construct of 5'-constant+N30+3'-constant type
  structure. returns only the variable region"

  [s]
  (-> (str prime5-const 
           "(\\w{30})"
           prime3-const)
      re-pattern
      (re-find s)
      second))



(defn get-seqs [infastq1 infastq2]
  (map (fn [r1 r2]
         (let [[nm1 inseq1] r1
               [nm2 inseq2] r2]
           (if (forward? inseq1)
             [nm1 inseq1]
             [nm2 inseq2])))
       (partition-all 4 (io/read-lines infastq1))
       (partition-all 4 (io/read-lines infastq2))))

(defn get-freqs []
  (r/fold (fn ([] {})
            ([l r] (merge-with conj l r)))
          (fn [m [nm fseq]]
            (let [nm-map (select-keys (hts-name-parts nm)
                                      [:flowcell-lane
                                       :xcoord :ycoord :pair])
                  k [(var-seeker fseq) (nm-map :flowcell-lane)]
                  nm-list (get m k [])]
              (if (first k) (assoc m k (conj nm-list nm-map)) m)))
          (get-seqs S15-r1 S15-r2)))

(comment
  (def var-region4-freqs (group-by (fn [[[n30 lane] _]]
                                     lane)
                                   (get-freqs)))
  (def round4-n30-seqs (map ffirst (var-region4-freqs "1")))
  (def round9-n30-seqs (map ffirst (var-region4-freqs "2")))

  (def bob-cluster
    (future 
      (let [bobs (robustness/mutant-neighbor variable-region)]
        (doall
         (io/with-out-writer "/home/peis/TMP/bob-cluster.fna"
           (doseq [[r1 r2] (map vector (read-fastq S15-r1) (read-fastq S15-r2))
                   bob bobs
                   :let [iseq (second r1)
                         jseq (second r2)
                         re (re-pattern bob)
                         x (re-find re iseq)
                         y (re-find re jseq)]
                   :when (or x y)]
             (cond 
              x (dorun (map println r1))
              y (dorun (map println r2)))))))))

;;;frequency of the number of times an N30 region appears
;;;the map is k=#appearances v=#species appearing k times
;;;last number shows the number of nils
  (frequencies (map #(-> % second count) (var-region4-freqs "1")))
;;;round 4 {1 3107618, 2 12011, 3 46, 5 1, 2358643 1}
;;;round 9 {1 3114378, 2 12479, 3 35, 2348177 1}

  (map #(relative-entropy (probs (frequencies %))
                          {\A 0.25 \C 0.25 \G 0.25 \T 0.25})
       (transpose (take 10 round4-n30-seqs)))

  (defn create-pseudoseq
    "make up seqs for a pseudogene and prints to file f ie
  /home/peis/TMP/bob-pseudoseqs3.fna"

    [f]
    (let [bob-nums (atom 0)
          fun (fn [n L]
                (let [P (probs 1 [\A \C \G \T])]
                  (->> #(markov-step P)
                       repeatedly
                       (partition-all L )
                       (map #(apply str %) ) 
                       distinct
                       (take n)
                       )))]
      (io/with-out-writer f 
        (doseq [i (fun 10000 30)
                j (conj (n-mutant-neighbor i 2) i)]
          (do (println (str ">bob-" @bob-nums)) 
              (println (str (str/tail 10 prime5-const)
                            j
                            (str/take 10 prime3-const)))
              (swap! bob-nums inc))))))

;;;pseudogenome
  (let [n 8 ;number files to split gnome into
        pseudogene (read-seqs "/home/peis/TMP/bob-pseudoseqs3.fna")
        c (/ 40059999 n)                      ;last name + 1
        ] 
    (loop [i pseudogene
           j (atom 0)
           outfile (for [i (range 1 (inc n))]
                     [(str "/home/peis/TMP/bob-pseudogene3-" i ".fna") ;filename
                      (str ">bob|pseudogene3." i)])] ;gene name
      (when (seq i)
        (io/with-out-writer (-> outfile first first)
          (println (-> outfile first second))
          (doseq [k i
                  :while (< @j c)]
            (do (print k)
                (swap! j inc))))
        (recur (drop @j i)
               (atom 0)
               (rest outfile)))))

  (defn sam-reader [f]
    (let [sam (iota/seq f)
          aligned (r/remove #(re-find #"^\@\w+" %) sam)
          ]
      (r/map #(str/split #"\t" %) aligned)))

  (defn sam-string [sam-vec]
    (apply hash-map
           (interleave [:qname :flag :rname :pos :mapq :cigar :rnext :pnext :tlen :seq :qual]
                       sam-vec)))

  (defn cigar-filtered

    [f]
    (let [sam-aligned (sam-reader f)]
      (r/filter (fn [line]
                  (let [{cigar :cigar} (sam-string line)
                        cigar (re-find #"(\d+)S\d+M(\d+)S" (nth line 5))
                        soft5 (if-not (nil? cigar) (Integer/parseInt (cigar 1)) 0)
                        soft3 (if-not (nil? cigar) (Integer/parseInt (cigar 2)) 0)] 
                    (and (>= soft5 10)
                         (>= soft3 10))))
                sam-aligned)))

;;;filter according to cigar
  (def foo
    (let [bins (partition 2 1 (map #(* 4006 50 %) (range 10001)))
          
          find-bin (fn [x] 
                     (-> (filter #(<=_< (first %) x (second %)) bins) first vec))]
      (r/fold (fn ([] {})
                ([l r] (merge-with + l r)))
              (fn ([] {})
                ([m sam-line]
                   (let [{pos :pos rname :rname qname :qname} (sam-string sam-line)
                         pos (pos->pseudo-num pos rname)
                         bin (find-bin pos)
                         cur-val (get m bin 0)]
                     (assoc m bin (inc cur-val)))))
                              ;;;filter according to cigar
              (sam-reader "/home/peis/TMP/test.sort.sam"))))
  (def bar
    (map (fn [[pos-vec seed-seq]] 
           (let [bins (bins-of-interest pos-vec)]
             (r/fold (fn ([] [])
                       ([l r] (concat l r)))
                     (fn ([] [])
                       ([l line]
                          (let [{qname :qname seq :seq} (sam-string line)
                                lane (->> (str/split #":" qname) (drop 3) first (Integer/parseInt))
                                alignment (smith-waterman/sw seed-seq seq :global true
                                                             :match-weight 1
                                                             :gap-weight -2)
                                score (ffirst alignment)]
                            (if (pos? score)
                              (conj l 
                                    [qname
                                     (cond
                                      (= lane 1) ":round4"
                                      (= lane 2) ":round9")
                                     (str score)
                                     (-> alignment first third)])
                              l))))
                     bins)))
         (map vector
              (map first (filter #(>= (second %) 10000) foo))
              ["AACGTACACTCCAAAGTTTAAAGATGAATGGTTGTGGCAGTCATTCTATA"
               "AACGTACACTAGATGTCTCGTCTAGTCTCTTTTATCGAGTTCATTCTATA"
               "AACGTACACTAGACGTCTCTGTCCTTGTTGTTACTGCGTGTCATTCTATA"])))

;;;separate out a list of hits
  (pmap (fn [[pos-vec seed-seq]] 
          (let [bins (bins-of-interest pos-vec)]
            (io/with-out-writer (str "/home/peis/TMP/bins-of-interest-" (first pos-vec) ".txt")
              (println (apply str (apply str (repeat 14 "ref")) "\t:round0\t0\t" seed-seq))
              (doseq [line bins
                      :let [{qname :qname seq :seq} (sam-string line)
                            lane (->> (str/split #":" qname) (drop 3) first (Integer/parseInt))
                            alignment (smith-waterman/sw seed-seq seq :global true
                                                         :match-weight 1
                                                         :gap-weight -2)
                            score (ffirst alignment)]
                      :when (pos? score)]
                (-> [[qname
                      (cond
                       (= lane 1) ":round4"
                       (= lane 2) ":round9")
                      (str score)
                      (-> alignment first third)]]
                    (Csv/write-csv :delimiter "\t") 
                    print)))))
        (map vector
             (map first (filter #(>= (second %) 10000) foo))
             ["AACGTACACTCCAAAGTTTAAAGATGAATGGTTGTGGCAGTCATTCTATA"
              "AACGTACACTAGATGTCTCGTCTAGTCTCTTTTATCGAGTTCATTCTATA"
              "AACGTACACTAGACGTCTCTGTCCTTGTTGTTACTGCGTGTCATTCTATA"]))

;;;make fastas from the text files containing score output after
;;;filtering for seqs that align with positive score to a seed seq
  (let [f (io/read-lines "/home/peis/TMP/bins-of-interest-881520300.txt")
        c (atom 0)]
    (io/with-out-writer "/home/peis/TMP/round4-881520300.fna"
      (doseq [i f
              :let [[qname round score seq] (str/split #"\t" i)]
              :when (= round ":round4")]
        (do (println (str ">" qname "-" @c) round)
            (println (str/replace-re #"-" "" seq))
            (swap! c inc)))))

;;;removal of parasites after examining common sequences appearing in
;;;the clustalw alignment we noticed that a lot of sequences start
;;;with the 5' primer TACACT but don't end with the 3' primer. This
;;;lead to the discovery of a lot of parasites which contain the
;;;ATGTCT motif
  (let [fs (fs/re-directory-files "/home/peis/TMP" "round*fna")]
    (doseq [f fs]
      (io/with-out-writer (fs/replace-type f ".no-parasite.fna")
        (doseq [i (->> (read-seqs f :info :both)
                       (remove #(re-find #"TACACT.{1,10}ATGTCT" (second %))))                         
                j i]
          (println j)))))

;;;diversity of files
  (let [fnas (fs/re-directory-files "/home/peis/TMP" #"round.*.fna") ] 
    (sort-by first 
             (for [fna fnas
                   :let [f (read-seqs fna)]]
               [fna (shannon-entropy (probs 1 f))])))

  (defn unblock-clustalw [f]
    (let [inseqs (read-aln2 f :info :both)]
      (io/with-out-writer (fs/replace-type f ".unblock.aln")
        (println "CLUSTAL 2.1 multiple sequence alignment\n\n")
        (doseq [i inseqs] (println (str/join "\t" i))))))

;;;output the first 2000 N30 sequences that appear exactly 3 times
;;;from round 9. the name is simply the xy coord
  (io/with-out-writer "/home/peis/TMP/foo9.txt"
    (doseq [i (->> (var-region4-freqs "2")
                   (filter #(= (-> % second count) 3) )
                   (take 2000))
            :let [m (-> i second first)
                  outname (str ">" (m :xcoord) ":" (m :ycoord))]]
      (do (println outname)
          (println (ffirst i)))))

;;;difference between variable region of paired sequences
  (def pairs
    (frequencies
     (apply map #(levenshtein (str/drop 6 (second %1))
                              (reverse-compliment (str/drop 6 (second %2))))
            (map read-fastq [S15-r1 S15-r2]))))


;;;filtering S15 mate pair 1 for sequences of potentially valid
;;;lenght. the total length of S15 regulatory RNA is 57 so we in
;;;essence allow a 10nt insertion/deletion in the length. this filter
;;;should remove any parasites in the data
  (let [fastq (read-fastq S15-r1)
        outfile "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.filtered.fastq"]
    (io/with-out-writer outfile
      (doseq[i (filter #(re-find #"ACT.{45,65}ATGTCT" (second %))
                       fastq)
             j i]
        (println j))))

;;;from the fastq file create a new file that contains only the
;;;relevant part of the sequence
  (let [fastq (iota/seq "/home/peis/S15SELEXHTSdata/index46_TCCCGA_combined_R1_001.filtered.fastq")]
    (io/with-out-writer "/home/peis/TMP/foo.txt"
      (doseq [i (r/fold (fn ([] [])
                          ([l r] (concat l r)))
                        (fn ([] [])
                          ([V x]
                             (let [y (re-find #"ACT.{45,65}ATG" x)]
                               (if y (conj V y) V))))
                        fastq)]
        (println i))))

;;;attempt to process 1000 lines of seq and compare using levenshtein
;;;distance. this method is infeasible due to the computational
;;;restraints
  (let [fastq (iota/seq "/home/peis/TMP/foo.txt")
        d (loop [i (into [] (r/take 1000 fastq))
                 M {}]
            (if-not (empty? i)
              (let [x (first i)
                    y (re-find #"ACT.{20,40}TCATT" x)
                    ;;find parts of the seq which match the regex
                    ;;pattern and then compare the given seq to the ith
                    ;;seq and store if it's within a distance of 10
                    data (r/fold 25
                                 concat
                                 (fn ([] [])
                                   ([m all]
                                      (let [z (re-find #"ACT.{20,40}TCATT" all)]
                                        (if (<= (levenshtein y z) 10)
                                          (conj m z)
                                          m))))
                                 i)]
                (recur (vec (remove (set data) (rest i)));remove found seqs
                       (if (>= (count data) 2);>=2 seqs to be unique cluster
                         (assoc M y data)
                         M)))
              M))]
    ;;print data to a file
    (io/with-out-writer "/home/peis/TMP/foo2.txt"
      (println d)))

;;;make fasta file or seqs and name with arbitrary number. the seqs in
;;;the new file are used as blast query sequences. they only have
;;;ACT+N30variable+TCATT construct
  (let [foo (io/read-lines "/home/peis/TMP/foo.txt")
        c (atom 0)]
    (io/with-out-writer "/home/peis/TMP/myblasttest.fna"
      (doseq [[i s] (->> foo 
                         (map-indexed (fn [i s] [i (re-find #"ACT.{20,40}TCATT" s)]))
                         (remove #(nil? (second %))))
              :while (< @c 1000)]
        (println (str ">" i))
        (println s)
        (swap! c inc))))

;;;blast the query file against the database of HTS reads
  (let [blast (blastn "/home/peis/TMP/myblasttest.fna"
                      :out (str "/home/peis/TMP/myblasttest." i ".out")
                      :blastdb (str "/home/peis/S15SELEXHTSblastdb/s15SELEXHTSfiltered." i)
                      :word-size 5
                      :misc "-reward 1 -penalty -1 -gapopen 3 -gapextend 2 -num_threads 15 -dbsize 590400800"
                      :fmt "10 qseqid qstart qend evalue nident sseqid sstart send")
        blast-result (io/read-lines "/home/peis/TMP/myblasttest.out")
        ]
    (io/with-out-writer "/home/peis/TMP/myblasttest.filtered.out"
      (doseq [i (filter (fn [result]
                          (let [[qseid qstart qend evalue nident sseqid start send]
                                (first (csv/parse-csv result))]
                            (>= (Integer/parseInt nident) 28)))
                        blast-result)]
        (println i))))

;;;remade the blastdb by splitting into 10 subparts. now blast occurs
;;;for the seq across 10 dbs. the results are unioned together to
;;;create the final output. Effective dbsize is added to reflect the
;;;size of full database so that evalues do not seem more significant
  (future
    (time
     (let [outfiles
           (doall
            (for [i (range 10)]
              (blastn "/home/peis/TMP/myblasttest.fna" :out (str "/home/peis/TMP/myblasttest." i ".out") :blastdb (str "/home/peis/S15SELEXHTSblastdb/s15SELEXHTSfiltered." i) :word-size 5 :misc "-reward 1 -penalty -1 -gapopen 3 -gapextend 2 -num_threads 15 -dbsize 590400800" :fmt "10 qseqid qstart qend evalue nident sseqid sstart send")))]
       (io/with-out-writer "/home/peis/TMP/myblasttest2.out"
         (doseq [i outfiles
                 j (io/read-lines i)]
           (println j)))))))

(defn combin-blast

  [directory prefix outcsv]
  (let [files (fs/re-directory-files directory
                                     (re-pattern (str prefix "\\d+.out")))]
    (io/with-out-writer outcsv
      (doseq [i (mapcat io/read-lines files)]
        (println i)))))

(defn blast-hts-wrapper

  [in out blastdb]
  (blastn in :out out
          :blastdb blastdb
          :word-size 5
          :misc "-reward 1 -penalty -1 -gapopen 3 -gapextend 2 -num_threads 25 -dbsize 590400800"
          :fmt "10 qseqid qstart qend evalue nident sseqid sstart send"))

(defn filter-blast-results

  [blast-result & {:keys [thr] :or {thr 28}}]
  (filter (fn [result]
            (let [[qseid qstart qend evalue nident sseqid start send]
                  (str/split #"," result)]
              (>= (Integer/parseInt nident) thr)))
          (io/read-lines blast-result)))

(defn main-blast-driver

  [& args]
  (let [[opts _ usage] (cli args
                            ["-f" "--file" "queries to blast" :default nil]
                            ["-o" "--outfile" "blast output" :default nil]
                            ["-db" "--blastdb" "blastdb to search" :default nil])
        in (opts :file)
        out (opts :outfile)
        blastdb (opts :blastdb)]
    (if (or (nil? args)
              (opts :help)
              (some nil? args))
      (print usage)
      (do (blast-hts-wrapper in out blastdb)
          (io/with-out-writer (fs/replace-type out ".filtered.out")
            (doseq [i (filter-blast-results)]))))))

(defn pbs-blast-driver

  [& args]
  (let [[opts _ usage] (cli args
                            ["-f" "--file" "queries to blast" :default nil]
                            ["-o" "--outfile" "blast output prefix" :default nil]
                            ["-db" "--blastdb" "blastdb to search" :default nil]
                            ["-di" "--pdir" "dir in which to produce pbs files"
                             :default nil]
                            ["-p" "--pbs" "prefix for pbs files to produce"
                             :default nil])
        in (opts :file)
        out (opts :outfile)
        blastdb (opts :blastdb)]
    (if (or (nil? args)
              (opts :help)
              (some nil? args))
      (print usage)
      (doseq [i (range 10)]
       (let [template {:shell "#!/bin/bash"
                       :resource (str "#PBS -l mem=5gb,nodes=1:ppn=16" ",walltime=100:00:00")
                       :working-dir "#PBS -d /home/peis/bin/gaisr/"
                       :module "module add ncbiblast"
                       :command (str
                                 "lein run -m ht-data-analysis/main-blast-driver"
                                 " -f " "\"" in "\""
                                 " -o " (str out "." i ".out")
                                 " -db " (str blastdb "." i)
                                 )}
             outpbs (str (fs/join (opts :pdir) (opts :pbs)) "." i ".pbs")]
         (println opts)
         (io/with-out-writer outpbs
           (println (template :shell))
           (println (template :resource))
           (println (template :working-dir))
           (println (template :module) "\n")
           (println (template :command)))
         outpbs)))))

;;;round4
;;;"0.12884196903489223,0.0864283103545948,0.08472256052753267,0.06842485783825393,0.07770744199272558,0.0704715728434399,0.06465902210503963,0.05802716246662663,0.05628170153091212,0.05188434425461899,0.05218078594389408,0.05127928523073287,0.05096776039074802,0.04714972645197674,0.05041856378523632,0.04789509265751778,0.04591946697771544,0.04581962196081607,0.045535841155522144,0.042609712065177247,0.04277477532893913,0.04234552418236616,0.04123868679683079,0.03654822678546027,0.03653470224714478,0.03335048624958892,0.02743214113274514,0.025968903040796174,0.027990975906299106,0.02430884248707646"


;;;round9
;;;"0.1288781804324984,0.08600645802736594,0.08499500423555527,0.06884838430878303,0.07797702107165977,0.07074144517907104,0.06471124292304158,0.05879414032082503,0.05553673786924779,0.051905318792766886,0.05269249255461443,0.05100304499197318,0.0509856171367582,0.047407446335724426,0.050893133752294484,0.048086870055667225,0.04591902514987413,0.04542840915439787,0.04538885931554702,0.042629145167125335,0.04280327884163082,0.04247415019397602,0.041341748393219344,0.03684887253126677,0.03645624038032161,0.03314059262066983,0.02719336571841155,0.02622975012873685,0.028008647415719373,0.024261784639797654"

;;;jsd comparing round 4 to round 9
;;;"1.5838557975681345E-6,1.5819957371086086E-6,8.01549332228829E-7,5.8944928488908436E-6,1.950699645866924E-6,9.657815317308668E-7,1.4622032928925913E-6,1.280304688027644E-6,5.247252771242413E-6,1.168632599337808E-6,3.919911398426741E-6,1.1887391760298223E-6,5.562152249357434E-6,1.5126484241453865E-6,4.132107187430624E-6,6.74198074636221E-6,3.1095422971483765E-6,2.0884556368633332E-6,1.2504472252076895E-6,5.785441678858485E-6,4.836221526700012E-6,1.0784143811082086E-5,9.897764456182008E-6,7.986904609809452E-6,2.6397845078018737E-6,2.14540363846323E-6,4.442589476412326E-6,4.6044841497616754E-6,4.401917250023969E-6,3.040629339682118E-6"

(comment
  EAS139	the unique instrument name
  136	the run id
  FC706VJ	the flowcell id
  2	flowcell lane
  2104	tile number within the flowcell lane
  15343	'x'-coordinate of the cluster within the tile
  197393	'y'-coordinate of the cluster within the tile
  1	the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
  Y	Y if the read fails filter (read is bad), N otherwise
  18	0 when none of the control bits are on, otherwise it is an even number
  ATCACG	index sequence

  )


